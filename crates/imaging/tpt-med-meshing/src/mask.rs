//! Segmentation masks on the image voxel grid, from DICOM or NIfTI sources.

use tpt_med_dicom::DicomSeries;
use tpt_med_geometry::{ras_to_lps, Vec3};
use tpt_med_nifti::NiftiVolume;

/// A segmentation on the CT voxel grid: solid/empty per voxel, with the HU
/// values and patient-space geometry needed for meshing.
#[derive(Debug, Clone)]
pub struct SegmentationMask {
    /// Volume dimensions `(nx, ny, nz)` in voxels.
    pub dims: (usize, usize, usize),
    /// Patient-space centre of voxel (0,0,0) (mm).
    pub origin: Vec3,
    /// Image basis vectors (unit): row direction (+x voxels).
    pub row_dir: Vec3,
    /// Image basis vectors (unit): column direction (+y voxels).
    pub col_dir: Vec3,
    /// Image basis vectors (unit): slice direction (+z voxels).
    pub slice_dir: Vec3,
    /// Voxel pitch along `(row, col, slice)` axes in mm.
    pub spacing: (f64, f64, f64),
    /// Solid flag per voxel, index order `(z * ny + y) * nx + x`.
    pub voxels: Vec<bool>,
    /// HU per voxel, same indexing as [`Self::voxels`].
    pub hu: Vec<f64>,
}

impl SegmentationMask {
    /// Thresholds a CT series into a bone mask: voxels with
    /// `HU >= min_hu` are solid. Uses
    /// [`HounsfieldMapper::DEFAULT_BONE_THRESHOLD_HU`](tpt_med_dicom::HounsfieldMapper::DEFAULT_BONE_THRESHOLD_HU)
    /// when a data-driven default is acceptable.
    pub fn threshold_hu(series: &DicomSeries, min_hu: f64) -> Self {
        let (nx, ny, nz) = series.dims();
        let mut voxels = vec![false; nx * ny * nz];
        let mut hu = vec![f64::NAN; nx * ny * nz];

        let first = &series.slices[0];
        let frame = first.frame();
        let row_dir = frame.row_dir.normalize();
        let col_dir = frame.col_dir.normalize();
        let slice_dir = frame.normal();

        // z pitch: mean spacing between consecutive slice positions.
        let z_pitch = if series.slices.len() > 1 {
            let mut sum = 0.0;
            for w in series.slices.windows(2) {
                sum += (w[1].position - w[0].position).dot(slice_dir).abs();
            }
            sum / (series.slices.len() - 1) as f64
        } else {
            series.slice_thickness
        };

        for (z, slice) in series.slices.iter().enumerate() {
            for (idx, &stored) in slice.pixel_data.iter().enumerate() {
                let h = stored as f64 * slice.rescale_slope + slice.rescale_intercept;
                let v = z * nx * ny + idx;
                hu[v] = h;
                voxels[v] = h >= min_hu;
            }
        }

        Self {
            dims: (nx, ny, nz),
            origin: first.position,
            row_dir,
            col_dir,
            slice_dir,
            spacing: (series.pixel_spacing.1, series.pixel_spacing.0, z_pitch),
            voxels,
            hu,
        }
    }

    /// Thresholds a NIfTI volume into a bone mask: voxels whose value is
    /// `>= min_hu` are solid — the same rule and the same default
    /// (`HounsfieldMapper::DEFAULT_BONE_THRESHOLD_HU`) as
    /// [`Self::threshold_hu`], because a CT NIfTI export's `values` are HU
    /// when its writer baked DICOM's `RescaleSlope`/`RescaleIntercept` into
    /// `scl_slope`/`scl_inter`, as `dcm2niix` does.
    ///
    /// NIfTI geometry is **RAS**; origin and direction columns are converted
    /// through [`ras_to_lps`] so the mask lands in the same **LPS** patient
    /// frame [`Self::threshold_hu`] produces — a mask means the same thing
    /// regardless of which format it was read from. The source volume must
    /// be well-formed (as the parser guarantees): `values` covers exactly
    /// `dims.0 * dims.1 * dims.2` voxels in x-fastest order.
    pub fn threshold_nifti(volume: &NiftiVolume, min_hu: f64) -> Self {
        let (nx, ny, nz) = volume.dims;
        debug_assert_eq!(volume.values.len(), nx * ny * nz);
        let voxels: Vec<bool> = volume.values.iter().map(|&h| h >= min_hu).collect();

        // RAS → LPS flips x and y, so each direction column maps through
        // the same involution the origin does (normalised for parity with
        // the DICOM path's direction-cosine handling).
        let row_dir = ras_to_lps(volume.rotation.col(0)).normalize();
        let col_dir = ras_to_lps(volume.rotation.col(1)).normalize();
        let slice_dir = ras_to_lps(volume.rotation.col(2)).normalize();

        Self {
            dims: volume.dims,
            origin: ras_to_lps(volume.origin),
            row_dir,
            col_dir,
            slice_dir,
            spacing: volume.voxel_spacing,
            voxels,
            hu: volume.values.clone(),
        }
    }

    /// Number of solid voxels.
    pub fn solid_count(&self) -> usize {
        self.voxels.iter().filter(|&&v| v).count()
    }

    /// Flat index of a voxel, bounds-checked.
    pub fn index(&self, x: usize, y: usize, z: usize) -> Option<usize> {
        let (nx, ny, _) = self.dims;
        if x < nx && y < ny && z < self.dims.2 {
            Some((z * ny + y) * nx + x)
        } else {
            None
        }
    }

    /// Solid flag at (x, y, z); out-of-bounds is empty.
    pub fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        self.index(x, y, z).is_some_and(|i| self.voxels[i])
    }

    /// Patient-space centre of the voxel at (x, y, z).
    pub fn voxel_center(&self, x: usize, y: usize, z: usize) -> Vec3 {
        self.origin
            + self.row_dir * (x as f64 * self.spacing.0)
            + self.col_dir * (y as f64 * self.spacing.1)
            + self.slice_dir * (z as f64 * self.spacing.2)
    }

    /// Patient-space position of grid node `(i, j, k)` — the corner between
    /// voxels `(i-1..i, j-1..j, k-1..k)`.
    pub fn node_position(&self, i: usize, j: usize, k: usize) -> Vec3 {
        self.voxel_center(0, 0, 0)
            - (self.row_dir * self.spacing.0
                + self.col_dir * self.spacing.1
                + self.slice_dir * self.spacing.2)
                * 0.5
            + self.row_dir * (i as f64 * self.spacing.0)
            + self.col_dir * (j as f64 * self.spacing.1)
            + self.slice_dir * (k as f64 * self.spacing.2)
    }

    /// Mean HU over solid voxels (diagnostic helper).
    pub fn mean_solid_hu(&self) -> Option<f64> {
        let mut sum = 0.0;
        let mut n = 0usize;
        for (v, &solid) in self.voxels.iter().enumerate() {
            if solid {
                sum += self.hu[v];
                n += 1;
            }
        }
        (n > 0).then_some(sum / n as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_med_dicom::synthetic;

    #[test]
    fn threshold_segments_core_and_shell() {
        let series = synthetic::femur_phantom(16, 16, 8);
        let mask = SegmentationMask::threshold_hu(
            &series.parse().unwrap(),
            tpt_med_dicom::HounsfieldMapper::DEFAULT_BONE_THRESHOLD_HU,
        );
        assert_eq!(mask.dims, (16, 16, 8));
        // Bone (>=200 HU): cortical shell (700) + trabecular core (150 < 200
        // → NOT bone at the default threshold).
        let bone = mask.solid_count();
        assert!(bone > 0 && bone < 16 * 16 * 8, "bone={bone}");
        // Centre voxel (150 HU) is not bone at 200 HU threshold.
        assert!(!mask.is_solid(8, 8, 4));
        // Shell voxel at r≈5 from centre is bone.
        assert!(mask.is_solid(8 + 5, 8, 4));
        // Soft-tissue corner is not bone.
        assert!(!mask.is_solid(0, 0, 0));
        assert!(mask.mean_solid_hu().unwrap() >= 200.0);
    }

    #[test]
    fn node_and_voxel_positions() {
        let series = synthetic::femur_phantom(8, 8, 4);
        let mask = SegmentationMask::threshold_hu(&series.parse().unwrap(), 200.0);
        // Node (0,0,0) sits half a voxel below/left/behind the first voxel.
        let n0 = mask.node_position(0, 0, 0);
        let v0 = mask.voxel_center(0, 0, 0);
        assert!((v0.x - n0.x - 0.5).abs() < 1e-12);
        assert!((v0.y - n0.y - 0.5).abs() < 1e-12);
        assert!((v0.z - n0.z - 0.5).abs() < 1e-12);
        // Node grid advances by spacing.
        let n1 = mask.node_position(1, 1, 1);
        assert!((n1.x - n0.x - 1.0).abs() < 1e-12);
        assert!((n1.y - n0.y - 1.0).abs() < 1e-12);
        assert!((n1.z - n0.z - 1.0).abs() < 1e-12);
    }

    /// A 4x3x2 volume with spacing (1, 2, 3) mm and origin (10, 20, 30)
    /// RAS: the `i < 2` half is bone (800 HU), the rest soft tissue (50 HU).
    fn synthetic_bone_volume() -> NiftiVolume {
        let (nx, ny, nz) = (4usize, 3, 2);
        let mut values = vec![50.0; nx * ny * nz];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..2 {
                    values[i + nx * (j + ny * k)] = 800.0;
                }
            }
        }
        let bytes = tpt_med_nifti::synthetic::SyntheticNiftiBuilder::new(4, 3, 2, (1.0, 2.0, 3.0))
            .with_origin(Vec3::new(10.0, 20.0, 30.0))
            .with_values_f32(&values)
            .build();
        NiftiVolume::parse_bytes(&bytes).expect("synthetic NIfTI parses")
    }

    #[test]
    fn threshold_nifti_segments_and_converts_ras_to_lps() {
        let vol = synthetic_bone_volume();
        let mask = SegmentationMask::threshold_nifti(&vol, 200.0);
        assert_eq!(mask.dims, (4, 3, 2));
        // The i < 2 half is solid: 2 * 3 * 2 voxels, boundary included.
        assert_eq!(mask.solid_count(), 12);
        assert!(mask.is_solid(0, 0, 0));
        assert!(mask.is_solid(1, 0, 0));
        assert!(!mask.is_solid(2, 0, 0));
        assert!(!mask.is_solid(3, 2, 1));
        assert!(mask.mean_solid_hu().unwrap() >= 200.0);

        // RAS (10, 20, 30) → LPS (-10, -20, 30); the +x/+y columns flip
        // and +z does not — the frame `threshold_hu`'s DICOM path produces.
        assert!((mask.origin - Vec3::new(-10.0, -20.0, 30.0)).norm() < 1e-6);
        assert!((mask.row_dir - Vec3::new(-1.0, 0.0, 0.0)).norm() < 1e-9);
        assert!((mask.col_dir - Vec3::new(0.0, -1.0, 0.0)).norm() < 1e-9);
        assert!((mask.slice_dir - Vec3::new(0.0, 0.0, 1.0)).norm() < 1e-9);
        assert!((mask.spacing.0 - 1.0).abs() < 1e-6);
        assert!((mask.spacing.1 - 2.0).abs() < 1e-6);
        assert!((mask.spacing.2 - 3.0).abs() < 1e-6);

        // The mask's voxel centres are the volume's voxel positions in LPS.
        let p_ras = vol.voxel_position(1, 1, 1);
        let p_lps = mask.voxel_center(1, 1, 1);
        assert!((ras_to_lps(p_ras) - p_lps).norm() < 1e-6);
    }

    #[test]
    fn nifti_mask_meshes_end_to_end() {
        let vol = synthetic_bone_volume();
        let mask = SegmentationMask::threshold_nifti(&vol, 200.0);
        let mesh = crate::MedicalMesher::default()
            .voxels_to_hex_mesh(&mask)
            .expect("a NIfTI-sourced mask meshes like a DICOM one");
        assert_eq!(mesh.elements.len(), mask.solid_count());
    }
}
