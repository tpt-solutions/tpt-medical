//! Binary segmentation masks on the CT voxel grid.

use tpt_med_dicom::DicomSeries;
use tpt_med_geometry::Vec3;

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
    /// [`HounsfieldMapper::DEFAULT_BONE_THRESHOLD_HU`] when a data-driven
    /// default is acceptable.
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
}
