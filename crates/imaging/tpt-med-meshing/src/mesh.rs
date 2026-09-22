//! Voxel grid → hexahedral element mesh conversion.

use std::collections::HashMap;

use crate::mask::SegmentationMask;
use crate::MeshError;
use tpt_med_dicom::{BoneRegion, HounsfieldMapper};
use tpt_med_geometry::Vec3;

/// Material assignment per element.
#[derive(Debug, Clone, Copy)]
pub struct ElementMaterial {
    /// Mean HU over the element's voxels.
    pub mean_hu: f64,
    /// Apparent density (g/cm³) from the HU correlation.
    pub density: f64,
    /// Young's modulus (MPa) from the region-specific density law.
    pub youngs_modulus: f64,
    /// Which bone-region law was applied.
    pub region: BoneRegion,
    /// Poisson's ratio assigned to the element.
    pub poissons_ratio: f64,
}

/// A structured hexahedral mesh: one element per solid voxel, shared corner
/// nodes, per-element HU-derived materials, in patient coordinates (mm).
#[derive(Debug, Clone)]
pub struct VoxelHexMesh {
    /// Node positions (mm, patient coordinates).
    pub nodes: Vec<Vec3>,
    /// Hex connectivity, 8 node indices per element in the order
    /// `(±x, ±y, ±z)` corners: `[000,100,110,010,001,101,111,011]` with
    /// bit order `(i,j,k)` relative to the voxel's minimum corner.
    pub elements: Vec<[u32; 8]>,
    /// Material properties per element (same order as [`Self::elements`]).
    pub materials: Vec<ElementMaterial>,
}

/// Voxel-to-hex mesher with tunable material assignment.
#[derive(Debug, Clone)]
pub struct MedicalMesher {
    /// Elements with density ≥ this (g/cm³) use the cortical modulus law;
    /// below it, the trabecular law. Default 1.3 g/cm³ (≈300 HU), a common
    /// CT-FEM heuristic.
    pub region_split_density: f64,
    /// Poisson's ratio for bone elements. Default 0.30.
    pub poissons_ratio: f64,
}

impl Default for MedicalMesher {
    fn default() -> Self {
        Self {
            region_split_density: 1.3,
            poissons_ratio: 0.30,
        }
    }
}

impl MedicalMesher {
    /// Converts a segmentation mask into a hex mesh. Returns
    /// [`MeshError::EmptyMask`] when no voxel passes the threshold.
    pub fn voxels_to_hex_mesh(&self, mask: &SegmentationMask) -> Result<VoxelHexMesh, MeshError> {
        let (nx, ny, nz) = mask.dims;
        // Node grid: (nx+1)(ny+1)(nz+1); index = (i*(ny+1)+j)*(nz+1)+k
        let stride_i = (ny + 1) * (nz + 1);
        let stride_j = nz + 1;
        let node_idx = |i: usize, j: usize, k: usize| i * stride_i + j * stride_j + k;

        // First pass: which nodes are referenced (compact output).
        let mut referenced = vec![false; (nx + 1) * (ny + 1) * (nz + 1)];
        let mut elements = Vec::new();
        let mut materials = Vec::new();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if !mask.is_solid(x, y, z) {
                        continue;
                    }
                    let corners = [
                        node_idx(x, y, z),
                        node_idx(x + 1, y, z),
                        node_idx(x + 1, y + 1, z),
                        node_idx(x, y + 1, z),
                        node_idx(x, y, z + 1),
                        node_idx(x + 1, y, z + 1),
                        node_idx(x + 1, y + 1, z + 1),
                        node_idx(x, y + 1, z + 1),
                    ];
                    for &c in &corners {
                        referenced[c] = true;
                    }
                    // Mean HU over the voxel (one voxel = one element here).
                    let mean_hu = mask.index(x, y, z).map(|i| mask.hu[i]).unwrap_or(f64::NAN);
                    let density = HounsfieldMapper::hu_to_density(mean_hu);
                    let region = if density.value() >= self.region_split_density {
                        BoneRegion::Cortical
                    } else {
                        BoneRegion::Trabecular
                    };
                    let modulus = HounsfieldMapper::density_to_youngs_modulus(density, region);
                    elements.push(corners.map(|c| c as u32));
                    materials.push(ElementMaterial {
                        mean_hu,
                        density: density.value(),
                        youngs_modulus: modulus.to_mpa(),
                        region,
                        poissons_ratio: self.poissons_ratio,
                    });
                }
            }
        }
        if elements.is_empty() {
            return Err(MeshError::EmptyMask);
        }

        // Compact referenced nodes into the output arrays.
        let mut remap: HashMap<usize, u32> = HashMap::new();
        let mut nodes = Vec::new();
        for g in 0..referenced.len() {
            if referenced[g] {
                let (i, j, k) = (g / stride_i, (g / stride_j) % (ny + 1), g % (nz + 1));
                remap.insert(g, nodes.len() as u32);
                nodes.push(mask.node_position(i, j, k));
            }
        }
        for el in &mut elements {
            for c in el.iter_mut() {
                *c = remap[&(*c as usize)];
            }
        }

        Ok(VoxelHexMesh {
            nodes,
            elements,
            materials,
        })
    }
}

impl VoxelHexMesh {
    /// Maximum Young's modulus in the mesh (MPa).
    pub fn max_modulus(&self) -> f64 {
        self.materials
            .iter()
            .map(|m| m.youngs_modulus)
            .fold(0.0, f64::max)
    }

    /// Mean Young's modulus (MPa).
    pub fn mean_modulus(&self) -> f64 {
        if self.elements.is_empty() {
            return 0.0;
        }
        self.materials.iter().map(|m| m.youngs_modulus).sum::<f64>() / self.elements.len() as f64
    }

    /// Mesh bounding box in patient coordinates.
    pub fn bounds(&self) -> (Vec3, Vec3) {
        let mut min = self.nodes[0];
        let mut max = self.nodes[0];
        for n in &self.nodes[1..] {
            min = Vec3::new(min.x.min(n.x), min.y.min(n.y), min.z.min(n.z));
            max = Vec3::new(max.x.max(n.x), max.y.max(n.y), max.z.max(n.z));
        }
        (min, max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SegmentationMask;
    use tpt_med_dicom::synthetic;

    #[test]
    fn phantom_meshes_with_expected_topology() {
        let series = synthetic::femur_phantom(12, 12, 6);
        let mask = SegmentationMask::threshold_hu(&series.parse().unwrap(), 200.0);
        let mesher = MedicalMesher::default();
        let mesh = mesher.voxels_to_hex_mesh(&mask).expect("meshes");

        assert_eq!(mesh.elements.len(), mask.solid_count());
        // Shared corner nodes: strictly fewer nodes than 8×elements.
        assert!(mesh.nodes.len() < 8 * mesh.elements.len());

        // All cortical elements in the shell (700 HU → ρ=1.7 ≥ 1.3).
        assert!(mesh
            .materials
            .iter()
            .all(|m| m.region == BoneRegion::Cortical));
        // Cortical modulus at 700 HU: 10500·1.7² = 30345 MPa.
        assert!((mesh.mean_modulus() - 10500.0 * 1.7f64.powi(2)).abs() < 1e-6);

        let (min, max) = mesh.bounds();
        assert!(max.x > min.x);
    }

    #[test]
    fn empty_mask_errors() {
        let series = synthetic::SyntheticCtBuilder::new()
            .cols(4)
            .rows(4)
            .slices(2)
            .hu_fn(|_, _, _| -1000.0) // pure air
            .build("1.2.3");
        let mask = SegmentationMask::threshold_hu(&series.parse().unwrap(), 200.0);
        let err = MedicalMesher::default()
            .voxels_to_hex_mesh(&mask)
            .unwrap_err();
        assert!(matches!(err, MeshError::EmptyMask));
    }

    #[test]
    fn trabecular_only_series_uses_trabecular_law() {
        // 300 HU → ρ = 1.3 g/cm³ exactly at split → cortical; 250 HU →
        // ρ = 1.25 < 1.3 → trabecular.
        let series = synthetic::SyntheticCtBuilder::new()
            .cols(4)
            .rows(4)
            .slices(2)
            .hu_fn(|_, _, _| 250.0)
            .build("1.2.3.4");
        let parsed: Vec<_> = series
            .slices
            .iter()
            .map(|s| tpt_med_dicom::DicomParser::parse_bytes(&s.bytes).unwrap())
            .collect();
        let ct = tpt_med_dicom::DicomSeries::from_slices(parsed).unwrap();
        let mask = SegmentationMask::threshold_hu(&ct, 200.0);
        let mesh = MedicalMesher::default().voxels_to_hex_mesh(&mask).unwrap();
        assert!(mesh
            .materials
            .iter()
            .all(|m| m.region == BoneRegion::Trabecular));
        let e = 6850.0 * 1.25f64.powf(1.49);
        assert!((mesh.mean_modulus() - e).abs() < 1e-6);
    }
}
