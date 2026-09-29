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
        self.voxels_to_hex_mesh_with_overrides(mask, &HashMap::new())
    }

    /// Like [`Self::voxels_to_hex_mesh`], but with **per-voxel modulus
    /// overrides**: for every `(x, y, z)` key present in `overrides`, the
    /// element's Young's modulus (MPa) is taken from the map instead of the
    /// HU correlation — the direct hook for QCT-calibrated or
    /// region-specific material assignment. Unlisted voxels keep the
    /// default law; `mean_hu`/`density` stay informational either way.
    pub fn voxels_to_hex_mesh_with_overrides(
        &self,
        mask: &SegmentationMask,
        overrides: &HashMap<(usize, usize, usize), f64>,
    ) -> Result<VoxelHexMesh, MeshError> {
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
                    let modulus_mpa = match overrides.get(&(x, y, z)) {
                        Some(&e) => e,
                        None => {
                            HounsfieldMapper::density_to_youngs_modulus(density, region).to_mpa()
                        }
                    };
                    elements.push(corners.map(|c| c as u32));
                    materials.push(ElementMaterial {
                        mean_hu,
                        density: density.value(),
                        youngs_modulus: modulus_mpa,
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
    /// Welds nodes closer than `tolerance` (mm) into single vertices and
    /// remaps connectivity, returning the number of nodes removed.
    ///
    /// Grid-corner construction never produces coincident nodes, but two
    /// *disconnected* components can end up within welding distance after
    /// Laplacian smoothing, and imports (e.g. a re-meshed half after an
    /// osteotomy) can carry duplicate positions outright. Welding is
    /// opt-in because it changes the topology contract of the mesh.
    ///
    /// The merge is position-hash based: nodes hash to a uniform grid of
    /// cell size `tolerance`, and any two nodes in the same or adjacent
    /// cell within `tolerance` merge into the lower index. Merged
    /// positions average their members.
    pub fn weld_nodes(&mut self, tolerance: f64) -> usize {
        assert!(tolerance >= 0.0, "tolerance must be non-negative");
        if tolerance == 0.0 || self.nodes.is_empty() {
            return 0;
        }
        let inv = 1.0 / tolerance;
        let cell = |v: f64| (v * inv).floor() as i64;

        struct Cand {
            index: u32,
            x: f64,
            y: f64,
            z: f64,
        }
        let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<Cand>> =
            std::collections::HashMap::new();
        let mut remap: Vec<u32> = (0..self.nodes.len() as u32).collect();
        let mut removed = 0usize;

        for (i, node) in self.nodes.iter().enumerate() {
            let key = (cell(node.x), cell(node.y), cell(node.z));
            let mut target: Option<u32> = None;
            'search: for dx in [-1i64, 0, 1] {
                for dy in [-1, 0, 1] {
                    for dz in [-1, 0, 1] {
                        if let Some(cands) = grid.get(&(key.0 + dx, key.1 + dy, key.2 + dz)) {
                            for c in cands {
                                let d2 = (node.x - c.x).powi(2)
                                    + (node.y - c.y).powi(2)
                                    + (node.z - c.z).powi(2);
                                if d2 <= tolerance * tolerance {
                                    target = Some(c.index);
                                    break 'search;
                                }
                            }
                        }
                    }
                }
            }
            match target {
                Some(t) => {
                    remap[i] = t;
                    removed += 1;
                }
                None => {
                    grid.entry(key).or_default().push(Cand {
                        index: i as u32,
                        x: node.x,
                        y: node.y,
                        z: node.z,
                    });
                }
            }
        }

        if removed == 0 {
            return 0;
        }

        // Path-compress: remap may point at a node that itself was merged.
        for i in 0..remap.len() {
            let mut r = remap[i];
            while remap[r as usize] != r {
                r = remap[r as usize];
            }
            remap[i] = r;
        }

        // Compact: one slot per representative; average merged positions.
        let mut keep_of: HashMap<u32, u32> = HashMap::new();
        let mut members: HashMap<u32, Vec<usize>> = HashMap::new();
        let mut new_nodes: Vec<Vec3> = Vec::new();
        for (i, &r) in remap.iter().enumerate() {
            match keep_of.get(&r) {
                None => {
                    keep_of.insert(r, new_nodes.len() as u32);
                    members.insert(r, vec![i]);
                    new_nodes.push(self.nodes[i]);
                }
                Some(_) => members.get_mut(&r).expect("seeded").push(i),
            }
        }
        for (r, idxs) in &members {
            if idxs.len() > 1 {
                let k = *keep_of.get(r).expect("seeded");
                let inv_n = 1.0 / idxs.len() as f64;
                let mut sum = Vec3::ZERO;
                for &i in idxs {
                    sum += self.nodes[i];
                }
                new_nodes[k as usize] = sum * inv_n;
            }
        }

        for el in &mut self.elements {
            for c in el.iter_mut() {
                *c = remap[*c as usize];
            }
        }
        let removed_total = self.nodes.len() - new_nodes.len();
        self.nodes = new_nodes;
        removed_total
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
    fn modulus_overrides_replace_the_correlation() {
        let series = synthetic::femur_phantom(8, 8, 4);
        let mask = SegmentationMask::threshold_hu(&series.parse().unwrap(), 200.0);
        let mut solid = None;
        'find: for z in 0..mask.dims.2 {
            for y in 0..mask.dims.1 {
                for x in 0..mask.dims.0 {
                    if mask.is_solid(x, y, z) {
                        solid = Some((x, y, z));
                        break 'find;
                    }
                }
            }
        }
        let (x, y, z) = solid.expect("phantom has bone");
        let mut overrides = HashMap::new();
        overrides.insert((x, y, z), 1234.5);
        let mesh = MedicalMesher::default()
            .voxels_to_hex_mesh_with_overrides(&mask, &overrides)
            .unwrap();
        let overridden = mesh
            .materials
            .iter()
            .filter(|m| (m.youngs_modulus - 1234.5).abs() < 1e-9)
            .count();
        assert_eq!(overridden, 1, "exactly one element overridden");
        let hit = mesh
            .materials
            .iter()
            .find(|m| (m.youngs_modulus - 1234.5).abs() < 1e-9)
            .expect("override applied");
        assert_eq!(hit.mean_hu, 700.0, "HU stays informational");
    }

    #[test]
    fn weld_merges_coincident_nodes_across_components() {
        let mat = || ElementMaterial {
            mean_hu: 0.0,
            density: 1.0,
            youngs_modulus: 1.0,
            region: BoneRegion::Cortical,
            poissons_ratio: 0.3,
        };
        let mut mesh = VoxelHexMesh {
            nodes: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 1.0),
                Vec3::new(1.0, 1.0, 1.0),
                Vec3::new(0.0, 1.0, 1.0),
                // A second cube whose bottom face duplicates the first
                // cube's top face within welding tolerance (a stale import).
                Vec3::new(0.0, 0.0, 1.0 + 1e-9),
                Vec3::new(1.0, 0.0, 1.0 + 1e-9),
                Vec3::new(1.0, 1.0, 1.0 + 1e-9),
                Vec3::new(0.0, 1.0, 1.0 + 1e-9),
            ],
            elements: vec![[0, 1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 4, 5, 6, 7]],
            materials: vec![mat(), mat()],
        };
        let removed = mesh.weld_nodes(1e-6);
        assert_eq!(removed, 4, "the 4 duplicated shared-face nodes merge");
        assert_eq!(mesh.nodes.len(), 8);
        assert_eq!(mesh.elements[0][4], mesh.elements[1][4]);
        assert_eq!(mesh.elements[0][7], mesh.elements[1][7]);
        for el in &mesh.elements {
            for c in el {
                assert!((*c as usize) < mesh.nodes.len());
            }
        }
    }

    #[test]
    fn weld_with_zero_tolerance_is_noop() {
        let mut mesh = VoxelHexMesh {
            nodes: vec![Vec3::ZERO, Vec3::ZERO],
            elements: vec![[0, 0, 0, 0, 0, 0, 1, 1]],
            materials: vec![ElementMaterial {
                mean_hu: 0.0,
                density: 1.0,
                youngs_modulus: 1.0,
                region: BoneRegion::Cortical,
                poissons_ratio: 0.3,
            }],
        };
        assert_eq!(mesh.weld_nodes(0.0), 0);
        assert_eq!(mesh.nodes.len(), 2);
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
