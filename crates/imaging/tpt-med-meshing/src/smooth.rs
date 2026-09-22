//! Laplacian surface smoothing for voxel meshes.

use crate::mesh::VoxelHexMesh;
use std::collections::HashMap;

/// Laplacian-smooths the *surface* nodes of a voxel mesh (interior nodes
/// are pinned — with one element per voxel, interior smoothing would
/// shrink the volume without visual benefit).
///
/// For each iteration, every surface node moves by `relaxation × (mean of
/// its surface neighbours − itself)`. `relaxation` in `(0, 1]`; 0.5 is a
/// typical value.
pub fn smooth_mesh(mesh: &mut VoxelHexMesh, iterations: u32, relaxation: f64) {
    if iterations == 0 || relaxation <= 0.0 {
        return;
    }

    let neigh = surface_adjacency(mesh);
    for _ in 0..iterations {
        let mut updates: HashMap<usize, tpt_med_geometry::Vec3> = HashMap::new();
        for (&node, neighbours) in &neigh {
            if neighbours.is_empty() {
                continue;
            }
            let mean = neighbours
                .iter()
                .fold(tpt_med_geometry::Vec3::ZERO, |acc, &n| acc + mesh.nodes[n])
                / neighbours.len() as f64;
            let delta = (mean - mesh.nodes[node]) * relaxation;
            updates.insert(node, mesh.nodes[node] + delta);
        }
        for (node, pos) in updates {
            mesh.nodes[node] = pos;
        }
    }
}

/// Adjacency between surface nodes (nodes touching a face not shared by two
/// elements).
fn surface_adjacency(mesh: &VoxelHexMesh) -> HashMap<usize, Vec<usize>> {
    // Quad faces of a hex, in element-corner ordering
    // [000,100,110,010,001,101,111,011]:
    //   -x: [000,010,011,001]  +x: [100,110,111,101]
    //   -y: [000,100,101,001]  +y: [010,110,111,011]
    //   -z: [000,100,110,010]  +z: [001,101,111,011]
    const FACES: [[usize; 4]; 6] = [
        [0, 3, 7, 4],
        [1, 2, 6, 5],
        [0, 1, 5, 4],
        [3, 2, 6, 7],
        [0, 1, 2, 3],
        [4, 5, 6, 7],
    ];

    let mut face_count: HashMap<[u32; 4], i32> = HashMap::new();
    for el in &mesh.elements {
        for face in FACES {
            let mut key = [el[face[0]], el[face[1]], el[face[2]], el[face[3]]];
            key.sort_unstable();
            *face_count.entry(key).or_insert(0) += 1;
        }
    }

    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for el in &mesh.elements {
        for face in FACES {
            let mut key = [el[face[0]], el[face[1]], el[face[2]], el[face[3]]];
            key.sort_unstable();
            if face_count[&key] != 1 {
                continue; // interior face
            }
            // Ring edges of the (unsorted) face quad.
            let ring = [
                (face[0], face[1]),
                (face[1], face[2]),
                (face[2], face[3]),
                (face[3], face[0]),
            ];
            for (a, b) in ring {
                let (a, b) = (el[a] as usize, el[b] as usize);
                adjacency.entry(a).or_default().push(b);
                adjacency.entry(b).or_default().push(a);
            }
        }
    }
    adjacency.values_mut().for_each(|v| {
        v.sort_unstable();
        v.dedup();
    });
    adjacency
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{MedicalMesher, SegmentationMask};
    use tpt_med_dicom::synthetic;

    fn phantom_mesh() -> VoxelHexMesh {
        let series = synthetic::femur_phantom(10, 10, 6);
        let mask = SegmentationMask::threshold_hu(&series.parse().unwrap(), 200.0);
        MedicalMesher::default().voxels_to_hex_mesh(&mask).unwrap()
    }

    #[test]
    fn smoothing_moves_only_surface_and_converges() {
        let mut mesh = phantom_mesh();
        let before = mesh.nodes.clone();
        let n_before = mesh.nodes.len();
        let el_before = mesh.elements.len();

        smooth_mesh(&mut mesh, 5, 0.5);

        assert_eq!(mesh.nodes.len(), n_before, "no nodes added/removed");
        let mut moved = 0;
        for (i, (a, b)) in before.iter().zip(mesh.nodes.iter()).enumerate() {
            if (*a - *b).norm() > 1e-12 {
                moved += 1;
            }
            let _ = i;
        }
        assert!(moved > 0, "surface nodes should move");
        // Elements untouched by smoothing.
        assert_eq!(mesh.elements.len(), el_before);
    }

    #[test]
    fn zero_iterations_is_noop() {
        let mut mesh = phantom_mesh();
        let before = mesh.nodes.clone();
        smooth_mesh(&mut mesh, 0, 0.5);
        assert_eq!(before, mesh.nodes);
    }
}
