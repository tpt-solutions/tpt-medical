//! CSV export/import for voxel hex meshes.
//!
//! Format (stable, versioned in the header):
//!
//! ```text
//! # tpt-medical voxel hex mesh v1
//! # units: mm, MPa
//! nodes,<count>
//! node,<id>,<x>,<y>,<z>
//! elements,<count>
//! hex,<id>,<n0>,...,<n7>,<youngs_modulus_mpa>,<density_gcm3>
//! ```

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;

use crate::mesh::{ElementMaterial, VoxelHexMesh};
use crate::{MeshError, Result};
use tpt_med_dicom::BoneRegion;
use tpt_med_geometry::Vec3;

impl VoxelHexMesh {
    /// Writes the mesh to a CSV file.
    pub fn write_csv(&self, path: &Path) -> Result<()> {
        let file = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(file);
        self.write_csv_to(&mut w)?;
        w.flush()?;
        Ok(())
    }

    /// Writes the mesh CSV into any writer.
    pub fn write_csv_to(&self, w: &mut impl std::io::Write) -> Result<()> {
        writeln!(w, "# tpt-medical voxel hex mesh v1")?;
        writeln!(w, "# units: mm, MPa")?;
        writeln!(w, "nodes,{}", self.nodes.len())?;
        for (i, n) in self.nodes.iter().enumerate() {
            writeln!(w, "node,{i},{:.6},{:.6},{:.6}", n.x, n.y, n.z)?;
        }
        writeln!(w, "elements,{}", self.elements.len())?;
        for (i, el) in self.elements.iter().enumerate() {
            let m = &self.materials[i];
            let mut line = String::with_capacity(96);
            let _ = write!(line, "hex,{i}");
            for &c in el {
                let _ = write!(line, ",{c}");
            }
            let _ = write!(line, ",{:.3},{:.6}", m.youngs_modulus, m.density);
            writeln!(w, "{line}")?;
        }
        Ok(())
    }

    /// Parses back a mesh written by [`Self::write_csv`]. Round-trip
    /// fidelity is exact to the printed precision (6 decimals for
    /// coordinates, 3 for moduli, 6 for density). Mean HU and Poisson's
    /// ratio are not exported, so parsed materials carry `NAN` mean HU and
    /// the workspace-default Poisson ratio; bone region is re-inferred from
    /// the density with the default mesher split.
    pub fn parse_csv(text: &str) -> Result<Self> {
        let mut nodes = Vec::new();
        let mut elements = Vec::new();
        let mut materials = Vec::new();
        let mut phase = 0u8; // 0 = nodes, 1 = elements
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split(',').collect();
            match parts[0] {
                "nodes" => {
                    let n: usize = parts
                        .get(1)
                        .ok_or_else(|| bad(line))?
                        .parse()
                        .map_err(|e| MeshError::Inconsistent(format!("nodes count: {e}")))?;
                    nodes.reserve(n);
                    phase = 0;
                }
                "elements" => {
                    let n: usize = parts
                        .get(1)
                        .ok_or_else(|| bad(line))?
                        .parse()
                        .map_err(|e| MeshError::Inconsistent(format!("elements count: {e}")))?;
                    elements.reserve(n);
                    materials.reserve(n);
                    phase = 1;
                }
                "node" => {
                    if phase != 0 || parts.len() < 5 {
                        return Err(bad(line));
                    }
                    let x: f64 = parts[2].parse().map_err(|_| bad(line))?;
                    let y: f64 = parts[3].parse().map_err(|_| bad(line))?;
                    let z: f64 = parts[4].parse().map_err(|_| bad(line))?;
                    nodes.push(Vec3::new(x, y, z));
                }
                "hex" => {
                    if phase != 1 || parts.len() < 12 {
                        return Err(bad(line));
                    }
                    let mut el = [0u32; 8];
                    for (i, p) in parts[2..10].iter().enumerate() {
                        el[i] = p.parse().map_err(|_| bad(line))?;
                    }
                    let modulus: f64 = parts[10].parse().map_err(|_| bad(line))?;
                    let density: f64 = parts[11].parse().map_err(|_| bad(line))?;
                    let region = if density >= 1.3 {
                        BoneRegion::Cortical
                    } else {
                        BoneRegion::Trabecular
                    };
                    elements.push(el);
                    materials.push(ElementMaterial {
                        mean_hu: f64::NAN,
                        density,
                        youngs_modulus: modulus,
                        region,
                        poissons_ratio: 0.30,
                    });
                }
                other => {
                    return Err(MeshError::Inconsistent(format!("unknown record: {other}")));
                }
            }
        }
        if nodes.is_empty() || elements.is_empty() {
            return Err(MeshError::Inconsistent(
                "missing node or element records".into(),
            ));
        }
        Ok(Self {
            nodes,
            elements,
            materials,
        })
    }
}

fn bad(line: &str) -> MeshError {
    MeshError::Inconsistent(format!("malformed record: {line}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{MedicalMesher, SegmentationMask};
    use tpt_med_dicom::synthetic;

    #[test]
    fn csv_roundtrip() {
        let series = synthetic::femur_phantom(8, 8, 4);
        let mask = SegmentationMask::threshold_hu(&series.parse().unwrap(), 200.0);
        let mesh = MedicalMesher::default().voxels_to_hex_mesh(&mask).unwrap();

        let mut buf = Vec::new();
        mesh.write_csv_to(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        let parsed = VoxelHexMesh::parse_csv(&text).unwrap();
        assert_eq!(parsed.nodes.len(), mesh.nodes.len());
        assert_eq!(parsed.elements.len(), mesh.elements.len());
        // Node coordinates round-trip to 6 decimals.
        for (a, b) in mesh.nodes.iter().zip(&parsed.nodes) {
            assert!((*a - *b).norm() < 1e-5, "{a:?} vs {b:?}");
        }
        // Moduli round-trip to 3 decimals.
        for (a, b) in mesh.materials.iter().zip(&parsed.materials) {
            assert!((a.youngs_modulus - b.youngs_modulus).abs() < 1e-2);
            assert!((a.density - b.density).abs() < 1e-5);
        }
    }
}
