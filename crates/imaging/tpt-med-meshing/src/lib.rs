//! Threshold segmentation and voxel-to-hexahedral meshing.
//!
//! The pipeline: a [`DicomSeries`](tpt_med_dicom::DicomSeries) is thresholded
//! into a [`SegmentationMask`] (per-voxel bone / not-bone), and
//! [`MedicalMesher::voxels_to_hex_mesh`] converts the mask into a structured
//! hexahedral mesh with one element per solid voxel and HU-derived elastic
//! moduli per element. Laplacian smoothing relaxes the stair-step surface
//! for visualization.
//!
//! This "voxel-to-hex" strategy avoids GPL meshing tools entirely and is
//! the standard approach for CT-based biomechanical models.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod csv;
pub mod mask;
pub mod mesh;
pub mod smooth;

pub use mask::SegmentationMask;
pub use mesh::{ElementMaterial, MedicalMesher, VoxelHexMesh};
pub use smooth::smooth_mesh;

/// Crate-level error type.
#[derive(Debug)]
pub enum MeshError {
    /// I/O failure while writing an export.
    Io(std::io::Error),
    /// The mask has no solid voxels to mesh.
    EmptyMask,
    /// Inconsistent mask/mesh state.
    Inconsistent(String),
}

impl core::fmt::Display for MeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            MeshError::Io(e) => write!(f, "io error: {e}"),
            MeshError::EmptyMask => write!(f, "segmentation mask has no solid voxels"),
            MeshError::Inconsistent(why) => write!(f, "inconsistent mesh state: {why}"),
        }
    }
}

impl std::error::Error for MeshError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            MeshError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for MeshError {
    fn from(e: std::io::Error) -> Self {
        MeshError::Io(e)
    }
}

/// Crate result alias.
pub type Result<T> = core::result::Result<T, MeshError>;
