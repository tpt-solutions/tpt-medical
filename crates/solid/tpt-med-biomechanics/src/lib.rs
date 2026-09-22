//! Voxel hexahedral FEM solver core for patient-specific biomechanics.
//!
//! Displacement-based small-strain finite elements:
//!
//! - 8-node trilinear hexahedra (Q1) with 2×2×2 Gauss quadrature — the
//!   natural match for voxel meshes produced by `tpt-med-meshing`.
//! - Isotropic linear-elastic materials with per-element properties from
//!   the HU → density → modulus pipeline.
//! - Global assembly into CSR, Jacobi-preconditioned conjugate gradient.
//! - Post-processing: element strains/stresses, von Mises, principal
//!   stresses (analytic symmetric eigenvalues from `tpt-med-geometry`).
//!
//! Nonlinear (hyperelastic, large deformation, contact) capability is the
//! documented upgrade path via `tpt-fem` (see RFC 0002); this crate
//! deliberately keeps the zero-dependency, WASM-friendly linear core that
//! covers screening-level patient-specific stress analysis.
//!
//! # Examples
//!
//! ```no_run
//! use tpt_med_biomechanics::{BiomechanicsModel, BoundaryConditions};
//!
//! # let mesh = tpt_med_meshing::VoxelHexMesh {
//! #     nodes: vec![],
//! #     elements: vec![],
//! #     materials: vec![],
//! # };
//! let model = BiomechanicsModel::from_voxel_mesh(&mesh);
//! let bc = BoundaryConditions::default();
//! // Pin nodes via `bc.fix_nodes(..)`, add loads via `bc.add_force(..)`,
//! // then:
//! # let result = model.solve(&bc, 1e-8, 1).unwrap_or_else(|_| unreachable!());
//! println!("max von Mises: {} MPa", result.max_von_mises());
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod hex;
pub mod results;
pub mod solver;
pub mod sparse;

pub use hex::trilinear_hex_stiffness;
pub use results::{ElementStress, StressResult};
pub use solver::{BiomechanicsModel, BoundaryConditions, SolverError};

#[cfg(test)]
mod tests;
