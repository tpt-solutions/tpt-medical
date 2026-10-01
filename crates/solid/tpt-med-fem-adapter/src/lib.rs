//! 3-D nonlinear hyperelastic FEM assembly and unilateral contact, built on
//! the TPT substrate's `tpt-fem` primitives.
//!
//! This crate is the adapter `rfcs/0002-hyperelastic-tissue.md` and
//! `rfcs/0004-nitinol-superelasticity.md` pointed at and
//! `rfcs/0009-nonlinear-fem-substrate-adapter.md` scoped. The substrate at
//! 0.1.0 ships `Hex8` shape functions, a `CellType::Hex` mesh, hyperelastic
//! *stress functions*, a generic Newton driver and a COO/CSR assembly — but no
//! 3-D hyperelastic assembly, which is what this crate adds:
//!
//! | Piece | Provided by | This crate |
//! |---|---|---|
//! | `Hex8`/`Hex20`/`Hex27` shape functions, gradients, `J^-T` | `tpt-fem-element` | used as-is, via `mesh::Mesh<E>` |
//! | `2x2x2` / `3x3x3` tensor-product rules | `tpt-fem-quadrature` (via `tpt-fem-element`) | used as-is |
//! | `P = dW/dF` per quadrature point | `tpt-med-tissue` | used as-is |
//! | Internal force `int B^T P dV` | — | [`assembly::internal_force`] |
//! | Tangent stiffness `int B^T A B dV` | — | [`assembly::tangent_stiffness`] |
//! | Newton driver, Dirichlet condensation | `tpt-fem-solve` | wired in [`solver::solve_static`] |
//! | Load path (proportional stepping + cutback) | — | [`loadpath::solve_load_path`] |
//! | Global sparse assembly + solve | `tpt-fem-sparse` | wired in |
//! | Unilateral penalty, node pairing | `tpt-fem-contact` | extended to the nonlinear case in [`contact`] |
//! | Friction | — (substrate has none at 0.1.0) | regularized Coulomb layer in [`friction`] |
//!
//! Units follow the workspace convention: lengths in mm, stresses in MPa, so
//! the assembled stiffness is in MPa and a penalty in the same units is a
//! stiffness too (force per unit area per unit length, per node).
//!
//! # Verification
//!
//! Code verification against the closed form the in-house core already
//! verifies against (`mu (lambda - lambda^-2)`, the same nominal uniaxial
//! stress `tpt-med-tissue`'s `substrate-cross-check` compares to
//! `tpt-fem-hyperelastic`), plus a patch-test identity that is exact for a
//! uniform deformation, tangent cross-checks between the two tangent
//! strategies, and mesh refinement. See this crate's README for the table.
//!
//! Uniaxial tension of a cube on its symmetry planes, prescribed on the top
//! face and free on the sides, solved to equilibrium and reduced to a nominal
//! stress.
//!
//! ```
//! use tpt_med_fem_adapter::{hex_box, internal_force, solve_static, AssemblyOptions, SolveOptions};
//! use tpt_med_tissue::{NeoHookeanParams, TissueModel};
//!
//! let l = 10.0;
//! let lam = 1.3;
//! let mesh = hex_box(4, 4, 4, l, l, l).expect("box");
//! // `d1` is the in-house volumetric penalty parameter: a *small* value is a
//! // *stiff* penalty (see this crate's README, "Conventions").
//! let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.49, d1: 0.5 });
//!
//! let mut dirichlet: Vec<(usize, f64)> = Vec::new();
//! for n in mesh.face_nodes(1, false) {
//!     for c in 0..3 {
//!         dirichlet.push((mesh.dof(n, c), 0.0));
//!     }
//! }
//! for n in mesh.face_nodes(2, false) {
//!     dirichlet.push((mesh.dof(n, 2), 0.0));
//! }
//! for n in mesh.face_nodes(1, true) {
//!     dirichlet.push((mesh.dof(n, 1), (lam - 1.0) * l));
//! }
//!
//! let opts = SolveOptions {
//!     assembly: AssemblyOptions::with_quadrature_order(3),
//!     ..SolveOptions::default()
//! };
//! let result = solve_static(
//!     &mesh,
//!     &model,
//!     &vec![0.0; mesh.dof_count()],
//!     &dirichlet,
//!     &opts,
//!     None,
//! )
//! .expect("converges");
//!
//! // Nominal axial stress from the reaction on the prescribed top face.
//! let internal = internal_force(&mesh, &model, &result.displacement, &opts.assembly)
//!     .expect("assembly");
//! let reaction: f64 = mesh
//!     .face_nodes(1, true)
//!     .iter()
//!     .map(|&n| internal[mesh.dof(n, 1)])
//!     .sum();
//! let nominal = reaction / (l * l);
//!
//! // The incompressible Neo-Hookean closed form the in-house core also uses.
//! let expected = 2.0 * 0.49 * (lam - lam.powi(-2));
//! assert!(
//!     (nominal - expected).abs() / expected < 0.03,
//!     "{nominal} vs {expected}"
//! );
//! ```
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod assembly;
pub mod contact;
pub mod friction;
pub mod loadpath;
pub mod mesh;
pub mod mixed;
pub mod solver;

pub use assembly::{
    coo_max_abs_diff, element_deformation_gradient, internal_force, material_tangent,
    tangent_stiffness, tangent_stiffness_numerical, AssemblyOptions, Constitutive, FnModel,
    Tensor4,
};
pub use contact::{
    Constraint as ContactConstraint, ContactCandidate, ContactError, ContactPairing,
};
pub use friction::{friction_terms, FrictionConfig, FrictionError, FrictionState, FrictionTerms};
pub use loadpath::{solve_load_path, LoadPath, LoadPathError, LoadPathOptions, LoadStep};
pub use mesh::{
    hex_box, hex_box_of, tet_box, tet_box_of, ElementFamily, Hex8Mesh, HexFamily, Mesh, MeshError,
    TetFamily,
};
pub use mixed::{solve_mixed_static, MixedOptions, MixedSolveResult};
pub use solver::{
    residual, solve_static, ContactConfig, ContactSummary, SolveError, SolveOptions, SolveResult,
};

#[cfg(test)]
mod tests;
