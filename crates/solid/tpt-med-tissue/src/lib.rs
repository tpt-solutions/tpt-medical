//! Hyperelastic soft-tissue constitutive models.
//!
//! Strain-energy functions over the deformation gradient `F`, with analytic
//! first Piola–Kirchhoff stress and a finite-difference reference used to
//! verify every analytic derivative in CI (ASME V&V 40 code verification).
//!
//! Models: Neo-Hookean, Mooney–Rivlin, Yeoh, compressible-form Ogden
//! (energy), and Holzapfel–Gasser–Ogden (arterial wall, see RFC 0002).
//!
//! Units: stresses in MPa when moduli are given in MPa (the workspace
//! convention).
//!
//! # Examples
//!
//! ```
//! use tpt_med_geometry::{Mat3, Vec3};
//! use tpt_med_tissue::{TissueModel, NeoHookeanParams};
//!
//! let model = TissueModel::NeoHookean(NeoHookeanParams {
//!     c10: 0.49, // ~1.17 kPa shear modulus in MPa-ish scale
//!     d1: 2.0,
//! });
//! // Uniaxial stretch λ = 1.5 with incompressible lateral contraction.
//! let lam = 1.5f64;
//! let f = Mat3::from_rows(
//!     Vec3::new(lam, 0.0, 0.0),
//!     Vec3::new(0.0, 1.0 / lam.sqrt(), 0.0),
//!     Vec3::new(0.0, 0.0, 1.0 / lam.sqrt()),
//! );
//! let p = model.first_piola(&f);
//! let pn = model.first_piola_numerical(&f);
//! // Analytic stress is locked to the finite-difference reference.
//! assert!((p.at(0, 0) - pn.at(0, 0)).abs() < 1e-6);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod hgo;
pub mod models;

pub use hgo::HgoParams;
pub use models::{
    MaterialTangent, MooneyRivlinParams, NeoHookeanParams, OgdenParams, PlaneCondition,
    PlaneSolution, ReducedPlaneModel, SoftTissueMaterial, TissueModel, YeohParams,
};

#[cfg(test)]
mod tests;

#[cfg(all(test, feature = "substrate-cross-check"))]
mod substrate_cross_check;
