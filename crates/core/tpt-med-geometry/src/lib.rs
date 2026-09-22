//! Geometry primitives and anatomical coordinate transforms.
//!
//! Self-contained `std`-only linear algebra for medical simulation:
//! [`Vec3`]/[`Mat3`] fixed-size types, [`Plane`] and [`Aabb`] primitives,
//! analytic eigenvalues of symmetric 3×3 matrices (principal stresses), and
//! the DICOM (LPS) ↔ image research (RAS) coordinate conversions.
//!
//! All angles are radians, lengths are caller-defined (the workspace
//! convention is millimetres via `tpt-med-units`).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aabb;
pub mod coords;
pub mod mat3;
pub mod plane;
pub mod vec3;

pub use aabb::Aabb;
pub use coords::{lps_to_ras, ras_to_lps, CoordinateSystem, ImageFrame};
pub use mat3::Mat3;
pub use plane::Plane;
pub use vec3::Vec3;

/// Tolerance used across geometry predicates.
pub const EPS: f64 = 1.0e-12;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reexports_usable() {
        let v = Vec3::new(1.0, 2.0, 2.0);
        assert!((v.norm() - 3.0).abs() < EPS);
        let m = Mat3::IDENTITY;
        assert_eq!(m * v, v);
    }
}
