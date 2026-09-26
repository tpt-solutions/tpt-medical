//! qform quaternion → rotation matrix, per `nifti1.h`'s documented formula.

use tpt_med_geometry::Mat3;

/// Builds the unit-column rotation matrix from the qform quaternion
/// components `(b, c, d)` (the stored components; `a` is derived).
///
/// `a = sqrt(max(0, 1 - b² - c² - d²))` — clamped rather than left to go
/// negative under a non-unit `(b,c,d)`, per the module docs' documented
/// "trust the header, verify what's cheap to verify" posture (see RFC 0006,
/// Drawbacks).
///
/// The returned matrix's columns are the voxel i/j/k axis directions in RAS,
/// each unit length — voxel spacing is applied separately by the caller
/// (`pixdim[1..4]`, with the k column additionally negated when `qfac` is
/// `-1`), matching how [`tpt_med_geometry::ImageFrame`] keeps direction
/// cosines and spacing apart rather than baking one into the other.
pub fn quaternion_to_rotation(b: f64, c: f64, d: f64) -> Mat3 {
    let a = (1.0 - b * b - c * c - d * d).max(0.0).sqrt();

    // Standard unit-quaternion -> rotation matrix (nifti1.h's own formula,
    // identical to the general quaternion-to-matrix identity).
    Mat3::from_array([
        a * a + b * b - c * c - d * d,
        2.0 * b * c - 2.0 * a * d,
        2.0 * b * d + 2.0 * a * c,
        2.0 * b * c + 2.0 * a * d,
        a * a + c * c - b * b - d * d,
        2.0 * c * d - 2.0 * a * b,
        2.0 * b * d - 2.0 * a * c,
        2.0 * c * d + 2.0 * a * b,
        a * a + d * d - b * b - c * c,
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_quaternion_is_identity_rotation() {
        let r = quaternion_to_rotation(0.0, 0.0, 0.0);
        assert_eq!(r, Mat3::IDENTITY);
    }

    /// Hand-derived: b=0, c=0, d=1 => a=0. Substituting into the formula
    /// above gives diag(-1, -1, 1), a 180-degree rotation about the k-axis
    /// (x -> -x, y -> -y, z -> z). Worked by hand for this RFC's review, not
    /// copied from an external tool -- see rfcs/0006-nifti-ingestion.md.
    #[test]
    fn hand_derived_180_degree_rotation_about_k() {
        let r = quaternion_to_rotation(0.0, 0.0, 1.0);
        assert_eq!(r, Mat3::diagonal([-1.0, -1.0, 1.0]));
    }

    /// Property check independent of any single worked example: for a swept
    /// range of *unit* quaternions, the resulting matrix must be a proper
    /// rotation (orthonormal, determinant +1). This catches an algebra
    /// transcription error that a single hand-checked case could miss.
    #[test]
    fn unit_quaternions_produce_orthonormal_rotations() {
        let cases: [(f64, [f64; 3]); 5] = [
            (0.3, [1.0, 0.0, 0.0]),
            (1.1, [0.0, 1.0, 0.0]),
            (2.0, [0.0, 0.0, 1.0]),
            (0.7, [1.0, 1.0, 1.0]),
            (2.5, [1.0, -1.0, 0.5]),
        ];
        for &(theta, axis) in &cases {
            let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
            let (bx, by, bz) = (axis[0] / norm, axis[1] / norm, axis[2] / norm);
            let s = (theta / 2.0).sin();
            let (b, c, d) = (bx * s, by * s, bz * s);

            let r = quaternion_to_rotation(b, c, d);
            let rt = r.transpose();
            let identity = r.mul_mat(&rt);
            for i in 0..3 {
                for j in 0..3 {
                    let want = if i == j { 1.0 } else { 0.0 };
                    assert!(
                        (identity.at(i, j) - want).abs() < 1e-10,
                        "R*R^T not identity at ({i},{j}): {}",
                        identity.at(i, j)
                    );
                }
            }
            assert!(
                (r.det() - 1.0).abs() < 1e-9,
                "determinant {} != 1 for theta={theta} axis={axis:?}",
                r.det()
            );
        }
    }

    /// 90-degree rotation about the k-axis (b=c=0, d=sin(45deg)): x -> y,
    /// y -> -x, z -> z. A second hand-derived case at a different angle
    /// than the 180-degree one above.
    #[test]
    fn hand_derived_90_degree_rotation_about_k() {
        let s = std::f64::consts::FRAC_1_SQRT_2; // sin(45 deg) = cos(45 deg)
        let r = quaternion_to_rotation(0.0, 0.0, s);
        // a = sqrt(1 - s^2) = s, so a == d here.
        let want = Mat3::from_array([0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]);
        for i in 0..3 {
            for j in 0..3 {
                assert!((r.at(i, j) - want.at(i, j)).abs() < 1e-10);
            }
        }
    }
}
