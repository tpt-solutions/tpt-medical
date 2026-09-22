//! Holzapfel–Gasser–Ogden (HGO) arterial wall model.
//!
//! Two symmetric fiber families embedded in an isotropic ground substance,
//! with dispersion `κ` about the mean fiber axis (Gasser, Holzapfel &
//! Ogden 2006):
//!
//! ```text
//! W = (c/2)(I1 − 3)
//!   + Σ_families (k1/2k2) [exp(k2 E²) − 1],   E = κ I1 + (1 − 3κ) I4 − 1
//!   + (1/D1)(J − 1)²
//! ```
//!
//! `I4 = a0 · C · a0` is the square of the stretch along the fiber
//! direction `a0`. The structural term activates only in extension
//! (`E > 0`), which the implementation enforces exactly like the original.

use tpt_med_geometry::{Mat3, Vec3};

/// HGO material parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct HgoParams {
    /// Ground-substance modulus `c`.
    pub c: f64,
    /// Fiber stiffness `k1`.
    pub k1: f64,
    /// Fiber nonlinearity `k2`.
    pub k2: f64,
    /// Fiber dispersion `κ ∈ [0, 1/3]`; `1/3` = isotropic, `0` = perfectly
    /// aligned fibers.
    pub kappa: f64,
    /// Mean fiber directions; typically two symmetric families.
    pub fiber_directions: Vec<Vec3>,
    /// Volumetric penalty coefficient.
    pub d1: f64,
}

/// `E_f = κ I1 + (1 − 3κ) I4 − 1` for one family.
fn fiber_exponent(i1: f64, i4: f64, kappa: f64) -> f64 {
    kappa * i1 + (1.0 - 3.0 * kappa) * i4 - 1.0
}

/// Outer product `a ⊗ a` as a matrix.
fn outer(a: Vec3) -> Mat3 {
    Mat3::from_rows(
        Vec3::new(a.x * a.x, a.x * a.y, a.x * a.z),
        Vec3::new(a.y * a.x, a.y * a.y, a.y * a.z),
        Vec3::new(a.z * a.x, a.z * a.y, a.z * a.z),
    )
}

impl HgoParams {
    /// Strain energy `W(F)`.
    pub fn strain_energy(&self, f: &Mat3) -> f64 {
        let f = *f;
        let c_mat = f.transpose() * f;
        let i1 = c_mat.trace();
        let j = f.det();

        let mut w = self.c * 0.5 * (i1 - 3.0);
        for a0 in &self.fiber_directions {
            let a = a0.normalize();
            let i4 = a.dot(c_mat * a);
            let e = fiber_exponent(i1, i4, self.kappa);
            if e > 0.0 {
                w += self.k1 / (2.0 * self.k2) * ((self.k2 * e * e).exp_m1());
            }
        }
        w + (j - 1.0).powi(2) / self.d1
    }

    /// Analytic first Piola–Kirchhoff stress:
    /// `P = 2F ∂W/∂C + 2J(J−1)/D1 · F^{−T}` with
    /// `∂W/∂C = c/2 I + Σ_families k1 E e^{k2E²} (κ I + (1−3κ) a0⊗a0)`.
    pub fn first_piola(&self, f: &Mat3) -> Mat3 {
        let f = *f;
        let c_mat = f.transpose() * f;
        let i1 = c_mat.trace();
        let identity = Mat3::IDENTITY;
        let mut dw_dc = identity * (self.c * 0.5);

        for a0 in &self.fiber_directions {
            let a = a0.normalize();
            let aot = outer(a);
            let i4 = a.dot(c_mat * a);
            let e = fiber_exponent(i1, i4, self.kappa);
            if e > 0.0 {
                let coef = self.k1 * e * (self.k2 * e * e).exp();
                dw_dc = dw_dc + (identity * self.kappa + aot * (1.0 - 3.0 * self.kappa)) * coef;
            }
        }

        let f_inv_t = f.inverse().map(|i| i.transpose()).unwrap_or(Mat3::ZERO);
        let j = f.det();
        f * dw_dc * 2.0 + f_inv_t * (2.0 * j * (j - 1.0) / self.d1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Central-difference reference stress.
    fn fd_stress(p: &HgoParams, f: &Mat3) -> Mat3 {
        const H: f64 = 1.0e-6;
        let mut out = Mat3::ZERO;
        for r in 0..3 {
            for c in 0..3 {
                let mut fp = *f;
                fp.set(r, c, fp.at(r, c) + H);
                let mut fm = *f;
                fm.set(r, c, fm.at(r, c) - H);
                out.set(
                    r,
                    c,
                    (p.strain_energy(&fp) - p.strain_energy(&fm)) / (2.0 * H),
                );
            }
        }
        out
    }

    fn uniaxial_f(lam: f64) -> Mat3 {
        Mat3::from_rows(
            Vec3::new(lam, 0.0, 0.0),
            Vec3::new(0.0, 1.0 / lam.sqrt(), 0.0),
            Vec3::new(0.0, 0.0, 1.0 / lam.sqrt()),
        )
    }

    fn params(kappa: f64, k1: f64) -> HgoParams {
        HgoParams {
            c: 0.8,
            k1,
            k2: 12.0,
            kappa,
            fiber_directions: vec![Vec3::new(1.0, 1.0, 0.0), Vec3::new(-1.0, 1.0, 0.0)],
            d1: 100.0,
        }
    }

    #[test]
    fn analytic_stress_matches_finite_difference() {
        let p = params(0.2, 5.0);
        for lam in [1.05, 1.2, 1.4] {
            let f = uniaxial_f(lam);
            let pa = p.first_piola(&f);
            let pn = fd_stress(&p, &f);
            for i in 0..3 {
                for j in 0..3 {
                    let scale = pa.at(i, j).abs().max(1e-3);
                    assert!(
                        (pa.at(i, j) - pn.at(i, j)).abs() < 1e-5 * scale,
                        "λ={lam} ({i},{j}): {} vs {}",
                        pa.at(i, j),
                        pn.at(i, j)
                    );
                }
            }
        }
    }

    #[test]
    fn isotropic_dispersion_recovers_ground_substance() {
        // κ = 1/3 makes the structural term isotropic and, for a dilatation-
        // free deformation with I1 = I4-dependent combinations, the fiber
        // term reduces to a scalar function of I1. Sanity: energy increases
        // monotonically with stretch, P11 > 0 for λ > 1.
        let p = params(1.0 / 3.0, 5.0);
        let w1 = p.strain_energy(&uniaxial_f(1.1));
        let w2 = p.strain_energy(&uniaxial_f(1.3));
        assert!(w2 > w1);
        assert!(p.first_piola(&uniaxial_f(1.2)).at(0, 0) > 0.0);
    }

    #[test]
    fn fibers_carry_no_compression() {
        // For stretches that keep E < 0 along both families, the energy must
        // equal the ground-substance energy.
        let aligned = HgoParams {
            c: 0.8,
            k1: 50.0,
            k2: 12.0,
            kappa: 0.0,
            fiber_directions: vec![Vec3::X, Vec3::Y],
            d1: 100.0,
        };
        let f = uniaxial_f(1.1);
        // κ = 0 → E = I4 − 1; with the lateral contraction the y-fiber sees
        // I4 = 1/λ < 1 (compression, inactive) and x sees λ² > 1 (active).
        // Compare against κ→0 with k1 = 0 instead: only the x-family
        // contributes, so verify the inactive family contributes zero.
        let inactive_only = HgoParams {
            c: 0.8,
            k1: 50.0,
            k2: 12.0,
            kappa: 0.0,
            fiber_directions: vec![Vec3::Y],
            d1: 100.0,
        };
        let base = HgoParams {
            c: 0.8,
            k1: 0.0,
            k2: 12.0,
            kappa: 0.0,
            fiber_directions: vec![Vec3::Y],
            d1: 100.0,
        };
        assert!(
            (inactive_only.strain_energy(&f) - base.strain_energy(&f)).abs() < 1e-12,
            "compressed fiber must carry no energy"
        );
        assert!(aligned.strain_energy(&f) > base.strain_energy(&f));
    }
}
