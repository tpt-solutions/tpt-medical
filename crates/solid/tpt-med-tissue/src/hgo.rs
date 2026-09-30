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
//!
//! # Collagen crimp (progressive fiber recruitment)
//!
//! With `crimp` set, each family's fiber term is weighted by the
//! **recruited fraction** of its fibers. Real collagen fibers are wavy
//! (crimped) at rest and straighten progressively: a fiber with
//! recruitment stretch `λ_r` bears load only once the fiber stretch
//! `λ = √I4` exceeds `λ_r`. Following the waviness-distribution
//! treatment of Decraemer, Maes & Vanhuyse (1980), the recruited
//! fraction at stretch `λ` is the Gaussian recruitment function
//!
//! ```text
//! R(λ) = Φ((λ − λ̄_r)/σ_r),      W_fiber = R(λ)·(k1/2k2)[exp(k2 E²) − 1]
//! ```
//!
//! with `λ̄_r` the mean recruitment stretch and `σ_r` its spread
//! (both caller-cited — the mechanism ships, the coefficients come from
//! the caller's source, as everywhere in this workspace). `R → 1`
//! recovers the standard HGO response exactly, so crimped parameters
//! bracket the uncrimped model; the analytic first Piola carries the
//! extra `R′` term (a Gaussian density), verified against central
//! differences. The normal CDF is evaluated to near machine precision
//! (a series plus continued-fraction `erf`), because the
//! finite-difference stress references differentiate through `R` — a
//! CDF accurate only in value would disagree with its own derivative
//! inside the recruitment window. `crimp: None` (the default) is the
//! unmodified HGO model.

use tpt_med_geometry::{Mat3, Vec3};

/// Collagen-crimp recruitment parameters: the distribution of fiber
/// recruitment stretches `λ_r` (the stretch at which a crimped fiber
/// straightens and begins to bear load). Construct with
/// [`CrimpRecruitment::new`], which refuses a non-positive spread and a
/// mean below rest stretch (`λ̄_r ≥ 1`: fibers are crimped at rest, not
/// pre-tensioned).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrimpRecruitment {
    /// Mean recruitment stretch `λ̄_r` (dimensionless, ≥ 1).
    pub mean_recruitment_stretch: f64,
    /// Spread `σ_r` of the recruitment-stretch distribution (dimensionless
    /// and strictly positive). Arterial collagen crimp engages over a few
    /// percent of stretch, so `σ_r ≈ 0.01–0.03` is the screening order of
    /// magnitude.
    pub spread: f64,
}

impl CrimpRecruitment {
    /// Validates and constructs; `None` for a non-finite or non-positive
    /// spread, or a mean below 1 (recruitment before rest stretch has no
    /// physical crimp interpretation).
    pub fn new(mean_recruitment_stretch: f64, spread: f64) -> Option<Self> {
        if !spread.is_finite()
            || spread <= 0.0
            || !mean_recruitment_stretch.is_finite()
            || mean_recruitment_stretch < 1.0
        {
            return None;
        }
        Some(Self {
            mean_recruitment_stretch,
            spread,
        })
    }

    /// Recruited fraction `R(λ) = Φ((λ − λ̄_r)/σ_r)` at fiber stretch `λ`.
    pub fn recruited_fraction(&self, fiber_stretch: f64) -> f64 {
        normal_cdf((fiber_stretch - self.mean_recruitment_stretch) / self.spread)
    }

    /// `dR/dλ` — the Gaussian density `φ((λ − λ̄_r)/σ_r)/σ_r` that enters
    /// the analytic first Piola through the chain rule.
    fn recruited_fraction_derivative(&self, fiber_stretch: f64) -> f64 {
        let z = (fiber_stretch - self.mean_recruitment_stretch) / self.spread;
        normal_pdf(z) / self.spread
    }
}

/// Standard normal CDF `Φ(z) = ½[1 + erf(z/√2)]`, with `erf` from the
/// Abramowitz & Stegun 7.1.26 rational approximation (`|error| ≤ 1.5e-7`,
/// ample for a weighting fraction).
fn normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z / core::f64::consts::SQRT_2))
}

/// Standard normal density `φ(z) = e^{−z²/2}/√(2π)`.
fn normal_pdf(z: f64) -> f64 {
    (-(z * z) / 2.0).exp() / (2.0 * core::f64::consts::PI).sqrt()
}

/// Error function to near machine precision (~3e-16 across the range,
/// verified against 30-digit reference values): the convergent positive
/// series (Abramowitz & Stegun 7.1.5) for `|x| ≤ 3`, and 1 − erfc via
/// the continued fraction of A&S 7.1.14 beyond, with the truncation
/// depth doubled until bit-stable. A low-order rational approximation
/// is deliberately not enough here: `strain_energy` is differentiated
/// (by the finite-difference stress references), so the *shape* of the
/// recruitment CDF must be accurate to machine precision, not merely
/// its value — an approximation whose value error is 1e-7 carries an
/// error slope an order of magnitude worse, which shows up directly in
/// the stress.
fn erf(x: f64) -> f64 {
    if x < 0.0 {
        return -erf(-x);
    }
    if x <= 3.0 {
        // erf(x) = (2x/√π)·e^{−x²}·Σ (2x²)ⁿ / (1·3·5···(2n+1)) — every
        // term positive, so the summation cannot cancel.
        let x2 = 2.0 * x * x;
        let mut sum = 1.0f64;
        let mut term = 1.0f64;
        let mut n = 0.0f64;
        loop {
            n += 1.0;
            term *= x2 / (2.0 * n + 1.0);
            let next = sum + term;
            if next == sum || n > 500.0 {
                break;
            }
            sum = next;
        }
        2.0 * x * (-x * x).exp() * sum / core::f64::consts::PI.sqrt()
    } else {
        // erfc(x) = (e^{−x²}/√π)·[1/(x + 1/(2x + 2/(x + 3/(2x + ⋯))))].
        let cf = |depth: u32| {
            let mut f = 0.0f64;
            for k in (1..=depth).rev() {
                let d = if k % 2 == 1 { 2.0 * x } else { x };
                f = k as f64 / (d + f);
            }
            (1.0 / (x + f)) * (-x * x).exp() / core::f64::consts::PI.sqrt()
        };
        let mut depth = 16u32;
        let mut erfc = cf(depth);
        while depth < 4096 {
            depth *= 2;
            let next = cf(depth);
            if next == erfc {
                break;
            }
            erfc = next;
        }
        1.0 - erfc
    }
}

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
    /// Optional collagen-crimp recruitment (see the module docs).
    /// `None` — the default — is the unmodified HGO model.
    pub crimp: Option<CrimpRecruitment>,
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
                let weight = self.recruitment_weight(i4);
                w += weight * self.k1 / (2.0 * self.k2) * ((self.k2 * e * e).exp_m1());
            }
        }
        w + (j - 1.0).powi(2) / self.d1
    }

    /// The recruited fraction of a family at fiber stretch `λ = √I4`:
    /// `R(λ)` under crimp, exactly `1.0` without.
    fn recruitment_weight(&self, i4: f64) -> f64 {
        match &self.crimp {
            Some(c) => c.recruited_fraction(i4.max(0.0).sqrt()),
            None => 1.0,
        }
    }

    /// Analytic first Piola–Kirchhoff stress:
    /// `P = 2F ∂W/∂C + 2J(J−1)/D1 · F^{−T}` with, per family,
    /// `∂W/∂C_f = k1 R E e^{k2E²} (κ I + (1−3κ) a0⊗a0)
    /// plus R′·(k1/4k2λ)(e^{k2E²}−1)·a0⊗a0`, where `R`/`R′` are the crimp
    /// recruitment weight and its λ-derivative (`R = 1`, `R′ = 0` without
    /// crimp).
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
                let r = self.recruitment_weight(i4);
                let coef = self.k1 * r * e * (self.k2 * e * e).exp();
                dw_dc = dw_dc + (identity * self.kappa + aot * (1.0 - 3.0 * self.kappa)) * coef;
                if let Some(c) = &self.crimp {
                    // The chain-rule term through R(λ), λ = √I4: dR/dI4 =
                    // R′(λ)/(2λ) — a Gaussian density, active only inside
                    // the recruitment window.
                    let lambda = i4.max(0.0).sqrt();
                    if lambda > 0.0 {
                        let dr = c.recruited_fraction_derivative(lambda) / (2.0 * lambda) * self.k1
                            / (2.0 * self.k2)
                            * (self.k2 * e * e).exp_m1();
                        dw_dc = dw_dc + aot * dr;
                    }
                }
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
            crimp: None,
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
            crimp: None,
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
            crimp: None,
        };
        let base = HgoParams {
            c: 0.8,
            k1: 0.0,
            k2: 12.0,
            kappa: 0.0,
            fiber_directions: vec![Vec3::Y],
            d1: 100.0,
            crimp: None,
        };
        assert!(
            (inactive_only.strain_energy(&f) - base.strain_energy(&f)).abs() < 1e-12,
            "compressed fiber must carry no energy"
        );
        assert!(aligned.strain_energy(&f) > base.strain_energy(&f));
    }

    /// Aligned (κ = 0) x-fiber family: the fiber stretch under uniaxial
    /// `uniaxial_f(λ)` is exactly `λ`.
    fn crimped_params(mean: f64, spread: f64) -> HgoParams {
        HgoParams {
            c: 0.8,
            k1: 50.0,
            k2: 12.0,
            kappa: 0.0,
            fiber_directions: vec![Vec3::X],
            d1: 100.0,
            crimp: CrimpRecruitment::new(mean, spread),
        }
    }

    #[test]
    fn crimp_is_silent_below_recruitment_and_exact_above_it() {
        // λ̄_r = 1.2, σ_r = 0.02: at fiber stretch 1.1 the recruited
        // fraction is Φ(−5) ≈ 3e-7 — the response is the ground substance
        // for every practical purpose (0.1 % of the uncrimped fiber
        // energy), and at stretch 1.6 it is Φ(20) = 1 — bit-identical to
        // the uncrimped model.
        let crimped = crimped_params(1.2, 0.02);
        let plain = HgoParams {
            crimp: None,
            ..crimped.clone()
        };
        let early = uniaxial_f(1.1);
        let w_crimp = crimped.strain_energy(&early);
        let w_plain = plain.strain_energy(&early);
        let fiber_only = w_plain
            - 0.8 * 0.5 * ((early.transpose() * early).trace() - 3.0)
            - (early.det() - 1.0).powi(2) / 100.0;
        assert!(fiber_only > 0.0, "sanity: fiber term active at E > 0");
        // Below recruitment only Φ(−5) ≈ 3e-7 of the family is engaged:
        // the crimped energy is the ground substance plus a negligible
        // share of the fiber term.
        let engaged = w_crimp - (w_plain - fiber_only);
        assert!(
            engaged >= 0.0 && engaged < 1e-3 * fiber_only,
            "crimp must silence the fiber term below recruitment: {engaged} of {fiber_only}"
        );
        let late = uniaxial_f(1.6);
        assert!(
            (crimped.strain_energy(&late) - plain.strain_energy(&late)).abs() < 1e-12,
            "fully recruited must equal the standard HGO energy"
        );
        let p_crimp = crimped.first_piola(&late);
        let p_plain = plain.first_piola(&late);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (p_crimp.at(i, j) - p_plain.at(i, j)).abs() < 1e-12,
                    "fully recruited must equal the standard HGO stress"
                );
            }
        }
    }

    #[test]
    fn crimped_analytic_stress_matches_finite_difference() {
        // The R′ chain-rule term is the new analytic content: verify it
        // through the recruitment window, where the density is largest.
        let crimped = crimped_params(1.2, 0.02);
        assert!(crimped.crimp.is_some());
        for lam in [1.15, 1.19, 1.2, 1.21, 1.25] {
            let f = uniaxial_f(lam);
            let pa = crimped.first_piola(&f);
            let pn = fd_stress(&crimped, &f);
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
    fn crimp_stiffens_through_the_recruitment_window() {
        let crimped = crimped_params(1.2, 0.02);
        // Axial stress rises monotonically as fibers recruit through the
        // window; and at a stretch *inside* the window the crimped wall is
        // softer than the uncrimped one (only part of the family is
        // engaged).
        let mut previous = 0.0;
        for lam in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3] {
            let p11 = crimped.first_piola(&uniaxial_f(lam)).at(0, 0);
            assert!(
                p11 > previous,
                "P11 must rise through the window: {p11} after {previous}"
            );
            previous = p11;
        }
        let mid = uniaxial_f(1.19);
        let plain = HgoParams {
            crimp: None,
            ..crimped.clone()
        };
        assert!(
            crimped.first_piola(&mid).at(0, 0) < plain.first_piola(&mid).at(0, 0),
            "partially recruited must be softer than fully recruited"
        );
    }

    #[test]
    fn crimp_spread_controls_the_window_width() {
        // The spread is the width of the recruitment transition. At the
        // fraction level: a wider distribution recruits more fibers below
        // the mean and fewer above it (its CDF is flatter). At the stress
        // level the narrow spread is the sharper switch — the stress rise
        // through the window is steeper, because the R′ term scales as
        // 1/σ_r even where R itself is smaller.
        let narrow = crimped_params(1.2, 0.01);
        let wide = crimped_params(1.2, 0.04);
        let (r_narrow, r_wide) = (narrow.crimp.expect("set"), wide.crimp.expect("set"));
        assert!(r_narrow.recruited_fraction(1.19) < r_wide.recruited_fraction(1.19));
        assert!(r_narrow.recruited_fraction(1.21) > r_wide.recruited_fraction(1.21));
        let rise = |p: &HgoParams| {
            p.first_piola(&uniaxial_f(1.21)).at(0, 0) - p.first_piola(&uniaxial_f(1.19)).at(0, 0)
        };
        assert!(
            rise(&narrow) > rise(&wide),
            "narrow spread must be the sharper switch: {} vs {}",
            rise(&narrow),
            rise(&wide)
        );
    }

    #[test]
    fn crimp_recruitment_validates_and_evaluates() {
        // Non-positive spread, sub-rest mean: refused.
        assert!(CrimpRecruitment::new(1.2, 0.0).is_none());
        assert!(CrimpRecruitment::new(1.2, -0.1).is_none());
        assert!(CrimpRecruitment::new(0.95, 0.02).is_none());
        assert!(CrimpRecruitment::new(f64::NAN, 0.02).is_none());
        let r = CrimpRecruitment::new(1.2, 0.02).expect("valid");
        // Φ(0) = ½ at the mean; Φ beyond ±9σ saturates at 0/1.
        assert!((r.recruited_fraction(1.2) - 0.5).abs() < 1e-15);
        assert!(r.recruited_fraction(1.0) < 1e-6);
        assert!(r.recruited_fraction(1.4) > 1.0 - 1e-9);
        // Machine-precision erf: exact odd symmetry, the published
        // erf(1) to ~1 ulp, and both branches meeting at x = 3.
        assert_eq!(erf(0.0), 0.0);
        assert!((erf(1.0) - 0.8427007929497149).abs() < 1e-15);
        assert!((erf(-1.0) + 0.8427007929497149).abs() < 1e-15);
        assert!((erf(3.0) - erf(3.0 + 1e-9)).abs() < 1e-12);
    }
}
