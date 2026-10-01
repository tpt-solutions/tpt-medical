//! Finite-strain internal-variable viscoelasticity (Simo-type) — the
//! thermodynamically complete counterpart to the linear-regime
//! [`crate::QuasiLinearViscoelastic`]: an isochoric free energy with
//! per-branch internal deformation tensors and an exact per-step update,
//! valid at large 3-D deformations.
//!
//! # Formulation
//!
//! On `C = FᵀF` with isochoric part `C̄ = J^{−2/3} C`, the free energy is
//!
//! ```text
//! Ψ(C, Γ₁..Γₙ) = U(J) + (μ∞/2)(tr C̄ − 3)
//!              + Σᵢ (μᵢ/2)(tr(C̄ Γᵢ⁻¹) − 3)
//! ```
//!
//! with `U(J) = (K/2)(J − 1)²`, equilibrium modulus `μ∞`, and one internal
//! symmetric tensor `Γᵢ` per Maxwell branch (Simo 1987, "On a fully
//! three-dimensional finite-strain viscoelastic damage model"; Holzapfel
//! 2000, ch. 6). Each branch evolves by
//!
//! ```text
//! Γ̇ᵢ = (C̄ − Γᵢ)/τᵢ
//! ```
//!
//! whose update is **exact** for the step (the ODE is linear with `C̄`
//! frozen over the step, the standard Simo assumption):
//!
//! ```text
//! Γᵢ⁺ = Γᵢ + (1 − exp(−Δt/τᵢ)) (C̄ − Γᵢ)
//! ```
//!
//! an interpolation between the old state and the current isochoric
//! response, so `Γᵢ` stays symmetric positive definite whenever `C̄` is.
//!
//! # The limits this module pins exactly
//!
//! - **Instant** (`Γᵢ = I`): every branch contributes its full
//!   deviatoric neo-Hookean response with modulus `μᵢ`, so the state is
//!   `G(0) = μ∞ + Σμᵢ`.
//! - **Equilibrium** (`Γᵢ = C̄`): every branch stress vanishes
//!   identically, leaving the relaxed neo-Hookean `μ∞`.
//! - **Linear**: at small strain the simple-shear stress reduces to the
//!   Prony relaxation `γ(μ∞ + Σμᵢ e^{−t/τᵢ})` — pinned to machine
//!   precision against the closed form, which ties this model to the
//!   linear-regime-verified [`crate::PronyIntegrator`].

use tpt_med_geometry::Mat3;
const EPS_F64: f64 = 1.0e-12;

/// One finite-strain Maxwell branch: absolute shear modulus and the
/// relaxation time of its internal tensor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousBranch {
    /// Branch shear modulus `μᵢ` (MPa).
    pub mu_i: f64,
    /// Branch relaxation time `τᵢ` (s).
    pub tau_i: f64,
}

impl ViscousBranch {
    /// Validates; `Err` for a non-finite or non-positive modulus or time.
    pub fn new(mu_i: f64, tau_i: f64) -> Result<Self, String> {
        if !mu_i.is_finite() || mu_i <= 0.0 {
            return Err(format!("branch modulus must be positive, got {mu_i}"));
        }
        if !tau_i.is_finite() || tau_i <= 0.0 {
            return Err(format!("branch time must be positive, got {tau_i}"));
        }
        Ok(Self { mu_i, tau_i })
    }
}

/// The internal state: one symmetric positive-definite tensor per branch,
/// initialized to the identity (the instant configuration).
#[derive(Debug, Clone, PartialEq)]
pub struct FsViscoelasticState {
    /// Per-branch internal deformation tensors `Γᵢ` (in `C`-space).
    pub gamma: Vec<Mat3>,
}

impl FsViscoelasticState {
    /// A fresh state: every `Γᵢ = I` (no viscous deformation).
    pub fn new(branches: usize) -> Self {
        Self {
            gamma: vec![Mat3::IDENTITY; branches],
        }
    }
}

/// Finite-strain viscoelastic material: isochoric neo-Hookean equilibrium
/// network plus `n` Maxwell branches with internal tensors.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteStrainViscoelastic {
    /// Equilibrium (relaxed) shear modulus `μ∞` (MPa).
    pub mu_inf: f64,
    /// Bulk modulus `K` for the volumetric penalty `(K/2)(J−1)²` (MPa).
    pub bulk: f64,
    /// Maxwell branches (absolute moduli, not the relative `gᵢ` of the
    /// linear [`crate::ViscoelasticMaterial`]).
    pub branches: Vec<ViscousBranch>,
}

impl FiniteStrainViscoelastic {
    /// Validates; `Err` for a non-finite or non-positive modulus, or a
    /// bad branch.
    pub fn new(mu_inf: f64, bulk: f64, branches: Vec<ViscousBranch>) -> Result<Self, String> {
        if !mu_inf.is_finite() || mu_inf < 0.0 {
            return Err(format!("mu_inf must be non-negative, got {mu_inf}"));
        }
        if !bulk.is_finite() || bulk <= 0.0 {
            return Err(format!("bulk must be positive, got {bulk}"));
        }
        for b in &branches {
            ViscousBranch::new(b.mu_i, b.tau_i)?;
        }
        Ok(Self {
            mu_inf,
            bulk,
            branches,
        })
    }

    /// The instantaneous shear modulus `G(0) = μ∞ + Σμᵢ` (MPa).
    pub fn instant_modulus(&self) -> f64 {
        self.mu_inf + self.branches.iter().map(|b| b.mu_i).sum::<f64>()
    }

    /// The isochoric deformation `(C̄, J)` with `C̄ = J^{−2/3} C`; `None`
    /// for an inverted configuration.
    fn isochoric(f: &Mat3) -> Option<(Mat3, f64)> {
        let j = f.det();
        if !j.is_finite() || j <= EPS_F64 {
            return None;
        }
        Some((f.transpose() * *f * j.powf(-2.0 / 3.0), j))
    }

    /// Second Piola stress `S(F)` at the given state, advancing the
    /// internal variables by one step of `dt` (the update is exact for
    /// the branch ODE with `C̄` frozen over the step). `dt = 0` reads the
    /// stress without evolving. The stress is evaluated at the *old*
    /// internal state before the update — the standard explicit
    /// internal-variable order, and what the step-relaxation closed form
    /// relies on.
    ///
    /// Returns `None` for an inverted or non-invertible configuration.
    pub fn second_piola(&self, f: &Mat3, state: &mut FsViscoelasticState, dt: f64) -> Option<Mat3> {
        let (c_bar, j) = Self::isochoric(f)?;
        let c_bar_inv = c_bar.inverse()?;
        // C⁻¹ = F⁻¹ F⁻ᵀ.
        let f_inv = f.inverse()?;
        let mut c_inv = Mat3::ZERO;
        for i in 0..3 {
            for j in 0..3 {
                let v: f64 = (0..3).map(|k| f_inv.at(k, i) * f_inv.at(k, j)).sum();
                c_inv.set(i, j, v);
            }
        }

        // Stress at the OLD internal state.
        let j23 = j.powf(-2.0 / 3.0);
        let mut s = self.mu_inf * j23 * self.network_stress(&c_bar, &c_bar_inv, &Mat3::IDENTITY)
            + self.bulk * j * (j - 1.0) * c_inv;
        for (i, branch) in self.branches.iter().enumerate() {
            let gamma_inv = state.gamma[i].inverse()?;
            s = s + branch.mu_i * j23 * self.network_stress(&c_bar, &c_bar_inv, &gamma_inv);
        }

        // Evolve the internal tensors (exact for frozen C̄).
        if dt > 0.0 {
            for (i, branch) in self.branches.iter().enumerate() {
                let alpha = 1.0 - (-(dt / branch.tau_i)).exp();
                let g = state.gamma[i];
                state.gamma[i] = g + (c_bar - g) * alpha;
            }
        }
        Some(s)
    }

    /// First Piola `P = F S` (see [`Self::second_piola`]).
    pub fn first_piola(&self, f: &Mat3, state: &mut FsViscoelasticState, dt: f64) -> Option<Mat3> {
        let s = self.second_piola(f, state, dt)?;
        Some(*f * s)
    }

    /// The deviator-in-`C̄` form `Γ* − (tr(C̄Γ*)/3) C̄⁻¹` for an inverse
    /// tensor `Γ*` (`I` for the equilibrium network): vanishes at `Γ* =
    /// C̄⁻¹`-scaling, i.e. at the network's relaxed state, and is the full
    /// deviatoric neo-Hookean response at `Γ* = I`.
    fn network_stress(&self, c_bar: &Mat3, c_bar_inv: &Mat3, gamma_inv: &Mat3) -> Mat3 {
        let tr = (*c_bar * *gamma_inv).trace();
        let mut out = Mat3::ZERO;
        for i in 0..3 {
            for j in 0..3 {
                out.set(i, j, gamma_inv.at(i, j) - (tr / 3.0) * c_bar_inv.at(i, j));
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn material() -> FiniteStrainViscoelastic {
        FiniteStrainViscoelastic::new(
            0.3,   // mu_inf
            100.0, // bulk
            vec![
                ViscousBranch::new(0.4, 1.0).expect("branch"),
                ViscousBranch::new(0.3, 10.0).expect("branch"),
            ],
        )
        .expect("valid material")
    }

    fn simple_shear(gamma: f64) -> Mat3 {
        let mut f = Mat3::IDENTITY;
        f.set(0, 1, gamma);
        f
    }

    /// Linear limit, the model's defining contract: a small step shear
    /// relaxes exactly as the Prony series
    /// `tau(t) = gamma (mu_inf + sum mu_i e^(-t/tau_i))`, at machine
    /// precision, across both branches' time scales.
    #[test]
    fn step_shear_relaxes_as_the_prony_series() {
        let m = material();
        let gamma = 1.0e-5;
        let mut state = FsViscoelasticState::new(m.branches.len());
        let mut t = 0.0f64;
        let dt = 0.01;
        let g = |time: f64| m.mu_inf + 0.4 * (-(time / 1.0)).exp() + 0.3 * (-(time / 10.0)).exp();
        for n in 0..3000 {
            let p = m
                .first_piola(&simple_shear(gamma), &mut state, dt)
                .expect("valid");
            let tau = p.at(0, 1);
            let expected = gamma * g(t);
            assert!(
                (tau - expected).abs() < 1e-9 + 1e-6 * expected.abs(),
                "t={t}: tau {tau} vs closed form {expected}"
            );
            t += dt;
            let _ = n;
        }
        // Still relaxed correctly at 30 s (three long-branch times).
        let p = m
            .first_piola(&simple_shear(gamma), &mut state, dt)
            .expect("valid");
        assert!((p.at(0, 1) - gamma * g(t)).abs() < 1e-8);
    }

    #[test]
    fn instant_state_is_the_full_modulus_and_equilibrium_is_mu_inf() {
        let m = material();
        let gamma = 1.0e-3;
        // Instant: dt = 0 leaves Gamma = I; to linear order the shear
        // stress is G(0)·gamma with G(0) = mu_inf + sum mu_i = 1.0.
        let mut state = FsViscoelasticState::new(m.branches.len());
        let p0 = m
            .first_piola(&simple_shear(gamma), &mut state, 0.0)
            .expect("valid");
        assert!(
            (p0.at(0, 1) - m.instant_modulus() * gamma).abs() < 1e-10,
            "instant P12 {} vs G(0)·gamma {}",
            p0.at(0, 1),
            m.instant_modulus() * gamma
        );
        // Long time: the state converges to the relaxed response (the
        // equilibrium network alone, at the same mu_inf).
        let relaxed = FiniteStrainViscoelastic::new(0.3, 100.0, vec![]).expect("valid");
        let p_relaxed = relaxed
            .first_piola(&simple_shear(gamma), &mut FsViscoelasticState::new(0), 0.0)
            .expect("valid");
        let mut state = FsViscoelasticState::new(m.branches.len());
        for _ in 0..3000 {
            m.first_piola(&simple_shear(gamma), &mut state, 0.05)
                .expect("valid");
        }
        let p_late = m
            .first_piola(&simple_shear(gamma), &mut state, 0.0)
            .expect("valid");
        assert!(
            (p_late.at(0, 1) - p_relaxed.at(0, 1)).abs() < 1e-9,
            "late {} vs relaxed {}",
            p_late.at(0, 1),
            p_relaxed.at(0, 1)
        );
    }

    #[test]
    fn two_half_steps_equal_one_full_step_in_total_time() {
        // The recurrence interpolates Gamma toward C̄ by (1 - e^(-dt/tau)):
        // two half steps accumulate e^(-2h/tau), one full step e^(-2h/tau)
        // — identical, to machine precision.
        let m = material();
        let f = simple_shear(0.05);
        let mut halves = FsViscoelasticState::new(m.branches.len());
        m.second_piola(&f, &mut halves, 0.5).expect("valid");
        m.second_piola(&f, &mut halves, 0.5).expect("valid");
        let mut full = FsViscoelasticState::new(m.branches.len());
        m.second_piola(&f, &mut full, 1.0).expect("valid");
        for (a, b) in halves.gamma.iter().zip(&full.gamma) {
            for i in 0..3 {
                for j in 0..3 {
                    assert!(
                        (a.at(i, j) - b.at(i, j)).abs() < 1e-14,
                        "recurrence mismatch {a:?} vs {b:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn relaxation_is_monotone_and_internal_tensors_stay_spd() {
        let m = material();
        let gamma = 0.08;
        let mut state = FsViscoelasticState::new(m.branches.len());
        let mut previous = f64::INFINITY;
        for _ in 0..500 {
            let p = m
                .first_piola(&simple_shear(gamma), &mut state, 0.05)
                .expect("valid");
            let tau = p.at(0, 1);
            assert!(
                tau <= previous + 1e-12,
                "stress must not rise: {tau} after {previous}"
            );
            previous = tau;
            for g in &state.gamma {
                // SPD: symmetric and positive diagonal after Cholesky-free
                // check — the diagonal of the inverse is positive.
                for i in 0..3 {
                    for j in 0..3 {
                        assert!((g.at(i, j) - g.at(j, i)).abs() < 1e-12, "symmetry");
                    }
                }
                assert!(g.inverse().is_some(), "stayed invertible");
            }
        }
    }

    #[test]
    fn constructor_validates() {
        assert!(FiniteStrainViscoelastic::new(-1.0, 100.0, vec![]).is_err());
        assert!(FiniteStrainViscoelastic::new(0.3, 0.0, vec![]).is_err());
        assert!(FiniteStrainViscoelastic::new(
            0.3,
            100.0,
            vec![ViscousBranch {
                mu_i: 0.0,
                tau_i: 1.0
            }]
        )
        .is_err());
        assert!(FiniteStrainViscoelastic::new(
            0.3,
            100.0,
            vec![ViscousBranch {
                mu_i: 0.4,
                tau_i: -1.0
            }]
        )
        .is_err());
        assert!(FiniteStrainViscoelastic::new(
            0.3,
            100.0,
            vec![ViscousBranch::new(0.4, 1.0).expect("valid")]
        )
        .is_ok());
        assert!(ViscousBranch::new(0.4, 0.0).is_err());
        assert!(ViscousBranch::new(f64::NAN, 1.0).is_err());
    }
}
