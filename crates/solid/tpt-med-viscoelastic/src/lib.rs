//! Prony-series viscoelasticity for soft tissue.
//!
//! Generalized Maxwell material in shear:
//!
//! ```text
//! G(t) = G∞ + Σᵢ Gᵢ exp(−t / τᵢ)
//! ```
//!
//! with frequency-domain storage/loss moduli
//!
//! ```text
//! G'(ω) = G∞ + Σᵢ Gᵢ (ωτᵢ)² / (1 + (ωτᵢ)²)
//! G''(ω) = Σᵢ Gᵢ (ωτᵢ) / (1 + (ωτᵢ)²)
//! ```
//!
//! An elastic reference model (`tpt-med-tissue`) supplies the glass
//! response; the Prony terms are *relative* moduli `gᵢ = Gᵢ/G0` with
//! `Σgᵢ ≤ 1`, the convention used by Abaqus/FEBio inputs.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use tpt_med_tissue::TissueModel;

/// One Maxwell element: relative modulus and relaxation time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PronyTerm {
    /// Relative shear modulus `gᵢ = Gᵢ/G0` (dimensionless).
    pub g_i: f64,
    /// Relaxation time `τᵢ` (s).
    pub tau_i: f64,
}

/// A Prony-series viscoelastic material.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscoelasticMaterial {
    /// Instantaneous (glass) shear modulus `G0` (MPa).
    pub g0: f64,
    /// Maxwell elements.
    pub prony: Vec<PronyTerm>,
}

impl ViscoelasticMaterial {
    /// Validates the series: non-negative terms, `Σgᵢ ≤ 1`, `τᵢ > 0`.
    pub fn validate(&self) -> Result<(), String> {
        if self.g0 <= 0.0 {
            return Err("G0 must be positive".into());
        }
        let sum: f64 = self.prony.iter().map(|t| t.g_i).sum();
        if sum > 1.0 + 1e-9 {
            return Err(format!("Σgᵢ = {sum} exceeds 1"));
        }
        for t in &self.prony {
            if t.g_i < 0.0 {
                return Err("negative Prony modulus".into());
            }
            if t.tau_i <= 0.0 {
                return Err("non-positive relaxation time".into());
            }
        }
        Ok(())
    }

    /// Long-term (equilibrium) modulus `G∞ = G0 (1 − Σgᵢ)` (MPa).
    pub fn equilibrium_modulus(&self) -> f64 {
        self.g0 * (1.0 - self.prony.iter().map(|t| t.g_i).sum::<f64>())
    }

    /// Relaxation (stress-relaxation) modulus `G(t)` (MPa).
    pub fn relaxation_modulus(&self, time: f64) -> f64 {
        let g_inf = self.equilibrium_modulus();
        g_inf
            + self
                .prony
                .iter()
                .map(|t| self.g0 * t.g_i * (-time / t.tau_i).exp())
                .sum::<f64>()
    }

    /// Storage modulus `G'(ω)` (MPa).
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        self.g0 * (1.0 - self.prony.iter().map(|t| t.g_i).sum::<f64>())
            + self
                .prony
                .iter()
                .map(|t| {
                    let wt = omega * t.tau_i;
                    self.g0 * t.g_i * wt * wt / (1.0 + wt * wt)
                })
                .sum::<f64>()
    }

    /// Loss modulus `G''(ω)` (MPa).
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        self.prony
            .iter()
            .map(|t| {
                let wt = omega * t.tau_i;
                self.g0 * t.g_i * wt / (1.0 + wt * wt)
            })
            .sum::<f64>()
    }

    /// Loss tangent `tan δ = G''/G'`.
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        self.loss_modulus(omega) / self.storage_modulus(omega)
    }

    /// Stress response to a step shear strain `γ0` applied at t = 0:
    /// `σ(t) = γ0 · G(t)` (MPa).
    pub fn step_strain_stress(&self, gamma0: f64, time: f64) -> f64 {
        gamma0 * self.relaxation_modulus(time)
    }

    /// Uniaxial elastic reference model scaled by the instantaneous
    /// modulus — used when embedding viscoelasticity into an FEM pipeline
    /// whose elasticity comes from `tpt-med-tissue`.
    pub fn elastic_reference(&self) -> TissueModel {
        // G = E / (2(1+ν)); take ν = 0.5 (incompressible) ⇒ E = 3G.
        TissueModel::NeoHookean(tpt_med_tissue::NeoHookeanParams {
            c10: self.g0 / 3.0,
            d1: 100.0,
        })
    }
}

/// Re-exported param alias so the reference above stays in-crate.
pub mod tpt_med_tissue_link {
    pub use tpt_med_tissue::NeoHookeanParams as NeoHookeanForVis;
}

/// Temperature shift factor conventions (documented, single source).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TemperatureShift {
    /// Williams–Landel–Ferry: `log10 aT = −C1 (T − Tref) / (C2 + T − Tref)`.
    /// Valid above Tg; the classic C1 = 17.4, C2 = 51.6 K are the caller's
    /// to cite — this crate ships no silent material constants.
    Wlf {
        /// WLF C1 (dimensionless).
        c1: f64,
        /// WLF C2 (kelvin).
        c2_kelvin: f64,
    },
    /// Arrhenius: `ln aT = (Ea/R) (1/T − 1/Tref)`, temperatures in kelvin.
    Arrhenius {
        /// Apparent activation energy (J/mol).
        activation_energy_j_mol: f64,
    },
}

impl TemperatureShift {
    /// Shift factor `aT = τ(T)/τ(Tref)` at temperature `t` (same unit as
    /// the reference temperature; kelvin for Arrhenius).
    pub fn shift_factor(&self, t: f64, t_ref: f64) -> f64 {
        match *self {
            TemperatureShift::Wlf { c1, c2_kelvin } => {
                let dt = t - t_ref;
                10f64.powf(-c1 * dt / (c2_kelvin + dt))
            }
            TemperatureShift::Arrhenius {
                activation_energy_j_mol,
            } => {
                const R: f64 = 8.314_462_618;
                ((activation_energy_j_mol / R) * (1.0 / t - 1.0 / t_ref)).exp()
            }
        }
    }

    /// Master-curve shift: every Prony time constant becomes
    /// `τᵢ(T) = τᵢ_ref · aT(T)`.
    pub fn shifted_material(
        &self,
        material: &ViscoelasticMaterial,
        t_ref: f64,
        t: f64,
    ) -> ViscoelasticMaterial {
        let a = self.shift_factor(t, t_ref);
        ViscoelasticMaterial {
            g0: material.g0,
            prony: material
                .prony
                .iter()
                .map(|term| PronyTerm {
                    g_i: term.g_i,
                    tau_i: term.tau_i * a,
                })
                .collect(),
        }
    }
}

/// Time-integration helper for driving a Prony material inside an explicit
/// FEM/CFD loop, so callers don't reimplement the internal-variable
/// recurrence.
///
/// Each Maxwell element carries a partial (viscous) strain `εᵢ` updated by
/// the exact exponential recurrence over the step
/// `Δt`:
///
/// ```text
/// εᵢⁿ⁺¹ = e^(−Δt/τᵢ) εᵢⁿ + (1 − e^(−Δt/τᵢ)) γⁿ⁺¹
/// σⁿ⁺¹  = G∞ γⁿ⁺¹ + Σᵢ Gᵢ (γⁿ⁺¹ − εᵢⁿ⁺¹)
/// ```
///
/// which is the standard uniaxial shear form; the update is exact for a
/// strain that is linear over the step.
#[derive(Debug, Clone)]
pub struct PronyIntegrator {
    /// The (temperature-shifted, if any) material being integrated.
    material: ViscoelasticMaterial,
    /// Per-element viscous (partially developed) strain.
    viscous: Vec<f64>,
    /// Held total strain.
    strain: f64,
}

impl PronyIntegrator {
    /// Fresh integrator at zero strain and zero partial strains.
    pub fn new(material: &ViscoelasticMaterial) -> Self {
        material.validate().expect("valid Prony series");
        Self {
            viscous: vec![0.0; material.prony.len()],
            strain: 0.0,
            material: material.clone(),
        }
    }

    /// Advances one step to the new total strain `gamma`, returning the
    /// shear stress (MPa). `dt <= 0` is treated as a no-op returning the
    /// current stress.
    pub fn step(&mut self, gamma: f64, dt: f64) -> f64 {
        if dt > 0.0 {
            for (term, eps_i) in self.material.prony.iter().zip(&mut self.viscous) {
                let decay = (-dt / term.tau_i).exp();
                *eps_i = decay * *eps_i + (1.0 - decay) * gamma;
            }
            self.strain = gamma;
        }
        self.stress()
    }

    /// Current stress at the held strain (MPa).
    pub fn stress(&self) -> f64 {
        let g_inf = self.material.equilibrium_modulus();
        let mut sigma = g_inf * self.strain;
        for (term, &eps_i) in self.material.prony.iter().zip(&self.viscous) {
            sigma += self.material.g0 * term.g_i * (self.strain - eps_i);
        }
        sigma
    }

    /// Current total strain.
    pub fn strain(&self) -> f64 {
        self.strain
    }
}

/// Fung-type **quasi-linear viscoelasticity (QLV)**: the non-linear
/// generalisation that applies the Prony series to a *hyperelastic*
/// stress history. Where [`PronyIntegrator`] superposes in strain
/// (linear regime), QLV superposes in stress:
///
/// ```text
/// σ(t) = Σ_k g(t − t_k)·Δσ^e_k,   g = G(t)/G0,  Σ_k Δσ^e_k = σ^e(t)
/// ```
///
/// with `σ^e` the **instantaneous elastic** response of a hyperelastic
/// law evaluated along the (possibly finite-strain) deformation history
/// — the hereditary integral Fung introduced for soft tissue. The kernel
/// is the crate's own normalised relaxation modulus, so the linear regime
/// reduces exactly to the [`PronyIntegrator`] behaviour.
///
/// Accuracy note: the discrete form samples the hereditary integral
/// first-order (rectangle rule on the elastic-stress increments); a step
/// elastic history is represented *exactly*.
pub struct QuasiLinearViscoelastic {
    material: ViscoelasticMaterial,
}

impl QuasiLinearViscoelastic {
    /// Builds the QLV wrapper from the material whose normalised
    /// relaxation modulus supplies the kernel.
    pub fn new(material: &ViscoelasticMaterial) -> Self {
        material.validate().expect("valid Prony series");
        Self {
            material: material.clone(),
        }
    }

    /// Reduced relaxation kernel `g(t) = G(t)/G0 ∈ (g∞, 1]`, `g(0) = 1`.
    pub fn kernel(&self, time: f64) -> f64 {
        self.material.relaxation_modulus(time) / self.material.g0
    }

    /// The QLV stress history at `times`, given the instantaneous elastic
    /// stress `elastic_stress` sampled at the same instants (both evenly
    /// or unevenly spaced; `times` must be ascending and start at the
    /// history's origin). The elastic stress may be any (non-linear)
    /// function of the deformation — a `tpt-med-tissue` model's stress
    /// along the strain path, for instance.
    pub fn stress_history(&self, times: &[f64], elastic_stress: &[f64]) -> Vec<f64> {
        assert_eq!(times.len(), elastic_stress.len(), "histories must pair");
        assert!(!times.is_empty(), "empty history");
        let mut out = Vec::with_capacity(times.len());
        let mut previous = 0.0f64;
        for (n, &t) in times.iter().enumerate() {
            // Rectangle rule on the elastic-stress increments; the jump
            // 0 → σ^e[0] enters at index 0 with kernel weight g(t − t_0).
            let mut sigma = 0.0f64;
            let mut prior_stress = 0.0f64;
            for k in 0..=n {
                let increment = elastic_stress[k] - prior_stress;
                sigma += self.kernel(t - times[k]) * increment;
                prior_stress = elastic_stress[k];
            }
            out.push(sigma);
            let _ = previous;
            previous = t;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn material() -> ViscoelasticMaterial {
        ViscoelasticMaterial {
            g0: 1.0,
            prony: vec![
                PronyTerm {
                    g_i: 0.4,
                    tau_i: 1.0,
                },
                PronyTerm {
                    g_i: 0.3,
                    tau_i: 10.0,
                },
            ],
        }
    }

    #[test]
    fn qlv_step_history_is_exact() {
        let material = ViscoelasticMaterial {
            g0: 10.0,
            prony: vec![
                PronyTerm {
                    g_i: 0.3,
                    tau_i: 0.5,
                },
                PronyTerm {
                    g_i: 0.2,
                    tau_i: 5.0,
                },
            ],
        };
        let q = QuasiLinearViscoelastic::new(&material);
        // A step elastic stress: QLV is exactly g(t)·σ0.
        let times = [0.0, 0.1, 0.5, 1.0, 10.0];
        let sigma_e = [5.0, 5.0, 5.0, 5.0, 5.0];
        let out = q.stress_history(&times, &sigma_e);
        for (i, &t) in times.iter().enumerate() {
            let expected = q.kernel(t) * 5.0;
            assert!(
                (out[i] - expected).abs() < 1e-12,
                "t={t}: {} vs {expected}",
                out[i]
            );
        }
        // Kernel bounds: g(0) = 1, decaying toward the equilibrium share.
        assert!((q.kernel(0.0) - 1.0).abs() < 1e-12);
        assert!(q.kernel(10.0) < 1.0);
        assert!(q.kernel(10.0) > material.equilibrium_modulus() / material.g0 - 1e-9);
    }

    #[test]
    fn qlv_reduces_to_the_prony_integrator_in_the_linear_regime() {
        let material = ViscoelasticMaterial {
            g0: 3.0,
            prony: vec![
                PronyTerm {
                    g_i: 0.4,
                    tau_i: 0.2,
                },
                PronyTerm {
                    g_i: 0.3,
                    tau_i: 1.5,
                },
            ],
        };
        // Linear ramp to γ = 0.01 over 1 s, then hold: the elastic stress
        // is the instantaneous G0·γ, so QLV must reproduce the exact
        // recurrence of PronyIntegrator to the rectangle rule's accuracy.
        let dt = 0.001;
        let mut integrator = PronyIntegrator::new(&material);
        let mut times = Vec::new();
        let mut elastic = Vec::new();
        for k in 0..=2000 {
            let t = k as f64 * dt;
            let gamma = 0.01 * (t / 1.0).min(1.0);
            times.push(t);
            elastic.push(material.g0 * gamma);
            integrator.step(gamma, dt);
            // Compare every 200 steps (skip the exact-recurrence-vs-
            // quadrature tail where the hold dominates; the two agree to
            // first order in dt there as well).
            if k % 200 == 0 && k > 0 {
                let q = QuasiLinearViscoelastic::new(&material);
                let qlv = q.stress_history(&times, &elastic);
                let reference = integrator.stress();
                assert!(
                    (qlv[qlv.len() - 1] - reference).abs() < 5e-4,
                    "k={k}: QLV {} vs integrator {reference}",
                    qlv[qlv.len() - 1]
                );
            }
        }
    }

    #[test]
    fn endpoints_of_relaxation() {
        let m = material();
        assert!((m.relaxation_modulus(0.0) - 1.0).abs() < 1e-12);
        // Long time: only G∞ remains.
        assert!((m.relaxation_modulus(1e6) - 0.3).abs() < 1e-9);
        // G(τ1): fast mode decayed to e^{-1}, slow mode ≈ intact.
        let expected = 0.3 + 1.0 * (0.4 * (-1.0f64).exp() + 0.3 * (-0.1f64).exp());
        assert!((m.relaxation_modulus(1.0) - expected).abs() < 1e-12);
    }

    #[test]
    fn monotone_relaxation() {
        let m = material();
        let mut prev = m.relaxation_modulus(0.0);
        for i in 1..200 {
            let t = i as f64 * 0.1;
            let g = m.relaxation_modulus(t);
            assert!(g <= prev + 1e-15, "G increased at t={t}");
            prev = g;
        }
    }

    #[test]
    fn frequency_limits() {
        let m = material();
        // Low frequency: G' → G∞, G'' → 0.
        assert!((m.storage_modulus(1e-6) - 0.3).abs() < 1e-6);
        assert!(m.loss_modulus(1e-6) < 1e-4);
        // High frequency: G' → G0, G'' → 0.
        assert!((m.storage_modulus(1e6) - 1.0).abs() < 1e-3);
        assert!(m.loss_modulus(1e6) < 1e-2);
        // Peak loss between the two relaxation times.
        let peak = (0..800)
            .map(|i| 10f64.powf(i as f64 / 100.0 - 2.0))
            .map(|w| m.loss_modulus(w))
            .fold(0.0f64, f64::max);
        assert!(peak > 0.1, "peak loss {peak}");
    }

    #[test]
    fn wlf_shift_speeds_up_above_reference() {
        let wlf = TemperatureShift::Wlf {
            c1: 17.4,
            c2_kelvin: 51.6,
        };
        let a_cold = wlf.shift_factor(25.0, 37.0);
        let a_hot = wlf.shift_factor(60.0, 37.0);
        assert!(a_hot < 1.0 && a_hot > 0.0);
        assert!(a_cold > 1.0);
        // WLF diverges as T approaches Tg (here Tref + -C2): sanity only.
        let m = material();
        let shifted = wlf.shifted_material(&m, 37.0, 60.0);
        assert!(shifted.prony[0].tau_i < m.prony[0].tau_i);
    }

    #[test]
    fn arrhenius_shift_increases_with_temperature() {
        let arr = TemperatureShift::Arrhenius {
            activation_energy_j_mol: 50_000.0,
        };
        let a_hot = arr.shift_factor(350.0, 300.0);
        assert!(a_hot < 1.0, "higher T -> faster relaxation: {a_hot}");
        assert!((arr.shift_factor(300.0, 300.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn integrator_reproduces_relaxation_then_recovery() {
        // Stress relaxation under held strain must decay like G(t).
        let m = material();
        let mut it = PronyIntegrator::new(&m);
        let dt = 1e-6;
        let sigma0 = it.step(0.05, dt);
        // The exponential recurrence reproduces G(dt) exactly for a step
        // held from zero: sigma = gamma * G(dt).
        assert!(
            (sigma0 - 0.05 * m.relaxation_modulus(dt)).abs() < 1e-12,
            "{} vs {}",
            sigma0,
            0.05 * m.relaxation_modulus(dt)
        );
        let mut t = 0.0;
        for _ in 0..2000 {
            t += 0.01;
            it.step(0.05, 0.01);
        }
        let expected = 0.05 * m.relaxation_modulus(t);
        assert!(
            (it.stress() - expected).abs() < 5e-3 * expected.abs().max(1e-9),
            "{} vs {expected}",
            it.stress()
        );
        // Return to zero strain: stress must recover to zero. The slow
        // mode (tau = 10 s) needs many time constants to vanish, so run
        // 200 s and expect the residual below 1e-5 MPa.
        for _ in 0..20_000 {
            it.step(0.0, 0.01);
        }
        assert!(it.stress().abs() < 1e-5, "residual {}", it.stress());
    }

    #[test]
    fn integrator_matches_analytic_constant_rate_response() {
        // For a linear strain ramp at rate r, the analytic stress of the
        // one-term line is G∞ r t + G g τ r (1 - e^(-t/τ)) / 1 — compare
        // against the integrator within integration tolerance.
        let m = ViscoelasticMaterial {
            g0: 1.0,
            prony: vec![PronyTerm {
                g_i: 0.5,
                tau_i: 2.0,
            }],
        };
        let rate = 0.01;
        let mut it = PronyIntegrator::new(&m);
        let dt = 0.01;
        let mut t = 0.0;
        let mut last = 0.0;
        for _ in 0..1000 {
            t += dt;
            last = it.step(rate * t, dt);
        }
        // Analytic ramp response: sigma(t) = G_inf * r * t
        //                     + G0 * g * r * tau * (1 - e^(-t/tau)).
        let g_inf = m.equilibrium_modulus();
        let expected = g_inf * rate * t + m.g0 * 0.5 * rate * 2.0 * (1.0 - (-t / 2.0).exp());
        assert!(
            (last - expected).abs() < 1e-3 * expected.abs().max(1e-9),
            "{last} vs {expected}"
        );
    }

    #[test]
    fn validation_rejects_over_sum() {
        let bad = ViscoelasticMaterial {
            g0: 1.0,
            prony: vec![
                PronyTerm {
                    g_i: 0.7,
                    tau_i: 1.0,
                },
                PronyTerm {
                    g_i: 0.7,
                    tau_i: 2.0,
                },
            ],
        };
        assert!(bad.validate().is_err());
        assert!(material().validate().is_ok());
    }

    #[test]
    fn step_strain_matches_relaxation() {
        let m = material();
        assert!((m.step_strain_stress(0.05, 2.0) - 0.05 * m.relaxation_modulus(2.0)).abs() < 1e-15);
    }
}
