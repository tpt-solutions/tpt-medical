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
