//! Blood rheology models.

/// Shear-rate-dependent blood viscosity models (μ in Pa·s, γ̇ in s⁻¹).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BloodModel {
    /// Constant viscosity.
    Newtonian {
        /// Dynamic viscosity (Pa·s); ~0.0035 Pa·s for normal hematocrit.
        viscosity: f64,
    },
    /// Carreau–Yasuda shear-thinning model:
    /// `μ = μ∞ + (μ0 − μ∞) [1 + (λγ̇)^a]^((n−1)/a)`.
    CarreauYasuda {
        /// Zero-shear viscosity (Pa·s).
        mu0: f64,
        /// Infinite-shear viscosity (Pa·s).
        mu_inf: f64,
        /// Time constant λ (s).
        lambda: f64,
        /// Transition sharpness `a`.
        a: f64,
        /// Power-law index `n` (< 1 for shear thinning).
        n: f64,
    },
    /// Casson yield model: `√τ = √τ_y + √(μ γ̇)`.
    Casson {
        /// Yield stress τ_y (Pa); ~0.005 Pa for blood.
        tau_y: f64,
        /// Casson viscosity slope (Pa·s).
        mu: f64,
    },
}

impl BloodModel {
    /// Classic 3.5 cP Newtonian blood.
    pub const NEWTONIAN_BLOOD: BloodModel = BloodModel::Newtonian { viscosity: 0.0035 };

    /// Human blood Carreau–Yasuda parameter set (literature screening
    /// values: μ0 = 0.056 Pa·s, μ∞ = 0.0035 Pa·s, λ = 11.5 s, a = 1.9,
    /// n = 0.22).
    pub const CARREAU_YASUDA_BLOOD: BloodModel = BloodModel::CarreauYasuda {
        mu0: 0.056,
        mu_inf: 0.0035,
        lambda: 11.5,
        a: 1.9,
        n: 0.22,
    };

    /// Apparent viscosity at a local shear rate.
    pub fn viscosity(&self, shear_rate: f64) -> f64 {
        match *self {
            BloodModel::Newtonian { viscosity } => viscosity,
            BloodModel::CarreauYasuda {
                mu0,
                mu_inf,
                lambda,
                a,
                n,
            } => {
                let g = shear_rate.max(1e-6);
                mu_inf + (mu0 - mu_inf) * (1.0 + (lambda * g).powf(a)).powf((n - 1.0) / a)
            }
            BloodModel::Casson { tau_y, mu } => {
                // Solve √τ = √τy + √(μ γ̇); τ = μ_app γ̇ gives
                // μ_app = τ_y/γ̇ + μ + 2√(τ_y μ / γ̇).
                let g = shear_rate.max(1e-6);
                tau_y / g + mu + 2.0 * (tau_y * mu / g).sqrt()
            }
        }
    }

    /// Plasma reference viscosity used in relative-viscosity reporting.
    pub const PLASMA: f64 = 0.0012;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newtonian_is_constant() {
        let m = BloodModel::NEWTONIAN_BLOOD;
        for g in [0.1, 10.0, 1000.0] {
            assert!((m.viscosity(g) - 0.0035).abs() < 1e-12);
        }
    }

    #[test]
    fn carreau_yasuda_shear_thins() {
        let m = BloodModel::CARREAU_YASUDA_BLOOD;
        let low = m.viscosity(1.0);
        let mid = m.viscosity(100.0);
        let high = m.viscosity(10_000.0);
        assert!(low > mid && mid > high, "{low} {mid} {high}");
        // High-shear plateau near μ∞.
        assert!((high - 0.0035).abs() < 0.002);
        // Low-shear viscosity well above μ∞.
        assert!(low > 0.01, "low-shear μ = {low}");
    }

    #[test]
    fn casson_limits() {
        let m = BloodModel::Casson {
            tau_y: 0.005,
            mu: 0.0035,
        };
        // High shear: apparent viscosity → μ.
        assert!((m.viscosity(1e6) - 0.0035).abs() < 1e-4);
        // Low shear: yield term dominates → large apparent viscosity.
        assert!(m.viscosity(1e-3) > 5.0);
    }
}
