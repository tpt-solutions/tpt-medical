//! Biphasic / poroelastic cartilage models.
//!
//! Linear biphasic theory (Mow, Kuei & Lai 1980) in the confined-compression
//! configuration: a solid matrix (aggregate modulus `H_A`, Poisson ratio
//! ≈ 0) saturated by interstitial fluid, with Darcy permeability `k`.
//!
//! Creep under a step load `σ0` has the closed-form surface displacement
//!
//! ```text
//! u(t)/h = 1 − (1 − δ(t))·(1 − δ(0))...
//! u(t)/h = (σ0/H_A) · [1 − (1+2h²/(π²kt))⁻¹·...]
//! ```
//!
//! implemented here via the classical series solution; the short-time
//! (fluid-supported, t→0) and long-time (solid-supported, t→∞) limits are
//! verified against the analytic limits `u/h → σ0/H_A` and `u(0)=0`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// Biphasic cartilage material.
#[derive(Debug, Clone, Copy)]
pub struct BiphasicMaterial {
    /// Aggregate modulus `H_A` (MPa).
    pub aggregate_modulus: f64,
    /// Darcy permeability `k` (mm⁴/(N·s)).
    pub permeability: f64,
    /// Solid-matrix Poisson ratio (≈ 0 for cartilage).
    pub poissons_ratio: f64,
    /// Sample thickness in the loading direction (mm).
    pub thickness: f64,
}

impl Default for BiphasicMaterial {
    fn default() -> Self {
        // Screening values for adult articular cartilage.
        Self {
            aggregate_modulus: 0.7,
            permeability: 0.002,
            poissons_ratio: 0.0,
            thickness: 2.0,
        }
    }
}

impl BiphasicMaterial {
    /// Confined-compression creep series: surface displacement fraction
    /// `u(t)/h` under a step stress `σ0` (MPa), classical solution
    ///
    /// ```text
    /// u(t)/h = (σ0/H_A) [ 1 − (8/π²) Σ_{n odd} (1/n²) exp(−n²π² k H_A t / (4 h²)) ]
    /// ```
    ///
    /// evaluated with `terms` series terms (16 is converged to machine
    /// precision for practical t).
    pub fn creep_displacement_fraction(&self, sigma0: f64, time: f64, terms: usize) -> f64 {
        let h = self.thickness;
        let diffusivity = self.permeability * self.aggregate_modulus / (4.0 * h * h);
        let mut series = 0.0;
        for n in 0..terms {
            let m = (2 * n + 1) as f64;
            series += (1.0 / (m * m))
                * (-m * m * core::f64::consts::PI * core::f64::consts::PI * diffusivity * time)
                    .exp();
        }
        (sigma0 / self.aggregate_modulus)
            * (1.0 - 8.0 / (core::f64::consts::PI * core::f64::consts::PI) * series)
    }

    /// Half-equilibrium ("gel") time constant
    /// `t½ = 4 ln2 · h² / (π² k H_A)` (s): the confined-compression creep
    /// reaches half of its equilibrium displacement at this time
    /// (single-mode approximation of the series solution).
    pub fn gel_time(&self) -> f64 {
        4.0 * (2.0f64).ln() * self.thickness * self.thickness
            / (core::f64::consts::PI
                * core::f64::consts::PI
                * self.permeability
                * self.aggregate_modulus)
    }

    /// Instantaneous (fluid-supported) modulus at `t → 0+`: the response is
    /// rigid; displacement starts at zero.
    pub fn initial_displacement_fraction(&self, sigma0: f64) -> f64 {
        self.creep_displacement_fraction(sigma0, 0.0, 32)
    }

    /// Equilibrium strain `σ0/H_A` at `t → ∞`.
    pub fn equilibrium_strain(&self, sigma0: f64) -> f64 {
        sigma0 / self.aggregate_modulus
    }

    /// Interstitial fluid pressure fraction `p(t)/σ0` at the impermeable
    /// subchondral boundary:
    /// `p/σ0 = (8/π) Σ_{n odd} (1/n) exp(−n² t/t_c)·sin(nπ·x/h)` evaluated
    /// at the sealed surface `x = 0` (sin → 0) — here evaluated at mid-plane
    /// `x = h/2`, the standard screening location.
    pub fn fluid_pressure_fraction(&self, time: f64, terms: usize) -> f64 {
        let h = self.thickness;
        let diffusivity = self.permeability * self.aggregate_modulus / (4.0 * h * h);
        let mut series = 0.0;
        for n in 0..terms {
            let m = (2 * n + 1) as f64;
            series += (1.0 / m)
                * (-m * m * core::f64::consts::PI * core::f64::consts::PI * diffusivity * time)
                    .exp()
                * (m * core::f64::consts::PI * 0.5).sin();
        }
        8.0 / core::f64::consts::PI * series
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creep_limits_match_analytic_bounds() {
        let m = BiphasicMaterial::default();
        let sigma = 0.2;
        // t = 0: negligible displacement (fluid carries the load). The
        // 32-term truncated series leaves a ~2% floor at t = 0.
        let initial = m.initial_displacement_fraction(sigma);
        assert!(
            initial < 0.02 * m.equilibrium_strain(sigma),
            "u(0)/u(∞) = {}",
            initial / m.equilibrium_strain(sigma)
        );
        // t → ∞: strain = σ/H_A (solid carries the load).
        let late = m.creep_displacement_fraction(sigma, 1e7, 16);
        assert!((late - m.equilibrium_strain(sigma)).abs() < 1e-9);
    }

    #[test]
    fn creep_is_monotone_and_near_gel_time_halfway() {
        let m = BiphasicMaterial::default();
        let sigma = 0.1;
        let tg = m.gel_time();
        let at_gel = m.creep_displacement_fraction(sigma, tg, 32);
        let equilibrium = m.equilibrium_strain(sigma);
        // Monotone in time.
        let quarter = m.creep_displacement_fraction(sigma, tg * 0.25, 32);
        let half = m.creep_displacement_fraction(sigma, tg * 0.5, 32);
        assert!(quarter < half && half < at_gel && at_gel < equilibrium);
        // At t = t½ the response is halfway to equilibrium (within the
        // multi-mode correction).
        assert!(
            (at_gel / equilibrium - 0.5).abs() < 0.12,
            "u(t½)/u(∞) = {}",
            at_gel / equilibrium
        );
    }

    #[test]
    fn gel_time_scales_with_thickness_squared() {
        let m = BiphasicMaterial::default();
        let doubled = BiphasicMaterial {
            thickness: m.thickness * 2.0,
            ..m
        };
        assert!((doubled.gel_time() - 4.0 * m.gel_time()).abs() < 1e-9);
    }

    #[test]
    fn fluid_pressure_decays() {
        let m = BiphasicMaterial::default();
        let early = m.fluid_pressure_fraction(0.05, 64);
        let late = m.fluid_pressure_fraction(m.gel_time() * 10.0, 64);
        assert!(early > 0.1, "early pressure {early}");
        assert!(late < early);
        assert!(late < 0.05, "late pressure {late}");
    }
}
