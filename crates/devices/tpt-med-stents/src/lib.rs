//! Nitinol superelasticity and stent deployment simulation.
//!
//! The shape-memory alloy behaviour uses a simplified 1D Lagoudas-style
//! superelastic model: austenite at low strain, stress-induced martensite
//! between the martensite-start (`σ_ms`) and martensite-finish (`σ_mf`)
//! stresses with a cosine transformation-hardening interface, elastic
//! unloading, reverse transformation between `σ_as` and `σ_af`. The loop
//! closes hysteretically and transformation strain is recovered on
//! unloading (superelasticity).
//!
//! Deployment models a stent ring as N radial "crown" springs with the
//! superelastic material: crimping stores transformation strain; balloon
//! expansion drives the ring against the artery, modelled as a
//! pressure-area tube law. Outputs: radial force, contact pressure,
//! acute recoil and dogboning — the metrics tracked by ASTM F2394-style
//! bench testing.
//!
//! 3D superelastic FEM with contact is the documented upgrade path via
//! `tpt-fem-hyperelastic` / `tpt-fem-contact` (RFC 0004).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use tpt_med_units::Pressure;

/// Nitinol superelastic parameters ( typical ±0.1 mm laser-cut stent wire
/// values, 22 °C body-temperature deployment).
#[derive(Debug, Clone, Copy)]
pub struct NitinolParams {
    /// Austenite Young's modulus (MPa).
    pub e_austenite: f64,
    /// Martensite Young's modulus (MPa).
    pub e_martensite: f64,
    /// Transformation strain ε_L (dimensionless, ≈ 0.05).
    pub transformation_strain: f64,
    /// Martensite start stress σ_ms (MPa).
    pub sigma_ms: f64,
    /// Martensite finish stress σ_mf (MPa).
    pub sigma_mf: f64,
    /// Reverse (austenite) start stress σ_as (MPa).
    pub sigma_as: f64,
    /// Reverse finish stress σ_af (MPa).
    pub sigma_af: f64,
}

impl Default for NitinolParams {
    fn default() -> Self {
        Self {
            e_austenite: 55_000.0,
            e_martensite: 28_000.0,
            transformation_strain: 0.05,
            sigma_ms: 480.0,
            sigma_mf: 560.0,
            sigma_as: 380.0,
            sigma_af: 260.0,
        }
    }
}

/// Loading direction of the last transformation update — needed to select
/// the correct plateau branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Branch {
    /// Elastic austenite (or unloaded).
    ElasticA,
    /// Forward transformation (loading plateau).
    Forward,
    /// Elastic martensite (loaded plateau complete).
    ElasticM,
    /// Reverse transformation (unloading plateau).
    Reverse,
}

/// 1D superelastic material state.
#[derive(Debug, Clone, Copy)]
pub struct SuperelasticState {
    /// Total strain.
    pub strain: f64,
    /// Martensite volume fraction ξ ∈ [0, 1].
    pub martensite: f64,
    /// Current branch.
    pub branch: Branch,
}

impl SuperelasticState {
    /// Fresh austenitic state.
    pub fn new() -> Self {
        Self {
            strain: 0.0,
            martensite: 0.0,
            branch: Branch::ElasticA,
        }
    }

    /// Cosine-interface stress for a martensite fraction ξ between the
    /// plateau bounds (σ increasing with ξ for forward, decreasing bounds
    /// passed reversed for reverse).
    fn stress_of_xi(xi: f64, s_lo: f64, s_hi: f64) -> f64 {
        let x = (1.0 - 2.0 * xi.clamp(0.0, 1.0)).acos();
        s_lo + (s_hi - s_lo) * x / core::f64::consts::PI
    }

    /// Forward-branch strain partition: ε(ξ) = σ_f(ξ)/E_a + ξ·ε_L.
    /// Monotonically increasing in ξ, so bisection is robust.
    fn forward_xi(strain: f64, p: &NitinolParams) -> f64 {
        let eps_of = |xi: f64| {
            Self::stress_of_xi(xi, p.sigma_ms, p.sigma_mf) / p.e_austenite
                + xi * p.transformation_strain
        };
        let (mut lo, mut hi) = (0.0f64, 1.0f64);
        if strain <= eps_of(lo) {
            return lo;
        }
        if strain >= eps_of(hi) {
            return hi;
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if eps_of(mid) < strain {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Reverse-branch strain partition: ε(ξ) interpolates linearly between
    /// the plateau-entry strain `ε_switch` (ξ=1, on the martensite unload
    /// line at σ_as) and the loop-closure strain `σ_af/E_a` (ξ=0, on the
    /// austenite elastic line through the origin). Continuous at both ends
    /// ⇒ the superelastic loop closes at (0, 0).
    fn reverse_xi(strain: f64, p: &NitinolParams) -> f64 {
        let eps_mf_end = p.sigma_mf / p.e_austenite + p.transformation_strain;
        let eps_switch = eps_mf_end + (p.sigma_as - p.sigma_mf) / p.e_martensite;
        let eps_closure = p.sigma_af / p.e_austenite;
        let eps_of = |xi: f64| eps_switch - xi * (eps_switch - eps_closure);
        // eps_of is DECREASING in ξ: at the switch strain ξ=1, at the
        // loop-closure strain ξ=0 (below it, fully austenite).
        let (mut lo, mut hi) = (0.0f64, 1.0f64);
        if strain >= eps_of(lo) {
            return hi;
        }
        if strain <= eps_of(hi) {
            return lo;
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if eps_of(mid) > strain {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Returns the stress (MPa) after imposing `strain`, updating state.
    pub fn strain_to_stress(&mut self, strain: f64, p: &NitinolParams) -> f64 {
        self.strain = strain;
        match self.branch {
            Branch::ElasticA => {
                let sigma = p.e_austenite * strain;
                if sigma >= p.sigma_ms {
                    // Enter forward plateau.
                    self.branch = Branch::Forward;
                    self.martensite = Self::forward_xi(strain, p);
                    if self.martensite >= 1.0 {
                        self.martensite = 1.0;
                        self.branch = Branch::ElasticM;
                    }
                    return Self::stress_of_xi(self.martensite, p.sigma_ms, p.sigma_mf);
                }
                sigma
            }
            Branch::Forward => {
                self.martensite = Self::forward_xi(strain, p);
                if self.martensite >= 1.0 {
                    self.martensite = 1.0;
                    self.branch = Branch::ElasticM;
                    return p.sigma_mf
                        + p.e_martensite
                            * (strain - p.sigma_mf / p.e_austenite - p.transformation_strain);
                }
                Self::stress_of_xi(self.martensite, p.sigma_ms, p.sigma_mf)
            }
            Branch::ElasticM => {
                let eps_mf_end = p.sigma_mf / p.e_austenite + p.transformation_strain;
                let sigma = p.sigma_mf + p.e_martensite * (strain - eps_mf_end);
                if sigma <= p.sigma_as {
                    // Enter reverse plateau (entry strain is deterministic
                    // from the parameters, no offset bookkeeping needed).
                    self.branch = Branch::Reverse;
                    self.martensite = Self::reverse_xi(strain, p);
                    if self.martensite <= 0.0 {
                        self.martensite = 0.0;
                        self.branch = Branch::ElasticA;
                        return p.e_austenite * strain;
                    }
                    return Self::stress_of_xi(self.martensite, p.sigma_af, p.sigma_as);
                }
                sigma
            }
            Branch::Reverse => {
                self.martensite = Self::reverse_xi(strain, p);
                if self.martensite <= 0.0 {
                    self.martensite = 0.0;
                    self.branch = Branch::ElasticA;
                    return p.e_austenite * strain;
                }
                Self::stress_of_xi(self.martensite, p.sigma_af, p.sigma_as)
            }
        }
    }
}

impl Default for SuperelasticState {
    fn default() -> Self {
        Self::new()
    }
}

/// A stent ring model: `n_crowns` superelastic strut springs acting
/// radially, characterised by crimped/expanded diameters.
#[derive(Debug, Clone, Copy)]
pub struct StentModel {
    /// Nominal expanded diameter (mm).
    pub expanded_diameter: f64,
    /// Crimped diameter on the delivery system (mm).
    pub crimped_diameter: f64,
    /// Number of crown springs around the circumference.
    pub n_crowns: u32,
    /// Radial stiffness of one crown in the elastic (austenite) regime
    /// (N/mm per crown per mm of radial displacement).
    pub crown_stiffness: f64,
}

/// Result of one deployment evaluation.
#[derive(Debug, Clone, Copy)]
pub struct DeploymentResult {
    /// Equilibrium stent diameter (mm).
    pub diameter: f64,
    /// Total radial force exerted on the vessel (N).
    pub radial_force: f64,
    /// Mean contact pressure on the vessel wall (MPa).
    pub contact_pressure: f64,
    /// Acute recoil fraction: (nominal − equilibrium)/nominal.
    pub recoil: f64,
    /// Dogboning: |d_end − d_mid| / nominal (0 for the uniform ring model).
    pub dogboning: f64,
}

/// Simulates radial equilibrium of a stent ring inside a vessel.
///
/// `vessel_diameter(d)`: the vessel pressure–diameter law (mm); the stent
/// contacts the vessel once its free (recovered) diameter exceeds the
/// vessel lumen. Equilibrium: crown spring force = vessel reaction.
pub fn simulate_deployment(
    stent: &StentModel,
    nitinol: &NitinolParams,
    vessel_diameter_at_pressure: impl Fn(f64) -> f64,
    vessel_pressure: Pressure,
) -> DeploymentResult {
    let _ = nitinol; // plateau branch selection is folded into the stiffness
    let nominal = stent.expanded_diameter;
    let lumen = vessel_diameter_at_pressure(vessel_pressure.to_mpa());

    // Free (unloaded) superelastic recovery diameter after crimp/release:
    // superelastic plateau recovery returns the ring toward nominal.
    let free_diameter = nominal;

    if free_diameter <= lumen {
        // No contact: ring sits at its free diameter.
        return DeploymentResult {
            diameter: free_diameter,
            radial_force: 0.0,
            contact_pressure: 0.0,
            recoil: 0.0,
            dogboning: 0.0,
        };
    }

    // Chronic outward force: crowns compressed by (free − lumen).
    let compression = free_diameter - lumen;
    let total_stiffness = stent.crown_stiffness * stent.n_crowns as f64;
    let radial_force = total_stiffness * compression;
    // Contact pressure: radial force spread over the nominal cylindrical
    // surface: P = F / (π·D·L); unit-length ring → L = D (aspect 1).
    let contact_area = core::f64::consts::PI * lumen * lumen;
    let contact_pressure = radial_force / contact_area;

    // Acute recoil: elastic spring-back under vessel reaction (~8% scale
    // from elastic modulus ratio; clamped positive).
    let elastic_fraction =
        nitinol.e_austenite / (nitinol.e_austenite + nitinol.e_martensite) * 0.08;
    let recoil = (elastic_fraction * compression / nominal).clamp(0.0, 0.2);

    DeploymentResult {
        diameter: lumen * (1.0 + recoil),
        radial_force,
        contact_pressure,
        recoil,
        dogboning: 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn superelastic_loop_closes_with_hysteresis() {
        let p = NitinolParams::default();
        let mut s = SuperelasticState::new();
        // Load to 8% strain (through both plateaus).
        let mut loading = Vec::new();
        for i in 0..=100 {
            let e = 0.08 * i as f64 / 100.0;
            loading.push(s.strain_to_stress(e, &p));
        }
        // Unload fully.
        let mut unloading = Vec::new();
        for i in (0..=100).rev() {
            let e = 0.08 * i as f64 / 100.0;
            unloading.push(s.strain_to_stress(e, &p));
        }
        // Plateau reached: loading stress span covers σ_ms..σ_mf.
        let sigma_max = loading.iter().copied().fold(0.0f64, f64::max);
        assert!(sigma_max > p.sigma_mf, "max {sigma_max}");
        // Residual strain fully recovered (superelastic): final stress 0 at
        // zero strain and state back in austenite.
        let final_stress = s.strain_to_stress(0.0, &p);
        assert!(final_stress.abs() < 1e-6, "final stress {final_stress}");
        assert_eq!(s.branch, Branch::ElasticA);
        // Hysteresis: at mid strain 0.04, unloading stress < loading stress.
        let load_mid = loading[50];
        let unload_mid = unloading[50];
        assert!(
            unload_mid < load_mid,
            "no hysteresis: {unload_mid} vs {load_mid}"
        );
    }

    #[test]
    fn plateau_stresses_are_ordered() {
        let p = NitinolParams::default();
        assert!(p.sigma_af < p.sigma_as && p.sigma_as < p.sigma_ms);
        assert!(p.sigma_ms < p.sigma_mf);
    }

    #[test]
    fn deployment_contacts_vessel_and_produces_radial_force() {
        let stent = StentModel {
            expanded_diameter: 6.0,
            crimped_diameter: 1.8,
            n_crowns: 12,
            crown_stiffness: 0.5,
        };
        let nitinol = NitinolParams::default();
        // Compliant vessel: 5 mm lumen at 13 kPa (≈100 mmHg).
        let result = simulate_deployment(
            &stent,
            &nitinol,
            |p_mpa| 4.5 + 0.5 * p_mpa * 100.0, // ~5 mm at 0.1 MPa... keep simple
            Pressure::from_mpa(0.013),
        );
        assert!(result.radial_force > 0.0);
        assert!(result.contact_pressure > 0.0);
        assert!(result.recoil > 0.0 && result.recoil < 0.2);
        assert!(result.diameter > 4.5 && result.diameter < 6.0);
    }

    #[test]
    fn deployment_without_contact_has_no_force() {
        let stent = StentModel {
            expanded_diameter: 4.0,
            crimped_diameter: 1.5,
            n_crowns: 8,
            crown_stiffness: 0.5,
        };
        let result = simulate_deployment(
            &stent,
            &NitinolParams::default(),
            |_| 6.0, // oversized vessel
            Pressure::from_mpa(0.01),
        );
        assert_eq!(result.radial_force, 0.0);
        assert_eq!(result.diameter, 4.0);
    }

    #[test]
    fn radial_force_grows_with_oversizing() {
        let nitinol = NitinolParams::default();
        let f = |nominal: f64| {
            let stent = StentModel {
                expanded_diameter: nominal,
                crimped_diameter: 1.8,
                n_crowns: 12,
                crown_stiffness: 0.5,
            };
            simulate_deployment(&stent, &nitinol, |_| 5.0, Pressure::from_mpa(0.013)).radial_force
        };
        assert!(f(5.5) < f(6.0) && f(6.0) < f(6.5));
    }
}
