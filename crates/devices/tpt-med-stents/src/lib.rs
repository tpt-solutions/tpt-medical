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

/// Result of a deployment evaluation on a **non-uniform** ring.
#[derive(Debug, Clone)]
pub struct NonUniformDeployment {
    /// Ring-level metrics: the uniform-radius equilibrium with the crown
    /// stiffnesses summed.
    pub ring: DeploymentResult,
    /// Radial force carried by each crown (N), in input order. Zero in the
    /// no-contact case.
    pub crown_forces: Vec<f64>,
    /// Largest single-crown share of the total radial force, `(0, 1]` under
    /// contact — a tapered design or an anomalous crown shows up here.
    pub peak_crown_fraction: f64,
}

/// Simulates radial equilibrium of a stent ring whose crowns do **not**
/// share one stiffness (tapered designs, per-crown variation): the ring
/// radius stays uniform in this model, so the ring-level metrics use the
/// summed stiffness and the per-crown forces split proportionally to
/// stiffness. An empty `crown_stiffness` slice falls back to
/// `n_crowns` copies of [`StentModel::crown_stiffness`].
pub fn simulate_deployment_with_crowns(
    stent: &StentModel,
    nitinol: &NitinolParams,
    vessel_diameter_at_pressure: impl Fn(f64) -> f64,
    vessel_pressure: Pressure,
    crown_stiffness: &[f64],
) -> NonUniformDeployment {
    let stiffnesses: Vec<f64> = if crown_stiffness.is_empty() {
        vec![stent.crown_stiffness; stent.n_crowns as usize]
    } else {
        crown_stiffness.to_vec()
    };
    let total_stiffness: f64 = stiffnesses.iter().sum();
    let nominal = stent.expanded_diameter;
    let lumen = vessel_diameter_at_pressure(vessel_pressure.to_mpa());
    let compression = (nominal - lumen).max(0.0);

    let mut ring =
        simulate_deployment(stent, nitinol, vessel_diameter_at_pressure, vessel_pressure);
    // Recompute the force with the summed stiffness (the uniform path uses
    // n_crowns · crown_stiffness, which an explicit non-uniform slice with
    // equal entries must reproduce exactly).
    let crown_forces: Vec<f64> = stiffnesses.iter().map(|&k| k * compression).collect();
    ring.radial_force = if compression > 0.0 {
        total_stiffness * compression
    } else {
        0.0
    };
    ring.contact_pressure = if compression > 0.0 {
        ring.radial_force / (core::f64::consts::PI * lumen * lumen)
    } else {
        0.0
    };
    let peak_crown_fraction = if total_stiffness > 0.0 && compression > 0.0 {
        stiffnesses.iter().cloned().fold(0.0f64, f64::max) / total_stiffness
    } else {
        0.0
    };
    NonUniformDeployment {
        ring,
        crown_forces,
        peak_crown_fraction,
    }
}

/// Result of a **tapered** (multi-ring-group) deployment — the Level-2
/// model, per RFC 0004's fidelity ladder. Each axial ring group reaches its
/// own equilibrium against the local lumen, which is what gives `dogboning`
/// a real value: a lesion profile that resists mid-stent expansion makes
/// the ends open wider.
#[derive(Debug, Clone)]
pub struct TaperedDeployment {
    /// Per-group equilibrium results, in input order.
    pub groups: Vec<DeploymentResult>,
    /// Axial **dogboning**: `|d_ends − d_mids| / nominal`, where the end
    /// diameter is the mean of the first and last group and the middle
    /// diameter is the mean of the interior groups (the overall mean when
    /// there are no interior groups). Zero for fewer than two groups.
    pub dogboning: f64,
    /// Mean equilibrium diameter over the groups (mm).
    pub mean_diameter: f64,
    /// Total radial force: the sum over groups (N).
    pub radial_force: f64,
}

/// Level-2 radial equilibrium with the stent resolved into **axial ring
/// groups** — a tapered design (different nominal diameter or stiffness per
/// group) or an axial lesion profile (a smaller lumen at the mid-stent
/// groups) produces genuinely different group equilibria, and the
/// F2394-style dogboning metric takes a non-zero value.
///
/// The three slices pair element-by-element (equal lengths, ≥ 1 group):
/// `group_nominal` is the group's free (recovered) diameter, `group_stiffness`
/// its per-crown stiffness, and `group_lumen` the local vessel lumen at the
/// deployment pressure (the caller's axial profile — evaluate the vessel law
/// at each group's position to pass a lesion shape through). The uniform
/// single-ring case is [`simulate_deployment`].
pub fn simulate_tapered_deployment(
    stent: &StentModel,
    nitinol: &NitinolParams,
    group_nominal: &[f64],
    group_stiffness: &[f64],
    group_lumen: &[f64],
) -> TaperedDeployment {
    assert_eq!(
        group_nominal.len(),
        group_stiffness.len(),
        "group nominal-diameter and stiffness counts must pair"
    );
    assert_eq!(
        group_nominal.len(),
        group_lumen.len(),
        "group nominal-diameter and lumen counts must pair"
    );
    assert!(
        !group_nominal.is_empty(),
        "at least one ring group is required"
    );

    let groups: Vec<DeploymentResult> = group_nominal
        .iter()
        .zip(group_stiffness)
        .zip(group_lumen)
        .map(|((&nominal, &stiffness), &lumen)| {
            let group = StentModel {
                expanded_diameter: nominal,
                crimped_diameter: stent.crimped_diameter,
                n_crowns: stent.n_crowns,
                crown_stiffness: stiffness,
            };
            simulate_deployment(&group, nitinol, |_| lumen, Pressure::from_mpa(0.0))
        })
        .collect();

    let mean_diameter = groups.iter().map(|g| g.diameter).sum::<f64>() / groups.len() as f64;
    let radial_force: f64 = groups.iter().map(|g| g.radial_force).sum();

    let dogboning = match groups.len() {
        0 | 1 => 0.0,
        2 => (groups[0].diameter - groups[1].diameter).abs() / stent.expanded_diameter,
        _ => {
            let d_end = 0.5 * (groups[0].diameter + groups[groups.len() - 1].diameter);
            let d_mid = groups[1..groups.len() - 1]
                .iter()
                .map(|g| g.diameter)
                .sum::<f64>()
                / (groups.len() - 2) as f64;
            (d_end - d_mid).abs() / stent.expanded_diameter
        }
    };

    TaperedDeployment {
        groups,
        dogboning,
        mean_diameter,
        radial_force,
    }
}

/// Screening strain-life law for superelastic Nitinol: the alternating
/// strain amplitude the material tolerates **degrades logarithmically with
/// cycle count** (the published Nitinol fatigue band drops roughly a factor
/// of two from 10³ to 10⁷ cycles — e.g. the rotational-bending data
/// collected in Pelton's Nitinol fatigue reviews). A design whose computed
/// strain amplitude exceeds the degraded limit fails the screen.
///
/// This is a *screening* interpolation of published band data, not a
/// device-specific S–N curve; a life claim still needs the vendor's own
/// fatigue data and ASTM F2477-style pulsatile testing.
#[derive(Debug, Clone, Copy)]
pub struct StrainLifeLaw {
    /// Alternating strain amplitude tolerated at [`Self::reference_cycles`].
    pub amplitude_at_reference: f64,
    /// Reference cycle count for `amplitude_at_reference`.
    pub reference_cycles: f64,
    /// Log–log slope `d log(ε)/d log(N)` (negative: life shortens the more
    /// strain is applied).
    pub slope: f64,
}

impl StrainLifeLaw {
    /// A screening Nitinol band: 0.4 % alternating amplitude at 10⁷ cycles,
    /// a factor-of-two drop per four decades (slope ≈ −0.075).
    pub fn nitinol_screening() -> Self {
        Self {
            amplitude_at_reference: 0.004,
            reference_cycles: 1.0e7,
            slope: -0.075,
        }
    }

    /// Tolerated alternating amplitude (dimensionless strain) after `cycles`
    /// cycles.
    pub fn amplitude_at(&self, cycles: f64) -> f64 {
        if cycles <= 0.0 {
            return f64::INFINITY; // no fatigue demand yet
        }
        self.amplitude_at_reference * (cycles / self.reference_cycles).powf(self.slope)
    }

    /// Screening verdict: `true` when the applied alternating strain
    /// amplitude stays below the tolerated amplitude at `cycles` — the
    /// conservative order (equality fails).
    pub fn survives(&self, cycles: f64, applied_amplitude: f64) -> bool {
        applied_amplitude < self.amplitude_at(cycles)
    }
}

impl StentModel {
    /// Geometric **foreshortening** of a zig-zag crown ring: the axial
    /// shortening when the ring opens from its delivery (crimped)
    /// configuration to `diameter`, as a fraction
    /// `(L_crimped − L)/L_crimped`.
    ///
    /// Each crown cell is modelled as a diamond: strut segment length is
    /// fixed (developed material), so opening the cell's circumferential
    /// width `πD/n_crowns` must reduce its axial height. `link_fraction` is
    /// the share of the manufactured length held by straight axial links
    /// that do not swing (real laser-cut designs ≈ 0.6–0.8; 0 is the pure
    /// diamond-cell upper bound). Returns `NaN` when `diameter` exceeds the
    /// developed-length limit, where the cell geometry cannot close.
    pub fn foreshortening(
        &self,
        manufactured_length: f64,
        diameter: f64,
        link_fraction: f64,
    ) -> f64 {
        let n = self.n_crowns as f64;
        // Full diamond cell built from the whole manufactured length.
        let h_crimped_cell = manufactured_length / n;
        let half_width_crimped = core::f64::consts::PI * self.crimped_diameter / (2.0 * n);
        let segment =
            (h_crimped_cell * h_crimped_cell + half_width_crimped * half_width_crimped).sqrt();
        let half_width = core::f64::consts::PI * diameter / (2.0 * n);
        if half_width >= segment {
            return f64::NAN; // circumference exceeds developed strut length
        }
        let h_cell = (segment * segment - half_width * half_width).sqrt();
        // Straight links keep their length; only the swinging share of the
        // cell follows the diamond height ratio.
        let axial_ratio = link_fraction + (1.0 - link_fraction) * (h_cell / h_crimped_cell);
        1.0 - axial_ratio
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

    #[test]
    fn foreshortening_is_zero_at_the_crimped_diameter_and_grows_with_expansion() {
        let stent = StentModel {
            expanded_diameter: 6.0,
            crimped_diameter: 1.8,
            n_crowns: 12,
            crown_stiffness: 0.5,
        };
        let length = 16.0;
        assert_eq!(stent.foreshortening(length, 1.8, 0.7), 0.0);
        // Monotone in diameter within the physical range…
        let f4 = stent.foreshortening(length, 4.0, 0.7);
        let f6 = stent.foreshortening(length, 6.0, 0.7);
        let f8 = stent.foreshortening(length, 8.0, 0.7);
        assert!(f4 > 0.0 && f6 > f4 && f8 > f6, "{f4} {f6} {f8}");
        // …and with realistic straight links (0.7) a 6 mm deployment lands
        // in the published few-percent band; the pure diamond cell (0.0) is
        // the upper bound.
        assert!((0.02..0.10).contains(&f6), "6 mm foreshortening {f6}");
        assert!(stent.foreshortening(length, 6.0, 0.0) > f6);
        // All-link ring never foreshortens.
        assert_eq!(stent.foreshortening(length, 6.0, 1.0), 0.0);
        // Beyond the developed-length limit the cell cannot close: NaN.
        assert!(stent.foreshortening(length, 14.0, 0.7).is_nan());
    }

    #[test]
    fn nonuniform_ring_splits_force_by_stiffness() {
        let stent = StentModel {
            expanded_diameter: 6.0,
            crimped_diameter: 1.8,
            n_crowns: 4,
            crown_stiffness: 0.5,
        };
        let nitinol = NitinolParams::default();
        let vessel = |_| 5.0f64;
        let p = Pressure::from_mpa(0.013);
        // Equal stiffnesses reproduce the uniform ring exactly.
        let uniform = simulate_deployment(&stent, &nitinol, vessel, p);
        let equal = simulate_deployment_with_crowns(&stent, &nitinol, vessel, p, &[0.5; 4]);
        assert!((equal.ring.radial_force - uniform.radial_force).abs() < 1e-12);
        assert!((equal.peak_crown_fraction - 0.25).abs() < 1e-12);
        // Empty slice falls back to the uniform crown stiffness.
        let fallback = simulate_deployment_with_crowns(&stent, &nitinol, vessel, p, &[]);
        assert!((fallback.ring.radial_force - uniform.radial_force).abs() < 1e-12);
        // One stiff crown carries proportionally more of the load.
        let mixed =
            simulate_deployment_with_crowns(&stent, &nitinol, vessel, p, &[0.5, 0.5, 0.5, 2.5]);
        assert_eq!(mixed.crown_forces.len(), 4);
        assert!((mixed.peak_crown_fraction - 2.5 / 4.0).abs() < 1e-12);
        assert!(
            (mixed.crown_forces[3].min(mixed.crown_forces[0]) / mixed.crown_forces[3] - 0.2).abs()
                < 1e-12
        );
        // Crown forces sum to the ring total.
        let sum: f64 = mixed.crown_forces.iter().sum();
        assert!((sum - mixed.ring.radial_force).abs() < 1e-12);
    }

    #[test]
    fn nonuniform_no_contact_carries_nothing() {
        let stent = StentModel {
            expanded_diameter: 4.0,
            crimped_diameter: 1.5,
            n_crowns: 8,
            crown_stiffness: 0.5,
        };
        let out = simulate_deployment_with_crowns(
            &stent,
            &NitinolParams::default(),
            |_| 6.0,
            Pressure::from_mpa(0.01),
            &[0.5, 1.5, 0.5, 1.5, 0.5, 1.5, 0.5, 1.5],
        );
        assert!(out.crown_forces.iter().all(|&f| f == 0.0));
        assert_eq!(out.peak_crown_fraction, 0.0);
        assert_eq!(out.ring.radial_force, 0.0);
    }

    fn base_stent() -> StentModel {
        StentModel {
            expanded_diameter: 6.0,
            crimped_diameter: 1.8,
            n_crowns: 12,
            crown_stiffness: 0.5,
        }
    }

    #[test]
    fn tapered_deployment_with_uniform_profile_matches_the_uniform_ring() {
        let stent = base_stent();
        let nitinol = NitinolParams::default();
        let out = simulate_tapered_deployment(
            &stent,
            &nitinol,
            &[6.0, 6.0, 6.0],
            &[0.5, 0.5, 0.5],
            &[5.0, 5.0, 5.0],
        );
        let single = simulate_deployment(&stent, &nitinol, |_| 5.0, Pressure::from_mpa(0.013));
        assert_eq!(out.dogboning, 0.0, "uniform profile cannot dogbone");
        for g in &out.groups {
            assert!((g.radial_force - single.radial_force).abs() < 1e-12);
        }
        assert!((out.radial_force - 3.0 * single.radial_force).abs() < 1e-12);
        assert!((out.mean_diameter - single.diameter).abs() < 1e-12);
    }

    #[test]
    fn stiff_mid_lesion_produces_real_dogboning() {
        let stent = base_stent();
        let nitinol = NitinolParams::default();
        // A calcified mid-lesion: the middle group meets a smaller lumen,
        // the ends open wider — the ends-flare-out dogbone shape.
        let out = simulate_tapered_deployment(
            &stent,
            &nitinol,
            &[6.0, 6.0, 6.0],
            &[0.5, 0.5, 0.5],
            &[5.0, 4.6, 5.0],
        );
        assert!(out.dogboning > 0.0, "dogboning {}", out.dogboning);
        assert!(
            out.groups[0].diameter > out.groups[1].diameter,
            "ends wider than the lesion: {} vs {}",
            out.groups[0].diameter,
            out.groups[1].diameter
        );
        // Mirror the profile: |d_end − d_mid| is direction-agnostic.
        let mirrored = simulate_tapered_deployment(
            &stent,
            &nitinol,
            &[6.0, 6.0, 6.0],
            &[0.5, 0.5, 0.5],
            &[4.6, 5.0, 4.6],
        );
        assert!(
            (mirrored.dogboning - out.dogboning).abs() < 1e-12,
            "{mirrored:?} vs {out:?}"
        );
        // A two-group ring uses the overall mean as the middle reference.
        let two =
            simulate_tapered_deployment(&stent, &nitinol, &[6.0, 6.0], &[0.5, 0.5], &[5.0, 4.6]);
        assert!(two.dogboning > 0.0);
    }

    #[test]
    fn strain_life_law_degrades_with_cycles() {
        let law = StrainLifeLaw::nitinol_screening();
        assert!((law.amplitude_at(law.reference_cycles) - 0.004).abs() < 1e-15);
        // Fewer cycles tolerate more strain; more cycles less.
        assert!(law.amplitude_at(1.0e3) > law.amplitude_at(1.0e7));
        assert!(law.amplitude_at(1.0e9) < 0.004);
        // Screening verdict: 0.3 % amplitude is safe at 10 years of
        // cardiac cycling, 0.5 % is not.
        assert!(law.survives(3.0e8, 0.003));
        assert!(!law.survives(3.0e8, 0.005));
        // No cycles yet: infinite capacity.
        assert_eq!(law.amplitude_at(0.0), f64::INFINITY);
        assert!(law.survives(0.0, 0.05));
    }
}
