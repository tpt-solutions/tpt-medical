//! Implant–bone micromotion and stress-shielding analysis for
//! orthopedic devices.
//!
//! - **Micromotion**: relative interface displacement between implant and
//!   bone under load, computed with a Winkler-foundation interface model
//!   (bone as distributed springs, implant as rigid punch). Osseointegration
//!   screening uses the classic threshold: interface motion < 150 µm favours
//!   bone fixation; the stricter < 50 µm band favours primary stability for
//!   cementless press-fit components (workspace convention).
//! - **Stress shielding**: per-zone ratio of bone strain energy with the
//!   implant in place versus the intact baseline (Engh-style index).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use tpt_med_geometry::Vec3;
use tpt_med_units::{Force, Length};

/// Interface foundation parameters.
#[derive(Debug, Clone, Copy)]
pub struct InterfaceModel {
    /// Bone foundation stiffness per unit area (N/mm³) — cortical ~2–20,
    /// trabecular ~0.2–2.
    pub foundation_stiffness: f64,
    /// Effective contact area of the implant–bone interface (mm²).
    pub contact_area: f64,
    /// Friction coefficient of the interface (press-fit titanium ~0.4–0.6).
    pub friction: f64,
}

/// Osseointegration risk classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RiskLevel {
    /// Motion below the primary-stability band.
    Low,
    /// Motion within the osseointegration band.
    Moderate,
    /// Motion above the fibrous-tissue threshold.
    High,
}

/// Result of a micromotion analysis.
#[derive(Debug, Clone)]
pub struct MicromotionResult {
    /// Maximum relative interface displacement (mm).
    pub max_micromotion: f64,
    /// Per-zone displacements (mm), in input zone order.
    pub zone_micromotion: Vec<f64>,
    /// Osseointegration risk classification.
    pub risk: RiskLevel,
}

impl MicromotionResult {
    /// Thresholds: < 0.05 mm low, < 0.15 mm moderate, else high (mm).
    pub fn classify(max_micromotion: f64) -> RiskLevel {
        if max_micromotion < 0.050 {
            RiskLevel::Low
        } else if max_micromotion < 0.150 {
            RiskLevel::Moderate
        } else {
            RiskLevel::High
        }
    }
}

/// Winkler-foundation interface micromotion under a joint reaction force.
///
/// The load introduces an interface shear `τ = μ·σ + F_tangential/A`; the
/// relative motion follows `δ = τ / k` per zone. Zones with zero contact
/// area are skipped.
pub fn micromotion_analysis(
    interface: &InterfaceModel,
    joint_reaction: Force,
    tangential_fraction: f64,
    zone_areas: &[f64],
) -> MicromotionResult {
    let shear_load = joint_reaction.to_n() * tangential_fraction.clamp(0.0, 1.0);
    let mut zone_micromotion = Vec::with_capacity(zone_areas.len());
    let mut max = 0.0f64;
    for &area in zone_areas {
        if area <= 0.0 {
            zone_micromotion.push(f64::NAN);
            continue;
        }
        // Normal pressure from the total reaction over the zone.
        let pressure = joint_reaction.to_n() / area;
        // Interface shear resistance: friction + mechanical interlock.
        let shear_stress = shear_load / area + interface.friction * pressure * 0.1;
        let delta = shear_stress / interface.foundation_stiffness / 1.0e3; // N/mm³ → mm
        zone_micromotion.push(delta);
        max = max.max(delta);
    }
    MicromotionResult {
        max_micromotion: max,
        zone_micromotion,
        risk: MicromotionResult::classify(max),
    }
}

/// Stress-shielding assessment per Gruen-like zone.
#[derive(Debug, Clone)]
pub struct StressShieldingResult {
    /// Zone strain-energy densities with implant (MPa = mJ/mm³).
    pub with_implant: Vec<f64>,
    /// Zone strain-energy densities of the intact baseline (MPa).
    pub intact: Vec<f64>,
    /// Shielding ratio per zone: 1 − (SED_implant / SED_intact), clamped to
    /// [0, 1]; higher = more shielded.
    pub shielding_index: Vec<f64>,
}

impl StressShieldingResult {
    /// Mean shielding index over zones.
    pub fn mean_index(&self) -> f64 {
        if self.shielding_index.is_empty() {
            return 0.0;
        }
        self.shielding_index.iter().sum::<f64>() / self.shielding_index.len() as f64
    }

    /// True if any proximal zone is severely shielded (> 0.7) — the classic
    /// resorption-risk flag.
    pub fn has_resorption_risk(&self) -> bool {
        self.shielding_index.iter().any(|&s| s > 0.7)
    }
}

/// Computes the shielding index from paired zone energies.
pub fn stress_shielding_analysis(
    intact_sed: &[f64],
    implanted_sed: &[f64],
) -> StressShieldingResult {
    assert_eq!(intact_sed.len(), implanted_sed.len(), "zone count mismatch");
    let shielding_index = intact_sed
        .iter()
        .zip(implanted_sed)
        .map(|(&base, &impl_)| {
            if base <= 0.0 {
                0.0
            } else {
                (1.0 - impl_ / base).clamp(0.0, 1.0)
            }
        })
        .collect();
    StressShieldingResult {
        with_implant: implanted_sed.to_vec(),
        intact: intact_sed.to_vec(),
        shielding_index,
    }
}

/// Convenience: micromotion in micrometres for reporting.
pub fn micromotion_um(result: &MicromotionResult) -> Length {
    Length::from_mm(result.max_micromotion * 1.0e3)
}

/// Effective interface stiffness when the implant itself is **compliant**:
/// the Winkler foundation in series with the implant's coating/stem
/// stiffness, so the micromotion the bone sees includes the implant's own
/// elastic deformation — a screening step from rigid-punch toward
/// compliant-implant evaluation (full continuum coupling is the
/// fem-adapter's job).
///
/// ```text
/// 1/k_eff = 1/k_foundation + t/(E·A_share)
/// ```
///
/// with `t` the implant's load-bearing thickness and `E` its modulus
/// (Ti-alloy stem ~110 GPa; a porous-coating or cement mantle drops the
/// effective modulus by an order of magnitude, which is exactly the case
/// where the rigid-punch assumption bites).
#[derive(Debug, Clone, Copy)]
pub struct CompliantImplant {
    /// Foundation stiffness of the bone bed (N/mm³).
    pub foundation_stiffness: f64,
    /// Implant modulus at the interface (MPa).
    pub implant_modulus_mpa: f64,
    /// Implant load-bearing thickness at the interface (mm).
    pub implant_thickness_mm: f64,
}

impl CompliantImplant {
    /// Effective foundation stiffness (N/mm³) seen by the micromotion
    /// model: the series combination of the bone bed and the implant's
    /// interface layer. A rigid punch is the `implant_modulus → ∞` limit,
    /// where this returns `foundation_stiffness` unchanged.
    pub fn effective_stiffness(&self) -> f64 {
        let k_implant = self.implant_modulus_mpa / self.implant_thickness_mm.max(1e-9);
        // MPa = N/mm²; a foundation stiffness is N/mm³ = N/mm²/mm, so the
        // implant layer stiffness k_i = E/t converts directly.
        let k_eff = 1.0 / (1.0 / self.foundation_stiffness + 1.0 / k_implant);
        let _ = k_implant;
        k_eff
    }

    /// An [`InterfaceModel`] with the foundation stiffness replaced by the
    /// series effective value — ready for [`micromotion_analysis`] or
    /// [`micromotion_over_cycle`].
    pub fn interface_model(&self, contact_area: f64, friction: f64) -> InterfaceModel {
        InterfaceModel {
            foundation_stiffness: self.effective_stiffness(),
            contact_area,
            friction,
        }
    }
}

/// Result of a cyclic (gait) micromotion analysis.
#[derive(Debug, Clone)]
pub struct CyclicMicromotionResult {
    /// Maximum interface displacement over the cycle (mm) — the peak
    /// instantaneous value, classified as in the static case.
    pub peak_micromotion: f64,
    /// Per-zone motion *amplitude* over the cycle, `max − min` (mm): the
    /// relative interface motion that cyclic loading actually imposes,
    /// which is what fibrous-tissue screening keys on.
    pub zone_amplitude: Vec<f64>,
    /// Largest per-zone amplitude (mm).
    pub peak_amplitude: f64,
    /// Osseointegration risk classification, on the peak.
    pub risk: RiskLevel,
}

/// Micromotion accumulated over one load cycle (e.g. a gait cycle) rather
/// than at a single static load: the interface sees every load sample in
/// turn, and both the peak motion and the per-zone motion amplitude are
/// reported. Load samples are typically produced by
/// [`GaitCycle::iso_double_hump`].
pub fn micromotion_over_cycle(
    interface: &InterfaceModel,
    loads: &[Force],
    tangential_fraction: f64,
    zone_areas: &[f64],
) -> CyclicMicromotionResult {
    let mut peak = 0.0f64;
    let mut zone_min = vec![f64::INFINITY; zone_areas.len()];
    let mut zone_max = vec![f64::NEG_INFINITY; zone_areas.len()];
    for load in loads {
        let r = micromotion_analysis(interface, *load, tangential_fraction, zone_areas);
        peak = peak.max(r.max_micromotion);
        for ((zmin, zmax), z) in zone_min
            .iter_mut()
            .zip(zone_max.iter_mut())
            .zip(r.zone_micromotion)
        {
            *zmin = zmin.min(z);
            *zmax = zmax.max(z);
        }
    }
    let zone_amplitude: Vec<f64> = zone_min
        .iter()
        .zip(&zone_max)
        .map(|(lo, hi)| {
            if *lo == f64::INFINITY {
                f64::NAN // skipped (zero-area) zone
            } else {
                hi - lo
            }
        })
        .collect();
    let peak_amplitude = zone_amplitude
        .iter()
        .filter(|a| !a.is_nan())
        .fold(0.0f64, |a, &x| a.max(x));
    CyclicMicromotionResult {
        peak_micromotion: peak,
        peak_amplitude,
        zone_amplitude,
        risk: MicromotionResult::classify(peak),
    }
}

/// A canonical axial gait loading profile, in the ISO 14243 (knee-wear
/// testing) *shape*: a double hump — heel-strike peak, stance trough,
/// push-off peak — normalised so the peak load is a multiple of body
/// weight. This is a screening waveform, not a patient's measured gait;
/// the standard's full force/alignment tables are not reproduced.
#[derive(Debug, Clone, Copy)]
pub struct GaitCycle {
    /// Peak load as a multiple of body weight (ISO 14243-style: ≈ 2.6).
    pub peak_bw: f64,
    /// Cycle length (s). Only affects the reported sample times.
    pub cycle_seconds: f64,
    /// Body weight the profile scales to (N).
    pub body_weight_n: f64,
}

impl GaitCycle {
    /// The ISO 14243-style double hump: peak ≈ 2.6 × body weight.
    pub fn iso_double_hump(body_weight_n: f64) -> Self {
        Self {
            peak_bw: 2.6,
            cycle_seconds: 1.0,
            body_weight_n,
        }
    }

    /// Load samples over one cycle. `t ∈ [0, 1)` is cycle fraction: the
    /// profile is two raised-cosine humps at 15 % and 45 % cycle fraction,
    /// with a stance trough between them and swing unload after.
    pub fn samples(&self, n: usize) -> Vec<(f64, Force)> {
        let peak = self.peak_bw * self.body_weight_n;
        let trough = 0.3 * peak;
        let mut v = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / n as f64;
            let load = hump(t, 0.15, 0.18, peak)
                .max(hump(t, 0.45, 0.18, peak))
                .max(trough * 0.5);
            v.push((t * self.cycle_seconds, Force::from_n(load)));
        }
        v
    }

    /// Force samples only (convenience for [`micromotion_over_cycle`]).
    pub fn load_samples(&self, n: usize) -> Vec<Force> {
        self.samples(n).into_iter().map(|(_, f)| f).collect()
    }
}

// Raised-cosine hump centred at `centre` with half-width `width`, peak `peak`.
fn hump(t: f64, centre: f64, width: f64, peak: f64) -> f64 {
    let x = (t - centre) / width;
    if x.abs() < 1.0 {
        peak * 0.5 * (1.0 + (core::f64::consts::PI * x).cos())
    } else {
        0.0
    }
}

/// Screening migration law: the time-dependent consequence of interface
/// micromotion. Per-cycle migration accrues proportionally to the motion
/// **amplitude** above a stability threshold, with the rate decaying
/// exponentially as the implant beds in:
///
/// ```text
/// dx/dN = k·(δ_amp − δ_th)·e^{−x/x_bed}   ⇒   x(N) = x_bed·ln(1 + k·(δ_amp − δ_th)·N / x_bed)
/// ```
///
/// which reproduces the classic logarithmic migration curve seen in
/// radiostereometric analysis (RSA): rapid early bedding-in, then a slow
/// creep whose velocity is the at-risk discriminator.
#[derive(Debug, Clone, Copy)]
pub struct MigrationModel {
    /// Per-cycle rate constant (mm per cycle, per mm of excess amplitude).
    pub rate: f64,
    /// Motion-amplitude stability threshold (mm): below it, no migration
    /// accrues (the interface is in the osseointegration band).
    pub threshold: f64,
    /// Bedding-in length scale (mm) — the exponential decay of the rate
    /// with accumulated migration.
    pub bedding_in: f64,
}

impl MigrationModel {
    /// Instantaneous per-cycle migration rate (mm/cycle) at a given motion
    /// amplitude and accumulated migration.
    pub fn rate_at(&self, amplitude: f64, migration: f64) -> f64 {
        if amplitude <= self.threshold {
            return 0.0;
        }
        let decay = (-(migration / self.bedding_in)).exp();
        self.rate * (amplitude - self.threshold) * decay
    }

    /// Closed-form cumulative migration (mm) after `cycles` cycles at a
    /// constant motion amplitude.
    pub fn cumulative(&self, amplitude: f64, cycles: f64) -> f64 {
        if amplitude <= self.threshold || cycles <= 0.0 {
            return 0.0;
        }
        let a = self.rate * (amplitude - self.threshold);
        self.bedding_in * (1.0 + a * cycles / self.bedding_in).ln()
    }

    /// Late migration velocity (mm/year) at a given accumulated migration —
    /// the RSA discriminator: sustained velocity after the bedding-in year
    /// flags at-risk fixation.
    pub fn velocity_per_year(&self, amplitude: f64, migration: f64, cycles_per_year: f64) -> f64 {
        self.rate_at(amplitude, migration) * cycles_per_year
    }

    /// RSA-style stability screening: `true` when the current migration
    /// velocity exceeds 0.2 mm/year (continued migration after the first
    /// year) — the at-risk band in RSA follow-up practice.
    pub fn is_at_risk(&self, amplitude: f64, migration: f64, cycles_per_year: f64) -> bool {
        self.velocity_per_year(amplitude, migration, cycles_per_year) > 0.2
    }
}

/// Gruen zones for femoral stem fixation assessment (Gruen, McNeice &
/// Amstutz 1979): seven periprosthetic regions — three lateral (1–3,
/// proximal to distal), the distal tip (4), and three medial (5–7, distal
/// to proximal). The classifier is geometric: a point is projected onto
/// the stem axis and, off-axis, assigned by side and axial third.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GruenZone {
    /// Lateral, proximal third.
    Z1,
    /// Lateral, middle third.
    Z2,
    /// Lateral, distal third.
    Z3,
    /// Distal tip.
    Z4,
    /// Medial, distal third.
    Z5,
    /// Medial, middle third.
    Z6,
    /// Medial, proximal third.
    Z7,
}

impl GruenZone {
    /// Stable key for reports and audit records.
    pub fn key(self) -> &'static str {
        match self {
            GruenZone::Z1 => "gruen-1",
            GruenZone::Z2 => "gruen-2",
            GruenZone::Z3 => "gruen-3",
            GruenZone::Z4 => "gruen-4",
            GruenZone::Z5 => "gruen-5",
            GruenZone::Z6 => "gruen-6",
            GruenZone::Z7 => "gruen-7",
        }
    }
}

/// Classifies a point into its Gruen zone.
///
/// * `stem_proximal` / `stem_tip`: axis end points (proximal shoulder and
///   distal tip of the stem, patient coordinates).
/// * `mediolateral`: a direction pointing toward the *lateral* side at the
///   stem (e.g. from the medial cortex toward the greater trochanter).
/// * `point`: the centroid of the evaluated region.
/// * `tip_band`: fraction of stem length adjacent to the tip that counts
///   as zone 4 (typical 0.1–0.15).
///
/// On-axis points (within `tip_band` of the tip) are zone 4; otherwise the
/// lateral/medial side is decided by the mediolateral projection and the
/// axial position by thirds of the stem length.
pub fn gruen_zone(
    stem_proximal: Vec3,
    stem_tip: Vec3,
    mediolateral: Vec3,
    point: Vec3,
    tip_band: f64,
) -> GruenZone {
    let axis = stem_tip - stem_proximal;
    let length = axis.norm().max(1e-9);
    let axial = axis / length;
    let ml = mediolateral.normalize();

    let rel = point - stem_proximal;
    let t = (rel.dot(axial) / length).clamp(0.0, 1.0); // 0 = proximal, 1 = tip
    let side = rel.dot(ml);

    // Tip band wins over side classification.
    if t >= 1.0 - tip_band {
        return GruenZone::Z4;
    }
    let lateral = side >= 0.0;
    let third = t < 1.0 / 3.0;
    let middle = t < 2.0 / 3.0;
    if lateral {
        if third {
            GruenZone::Z1
        } else if middle {
            GruenZone::Z2
        } else {
            GruenZone::Z3
        }
    } else if third {
        GruenZone::Z7
    } else if middle {
        GruenZone::Z6
    } else {
        GruenZone::Z5
    }
}

/// All seven zones in canonical reporting order (1..=7).
pub fn gruen_zones_in_order() -> [GruenZone; 7] {
    [
        GruenZone::Z1,
        GruenZone::Z2,
        GruenZone::Z3,
        GruenZone::Z4,
        GruenZone::Z5,
        GruenZone::Z6,
        GruenZone::Z7,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interface() -> InterfaceModel {
        InterfaceModel {
            foundation_stiffness: 5.0,
            contact_area: 800.0,
            friction: 0.5,
        }
    }

    #[test]
    fn gruen_zones_classify_by_side_and_third() {
        // Stem along z from proximal (0,0,10) to tip (0,0,0); lateral = +x.
        let prox = Vec3::new(0.0, 0.0, 10.0);
        let tip = Vec3::new(0.0, 0.0, 0.0);
        let ml = Vec3::new(1.0, 0.0, 0.0);
        let cases = [
            (Vec3::new(1.0, 0.0, 9.0), GruenZone::Z1), // lateral proximal
            (Vec3::new(1.0, 0.0, 5.0), GruenZone::Z2), // lateral middle
            (Vec3::new(1.0, 0.0, 1.5), GruenZone::Z3), // lateral distal
            (Vec3::new(0.0, 0.0, 0.5), GruenZone::Z4), // tip band
            (Vec3::new(-1.0, 0.0, 1.5), GruenZone::Z5), // medial distal
            (Vec3::new(-1.0, 0.0, 5.0), GruenZone::Z6), // medial middle
            (Vec3::new(-1.0, 0.0, 9.0), GruenZone::Z7), // medial proximal
        ];
        for (point, want) in cases {
            assert_eq!(gruen_zone(prox, tip, ml, point, 0.1), want);
        }
    }

    #[test]
    fn gruen_tip_band_wins_over_side() {
        let prox = Vec3::new(0.0, 0.0, 10.0);
        let tip = Vec3::new(0.0, 0.0, 0.0);
        let ml = Vec3::new(1.0, 0.0, 0.0);
        // Lateral point inside the 15 % tip band: zone 4, not 3.
        assert_eq!(
            gruen_zone(prox, tip, ml, Vec3::new(1.0, 0.0, 1.0), 0.15),
            GruenZone::Z4
        );
    }

    #[test]
    fn gruen_zone_keys_are_distinct() {
        let keys: Vec<_> = gruen_zones_in_order()
            .iter()
            .map(|z| z.key().to_string())
            .collect();
        assert_eq!(keys.len(), 7);
        let set: std::collections::BTreeSet<_> = keys.iter().collect();
        assert_eq!(set.len(), 7);
        assert_eq!(keys[0], "gruen-1");
        assert_eq!(keys[6], "gruen-7");
    }

    #[test]
    fn micromotion_thresholds_classify() {
        assert_eq!(MicromotionResult::classify(0.03), RiskLevel::Low);
        assert_eq!(MicromotionResult::classify(0.08), RiskLevel::Moderate);
        assert_eq!(MicromotionResult::classify(0.5), RiskLevel::High);
    }

    #[test]
    fn stiffer_foundation_reduces_micromotion() {
        let load = Force::from_n(2000.0);
        let soft = micromotion_analysis(&interface(), load, 0.3, &[800.0]);
        let stiff = micromotion_analysis(
            &InterfaceModel {
                foundation_stiffness: 20.0,
                ..interface()
            },
            load,
            0.3,
            &[800.0],
        );
        assert!(stiff.max_micromotion < soft.max_micromotion);
        assert!(soft.max_micromotion > 0.0);
        assert_eq!(soft.zone_micromotion.len(), 1);
    }

    #[test]
    fn micromotion_scales_with_load() {
        let small = micromotion_analysis(&interface(), Force::from_n(500.0), 0.3, &[800.0]);
        let large = micromotion_analysis(&interface(), Force::from_n(4000.0), 0.3, &[800.0]);
        assert!(large.max_micromotion > small.max_micromotion);
    }

    #[test]
    fn stress_shielding_index_bounded_and_flagged() {
        // Stem offloads proximal zones severely.
        let intact = [10.0, 8.0, 5.0, 3.0, 2.0, 2.0, 2.0];
        let implanted = [1.0, 2.0, 3.0, 2.8, 2.0, 2.0, 2.0];
        let result = stress_shielding_analysis(&intact, &implanted);
        assert!(result.shielding_index[0] > 0.85);
        assert!(result.shielding_index[6] < 0.05);
        assert!(result.has_resorption_risk());
        assert!(result.mean_index() > 0.0 && result.mean_index() < 1.0);
    }

    #[test]
    fn no_shielding_without_implant_effect() {
        let intact = [5.0, 5.0, 5.0];
        let result = stress_shielding_analysis(&intact, &intact);
        assert!(result.shielding_index.iter().all(|&s| s == 0.0));
        assert!(!result.has_resorption_risk());
    }

    #[test]
    fn micrometre_reporting() {
        let r = micromotion_analysis(&interface(), Force::from_n(1500.0), 0.25, &[800.0]);
        let um = micromotion_um(&r);
        assert!((um.to_mm() - r.max_micromotion * 1.0e3).abs() < 1e-9);
    }

    #[test]
    fn compliant_implant_softens_the_interface_in_series() {
        // Rigid limit: a very stiff implant returns the bare foundation.
        let stiff = CompliantImplant {
            foundation_stiffness: 5.0,
            implant_modulus_mpa: 110_000.0,
            implant_thickness_mm: 5.0,
        };
        let k_stiff = stiff.effective_stiffness();
        assert!((k_stiff - 5.0).abs() < 0.01, "rigid limit {k_stiff}");

        // A compliant coating (1 GPa, 2 mm) softens the interface.
        let coated = CompliantImplant {
            foundation_stiffness: 5.0,
            implant_modulus_mpa: 1_000.0,
            implant_thickness_mm: 2.0,
        };
        let k_coated = coated.effective_stiffness();
        assert!(k_coated < 5.0, "coating must soften: {k_coated}");
        // Series arithmetic hand-checked: 1/(1/5 + 2/1000) = 1/0.202 = 4.95.
        assert!((k_coated - 4.9505).abs() < 1e-3, "{k_coated}");
        // And the softened interface raises micromotion.
        let load = Force::from_n(2000.0);
        let rigid = micromotion_analysis(&stiff.interface_model(800.0, 0.5), load, 0.3, &[800.0]);
        let compliant =
            micromotion_analysis(&coated.interface_model(800.0, 0.5), load, 0.3, &[800.0]);
        assert!(
            compliant.max_micromotion > rigid.max_micromotion,
            "compliant implant must increase micromotion"
        );
    }

    #[test]
    fn gait_profile_is_a_double_hump() {
        let gait = GaitCycle::iso_double_hump(750.0);
        let samples = gait.samples(1000);
        assert_eq!(samples.len(), 1000);
        let peak = samples.iter().map(|(_, f)| f.to_n()).fold(0.0, f64::max);
        assert!((peak - 2.6 * 750.0).abs() < 1e-6, "peak {peak}");
        // Two distinct humps: local maxima in the stance window, separated
        // by a trough.
        let loads: Vec<f64> = samples.iter().map(|(_, f)| f.to_n()).collect();
        let humps_at = |frac: f64| loads[(frac * 1000.0) as usize];
        assert!(humps_at(0.15) > 2.0 * 750.0, "heel-strike hump");
        assert!(humps_at(0.45) > 2.0 * 750.0, "push-off hump");
        assert!(humps_at(0.30) < humps_at(0.15), "stance trough between");
        assert!(humps_at(0.85) < 0.3 * peak, "swing unloads");
        // Periodicity: loading at cycle end returns to unload.
        assert!(loads[999] < 0.3 * peak);
    }

    #[test]
    fn cyclic_micromotion_reports_peak_and_amplitude() {
        let gait = GaitCycle::iso_double_hump(750.0);
        let loads = gait.load_samples(200);
        let cyclic = micromotion_over_cycle(&interface(), &loads, 0.3, &[800.0, 400.0]);
        assert_eq!(cyclic.zone_amplitude.len(), 2);
        // The peak equals the static analysis at the peak load…
        let peak_load = loads.iter().copied().fold(0.0f64, |a, f| a.max(f.to_n()));
        let static_at_peak =
            micromotion_analysis(&interface(), Force::from_n(peak_load), 0.3, &[800.0, 400.0]);
        assert!(
            (cyclic.peak_micromotion - static_at_peak.max_micromotion).abs() < 1e-12,
            "{:.6} vs {:.6}",
            cyclic.peak_micromotion,
            static_at_peak.max_micromotion
        );
        // …and every zone has positive amplitude under cyclic loading…
        assert!(cyclic.zone_amplitude.iter().all(|a| *a > 0.0));
        assert!(cyclic.peak_amplitude >= cyclic.peak_micromotion * 0.1);
        // …and the risk class comes from the peak.
        assert_eq!(
            cyclic.risk,
            MicromotionResult::classify(cyclic.peak_micromotion)
        );
    }

    #[test]
    fn migration_closed_form_matches_numerical_integration() {
        let model = MigrationModel {
            rate: 1.0e-4,
            threshold: 0.05,
            bedding_in: 0.4,
        };
        let amplitude = 0.3;
        // Numeric integration of rate_at vs the closed form.
        let mut x = 0.0;
        let dt = 0.05f64;
        let cycles = 5000.0;
        let mut n = 0.0;
        while n < cycles {
            let h = dt.min(cycles - n);
            x += model.rate_at(amplitude, x) * h;
            n += h;
        }
        let closed = model.cumulative(amplitude, cycles);
        assert!(
            (x - closed).abs() < 0.02 * closed,
            "numeric {x} vs closed form {closed}"
        );
        // Below the threshold amplitude: no migration ever.
        assert_eq!(model.cumulative(0.05, 1.0e7), 0.0);
        assert_eq!(model.rate_at(0.04, 0.0), 0.0);
    }

    #[test]
    fn migration_velocity_decays_and_flags_risk() {
        let model = MigrationModel {
            rate: 1.0e-4,
            threshold: 0.05,
            bedding_in: 0.4,
        };
        let amplitude = 0.3;
        let cycles_per_year = 1.0e6;
        let v_new = model.velocity_per_year(amplitude, 0.0, cycles_per_year);
        let v_late = model.velocity_per_year(
            amplitude,
            model.cumulative(amplitude, 2.0e6),
            cycles_per_year,
        );
        assert!(v_new > 0.2, "fresh implant migrating fast: {v_new}");
        assert!(
            v_late < v_new * 0.5,
            "bedding-in must slow migration: {v_new} -> {v_late}"
        );
        // The at-risk flag follows the 0.2 mm/year boundary.
        assert!(model.is_at_risk(amplitude, 0.0, cycles_per_year));
        let stable = MigrationModel {
            threshold: 0.4, // amplitude below threshold → no migration at all
            ..model
        };
        assert!(!stable.is_at_risk(amplitude, 0.0, cycles_per_year));
        // Logarithmic growth: late doubling of the cycle count adds far
        // less than the first doubling did, and successive late increments
        // are nearly equal (the pure-log limit).
        let x1 = model.cumulative(amplitude, 1.0e6);
        let x2 = model.cumulative(amplitude, 2.0e6);
        let x4 = model.cumulative(amplitude, 4.0e6);
        assert!(
            x2 - x1 < 0.5 * x1,
            "late growth slower than early: {x1} {x2} {x4}"
        );
        let d1 = x2 - x1;
        let d2 = x4 - x2;
        assert!(
            (d2 - d1).abs() < 0.02 * d1,
            "log curve: nearly-equal late increments {d1} {d2}"
        );
    }
}
