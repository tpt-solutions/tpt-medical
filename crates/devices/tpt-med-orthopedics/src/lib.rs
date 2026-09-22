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
}
