//! Linear elastic bone mechanics, HU-based property assignment, and
//! Wolff's-law density remodeling.
//!
//! Modulus correlations live in `tpt-med-dicom` (`HounsfieldMapper`); this
//! crate adds the structural material descriptions, anisotropy
//! classification, and the time-domain remodeling law used for stress
//! shielding assessment (see `tpt-med-orthopedics`).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use tpt_med_core::BoneType;
use tpt_med_dicom::{BoneRegion, HounsfieldMapper};
use tpt_med_geometry::Vec3;
use tpt_med_units::{Density, Modulus};

/// Bone tissue classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TissueClass {
    /// Compact cortical bone.
    Cortical,
    /// Spongy trabecular bone.
    Trabecular,
}

/// Degree and direction of elastic anisotropy.
#[derive(Debug, Clone, PartialEq)]
pub enum Anisotropy {
    /// Isotropic (single E, ν).
    Isotropic,
    /// Transversely isotropic about an axis (e.g. along the femoral shaft).
    TransverselyIsotropic {
        /// Axis of symmetry.
        axis: Vec3,
        /// Longitudinal modulus (MPa).
        e_long: f64,
        /// Transverse modulus (MPa).
        e_trans: f64,
    },
    /// Orthotropic with three moduli along three axes.
    Orthotropic {
        /// Principal axes.
        axes: [Vec3; 3],
        /// Moduli along each axis (MPa).
        e: [f64; 3],
    },
}

/// A bone material description.
#[derive(Debug, Clone, PartialEq)]
pub struct BoneMaterial {
    /// Which bone this describes.
    pub bone_type: BoneType,
    /// Tissue class the properties refer to.
    pub tissue_class: TissueClass,
    /// Young's modulus (MPa).
    pub youngs_modulus: f64,
    /// Poisson's ratio (≈0.3 for bone).
    pub poissons_ratio: f64,
    /// Yield stress (MPa), 0 when unknown.
    pub yield_stress: f64,
    /// Ultimate stress (MPa), 0 when unknown.
    pub ultimate_stress: f64,
    /// Anisotropy classification.
    pub anisotropy: Anisotropy,
}

/// A density → modulus relation, so a study can supply a calibrated law
/// (e.g. fit from a QCT calibration phantom via
/// `tpt-med-dicom::QctCalibration`) instead of the default power law.
///
/// Implementors receive the apparent density and the bone region the
/// meshing pipeline classified the voxel into; the default [`PowerLaw`]
/// reproduces the Morgan–Keaveny-style correlations.
pub trait ModulusLaw {
    /// Modulus (MPa) at the given apparent density.
    fn modulus(&self, density: Density, region: BoneRegion) -> Modulus;
}

/// The default HU-correlation power law (cortical `10500 ρ²`,
/// trabecular `6850 ρ^1.49`, MPa), as used by
/// [`BoneMaterial::from_hu`].
#[derive(Debug, Clone, Copy, Default)]
pub struct PowerLaw;

impl ModulusLaw for PowerLaw {
    fn modulus(&self, density: Density, region: BoneRegion) -> Modulus {
        HounsfieldMapper::density_to_youngs_modulus(density, region)
    }
}

impl<F> ModulusLaw for F
where
    F: Fn(Density, BoneRegion) -> Modulus,
{
    fn modulus(&self, density: Density, region: BoneRegion) -> Modulus {
        self(density, region)
    }
}

impl BoneMaterial {
    /// Reference cortical femur properties (literature screening values:
    /// E ≈ 17 GPa, ν = 0.3, σy ≈ 110 MPa, σu ≈ 130 MPa).
    pub fn cortical_reference(bone_type: BoneType) -> Self {
        Self {
            bone_type,
            tissue_class: TissueClass::Cortical,
            youngs_modulus: 17_000.0,
            poissons_ratio: 0.3,
            yield_stress: 110.0,
            ultimate_stress: 130.0,
            anisotropy: Anisotropy::Isotropic,
        }
    }

    /// Reference trabecular femur properties (E ≈ 0.7 GPa at 0.3 g/cm³
    /// apparent density).
    pub fn trabecular_reference(bone_type: BoneType) -> Self {
        Self {
            bone_type,
            tissue_class: TissueClass::Trabecular,
            youngs_modulus: 700.0,
            poissons_ratio: 0.3,
            yield_stress: 10.0,
            ultimate_stress: 15.0,
            anisotropy: Anisotropy::Isotropic,
        }
    }

    /// Builds properties from a CT Hounsfield value using the
    /// [`HounsfieldMapper`] correlations and a density-based tissue split.
    pub fn from_hu(hu: f64, bone_type: BoneType, poissons_ratio: f64) -> Self {
        Self::from_hu_with_law(hu, bone_type, poissons_ratio, &PowerLaw)
    }

    /// Like [`Self::from_hu`], but evaluates the caller's density → modulus
    /// law instead of the default power law — the calibration hook for
    /// studies with a phantom-fitted relation.
    pub fn from_hu_with_law(
        hu: f64,
        bone_type: BoneType,
        poissons_ratio: f64,
        law: &impl ModulusLaw,
    ) -> Self {
        let density = HounsfieldMapper::hu_to_density(hu);
        let region = if density.value() >= 1.3 {
            BoneRegion::Cortical
        } else {
            BoneRegion::Trabecular
        };
        let e = law.modulus(density, region).to_mpa();
        Self {
            bone_type,
            tissue_class: match region {
                BoneRegion::Cortical => TissueClass::Cortical,
                BoneRegion::Trabecular => TissueClass::Trabecular,
            },
            youngs_modulus: e,
            poissons_ratio,
            yield_stress: 0.0,
            ultimate_stress: 0.0,
            anisotropy: Anisotropy::Isotropic,
        }
    }

    /// Transverse modulus for transversely isotropic classification; falls
    /// back to `youngs_modulus` for isotropic materials.
    pub fn effective_modulus(&self, direction: Vec3) -> f64 {
        match &self.anisotropy {
            Anisotropy::Isotropic => self.youngs_modulus,
            Anisotropy::TransverselyIsotropic {
                axis,
                e_long,
                e_trans,
            } => {
                let a = axis.normalize();
                let d = direction.normalize();
                let along = a.dot(d).powi(2);
                e_long * along + e_trans * (1.0 - along)
            }
            Anisotropy::Orthotropic { axes, e } => {
                let d = direction.normalize();
                axes.iter()
                    .zip(e)
                    .map(|(a, &modulus)| modulus * a.dot(d).powi(2))
                    .sum()
            }
        }
    }
}

/// Mechanical stimulus measures driving remodeling (Wolff's law).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RemodelingStimulus {
    /// Strain energy density (mJ/mm³ = MPa).
    StrainEnergyDensity,
    /// Maximum principal strain (dimensionless, microstrain/1000).
    PrincipalStrain,
    /// Damage accumulation rate.
    DamageAccumulation,
}

/// Wolff's-law remodeling: density adapts toward the mechanical stimulus
/// with a dead zone ("lazy zone") around the reference stimulus.
#[derive(Debug, Clone, Copy)]
pub struct BoneRemodelingModel {
    /// Stimulus measure used.
    pub stimulus: RemodelingStimulus,
    /// Maximum apposition rate (g/cm³ per day) above the lazy zone.
    pub apposition_rate: f64,
    /// Maximum resorption rate (g/cm³ per day) below the lazy zone.
    pub resorption_rate: f64,
    /// Lazy zone bounds `(lower, upper)` in stimulus units; no remodeling
    /// happens inside.
    pub lazy_zone: (f64, f64),
    /// Physiologically viable density bounds (g/cm³).
    pub viable_density: (f64, f64),
}

impl Default for BoneRemodelingModel {
    fn default() -> Self {
        // Frost-style mechanostat screening parameters: reference SED
        // stimulus ~0.004 mJ/mm³ with a ±35% lazy zone.
        Self {
            stimulus: RemodelingStimulus::StrainEnergyDensity,
            apposition_rate: 0.003,
            resorption_rate: 0.002,
            lazy_zone: (0.0026, 0.0054),
            viable_density: (0.02, 2.0),
        }
    }
}

impl BoneRemodelingModel {
    /// Spatial remodeling: drives a per-voxel density field from a solved
    /// stimulus field (e.g. strain energy density per element from a
    /// `tpt-med-biomechanics` result), advancing `dt_days` in lockstep.
    /// Input slices are parallel (`densities[i]` pairs with `stimuli[i]`).
    /// Viable-range clamping is applied per voxel.
    pub fn remodel_field(
        &self,
        densities: &[Density],
        stimuli: &[f64],
        dt_days: f64,
        viable: (f64, f64),
    ) -> Vec<Density> {
        assert_eq!(
            densities.len(),
            stimuli.len(),
            "density and stimulus fields must pair"
        );
        densities
            .iter()
            .zip(stimuli)
            .map(|(&rho, &stimulus)| {
                crate::clamp_viable(self.update_density(rho, stimulus, dt_days), viable)
            })
            .collect()
    }

    /// One-day density update for a given stimulus level (g/cm³).
    ///
    /// Response scales linearly with normalized over-/under-stimulus:
    /// `Δρ = rate · (S/S_ref − 1)` where `S_ref` is the nearer lazy-zone
    /// bound, clamped to the viable range.
    pub fn update_density(&self, current_density: Density, stimulus: f64, dt_days: f64) -> Density {
        let rho = current_density.value();
        let (lower, upper) = self.lazy_zone;
        let delta = if stimulus > upper {
            let over = stimulus / upper - 1.0;
            self.apposition_rate * over.min(3.0) * dt_days
        } else if stimulus < lower {
            let under = 1.0 - stimulus / lower;
            -self.resorption_rate * under.min(3.0) * dt_days
        } else {
            0.0
        };
        Density::from_gcm3(rho + delta)
    }

    /// Equilibrium prediction: repeatedly apply the update for `days`,
    /// returning the final density (screening helper).
    pub fn simulate_days(&self, start: Density, stimulus: f64, days: u32) -> Density {
        let mut rho = start;
        for _ in 0..days {
            rho = self.update_density(rho, stimulus, 1.0);
        }
        rho
    }
}

/// Viable-density clamp used by remodeling pipelines.
pub fn clamp_viable(density: Density, viable: (f64, f64)) -> Density {
    Density::from_gcm3(density.value().clamp(viable.0, viable.1))
}

/// Reference modulus accessor kept symbolically exported for the
/// regulatory layer's material tables.
pub fn reference_modulus(tissue: TissueClass) -> Modulus {
    match tissue {
        TissueClass::Cortical => Modulus::from_mpa(17_000.0),
        TissueClass::Trabecular => Modulus::from_mpa(700.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calibrated_law_replaces_power_law() {
        // A phantom-fitted law: 20 GPa at 1.7 g/cm3 regardless of region.
        let calibrated = |density: Density, _region: BoneRegion| {
            tpt_med_units::Modulus::from_mpa(20_000.0 * density.value() / 1.7)
        };
        let m = BoneMaterial::from_hu_with_law(700.0, BoneType::Femur, 0.3, &calibrated);
        assert!((m.youngs_modulus - 20_000.0).abs() < 1e-9);
        // Default law gives 10500 * 1.7^2 = 30345 at 700 HU.
        let default = BoneMaterial::from_hu(700.0, BoneType::Femur, 0.3);
        assert!((default.youngs_modulus - 30_345.0).abs() < 1e-6);
        // PowerLaw struct reproduces the default exactly.
        let pl = BoneMaterial::from_hu_with_law(700.0, BoneType::Femur, 0.3, &PowerLaw);
        assert_eq!(pl.youngs_modulus, default.youngs_modulus);
    }

    #[test]
    fn remodel_field_drives_each_voxel_independently() {
        let model = BoneRemodelingModel::default();
        let densities = [
            Density::from_gcm3(1.2),
            Density::from_gcm3(1.2),
            Density::from_gcm3(1.2),
        ];
        let stimuli = [0.02, 0.004, 0.0005]; // over / lazy / under
        let out = model.remodel_field(&densities, &stimuli, 30.0, (0.02, 2.0));
        assert_eq!(out.len(), 3);
        assert!(out[0].value() > 1.2, "over-stimulated gains density");
        assert_eq!(out[1].value(), 1.2, "lazy zone unchanged");
        assert!(out[2].value() < 1.2, "under-stimulated loses density");
        // Clamped at the viable ceiling.
        let hot = model.remodel_field(&[Density::from_gcm3(1.9)], &[1.0], 400.0, (0.02, 2.0));
        assert_eq!(hot[0].value(), 2.0);
    }

    #[test]
    fn hu_assignment_matches_dicom_correlations() {
        let m = BoneMaterial::from_hu(700.0, BoneType::Femur, 0.3);
        assert_eq!(m.tissue_class, TissueClass::Cortical);
        assert!((m.youngs_modulus - 10_500.0 * 1.7f64.powi(2)).abs() < 1e-6);

        let t = BoneMaterial::from_hu(250.0, BoneType::Vertebra { level: "L4".into() }, 0.3);
        assert_eq!(t.tissue_class, TissueClass::Trabecular);
        assert!((t.youngs_modulus - 6850.0 * 1.25f64.powf(1.49)).abs() < 1e-6);
    }

    #[test]
    fn remodeling_respects_lazy_zone() {
        let model = BoneRemodelingModel::default();
        let rho = Density::from_gcm3(1.2);
        // Inside the lazy zone: no change.
        let inside = model.update_density(rho, 0.004, 30.0);
        assert_eq!(inside.value(), 1.2);
        // Over-stimulated: gains density.
        let over = model.update_density(rho, 0.02, 30.0);
        assert!(over.value() > rho.value());
        // Under-stimulated: loses density.
        let under = model.update_density(rho, 0.0005, 30.0);
        assert!(under.value() < rho.value());
    }

    #[test]
    fn stress_shielding_scenario_remodels_down() {
        // A stiff implant shields bone (low stimulus) → density loss over
        // months; apposition saturates the other direction.
        let model = BoneRemodelingModel::default();
        let shielded = model.simulate_days(Density::from_gcm3(1.2), 0.0001, 180);
        assert!(shielded.value() < 1.2);
        assert!(shielded.value() >= 1.2 - model.resorption_rate * 3.0 * 180.0);
        let over = model.simulate_days(Density::from_gcm3(1.2), 0.05, 180);
        assert!(over.value() > 1.2);
    }

    #[test]
    fn anisotropic_effective_modulus_interpolates() {
        let m = BoneMaterial {
            anisotropy: Anisotropy::TransverselyIsotropic {
                axis: Vec3::Z,
                e_long: 17_000.0,
                e_trans: 11_000.0,
            },
            ..BoneMaterial::cortical_reference(BoneType::Femur)
        };
        assert!((m.effective_modulus(Vec3::Z) - 17_000.0).abs() < 1e-6);
        assert!((m.effective_modulus(Vec3::X) - 11_000.0).abs() < 1e-6);
        let mid = m.effective_modulus(Vec3::new(1.0, 0.0, 1.0));
        assert!((11_000.0..17_000.0).contains(&mid));
    }

    #[test]
    fn viable_clamp_bounds() {
        let rho = clamp_viable(Density::from_gcm3(5.0), (0.02, 2.0));
        assert_eq!(rho.value(), 2.0);
    }
}
