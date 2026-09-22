//! Hounsfield Unit (HU) → material property mapping.
//!
//! CT numbers map to apparent bone density and then to elastic moduli via
//! the correlations used throughout the CT-based FEM literature (the exact
//! relations are specified in the tpt-medical design document; they follow
//! the Morgan–Keaveny style power laws).

use tpt_med_units::{Density, Modulus};

/// Which bone tissue a density/modulus correlation applies to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoneRegion {
    /// Cortical (compact) bone shell.
    Cortical,
    /// Trabecular (cancellous) bone core.
    Trabecular,
}

/// HU to material property conversions.
#[derive(Debug, Clone, Copy)]
pub struct HounsfieldMapper;

impl HounsfieldMapper {
    /// Converts HU to apparent density in g/cm³ via the linear CT
    /// approximation: water (0 HU) → 1.0, air (−1000 HU) → 0.0.
    ///
    /// This is the first-order screening relation; quantitative CT (QCT)
    /// calibration against a phantom replaces it in regulated pipelines.
    pub fn hu_to_density(hu: f64) -> Density {
        Density::from_gcm3((hu + 1000.0) / 1000.0)
    }

    /// Apparent density (g/cm³) → Young's modulus (MPa).
    ///
    /// - Cortical: `E = 10500 · ρ^2.0`
    /// - Trabecular: `E = 6850 · ρ^1.49`
    ///
    /// Negative/zero densities clamp to zero modulus (air voxels).
    pub fn density_to_youngs_modulus(density: Density, region: BoneRegion) -> Modulus {
        let rho = density.value().max(0.0);
        let e = match region {
            BoneRegion::Cortical => 10_500.0 * rho.powf(2.0),
            BoneRegion::Trabecular => 6_850.0 * rho.powf(1.49),
        };
        Modulus::from_mpa(e)
    }

    /// Convenience: HU → Young's modulus for a region.
    pub fn hu_to_youngs_modulus(hu: f64, region: BoneRegion) -> Modulus {
        Self::density_to_youngs_modulus(Self::hu_to_density(hu), region)
    }

    /// Default HU threshold separating bone from soft tissue. Values in the
    /// literature range 130–300 HU; 200 HU is a common middle choice for
    /// appendicular CT.
    pub const DEFAULT_BONE_THRESHOLD_HU: f64 = 200.0;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs().max(1.0)
    }

    #[test]
    fn hu_density_landmarks() {
        assert!(close(
            HounsfieldMapper::hu_to_density(-1000.0).value(),
            0.0,
            1e-9
        ));
        assert!(close(
            HounsfieldMapper::hu_to_density(0.0).value(),
            1.0,
            1e-9
        ));
        assert!(close(
            HounsfieldMapper::hu_to_density(1000.0).value(),
            2.0,
            1e-9
        ));
    }

    #[test]
    fn cortical_modulus_law() {
        // E = 10500 · ρ²  at ρ = 2 g/cm³ → 42 GPa (upper cortical range)
        let e = HounsfieldMapper::density_to_youngs_modulus(
            Density::from_gcm3(2.0),
            BoneRegion::Cortical,
        );
        assert!(close(e.to_mpa(), 42_000.0, 1e-9));
    }

    #[test]
    fn trabecular_modulus_law() {
        // E = 6850 · ρ^1.49 at ρ = 0.5 g/cm³
        let e = HounsfieldMapper::density_to_youngs_modulus(
            Density::from_gcm3(0.5),
            BoneRegion::Trabecular,
        );
        assert!(close(e.to_mpa(), 6850.0 * 0.5f64.powf(1.49), 1e-9));
    }

    #[test]
    fn air_clamps_to_zero() {
        let e = HounsfieldMapper::density_to_youngs_modulus(
            Density::from_gcm3(-0.5),
            BoneRegion::Cortical,
        );
        assert_eq!(e.to_mpa(), 0.0);
    }

    #[test]
    fn hu_to_modulus_end_to_end() {
        // 700 HU → ρ = 1.7 → cortical E = 10500 · 1.7² = 30345
        let e = HounsfieldMapper::hu_to_youngs_modulus(700.0, BoneRegion::Cortical);
        assert!(close(e.to_mpa(), 10500.0 * 1.7f64.powi(2), 1e-9));
    }
}
