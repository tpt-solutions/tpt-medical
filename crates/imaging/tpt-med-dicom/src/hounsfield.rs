//! Hounsfield Unit (HU) → material property mapping.
//!
//! CT numbers map to apparent bone density and then to elastic moduli via
//! the correlations used throughout the CT-based FEM literature (the exact
//! relations are specified in the tpt-medical design document; they follow
//! the Morgan–Keaveny style power laws).

use crate::error::{DicomError, Result};
use tpt_med_units::{Density, Modulus};

/// A linear HU → density calibration, fitted from calibration-phantom
/// measurements rather than assumed.
///
/// [`HounsfieldMapper::hu_to_density`] uses a fixed two-point line — water
/// (0 HU) is 1.0 g/cm³, air (−1000 HU) is 0.0 — which is a screening
/// estimate, not a measurement: it ignores the scanner's actual calibration
/// (kVp, reconstruction kernel, day-to-day drift). Quantitative CT (QCT)
/// replaces that assumption with a line fitted through *measured* points: the
/// HU each rod of a calibration phantom reads in this scan, paired with that
/// rod's manufacturer-specified value, fitted by ordinary least squares.
///
/// # Units are the caller's responsibility
///
/// This type does not know or assume what physical quantity its fitted
/// points are in — it is a linear regression, nothing more. That matters
/// because the two conventions calibration phantoms actually use are
/// physically different quantities:
///
/// - **Apparent (whole-tissue) density**, g/cm³ — what
///   [`HounsfieldMapper::density_to_youngs_modulus`] expects. A phantom
///   rod's value is directly usable here only if the manufacturer specifies
///   it as equivalent apparent density.
/// - **Bone mineral density (BMD)**, typically mg/cm³ K₂HPO₄- or
///   CaHA-equivalent — what most clinical QCT phantoms (Mindways QCT Pro,
///   CIRS/Image Analysis) actually report. BMD is mineral concentration, not
///   whole-tissue density, and converting it to apparent density needs a
///   documented, protocol-specific relation this crate does not supply
///   (the literature has more than one, and picking one silently would be
///   exactly the "confidently wrong" number this crate elsewhere refuses to
///   produce). If your phantom reports BMD, convert to apparent density
///   *before* calling [`Self::fit`], using whatever relation your V&V
///   protocol documents — do not feed BMD values into
///   [`Self::hu_to_apparent_density`] and treat the result as density.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QctCalibration {
    slope: f64,
    intercept: f64,
}

impl QctCalibration {
    /// The fixed two-point line [`HounsfieldMapper::hu_to_density`] uses:
    /// water (0 HU) → 1.0 g/cm³, air (−1000 HU) → 0.0 g/cm³. Exposed so a
    /// caller can compare a real phantom fit against this screening default,
    /// or use the same [`Self::hu_to_apparent_density`] path for both.
    pub fn screening_default() -> Self {
        Self {
            slope: 0.001,
            intercept: 1.0,
        }
    }

    /// Fits a calibration line through `points` (`(HU, known_value)` pairs,
    /// one per calibration-phantom rod) by ordinary least squares.
    ///
    /// Requires at least two points with a non-degenerate spread of HU
    /// values — a real phantom has multiple rods at different densities
    /// precisely so this fit is over-determined, not just two points drawn
    /// through noise.
    pub fn fit(points: &[(f64, f64)]) -> Result<Self> {
        if points.len() < 2 {
            return Err(DicomError::Calibration(format!(
                "at least 2 calibration points are required, got {}",
                points.len()
            )));
        }
        if points
            .iter()
            .any(|(hu, v)| !hu.is_finite() || !v.is_finite())
        {
            return Err(DicomError::Calibration(
                "calibration points must be finite (no NaN/infinite HU or value)".into(),
            ));
        }

        let n = points.len() as f64;
        let mean_hu = points.iter().map(|(hu, _)| hu).sum::<f64>() / n;
        let mean_v = points.iter().map(|(_, v)| v).sum::<f64>() / n;

        let mut cov = 0.0;
        let mut var_hu = 0.0;
        for &(hu, v) in points {
            let dhu = hu - mean_hu;
            cov += dhu * (v - mean_v);
            var_hu += dhu * dhu;
        }

        // A near-zero variance means every point sits at (nearly) the same
        // HU — physically, every rod read the same number, which cannot
        // happen on a real phantom and cannot determine a slope. Reject
        // rather than divide by a near-zero denominator into a wild slope.
        const MIN_HU_VARIANCE: f64 = 1e-6;
        if var_hu < MIN_HU_VARIANCE {
            return Err(DicomError::Calibration(
                "calibration points have no meaningful spread of HU values \
                 (all rods read approximately the same HU)"
                    .into(),
            ));
        }

        let slope = cov / var_hu;
        let intercept = mean_v - slope * mean_hu;
        Ok(Self { slope, intercept })
    }

    /// The fitted slope, in calibration-value units per HU.
    pub fn slope(&self) -> f64 {
        self.slope
    }

    /// The fitted intercept: the calibrated value at 0 HU.
    pub fn intercept(&self) -> f64 {
        self.intercept
    }

    /// Evaluates the fitted line at `hu`, in whatever unit the calibration
    /// points were given in.
    pub fn evaluate(&self, hu: f64) -> f64 {
        self.intercept + self.slope * hu
    }

    /// Evaluates the fitted line as an apparent density in g/cm³.
    ///
    /// Only meaningful if this calibration was fitted from points already in
    /// g/cm³ apparent density — see the type-level docs on the BMD-vs-density
    /// distinction before calling this on a clinical QCT phantom's raw rod
    /// values.
    pub fn hu_to_apparent_density(&self, hu: f64) -> Density {
        Density::from_gcm3(self.evaluate(hu))
    }
}

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
    /// This is the first-order screening relation — exactly
    /// [`QctCalibration::screening_default`] — not a measurement; a real
    /// [`QctCalibration`] fitted from a calibration phantom's rods replaces
    /// it in regulated pipelines. See that type's docs for the units caveat
    /// before treating a calibrated result as interchangeable with this one.
    pub fn hu_to_density(hu: f64) -> Density {
        QctCalibration::screening_default().hu_to_apparent_density(hu)
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

    /// Convenience: HU → Young's modulus for a region, via a fitted
    /// [`QctCalibration`] rather than the fixed screening line. The
    /// calibration's units caveat applies here exactly as it does to
    /// [`QctCalibration::hu_to_apparent_density`].
    pub fn hu_to_youngs_modulus_calibrated(
        hu: f64,
        calibration: &QctCalibration,
        region: BoneRegion,
    ) -> Modulus {
        Self::density_to_youngs_modulus(calibration.hu_to_apparent_density(hu), region)
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

    #[test]
    fn screening_default_matches_the_fixed_hu_to_density_line() {
        // hu_to_density is now defined in terms of screening_default(); this
        // pins that the refactor didn't change its behaviour.
        let cal = QctCalibration::screening_default();
        for hu in [-1000.0, -500.0, 0.0, 500.0, 1000.0, 2500.0] {
            assert!(close(
                cal.hu_to_apparent_density(hu).value(),
                HounsfieldMapper::hu_to_density(hu).value(),
                1e-12
            ));
        }
    }

    #[test]
    fn fit_recovers_an_exact_line() {
        // Three points on y = 0.0008x + 1.05 exactly -- a noiseless
        // "phantom" the least-squares fit should reproduce exactly.
        let slope = 0.0008;
        let intercept = 1.05;
        let points: Vec<(f64, f64)> = [-200.0, 100.0, 800.0]
            .iter()
            .map(|&hu| (hu, intercept + slope * hu))
            .collect();
        let cal = QctCalibration::fit(&points).expect("fits");
        assert!(close(cal.slope(), slope, 1e-9));
        assert!(close(cal.intercept(), intercept, 1e-9));
    }

    #[test]
    fn fit_minimises_squared_error_on_noisy_points() {
        // Four rods roughly on y = 0.001x + 1.0 with symmetric +/- noise;
        // OLS should recover very close to the true line since the noise
        // cancels by construction.
        let points = [
            (-200.0, 1.0 + 0.001 * -200.0 - 0.01),
            (-200.0, 1.0 + 0.001 * -200.0 + 0.01),
            (600.0, 1.0 + 0.001 * 600.0 - 0.01),
            (600.0, 1.0 + 0.001 * 600.0 + 0.01),
        ];
        let cal = QctCalibration::fit(&points).expect("fits");
        assert!(close(cal.slope(), 0.001, 1e-6));
        assert!(close(cal.intercept(), 1.0, 1e-6));
    }

    #[test]
    fn fit_rejects_too_few_points() {
        let err = QctCalibration::fit(&[(0.0, 1.0)]).unwrap_err();
        assert!(matches!(err, DicomError::Calibration(_)));
    }

    #[test]
    fn fit_rejects_empty_points() {
        let err = QctCalibration::fit(&[]).unwrap_err();
        assert!(matches!(err, DicomError::Calibration(_)));
    }

    #[test]
    fn fit_rejects_degenerate_hu_spread() {
        // Every "rod" reads the same HU -- no slope is determinable.
        let err = QctCalibration::fit(&[(100.0, 1.0), (100.0, 1.5), (100.0, 2.0)]).unwrap_err();
        assert!(matches!(err, DicomError::Calibration(_)));
    }

    #[test]
    fn fit_rejects_non_finite_points() {
        assert!(QctCalibration::fit(&[(0.0, 1.0), (f64::NAN, 2.0)]).is_err());
        assert!(QctCalibration::fit(&[(0.0, 1.0), (f64::INFINITY, 2.0)]).is_err());
    }

    #[test]
    fn calibrated_modulus_uses_the_fitted_density() {
        let cal = QctCalibration::fit(&[(-200.0, 0.7), (800.0, 1.9)]).expect("fits");
        let expected_density = cal.hu_to_apparent_density(300.0);
        let expected =
            HounsfieldMapper::density_to_youngs_modulus(expected_density, BoneRegion::Cortical);
        let got =
            HounsfieldMapper::hu_to_youngs_modulus_calibrated(300.0, &cal, BoneRegion::Cortical);
        assert!(close(got.to_mpa(), expected.to_mpa(), 1e-12));
    }
}
