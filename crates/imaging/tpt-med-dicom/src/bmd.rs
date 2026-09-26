//! BMD → apparent-density conversion (`rfcs/0007-bmd-apparent-density-conversion.md`).
//!
//! Most clinical QCT calibration phantoms (Mindways QCT Pro, CIRS/Image
//! Analysis, the European Forearm Phantom) report rod values as bone
//! mineral density (BMD) — mg/cm³, K₂HPO₄- or hydroxyapatite(CaHA)-
//! equivalent — not as the apparent (whole-tissue) density
//! [`crate::HounsfieldMapper::density_to_youngs_modulus`] expects. Converting
//! one to the other is a two-stage, literature-published, protocol- and
//! site-dependent relation: BMD → ash density, then ash density → apparent
//! density via the tissue's ash fraction. Neither stage ships a built-in
//! default in this crate — see the RFC's "No built-in default relation" for
//! why that is a deliberate, structural choice rather than an omission: both
//! [`BmdToAshDensity::new`] and [`AshFraction::new`] refuse to construct
//! without a non-empty citation, so a caller cannot end up with a
//! provenance-free apparent density from this path.

use crate::error::{DicomError, Result};
use tpt_med_units::Density;

/// Which mineral-equivalent convention a phantom's BMD value uses. Two rods
/// of the same physical composition read different numbers under the two
/// conventions, so this selects which published BMD→ash-density relation is
/// valid to apply — it is metadata carried alongside a fit, not itself a
/// unit conversion between the two conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BmdConvention {
    /// K₂HPO₄-equivalent (e.g. Mindways QCT Pro solid phantom).
    K2Hpo4Equivalent,
    /// Hydroxyapatite/CaHA-equivalent (e.g. CIRS, European Forearm Phantom).
    HydroxyapatiteEquivalent,
}

fn require_finite(value: f64, what: &str) -> Result<()> {
    if !value.is_finite() {
        return Err(DicomError::Calibration(format!("{what} must be finite")));
    }
    Ok(())
}

fn require_citation(source: &str) -> Result<()> {
    if source.trim().is_empty() {
        return Err(DicomError::Calibration(
            "a citation (source) is required and cannot be empty — see \
             rfcs/0007-bmd-apparent-density-conversion.md \"No built-in \
             default relation\""
                .into(),
        ));
    }
    Ok(())
}

/// A published (or site-fitted) linear relation from a phantom's BMD
/// (mg/cm³, in [`BmdConvention`]'s units) to ash density (g/cm³).
///
/// Cannot be constructed without a `source` citation.
#[derive(Debug, Clone, PartialEq)]
pub struct BmdToAshDensity {
    convention: BmdConvention,
    slope: f64,
    intercept: f64,
    source: String,
}

impl BmdToAshDensity {
    /// `slope`/`intercept` in `ash_density_g_cm3 = intercept + slope * bmd_mg_cm3`.
    /// `source` (paper/table/protocol citation) must be non-empty.
    pub fn new(
        convention: BmdConvention,
        slope: f64,
        intercept: f64,
        source: impl Into<String>,
    ) -> Result<Self> {
        require_finite(slope, "slope")?;
        require_finite(intercept, "intercept")?;
        let source = source.into();
        require_citation(&source)?;
        Ok(Self {
            convention,
            slope,
            intercept,
            source,
        })
    }

    /// The mineral-equivalent convention this relation was fitted under.
    pub fn convention(&self) -> BmdConvention {
        self.convention
    }

    /// The relation's citation.
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Evaluates the fitted line, returning an ash density in g/cm³.
    pub fn ash_density(&self, bmd_mg_cm3: f64) -> Density {
        Density::from_gcm3(self.intercept + self.slope * bmd_mg_cm3)
    }
}

/// A published (or site-fitted) ash fraction: the dry-mass fraction of bone
/// that is mineral, converting ash density to apparent (whole-tissue)
/// density: `apparent = ash / ash_fraction`.
///
/// Also requires a `source` citation, for the same reason as
/// [`BmdToAshDensity`].
#[derive(Debug, Clone, PartialEq)]
pub struct AshFraction {
    fraction: f64,
    source: String,
}

impl AshFraction {
    /// `fraction` must be in `(0.0, 1.0]` — an ash fraction is a mass
    /// fraction, so it cannot be zero, negative, or exceed 1.
    pub fn new(fraction: f64, source: impl Into<String>) -> Result<Self> {
        require_finite(fraction, "fraction")?;
        if !(fraction > 0.0 && fraction <= 1.0) {
            return Err(DicomError::Calibration(format!(
                "ash fraction must be in (0.0, 1.0], got {fraction}"
            )));
        }
        let source = source.into();
        require_citation(&source)?;
        Ok(Self { fraction, source })
    }

    /// The fraction value.
    pub fn fraction(&self) -> f64 {
        self.fraction
    }

    /// The relation's citation.
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Converts an ash density to apparent (whole-tissue) density.
    pub fn apparent_density(&self, ash_density: Density) -> Density {
        Density::from_gcm3(ash_density.value() / self.fraction)
    }
}

/// The composed two-stage conversion, carrying both citations together so a
/// caller (and, downstream, an audit trail) has one object naming the full
/// provenance of a BMD-derived apparent density.
#[derive(Debug, Clone, PartialEq)]
pub struct BmdToApparentDensity {
    bmd_to_ash: BmdToAshDensity,
    ash_fraction: AshFraction,
}

impl BmdToApparentDensity {
    /// Composes a BMD→ash relation with an ash→apparent fraction.
    pub fn new(bmd_to_ash: BmdToAshDensity, ash_fraction: AshFraction) -> Self {
        Self {
            bmd_to_ash,
            ash_fraction,
        }
    }

    /// The BMD-convention this conversion accepts.
    pub fn convention(&self) -> BmdConvention {
        self.bmd_to_ash.convention()
    }

    /// Converts a BMD value (mg/cm³, in `self.convention()`'s units) to
    /// apparent (whole-tissue) density in g/cm³.
    pub fn apparent_density(&self, bmd_mg_cm3: f64) -> Density {
        self.ash_fraction
            .apparent_density(self.bmd_to_ash.ash_density(bmd_mg_cm3))
    }

    /// Both citations, concatenated for a report/audit line.
    pub fn provenance(&self) -> String {
        format!(
            "BMD->ash: {}; ash->apparent: {}",
            self.bmd_to_ash.source(),
            self.ash_fraction.source()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs().max(1.0)
    }

    #[test]
    fn bmd_to_ash_density_evaluates_the_line() {
        let rel =
            BmdToAshDensity::new(BmdConvention::K2Hpo4Equivalent, 0.001, 0.05, "test fixture")
                .expect("constructs");
        assert!(close(
            rel.ash_density(200.0).value(),
            0.05 + 0.001 * 200.0,
            1e-12
        ));
    }

    #[test]
    fn bmd_to_ash_density_rejects_empty_source() {
        let err =
            BmdToAshDensity::new(BmdConvention::K2Hpo4Equivalent, 0.001, 0.05, "").unwrap_err();
        assert!(matches!(err, DicomError::Calibration(_)));
    }

    #[test]
    fn bmd_to_ash_density_rejects_whitespace_only_source() {
        let err =
            BmdToAshDensity::new(BmdConvention::K2Hpo4Equivalent, 0.001, 0.05, "   ").unwrap_err();
        assert!(matches!(err, DicomError::Calibration(_)));
    }

    #[test]
    fn bmd_to_ash_density_rejects_non_finite_coefficients() {
        assert!(
            BmdToAshDensity::new(BmdConvention::K2Hpo4Equivalent, f64::NAN, 0.05, "test").is_err()
        );
        assert!(BmdToAshDensity::new(
            BmdConvention::K2Hpo4Equivalent,
            0.001,
            f64::INFINITY,
            "test"
        )
        .is_err());
    }

    #[test]
    fn ash_fraction_converts_correctly() {
        let f = AshFraction::new(0.6, "test fixture").expect("constructs");
        let apparent = f.apparent_density(Density::from_gcm3(0.6));
        assert!(close(apparent.value(), 1.0, 1e-12));
    }

    #[test]
    fn ash_fraction_rejects_out_of_range() {
        assert!(AshFraction::new(0.0, "test").is_err());
        assert!(AshFraction::new(-0.1, "test").is_err());
        assert!(AshFraction::new(1.1, "test").is_err());
    }

    #[test]
    fn ash_fraction_rejects_empty_source() {
        let err = AshFraction::new(0.6, "").unwrap_err();
        assert!(matches!(err, DicomError::Calibration(_)));
    }

    #[test]
    fn ash_fraction_rejects_non_finite() {
        assert!(AshFraction::new(f64::NAN, "test").is_err());
    }

    #[test]
    fn composed_conversion_chains_both_stages() {
        let bmd_to_ash = BmdToAshDensity::new(
            BmdConvention::HydroxyapatiteEquivalent,
            0.0012,
            0.02,
            "paper A",
        )
        .expect("constructs");
        let ash_fraction = AshFraction::new(0.6, "paper B").expect("constructs");
        let conversion = BmdToApparentDensity::new(bmd_to_ash, ash_fraction);

        let bmd = 300.0;
        let expected_ash = 0.02 + 0.0012 * bmd;
        let expected_apparent = expected_ash / 0.6;

        assert_eq!(
            conversion.convention(),
            BmdConvention::HydroxyapatiteEquivalent
        );
        assert!(close(
            conversion.apparent_density(bmd).value(),
            expected_apparent,
            1e-12
        ));
        assert_eq!(
            conversion.provenance(),
            "BMD->ash: paper A; ash->apparent: paper B"
        );
    }
}
