//! Grid-convergence reporting as first-class
//! `CalculationVerification` evidence (ASME V&V 40).
//!
//! A refinement study solves the same question of interest at several mesh
//! densities and reports the observed convergence order against the
//! expected asymptotic order of the discretisation (2nd order for the Q1
//! hex core in the energy norm and displacement quantities in the
//! pre-asymptotic-to-asymptotic range).

use crate::SolverError;

/// One refinement level.
#[derive(Debug, Clone, Copy)]
pub struct Level {
    /// Characteristic mesh size `h` (any consistent measure — element edge
    /// length, DOF count^(-1/3), …).
    pub h: f64,
    /// The quantity of interest computed at this level.
    pub value: f64,
}

/// Result of a refinement study.
#[derive(Debug, Clone)]
pub struct ConvergenceReport {
    /// Input levels, in submission order.
    pub levels: Vec<Level>,
    /// Reference value the errors are measured against (an analytic value,
    /// or the finest-mesh result when none exists — recorded as such).
    pub reference: f64,
    /// True when `reference` is an analytic value rather than the
    /// finest-mesh extrapolation.
    pub reference_is_analytic: bool,
    /// Relative errors |value − reference| / |reference| per level.
    pub relative_errors: Vec<f64>,
    /// Observed convergence order between successive levels,
    /// `p = ln(eᵢ/eᵢ₊₁) / ln(hᵢ₊₁/hᵢ)`; empty when fewer than two usable
    /// error pairs exist.
    pub observed_order: Vec<f64>,
}

/// Runs a refinement study from `(h, value)` pairs.
///
/// `expected_order` is only documentation on the report's consumer side
/// (the caller compares `observed_order` against it); the function computes
/// the measured order and never fails on non-converging data — it reports
/// what the numbers say, including negative orders for oscillating series.
pub fn convergence_study(
    mut levels: Vec<Level>,
    reference: Option<f64>,
) -> Result<ConvergenceReport, SolverError> {
    if levels.len() < 2 {
        return Err(SolverError::Invalid(
            "a refinement study needs at least two levels".into(),
        ));
    }
    levels.sort_by(|a, b| a.h.partial_cmp(&b.h).unwrap_or(core::cmp::Ordering::Equal));
    for l in &levels {
        if !(l.h.is_finite() && l.h > 0.0) || !l.value.is_finite() {
            return Err(SolverError::Invalid(
                "levels must carry finite positive h and finite values".into(),
            ));
        }
    }

    let (reference, reference_is_analytic) = match reference {
        Some(v) if v.is_finite() => (v, true),
        Some(_) => {
            return Err(SolverError::Invalid(
                "reference value must be finite when given".into(),
            ))
        }
        // Sorted ascending by h: the FIRST level is the finest mesh.
        None => (levels.first().expect("non-empty").value, false),
    };

    let relative_errors: Vec<f64> = levels
        .iter()
        .map(|l| (l.value - reference).abs() / reference.abs().max(1e-300))
        .collect();

    // Levels are sorted ascending by h, so consecutive pairs run fine ->
    // coarse: the observed order is
    // p = ln(e_coarse / e_fine) / ln(h_coarse / h_fine), positive when the
    // error shrinks with refinement.
    let mut observed_order = Vec::new();
    for i in 1..relative_errors.len() {
        let (e_fine, e_coarse) = (relative_errors[i - 1], relative_errors[i]);
        let h_ratio = levels[i].h / levels[i - 1].h;
        if e_fine > 0.0 && e_coarse > 0.0 && h_ratio > 1.0 {
            observed_order.push((e_coarse / e_fine).ln() / h_ratio.ln());
        } else {
            observed_order.push(f64::NAN);
        }
    }

    Ok(ConvergenceReport {
        levels,
        reference,
        reference_is_analytic,
        relative_errors,
        observed_order,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn second_order_series_reports_order_two() {
        // h = 1, 1/2, 1/4 with 2nd-order errors; levels are reported sorted
        // ascending by h (finest first), errors [0.005, 0.02, 0.08].
        let levels = vec![
            Level {
                h: 1.0,
                value: 1.0 + 0.08,
            },
            Level {
                h: 0.5,
                value: 1.0 + 0.02,
            },
            Level {
                h: 0.25,
                value: 1.0 + 0.005,
            },
        ];
        let report = convergence_study(levels, Some(1.0)).expect("study");
        assert!(report.reference_is_analytic);
        assert!((report.relative_errors[0] - 0.005).abs() < 1e-12);
        assert!((report.relative_errors[2] - 0.08).abs() < 1e-12);
        for p in &report.observed_order {
            assert!((p - 2.0).abs() < 1e-9, "order {p}");
        }
    }

    #[test]
    fn finest_mesh_reference_is_recorded_as_non_analytic() {
        let levels = vec![
            Level {
                h: 1.0,
                value: 1.10,
            },
            Level {
                h: 0.5,
                value: 1.02,
            },
            Level {
                h: 0.25,
                value: 1.00,
            },
        ];
        let report = convergence_study(levels, None).expect("study");
        assert!(!report.reference_is_analytic);
        assert_eq!(
            report.reference, 1.00,
            "finest (smallest h) is the reference"
        );
        // finest level error is zero by construction
        assert_eq!(*report.relative_errors.first().expect("levels"), 0.0);
    }

    #[test]
    fn rejects_single_level_and_bad_input() {
        assert!(convergence_study(vec![Level { h: 1.0, value: 1.0 }], None).is_err());
        assert!(convergence_study(
            vec![
                Level { h: 1.0, value: 1.0 },
                Level {
                    h: -1.0,
                    value: 1.0
                }
            ],
            None
        )
        .is_err());
    }
}
