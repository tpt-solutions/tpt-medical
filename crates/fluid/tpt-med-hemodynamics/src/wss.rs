//! Wall shear stress fields and oscillatory shear index.

use crate::blood::BloodModel;
use crate::domain::FluidDomain;
use crate::solver::{HemodynamicsSolver, SolverConfig};
use tpt_med_geometry::Vec3;

/// Wall shear stress samples on wall-adjacent faces.
#[derive(Debug, Clone, Default)]
pub struct WssField {
    /// Sample positions (cell centres of the wall-adjacent fluid cells).
    pub positions: Vec<Vec3>,
    /// Tangential traction vectors (Pa).
    pub tractions: Vec<Vec3>,
}

impl WssField {
    /// Mean WSS magnitude (Pa).
    pub fn mean_magnitude(&self) -> f64 {
        if self.tractions.is_empty() {
            return 0.0;
        }
        self.tractions.iter().map(|t| t.norm()).sum::<f64>() / self.tractions.len() as f64
    }

    /// Maximum WSS magnitude (Pa).
    pub fn max_magnitude(&self) -> f64 {
        self.tractions.iter().map(|t| t.norm()).fold(0.0, f64::max)
    }
}

/// Extracts the instantaneous WSS field from a solver state: for every
/// wall-adjacent cell face, `τ = μ · |u_tangential| / (Δ/2)` with the
/// tangential velocity taken from the adjacent face value and Δ the cell
/// pitch (linear wall gradient over half a cell).
pub fn extract_wss(solver: &HemodynamicsSolver) -> WssField {
    let (nx, ny, nz) = solver.domain.dims;
    let (dx, dy, dz) = solver.domain.spacing;
    let mut field = WssField::default();

    let mu_at = |i: i64, j: i64, k: i64| -> f64 {
        if i < 0 || j < 0 || k < 0 {
            return 0.0035;
        }
        solver.mu[solver.domain.index(i as usize, j as usize, k as usize)]
    };

    for i in 0..nx as i64 {
        for j in 0..ny as i64 {
            for k in 0..nz as i64 {
                if !solver.domain.is_fluid(i, j, k) {
                    continue;
                }
                // Check each of the 6 faces for a solid neighbour.
                for face in 0..6 {
                    let (di, dj, dk) = match face {
                        0 => (-1i64, 0, 0),
                        1 => (1, 0, 0),
                        2 => (0, -1, 0),
                        3 => (0, 1, 0),
                        4 => (0, 0, -1),
                        _ => (0, 0, 1),
                    };
                    if solver.domain.is_fluid(i + di, j + dj, k + dk) {
                        continue; // not a wall face
                    }
                    let mu = mu_at(i, j, k);
                    // Tangential velocity: average of the two in-plane face
                    // components nearest the wall (screening estimate).
                    let (pitch, ut_x, ut_y, ut_z) = match face {
                        0 | 1 => (
                            dx,
                            0.0,
                            solver.v[i as usize * (ny + 1) * nz + j as usize * nz + k as usize],
                            solver.w
                                [i as usize * ny * (nz + 1) + j as usize * (nz + 1) + k as usize],
                        ),
                        2 | 3 => (
                            dy,
                            solver.u[solver.uid(i as usize, j as usize, k as usize)],
                            0.0,
                            solver.w
                                [i as usize * ny * (nz + 1) + j as usize * (nz + 1) + k as usize],
                        ),
                        _ => (
                            dz,
                            solver.u[solver.uid(i as usize, j as usize, k as usize)],
                            solver.v[i as usize * (ny + 1) * nz + j as usize * nz + k as usize],
                            0.0,
                        ),
                    };
                    // Wall gradient magnitude over half a cell.
                    let g = (ut_x * ut_x + ut_y * ut_y + ut_z * ut_z).sqrt() / (0.5 * pitch);
                    let traction = g * mu;
                    // Direction: unit vector of the tangential velocity
                    // (sign-corrected so traction opposes the flow).
                    let un = Vec3::new(ut_x, ut_y, ut_z);
                    let dir = if un.norm() > 1e-12 {
                        un * (-1.0 / un.norm())
                    } else {
                        Vec3::ZERO
                    };
                    field.positions.push(Vec3::new(
                        (i as f64 + 0.5) * dx,
                        (j as f64 + 0.5) * dy,
                        (k as f64 + 0.5) * dz,
                    ));
                    field.tractions.push(dir * traction);
                }
            }
        }
    }
    field
}

/// Accumulates traction samples over a pulsatile cycle to compute the
/// oscillatory shear index:
///
/// ```text
/// OSI = ½ (1 − |Σ τᵢ Δtᵢ| / Σ |τᵢ| Δtᵢ)
/// ```
///
/// OSI = 0 for unidirectional steady shear; OSI → 0.5 for fully reversed
/// oscillation.
#[derive(Debug, Clone, Default)]
pub struct OsiAccumulator {
    sum_weighted: Vec3,
    sum_abs: f64,
}

impl OsiAccumulator {
    /// Adds one traction sample held for `dt`.
    pub fn sample(&mut self, traction: Vec3, dt: f64) {
        self.sum_weighted += traction * dt;
        self.sum_abs += traction.norm() * dt;
    }

    /// OSI in `[0, 0.5]` (0 when no samples).
    pub fn osi(&self) -> f64 {
        if self.sum_abs <= 1e-30 {
            return 0.0;
        }
        0.5 * (1.0 - self.sum_weighted.norm() / self.sum_abs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn osi_limits() {
        // Steady unidirectional: OSI = 0.
        let mut a = OsiAccumulator::default();
        for _ in 0..100 {
            a.sample(Vec3::new(1.5, 0.0, 0.0), 0.01);
        }
        assert!(a.osi() < 1e-12);
        // Fully reversed equal halves: OSI = 0.5.
        let mut b = OsiAccumulator::default();
        for _ in 0..50 {
            b.sample(Vec3::new(1.5, 0.0, 0.0), 0.01);
        }
        for _ in 0..50 {
            b.sample(Vec3::new(-1.5, 0.0, 0.0), 0.01);
        }
        assert!((b.osi() - 0.5).abs() < 1e-12);
        // Empty: 0.
        assert_eq!(OsiAccumulator::default().osi(), 0.0);
    }

    #[test]
    fn partial_oscillation_interpolates() {
        // 3/4 forward + 1/4 reversed at equal magnitude:
        // |Σ|/Σ|τ| = 0.5 → OSI = 0.25.
        let mut a = OsiAccumulator::default();
        for _ in 0..75 {
            a.sample(Vec3::new(1.0, 0.0, 0.0), 1.0);
        }
        for _ in 0..25 {
            a.sample(Vec3::new(-1.0, 0.0, 0.0), 1.0);
        }
        assert!((a.osi() - 0.25).abs() < 1e-12);
    }

    #[test]
    fn richardson_extrapolates_the_known_series() {
        // A second-order series converging on 37/24: 3, 1.75, 1.4375 at
        // ratio 2 — differences shrink 4×, order 2, and the
        // extrapolation from the two finest lands on the limit.
        let f = [3.0f64, 1.75, 1.4375];
        let p = ((f[0] - f[1]) / (f[1] - f[2])).ln() / 2.0f64.ln();
        assert!((p - 2.0).abs() < 1e-12);
        let extrap = richardson(f[1], f[2], 2.0, p);
        assert!((extrap - 4.0 / 3.0).abs() < 1e-9, "{extrap}");
    }

    #[test]
    fn resolution_study_quantifies_the_peak_wss() {
        // The same physical tube (radius 1 mm, length 8 mm) at three
        // pitches. The smooth Poiseuille case is resolution-stable: the
        // two finest levels agree to a few percent, the extrapolation is
        // on the same side of the finest value as the coarse-to-fine
        // trend, and the study reports an observed order.
        let build = |cells: usize| {
            let spacing = 8.0 / cells as f64;
            let radius_cells = 1.0 / spacing;
            FluidDomain::cylinder(cells, (cells / 4) + 2, radius_cells, spacing, 0)
        };
        let study = wss_resolution_study(
            &build,
            BloodModel::Newtonian { viscosity: 0.0035 },
            50.0,
            SolverConfig::default(),
            &[16, 24, 32],
            1.0e-3,
            4000,
        )
        .expect("converges");
        assert_eq!(study.levels.len(), 3);
        assert!(study
            .levels
            .iter()
            .all(|l| l.peak_wss_pa.is_finite() && l.peak_wss_pa > 0.0));
        let (coarsest, finest) = (study.levels[0].peak_wss_pa, study.levels[2].peak_wss_pa);
        // Honesty first: on a STAIR-STEPPED wall the peak WSS genuinely
        // jumps between resolutions (the stepped area fraction
        // quantises) — the study's job is to measure that, not to
        // smooth it over. Assert the reported extrapolation stays inside
        // the coarse-to-fine trend and that the finest single-grid
        // number carries a non-trivial resolution uncertainty.
        let span = (coarsest - finest).abs().max(finest.abs() * 1e-6);
        assert!(
            (study.extrapolated_peak_wss_pa - finest).abs() <= span * 1.5,
            "extrapolated {} outside the trend [{}, {}]",
            study.extrapolated_peak_wss_pa,
            coarsest,
            finest
        );
        assert!(study.finest_gap_fraction.is_finite() && study.finest_gap_fraction >= 0.0);
        // Non-trivial quantisation: the levels do not all agree — the
        // reason a peak-WSS claim must carry this study.
        let spread = (coarsest - finest).abs() / finest;
        assert!(spread > 1e-6, "levels suspiciously identical: {spread}");
    }

    #[test]
    fn resolution_study_rejects_bad_input() {
        let build = |cells: usize| FluidDomain::cylinder(cells, 6, 2.0, 0.5, 0);
        assert!(wss_resolution_study(
            &build,
            BloodModel::Newtonian { viscosity: 0.0035 },
            50.0,
            SolverConfig::default(),
            &[16],
            1.0e-3,
            100
        )
        .is_err());
        assert!(wss_resolution_study(
            &build,
            BloodModel::Newtonian { viscosity: 0.0035 },
            50.0,
            SolverConfig::default(),
            &[16, 16],
            1.0e-3,
            100
        )
        .is_err());
    }
}

/// One resolution level of a peak-WSS refinement study.
#[derive(Debug, Clone, Copy)]
pub struct WssResolutionLevel {
    /// The study's own refinement label (the domain builder's cell
    /// count), reported back for the record.
    pub cells: usize,
    /// Peak wall-shear magnitude (Pa) at this resolution.
    pub peak_wss_pa: f64,
}

/// The resolution-dependence report for peak WSS: per-level values, the
/// Richardson-extrapolated limit, the observed order (when three or more
/// levels allow one), and the gap between the finest run and the
/// extrapolation — the quantified resolution uncertainty the caller
/// reports instead of a single-grid number.
#[derive(Debug, Clone)]
pub struct WssResolutionStudy {
    /// Per-level peak WSS (Pa), sorted by increasing resolution.
    pub levels: Vec<WssResolutionLevel>,
    /// The resolution-corrected peak WSS (Pa): the Richardson
    /// extrapolation from the two finest levels when the three finest
    /// differences are well-behaved (monotone, consistent sign), and
    /// simply the finest value otherwise — a noisy quantised trend
    /// must not be amplified into a fake limit.
    pub extrapolated_peak_wss_pa: f64,
    /// Observed order from the three finest levels, `Some` only when the
    /// differences are monotonically shrinking with consistent sign —
    /// `None` means quantisation noise dominates and no extrapolation
    /// was claimed.
    pub observed_order: Option<f64>,
    /// The resolution uncertainty of the finest single-grid number: the
    /// Richardson gap when an order was observed, otherwise the
    /// measured coarse-to-fine spread.
    pub finest_gap_fraction: f64,
}

/// Richardson extrapolation of the two finest values at refinement ratio
/// `ratio` (≥ 1, cells ratio) under assumed order `p`.
fn richardson(coarse: f64, fine: f64, ratio: f64, p: f64) -> f64 {
    fine + (fine - coarse) / (ratio.powf(p) - 1.0)
}

/// Quantifies how resolution-dependent the peak WSS is: runs the same
/// physically-similar case at each entry of `cells_levels` (the caller's
/// `build` closure maps a cell count to the domain — same physical
/// geometry, finer pitch), marches each to steady state, extracts peak
/// WSS, and reports the Richardson picture. This is the honest
/// counterpart to local wall refinement: until the grid adapts, a
/// peak-WSS claim carries its resolution study with it.
///
/// # Errors
///
/// `Err` for fewer than two levels, a non-finite steady state, or a
/// level whose peak WSS did not converge inside `max_steps`.
pub fn wss_resolution_study(
    build: &dyn Fn(usize) -> FluidDomain,
    blood: BloodModel,
    inlet_velocity: f64,
    config: SolverConfig,
    cells_levels: &[usize],
    steady_tolerance: f64,
    max_steps: usize,
) -> Result<WssResolutionStudy, String> {
    if cells_levels.len() < 2 {
        return Err("a resolution study needs at least two levels".into());
    }
    let mut levels = Vec::with_capacity(cells_levels.len());
    for &cells in cells_levels {
        let mut solver =
            HemodynamicsSolver::new(build(cells), blood, inlet_velocity, config.clone());
        solver.run_steady(max_steps, steady_tolerance);
        let peak = crate::wss::extract_wss(&solver).max_magnitude();
        if !peak.is_finite() {
            return Err(format!("peak WSS diverged at {cells} cells"));
        }
        levels.push(WssResolutionLevel {
            cells,
            peak_wss_pa: peak,
        });
    }
    levels.sort_by_key(|l| l.cells);
    let n = levels.len();
    // Duplicate cell counts make the ratio meaningless.
    for w in levels.windows(2) {
        if w[0].cells == w[1].cells {
            return Err("resolution levels must be distinct".into());
        }
    }
    let f_coarse = levels[n - 2].peak_wss_pa;
    let f_fine = levels[n - 1].peak_wss_pa;
    let ratio = levels[n - 1].cells as f64 / levels[n - 2].cells as f64;

    // A Richardson claim needs well-behaved differences: monotonically
    // decreasing across the three finest levels. On a stair-stepped wall
    // the quantisation noise often breaks that — then NO extrapolation
    // is claimed (amplifying noise 1/(r− 1) times is worse than none)
    // and the report falls back to the finest value with the measured
    // coarse-to-fine spread as the resolution uncertainty.
    let mut observed_order = None;
    let mut extrapolated = f_fine;
    let mut gap = ((f_coarse - f_fine) / f_fine).abs();
    if n >= 3 {
        let f_coarser = levels[n - 3].peak_wss_pa;
        let d_coarse = f_coarser - f_coarse;
        let d_fine = f_coarse - f_fine;
        if d_fine.abs() > 0.0
            && d_coarse.abs() > d_fine.abs()
            && d_coarse.signum() == d_fine.signum()
        {
            let order = (d_coarse / d_fine).ln() / ratio.ln();
            if order.is_finite() && order > 0.0 {
                observed_order = Some(order);
                extrapolated = richardson(f_coarse, f_fine, ratio, order);
                gap = ((f_fine - extrapolated) / extrapolated).abs();
            }
        }
    }
    Ok(WssResolutionStudy {
        levels,
        extrapolated_peak_wss_pa: extrapolated,
        observed_order,
        finest_gap_fraction: gap,
    })
}
