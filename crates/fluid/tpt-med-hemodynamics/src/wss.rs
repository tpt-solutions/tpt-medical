//! Wall shear stress fields and oscillatory shear index.

use crate::solver::HemodynamicsSolver;
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
}
