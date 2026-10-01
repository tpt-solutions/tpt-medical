//! Passive scalar / temperature transport on the MAC grid — the fluid
//! side of conjugate heat transfer (the solid side, and a two-way coupled
//! interface, are the documented next slices).
//!
//! The scalar (blood temperature, contrast concentration) lives
//! cell-centered on the fluid mask; velocities are the solver's face
//! fields. Two operators, both in conservative flux form:
//!
//! - **Advection**: first-order upwind on the face velocities. Divergence-
//!   free face fields therefore conserve the total scalar exactly
//!   (telescoping fluxes), which the verification suite pins to machine
//!   precision on a discretely divergence-free field.
//! - **Diffusion**: centered 7-point Laplacian with thermal diffusivity
//!   `κ` (mm²/s; blood ≈ 0.11–0.12 mm²/s). Wall treatment per
//!   [`ScalarWall`]: `Fixed(tw)` — an isothermal wall, whose flux uses
//!   the half-cell distance to the wall temperature; `Insulated` — a
//!   zero-flux (adiabatic) wall.
//!
//! Stability: explicit update, so [`stable_time_step`] bounds `dt` by the
//! stricter of the advective (`s/u_max`) and diffusive (`s²/(6κ)` in
//! three dimensions) limits, times a safety factor.
//!
//! Wall faces (solid-adjacent or domain boundary) carry zero advective
//! flux when the neighbor is a no-flow solid; inflows through domain
//! boundaries carry the wall temperature under `Fixed` and the
//! zero-gradient interior value under `Insulated`.

use crate::domain::FluidDomain;

/// Wall thermal boundary condition for the scalar.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScalarWall {
    /// Isothermal wall at the given scalar value (e.g. `T_wall` in K or
    /// °C — the module is unit-agnostic) on every wall face.
    Fixed(f64),
    /// Adiabatic and impermeable wall: zero flux of either kind.
    Insulated,
    /// Per-face fixed values, ordered `[x_low, x_high, y_low, y_high,
    /// z_low, z_high]`; a `None` face is insulated. This is what makes a
    /// 1-D analytic test expressible on a 3-D box (hot one end, cold the
    /// other, insulated sides).
    Faces([Option<f64>; 6]),
}

/// Wall-face identifiers for [`ScalarWall::Faces`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WallFace {
    /// Low-x boundary.
    XLow,
    /// High-x boundary.
    XHigh,
    /// Low-y boundary.
    YLow,
    /// High-y boundary.
    YHigh,
    /// Low-z boundary.
    ZLow,
    /// High-z boundary.
    ZHigh,
}

/// The fixed temperature of one wall face, or `None` when insulated.
fn wall_value(wall: &ScalarWall, face: WallFace) -> Option<f64> {
    match wall {
        ScalarWall::Fixed(v) => Some(*v),
        ScalarWall::Insulated => None,
        ScalarWall::Faces(arr) => match face {
            WallFace::XLow => arr[0],
            WallFace::XHigh => arr[1],
            WallFace::YLow => arr[2],
            WallFace::YHigh => arr[3],
            WallFace::ZLow => arr[4],
            WallFace::ZHigh => arr[5],
        },
    }
}

/// The face a boundary neighbor lies on.
fn face_of(di: i64, dj: i64, dk: i64) -> WallFace {
    if di == -1 {
        WallFace::XLow
    } else if di == 1 {
        WallFace::XHigh
    } else if dj == -1 {
        WallFace::YLow
    } else if dj == 1 {
        WallFace::YHigh
    } else if dk == -1 {
        WallFace::ZLow
    } else {
        WallFace::ZHigh
    }
}

/// Face index on the u (x-face) array.
fn uidx(domain: &FluidDomain, i: usize, j: usize, k: usize) -> usize {
    let (_, ny, nz) = domain.dims;
    (i * ny + j) * nz + k
}

/// Face index on the v (y-face) array.
fn vidx(domain: &FluidDomain, i: usize, j: usize, k: usize) -> usize {
    let (_, ny, nz) = domain.dims;
    i * (ny + 1) * nz + j * nz + k
}

/// Face index on the w (z-face) array.
fn widx(domain: &FluidDomain, i: usize, j: usize, k: usize) -> usize {
    let (_, ny, nz) = domain.dims;
    i * ny * (nz + 1) + j * (nz + 1) + k
}

/// The explicit time step at or below which the update is stable:
/// `0.45 × min(s/|u_max|, s²/(6κ))` over the axes (the advective
/// cell-transit limit and the three-dimensional diffusive limit).
pub fn stable_time_step(domain: &FluidDomain, kappa: f64, u_max: f64) -> f64 {
    let (sx, sy, sz) = domain.spacing;
    let s_min = sx.min(sy).min(sz);
    let advective = if u_max > 0.0 {
        s_min / u_max
    } else {
        f64::INFINITY
    };
    let diffusive = if kappa > 0.0 {
        s_min * s_min / (6.0 * kappa)
    } else {
        f64::INFINITY
    };
    0.45 * advective.min(diffusive)
}

/// Advances the cell-centered scalar `t` by one explicit step of
/// advection–diffusion. `u`, `v`, `w` are the solver's face velocity
/// fields; `kappa` the diffusivity (mm²/s); `wall` the thermal boundary
/// condition. Solid cells keep their values and never receive flux.
///
/// # Panics
///
/// Panics if the arrays' lengths do not match the domain's face/cell
/// layout (a caller bug, not a runtime condition).
pub fn step(
    domain: &FluidDomain,
    t: &[f64],
    u: &[f64],
    v: &[f64],
    w: &[f64],
    kappa: f64,
    dt: f64,
    wall: ScalarWall,
) -> Vec<f64> {
    let (nx, ny, nz) = domain.dims;
    let (sx, sy, sz) = domain.spacing;
    let cell = |i: usize, j: usize, k: usize| (i * ny + j) * nz + k;
    assert_eq!(t.len(), nx * ny * nz, "scalar length mismatch");
    assert_eq!(u.len(), (nx + 1) * ny * nz, "u length mismatch");
    assert_eq!(v.len(), nx * (ny + 1) * nz, "v length mismatch");
    assert_eq!(w.len(), nx * ny * (nz + 1), "w length mismatch");

    let mut out = t.to_vec();

    // ---- Conservative upwind advection ----
    // Face fluxes (velocity × upwind scalar), zero through no-flow solids.
    let mut fx = vec![0.0f64; (nx + 1) * ny * nz];
    for i in 0..=nx {
        for j in 0..ny {
            for k in 0..nz {
                let uf = u[uidx(domain, i, j, k)];
                let idx = uidx(domain, i, j, k);
                if uf == 0.0 {
                    continue;
                }
                let left_ok = i > 0 && domain.is_fluid(i as i64 - 1, j as i64, k as i64);
                let right_ok = i < nx && domain.is_fluid(i as i64, j as i64, k as i64);
                if !left_ok || !right_ok {
                    // Solid-adjacent faces carry no flow (no-slip walls).
                    // A domain-boundary inflow under a fixed wall carries
                    // the wall temperature; insulated boundaries are also
                    // impermeable (no advective flux).
                    let inflow = (i == 0 && uf > 0.0) || (i == nx && uf < 0.0);
                    if inflow {
                        let face = if i == 0 {
                            WallFace::XLow
                        } else {
                            WallFace::XHigh
                        };
                        if let Some(tw) = wall_value(&wall, face) {
                            fx[idx] = uf * tw;
                        }
                    }
                    continue;
                }
                let left = cell(i - 1, j, k);
                let right = cell(i, j, k);
                fx[idx] = if uf >= 0.0 {
                    uf * t[left]
                } else {
                    uf * t[right]
                };
            }
        }
    }
    let mut fy = vec![0.0f64; nx * (ny + 1) * nz];
    for i in 0..nx {
        for j in 0..=ny {
            for k in 0..nz {
                let vf = v[vidx(domain, i, j, k)];
                let idx = vidx(domain, i, j, k);
                if vf == 0.0 {
                    continue;
                }
                let lower_ok = j > 0 && domain.is_fluid(i as i64, j as i64 - 1, k as i64);
                let upper_ok = j < ny && domain.is_fluid(i as i64, j as i64, k as i64);
                if !lower_ok || !upper_ok {
                    let inflow = (j == 0 && vf > 0.0) || (j == ny && vf < 0.0);
                    if inflow {
                        let face = if j == 0 {
                            WallFace::YLow
                        } else {
                            WallFace::YHigh
                        };
                        if let Some(tw) = wall_value(&wall, face) {
                            fy[idx] = vf * tw;
                        }
                    }
                    continue;
                }
                let lower = cell(i, j - 1, k);
                let upper = cell(i, j, k);
                fy[idx] = if vf >= 0.0 {
                    vf * t[lower]
                } else {
                    vf * t[upper]
                };
            }
        }
    }
    let mut fz = vec![0.0f64; nx * ny * (nz + 1)];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..=nz {
                let wf = w[widx(domain, i, j, k)];
                let idx = widx(domain, i, j, k);
                if wf == 0.0 {
                    continue;
                }
                let back_ok = k > 0 && domain.is_fluid(i as i64, j as i64, k as i64 - 1);
                let front_ok = k < nz && domain.is_fluid(i as i64, j as i64, k as i64);
                if !back_ok || !front_ok {
                    let inflow = (k == 0 && wf > 0.0) || (k == nz && wf < 0.0);
                    if inflow {
                        let face = if k == 0 {
                            WallFace::ZLow
                        } else {
                            WallFace::ZHigh
                        };
                        if let Some(tw) = wall_value(&wall, face) {
                            fz[idx] = wf * tw;
                        }
                    }
                    continue;
                }
                let back = cell(i, j, k - 1);
                let front = cell(i, j, k);
                fz[idx] = if wf >= 0.0 {
                    wf * t[back]
                } else {
                    wf * t[front]
                };
            }
        }
    }
    // Flux-form cell update.
    let ay = sy * sz;
    let ax = sx * sz;
    let az = sx * sy;
    let volume = sx * sy * sz;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if !domain.is_fluid(i as i64, j as i64, k as i64) {
                    continue;
                }
                let c = cell(i, j, k);
                let net = (fx[uidx(domain, i + 1, j, k)] - fx[uidx(domain, i, j, k)]) * ay
                    + (fy[vidx(domain, i, j + 1, k)] - fy[vidx(domain, i, j, k)]) * ax
                    + (fz[widx(domain, i, j, k + 1)] - fz[widx(domain, i, j, k)]) * az;
                out[c] -= dt * net / volume;
            }
        }
    }

    // ---- Centered diffusion with wall treatment ----
    let mut lap = vec![0.0f64; nx * ny * nz];
    let neighbor = |i: i64, j: i64, k: i64| -> Option<usize> {
        if i < 0 || j < 0 || k < 0 {
            return None;
        }
        let (iu, ju, ku) = (i as usize, j as usize, k as usize);
        if iu >= nx || ju >= ny || ku >= nz {
            return None;
        }
        if domain.is_fluid(i, j, k) {
            Some(cell(iu, ju, ku))
        } else {
            None // solid or boundary: handled by the wall term
        }
    };
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if !domain.is_fluid(i as i64, j as i64, k as i64) {
                    continue;
                }
                let c = cell(i, j, k);
                // Per-axis second differences: each direction divides by
                // its own spacing squared.
                let mut sum_x = 0.0f64;
                let mut sum_y = 0.0f64;
                let mut sum_z = 0.0f64;
                let mut wall_x = 0.0f64;
                let mut wall_y = 0.0f64;
                let mut wall_z = 0.0f64;
                for (di, dj, dk) in [
                    (1i64, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0),
                    (0, 0, 1),
                    (0, 0, -1),
                ] {
                    let (ii, jj, kk) = (i as i64 + di, j as i64 + dj, k as i64 + dk);
                    if let Some(nc) = neighbor(ii, jj, kk) {
                        let d = t[nc] - t[c];
                        if di != 0 {
                            sum_x += d;
                        } else if dj != 0 {
                            sum_y += d;
                        } else {
                            sum_z += d;
                        }
                    } else {
                        // A wall face (solid or domain boundary) at its
                        // fixed temperature, half a cell from the center
                        // (factor 2); insulated faces contribute nothing.
                        if let Some(tw) = wall_value(&wall, face_of(di, dj, dk)) {
                            let d = 2.0 * (tw - t[c]);
                            if di != 0 {
                                wall_x += d;
                            } else if dj != 0 {
                                wall_y += d;
                            } else {
                                wall_z += d;
                            }
                        }
                    }
                    // Insulated walls contribute zero flux.
                }
                lap[cell(i, j, k)] = (sum_x + wall_x) / (sx * sx)
                    + (sum_y + wall_y) / (sy * sy)
                    + (sum_z + wall_z) / (sz * sz);
            }
        }
    }
    for i in 0..nx * ny * nz {
        out[i] += dt * kappa * lap[i];
    }
    out
}

/// Advances the cell-centered scalar by one explicit step of advection
/// (fluid cells only — the MAC velocities are zero in the solid) plus
/// diffusion **through both regions**, with the interface faces taking
/// the **harmonic-mean** conductivity `2 κ_f κ_s/(κ_f + κ_s)`.
///
/// The harmonic mean is not a refinement choice: for cell-centered finite
/// volumes it makes the interface face's resistance the exact series sum
/// of the two half-cell resistances, which is why the two-layer steady
/// test reproduces the analytic composite-wall solution to machine
/// precision (the classic argument for harmonic means at material
/// interfaces).
///
/// This is the two-way coupled fluid–solid boundary of the CHT item: the
/// fluid's wall temperature is no longer prescribed — it *is* the solid
/// cell's temperature, and the flux through the interface is continuous.
///
/// # Panics
///
/// Panics if the arrays' lengths do not match the domain layout.
#[allow(clippy::too_many_lines)]
pub fn conjugate_step(
    domain: &FluidDomain,
    t: &[f64],
    u: &[f64],
    v: &[f64],
    w: &[f64],
    kappa_fluid: f64,
    kappa_solid: f64,
    dt: f64,
    wall: &ScalarWall,
) -> Vec<f64> {
    let (nx, ny, nz) = domain.dims;
    let (sx, sy, sz) = domain.spacing;
    let cell = |i: usize, j: usize, k: usize| (i * ny + j) * nz + k;
    assert_eq!(t.len(), nx * ny * nz, "scalar length mismatch");
    assert_eq!(u.len(), (nx + 1) * ny * nz, "u length mismatch");
    assert_eq!(v.len(), nx * (ny + 1) * nz, "v length mismatch");
    assert_eq!(w.len(), nx * ny * (nz + 1), "w length mismatch");
    let kappa_cell = |i: usize, j: usize, k: usize| {
        if domain.is_fluid(i as i64, j as i64, k as i64) {
            kappa_fluid
        } else {
            kappa_solid
        }
    };
    // Harmonic mean at a face between two cells.
    let harmonic = |ka: f64, kb: f64| 2.0 * ka * kb / (ka + kb);

    let mut out = t.to_vec();

    // ---- Advection: fluid cells only (velocities are zero in solids) ----
    let mut fx = vec![0.0f64; (nx + 1) * ny * nz];
    for i in 0..=nx {
        for j in 0..ny {
            for k in 0..nz {
                let uf = u[uidx(domain, i, j, k)];
                let idx = uidx(domain, i, j, k);
                if uf == 0.0 {
                    continue;
                }
                let left_ok = i > 0 && domain.is_fluid(i as i64 - 1, j as i64, k as i64);
                let right_ok = i < nx && domain.is_fluid(i as i64, j as i64, k as i64);
                if !left_ok || !right_ok {
                    // No flow through solids; boundary inflows under a
                    // fixed wall face carry the wall value.
                    if (i == 0 && uf > 0.0) || (i == nx && uf < 0.0) {
                        let face = if i == 0 {
                            WallFace::XLow
                        } else {
                            WallFace::XHigh
                        };
                        if let Some(tw) = wall_value(wall, face) {
                            fx[idx] = uf * tw;
                        }
                    }
                    continue;
                }
                let left = cell(i - 1, j, k);
                let right = cell(i, j, k);
                fx[idx] = if uf >= 0.0 {
                    uf * t[left]
                } else {
                    uf * t[right]
                };
            }
        }
    }
    let mut fy = vec![0.0f64; nx * (ny + 1) * nz];
    for i in 0..nx {
        for j in 0..=ny {
            for k in 0..nz {
                let vf = v[vidx(domain, i, j, k)];
                let idx = vidx(domain, i, j, k);
                if vf == 0.0 {
                    continue;
                }
                let lower_ok = j > 0 && domain.is_fluid(i as i64, j as i64 - 1, k as i64);
                let upper_ok = j < ny && domain.is_fluid(i as i64, j as i64, k as i64);
                if !lower_ok || !upper_ok {
                    if (j == 0 && vf > 0.0) || (j == ny && vf < 0.0) {
                        let face = if j == 0 {
                            WallFace::YLow
                        } else {
                            WallFace::YHigh
                        };
                        if let Some(tw) = wall_value(wall, face) {
                            fy[idx] = vf * tw;
                        }
                    }
                    continue;
                }
                let lower = cell(i, j - 1, k);
                let upper = cell(i, j, k);
                fy[idx] = if vf >= 0.0 {
                    vf * t[lower]
                } else {
                    vf * t[upper]
                };
            }
        }
    }
    let mut fz = vec![0.0f64; nx * ny * (nz + 1)];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..=nz {
                let wf = w[widx(domain, i, j, k)];
                let idx = widx(domain, i, j, k);
                if wf == 0.0 {
                    continue;
                }
                let back_ok = k > 0 && domain.is_fluid(i as i64, j as i64, k as i64 - 1);
                let front_ok = k < nz && domain.is_fluid(i as i64, j as i64, k as i64);
                if !back_ok || !front_ok {
                    if (k == 0 && wf > 0.0) || (k == nz && wf < 0.0) {
                        let face = if k == 0 {
                            WallFace::ZLow
                        } else {
                            WallFace::ZHigh
                        };
                        if let Some(tw) = wall_value(wall, face) {
                            fz[idx] = wf * tw;
                        }
                    }
                    continue;
                }
                let back = cell(i, j, k - 1);
                let front = cell(i, j, k);
                fz[idx] = if wf >= 0.0 {
                    wf * t[back]
                } else {
                    wf * t[front]
                };
            }
        }
    }
    let ay = sy * sz;
    let ax = sx * sz;
    let az = sx * sy;
    let volume = sx * sy * sz;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if !domain.is_fluid(i as i64, j as i64, k as i64) {
                    continue;
                }
                let c = cell(i, j, k);
                let net = (fx[uidx(domain, i + 1, j, k)] - fx[uidx(domain, i, j, k)]) * ay
                    + (fy[vidx(domain, i, j + 1, k)] - fy[vidx(domain, i, j, k)]) * ax
                    + (fz[widx(domain, i, j, k + 1)] - fz[widx(domain, i, j, k)]) * az;
                out[c] -= dt * net / volume;
            }
        }
    }

    // ---- Diffusion through BOTH regions, harmonic-mean interfaces ----
    let mut out2 = out;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let c = cell(i, j, k);
                let kc = kappa_cell(i, j, k);
                let mut flux_sum = 0.0f64;
                for (di, dj, dk) in [
                    (1i64, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0),
                    (0, 0, 1),
                    (0, 0, -1),
                ] {
                    let (ii, jj, kk) = (i as i64 + di, j as i64 + dj, k as i64 + dk);
                    let in_bounds = ii >= 0
                        && jj >= 0
                        && kk >= 0
                        && (ii as usize) < nx
                        && (jj as usize) < ny
                        && (kk as usize) < nz;
                    let s2 = if di != 0 {
                        sx * sx
                    } else if dj != 0 {
                        sy * sy
                    } else {
                        sz * sz
                    };
                    if in_bounds {
                        let nc = cell(ii as usize, jj as usize, kk as usize);
                        // Face conductivity: harmonic mean across a
                        // material interface, own conductivity otherwise.
                        let kf = harmonic(kc, kappa_cell(ii as usize, jj as usize, kk as usize));
                        flux_sum += kf * (t[nc] - t[c]) / s2;
                    } else if let Some(tw) = wall_value(wall, face_of(di, dj, dk)) {
                        // A wall face at its fixed temperature, half a
                        // cell from the center (factor 2); insulated
                        // faces contribute nothing.
                        flux_sum += 2.0 * kc * (tw - t[c]) / s2;
                    }
                }
                // Unit volumetric heat capacity: diffusivity = conductivity.
                out2[c] = t[c] + dt * flux_sum;
            }
        }
    }
    out2
}

/// Volume-weighted mean scalar over fluid cells (the bulk temperature).
pub fn bulk(domain: &FluidDomain, t: &[f64]) -> f64 {
    let (nx, ny, nz) = domain.dims;
    let mut total = 0.0f64;
    let mut volume = 0.0f64;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if domain.is_fluid(i as i64, j as i64, k as i64) {
                    total += t[(i * ny + j) * nz + k];
                    volume += 1.0;
                }
            }
        }
    }
    if volume > 0.0 {
        total / volume
    } else {
        f64::NAN
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_fluid_box(n: usize) -> FluidDomain {
        FluidDomain::from_mask((n, n, n), (1.0, 1.0, 1.0), vec![true; n * n * n], 0, true)
    }

    fn zero_flow(domain: &FluidDomain) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let (nx, ny, nz) = domain.dims;
        (
            vec![0.0; (nx + 1) * ny * nz],
            vec![0.0; nx * (ny + 1) * nz],
            vec![0.0; nx * ny * (nz + 1)],
        )
    }

    /// Pure diffusion, insulated walls: `T = 1 + sin(pi x/L)` decays as
    /// `1 + sin(pi x/L) e^(-kappa_rate t)` — the Fourier mode is the
    /// exact eigenfunction of the centered Laplacian with zero-flux
    /// ends, and the asserted rate is the *discrete* eigenvalue
    /// `(2 - 2 cos(pi/N)) kappa`, which the centered operator reproduces
    /// to machine precision.
    #[test]
    fn insulated_diffusion_matches_discrete_fourier_decay() {
        let n = 32;
        let domain = all_fluid_box(n);
        let kappa = 0.12;
        let l = n as f64;
        // The cell-centered Neumann-Neumann eigenfunction for k = 1:
        // cos(pi (i + 1/2) / N).
        let mut t: Vec<f64> = (0..n * n * n)
            .map(|idx| {
                let i = idx / (n * n);
                1.0 + (core::f64::consts::PI * (i as f64 + 0.5) / l).cos()
            })
            .collect();
        let (u, v, w) = zero_flow(&domain);
        let rate = kappa * (2.0 - 2.0 * (core::f64::consts::PI / n as f64).cos());
        let dt = 0.4 / 6.0;
        let steps = 200;
        for _ in 0..steps {
            t = step(&domain, &t, &u, &v, &w, kappa, dt, ScalarWall::Insulated);
        }
        // The eigenmode decays by the discrete amplification factor
        // (1 - dt*rate) per step — exactly, since the mode is an exact
        // eigenvector of the operator.
        let expected =
            1.0 + (core::f64::consts::PI * 16.5 / l).cos() * (1.0 - dt * rate).powi(steps);
        let center = t[(16 * n + 8) * n + 8];
        assert!(
            (center - expected).abs() < 1e-9,
            "center {center} vs discrete-eigenfunction {expected}"
        );
    }

    /// Fixed (isothermal) walls at zero: the 3-D product mode
    /// `prod sin(pi (a + 1/2)/N)` is the exact eigenfunction of the
    /// half-cell Dirichlet closure on all three axes (the wall value
    /// lives at the faces, half a cell outside the first centers), with
    /// the summed discrete eigenvalue `3 kappa 4 sin^2(pi/2N)`.
    #[test]
    fn fixed_wall_diffusion_matches_discrete_fourier_decay() {
        let n = 32;
        let domain = all_fluid_box(n);
        let kappa = 0.12;
        let l = n as f64;
        let mode = |a: usize| (core::f64::consts::PI * (a as f64 + 0.5) / l).sin();
        let mut t: Vec<f64> = (0..n * n * n)
            .map(|idx| {
                let (i, j, k) = (idx / (n * n), (idx / n) % n, idx % n);
                mode(i) * mode(j) * mode(k)
            })
            .collect();
        let (u, v, w) = zero_flow(&domain);
        let rate = 3.0 * kappa * (2.0 - 2.0 * (core::f64::consts::PI / n as f64).cos());
        let dt = 0.4 / 6.0;
        let steps = 200;
        for _ in 0..steps {
            t = step(&domain, &t, &u, &v, &w, kappa, dt, ScalarWall::Fixed(0.0));
        }
        let expected = mode(16).powi(3) * (1.0 - dt * rate).powi(steps);
        let center = t[(16 * n + 16) * n + 16];
        assert!(
            (center - expected).abs() < 1e-9,
            "center {center} vs {expected}"
        );
    }

    /// Conservative advection: on a discretely divergence-free face field
    /// (built from a stream function), the flux-form update conserves the
    /// total scalar to machine precision - the telescoping-flux property
    /// the conservative form exists for.
    #[test]
    fn advective_transport_conserves_total_scalar_exactly() {
        let n = 16;
        let domain = all_fluid_box(n);
        let psi = |i: usize, j: usize, _k: usize| {
            (core::f64::consts::PI * i as f64 / n as f64).sin()
                * (core::f64::consts::PI * j as f64 / n as f64).sin()
        };
        let mut u = vec![0.0f64; (n + 1) * n * n];
        for i in 0..=n {
            for j in 0..n {
                for k in 0..n {
                    u[(i * n + j) * n + k] = (psi(i, j, k) - psi(i + 1, j, k)) * 5.0;
                }
            }
        }
        let mut v = vec![0.0f64; n * (n + 1) * n];
        for i in 0..n {
            for j in 0..=n {
                for k in 0..n {
                    v[(i * (n + 1) + j) * n + k] = (psi(i, j + 1, k) - psi(i, j, k)) * 5.0;
                }
            }
        }
        let w = vec![0.0f64; n * n * (n + 1)];
        let kappa = 0.0; // pure advection

        let mut t: Vec<f64> = (0..n * n * n)
            .map(|idx| {
                let (i, j, k) = (idx / (n * n), (idx / n) % n, idx % n);
                if (i as f64 - 8.0).powi(2) + (j as f64 - 8.0).powi(2) + (k as f64 - 8.0).powi(2)
                    < 9.0
                {
                    1.0
                } else {
                    0.0
                }
            })
            .collect();
        let initial: f64 = t.iter().sum();
        let dt = stable_time_step(&domain, kappa, 5.0);
        for _ in 0..50 {
            t = step(&domain, &t, &u, &v, &w, kappa, dt, ScalarWall::Insulated);
        }
        let final_total: f64 = t.iter().sum();
        assert!(
            (final_total - initial).abs() < 1e-9 * initial.abs().max(1.0),
            "total scalar drifted: {initial} -> {final_total}"
        );
    }

    /// Thermal entry behavior in a developed Poiseuille tube: with cold
    /// walls the bulk temperature decays monotonically toward the wall
    /// value; with insulated walls the advected heat is retained (no
    /// spurious wall losses). The exact Nusselt number is deliberately
    /// not asserted on a stair-step coarse grid.
    #[test]
    fn thermal_entry_behaves_physically() {
        let domain = FluidDomain::cylinder(16, 10, 3.5, 0.5, 0);
        let blood = crate::BloodModel::Newtonian { viscosity: 0.08 };
        let config = crate::solver::SolverConfig {
            dt: 5.0e-4,
            poisson_iterations: 400,
            include_convection: false,
            viscosity_relaxation: 0.2,
            density: 1.06e-3,
            pressure_solver: crate::solver::PressureSolver::Sor,
        };
        let mut solver = crate::solver::HemodynamicsSolver::new(domain.clone(), blood, 5.0, config);
        solver.run_steady(900, 1e-4);

        let (nx, ny, nz) = domain.dims;
        let init = || {
            let mut t = vec![0.0f64; nx * ny * nz];
            for i in 0..2 {
                for j in 0..ny {
                    for k in 0..nz {
                        if domain.is_fluid(i as i64, j as i64, k as i64) {
                            t[(i * ny + j) * nz + k] = 1.0;
                        }
                    }
                }
            }
            t
        };
        let kappa = 0.12;
        let dt = stable_time_step(&domain, kappa, 10.0);

        // Cold walls: bulk decays monotonically and substantially.
        let initial_bulk = bulk(&domain, &init());
        let cold = ScalarWall::Fixed(0.0);
        let mut t = init();
        let mut bulk_prev = initial_bulk;
        for _ in 0..200 {
            t = step(
                &domain, &t, &solver.u, &solver.v, &solver.w, kappa, dt, cold,
            );
            let b = bulk(&domain, &t);
            assert!(
                b <= bulk_prev + 1e-12,
                "bulk must not rise: {b} after {bulk_prev}"
            );
            bulk_prev = b;
        }
        // Cold walls strip the heat: most of it leaves the domain.
        assert!(
            bulk_prev < 0.6 * initial_bulk,
            "heat must be lost to the cold walls: {bulk_prev} vs {initial_bulk}"
        );

        // Insulated walls are also impermeable: the advected scalar is
        // conserved exactly (the flowing-flow conservation check).
        let insulated = ScalarWall::Insulated;
        let mut t = init();
        let mut bulk_ins = initial_bulk;
        for _ in 0..200 {
            t = step(
                &domain, &t, &solver.u, &solver.v, &solver.w, kappa, dt, insulated,
            );
            bulk_ins = bulk(&domain, &t);
        }
        assert!(
            (bulk_ins - initial_bulk).abs() < 1e-9,
            "insulated walls must conserve the scalar: {bulk_ins} vs {initial_bulk}"
        );
    }

    #[test]
    fn stable_time_step_bounds_both_limits() {
        let domain = all_fluid_box(8);
        let dt_d = stable_time_step(&domain, 0.12, 0.0);
        assert!((dt_d - 0.45 * (1.0 / (6.0 * 0.12))).abs() < 1e-12);
        let dt_a = stable_time_step(&domain, 0.0, 2.0);
        assert!((dt_a - 0.45 * 0.5).abs() < 1e-12);
        let dt_b = stable_time_step(&domain, 0.12, 2.0);
        assert!((dt_b - dt_d.min(dt_a)).abs() < 1e-12);
    }
}

/// The composite-wall steady state: a fluid half and a solid half,
/// hot wall on the left, cold on the right, insulated elsewhere,
/// run to equilibrium under `conjugate_step`. The harmonic-mean
/// interface makes the interface face's resistance the exact series
/// sum of the two half-cell resistances, so the discrete steady
/// state IS the analytic composite-wall solution — every cell
/// matches the resistance-chain potential to machine precision, and
/// the interface temperature equals the two-layer formula
/// `T_int = T_hot - (T_hot - T_cold) R_hot/(R_hot + R_cold)`.
#[test]
fn conjugate_two_layer_steady_state_matches_series_resistance() {
    let n = 32;
    let dims = (n, 2, 2);
    // Fluid for x < 16, solid for x >= 16.
    let mut mask = vec![false; n * 2 * 2];
    for i in 0..n {
        for j in 0..2 {
            for k in 0..2 {
                mask[(i * 2 + j) * 2 + k] = i < 16;
            }
        }
    }
    let domain = FluidDomain::from_mask(dims, (1.0, 1.0, 1.0), mask, 0, true);
    let (k_f, k_s) = (0.12, 1.0);
    let (t_hot, t_cold) = (1.0f64, 0.0f64);
    let wall = ScalarWall::Faces([Some(t_hot), Some(t_cold), None, None, None, None]);

    let zeros = (
        vec![0.0; (n + 1) * 2 * 2],
        vec![0.0; n * 3 * 2],
        vec![0.0; n * 2 * 3],
    );
    let mut t = vec![0.5f64; n * 2 * 2];
    // March to steady state: the slowest mode is the fluid slab's
    // diffusion time (2L_f/pi)^2/k_f ~ 860 s; 1.4e5 steps x 0.15 s
    // ~ 2.1e4 s is ~24 of those time constants, putting the residual
    // transient below the 1e-9 assertion.
    let dt = 0.9 / 6.0;
    for _ in 0..150_000 {
        t = conjugate_step(
            &domain, &t, &zeros.0, &zeros.1, &zeros.2, k_f, k_s, dt, &wall,
        );
    }

    // Analytic discrete steady state from the resistance chain: q =
    // (T_hot - T_cold)/R_total with the half-cell wall resistances
    // and the (exact) harmonic-mean interface resistance.
    let dx = 1.0;
    let r_half_f = dx / (2.0 * k_f);
    let r_half_s = dx / (2.0 * k_s);
    let k_int = 2.0 * k_f * k_s / (k_f + k_s);
    let r_int = dx / k_int;
    let r_total = r_half_f + 15.0 * (dx / k_f) + r_int + 15.0 * (dx / k_s) + r_half_s;
    let q = (t_hot - t_cold) / r_total;

    let mut expected = Vec::with_capacity(n);
    let mut potential = t_hot;
    // Wall half-cell first, then one resistance per cell-to-cell face.
    potential -= q * r_half_f;
    expected.push(potential);
    for i in 1..n {
        let r_step = if i == 16 {
            r_int
        } else if i < 16 {
            dx / k_f
        } else {
            dx / k_s
        };
        potential -= q * r_step;
        expected.push(potential);
    }
    for i in 0..n {
        for j in 0..2 {
            for k in 0..2 {
                let got = t[(i * 2 + j) * 2 + k];
                assert!(
                    (got - expected[i]).abs() < 1e-9,
                    "cell ({i},{j},{k}): {got} vs {}",
                    expected[i]
                );
            }
        }
    }
    // Interface temperature equals the two-layer formula.
    let r_fluid_side = r_half_f + 15.0 * (dx / k_f);
    let t_interface = t_hot - (t_hot - t_cold) * r_fluid_side / r_total;
    assert!(
        (t[(15 * 2) * 2] - t_interface).abs() < 1e-9,
        "interface cell {} vs {}",
        t[(15 * 2) * 2],
        t_interface
    );
}

/// Two-region conservation: with insulated outer walls, the fluid and
/// solid energies sum to a constant even as heat redistributes
/// through the interface.
#[test]
fn conjugate_conservation_across_regions() {
    let n = 16;
    let dims = (n, 2, 2);
    let mut mask = vec![false; n * 2 * 2];
    for i in 0..n {
        for j in 0..2 {
            for k in 0..2 {
                mask[(i * 2 + j) * 2 + k] = i < 8;
            }
        }
    }
    let domain = FluidDomain::from_mask(dims, (1.0, 1.0, 1.0), mask, 0, true);
    let (k_f, k_s) = (0.12, 1.0);
    let wall = ScalarWall::Insulated;
    let zeros = (
        vec![0.0; (n + 1) * 2 * 2],
        vec![0.0; n * 3 * 2],
        vec![0.0; n * 2 * 3],
    );
    // Hot fluid, cold solid.
    let mut t = vec![0.0f64; n * 2 * 2];
    for i in 0..n {
        for j in 0..2 {
            for k in 0..2 {
                t[(i * 2 + j) * 2 + k] = if i < 8 { 1.0 } else { 0.0 };
            }
        }
    }
    let initial: f64 = t.iter().sum();
    let dt = 0.9 / 6.0;
    // Conservation is checked along the way, not just at the end.
    // ~40 fluid-internal time constants: the equilibration assertion
    // below is a real 1e-3 of the initial unit contrast.
    for step_n in 0..12000 {
        t = conjugate_step(
            &domain, &t, &zeros.0, &zeros.1, &zeros.2, k_f, k_s, dt, &wall,
        );
        if step_n % 500 == 0 {
            let total: f64 = t.iter().sum();
            assert!(
                (total - initial).abs() < 1e-9,
                "two-region energy drifted at step {step_n}: {initial} -> {total}"
            );
        }
    }
    let final_total: f64 = t.iter().sum();
    assert!(
        (final_total - initial).abs() < 1e-9,
        "two-region energy drifted: {initial} -> {final_total}"
    );
    // Interface exchange time ~ C/G ~ 38 s; the slower fluid-internal
    // mode is ~54 s, so 900 s is ~17 of those: the spread is below
    // 1e-3 of the initial unit contrast.
    let spread = t.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - t.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(spread < 1e-3, "not equilibrated: spread {spread}");
}
