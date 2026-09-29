//! Projection-method Navier–Stokes solver on the staggered (MAC) grid.
//!
//! Algorithm per time step:
//! 1. Boundary application (inlet plug velocity, no-slip walls, zero-
//!    gradient outlet).
//! 2. Explicit sub-step: `u* = u + dt(−(u·∇)u + ν∇²u)` on active faces
//!    (both adjacent cells fluid); wall faces pinned to zero.
//! 3. Projection: solve `∇²φ = ∇·u*/dt` (Neumann walls, Dirichlet φ = 0
//!    outlet layer) with the configured [`PressureSolver`] — SOR
//!    (Gauss–Seidel with over-relaxation, the default) or Jacobi-
//!    preconditioned conjugate gradient; `u = u* − dt ∇φ`, `p += φ`.
//! 4. Non-Newtonian viscosity update from the local shear rate.

use crate::blood::BloodModel;
use crate::domain::FluidDomain;

/// Pressure-Poisson solver choice. The discrete operator is identical for
/// both: SOR sweeps it in place (cheap per iteration, slow convergence on
/// large grids); conjugate gradient iterates matrix-free with Jacobi
/// (diagonal) preconditioning and exits on a relative-residual criterion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PressureSolver {
    /// Red-black-style in-place Gauss–Seidel with over-relaxation (ω = 1.9).
    #[default]
    Sor,
    /// Preconditioned conjugate gradient — the larger-grid choice: its
    /// iteration count grows with √(condition number) rather than with the
    /// grid's graph diameter.
    ConjugateGradient,
}

/// Solver run configuration.
#[derive(Debug, Clone)]
pub struct SolverConfig {
    /// Physical time step (s).
    pub dt: f64,
    /// Iteration cap for the pressure solve (SOR sweeps or CG iterations).
    pub poisson_iterations: usize,
    /// Include convective term (disable for Stokes-like creeping flows and
    /// faster convergence).
    pub include_convection: bool,
    /// Under-relaxation for the non-Newtonian viscosity update.
    pub viscosity_relaxation: f64,
    /// Density of blood (g/mm³ = kg/m³ × 1e-6; 1.06e-3 g/mm³ = 1060 kg/m³).
    pub density: f64,
    /// Pressure-Poisson solve method.
    pub pressure_solver: PressureSolver,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            dt: 2.0e-4,
            poisson_iterations: 400,
            include_convection: false,
            viscosity_relaxation: 0.2,
            density: 1.06e-3,
            pressure_solver: PressureSolver::Sor,
        }
    }
}

/// Scalars reported after a steady run.
#[derive(Debug, Clone, Copy)]
pub struct SteadyStats {
    /// Marching steps performed.
    pub steps: usize,
    /// Maximum velocity magnitude (mm/s).
    pub max_velocity: f64,
    /// Flow through the inlet face (mm³/s).
    pub inlet_flow: f64,
    /// Flow through the outlet face (mm³/s).
    pub outlet_flow: f64,
    /// Mean pressure over inlet-plane fluid cells (internal units).
    pub mean_pressure_inlet: f64,
    /// Mean pressure over outlet-plane fluid cells (internal units).
    pub mean_pressure_outlet: f64,
    /// Pressure drop inlet→outlet (internal units).
    pub pressure_drop: f64,
}

/// MAC-grid hemodynamics solver. Flow runs along **x** (the low-x face is
/// the inlet, high-x the outlet); other axes are fully wall-bounded by the
/// mask.
pub struct HemodynamicsSolver {
    /// Domain.
    pub domain: FluidDomain,
    /// Blood model.
    pub blood: BloodModel,
    /// Configuration.
    pub config: SolverConfig,
    /// x-face velocities (mm/s), dims (nx+1, ny, nz).
    pub u: Vec<f64>,
    /// y-face velocities (mm/s), dims (nx, ny+1, nz).
    pub v: Vec<f64>,
    /// z-face velocities (mm/s), dims (nx, ny, nz+1).
    pub w: Vec<f64>,
    /// Accumulated pressure potential Π = (p/ρ)·dt-units (mm²/s² per
    /// step); its gradient enters the predictor, its increments come from
    /// the projection.
    pub p: Vec<f64>,
    /// Cell apparent viscosity (Pa·s).
    pub mu: Vec<f64>,
    /// Inlet speed (mm/s).
    pub inlet_velocity: f64,
    /// Iterations used by the most recent pressure solve (SOR sweeps or CG
    /// iterations) — the number to watch when cost accounting a run.
    pub last_pressure_solve_iterations: usize,
    /// Dirichlet value of the projection at the outlet layer, in the
    /// solver's internal mm²/s² units. Zero (the gauge) unless a
    /// Windkessel-driven run sets it.
    pub(crate) outlet_anchor: f64,
}

impl HemodynamicsSolver {
    /// Creates a solver over `domain` with the given blood model and inlet
    /// plug velocity.
    pub fn new(
        domain: FluidDomain,
        blood: BloodModel,
        inlet_velocity: f64,
        config: SolverConfig,
    ) -> Self {
        assert_eq!(domain.flow_axis, 0, "v0 solver drives flow along x only");
        let (nx, ny, nz) = domain.dims;
        let mu0 = match blood {
            BloodModel::Newtonian { viscosity } => viscosity,
            _ => 0.0035,
        };
        Self {
            u: vec![0.0; (nx + 1) * ny * nz],
            v: vec![0.0; nx * (ny + 1) * nz],
            w: vec![0.0; nx * ny * (nz + 1)],
            p: vec![0.0; nx * ny * nz],
            mu: vec![mu0; nx * ny * nz],
            domain,
            blood,
            config,
            inlet_velocity,
            last_pressure_solve_iterations: 0,
            outlet_anchor: 0.0,
        }
    }

    #[inline]
    pub(crate) fn uid(&self, i: usize, j: usize, k: usize) -> usize {
        let (_, ny, nz) = self.domain.dims;
        (i * ny + j) * nz + k
    }

    /// Active x-face: both adjacent cells are fluid.
    pub(crate) fn ux_active(&self, i: usize, j: usize, k: usize) -> bool {
        self.domain.is_fluid(i as i64, j as i64, k as i64)
            && self.domain.is_fluid(i as i64 - 1, j as i64, k as i64)
    }

    fn vy_active(&self, i: usize, j: usize, k: usize) -> bool {
        self.domain.is_fluid(i as i64, j as i64, k as i64)
            && self.domain.is_fluid(i as i64, j as i64 - 1, k as i64)
    }

    fn wz_active(&self, i: usize, j: usize, k: usize) -> bool {
        self.domain.is_fluid(i as i64, j as i64, k as i64)
            && self.domain.is_fluid(i as i64, j as i64, k as i64 - 1)
    }

    /// Applies BCs: plug inlet at x=0, zero-gradient outlet at x=nx,
    /// no-slip walls elsewhere.
    pub fn apply_boundary(&mut self) {
        let (nx, ny, nz) = self.domain.dims;
        // Inlet: x = 0 faces with fluid cell (0,j,k).
        for j in 0..ny {
            for k in 0..nz {
                let idx = self.uid(0, j, k);
                if self.domain.is_fluid(0, j as i64, k as i64) {
                    self.u[idx] = self.inlet_velocity;
                } else {
                    self.u[idx] = 0.0;
                }
            }
        }
        // Outlet: zero gradient: u[nx] = u[nx-1] where fluid at nx-1.
        for j in 0..ny {
            for k in 0..nz {
                let idx_out = self.uid(nx, j, k);
                let idx_in = self.uid(nx - 1, j, k);
                if self.domain.is_fluid((nx - 1) as i64, j as i64, k as i64) {
                    self.u[idx_out] = self.u[idx_in];
                } else {
                    self.u[idx_out] = 0.0;
                }
            }
        }
        // No-slip on all v/w faces and on interior u faces adjacent to
        // solids.
        let (_, ny1, _) = (nx, ny + 1, nz);
        let _ = ny1;
        for i in 0..nx {
            for j in 0..=ny {
                for k in 0..nz {
                    let idx = i * (ny + 1) * nz + j * nz + k;
                    self.v[idx] = if self.vy_active(i, j, k) {
                        self.v[idx]
                    } else {
                        0.0
                    };
                }
            }
        }
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..=nz {
                    let idx = i * ny * (nz + 1) + j * (nz + 1) + k;
                    self.w[idx] = if self.wz_active(i, j, k) {
                        self.w[idx]
                    } else {
                        0.0
                    };
                }
            }
        }
        for i in 1..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = self.uid(i, j, k);
                    self.u[idx] = if self.ux_active(i, j, k) {
                        self.u[idx]
                    } else {
                        0.0
                    };
                }
            }
        }
    }

    /// Explicit advection+viscous sub-step; returns max change.
    pub(crate) fn substep(&mut self) -> f64 {
        let (nx, ny, nz) = self.domain.dims;
        let (dx, dy, dz) = self.domain.spacing;
        let dt = self.config.dt;
        let rho = self.config.density;
        let mut mu_eff = 0.0;
        let mut count = 0;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if self.domain.mask[self.domain.index(i, j, k)] {
                        mu_eff += self.mu[self.domain.index(i, j, k)];
                        count += 1;
                    }
                }
            }
        }
        mu_eff /= count.max(1) as f64;
        let nu = mu_eff / rho; // mm²/s

        let u0 = self.u.clone();
        let v0 = self.v.clone();
        let w0 = self.w.clone();

        // x-faces
        for i in 1..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if !self.ux_active(i, j, k) {
                        continue;
                    }
                    let idx = self.uid(i, j, k);
                    let mut dudt = 0.0;
                    // Persistent pressure gradient (cell-centred Π).
                    dudt -= (self.p[self.domain.index(i, j, k)]
                        - self.p[self.domain.index(i - 1, j, k)])
                        / dx;
                    dudt += nu
                        * (u0[self.uid(i + 1, j, k)] - 2.0 * u0[idx] + u0[self.uid(i - 1, j, k)])
                        / (dx * dx);
                    if j > 0 && j + 1 < ny {
                        dudt += nu
                            * (u0[self.uid(i, j + 1, k)] - 2.0 * u0[idx]
                                + u0[self.uid(i, j - 1, k)])
                            / (dy * dy);
                    }
                    if k > 0 && k + 1 < nz {
                        dudt += nu
                            * (u0[self.uid(i, j, k + 1)] - 2.0 * u0[idx]
                                + u0[self.uid(i, j, k - 1)])
                            / (dz * dz);
                    }
                    // Convection via cell-centered interpolation (skipped in
                    // creeping-flow configuration).
                    if self.config.include_convection {
                        let uc = |a: usize, b: usize, c: usize| {
                            0.5 * (u0[self.uid(a, b, c)] + u0[self.uid(a + 1, b, c)])
                        };
                        let uu = uc(i, j, k);
                        let duudx =
                            (u0[self.uid(i + 1, j, k)] - u0[self.uid(i - 1, j, k)]) / (2.0 * dx);
                        dudt -= uu * duudx;
                    }
                    self.u[idx] = u0[idx] + dt * dudt;
                }
            }
        }
        // y-faces: viscous Laplacian only (screening simplification).
        for i in 0..nx {
            for j in 1..ny {
                for k in 0..nz {
                    if !self.vy_active(i, j, k) {
                        continue;
                    }
                    let idx = i * (ny + 1) * nz + j * nz + k;
                    let mut dvdt = 0.0;
                    dvdt -= (self.p[self.domain.index(i, j, k)]
                        - self.p[self.domain.index(i, j - 1, k)])
                        / dy;
                    if i > 0 && i + 1 < nx {
                        dvdt += nu
                            * (v0[idx + (ny + 1) * nz] - 2.0 * v0[idx] + v0[idx - (ny + 1) * nz])
                            / (dx * dx);
                    }
                    if j < ny && j >= 1 {
                        dvdt += nu * (v0[idx + nz] - 2.0 * v0[idx] + v0[idx - nz]) / (dy * dy);
                    }
                    self.v[idx] = v0[idx] + dt * dvdt;
                }
            }
        }
        // z-faces.
        for i in 0..nx {
            for j in 0..ny {
                for k in 1..nz {
                    if !self.wz_active(i, j, k) {
                        continue;
                    }
                    let idx = i * ny * (nz + 1) + j * (nz + 1) + k;
                    let mut dwdt = 0.0;
                    dwdt -= (self.p[self.domain.index(i, j, k)]
                        - self.p[self.domain.index(i, j, k - 1)])
                        / dz;
                    if i > 0 && i + 1 < nx {
                        dwdt += nu
                            * (w0[idx + ny * (nz + 1)] - 2.0 * w0[idx] + w0[idx - ny * (nz + 1)])
                            / (dx * dx);
                    }
                    if k < nz && k >= 1 {
                        dwdt += nu * (w0[idx + 1] - 2.0 * w0[idx] + w0[idx - 1]) / (dz * dz);
                    }
                    self.w[idx] = w0[idx] + dt * dwdt;
                }
            }
        }
        let mut max_change = 0.0f64;
        for (a, b) in u0.iter().zip(&self.u) {
            max_change = max_change.max((a - b).abs());
        }
        max_change
    }

    /// Face divergence per fluid cell (1/s).
    pub(crate) fn divergence(&self) -> Vec<f64> {
        let (nx, ny, nz) = self.domain.dims;
        let (dx, dy, dz) = self.domain.spacing;
        let mut div = vec![0.0; nx * ny * nz];
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if !self.domain.mask[self.domain.index(i, j, k)] {
                        continue;
                    }
                    let du = (self.u[self.uid(i + 1, j, k)] - self.u[self.uid(i, j, k)]) / dx;
                    let dv = (self.v[i * (ny + 1) * nz + (j + 1) * nz + k]
                        - self.v[i * (ny + 1) * nz + j * nz + k])
                        / dy;
                    let dw = (self.w[i * ny * (nz + 1) + j * (nz + 1) + (k + 1)]
                        - self.w[i * ny * (nz + 1) + j * (nz + 1) + k])
                        / dz;
                    div[self.domain.index(i, j, k)] = du + dv + dw;
                }
            }
        }
        div
    }

    /// Pressure-Poisson solve for φ with Neumann walls and a Dirichlet
    /// φ = 0 outlet layer, dispatching on
    /// [`SolverConfig::pressure_solver`]. Returns `(φ, iterations_used)`.
    pub(crate) fn solve_poisson(&self, rhs: &[f64], outlet_anchor: f64) -> (Vec<f64>, usize) {
        match self.config.pressure_solver {
            PressureSolver::Sor => self.solve_poisson_sor(rhs, outlet_anchor),
            PressureSolver::ConjugateGradient => self.solve_poisson_cg(rhs, outlet_anchor),
        }
    }

    /// Jacobi-preconditioned conjugate gradient on the same discrete
    /// operator the SOR sweeps drive to their fixed point — the SPD
    /// negative Laplacian `A = −∇²` (denom·φᵢ − Σ w·φ_nb), matrix-free, so
    /// the mask needs no matrix assembly. Solid cells and the Dirichlet
    /// outlet layer carry identity rows with zero right-hand side, so they
    /// stay pinned at 0 (exactly what the SOR sweep enforces in place).
    ///
    /// Exits when the ℓ2 residual falls below `1e-6` of the right-hand
    /// side's norm, or at the `poisson_iterations` cap.
    fn solve_poisson_cg(&self, rhs: &[f64], outlet_anchor: f64) -> (Vec<f64>, usize) {
        let n = rhs.len();
        // The SOR sweep's fixed point is `denom·φ − Σ w·φ_nb = −rhs`, i.e.
        // A = −∇²; CG must solve the same system, so the right-hand side
        // enters negated.
        let mut b: Vec<f64> = rhs.iter().map(|&v| -v).collect();
        let (nx, _ny, _nz) = self.domain.dims;
        // The outlet fluid layer is Dirichlet at `outlet_anchor`: its row
        // is the identity, and a Dirichlet value is NOT part of the
        // negated-Laplacian convention (only ∇² equations are negated), so
        // the rhs carries the anchor as-is.
        for idx in 0..n {
            let (i, _j, _k) = self.domain.coords(idx);
            if i + 1 == nx {
                b[idx] = outlet_anchor;
            }
        }

        // Start at the anchor level: the Dirichlet constant is then exact
        // from iteration zero and CG only solves the (small) gauge-
        // relative part — starting from zero would spend the whole
        // iteration budget climbing a huge constant.
        let mut x = vec![outlet_anchor; n];
        let mut r = vec![0.0; n];
        self.laplacian_apply(&x, &mut r);
        for (ri, &bi) in r.iter_mut().zip(&b) {
            *ri = bi - *ri;
        }
        let norm_b = r.iter().map(|&v| v * v).sum::<f64>().sqrt();
        if norm_b == 0.0 {
            return (x, 0);
        }
        let exit = 1e-6 * norm_b;

        let z = self.jacobi_precondition(&r);
        let mut p = z.clone();
        let mut rz: f64 = r.iter().zip(&z).map(|(&a, &b)| a * b).sum();
        let mut ap = vec![0.0; n];
        let mut used = 0usize;
        for _ in 0..self.config.poisson_iterations {
            used += 1;
            self.laplacian_apply(&p, &mut ap);
            let pap: f64 = p.iter().zip(&ap).map(|(&a, &b)| a * b).sum();
            if pap <= 0.0 || !pap.is_finite() {
                break; // non-SPD breakdown; keep the current iterate
            }
            let alpha = rz / pap;
            for ((xi, pi), &_api) in x.iter_mut().zip(&p).zip(&ap) {
                *xi += alpha * pi;
            }
            let mut r_norm2 = 0.0f64;
            for ((ri, _pi), &api) in r.iter_mut().zip(&p).zip(&ap) {
                *ri -= alpha * api;
                r_norm2 += *ri * *ri;
            }
            if r_norm2.sqrt() < exit {
                break;
            }
            let z_next = self.jacobi_precondition(&r);
            let rz_next: f64 = r.iter().zip(&z_next).map(|(&a, &b)| a * b).sum();
            let beta = rz_next / rz;
            for idx in 0..n {
                p[idx] = z_next[idx] + beta * p[idx];
            }
            rz = rz_next;
        }
        (x, used)
    }

    /// Diagonal (Jacobi) preconditioner: interior fluid cells use their
    /// stencil diagonal (including Neumann mirror weights); pinned rows
    /// have diagonal 1.
    pub(crate) fn jacobi_precondition(&self, r: &[f64]) -> Vec<f64> {
        let (nx, ny, nz) = self.domain.dims;
        let (dx, dy, dz) = self.domain.spacing;
        let denom = 2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz));
        r.iter()
            .enumerate()
            .map(|(idx, &v)| {
                let (i, j, k) = self.domain.coords(idx);
                if !self.domain.mask[idx] || i + 1 == nx {
                    v
                } else {
                    let mut diag = denom;
                    if i > 0 && !self.domain.mask[self.domain.index(i - 1, j, k)] {
                        diag += 1.0 / (dx * dx);
                    }
                    if !self.domain.mask[self.domain.index(i + 1, j, k)] {
                        diag += 1.0 / (dx * dx);
                    }
                    if j > 0 && !self.domain.mask[self.domain.index(i, j - 1, k)] {
                        diag += 1.0 / (dy * dy);
                    }
                    if j + 1 < ny && !self.domain.mask[self.domain.index(i, j + 1, k)] {
                        diag += 1.0 / (dy * dy);
                    }
                    if k > 0 && !self.domain.mask[self.domain.index(i, j, k - 1)] {
                        diag += 1.0 / (dz * dz);
                    }
                    if k + 1 < nz && !self.domain.mask[self.domain.index(i, j, k + 1)] {
                        diag += 1.0 / (dz * dz);
                    }
                    v / diag
                }
            })
            .collect()
    }

    /// Matrix-free application of the masked Laplacian: `out = A·x` with
    /// Neumann mirrors at wall faces and an explicit 0 beyond the Dirichlet
    /// outlet layer.
    pub(crate) fn laplacian_apply(&self, x: &[f64], out: &mut [f64]) {
        let (nx, ny, nz) = self.domain.dims;
        let (dx, dy, dz) = self.domain.spacing;
        let denom = 2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz));
        for idx in 0..x.len() {
            let (i, j, k) = self.domain.coords(idx);
            if !self.domain.mask[idx] || i + 1 == nx {
                out[idx] = x[idx]; // identity rows on pinned unknowns
                continue;
            }
            let m = &self.domain.mask;
            let xm = if i > 0 && m[self.domain.index(i - 1, j, k)] {
                x[self.domain.index(i - 1, j, k)]
            } else {
                x[idx]
            };
            let xp = if m[self.domain.index(i + 1, j, k)] {
                x[self.domain.index(i + 1, j, k)]
            } else {
                x[idx]
            };
            let ym = if j > 0 && m[self.domain.index(i, j - 1, k)] {
                x[self.domain.index(i, j - 1, k)]
            } else {
                x[idx]
            };
            let yp = if j + 1 < ny && m[self.domain.index(i, j + 1, k)] {
                x[self.domain.index(i, j + 1, k)]
            } else {
                x[idx]
            };
            let zm = if k > 0 && m[self.domain.index(i, j, k - 1)] {
                x[self.domain.index(i, j, k - 1)]
            } else {
                x[idx]
            };
            let zp = if k + 1 < nz && m[self.domain.index(i, j, k + 1)] {
                x[self.domain.index(i, j, k + 1)]
            } else {
                x[idx]
            };
            out[idx] = denom * x[idx]
                - (xm + xp) / (dx * dx)
                - (ym + yp) / (dy * dy)
                - (zm + zp) / (dz * dz);
        }
    }

    /// SOR (in-place Gauss–Seidel with over-relaxation) Poisson solve for
    /// φ with Neumann walls and a Dirichlet φ = 0 outlet layer. SOR at
    /// ω = 1.9 converges an order of magnitude faster than Jacobi on this
    /// grid, which is required to keep the post-projection divergence
    /// near zero. Returns `(φ, sweeps_used)`.
    fn solve_poisson_sor(&self, rhs: &[f64], outlet_anchor: f64) -> (Vec<f64>, usize) {
        let (nx, ny, nz) = self.domain.dims;
        let (dx, dy, dz) = self.domain.spacing;
        // Start at the anchor level (see solve_poisson_cg): the Dirichlet
        // constant is exact from sweep zero.
        let mut phi = vec![outlet_anchor; nx * ny * nz];
        let denom = 2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz));
        let omega = 1.9;
        let rhs_scale = rhs.iter().fold(0.0f64, |m, &v| m.max(v.abs())).max(1e-10);
        // Natural φ magnitude is rhs/denom; exit when per-sweep updates are
        // a 1e-4 fraction of it (post-projection divergence ≪ gradients).
        let residual_exit = 1e-4 * rhs_scale / denom;

        let mut sweeps = 0usize;
        for _ in 0..self.config.poisson_iterations {
            sweeps += 1;
            let mut max_update = 0.0f64;
            for i in 0..nx {
                for j in 0..ny {
                    for k in 0..nz {
                        let idx = self.domain.index(i, j, k);
                        if !self.domain.mask[idx] {
                            continue;
                        }
                        if i + 1 == nx {
                            phi[idx] = outlet_anchor; // Dirichlet outlet
                            continue;
                        }
                        // Neumann walls: mirror the centre value across
                        // non-fluid neighbours; the outlet layer is the
                        // Dirichlet φ = 0 anchor.
                        let xm = if i > 0 && self.domain.mask[self.domain.index(i - 1, j, k)] {
                            phi[self.domain.index(i - 1, j, k)]
                        } else {
                            phi[idx]
                        };
                        let xp = if i + 1 == nx {
                            0.0
                        } else if self.domain.mask[self.domain.index(i + 1, j, k)] {
                            phi[self.domain.index(i + 1, j, k)]
                        } else {
                            phi[idx]
                        };
                        let ym = if j > 0 && self.domain.mask[self.domain.index(i, j - 1, k)] {
                            phi[self.domain.index(i, j - 1, k)]
                        } else {
                            phi[idx]
                        };
                        let yp = if j + 1 < ny && self.domain.mask[self.domain.index(i, j + 1, k)] {
                            phi[self.domain.index(i, j + 1, k)]
                        } else {
                            phi[idx]
                        };
                        let zm = if k > 0 && self.domain.mask[self.domain.index(i, j, k - 1)] {
                            phi[self.domain.index(i, j, k - 1)]
                        } else {
                            phi[idx]
                        };
                        let zp = if k + 1 < nz && self.domain.mask[self.domain.index(i, j, k + 1)] {
                            phi[self.domain.index(i, j, k + 1)]
                        } else {
                            phi[idx]
                        };
                        let sum = xm / (dx * dx)
                            + xp / (dx * dx)
                            + ym / (dy * dy)
                            + yp / (dy * dy)
                            + zm / (dz * dz)
                            + zp / (dz * dz);
                        let target = (sum - rhs[idx]) / denom;
                        let update = omega * (target - phi[idx]);
                        phi[idx] += update;
                        max_update = max_update.max(update.abs());
                    }
                }
            }
            if max_update < residual_exit {
                break;
            }
        }
        (phi, sweeps)
    }

    /// One full step; returns max |velocity change| (mm/s).
    pub fn step(&mut self) -> f64 {
        self.apply_boundary();
        let change = self.substep();
        self.apply_boundary();

        let div = self.divergence();
        let dt = self.config.dt;
        let rhs: Vec<f64> = div.iter().map(|&d| d / dt).collect();
        let (phi, iterations) = self.solve_poisson(&rhs, self.outlet_anchor);
        self.last_pressure_solve_iterations = iterations;

        let (nx, ny, nz) = self.domain.dims;
        let (dx, dy, dz) = self.domain.spacing;
        // Pressure correction on faces.
        for i in 1..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if self.ux_active(i, j, k) {
                        let idx = self.uid(i, j, k);
                        self.u[idx] -= dt
                            * (phi[self.domain.index(i, j, k)]
                                - phi[self.domain.index(i - 1, j, k)])
                            / dx;
                    }
                }
            }
        }
        for i in 0..nx {
            for j in 1..ny {
                for k in 0..nz {
                    if self.vy_active(i, j, k) {
                        let idx = i * (ny + 1) * nz + j * nz + k;
                        self.v[idx] -= dt
                            * (phi[self.domain.index(i, j, k)]
                                - phi[self.domain.index(i, j - 1, k)])
                            / dy;
                    }
                }
            }
        }
        for i in 0..nx {
            for j in 0..ny {
                for k in 1..nz {
                    if self.wz_active(i, j, k) {
                        let idx = i * ny * (nz + 1) + j * (nz + 1) + k;
                        self.w[idx] -= dt
                            * (phi[self.domain.index(i, j, k)]
                                - phi[self.domain.index(i, j, k - 1)])
                            / dz;
                    }
                }
            }
        }
        // Zero-gradient outlet re-applied POST-correction so the outlet
        // face participates in mass balance (it has no adjacent correction
        // cell).
        let nx_out = nx;
        for j in 0..ny {
            for k in 0..nz {
                let idx_out = self.uid(nx_out, j, k);
                let idx_in = self.uid(nx_out - 1, j, k);
                self.u[idx_out] = self.u[idx_in];
            }
        }

        // Pressure accumulation, gauge-relative: with a Windkessel anchor
        // the Dirichlet level would otherwise be integrated once per step,
        // so the stored field carries only the correction relative to the
        // outlet reference (the absolute level is the boundary model's).
        for (pp, &ph) in self.p.iter_mut().zip(&phi) {
            *pp += ph - self.outlet_anchor;
        }
        // Non-Newtonian viscosity update.
        if !matches!(self.blood, BloodModel::Newtonian { .. }) {
            self.update_viscosity();
        }
        change
    }

    fn update_viscosity(&mut self) {
        let (nx, ny, nz) = self.domain.dims;
        let (dx, dy, dz) = self.domain.spacing;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = self.domain.index(i, j, k);
                    if !self.domain.mask[idx] {
                        continue;
                    }
                    // Velocity gradients from face values (central).
                    let dudx = (self.u[self.uid(i + 1, j, k)] - self.u[self.uid(i, j, k)]) / dx;
                    let dvdy = (self.v[i * (ny + 1) * nz + (j + 1) * nz + k]
                        - self.v[i * (ny + 1) * nz + j * nz + k])
                        / dy;
                    let dwdz = (self.w[i * ny * (nz + 1) + j * (nz + 1) + (k + 1)]
                        - self.w[i * ny * (nz + 1) + j * (nz + 1) + k])
                        / dz;
                    let dudy = 0.5
                        * (self.u[self.uid(i, core::cmp::min(j + 1, ny - 1), k)]
                            - self.u[self.uid(i, j.saturating_sub(1), k)])
                        / dy;
                    let g = (dudx * dudx
                        + dvdy * dvdy
                        + dwdz * dwdz
                        + 0.5 * (dudy * dudy + dvdy * dvdy))
                        .max(0.0)
                        .sqrt()
                        + 1e-6;
                    let mu_target = self.blood.viscosity(g);
                    let relax = self.config.viscosity_relaxation;
                    self.mu[idx] = self.mu[idx] * (1.0 - relax) + mu_target * relax;
                }
            }
        }
    }

    /// One full step with a **Windkessel-driven outlet** (`tpt-med-
    /// cardiovascular`'s `CoupledWindkessel`): the projection's Dirichlet
    /// anchor is set to the boundary model's current pressure — converted
    /// to the solver's internal mm²/s² units via `p·10⁶/ρ` — and after
    /// the step the measured outlet flow advances the 0-D model. Explicit
    /// staggered coupling: the imposed absolute pressure level acts
    /// through the pressure field, the flow feedback lags one step (the
    /// same stability posture as `CoupledWindkessel` itself; check
    /// `WindkesselModel::time_constant` against `dt` at the call site).
    ///
    /// After a coupled run, the reported pressures are **gauge-relative
    /// to the boundary model's outlet pressure**: the absolute outlet
    /// pressure *is* `wk.pressure()`, and absolute values elsewhere are
    /// the reported field plus that reference (converted with
    /// `p_mpa = Π·ρ/10⁶`).
    pub fn step_coupled(&mut self, wk: &mut tpt_med_cardiovascular::CoupledWindkessel) -> f64 {
        self.outlet_anchor = wk.pressure() * 1.0e6 / self.config.density;
        let change = self.step();
        let flow = self.outlet_flow();
        wk.advance(flow, self.config.dt);
        change
    }

    /// Face-integrated flow through the outlet plane (mm³/s).
    fn outlet_flow(&self) -> f64 {
        let (nx, ny, nz) = self.domain.dims;
        let (_dx, dy, dz) = self.domain.spacing;
        let cell_face = dy * dz;
        let mut outlet = 0.0;
        for j in 0..ny {
            for k in 0..nz {
                if self.domain.is_fluid((nx - 1) as i64, j as i64, k as i64) {
                    outlet += self.u[self.uid(nx, j, k)] * cell_face;
                }
            }
        }
        outlet
    }

    /// Marches to steady state; stops when the velocity change per step is
    /// below `tolerance` (mm/s) or `max_steps` is exhausted.
    pub fn run_steady(&mut self, max_steps: usize, tolerance: f64) -> SteadyStats {
        let mut steps = 0;
        for _ in 0..max_steps {
            let change = self.step();
            steps += 1;
            if change < tolerance {
                break;
            }
        }
        self.stats(steps)
    }

    /// Face-integrated flows and pressures.
    pub fn stats(&self, steps: usize) -> SteadyStats {
        let (nx, ny, nz) = self.domain.dims;
        let (_dx, dy, dz) = self.domain.spacing;
        let cell_face = dy * dz;
        let mut inlet = 0.0;
        for j in 0..ny {
            for k in 0..nz {
                if self.domain.is_fluid(0, j as i64, k as i64) {
                    inlet += self.u[self.uid(0, j, k)] * cell_face;
                }
            }
        }
        let mut outlet = 0.0;
        for j in 0..ny {
            for k in 0..nz {
                if self.domain.is_fluid((nx - 1) as i64, j as i64, k as i64) {
                    outlet += self.u[self.uid(nx, j, k)] * cell_face;
                }
            }
        }
        let mut max_v = 0.0f64;
        for &uv in &self.u {
            max_v = max_v.max(uv.abs());
        }
        let mut pin = 0.0;
        let mut nin = 0;
        for j in 0..ny {
            for k in 0..nz {
                if self.domain.is_fluid(0, j as i64, k as i64) {
                    pin += self.p[self.domain.index(0, j, k)];
                    nin += 1;
                }
            }
        }
        let mut pout = 0.0;
        let mut nout = 0;
        for j in 0..ny {
            for k in 0..nz {
                if self.domain.is_fluid((nx - 1) as i64, j as i64, k as i64) {
                    pout += self.p[self.domain.index(nx - 1, j, k)];
                    nout += 1;
                }
            }
        }
        SteadyStats {
            steps,
            max_velocity: max_v,
            inlet_flow: inlet,
            outlet_flow: outlet,
            mean_pressure_inlet: pin / nin.max(1) as f64,
            mean_pressure_outlet: pout / nout.max(1) as f64,
            pressure_drop: pin / nin.max(1) as f64 - pout / nout.max(1) as f64,
        }
    }

    /// Cell-centered axial velocity at plane `i` (mm/s) — verification
    /// helper.
    pub fn axial_velocity_profile(&self, i: usize) -> Vec<f64> {
        let (_, ny, nz) = self.domain.dims;
        let mut out = Vec::with_capacity(ny * nz);
        for j in 0..ny {
            for k in 0..nz {
                let fluid = self.domain.is_fluid(i as i64, j as i64, k as i64);
                if fluid {
                    let val = 0.5 * (self.u[self.uid(i, j, k)] + self.u[self.uid(i + 1, j, k)]);
                    out.push(val);
                }
            }
        }
        out
    }
}
