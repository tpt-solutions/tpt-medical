# Changelog

All notable changes to `tpt-med-hemodynamics` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- **`MembraneWall` / `step_coupled` / `MembraneWallState` (`fsi`) — wall
  compliance, the reduced fluid–structure coupling.** The vessel wall is
  an axisymmetric membrane: each axial section's radius relaxes
  first-order toward the pressure-set equilibrium
  `r_eq = r₀(1 + C·(p − p_ext))` (caller-cited fractional compliance and
  relaxation time), and the domain mask is rebuilt from the radii every
  step — the fluid feels the wall move, the wall feels the fluid's
  pressure. Section transmural pressure is read from the accumulated
  projection potential in the exact inverse of the Windkessel anchor's
  write convention (`p[MPa] = ΔΠ·ρ/1e6`), gauge-relative to the outlet
  layer, so only pressure differences drive the wall. Scope stated at
  the API: a fixed grid with a moving stair-stepped mask (no ALE, no
  immersed boundary), newly-wetted cells start at rest, and the
  relaxation time should keep per-step wall motion well inside a cell.
  Verified: a zero-compliance wall reproduces the fixed-domain run
  **bit for bit** (flow fields, pressure potential and mask identical);
  a compliant wall settles onto the pressure-set equilibrium radius
  — re-derived independently from the final flow state — section by
  section; the mask tracks the radii exactly; and the wall parameters
  are validated.
- `step_coupled(&mut CoupledWindkessel)`: **coupling to
  `tpt-med-cardiovascular` for a driven outlet boundary** — the
  projection's Dirichlet anchor is set to the boundary model's current
  pressure (converted `p·10⁶/ρ` into the solver's mm²/s² units) and the
  measured outlet flow advances the 0-D model each step (explicit
  staggered coupling; the usual `time_constant`-vs-`dt` check applies).
  Both Poisson solvers now start from the anchor level — the constant
  mode is exact from iteration zero, which a large anchor otherwise
  turns into a spurious outlet gradient — and the stored pressure field
  accumulates only the gauge-relative correction (reported pressures
  become relative to the boundary model's outlet pressure; the absolute
  level is `wk.pressure()`). Verified by an exact replay assertion: the
  flow history pushed through the boundary model's own RK4 reproduces
  the advanced state to machine precision, plus flow conservation.
  `tpt-med-cardiovascular` becomes a dependency (acyclic; it depends on
  nothing here).
- Crate README stating the screening-grade scope up front — laminar only, no
  body-fitted mesh, no turbulence model — and directing high-fidelity users to
  `tpt-sci-cfd-core` / `tpt-sci-hemodynamics`.

- `PressureSolver::{Sor, ConjugateGradient}` + `SolverConfig::pressure_solver`
  (+ `HemodynamicsSolver::last_pressure_solve_iterations`): the pressure
  projection can now run **Jacobi-preconditioned conjugate gradient** —
  matrix-free on the masked grid, pinned unknowns (solid cells, Dirichlet
  outlet layer) as identity rows — as an alternative to the default SOR
  sweeps. Verified against a manufactured solution of the discrete
  operator (which would catch a sign error: the SOR fixed point is the SPD
  negative Laplacian, so the rhs enters CG negated), and by the full
  Poiseuille march under CG (developed, concave, symmetric profile).
  On the verification tube CG settles in ~40 iterations per projection
  against SOR's 400-sweep cap. **The default remains SOR**, so existing
  golden datasets and benchmark baselines are unchanged; switching a
  production run to CG is a one-line config change and a golden re-run.

### Added
- **`heat` module: passive scalar / temperature transport on the MAC
  grid** — the fluid side of the conjugate-heat-transfer item.
  Conservative flux-form upwind advection plus centered 7-point
  diffusion with thermal diffusivity κ (mm²/s; blood ≈ 0.12), wall
  treatment per `ScalarWall` (`Fixed(t)`: isothermal wall with the
  half-cell flux; `Insulated`: adiabatic AND impermeable — no advective
  flux crosses it, which is what makes the insulated conservation test
  exact). `stable_time_step` bounds the explicit step by the stricter of
  the advective and three-dimensional diffusive limits. Verification:
  the insulated and fixed-wall discrete Fourier eigenmodes decay by
  exactly the discrete amplification factor (the initial conditions are
  the operator's exact eigenvectors — cell-centered Neumann
  `cos(π(i+½)/N)` and half-cell Dirichlet `sin(π(i+½)/N)` product modes);
  flux-form advection on a stream-function-built divergence-free field
  conserves the total scalar to machine precision; the thermal-entry
  scenario in a developed Poiseuille tube decays monotonically toward
  cold walls and retains advected heat under insulated walls (the exact
  Nusselt number is deliberately not asserted on a stair-step coarse
  grid). Remaining for full CHT: the solid-side conduction solver and a
  two-way coupled interface.
- **`conjugate_step`: the solid side, completing CHT** — diffusion
  through *both* regions in one step, with interface faces taking the
  **harmonic-mean** conductivity `2 κ_f κ_s/(κ_f + κ_s)`: for
  cell-centered finite volumes that makes the interface face's
  resistance the exact series sum of the two half-cell resistances, so
  the two-layer steady state IS the analytic composite-wall solution.
  Verified to machine precision across all 32 cells (fluid and solid)
  against the resistance-chain potentials, including the two-layer
  interface-temperature formula; two-region conservation holds exactly
  through the interface exchange, and the regions equilibrate. The
  fluid's wall temperature is no longer prescribed — it *is* the solid
  cell's temperature, and the interface flux is continuous. `ScalarWall`
  gained a `Faces([...])` variant (per-face fixed values, `None` =
  insulated) so 1-D analytic tests are expressible on a 3-D box.

### Planned
- Optional local wall refinement, so peak WSS at a geometric corner stops
  being resolution dependent. (Wall compliance is delivered — the reduced
  membrane coupling in `fsi`; a full FSI/ALE solve remains the substrate
  upgrade path.)


- Multigrid pressure solve (a level beyond the new CG option) if CG's
  √-condition-number scaling is ever insufficient on production grids.

### Notes
- **This crate is laminar only and must not be reported as if it were
  turbulent.** Adding a turbulence model is an RFC-scale change with its own
  verification programme.
- `SolverConfig::default()` has `include_convection = false`, because for
  creeping arterial flow the convective term is small and disabling it
  converges faster and more robustly. Turning it on changes the answer.
- The cell mask index is `(i*ny + j)*nz + k`, which is **not** the
  `(z·ny+y)·nx+x` ordering used by the imaging crates. Mixing them up produces
  a transposed domain that still runs.

## [0.1.0] - 2026-09-22

### Added
- **Projection-method Navier–Stokes** on a staggered MAC grid: explicit
  convection (optional) and viscous sub-step, then a pressure-Poisson
  projection enforcing incompressibility inside the voxel fluid mask.
- **Blood rheology** — `BloodModel::Newtonian`,
  `BloodModel::CarreauYasuda` (shear-thinning) and `BloodModel::Casson`
  (yield stress), with `viscosity(shear_rate)` and
  `NEWTONIAN_BLOOD` (0.0035 Pa·s) / `CARREAU_YASUDA_BLOOD` presets, plus
  under-relaxation (`viscosity_relaxation`) for stability.
- `PLASMA: f64` = 0.0012 Pa·s.
- **Voxel domains** — `FluidDomain::from_mask` for real lumen geometry taken
  from a segmentation, and `FluidDomain::cylinder` as the canonical
  verification domain. Accessors `index`, `is_fluid` and `fluid_stats`.
- `SolverConfig` — `dt` (2.0e-4 s), `poisson_iterations` (400),
  `include_convection` (false), `viscosity_relaxation` (0.2) and `density`
  (1.06e-3 g/mm³).
- `HemodynamicsSolver::new`, `::apply_boundary` (inlet plug velocity, Dirichlet
  pressure anchor at the outlet, no-slip walls elsewhere), `::step`,
  `::run_steady(max_steps, tolerance)`, `::stats` and
  `::axial_velocity_profile` for direct comparison with theory.
- `SteadyStats` — `steps`, `max_velocity`, `inlet_flow`, `outlet_flow`,
  `mean_pressure_inlet`, `mean_pressure_outlet`, `pressure_drop`.
- **Post-processing** — `extract_wss(&HemodynamicsSolver) -> WssField` with
  `positions` and `tractions` (Pa) and `mean_magnitude`/`max_magnitude`;
  `OsiAccumulator` with `sample(traction, dt)` and `osi()` for the oscillatory
  shear index over a pulsatile cycle.

### Verification
Pipe flow has exact solutions, so the checks are analytic:
- **Poiseuille profile** — `axial_velocity_profile` compared against the
  analytic paraboloid `u(r) = (Δp/4μL)(R² − r²)`; both curvature and peak value
  asserted.
- **Tube wall-shear law** — mean WSS compared against `τ = 4μQ/(πR³)`.
- **Mass conservation** — inlet and outlet flow agree after a steady run, which
  doubles as the incompressibility check.
- **Shear-thinning ordering** — the Carreau–Yasuda solution has strictly lower
  resistance than Newtonian at the same flow rate. A reversed comparison would
  expose a sign error that a magnitude-only check would miss.
- **Reynolds insensitivity** — the laminar solution is invariant under
  Reynolds scaling, catching an accidental unit error.
- Golden datasets `test-data/golden/fluid/carotid_bifurcation_cfd.json` and
  `aortic_aneurysm_flow.json`.

### Known limitations
- **Laminar only** — no turbulence model. Not for aortic aneurysm or
  high-Reynolds geometry.
- Stair-step walls; no body-fitted mesh or wall refinement, so peak WSS at a
  corner is resolution dependent.
- No conjugate heat transfer, no wall compliance, no FSI coupling.
- 2D-in-3D is supported; genuinely 2D cases should be treated as 3D slabs.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
