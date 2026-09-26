# tpt-med-hemodynamics

Incompressible Navier–Stokes CFD for blood flow on voxel domains, with
wall shear stress (WSS) and oscillatory shear index (OSI) post-processing.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--hemodynamics-orange)](https://crates.io/crates/tpt-med-hemodynamics)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--hemodynamics-blue)](https://docs.rs/tpt-med-hemodynamics)

| | |
|---|---|
| **Layer** | `fluid` |
| **Status** | Alpha, `0.1.0` |
| **Class** | Screening-grade, laminar only |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-geometry`](../../core/tpt-med-geometry), [`tpt-med-meshing`](../../imaging/tpt-med-meshing), [`tpt-med-units`](../../core/tpt-med-units) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Arterial screening questions — *is this plaque region high-shear? is this
stent underexpanded where the OSI is high?* — are answered by computing the
velocity field and looking at the wall. Doing that with a commercial CFD suite
means meshing a patient-specific lumen, a licence, and a cloud upload of the
patient's CT. This crate does it on the voxel grid you already have.

**Be clear about what this is.** It is a *screening-grade* solver:

- Staggered-grid (MAC) projection method: explicit advection and viscous
  sub-step, then an SOR (successive over-relaxation) pressure-Poisson
  projection enforcing incompressibility inside a voxel fluid mask.
- Stair-step walls — no body-fitted mesh, no boundary-layer refinement.
- **Laminar only.** No turbulence model. This is the appropriate regime for
  most arterial screening (Re < 2000) and is *not* appropriate for aortic
  aneurysm flow or any turbulent case.
- Unstructured, high-fidelity CFD is the documented upgrade path via
  `tpt-sci-cfd-core` / `tpt-sci-hemodynamics`, pinned in the workspace
  manifest.

The value is not that it replaces Fluent. It is that a surgeon can run it, in a
browser, in seconds, without a licence.

## Features

- **Projection-method Navier–Stokes** on a MAC grid: convection (optional),
  viscous diffusion, and an SOR (Gauss–Seidel with over-relaxation, ω = 1.9)
  pressure Poisson projection — an order of magnitude faster to converge than
  plain Jacobi on this grid.
- **Blood rheology** — `Newtonian`, `CarreauYasuda` (shear-thinning) and
  `Casson` (yield-stress), with `viscosity(shear_rate)` and under-relaxation
  for stability.
- **Voxel domains** — `FluidDomain::from_mask` for real lumen geometry taken
  from a segmentation, and `FluidDomain::cylinder` as the canonical
  verification domain.
- **Boundary conditions** — inlet plug velocity, Dirichlet pressure anchor at
  the outlet, Neumann (no-slip) walls everywhere else.
- **WSS extraction** — `extract_wss` computes the tangential traction on
  wall-adjacent faces with a linear wall gradient over half a cell, plus
  `mean_magnitude` / `max_magnitude`.
- **Oscillatory shear index** — `OsiAccumulator` integrates over a pulsatile
  cycle, which is the standard atheroma-relevant metric.
- **Diagnostics** — `SteadyStats` reports inlet vs. outlet flow (mass
  conservation), pressure drop, and max velocity; `axial_velocity_profile`
  exposes the centreline profile for direct comparison with theory.
- `PLASMA` constant (0.0012 Pa·s) for reference.

## Conventions

- Lengths **mm**, velocity **mm/s**, time **s**.
- Viscosity in **Pa·s**, stress in **Pa**, density in **g/mm³**
  (1.06e-3 g/mm³ = 1060 kg/m³).
- Cell mask index: `(i*ny + j)*nz + k` (note the ordering — it is *not* the
  meshing crate's `(z·ny+y)·nx+x`).
- `SolverConfig::default()`: `dt = 2.0e-4 s`, `poisson_iterations = 400`,
  `include_convection = false`, `viscosity_relaxation = 0.2`,
  `density = 1.06e-3`.
  Convection is **off by default**: for creeping arterial flow the term is
  small and disabling it converges faster and more robustly.

## Usage

```rust
use tpt_med_hemodynamics::{
    extract_wss, BloodModel, FluidDomain, HemodynamicsSolver, OsiAccumulator, SolverConfig,
};

fn main() {
    // A 4.5 mm-radius tube, 0.5 mm cells, 24 cells long along x.
    let domain = FluidDomain::cylinder(24, 12, 4.5, 0.5, 0);
    let (n_fluid, _) = domain.fluid_stats();

    let mut solver = HemodynamicsSolver::new(
        domain,
        BloodModel::CARREAU_YASUDA_BLOOD,
        30.0, // inlet plug velocity, mm/s
        SolverConfig::default(),
    );
    solver.apply_boundary();
    let stats = solver.run_steady(500, 1.0e-4);

    // Mass conservation is the first thing to check.
    assert!((stats.inlet_flow - stats.outlet_flow).abs() / stats.inlet_flow.max(1e-9) < 0.05);
    println!("{} fluid cells, max u = {:.2} mm/s", n_fluid, stats.max_velocity);
    println!("dp = {:.4}", stats.pressure_drop);

    // Wall shear stress on the lumen wall.
    let wss = extract_wss(&solver);
    println!("WSS: mean {:.2} Pa, max {:.2} Pa",
             wss.mean_magnitude(), wss.max_magnitude());

    // Accumulate the oscillatory shear index over a pulsatile cycle.
    let mut osi = OsiAccumulator::default();
    for _ in 0..200 {
        let traction = wss.tractions.first().copied().unwrap_or_default();
        osi.sample(traction, 0.005);
    }
    let index = osi.osi();
    assert!((0.0..=0.5).contains(&index));

    // Centreline profile, for direct comparison with Poiseuille theory.
    let profile = solver.axial_velocity_profile(12);
    assert!(profile[0] < profile[profile.len() / 2]);
}
```

## API Overview

| Item | Purpose |
|---|---|
| `BloodModel` | `Newtonian { viscosity }`, `CarreauYasuda { .. }`, `Casson { .. }` |
| `BloodModel::viscosity(shear_rate) -> f64` | Effective viscosity at a local shear rate (Pa·s) |
| `BloodModel::NEWTONIAN_BLOOD` | 0.0035 Pa·s — normal hematocrit |
| `BloodModel::CARREAU_YASUDA_BLOOD` | Shear-thinning, typical vessel-wall parameters |
| `PLASMA: f64` | 0.0012 Pa·s |
| `FluidDomain` | `dims`, `spacing`, `mask`, `flow_axis`, `inlet_low` |
| `FluidDomain::from_mask(dims, spacing, mask, flow_axis, inlet_low)` | Real lumen geometry from a segmentation |
| `FluidDomain::cylinder(n_axial, n_radius, radius_cells, spacing, flow_axis)` | Canonical verification geometry |
| `FluidDomain::{index, is_fluid, fluid_stats}` | Mask accessors; `(fluid_count, volume)` |
| `SolverConfig` | `dt`, `poisson_iterations`, `include_convection`, `viscosity_relaxation`, `density` |
| `HemodynamicsSolver::new(domain, blood, inlet_velocity, config)` | Construct over a domain |
| `HemodynamicsSolver::apply_boundary()` | Impose inlet velocity, outlet pressure, no-slip walls |
| `HemodynamicsSolver::step() -> f64` | One time step; returns the max velocity change |
| `HemodynamicsSolver::run_steady(max_steps, tolerance) -> SteadyStats` | March to steady state |
| `HemodynamicsSolver::stats(steps) -> SteadyStats` | Current scalars without marching |
| `HemodynamicsSolver::axial_velocity_profile(i) -> Vec<f64>` | Velocity along an axial station — for theory comparison |
| `SteadyStats` | `steps`, `max_velocity`, `inlet_flow`, `outlet_flow`, `mean_pressure_inlet`, `mean_pressure_outlet`, `pressure_drop` |
| `WssField` | `positions: Vec<Vec3>`, `tractions: Vec<Vec3>` (Pa) |
| `WssField::{mean_magnitude, max_magnitude}` | WSS statistics |
| `extract_wss(&HemodynamicsSolver) -> WssField` | Wall traction on wall-adjacent faces |
| `OsiAccumulator` | `sample(traction, dt)` per step; `osi()` returns the index |
| `tpt_med_geometry::Vec3` | Wall sample positions and traction vectors |

## Verification

This crate has some of the strongest verification in the workspace, because
pipe flow has exact solutions:

- **Poiseuille profile** — the centreline velocity profile from
  `axial_velocity_profile` is compared against the analytic paraboloid
  `u(r) = (Δp/4μL)(R² − r²)`; both curvature and peak value are asserted.
- **Tube wall-shear law** — mean WSS is compared against `τ = 4μQ/(πR³)`
  within a stated tolerance.
- **Mass conservation** — inlet and outlet flow rates must agree after a
  steady run, which doubles as the incompressibility check.
- **Shear-thinning ordering** — the Carreau–Yasuda solution has *lower*
  resistance than Newtonian at the same flow rate, asserted as a strict
  inequality. A reversed comparison here would expose a sign error that a
  magnitude-only check would miss.
- **Reynolds insensitivity** — the laminar solution is invariant under
  Reynolds scaling, asserted to catch an accidental unit error.
- Golden reference datasets:
  `test-data/golden/fluid/carotid_bifurcation_cfd.json` and
  `aortic_aneurysm_flow.json`.

## Known Limitations

- **Laminar only.** No turbulence model. Do not use for aortic aneurysm or
  high-Reynolds geometry, and do not report a laminar result as if it were
  turbulent.
- Stair-step walls; no body-fitted mesh and no wall refinement. Peak WSS at a
  geometric corner is therefore resolution-dependent.
- No conjugate heat transfer, no wall compliance, no FSI coupling. Vessel
  compliance is handled on the boundary side by
  [`tpt-med-cardiovascular`](../tpt-med-cardiovascular/README.md).
- 2D-in-3D is supported (a `cylinder` with `flow_axis = 2`), but genuinely 2D
  cases should be treated as 3D slabs.

## Related Crates

- [`tpt-med-cardiovascular`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/fluid/tpt-med-cardiovascular) — Windkessel boundary conditions, FFR, and the pulsatile waveforms that drive this solver.
- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) and [`tpt-med-meshing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-meshing) — supply the lumen segmentation mask.
- [`tpt-med-stents`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-stents) — stent deployment, the main downstream use case.
- [`tpt-med-geometry`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-geometry) — `Vec3` for WSS positions and tractions.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Numerical code
must ship verification against an analytical or published reference. Cite the
source of every rheological constant.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. Not a
diagnostic device; hemodynamics output is for research screening.
