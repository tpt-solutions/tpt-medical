# Changelog

All notable changes to `tpt-med-cardiovascular` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- `CoupledWindkessel`: a stateful boundary condition for lockstep coupling
  with a `tpt-med-hemodynamics` time step — the CFD outlet flow advances
  the 0-D model (RK4) and the returned pressure feeds back as the next
  step's outlet value. Explicit coupling; the stability ratio against
  `time_constant` is the caller's check.
- Crate README stating the lumped-model limitation explicitly: a Windkessel
  boundary condition reproduces the *global* impedance of the distal tree, not
  the shape or timing of a reflected wave, so it is the wrong tool for
  pulse-wave-velocity or augmentation-index work.

- `FourElementWindkessel`: the **four-element Windkessel** — the 3-element
  model plus an inertance `L` in the series branch (Stergiopoulos, Young &
  Westerhof 1999) — in its pressure-driven form (state `(p, Q)`, prescribed
  inlet pressure, RK4), which is the formulation where the inertial element
  produces genuinely new (second-order, ring-down) dynamics. Steady state is
  exact and independent of `L`. Four new tests: steady state (including
  `L`-independence), underdamped ring-down against the theoretical envelope
  decay, DC-gain/superposition of the whole transient, and the exact
  `(p, Q)` steady values.
- `WaterfallResistor`: the **vascular waterfall** (Starling-resistor)
  non-linear pressure–flow relation — flow is independent of downstream
  pressure once the vessel collapses below a critical closing pressure
  (Permutt & Bromberger-Barnea), the standard non-linear element of systemic
  and cerebral circulation modelling.

### Planned
- Waveform-based instantaneous-hyperbolic FFR, alongside the pressure-ratio
  definition implemented here.
- Patient-specific waveform fitting, rather than the fixed analytic shapes.

### Notes
- The FFR ischaemia threshold is `≤ 0.80`, and it is **inclusive**. The
  boundary is asserted explicitly, because an off-by-one here is a clinical
  misclassification.

## [0.1.0] - 2026-09-22

### Added
- `WindkesselModel { r_c, r_p, c, p_out }` — the 3-element Westerhof model
  obeying `C dp/dt = (1 + Rc/Rp)·Q − (p − p_out)/Rp`, with `two_element`
  implemented as the `Rc = 0` case rather than a separate code path, so the
  two cannot drift apart.
- Analytic quantities: `dp_dt(p, flow)`, `steady_state_pressure(flow)` and
  `time_constant()`, so a boundary condition can be sanity-checked without
  integrating anything.
- **Two integrators** — `step_rk4` (explicit RK4, the accuracy reference) and
  `simulate(p0, dt, steps, flow_fn)` (the semi-implicit scheme used online by
  the CFD boundary coupling; unconditionally stable, not the most accurate).
- `FractionalFlowReserve::calculate(p_distal, p_aortic)` and
  `::is_ischemic(ffr)`, using the accepted `≤ 0.80` cut-off for
  physiologically significant stenosis.
- `FlowWaveform::flow(t)` with `carotid_default()` and `coronary_default()`,
  which drive pulsatile CFD runs.
- Typed units throughout: resistances in MPa·s/mm³, compliance in mm³/MPa,
  pressures in MPa, flow in mm³/s.
- No dependencies; `#![forbid(unsafe_code)]`.

### Verification
- **Fixed point** — `dp_dt` is zero at `steady_state_pressure(flow)`, across a
  sweep of flow and both 2- and 3-element configurations.
- **Integrator vs. closed form** — the analytic `time_constant()` is compared
  against the time for the RK4 solution to decay back to within 1/e of
  equilibrium, so the integrator and the theory must agree rather than merely
  being self-consistent.
- **Two-element limit** — `Rc = 0` reproduces `two_element` exactly, element
  by element.
- **FFR boundary** — the classifier is asserted *at* 0.80 (ischaemic; the
  threshold is inclusive), just above and just below.
- **FFR degeneracies** — zero aortic pressure yields neither `NaN` nor a
  panic, and FFR > 1 is reported rather than silently clamped.
- **Waveform sanity** — both waveforms are strictly positive, bounded, return
  to a diastolic baseline within a cycle, and satisfy `flow(0) == flow(cycle)`
  so they are continuous.
- Golden dataset `test-data/golden/fluid/coronary_ffr.json`.

### Known limitations
- Lumped, not distributed: no wave-reflection shape or timing, no
  pulse-wave-velocity or augmentation-index work.
- Two- and three-element only.
- FFR is the pressure-ratio definition; instantaneous-hyperbolic FFR is not
  implemented.
- Waveforms are fixed analytic shapes, not patient-specific measurements.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
