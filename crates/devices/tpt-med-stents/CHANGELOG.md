# Changelog

All notable changes to `tpt-med-stents` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README naming the fidelity ladder from RFC 0004 — Level 1 ring model
  (shipped) → tapered ring groups → 3D superelastic FEM with frictional contact
  on `tpt-fem` — and stating that the default parameters are literature-typical
  starting points, **not vendor data**.
- `StentModel::foreshortening(manufactured_length, diameter, link_fraction)`:
  **geometric foreshortening** — the axial shortening when a zig-zag crown
  ring opens from its crimped configuration, from fixed strut segment
  lengths in a diamond-cell model. `link_fraction` (straight axial links,
  ≈0.6–0.8 for real laser-cut designs) calibrates the pure-geometry upper
  bound to the published few-percent band; `NaN` beyond the
  developed-length limit where the cell cannot close.
- `simulate_deployment_with_crowns` + `NonUniformDeployment`: **per-crown
  stiffness variation** for a non-uniform ring — ring-level metrics from the
  summed stiffness, with the radial force split proportionally to per-crown
  stiffness and the largest single-crown share reported (an empty slice
  falls back to the uniform model exactly).
- `simulate_tapered_deployment` / `TaperedDeployment`: **Level 2 of the
  RFC 0004 ladder** — the stent resolved into axial ring groups, each at
  its own equilibrium against a caller-supplied axial lumen profile. A
  stiff mid-lesion makes the ends open wider and `dogboning` takes a real
  value (`|d_ends − d_mids|/nominal`; the two ends are compared directly
  for a two-group ring), replacing the structural `0.0` the uniform ring
  reports. Group-level metrics reuse the Level-1 equilibrium unchanged, so
  a uniform profile reproduces `simulate_deployment` exactly.
- `StrainLifeLaw` (+ `nitinol_screening()`): **cyclic degradation of the
  strain capacity** — a screening log-log strain-life law (0.4 %
  alternating amplitude at 10⁷ cycles, a factor-of-two drop per four
  decades, per the published Nitinol fatigue band) with `amplitude_at(N)`
  and a conservative `survives(N, ε)` verdict. A screening interpolation
  of band data, not a device S–N curve: a life claim still needs vendor
  fatigue data and ASTM F2477-style pulsatile testing.
- Test groups: foreshortening (zero at the crimped diameter, monotone in
  diameter, few-percent band at realistic link fraction, upper bound at
  `link_fraction = 0`, NaN limit); non-uniform deployment (equal
  stiffnesses reproduce the uniform ring, proportional force split, forces
  summing to the ring total, no-contact carries nothing); tapered
  deployment (uniform-profile equivalence to the Level-1 ring, the
  stiff-mid-lesion dogbone with direction-agnostic magnitude and the
  two-group variant); and the strain-life law's degradation and screening
  verdicts.

### Planned
- **Level 3** — 3D superelastic FEM with frictional contact via
  `tpt-fem-hyperelastic` / `tpt-fem-contact`.
- Direct coupling to a `tpt-med-hemodynamics` solution in the same solve,
  rather than a prescribed vessel law.

### Notes
- `simulate_deployment` takes the vessel pressure–diameter law as a caller
  closure `Fn(f64 /*MPa*/) -> f64 /*mm*/`, so vessel compliance is the
  caller's model. This is a deliberate design choice, not a missing feature.
- `DeploymentResult.dogboning` remains `0.0` for the single uniform ring
  (one group has no ends-vs-middle reference); the Level-2
  `simulate_tapered_deployment` is where a real value comes from.
- `NitinolParams::default()` is a **literature-typical** parameter set. Feeding
  vendor data is the caller's responsibility, and the defaults are not a
  substitute for it.


## [0.1.0] - 2026-09-22

### Added
- **`NitinolParams`** — austenite and martensite Young's moduli, transformation
  strain `ε_L`, and the four transformation stresses `σ_ms`, `σ_mf`, `σ_as`,
  `σ_af`. `Default` gives typical ±0.1 mm laser-cut stent wire values at 22 °C
  body temperature (`E_a = 55 000` MPa, `E_m = 28 000` MPa, `ε_L = 0.05`,
  `σ_ms = 480`, `σ_mf = 560`, `σ_as = 380`, `σ_af = 260` MPa).
- **1D superelastic material** — `SuperelasticState` with `strain`, martensite
  fraction `ξ ∈ [0, 1]` and the active `Branch` (`ElasticA`, `Forward`,
  `ElasticM`, `Reverse`), and `strain_to_stress(strain, &NitinolParams)`.
  The branch is what makes the loop close: a simplified 1D Lagoudas-style model
  with a cosine transformation-hardening interface, so the same strain gives a
  different stress on loading and unloading.
- **Stent ring deployment** — `StentModel { expanded_diameter,
  crimped_diameter, n_crowns, crown_stiffness }` and
  `simulate_deployment(&stent, &nitinol, vessel_diameter_at_pressure,
  vessel_pressure)`, solving radial equilibrium for N crown springs driven
  against an artery modelled as a pressure–diameter tube law.
- **`DeploymentResult`** — `diameter`, `radial_force` (N), `contact_pressure`
  (MPa), acute `recoil` fraction (clamped to `[0, 0.2]`) and `dogboning`. The
  metric family tracked by ASTM F2394 / F2079-style bench testing.
- Typed `Pressure` from `tpt-med-units` for the deployment pressure.
- No dependencies beyond `tpt-med-units`; `#![forbid(unsafe_code)]`.

### Verification
- **Loop closure** — after loading to 8 % strain and unloading, the stress
  returns to < `1e-6` at zero strain and the state is back to `Branch::ElasticA`.
  A model that does not fully recover residual strain is not superelastic.
- **Hysteresis** — at 4 % strain the unloading stress is **strictly** below the
  loading stress. This is a strict inequality rather than a tolerance: a model
  that collapses the loop loses the physics it exists to capture, and an
  equality-only test would pass for a linear model.
- **Plateau coverage** — the peak loading stress exceeds `σ_mf`, confirming
  both transformation branches were traversed.
- **No-contact case** — an oversized stent in a large vessel returns
  `radial_force == 0.0`, asserted exactly. Silently reporting a positive force
  for a stent that never touches the artery is the worst failure mode available
  to this calculation.
- **Sign and clamp invariants** — `recoil ∈ [0, 0.2]`,
  `contact_pressure ≥ 0`, `diameter > 0`, all outputs finite.
- **Radial-force linearity** — with a fixed vessel law, doubling the crown
  count doubles the radial force, pinning the force assembly.
- Golden dataset `test-data/golden/devices/stent_expansion.json`. The ASTM
  F2394 radial-stiffness benchmark scaffold and its literature band are in
  place; the Level-3 FEM correlation remains pending.

### Known limitations
- **Level-1 ring model** — uniform ring and stent, no taper, no
  per-segment variation, no foreshortening, no 3D bending stiffness.
  `dogboning` is structurally `0.0`.
- No friction or contact mechanics; crowns are independent radial springs with
  smooth wall contact.
- 1D material: no multi-axial transformation, no Bauschinger effect.
- No fatigue or wire-fracture prediction.
- No coupling to a flow solution within the same solve.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
