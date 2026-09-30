# Changelog

All notable changes to `tpt-med-cartilage` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- `PermeabilityLaw` / `ConstantPermeability` /
  `StrainDependentPermeability`: the **strain-dependent permeability
  hook** — the mechanism half of the nonlinear-biphasic item. Callers
  supply a cited `k(J)` law (closure or type); the crate validates the
  evaluation and provides the equilibrium-compaction evaluation point.
  Constant-`k` remains the closed-form baseline: the creep series solve
  the constant-coefficient problem, so a strain-dependent `k` is
  evaluated per compaction zone for caller-driven stepping, not
  substituted into the series.
- `BiphasicMaterial::unconfined_equilibrium_modulus`: the **unconfined
  equilibrium limit** — `E_s = H_A(1+ν_s)(1−2ν_s)/(1−ν_s)`, equal to
  `H_A` at the cartilage default `ν_s = 0` (confined and unconfined
  long-time responses coincide) and softer for `ν_s > 0`. Both analytic
  limits (rigid instantaneous, this equilibrium) are exact and tested;
  the transient between them remains the classical Bessel-series
  solution, deliberately not reproduced from memory.
- `BiphasicMaterial::solid_shear_modulus`: the **shear half of the
  unconfined-shear boundary-condition item** — first-order biphasic shear
  produces no volumetric strain, so the interstitial fluid never
  pressurises and the response is the solid matrix at all times
  (`G = H_A(1−2ν_s)/(2(1−ν_s))`, permeability-independent; `None` at
  `ν_s ≥ ½`). A closed form rather than a solve, which is *why* there is
  no shear transient to implement.
- Crate README explaining why poroelastic rather than elastic: interstitial
  fluid transport, not solid elasticity, dominates the time-dependent
  response, and a purely elastic model gets the equilibrium roughly right and
  the timing completely wrong.

### Planned
- Unconfined compression (its analytic transient needs the classical Bessel-series coefficients; shear, the other half, is delivered above as a closed form).
- Nonlinear biphasic theory proper (the strain-dependent permeability
  hook is delivered; a full nonlinear solid matrix and a time-stepping
  solver for non-constant `k` remain).
- A lubrication/repulsion term for the contact interface, so the model can be
  driven by a contact solver rather than a prescribed step load.
- Fibrous-cartilage support (a fibre-reinforced solid matrix).

### Notes
- The **term count** is a runtime parameter, not a hard-coded truncation. The
  series is truncated conservatively, and increasing `terms` moves the result
  monotonically toward the analytic limit; this is asserted.
- Changing `BiphasicMaterial::default()` is a semver-minor change that alters
  results for any caller relying on it, and requires a V&V re-run.

## [0.1.0] - 2026-09-22

### Added
- `BiphasicMaterial { aggregate_modulus, permeability, poissons_ratio,
  thickness }` — linear biphasic theory (Mow, Kuei & Lai 1980) in the
  confined-compression configuration: a solid matrix with aggregate modulus
  `H_A` and Poisson ratio ≈ 0, saturated with fluid moving by Darcy's law.
- `Default` screening values for adult articular cartilage: `H_A = 0.7 MPa`,
  `k = 0.002 mm⁴/(N·s)`, `ν = 0`, `h = 2.0 mm`.
- `creep_displacement_fraction(sigma0, time, terms)` — the classical
  series solution, with a caller-controlled term count so convergence is an
  explicit input.
- `equilibrium_strain(sigma0) = σ₀ / H_A` — the asymptotic solid-supported
  limit, computed independently of the series so the two can be cross-checked.
- `initial_displacement_fraction(sigma0)` — the `t = 0` value.
- `fluid_pressure_fraction(time, terms)` — the fraction of the load still
  carried by interstitial fluid, the fluid-supported phase.
- `gel_time()` — the time to 50 % of equilibrium displacement, the inverse of
  the classical gel time.
- Pure functions, no allocation and no state, so the model embeds trivially in
  a contact solver or a browser loop.
- No dependencies; `#![forbid(unsafe_code)]`.

### Verification
Confined compression has closed-form limits, and every one is asserted:
- **Short time** — `u(t) → 0` as `t → 0`, with the fluid carrying the full
  load, so `fluid_pressure_fraction(0⁺) → 1`.
- **Long time** — `u(t) → h·σ₀/H_A` as `t → ∞`, compared against the
  independent `equilibrium_strain` computation rather than against itself.
- **Gel time** — at `gel_time()` the displacement fraction is 0.5, verified
  against the closed form rather than by bisection.
- **Monotonicity** — creep is non-decreasing in time and bounded above by
  equilibrium, across a sweep of `σ₀`, `k` and `h`.
- **Permeability sensitivity** — a lower `k` delays the response but does not
  change the equilibrium, which a wrong `H_A` implementation would break
  immediately.
- **Term-count convergence** — increasing `terms` moves monotonically toward
  the analytic limit, so the truncation is always conservative.

### Known limitations
- Linear theory only. Confined compression at 100 kPa is within the linear
  regime; impact loading is not.
- Confined compression only; unconfined compression and shear are not
  implemented.
- No nonlinear permeability–strain coupling and no lubrication term.
- 1D through-thickness; a full 3D poroelastic solve is out of scope.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
