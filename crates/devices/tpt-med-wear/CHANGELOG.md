# Changelog

All notable changes to `tpt-med-wear` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README stating that wear coefficients are **inputs, never defaults**,
  with the reasoning: coefficients in this domain span orders of magnitude
  between materials, bearing designs and test protocols, so a library-level
  default would be a fabricated number wearing a citation.

### Planned
- Wear-debris-induced damage feedback, so wear changes the contact geometry and
  pressures — without which the runaway that ends real implant life is not
  captured.
- A coupling helper to a contact solver, so pressures and sliding distances can
  be solved rather than supplied.
- Uncertainty propagation over the wear coefficient, which scatters over orders
  of magnitude between studies and which a defensible screening study should
  quantify.
- A run-in period and activity-level variation, so gait extrapolation is not
  strictly linear in cycle count.

### Notes
- Changing `WearLaw` or `WearResult` is a breaking change. Adding a new law
  variant is breaking for downstream exhaustive matches, even though it looks
  additive, and requires an RFC.
- `exceeds_iso14879_screen` takes the limit from the caller rather than
  hard-coding one, because acceptable wear is device- and material-specific.

## [0.1.0] - 2026-09-22

### Added
- **`WearLaw::Archard { k }`** — the classical law `V = k · F · s`, with `k` the
  specific wear rate in mm³/(N·m), absorbing everything the model does not
  represent.
- **`WearLaw::CrossLand { k, pressure_threshold }`** — the Archard-type
  pressure-threshold law `dh = K·(p − p₀)·ds`, with `p₀` in MPa. The threshold
  is the physically important feature: wear is **exactly zero** below it, which
  is what makes a well-functioning, low-contact-pressure bearing nearly
  wear-free while a mis-loaded one fails catastrophically rather than
  proportionally.
- `WearModel { law, gait_cycles }` and
  `::simulate_wear(&contact_pressures, &sliding_distances, area)`, taking
  **parallel per-zone arrays** so a real bearing with a non-uniform contact
  pattern is representable. `area` (mm²) converts volumetric to linear wear.
- `WearResult` — `volumetric_wear` (mm³), `linear_wear` (mm) and
  `wear_per_megacycle` (mm³/Mc), the last being the unit ISO 14879 expresses
  limits in.
- `WearModel::exceeds_iso14879_screen(&WearResult, limit_mm3_per_mc)`, so a
  screening verdict is one call rather than a hand calculation. The limit is
  caller-supplied; ~30 mm³/Mc is commonly quoted for UHMWPE tibial inserts.
- No dependencies, no allocation; trivially embeddable in a parameter sweep.
- `#![forbid(unsafe_code)]`.

### Verification
Both laws have exact algebraic forms, so the checks are hand-computed values
rather than stored snapshots:
- **Archard closed form** — compared against
  `k · Σ(μp_i · A · s_i · N)` computed independently, pinning the
  pressure-to-force conversion (`F = p · A`) and the metre/millimetre
  conversion in the same test.
- **Cross–Land threshold behaviour** — a zone *exactly at* `p₀` produces
  **zero** wear, and a zone just above it produces a positive amount that
  vanishes as `p → p₀⁺`. This is the defining property of the law.
- **Cross–Land above threshold** reduces to Archard with `k' = k(1 − p₀/p)`,
  asserted against the algebraic identity.
- **Unit conversions** — `gait_cycles == 1_000_000` is asserted to give
  `wear_per_megacycle == volumetric_wear`, with the mm³/mm and mm³/Mc units
  each pinned by an independent calculation.
- **Linear wear** is asserted to equal `volumetric_wear / area` exactly.
- **Multi-zone consistency** — splitting one zone into two with proportional
  areas reproduces the same total.
- **Monotonicity** — wear is non-decreasing in `gait_cycles`, contact pressure
  and sliding distance; a strict physical requirement and a cheap detector for
  sign errors.
- Golden dataset `test-data/golden/devices/knee_wear_10mcycles.json`.

### Known limitations
- No wear-debris-induced damage feedback, so the contact geometry does not
  evolve. Screening tool, not a life model.
- No lubricant or contact-mechanics solution: pressures and sliding distances
  are inputs, not solved quantities.
- Steady-state, single-condition extrapolation, linear in `gait_cycles`.
- No coefficient uncertainty propagation.
- Not a substitute for ASTM F2028 or ISO 14879 standard testing.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
