# Changelog

All notable changes to `tpt-med-units` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Documented per-quantity conversion tables in the crate README, including the
  rationale for keeping `Modulus` and `Pressure` as distinct types despite both
  being megapascals.

### Notes
- The workspace MSRV is **1.82** because the quantity constructors are `const
  fn` returning floats. Lowering it below 1.82 would require dropping `const`
  float arithmetic and is a semver-minor change.
- Adding a new quantity is a **semver-minor** change. Changing any
  canonical-unit convention, or an existing `from_*`/`to_*` conversion, is a
  **semver-major** change: every downstream solver's numbers move.

## [0.1.0] - 2026-09-22

### Added
- Ten quantity newtypes generated from a single internal macro, so behaviour
  (`ZERO`, `new`, `from_*`, `to_*`, `value`, `abs`, `max`, `min`,
  `is_finite`, `UNIT`, `Display`, `Add`, `Sub`, `Mul<f64>`, `Div<f64>`,
  `Div<Self>`) is identical across all of them.
- Canonical quantities: `Length` (mm), `Pressure` (MPa), `Force` (N),
  `Density` (g/cm³), `Time` (s), `Angle` (rad), `Viscosity` (Pa·s),
  `Velocity` (mm/s), `Modulus` (MPa), `FlowRate` (mm³/s).
- Explicit non-canonical conversions: `Length::{from_cm, to_cm, from_m, to_m}`,
  `Pressure::{from_pa, to_pascal, from_kpa, to_kpa, from_mmhg, to_mmhg}`,
  `Force::{from_kn, to_kn}`, `Time::{from_ms, to_ms, from_min}`,
  `Viscosity::{from_cp, to_cp}`,
  `FlowRate::{from_ml_per_min, to_ml_per_min}`.
- `Force::body_weights(multiple, body_weight)` — the single implementation of
  the "×N body weight" loading convention, so it is applied consistently
  across the stack.
- `const fn` constructors, enabling quantities in `const` items and static
  tables.
- `Display` renders with an explicit unit suffix (`"30.000000 mm"`) so audit
  trails and CSV exports are unambiguous.

### Verification
- Round-trip conversion tests for every non-trivial conversion, locked against
  the exact constant (120 mmHg → 15 998.7 Pa ± 1 Pa, 3 cm → 30 mm, and so on).
- `#![forbid(unsafe_code)]`; no dependencies beyond `std`.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
