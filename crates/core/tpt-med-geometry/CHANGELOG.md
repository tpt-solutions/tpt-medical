# Changelog

All notable changes to `tpt-med-geometry` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README documenting the right-handed column-vector convention, the
  LPS/RAS axis definitions, and the millimetre convention the workspace
  applies to these otherwise-unitless types.

### Notes
- `EPS` (`1.0e-12`) is the shared predicate tolerance and is **public API**.
  Changing it changes the result of every geometric predicate in the stack.
- Adding a variant to `CoordinateSystem` is semver-minor; changing the
  sign convention of `flip_xy`, `lps_to_ras` or `ras_to_lps` is semver-major
  and would silently mirror every model in the workspace.

## [0.1.0] - 2026-09-22

### Added
- `Vec3` — `const fn` constructors (`new`, `splat`) and the `X`/`Y`/`Z`/`ZERO`
  constants, plus dot and cross products, `normalize`, `length`, and
  componentwise `min`/`max`/`abs`.
- `Mat3` — `IDENTITY`, `ZERO`, `from_rows`, `from_cols`, `from_array`,
  `diagonal`, indexing via `at`/`set`/`row`/`col`/`to_array`, `mul_vec`,
  `mul_mat`, `transpose`, `det`, `inverse`, `rotation_axis_angle`, `scaling`,
  `trace`, and `symmetric_eigenvalues` returning `[f64; 3]`.
- `Plane` — `from_point_normal`, `from_three_points`, `signed_distance`, and
  `ray_intersection`.
- `Aabb` — `point`, `from_points`, `expand`, `union`, `center`, `extents`,
  `contains`, `intersects`; containment and intersection are boundary
  inclusive.
- `CoordinateSystem::{Lps, Ras}` and the anatomical coordinate transforms
  `lps_to_ras`, `ras_to_lps`, `flip_xy`.
- `ImageFrame` — DICOM `ImagePositionPatient` + `ImageOrientationPatient`
  with `normal`, `voxel_position`, and `rotation`.
- `EPS: f64` = `1.0e-12`, the shared tolerance for geometric predicates.
- No dependencies beyond `std`; `#![forbid(unsafe_code)]`.

### Verification
- Unit tests for every constructor, predicate and transform, including the
  LPS↔RAS round trip and image-frame normalisation.
- `Mat3::inverse` returns `Option` and yields `None` on singular input rather
  than a silent `NaN`.
- `symmetric_eigenvalues` is exercised on the stress tensors produced by
  `tpt-med-biomechanics`, so principal-stress regressions surface here.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
