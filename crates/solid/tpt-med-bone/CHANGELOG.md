# Changelog

All notable changes to `tpt-med-bone` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- `ModulusLaw` trait (+ `PowerLaw`, and a blanket impl for closures) and
  `BoneMaterial::from_hu_with_law`: **QCT calibration hook** — a study can
  evaluate a phantom-fitted density→modulus relation instead of the default
  power law.
- `BoneRemodelingModel::remodel_field`: **spatial remodeling** — drives a
  per-voxel density field from a solved stimulus field (e.g. SED per element
  from a `tpt-med-biomechanics` result), with per-voxel viable clamping.
- Crate README distinguishing the HU→modulus correlations (owned by
  `tpt-med-dicom`) from the structural description, anisotropy and remodeling
  law owned here, and stating explicitly that the reference constants are
  literature-typical population values rather than patient-calibrated ones.

- `ResorptionDeadline` with
  `BoneRemodelingModel::{update_density_with_deadline, remodel_field_with_deadline}`:
  the **disuse/resorption-deadline model** — per-voxel days of continuous
  disuse are tracked in a caller-held counter, reset by reloading, and past
  the deadline resorption runs at a configurable multiplier. Apposition
  ignores the deadline.
- `RateAugmentation`: **load-rate-dependent remodeling** (screening
  heuristic after Turner's loading-rule observations) — the stimulus is
  multiplied by a log-scaled, saturation-capped factor above a reference
  quasi-static rate. `update_density_with_deadline` composes with it at the
  call site.
- Four new tests: deadline gating (baseline before, ×multiplier after),
  counter reset on reload, 180-day shielded-voxel field run losing strictly
  more than the plain law, and rate augmentation (identity at/below
  reference, log growth, cap, lazy-zone escape).

### Notes
- Reference material constants are **screening values**. Changing any of them
  changes every downstream stress field, so it is a semver-minor change that
  requires an RFC, a citation, and a V&V re-run.
- Adding a variant to `TissueClass`, `Anisotropy` or `RemodelingStimulus` is
  breaking for downstream exhaustive matches.

## [0.1.0] - 2026-09-22

### Added
- `BoneMaterial` — `bone_type`, `tissue_class`, `youngs_modulus`,
  `poissons_ratio`, `yield_stress`, `ultimate_stress` and `anisotropy`.
- `BoneMaterial::cortical_reference(BoneType)` and
  `::trabecular_reference(BoneType)` — literature-typical materials per bone.
- `BoneMaterial::from_hu(hu, bone_type, poissons_ratio)` — HU to material,
  classifying cortical versus trabecular at 1.3 g/cm³ apparent density
  (≈300 HU) and applying the matching `tpt-med-dicom` power law automatically.
- `Anisotropy` — `Isotropic`,
  `TransverselyIsotropic { axis, e_long, e_trans }` and
  `Orthotropic { axes, e }`; `BoneMaterial::effective_modulus(direction)`
  performs a directional projection rather than returning a scalar, and falls
  back to `youngs_modulus` when isotropic.
- `TissueClass` — `Cortical`, `Trabecular`; `reference_modulus(TissueClass)`.
- **Wolff's-law remodeling** — `BoneRemodelingModel` with an explicit **lazy
  zone** (the mechanostat's dead band) around the reference stimulus,
  separate `apposition_rate` and `resorption_rate`, a response that is linear
  in the normalised over/under-stimulus
  (`Δρ = rate · (S/S_ref − 1)` with `S_ref` the nearer lazy-zone bound, a 3×
  overshoot cap), and a viability clamp. `Default` uses Frost-style
  screening parameters: reference SED ≈0.004 mJ/mm³, ±35 % lazy zone
  `(0.0026, 0.0054)`, apposition 0.003 and resorption 0.002 g/cm³/day.
- `update_density(current, stimulus, dt_days)` and
  `simulate_days(start, stimulus, days)`.
- `RemodelingStimulus` — `StrainEnergyDensity`, `PrincipalStrain`,
  `DamageAccumulation`, so the remodeler can be driven by whatever the solver
  actually produces.
- `clamp_viable(density, (low, high))` for density windows.
- `Density` and `Modulus` from `tpt-med-units` throughout, so the working unit
  set cannot drift.

### Verification
- HU → modulus is locked to the `tpt-med-dicom` power laws, and the
  cortical/trabecular split is asserted on both sides of 1.3 g/cm³ and exactly
  at the boundary.
- Reference materials are asserted to be ordered cortical > trabecular for
  every `BoneType`, which catches a transposed constant immediately.
- `effective_modulus` is pinned at the symmetry axis, perpendicular to it, and
  for the isotropic fallback, plus a rotational invariance check.
- Remodeling is verified in all three regimes — inside the lazy zone (no
  change), above (apposition), below (resorption) — with the
  apposition/resorption asymmetry and the viability clamp asserted, and
  `update_density` asserted monotone in the stimulus.
- Golden datasets `test-data/golden/solid/` and
  `test-data/golden/devices/hip_stem_micromotion.json`.

### Known limitations
- Reference constants are population-typical, not patient-calibrated.
- Remodeling is lumped and spatially uniform per call; per-element remodeling
  from a solved strain-energy field is a caller-side loop.
- No disuse or resorption-deadline model; the lazy zone is the only memory.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
