# Changelog

All notable changes to `tpt-med-electrophysiology` are documented here. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently**
on a six-week cadence. Workspace-level and cross-cutting changes are recorded
in the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only
what changes for consumers of this crate.

## [Unreleased]

### Planned
- Stage 2 (ECG/EGM forward problem via pseudo-bidomain lead-field
  projection) and Stage 3 (ablation screening), both kept at roadmap depth
  in `rfcs/0005-cardiac-electrophysiology.md` pending Stage 1 usage and,
  for Stage 3, clinical-data validation.
- Anisotropic (fiber-direction) conductivity — needs a fiber-field source
  (atlas or DTI derivation) this crate has no source for yet.
- Promotion path to `tpt-science`'s electrophysiology crate for ionic-model
  breadth beyond Mitchell-Schaeffer, once a real workflow needs it.

## [0.1.0] - 2026-09-27

Initial release: Stage 1 of `rfcs/0005-cardiac-electrophysiology.md`.

### Added
- **`MitchellSchaefferParams`** — validated Mitchell & Schaeffer (2003)
  two-current ionic kinetics (`dV/dt = J_in + J_out + J_stim`,
  `dh/dt` switching at `v_gate`). `new()` rejects non-finite or non-positive
  time constants and a `v_gate` outside `(0, 1)`.
  `human_ventricular_default()` is the commonly-reproduced generic-
  ventricular fit from the original paper's Table 1 — a screening default,
  not a validated patient- or species-specific one.
- **`MonodomainTissue`** — isotropic, homogeneous-diffusivity monodomain
  reaction-diffusion, built via `from_mask` directly from a
  `tpt_med_meshing::SegmentationMask`. No-flux (Neumann) boundary at the
  mask's own edge and at any tissue/non-tissue voxel interface, via a
  discrete Laplacian that only sums over neighbours that both exist on the
  grid and are solid tissue.
- **Explicit RK2 (Heun) time stepping** (`step`), with `max_stable_dt()`
  computing the smaller of a diffusion-CFL bound (accounting for
  anisotropic voxel spacing) and a reaction-stiffness bound
  (`0.2 * tau_in`). The reaction bound was added during development after
  a diffusion-only bound let an explicit step overshoot the upstroke's
  excitable manifold into `NaN` rather than saturating — see the RFC's
  Verification strategy and this crate's README.
- **`stimulate(region, amplitude)` / `activation_map()`** — sets the active
  stimulus for exactly the next `step()` call; `activation_map()` returns
  the first `V >= v_gate` crossing time per voxel.
- **`S1S2Protocol`** — a single-cell (0D, no diffusion) S1-S2 restitution
  protocol, `run()` returning `(diastolic_interval, apd)` pairs. Used as
  this crate's own code-verification fixture: no closed-form restitution
  curve exists for this model, so correctness is checked against the
  qualitative shape (APD increasing monotonically with DI) Mitchell &
  Schaeffer's own paper reports, not an analytic target.
- **`EpError`**: `EmptyMask`, `NonFiniteParameter`, `UnstableTimeStep`,
  `InvalidDiffusivity`, `RestitutionFailed`.
- New golden dataset
  `test-data/golden/electrophysiology/monodomain_restitution.json` — pins
  the exact `S1S2Protocol` output for a specific parameter set and
  protocol, a code-verification fixture (no external reference exists).

### Known limitations
- Isotropic conductivity only; no anisotropic (fiber-direction) support.
- `S1S2Protocol` is single-cell; not a tissue-level restitution behaviour.
- No Stage 2 (ECG) or Stage 3 (ablation screening).

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
