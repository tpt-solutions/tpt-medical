# Changelog

All notable changes to `tpt-med-viscoelastic` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README documenting the Abaqus/FEBio **relative** Prony convention
  (`gᵢ = Gᵢ/G₀`, `Σgᵢ ≤ 1`) and why a relative series is preferred: a
  glass-modulus change does not invalidate a fitted series.

### Planned
- Temperature shifting of relaxation times (WLF and Arrhenius relations), so
  `τᵢ` values can be generated rather than supplied.
- Non-linear hyperviscoelastic formulations, applying the Prony series to the
  hyperelastic energy in finite strain rather than in the linear-viscoelastic
  regime.
- A time-integration helper, so callers driving a finite-element inner loop do
  not each reimplement the recurrence.

### Notes
- The series is defined **in shear**; the volumetric response comes from the
  `tpt-med-tissue` reference model and is not itself rate-dependent. Making
  it rate-dependent is an RFC-scale change.
- `ViscoelasticMaterial::validate` returns `Result<(), String>` deliberately:
  the error text names the offending term, which a bare boolean cannot.

## [0.1.0] - 2026-09-22

### Added
- `PronyTerm { g_i, tau_i }` — one Maxwell element, with `g_i` a
  **relative, dimensionless** shear modulus and `tau_i` in seconds.
- `ViscoelasticMaterial { g0, prony }` — the glass (instantaneous) shear
  modulus `G0` in MPa plus the series.
- `::validate()` — enforces `G0 > 0`, `gᵢ ≥ 0`, `Σgᵢ ≤ 1` and `τᵢ > 0`,
  returning a descriptive error naming the offending term.
- `::equilibrium_modulus()` — `G∞ = G0(1 − Σgᵢ)`.
- `::relaxation_modulus(time)` — `G(t) = G∞ + Σ Gᵢ exp(−t/τᵢ)`.
- `::storage_modulus(omega)` — `G'(ω) = G∞ + Σ Gᵢ(ωτᵢ)²/(1+(ωτᵢ)²)`.
- `::loss_modulus(omega)` — `G''(ω) = Σ Gᵢ(ωτᵢ)/(1+(ωτᵢ)²)`.
- `::loss_tangent(omega)` — `G''/G'`, the scalar a damping requirement is
  usually written in.
- `::step_strain_stress(gamma0, time)` — the classic step-strain protocol, for
  comparison against a published relaxation curve.
- `::elastic_reference() -> TissueModel` — the `tpt-med-tissue` model supplying
  the glass response, so the hyperelastic implementation is not duplicated
  across crates.
- A documented re-export shim, `tpt_med_tissue_link`, mapping this crate's
  types onto `tpt-med-tissue` names for input-file compatibility.
- Follows the Abaqus/FEBio convention of **relative** Prony moduli, so a
  material defined here ports to a commercial solver without
  reparameterisation.

### Verification
Viscoelasticity has unusually strong analytic structure, and every limit is
asserted:
- **Short time** — `G(0) = G0`, `G'(∞) → G0`, `G''(∞) → 0`.
- **Long time** — `G(t) → G∞` as `t → ∞`, and `G'(0) = G∞`, `G''(0) = 0`.
- **Monotonicity** — `G(t)` strictly decreasing and bounded below by `G∞`;
  `G'(ω)` non-decreasing and bounded above by `G0`.
- **Step-strain stress** decays from `G0·γ₀` toward `G∞·γ₀`, monotone in time.
- **Loss tangent** is asserted equal to `G''/G'` and strictly positive for any
  non-degenerate series.
- **Linear scaling** — multiplying `g0` by `k` multiplies every modulus output
  by `k`, while the loss tangent is **invariant**, which is a clean test that
  the series is genuinely relative.
- **Validation** — `Σgᵢ > 1`, negative `gᵢ`, non-positive `τᵢ` and `G0 ≤ 0` are
  each rejected with their specific message.
- **Single-term reduction** — a one-term series is compared against the exact
  closed-form Maxwell relaxation `G₀ − (G₀−G∞)exp(−t/τ)`.

### Known limitations
- Shear only; the volumetric response is not rate-dependent.
- Linear viscoelasticity — large-strain and hyperviscoelastic formulations are
  not implemented.
- Temperature dependence is the caller's responsibility; `τᵢ` values are
  supplied, not shifted.
- No time-integration helper for a finite-element inner loop.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
