# Changelog

All notable changes to `tpt-med-tissue` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README stating the verification convention explicitly: closed-form
  verification compares **deviatoric** Cauchy stress `s = σ − (tr σ/3)I`, not
  full stress, because penalty formulations carry model-internal hydrostatic
  pressure at `J = 1` and the pressure-dependent components are therefore not
  unique.

### Planned
- Second-order tangent moduli per model, which the nonlinear `tpt-fem` upgrade
  path needs for Newton convergence.
- Plane-stress and reduced-order wrappers over the full 3×3 `F` interface.
- Fiber-family rotation in HGO (collagen crimp), and the two-family
  elastin/collagen parameterisation used in some literature.

### Notes
- Adding a `TissueModel` variant is a **breaking** change for any downstream
  `match` on the enum, so it is treated as semver-major despite being additive
  in appearance. New constitutive models require an RFC.
- Adding a field to a `*Params` struct is breaking; these are constructed with
  struct literals, not builders.
- Every new model must ship both a finite-difference derivative test and a
  closed-form verification case before it can be merged.

## [0.1.0] - 2026-09-22

### Added
- **Constitutive models**, each defined as a strain energy `W(F)` over the
  deformation gradient with an **analytic** first Piola–Kirchhoff stress
  `P = ∂W/∂F`:
  - `NeoHookean(NeoHookeanParams { c10, d1 })` — compressible penalty form.
  - `MooneyRivlin(MooneyRivlinParams { c10, c01, d1 })`.
  - `Yeoh(YeohParams { c1, c2, c3, d1 })`.
  - `Ogden(OgdenParams { mu, alpha, d1 })` — principal-stretch power series;
    `alpha` and `mu` are parallel vectors.
  - `HolzapfelGasserOgden(HgoParams { c, k1, k2, kappa, fiber_directions,
    d1 })` — two symmetric fiber families with dispersion `κ ∈ [0, 1/3]`
    (Gasser, Holzapfel & Ogden 2006); the structural term activates only in
    extension, as in the original.
- `TissueModel::{strain_energy, first_piola, first_piola_numerical}` — uniform
  dispatch with no generics or dynamic dispatch in the inner loop.
- `first_piola_numerical` — a central-difference reference implementation that
  serves as the **CI oracle** for every analytic derivative. If the analytic
  and numerical stresses disagree by more than `1e-6`, the build fails.
- Invariant helpers `invariant_i1`, `invariant_i2` and `principal_stretches`
  (singular values of `F`, so the Ogden form is real-valued without a signed
  convention).
- `SoftTissueMaterial { model, density, is_incompressible }` — a named
  material with physical metadata.
- `HgoParams` exposed directly for callers that want the fiber model without
  going through the enum.

### Verification
- Analytic vs. finite-difference stress for every model across a sweep of
  deformation gradients, locked to `1e-6`. This is the primary test in the
  crate: it catches a sign error, a missing invariant or a dropped `J` term.
- Uniaxial tension — Neo-Hookean, Mooney–Rivlin and Yeoh deviatoric Cauchy
  stress against the textbook incompressible solutions.
- **Ogden reduces to Neo-Hookean at `α = 2`**, an independent algebraic check.
- Volumetric branch under pure dilatation, verifying the `J` terms and the
  compressibility response.
- HGO: zero strain energy and zero stress at `F = I` (stress-free reference
  state), plus known analytic values for stretch along and transverse to the
  fiber.
- Energy objectivity asserted for all models.
- Golden dataset `test-data/golden/solid/arterial_wall_inflation.json`.

### Known limitations
- Operates on a full 3×3 `F`; no plane-stress or reduced-order wrappers.
- Derivatives are hand-written, not automatic. Keeping them in sync with `W` is
  enforced by the finite-difference test rather than by construction.
- HGO covers dispersion but not the full two-family collagen/elastin
  parameterisation.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
