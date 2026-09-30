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
- **Collagen crimp in HGO (`CrimpRecruitment` + `HgoParams::crimp`)** —
  the recruitment half of the fiber-rotation/crimp roadmap item. Real
  collagen fibers are wavy at rest and straighten progressively: with
  crimp set, each family's fiber term is weighted by the recruited
  fraction `R(λ) = Φ((λ − λ̄_r)/σ_r)` at fiber stretch `λ` (the Gaussian
  waviness-distribution treatment of Decraemer, Maes & Vanhuyse 1980).
  The caller cites `λ̄_r` and `σ_r` — the mechanism ships, the coefficients
  come from the caller's source, as everywhere in this workspace;
  `CrimpRecruitment::new` refuses a non-positive spread and a sub-rest
  mean. `crimp: None` — the default — is the unmodified HGO model,
  bit-identical; `R → 1` recovers it exactly, below recruitment the fiber
  term is silent. The analytic `HgoParams::first_piola` carries the
  chain-rule `R′` term, and the normal CDF is evaluated to near machine
  precision (convergent positive series plus the continued fraction of
  A&S 7.1.14): a 1e-7 rational approximation was tried first and
  **rejected** — the finite-difference stress reference differentiates
  through `R`, so a CDF accurate only in value disagrees with its own
  derivative inside the recruitment window (~2e-5 relative). Five new
  tests: silence below recruitment, exact recovery above it,
  analytic-vs-FD stress through the window, monotone stiffening, the
  spread ordering (fraction level and rise steepness), and constructor
  validation including the erf accuracy.
- **`ReducedPlaneModel` (+ `PlaneCondition`, `PlaneSolution`)**:
  plane-strain and plane-stress wrappers over the full 3×3 `F` interface —
  the caller supplies the four in-plane gradient components; plane strain
  pins `F₃₃ = 1`, plane stress solves `P₃₃ = 0` by bracketed bisection on
  the scalar `F₃₃` (`P₃₃` is monotone decreasing in `F₃₃` for every law in
  this crate, which is what makes the scalar bracket robust). Verified
  against the incompressible closed form `F₃₃ = 1/det F₂ₓ₂` and
  traction-freeness.
- Crate README stating the verification convention explicitly: closed-form
  verification compares **deviatoric** Cauchy stress `s = σ − (tr σ/3)I`, not
  full stress, because penalty formulations carry model-internal hydrostatic
  pressure at `J = 1` and the pressure-dependent components are therefore not
  unique.
- **`substrate-cross-check` cargo feature** (off by default), per
  `rfcs/0009-nonlinear-fem-substrate-adapter.md`: cross-checks the in-house
  closed-form uniaxial Neo-Hookean stress against
  `tpt-fem-hyperelastic::solve_hyperelastic_bar`'s independent 1-D bar
  Newton solve. Adds two optional dependencies
  (`tpt-fem-hyperelastic`, `tpt-fem-mesh`, both pinned `=0.1.0` in the root
  workspace manifest) and one test module; no production API, no change to
  the default build.
- **`TissueModel::volumetric_first_piola`** — the volumetric part of the first
  Piola, `d/dF [(J-1)^2/d1] = 2J(J-1)/d1 * F^-T`, split out from
  `first_piola`. Additive: `first_piola` is unchanged and still returns the
  fused deviatoric-plus-volumetric stress.
  - Added for `tpt-med-fem-adapter`'s selective reduced integration, which
    needs the deviatoric and volumetric responses integrated on different rules.
    The split is not recoverable from a single `Mat3` in general, so the law has
    to expose it.
  - One closed form serves all five variants, because they share the identical
    `(J-1)^2/d1` penalty. Returns zero for `J <= 0`, matching the guard in
    `first_piola`.
  - This is also the prerequisite for a mixed `u`-`p` formulation, where the
    pressure unknown is exactly this term.

- **`TissueModel::material_tangent` / `material_tangent_numerical` /
  `volumetric_tangent` (+ `MaterialTangent`)**: **second-order material
  tangents** `A[i][j](k,l) = ∂P_ij/∂F_kl` — the tensor a nonlinear Newton
  solve assembles element stiffness from. Analytic for the Neo-Hookean and
  Yeoh families (which share the `P_dev = 2β·q·G` structure; the `q′` term
  carries a second `β` because `∂s/∂F` has its own) plus the
  model-independent volumetric tangent of the shared `(J−1)²/d1` penalty;
  central differences through `first_piola` for Mooney–Rivlin, Ogden and
  HGO (same policy as `first_piola`'s analytic/numerical split; FD-of-FD
  round-off ~1e-4 relative is the consumer's tolerance for those three).
  Zero for `J <= 0`, matching the stress-side guards. Four new tests:
  analytic-vs-FD agreement, major symmetry `A_ij,kl = A_kl,ij` for all five
  models (the *minor* symmetry deliberately not asserted — `P` is not
  symmetric), volumetric tangent vs FD of `volumetric_first_piola`, and
  positivity/inverted-configuration guards.

### Added
- **`HgoParams::family_moduli` — the two-family elastin/collagen
  parameterisation**: per-family `(k1, k2)` overrides parallel to
  `fiber_directions` (a compliant elastin family alongside a stiff
  collagen family), replacing the shared `k1`/`k2` for the families they
  cover. `None` — the default — is the shared pair, bit-identical; a
  length mismatch is an assert-level API error, like the stent crate's
  paired-slice contracts. Because the fiber term sums over families, the
  two-family response is verified as the *exact* sum of the single-family
  responses (an identity, not a tolerance), alongside finite-difference
  stress agreement for the full combination (per-family moduli + crimp
  together), shared-pair equivalence, and the pairing panic. Crimp
  recruitment weights every family by the same `R(λ)` — per-family
  recruitment windows remain open.

### Planned
- (The fiber item is delivered: crimp recruitment above, per-family
  moduli in this release. Remaining as a small nuance: per-family
  recruitment windows, so an elastin family can be engaged before a
  crimp-gated collagen one.)

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
