# RFC 0002: Hyperelastic Tissue Models

- **Status:** Accepted
- **Started:** 2026-09-20
- **Crates:** `tpt-med-tissue`, `tpt-med-biomechanics`, `tpt-med-viscoelastic`

## Summary

Select the constitutive model set, stress formulations, and verification
policy for soft tissue, and define the boundary between the in-house linear
FEM core and the `tpt-fem` substrate for nonlinear analysis.

## Motivation

Soft tissue is nonlinear, nearly incompressible, and (for arteries)
fiber-reinforced. The model set must cover screening-grade clinical questions
with verifiable mathematics, while acknowledging that large-deformation
contact (stent deployment, joint replacement) needs substrate-grade solvers.

## Design

### Model set (implemented)

| Model | Energy | Primary use |
|---|---|---|
| Neo-Hookean | `C10(Ī1−3) + (J−1)²/D1` | fat, muscle screening |
| Mooney–Rivlin | `+ C01(Ī2−3)` | elastomers, skin |
| Yeoh | `Σcᵢ(Ī1−3)ⁱ` | rubber-like seals, tissue over wide strain |
| Ogden | principal-stretch power series | artery bulking data fits |
| HGO | ground substance + 2 dispersed fiber families | arterial wall |

All share the deviatoric–volumetric split with penalty incompressibility;
HGO uses the non-deviatoric I1/I4 form with extension-only fiber activation.

### Verification policy (ASME V&V 40 code verification)

1. Every analytic `∂W/∂F` is locked to a central-difference reference in CI
   (catches invariant/algebra errors independent of closed forms).
2. Closed-form verification uses the **deviatoric Cauchy stress** `s = σ −
   (tr σ/3)I`: penalty formulations carry model-internal hydrostatic
   pressure at J = 1, so pressure-dependent components are not unique.
   Verified: NH/MR/Yeoh uniaxial deviatorics vs textbook incompressible
   solutions; Ogden α=2 ↔ NH equivalence; volumetric branch under pure
   dilatation.
3. The FEM core's Q1 hexes are verified by exact uniaxial, patch, and
   rigid-translation tests, plus a cantilever vs Euler–Bernoulli comparison
   with a documented locking band (full 2×2×2 integration on Q1 hexes
   shear-locks; the band is a property of the discretisation, recorded in
   the golden dataset).

### Boundary with `tpt-fem`

The in-house core deliberately covers: structured voxel hexes, linear
elasticity, CG solving, WASM-safe `std`-only code. `tpt-fem` (pinned in the
workspace manifest) is the sanctioned path for: unstructured meshes,
large-deformation hyperelasticity with `tpt-fem-hyperelastic`, frictional
contact with `tpt-fem-contact`, and nonlinear solvers. Integration lands as
an adapter crate (`tpt-med-biomechanics` trait → `tpt-fem` model) behind a
cargo feature so default builds stay dependency-free; each substrate bump
re-runs the full V&V suite since credibility evidence is version-sensitive.

## Alternatives considered

- **Fung orthotropic** for myocardium: deferred until electrophysiology
  (RFC 0005) needs it; energy form is in the literature queue.
- **Compressible Ogden (Helfer et al.)**: the deviatoric form covers
  current needs; revisit if foam-like tissue appears.

## Unresolved questions

- Growth/remodeling coupling (arterial wall adaptation) — needs a
  thermodynamically consistent framework; candidate for a later RFC.
- Damage (contusion/tear) criteria for trauma screening.
