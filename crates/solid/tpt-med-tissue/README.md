# tpt-med-tissue

Hyperelastic soft-tissue constitutive models — Neo-Hookean, Mooney–Rivlin,
Yeoh, Ogden, and Holzapfel–Gasser–Ogden. Four of the five (Neo-Hookean,
Mooney–Rivlin, Yeoh, HGO) have **analytic** first Piola–Kirchhoff stress;
Ogden's principal-stretch form is evaluated **numerically** (finite
difference) instead, since its analytic derivative needs eigenvectors that
are not implemented here.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--tissue-orange)](https://crates.io/crates/tpt-med-tissue)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--tissue-blue)](https://docs.rs/tpt-med-tissue)

| | |
|---|---|
| **Layer** | `solid` |
| **Status** | Alpha, `0.1.0` |
| **Scope** | RFC 0002 — hyperelastic tissue |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-geometry`](../../core/tpt-med-geometry), [`tpt-med-units`](../../core/tpt-med-units) (declared, currently unused); optional `tpt-fem-hyperelastic`/`tpt-fem-mesh` behind `substrate-cross-check` |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Soft tissue is hyperelastic: it is load-bearing at large strain, and a linear
model is wrong in the regime that matters (a vessel at 40 % strain, a ligament
in a knee at 15 %). Every constitutive model here is defined as a strain energy
function `W(F)` over the deformation gradient `F`, with stress derived
analytically.

The analytic derivative matters more than it might seem. In a nonlinear FEM
inner loop, the tangent stiffness *is* the derivative of the residual — a
finite-difference stress reintroduces noise into Newton convergence. So this
crate computes `P = ∂W/∂F` in closed form, and separately exposes
`first_piola_numerical` as a **finite-difference reference used to verify
every analytic derivative in CI**. If the two disagree by more than `1e-6`,
the build fails.

## Features

- **Neo-Hookean** — `C10` with optional compressibility via `D1`.
- **Mooney–Rivlin** — two-term invariant polynomial (`C10`, `C01`).
- **Yeoh** — third-order, numerically well-behaved at large strain.
- **Ogden** — compressible-form, arbitrary `n` terms; at `α = 2` it reduces to
  the Neo-Hookean case, which is asserted as a verification test.
  `TissueModel::first_piola` for `Ogden` is evaluated via
  `first_piola_numerical` (central difference) rather than a closed form,
  because the principal-stretch derivative needs eigenvectors of `C` that are
  not implemented analytically.
- **Holzapfel–Gasser–Ogden (HGO)** — fibrous tissue with **fiber dispersion**,
  the standard model for arterial wall and tendon/meniscus (RFC 0002), plus
  optional **collagen crimp** (`CrimpRecruitment`: progressive fiber
  recruitment through a caller-cited Gaussian recruitment distribution) and
  **per-family moduli** (`family_moduli`: the two-family
  elastin/collagen parameterisation — a compliant elastin family alongside a
  stiff collagen one).
- **Invariant helpers** — `invariant_i1`, `invariant_i2`, `principal_stretches`.
- **Uniform dispatch** — the `TissueModel` enum, so a solver can hold a
  runtime-selected model without generics or dynamic dispatch overhead in the
  inner loop.
- **A finite-difference reference implementation** (`first_piola_numerical`)
  that doubles as a test oracle.
- **`substrate-cross-check` cargo feature** (off by default, test-only) —
  cross-checks the in-house closed-form uniaxial Neo-Hookean stress against
  `tpt-fem-hyperelastic`'s independently-implemented 1-D bar Newton solve.
  Adds no production API and no default-build dependency; see
  `rfcs/0009-nonlinear-fem-substrate-adapter.md` for why this is the
  deliberately narrow first slice of substrate integration rather than a
  full 3D adapter (which the substrate does not yet provide at 0.1.0).

- **Second-order material tangents** — `material_tangent` returns
  `A[i][j](k,l) = ∂P_ij/∂F_kl`: analytic for Neo-Hookean and Yeoh (plus the
  model-independent `volumetric_tangent` of the shared penalty), central
  differences for Mooney–Rivlin/Ogden/HGO. This is the tensor a nonlinear
  Newton stiffness assembly needs, verified against finite differences and
  by the major symmetry `A_ij,kl = A_kl,ij` (the minor symmetry does not
  hold — `P` is not symmetric).

## Conventions

- `F` is the deformation gradient as a `Mat3` (right-handed, column-vector
  convention per `tpt-med-geometry`).
- Strain energy is in **MPa** (energy density = stress), moduli in **MPa**.
- `principal_stretches` returns the singular values of `F` (non-negative, so
  Ogden's `α = 2` form is real-valued without a signed convention).

## Usage

```rust
use tpt_med_geometry::{Mat3, Vec3};
use tpt_med_tissue::{NeoHookeanParams, OgdenParams, TissueModel};

/// Incompressible uniaxial tension with free lateral contraction.
fn uniaxial(lam: f64) -> Mat3 {
    Mat3::from_rows(
        Vec3::new(lam, 0.0, 0.0),
        Vec3::new(0.0, 1.0 / lam.sqrt(), 0.0),
        Vec3::new(0.0, 0.0, 1.0 / lam.sqrt()),
    )
}

fn main() {
    let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.49, d1: 2.0 });
    let f = uniaxial(1.5);

    // Analytic stress, and the finite-difference reference it is locked to.
    let p = model.first_piola(&f);
    let pn = model.first_piola_numerical(&f);
    assert!((p.at(0, 0) - pn.at(0, 0)).abs() < 1e-6);

    // Strain energy is non-negative and grows with stretch.
    let w1 = model.strain_energy(&f);
    let w2 = model.strain_energy(&uniaxial(1.2));
    assert!(w1 > w2 && w1 > 0.0);

    // Ogden at alpha = 2 reduces to the Neo-Hookean response.
    let ogden = TissueModel::Ogden(OgdenParams { mu: vec![0.49], alpha: vec![2.0], d1: 2.0 });
    let po = ogden.first_piola(&f);
    assert!((po.at(0, 0) - p.at(0, 0)).abs() < 1e-3);
}
```

Fiber-reinforced arterial wall with dispersion:

```rust
use tpt_med_geometry::{Mat3, Vec3};
use tpt_med_tissue::HgoParams;

fn main() {
    // Two symmetric fiber families, circumferentially oriented (Gasser,
    // Holzapfel & Ogden 2006).
    let hgo = HgoParams {
        c: 0.03,                       // ground-substance modulus, MPa
        k1: 0.20,                      // fiber stiffness, MPa
        k2: 12.0,                      // fiber nonlinearity, MPa^-1
        kappa: 0.10,                   // dispersion in [0, 1/3]
        fiber_directions: vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
        ],
        d1: 2.0 / 0.03,                // volumetric penalty
    };
    let f = Mat3::IDENTITY;
    // Stress-free reference state: zero energy and zero stress at F = I.
    assert!(hgo.strain_energy(&f).abs() < 1e-12);
    assert!(hgo.first_piola(&f).at(0, 0).abs() < 1e-9);
}
```

## API Overview

| Item | Purpose |
|---|---|
| `TissueModel` | Enum dispatch: `NeoHookean`, `MooneyRivlin`, `Yeoh`, `Ogden`, `HolzapfelGasserOgden` |
| `TissueModel::strain_energy(&Mat3) -> f64` | `W(F)`, MPa |
| `TissueModel::first_piola(&Mat3) -> Mat3` | `P = ∂W/∂F`, MPa — **analytic** for Neo-Hookean, Mooney–Rivlin and Yeoh; **numerical** (finite difference) for Ogden and HolzapfelGasserOgden |
| `TissueModel::first_piola_numerical(&Mat3) -> Mat3` | Central-difference reference; the CI oracle |
| `NeoHookeanParams { c10, d1 }` | Compressible Neo-Hookean; `d1 = 0` for incompressible |
| `MooneyRivlinParams { c10, c01, d1 }` | Two-term invariant polynomial |
| `YeohParams { c1, c2, c3, d1 }` | Third-order model in `Ī1 − 3` |
| `OgdenParams { mu, alpha, d1 }` | Principal-stretch series; `mu`/`alpha` are parallel vectors |
| `HgoParams { c, k1, k2, kappa, fiber_directions, d1, crimp, family_moduli }` | Fiber-reinforced model; `kappa ∈ [0, 1/3]` (1/3 = isotropic, 0 = aligned); optional crimp recruitment and per-family `(k1, k2)` overrides |
| `HgoParams::{strain_energy, first_piola}` | The same interface for HGO, exposed as `TissueModel::HolzapfelGasserOgden` |
| `CrimpRecruitment::new(mean_recruitment_stretch, spread) -> Option<_>` | Collagen-crimp recruitment `R(λ) = Φ((λ − λ̄_r)/σ_r)`; refuses non-positive spread / sub-rest mean |
| `ReducedPlaneModel::new(model, PlaneCondition::{PlaneStrain, PlaneStress})` | In-plane 2×2 `F` in, out-of-plane stretch solved (`F₃₃ = 1`, or `P₃₃ = 0` by bisection) |
| `SoftTissueMaterial { model, density, is_incompressible }` | A named material: model plus physical metadata |
| `invariant_i1(&Mat3) -> f64` | First invariant `tr C` |
| `invariant_i2(&Mat3) -> f64` | Second invariant of `C` |
| `principal_stretches(&Mat3) -> [f64; 3]` | Singular values of `F` (non-negative) |

## Verification

Code verification per ASME V&V 40 — every claim is against a closed-form or
published result, not a stored snapshot:

- **Analytic vs. finite-difference stress** for every model and a sweep of
  deformation gradients, locked to `1e-6`. This is the primary test in the
  crate: it catches a sign error, a missing invariant, or a dropped `J` term.
- **Uniaxial tension** — Neo-Hookean, Mooney–Rivlin and Yeoh deviatoric Cauchy
  stress compared against the textbook incompressible solutions.
- **Ogden ≡ Neo-Hookean at `α = 2`** — an independent algebraic check that the
  Ogden series is implemented correctly.
- **Volumetric branch** under pure dilatation, verifying the `J` terms and
  the compressibility response.
- **HGO invariants** — zero strain energy and zero stress at `F = I`; a known
  analytic value for uniaxial stretch along and transverse to the fiber.
- **Energy objectivity and stress-free reference state** are asserted for all
  models.
- **Verification compares deviatoric Cauchy stress**, `s = σ − (tr σ/3)I`, not
  the full stress. Penalty formulations carry model-internal hydrostatic
  pressure at `J = 1`, so pressure-dependent components are not unique — this
  is stated explicitly in RFC 0002 and encoded in the tests.

- **Substrate cross-check** (`substrate-cross-check` feature, off by
  default): the in-house Neo-Hookean uniaxial nominal stress (derived from
  the deviatoric Cauchy stress difference `sigma11 - sigma22`, the
  physically correct axial true stress under traction-free lateral
  surfaces) agrees with `tpt-fem-hyperelastic::solve_hyperelastic_bar`'s
  independent Newton solve on a 1-D bar mesh to `1e-9` — two separately
  implemented codebases agreeing on the same incompressible-tension
  physics. See `rfcs/0009-nonlinear-fem-substrate-adapter.md`.

Golden reference dataset: `test-data/golden/solid/arterial_wall_inflation.json`.

## Known Limitations

- No automatic differentiation — the analytic derivatives are hand-written and
  therefore must be kept in sync with `W` (which is exactly what the
  finite-difference test enforces).
- The HGO implementation covers dispersion, collagen-crimp recruitment and
  per-family moduli; crimp recruitment weights every family by the same
  `R(λ)` — per-family recruitment windows are not modelled.
- **Ogden has no analytic stress.** `TissueModel::first_piola` falls back to
  `first_piola_numerical` for `Ogden` (and, via the same enum match arm,
  currently also for `HolzapfelGasserOgden`) — even though `HgoParams`
  itself has an analytic `first_piola` when called directly rather than
  through the `TissueModel` enum.

## Related Crates

- [`tpt-med-biomechanics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-biomechanics) — the linear solver; nonlinear FEM is the documented upgrade path (RFC 0002).
- [`tpt-med-viscoelastic`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-viscoelastic) — wraps these models with a Prony series via `elastic_reference()`.
- [`tpt-med-cartilage`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-cartilage) — biphasic/poroelastic cartilage, built on a solid matrix model.
- [`tpt-med-geometry`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-geometry) — `Mat3` and `Vec3` are the function signatures.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). New constitutive
models require an [RFC](../../../rfcs) and must ship
(1) a finite-difference derivative test, and (2) a closed-form verification
case. Cite the source of every published constant in a doc comment.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.
