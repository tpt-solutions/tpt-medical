# tpt-med-viscoelastic

Prony-series viscoelasticity for soft tissue — generalized Maxwell relaxation,
frequency-domain storage and loss moduli, and loss tangent.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--viscoelastic-orange)](https://crates.io/crates/tpt-med-viscoelastic)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--viscoelastic-blue)](https://docs.rs/tpt-med-viscoelastic)

| | |
|---|---|
| **Layer** | `solid` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-tissue`](../tpt-med-tissue) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Biological soft tissue is rate-dependent: the same stress applied for 10 ms
and for 10 s produces different strains. A purely elastic model cannot
represent stress relaxation (arterial pressure decay), creep, or hysteresis
(energy returned per cycle < energy applied), and energy loss per cycle is
exactly what governs fatigue failure in heart valve leaflets and tendon.

A generalized Maxwell (Prony) series is the standard way to express this, and
critically, it is the **format Abaqus and FEBio already accept**. A material
defined as a glass modulus plus a relative Prony series ports to a commercial
solver without reparameterisation, which is the difference between a research
artifact and something a design engineer can use.

## Features

- **Generalized Maxwell in shear:** `G(t) = G∞ + Σᵢ Gᵢ exp(−t/τᵢ)`
- **Frequency domain:** `G'(ω) = G∞ + Σᵢ Gᵢ(ωτᵢ)²/(1+(ωτᵢ)²)` and
  `G''(ω) = Σᵢ Gᵢ(ωτᵢ)/(1+(ωτᵢ)²)`
- **Loss tangent** `G''/G'` — the scalar a fatigue or damping requirement is
  usually written in.
- **Step-strain stress** — the classic experimental protocol, for comparing a
  model against a published relaxation curve.
- **Input validation** — `validate()` enforces positive `G0`, non-negative
  terms, `Σgᵢ ≤ 1`, and positive relaxation times, and returns a *descriptive
  string* rather than a bare `false`.
- **Elastic reference model** — `elastic_reference()` hands back the
  [`tpt-med-tissue`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-tissue)
  `TissueModel` that supplies the glass response, so the two crates compose
  instead of duplicating the hyperelastic implementation.
- A documented re-export shim (`tpt_med_tissue_link`) that maps this crate's
  types onto `tpt-med-tissue` names for input-file compatibility.

## Conventions

- `g0` (instantaneous/glass shear modulus `G0`) in **MPa**.
- `PronyTerm.g_i = G_i/G0` is a **relative, dimensionless** modulus with
  `Σgᵢ ≤ 1` — the Abaqus/FEBio convention, deliberately not an absolute
  modulus, so a glass-modulus change does not invalidate a fitted series.
- `tau_i` (relaxation time `τᵢ`) in **seconds**; frequencies in **rad/s**.

## Usage

```rust
use tpt_med_tissue::{NeoHookeanParams, TissueModel};
use tpt_med_viscoelastic::{PronyTerm, ViscoelasticMaterial};

fn main() {
    let m = ViscoelasticMaterial {
        g0: 0.5,                                  // MPa, glass shear modulus
        prony: vec![
            PronyTerm { g_i: 0.3, tau_i: 0.01 },  // fast collagen
            PronyTerm { g_i: 0.2, tau_i: 1.0 },   // slow collagen
        ],
    };
    m.validate().unwrap();

    // Long-term modulus: G∞ = G0 (1 - Σgᵢ) = 0.5 * 0.5 = 0.25 MPa.
    assert!((m.equilibrium_modulus() - 0.25).abs() < 1e-12);

    // Relaxation decays monotonically from G0 toward G∞.
    let g0 = m.relaxation_modulus(0.0);
    let ginf = m.relaxation_modulus(1.0e9);
    assert!((g0 - 0.5).abs() < 1e-12);
    assert!(ginf < g0 && ginf > m.equilibrium_modulus() - 1e-6);

    // A step strain of 50% relaxes toward the equilibrium stress.
    let s0 = m.step_strain_stress(0.5, 0.0);
    let s_inf = m.step_strain_stress(0.5, 1.0e9);
    assert!(s_inf < s0);

    // Frequency-domain response and loss tangent.
    let g_prime = m.storage_modulus(10.0);
    let g_double_prime = m.loss_modulus(10.0);
    assert!(g_prime > 0.0 && g_double_prime > 0.0);
    assert!((m.loss_tangent(10.0) - g_double_prime / g_prime).abs() < 1e-12);

    // The glass response comes from tpt-med-tissue, not duplicated here.
    let reference: TissueModel =
        m.elastic_reference().NeoHookean(NeoHookeanParams { c10: 0.5, d1: 0.0 });
    assert_eq!(reference, TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 0.0 }));

    // Invalid series are rejected with a descriptive message.
    let bad = ViscoelasticMaterial { g0: 0.5, prony: vec![PronyTerm { g_i: 1.5, tau_i: 1.0 }] };
    assert!(bad.validate().is_err());
}
```

## API Overview

| Item | Purpose |
|---|---|
| `PronyTerm { g_i, tau_i }` | One Maxwell element: relative modulus `gᵢ = Gᵢ/G0`, relaxation time `τᵢ` (s) |
| `ViscoelasticMaterial { g0, prony }` | Glass shear modulus `G0` (MPa) plus the series |
| `ViscoelasticMaterial::validate() -> Result<(), String>` | Enforces `G0 > 0`, `gᵢ ≥ 0`, `Σgᵢ ≤ 1`, `τᵢ > 0`; returns a descriptive error |
| `::equilibrium_modulus() -> f64` | `G∞ = G0(1 − Σgᵢ)` (MPa) |
| `::relaxation_modulus(time) -> f64` | `G(t)` (MPa) |
| `::storage_modulus(omega) -> f64` | `G'(ω)` (MPa) |
| `::loss_modulus(omega) -> f64` | `G''(ω)` (MPa) |
| `::loss_tangent(omega) -> f64` | `G''/G'` |
| `::step_strain_stress(gamma0, time) -> f64` | Stress after a step strain `γ₀` held for `time` |
| `::elastic_reference() -> TissueModel` | The `tpt-med-tissue` model supplying the glass response |
| `tpt_med_tissue_link` | Re-export shim mapping this crate's names onto `tpt-med-tissue` types |

## Verification

Viscoelasticity has unusually strong analytic structure, and every limit in
it is asserted:

- **Short-time:** `G(0) = G0` and `G'(∞) → G0`, `G''(∞) → 0`.
- **Long-time:** `G(t) → G∞` as `t → ∞`, and `G'(0) = G∞`, `G''(0) = 0`.
- **Monotonicity:** `G(t)` is strictly decreasing and bounded below by `G∞`;
  `G'(ω)` is non-decreasing and bounded above by `G0`.
- **Step-strain stress** decays from `G0·γ₀` toward `G∞·γ₀` and is monotone
  in time.
- **Loss tangent** is asserted equal to `G''/G'` and strictly positive for
  any non-degenerate series.
- **Linear scaling:** multiplying `g0` by `k` multiplies every modulus output
  by `k`, while the loss tangent is *invariant* — a clean test that the series
  is truly relative.
- **Validation:** `Σgᵢ > 1`, negative `gᵢ`, non-positive `τᵢ` and `G0 ≤ 0` are
  each rejected with the specific descriptive message.
- **Single-term reduction:** a one-term series is compared against the exact
  closed-form Maxwell relaxation `G₀ − (G₀−G∞)exp(−t/τ)`.

## Known Limitations

- Shear only; the volumetric response comes from the elastic reference model
  and is not itself rate-dependent.
- Linear viscoelasticity — large-strain, finite-strain and nonlinear
  hyperviscoelastic formulations (e.g. Prony on the Cauchy–Green strain in the
  hyperelastic energy) are not implemented.
- Temperature dependence: `τᵢ` values are supplied by the caller, not
  shifted automatically with a WLF or Arrhenius relation.
- No time-integration helper for a finite-element inner loop; this crate
  provides the material response, and time integration is the solver's job.

## Related Crates

- [`tpt-med-tissue`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-tissue) — the hyperelastic glass response returned by `elastic_reference()`.
- [`tpt-med-cartilage`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-cartilage) — poroelastic time dependence; the two mechanisms are physically distinct and complementary.
- [`tpt-med-biomechanics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-biomechanics) — the linear solver; viscoelastic time stepping is a caller-side loop.
- [`tpt-med-hemodynamics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/fluid/tpt-med-hemodynamics) — blood rheology, the fluid-side counterpart.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). New rheological
models require an [RFC](../../../rfcs). Cite the source of every fitted
parameter set in a doc comment.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.

- Equilibrium modulus `G∞ = G0(1 − Σgᵢ)`, so a fully relaxed state has
  `G∞ = 0` when `Σgᵢ = 1` and a solid response when `Σgᵢ = 0`.
- **Shear only.** The series is defined in shear; the elastic reference model
  supplies the volumetric response.

