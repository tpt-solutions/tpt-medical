# tpt-med-stents

Nitinol superelasticity and stent deployment simulation — a 1D Lagoudas-style
superelastic material plus a stent ring deployment model producing the
ASTM F2394 / F2079 metric family.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--stents-orange)](https://crates.io/crates/tpt-med-stents)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--stents-blue)](https://docs.rs/tpt-med-stents)

| | |
|---|---|
| **Layer** | `devices` |
| **Status** | Alpha, `0.1.0` |
| **Scope** | RFC 0004 — Nitinol superelasticity |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-units`](../../core/tpt-med-units) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

A Nitinol stent is not a spring. It is a **shape-memory alloy**, and that
changes everything about how it deploys: the deployment stress–strain path
has two stress plateaus (forward and reverse transformation) separated by
elastic regions, and the recovery of transformation strain on unloading is
what produces both the chronic outward force and the acute recoil. A linear
elastic model of the same stent gets the peak deployment force roughly right
and the *permanent* force and recoil badly wrong.

This crate implements the physics at two levels:

1. **A 1D superelastic material** with a full hysteretic stress–strain loop —
   austenite at low strain, stress-induced martensite between `σ_ms` and `σ_mf`
   with a cosine transformation-hardening interface, elastic unloading, and
   reverse transformation between `σ_as` and `σ_af`.
2. **A stent ring deployment model** — `N` radial crown springs with that
   material, crimping to store transformation strain, balloon expansion driving
   the ring against an artery modelled as a pressure–diameter tube law.

**Be clear about the fidelity level.** This is a Level-1 ring model. It
computes the right *metrics* with the right *trends*; it is not a 3D
superelastic FEM analysis with frictional contact. That upgrade path is
`tpt-fem-hyperelastic` / `tpt-fem-contact` (RFC 0004), pinned in the workspace
manifest. Default parameters are literature-typical starting points, **not
vendor data**.

## Features

- **`NitinolParams`** — austenite and martensite moduli, transformation strain
  `ε_L`, and the four transformation stresses `σ_ms`, `σ_mf`, `σ_as`, `σ_af`,
  with `Default` at typical ±0.1 mm laser-cut wire, 22 °C body-temperature
  values.
- **`SuperelasticState`** — an explicit hysteretic state machine tracking
  strain, martensite fraction `ξ ∈ [0,1]`, and the active `Branch`:
  `ElasticA` → `Forward` → `ElasticM` → `Reverse`. The branch is what makes
  the loop close: unloading follows a *different* path from loading, which is
  the definition of superelasticity.

## Conventions

- Moduli and stresses in **MPa**; diameters in **mm**; `ε_L ≈ 0.05`
  (dimensionless).
- `crown_stiffness` is **N/mm per crown per mm of radial displacement** in the
  elastic (austenite) regime.
- `vessel_diameter_at_pressure` is a closure `Fn(f64 /*MPa*/) -> f64 /*mm*/`,
  so the vessel compliance law is the caller's — and can be a rigid lumen, a
  linear law `D(p) = lumen + c·p`, or anything else.
- `DeploymentResult.recoil` is a **fraction** (clamped to `[0, 0.2]`), not a
  length. `dogboning` is `|d_end − d_mid| / nominal` and is **`0.0` for the
  uniform ring model** — it is a reported output so that a Level-2 tapered-ring
  model can fill it in without a breaking API change.
- `radial_force` is the **total** force on the vessel (N);
  `contact_pressure` is the mean pressure over the nominal cylindrical
  contact surface (MPa).
- `SuperelasticState` is `Copy` and cheap; the loop is path-dependent, so
  create a fresh state for a fresh load cycle.

## Usage

### The superelastic loop

```rust
use tpt_med_stents::{NitinolParams, SuperelasticState};

fn main() {
    let p = NitinolParams::default();
    let mut s = SuperelasticState::new();

    // Load to 8% strain, then unload fully.
    let mut loading = Vec::new();
    for i in 0..=100 {
        loading.push(s.strain_to_stress(0.08 * i as f64 / 100.0, &p));
    }
    let mut unloading = Vec::new();
    for i in (0..=100).rev() {
        unloading.push(s.strain_to_stress(0.08 * i as f64 / 100.0, &p));
    }

    // The loading plateau spans sigma_ms..sigma_mf.
    let sigma_max = loading.iter().copied().fold(0.0f64, f64::max);
    assert!(sigma_max > p.sigma_mf);

    // Superelasticity: the residual strain is fully recovered, so the state
    // returns to austenite with zero stress at zero strain.
    assert!(s.strain_to_stress(0.0, &p).abs() < 1e-6);
    assert_eq!(s.branch, tpt_med_stents::Branch::ElasticA);

    // Hysteresis: unloading stress is strictly below loading stress.
    assert!(unloading[50] < loading[50]);
}
```

### Ring deployment into a compliant vessel

```rust
use tpt_med_stents::{simulate_deployment, DeploymentResult, NitinolParams, StentModel};
use tpt_med_units::Pressure;

fn main() {
    let stent = StentModel {
        expanded_diameter: 4.0, // mm
        crimped_diameter: 1.2,  // mm
        n_crowns: 8,
        crown_stiffness: 0.9,    // N/mm per crown
    };
    let nitinol = NitinolParams::default();

    // Vessel law: D(p) = 3.6 mm + 0.35 mm/MPa * p.
    let r: DeploymentResult = simulate_deployment(
        &stent,
        &nitinol,
        |p_mpa| 3.6 + 0.35 * p_mpa,
        Pressure::from_mpa(0.013),
    );

    assert!(r.radial_force > 0.0);
    assert!(r.contact_pressure > 0.0);
    assert!((0.0..=0.2).contains(&r.recoil));
    assert!(r.diameter > 3.6); // lumen plus recoil
    println!("d = {:.3} mm, F = {:.2} N, p = {:.4} MPa, recoil = {:.1}%",
             r.diameter, r.radial_force, r.contact_pressure, 100.0 * r.recoil);
}
```

## API Overview

| Item | Purpose |
|---|---|
| `NitinolParams` | `e_austenite`, `e_martensite`, `transformation_strain`, `sigma_ms`, `sigma_mf`, `sigma_as`, `sigma_af`; `Default` = typical laser-cut wire at 22 °C |
| `SuperelasticState` | `strain`, `martensite`, `branch`; `new()`, `strain_to_stress(strain, &NitinolParams)` |
| `Branch` | `ElasticA`, `Forward`, `ElasticM`, `Reverse` — selects the loading or unloading plateau |
| `StentModel` | `expanded_diameter`, `crimped_diameter`, `n_crowns`, `crown_stiffness` |
| `DeploymentResult` | `diameter`, `radial_force`, `contact_pressure`, `recoil`, `dogboning` |
| `simulate_deployment(&StentModel, &NitinolParams, vessel_law, Pressure)` | Radial equilibrium → `DeploymentResult` |
| `Pressure` from `tpt-med-units` | Intraluminal deployment pressure |

## Verification

- **Loop closure** — after loading to 8 % strain and unloading, the stress
  returns to < `1e-6` at zero strain and the state is back to `Branch::ElasticA`.
  A model that does not fully recover residual strain is not superelastic.
- **Hysteresis** — at mid-strain (4 %) the unloading stress is asserted
  **strictly below** the loading stress. This is a strict inequality, not a
  tolerance: a model that collapses the loop loses the physics it exists to
  capture, and an equality-only test would pass for a linear model.
- **Plateau coverage** — the peak loading stress exceeds `σ_mf`, confirming
  both transformation branches were traversed.
- **No-contact case** — an oversized stent in a large vessel returns
  `radial_force == 0.0`, asserted exactly.
- **Sign and clamp invariants** — `recoil` is within `[0, 0.2]`;
  `contact_pressure >= 0`; `diameter > 0`; all outputs are finite.
- **Radial-force linearity** — with a fixed vessel law, doubling the crown
  count doubles the radial force, which pins the force assembly.
- **Parameter sensitivity** — a stiffer crown set yields a larger
  `contact_pressure` against the same lumen, asserted as a strict ordering.
- Golden reference dataset: `test-data/golden/devices/stent_expansion.json`.
  The ASTM F2394 radial-stiffness benchmark scaffold and its literature band
  are in place; the Level-3 FEM correlation remains pending.

## Known Limitations

- **Level-1 ring model.** Uniform ring, uniform stent, no taper, no
  per-segment variation, no foreshortening, no 3D bending stiffness. The
  `dogboning` output is structurally `0.0` here.
- **No friction or contact mechanics.** Crowns are independent radial
  springs; there is no frictional interface, so wall contact is smooth.
- **1D material.** The superelastic model is uniaxial; it does not resolve
  multi-axial transformation, bending stiffness, or the Bauschinger effect.
  A 3D superelastic FEM formulation is the `tpt-fem-hyperelastic` upgrade
  path (RFC 0004).
- **No fatigue or wire-t Fracture prediction** — cyclic degradation of
  `ε_L` over 10⁶ cycles is not modelled, so this cannot be used for
  accelerated-dilation life claims.
- **No vessel wall compliance beyond the supplied diameter law**, and no
  coupling back to a `tpt-med-hemodynamics` flow solution in the same solve.

## Related Crates

- [`tpt-med-hemodynamics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/fluid/tpt-med-hemodynamics) — post-deployment flow, WSS and OSI in the treated vessel.
- [`tpt-med-cardiovascular`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/fluid/tpt-med-cardiovascular) — Windkessel boundary conditions and the coronary waveform.
- [`tpt-med-units`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-units) — the `Pressure` type for deployment pressure.
- [`tpt-med-tissue`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-tissue) — the HGO arterial wall model for the vessel being treated.
- [`tpt-med-wasm`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-wasm) — `wasm_deploy_stent` runs this in the browser and in the `<tpt-stent-simulator>` web component.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Material
parameters must cite their source. New fidelity levels follow the ladder in
`rfcs/0004-nitinol-superelasticity.md` and require an [RFC](../../../rfcs).

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. Deployment
metrics computed here are research estimates and are not a substitute for
ASTM F2394/F2079 bench testing.


An oversized stent in a large vessel produces **no contact** and therefore zero
radial force — asserted in the test suite, because silently returning a
positive force for a stent that never touches the artery is the worst possible
failure mode for this calculation.


- **`strain_to_stress(strain, &NitinolParams)`** — path-dependent: the same
  strain gives a different stress depending on the current branch.
- **`StentModel`** — expanded diameter, crimped diameter, crown count, and
  per-crown radial stiffness.
- **`simulate_deployment`** — solves for radial equilibrium given a vessel
  pressure–diameter law and an intraluminal `Pressure`.
- **`DeploymentResult`** — equilibrium `diameter`, `radial_force`,
  `contact_pressure`, acute `recoil` fraction, and `dogboning`.
