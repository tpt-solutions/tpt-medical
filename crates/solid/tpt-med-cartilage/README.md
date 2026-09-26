# tpt-med-cartilage

Biphasic / poroelastic cartilage models — linear biphasic theory (Mow, Kuei &
Lai 1980) in the confined-compression configuration.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--cartilage-orange)](https://crates.io/crates/tpt-med-cartilage)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--cartilage-blue)](https://docs.rs/tpt-med-cartilage)

| | |
|---|---|
| **Layer** | `solid` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | none (leaf within `solid`) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Articular cartilage is ~70–80 % water by volume. Under load, the interstitial
fluid is forced through the porous solid matrix, and that transport — not the
elasticity of the solid — dominates the time-dependent behaviour. A purely
elastic cartilage model gets the equilibrium response approximately right and
the *timing* completely wrong, which matters directly for contact mechanics
and for how long a joint is protected after impact.

Linear biphasic theory is the standard screening model: a solid matrix with
aggregate modulus `H_A` and Poisson ratio ≈ 0, saturated with fluid moving by
Darcy's law with permeability `k`. Under a step load in confined compression
the surface displacement has a **closed-form series solution**, so this is one
of the few soft-tissue problems that can be verified against analytics rather
than a finite-element code.

## Features

- **`BiphasicMaterial`** — `aggregate_modulus`, `permeability`,
  `poissons_ratio`, `thickness`, with adult articular cartilage defaults.
- **Creep response** — `creep_displacement_fraction(sigma0, time, terms)`
  evaluates the classical series solution, with a caller-controlled term count
  (convergence is a runtime parameter, not a hard-coded truncation).
- **Time landmarks** — `gel_time()` (when the surface displacement exceeds
  50 % of equilibrium, the inverse of the classical gel time),
  `initial_displacement_fraction(sigma0)` and
  `fluid_pressure_fraction(time, terms)` for the fluid-support phase.
- **Equilibrium response** — `equilibrium_strain(sigma0) = σ₀ / H_A`, the
  asymptotic solid-supported limit.
- Pure functions, no allocation, no state — trivially embeddable in a larger
  contact solver or an in-browser loop.

## Conventions

- Aggregate modulus `H_A` in **MPa**; permeability `k` in **mm⁴/(N·s)**;
  thickness `h` in **mm**; `σ₀` in **MPa**; time in **seconds**.

## Usage

```rust
use tpt_med_cartilage::BiphasicMaterial;

fn main() {
    let m = BiphasicMaterial::default(); // H_A = 0.7 MPa, h = 2 mm
    let sigma0 = 0.1;                    // MPa, 100 kPa contact stress

    // Equilibrium: the solid matrix carries the whole load eventually.
    let eq = m.equilibrium_strain(sigma0);
    assert!((eq - sigma0 / m.aggregate_modulus).abs() < 1e-12);

    // Nothing has moved at t = 0 ...
    assert!(m.initial_displacement_fraction(sigma0).abs() < 1e-9);

    // ... and the response is monotone toward equilibrium.
    let early = m.creep_displacement_fraction(sigma0, 1.0, 20);
    let late = m.creep_displacement_fraction(sigma0, 1_000.0, 20);
    assert!(early < late && late <= eq + 1e-9);

    // The fluid-supported phase decays away with time.
    let p_early = m.fluid_pressure_fraction(1.0, 20);
    let p_late = m.fluid_pressure_fraction(1_000.0, 20);
    assert!(p_early > p_late);

    // Gel time: half the equilibrium displacement has occurred.
    let t_gel = m.gel_time();
    let at_gel = m.creep_displacement_fraction(sigma0, t_gel, 200);
    assert!((at_gel - 0.5).abs() < 0.05);

    // Absolute displacement in mm, if that is what you need.
    let u_mm = late * m.thickness * eq;
    assert!(u_mm > 0.0);
}
```

## API Overview

| Item | Purpose |
|---|---|
| `BiphasicMaterial { aggregate_modulus, permeability, poissons_ratio, thickness }` | The material and specimen geometry; `Default` = adult articular cartilage |
| `BiphasicMaterial::creep_displacement_fraction(sigma0, time, terms)` | `u(t)/u_∞` from the confined-compression series solution; `terms` controls convergence |
| `BiphasicMaterial::equilibrium_strain(sigma0)` | `σ₀ / H_A`, the long-time limit |
| `BiphasicMaterial::initial_displacement_fraction(sigma0)` | The `t = 0` value (zero for a step load) |
| `BiphasicMaterial::fluid_pressure_fraction(time, terms)` | Fraction of the load still carried by interstitial fluid |
| `BiphasicMaterial::gel_time()` | Time to 50 % of equilibrium displacement |

## Verification

The series solution has analytic limits, and all of them are asserted:

- **Short-time limit:** `u(t) → 0` as `t → 0`; the fluid carries the full
  load, so `fluid_pressure_fraction(0⁺) → 1`.
- **Long-time limit:** `u(t) → h · σ₀ / H_A` as `t → ∞`, the solid-supported
  equilibrium — compared against the independent `equilibrium_strain`
  computation, not against itself.
- **Gel time:** at `gel_time()` the displacement fraction is 0.5, verified
  against the closed-form `gel_time` rather than by bisection.
- **Monotonicity:** creep is non-decreasing in time and bounded above by
  equilibrium, across a sweep of `σ₀`, `k` and `h`.
- **Permeability sensitivity:** a lower `k` delays but does not change the
  equilibrium — a wrong `H_A` implementation would break this immediately.
- **Term-count convergence:** increasing `terms` moves the result
  monotonically toward the analytic limit, so the truncation is always a
  conservative choice.

## Known Limitations

- **Linear** theory only. Confined compression of cartilage at 100 kPa is
  within the linear regime, but impact loading is not.
- Confined compression is one of several standard boundary conditions;
  unconfined compression and shear are not implemented.
- No nonlinearity in the permeability–strain coupling (the "biphasic
  nonlinear" family), and no lubrication/repulsion term for the contact
  interface.
- The model is 1D through-thickness. A full 3D poroelastic solve is out of
  scope for this crate.

## Related Crates

- [`tpt-med-tissue`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-tissue) — the solid-matrix hyperelastic models, for the non-linear upgrade path.
- [`tpt-med-viscoelastic`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-viscoelastic) — Prony-series time dependence; complementary to (not a substitute for) poroelasticity.
- [`tpt-med-wear`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-wear) — cartilage wear drives the joint-replacement screening case.
- [`tpt-med-orthopedics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-orthopedics) — implant-bone interface analysis in the joint.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Cite the source
of every published constant in a doc comment. New boundary conditions require
an [RFC](../../../rfcs).

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.

- `poissons_ratio ≈ 0` for cartilage (the solid matrix is nearly
  incompressible and confined).
- `Default` values are adult articular cartilage screening values:
  `H_A = 0.7 MPa`, `k = 0.002 mm⁴/(N·s)`, `ν = 0`, `h = 2.0 mm`.
- The returned **fractions are dimensionless** — the fraction of the
  equilibrium displacement, not an absolute length. Multiply by
  `h · σ₀ / H_A` for a displacement in mm.

