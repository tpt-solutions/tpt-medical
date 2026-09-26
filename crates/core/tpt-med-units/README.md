# tpt-med-units

Type-safe unit system for biomedical simulation — millimetres, megapascals,
newtons, g/cm³ — with no implicit cross-unit arithmetic.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--units-orange)](https://crates.io/crates/tpt-med-units)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--units-blue)](https://docs.rs/tpt-med-units)

| | |
|---|---|
| **Layer** | `core` (leaf — no workspace dependencies) |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | none (`std` only) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Computational biomechanics constantly mixes three unit worlds: SI mechanics
(pascals, newtons), scanner geometry (millimetres), and clinical physiology
(mmHg, mL/min). A factor of 1000 between millimetres and metres, or between
mmHg and Pa, is silent in `f64` code and catastrophic in a femoral implant
sizing study. This crate makes that class of bug a **compile error**.

The design is deliberately minimal: `#[repr(transparent)]` newtypes over `f64`
with `const` constructors, `Add`/`Sub` for same-unit arithmetic, explicit
named conversions for everything else, and **no operator overloading across
different quantities** (there is no `Pressure * Length`; dimensional analysis
is left to the domain crates that actually know the physics).

## Features

- Ten quantity newtypes generated from one macro, so behaviour is identical
  across all of them.
- Canonical units fixed by convention (table below); every conversion is an
  explicit `const fn`.
- `ZERO` / `new` / `from_*` / `to_*` / `value` / `abs` / `max` / `min` /
  `is_finite` on every quantity.
- `UNIT` constant per quantity for display and log formatting.
- `Display` renders as `"12.345000 mm"`, so audit trails and CSV exports carry
  an unambiguous unit suffix.
- `const fn` constructors usable in `const` contexts and statics (this is why
  the workspace MSRV is 1.82).

## Conventions

| Quantity | Canonical unit | Constructor | Common conversions |
|---|---|---|---|
| `Length` | millimetre (`mm`) | `from_mm` / `to_mm` | `from_cm`/`to_cm`, `from_m`/`to_m` |
| `Pressure` | megapascal (`MPa`) | `from_mpa` / `to_mpa` | `from_pa`/`to_pascal`, `from_kpa`/`to_kpa`, `from_mmhg`/`to_mmhg` |
| `Force` | newton (`N`) | `from_n` / `to_n` | `from_kn`/`to_kn`, `body_weights(..)` |
| `Density` | g/cm³ | `from_gcm3` / `to_gcm3` | — |
| `Time` | second (`s`) | `from_s` / `to_s` | `from_ms`/`to_ms`, `from_min` |
| `Angle` | radian (`rad`) | `from_rad` / `to_rad` | — |

## Usage

```rust
use tpt_med_units::{Force, Length, Pressure};

fn main() {
    // Lesion length in millimetres, the workspace convention.
    let lesion = Length::from_mm(30.0);
    assert_eq!(lesion.to_cm(), 3.0);

    // Clinical units in, SI out.
    let systolic = Pressure::from_mmhg(120.0);
    assert!((systolic.to_pascal() - 15_998.7).abs() < 1.0);
    assert!((systolic.to_mpa() - 0.0159987).abs() < 1e-6);

    // ISO 7206-style body-weight loading convention.
    let bw = Force::from_n(80.0 * 9.81);
    let stance = Force::body_weights(3.0, bw);
    assert!((stance.to_n() - 3.0 * bw.to_n()).abs() < 1e-9);

    // Same-unit arithmetic is allowed; cross-unit is not.
    let doubled = lesion + lesion;
    assert_eq!(doubled.to_mm(), 60.0);

    // Display carries the unit.
    assert_eq!(format!("{}", lesion), "30.000000 mm");
}
```

`Force::body_weights` is the one piece of domain logic in this crate: it exists
in exactly one place so the "×N body weight" convention is applied
consistently across the stack.

## API Overview

| Item | Purpose |
|---|---|
| `Length`, `Pressure`, `Force`, `Density`, `Time`, `Angle`, `Viscosity`, `Velocity`, `Modulus`, `FlowRate` | The ten quantity newtypes |
| `::ZERO` | Zero-valued constant for every quantity |
| `::new(f64)`, `::from_*(f64)`, `::to_*() -> f64`, `::value() -> f64` | Construction and extraction |
| `::abs`, `::max`, `::min`, `::is_finite` | Shared numeric helpers |
| `::UNIT -> &'static str` | Unit symbol for display and exports |
| `Add`, `Sub`, `Mul<f64>`, `Div<f64>`, `Div<Self> -> f64` | Same-unit arithmetic |
| `Force::body_weights(f64, Force) -> Force` | Body-weight-multiple loading convention |

All conversions and constructors are `const fn`, so quantities can be built in
`const` items and static tables:

```rust
use tpt_med_units::{Length, Modulus, Pressure};

const CORTICAL_E: Modulus = Modulus::from_mpa(17_500.0);
const CORTICAL_YIELD: Pressure = Pressure::from_mpa(110.0);
const PITCH: Length = Length::from_mm(0.5);

const fn threshold() -> Pressure {
    CORTICAL_E / 10.0
}

fn main() {
    assert_eq!(threshold().to_mpa(), 1_750.0);
    assert_eq!(PITCH.UNIT, "mm");
    assert_eq!(CORTICAL_E.UNIT, "MPa");
}
```

## Verification

- Round-trip conversion tests for every non-trivial conversion, locked against
  the exact constant (`to_pascal` of 120 mmHg to ±1 Pa, and so on).
- The crate is exercised transitively by `tpt-med-biomechanics` (modulus and
  displacement), `tpt-med-hemodynamics` (velocity, viscosity, density),
  `tpt-med-stents` (pressure) and `tpt-med-orthopedics` (length), so a units
  regression fails the workspace test suite rather than a downstream report.

## Known Limitations

- **No dimensional algebra.** There is no `Length * Length -> Area` or
  `Force / Area -> Pressure`. Those expressions are the ones where unit
  mistakes actually cause clinical errors, and this crate deliberately stops
  at preventing *cross-quantity* arithmetic rather than checking dimensional
  consistency. Any crate combining quantities must do that check itself.
- **No uncertainty or tolerance tracking.** A `Length` carries a value and
  nothing else; the metrology needed for a patient-specific claim is not
  modelled.
- **No unit-aware vectors or matrices.** `Vec3` in `tpt-med-geometry` is
  unitless, so a componentwise `Length + Length` on `Vec3`s is not checked.
- **`is_finite` is opt-in.** Arithmetic operators do not validate their
  results, so a `NaN` propagates silently until something reads it. Callers
  that mix measured and computed quantities should check.
- **Conversion constants are `const fn` and not configurable.** Rounding
  constants for a specific unit system (for example, if a study used a
  non-standard mmHg definition) cannot be supplied without rebuilding.
- **`Display` uses six decimal places** regardless of magnitude, so a
  nanometre and a metre both print with the same number of digits. Convenient
  for logs, wrong for a paper.

## Related Crates

- [`tpt-med-geometry`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-geometry) — carries unitless `Vec3`/`Mat3`; the workspace convention is that geometry is fed millimetres.
- [`tpt-med-core`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-core) — patient demographics, including `body_weight_force()`.
- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) — produces `Density`/`Modulus` from Hounsfield Units.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Adding a new
quantity is a semver-minor change, and changing a canonical-unit convention
requires an RFC.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.

| `Viscosity` | Pa·s | `from_pas` / `to_pas` | `from_cp`/`to_cp` |
| `Velocity` | mm/s | `from_mms` / `to_mms` | — |
| `Modulus` | megapascal (`MPa`) | `from_mpa` / `to_mpa` | — |
| `FlowRate` | mm³/s | `from_mm3s` / `to_mm3s` | `from_ml_per_min`/`to_ml_per_min` |

> **Note on `Modulus` vs `Pressure`.** Both are megapascals but they are
> distinct types on purpose. A Young's modulus and a stress are physically
> different quantities, and conflating them is exactly the bug this crate
> exists to prevent. Convert explicitly with `Modulus::new(p.to_mpa())` if you
> truly mean the same number.
