# tpt-med-wear

Implant wear simulation — Archard and Archard-type Cross–Land laws with
gait-cycle extrapolation, for ISO 14879 / ASTM F2028-style screening.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--wear-orange)](https://crates.io/crates/tpt-med-wear)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--wear-blue)](https://docs.rs/tpt-med-wear)

| | |
|---|---|
| **Layer** | `devices` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-units`](../../core/tpt-med-units) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Polyethylene wear is the primary reason total knee replacements are revised.
The wear simulator in a hip simulator runs for five million cycles over
months of machine time; the design engineer needs an answer on Tuesday, and
the regulatory reviewer needs a number with a documented basis.

Two classical laws cover the engineering question:

- **Archard** — `V = k · F · s`, volumetric wear from contact load `F` and
  sliding distance `s`. Simple, and the `k` coefficient absorbs everything the
  model does not represent. It is the right law when you have measured
  specific wear rates.
- **Cross–Land** — per-cycle wear depth proportional to contact pressure
  *above a fatigue threshold* `p₀`: `dh = K·(p − p₀)·ds`. Below the threshold
  there is **no wear at all**, which is the physically important feature: it
  is what makes a well-functioning, low-contact-pressure bearing nearly
  wear-free, and a mis-loaded one fail catastrophically rather than
  proportionally.

Implement both, take the coefficients and the screening limit from the caller,
and the result is traceable to a published method rather than to a proprietary
black box.

## Features

- **`WearLaw::Archard { k }`** — the specific wear rate, in mm³/(N·m),
  multiplied by sliding metres.
- **`WearLaw::CrossLand { k, pressure_threshold }`** — pressure-threshold law
  with an explicit fatigue threshold `p₀` (MPa); wear is exactly zero below
  it.
- **`WearModel { law, gait_cycles }`** — law plus the number of cycles to
  simulate or extrapolate to.
- **Multi-zone input** — `simulate_wear` takes parallel arrays of per-zone
  contact pressures and per-cycle sliding distances, so a real bearing with a
  non-uniform contact pattern is representable.
- **`WearResult`** — `volumetric_wear` (mm³), `linear_wear` (mm, mean depth
  over the bearing area), and `wear_per_megacycle` (mm³/Mc), which is the
  unit ISO 14879 expresses limits in.
- **`exceeds_iso14879_screen`** — a direct comparison against a caller-supplied
  limit, so a screening verdict is one call rather than a hand calculation.
- No allocation, trivially embeddable in a parameter sweep.

- **Contact coupling with wear feedback** — `simulate_wear_with_contact`
  drives a `ContactSolver` (the built-in `WinklerContact` foundation, or an
  external solver implementing the trait) once per block: accumulated
  per-zone wear depths change the solved pressures, so the contact
  geometry evolves instead of staying prescribed. Worn zones shed load;
  contact loss past the penetration is flagged.

## Conventions

- Contact pressure in **MPa**; sliding distance in **mm per cycle**; bearing
  area in **mm²**.
- Volumetric wear in **mm³**; linear wear depth in **mm**;
  `wear_per_megacycle` in **mm³/Mc** (per million cycles).
- Wear coefficients: `k` in **mm³/(N·m)** for Archard, and the same units for
  Cross–Land `K`, with `p₀` in **MPa**.
- Wear coefficients span several orders of magnitude between implant materials,
  bearing designs and test protocols. They are **inputs, never defaults here** —
  a library-level default would be a fabricated number wearing a citation.

## Usage

```rust
use tpt_med_wear::{WearLaw, WearModel};

fn main() {
    // Archard, specific wear rate from a hip-simulator protocol.
    let archard = WearModel {
        law: WearLaw::Archard { k: 1.0e-8 }, // mm^3/(N*m)
        gait_cycles: 1_000_000,
    };
    // Two zones: 5 MPa / 3 mm per cycle, and 12 MPa / 1 mm per cycle.
    let r = archard.simulate_wear(&[5.0, 12.0], &[3.0, 1.0], 800.0);
    assert!(r.volumetric_wear > 0.0);
    assert!((r.linear_wear - r.volumetric_wear / 800.0).abs() < 1e-12);

    // The ISO 14879 screening verdict, against a caller-supplied limit.
    let limit = 30.0; // mm^3/Mc, commonly quoted for UHMWPE tibial inserts
    let fails = archard.exceeds_iso14879_screen(&r, limit);
    println!("{:.2} mm^3/Mc -> over limit: {}", r.wear_per_megacycle, fails);

    // Cross-Land: a threshold pressure, so a lightly loaded zone is wear-free.
    let crossland = WearModel {
        law: WearLaw::CrossLand { k: 1.0e-9, pressure_threshold: 8.0 },
        gait_cycles: 1_000_000,
    };
    // 5 MPa is below p0 = 8 MPa -> zero wear. 20 MPa is above -> wear.
    let below = crossland.simulate_wear(&[5.0], &[3.0], 800.0);
    let above = crossland.simulate_wear(&[20.0], &[3.0], 800.0);
    assert_eq!(below.volumetric_wear, 0.0);
    assert!(above.volumetric_wear > 0.0);
}
```

## API Overview

| Item | Purpose |
|---|---|
| `WearLaw::Archard { k }` | `V = k · F · s`; `k` is the specific wear rate in mm³/(N·m) |
| `WearLaw::CrossLand { k, pressure_threshold }` | `dh = K·(p − p₀)·ds`; `p₀` in MPa, wear is exactly zero at or below it |
| `WearModel { law, gait_cycles }` | Law selection plus the cycle count to simulate or extrapolate to |
| `WearModel::simulate_wear(&contact_pressures, &sliding_distances, area) -> WearResult` | Multi-zone accumulation over `gait_cycles`; `area` in mm² converts volumetric to linear wear |
| `WearResult` | `volumetric_wear` (mm³), `linear_wear` (mm), `wear_per_megacycle` (mm³/Mc) |
| `exceeds_iso14879_screen(&self, &WearResult, limit_mm3_per_mc) -> bool` | Screening verdict against a caller-supplied limit |

## Verification

Both laws have exact algebraic forms, so they are verified against
hand-computed values rather than snapshots:

- **Archard closed form** — `simulate_wear` is compared against
  `k · Σ(μp_i · A · s_i · N)` computed independently, to a tight relative
  tolerance. This pins the pressure-to-force conversion (`F = p · A`) and the
  metre/millimetre conversion in the same test.
- **Cross–Land threshold behaviour** — a zone *exactly at* `p₀` is asserted to
  produce **zero** wear, and a zone just above it produces a positive amount
  that vanishes as `p → p₀⁺`. This is the defining property of the law.
- **Cross-Land above threshold** reduces to the Archard form with
  `k' = k(1 − p₀/p)`, asserted against the algebraic identity.
- **Unit conversions** — one million cycles is asserted to produce
  `wear_per_megacycle == volumetric_wear` when `gait_cycles == 1_000_000`,
  and the mm³/mm and mm³/Mc units are each pinned by an independent
  calculation.
- **Linear wear** is asserted to equal `volumetric_wear / area` exactly.
- **Multi-zone vs. single-zone consistency** — splitting one zone into two
  with proportional areas reproduces the same total.
- **Monotonicity** — wear is non-decreasing in `gait_cycles`, contact pressure
  and sliding distance. This is a strict physical requirement and a cheap
  detector for sign errors.
- Golden reference dataset: `test-data/golden/devices/knee_wear_10mcycles.json`.

## Known Limitations

- **No wear-debris-induced damage feedback.** Wear does not change the contact
  geometry or pressures, so the runaway that ends real implant life is not
  captured. This is a screening tool, not a life model.
- **No lubricant or contact-mechanics solution.** Pressures and sliding
  distances are inputs, not solved quantities. A wear simulation coupled to a
  contact solver is a caller-side assembly.
- **Steady-state, single-cycle-condition extrapolation.** The law is linear in
  `gait_cycles`; a run-in period, varying activity level, or a change in gait
  over implant life is not represented.
- **No wear coefficient uncertainty propagation.** Coefficients in this domain
  scatter over orders of magnitude between studies; a defensible screening
  study should run a range, and this crate offers no help in doing that
  systematically.
- Not a substitute for ASTM F2028 hip-simulator or knee-simulator testing.

## Related Crates

- [`tpt-med-orthopedics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-orthopedics) — implant fixation and stress shielding, the other revision driver.
- [`tpt-med-cartilage`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-cartilage) — the articulating tissue whose degeneration drives the wear case.
- [`tpt-med-implant-sizing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/surgical/tpt-med-implant-sizing) — sizing, which determines the bearing area and hence the contact pressures.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Cite the source
of every wear coefficient and screening limit. New wear laws require an
[RFC](../../../rfcs).

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. Not a
substitute for ASTM F2028 or ISO 14879 standard testing.
