# tpt-med-orthopedics

Implant–bone micromotion and stress-shielding analysis for orthopedic devices.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--orthopedics-orange)](https://crates.io/crates/tpt-med-orthopedics)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--orthopedics-blue)](https://docs.rs/tpt-med-orthopedics)

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

Two numbers decide whether an uncemented orthopedic implant works, and both
are about the *interface* rather than the implant:

- **Micromotion** — how far the implant slides relative to the bone under
  load. Too much and the interface fibroses instead of osseointegrating; too
  little in a primary press-fit and there is no stimulus to bone ingrowth.
  The classic screening thresholds are **< 150 µm** for bone fixation and a
  stricter **< 50 µm** band for primary stability of a cementless press-fit
  component.
- **Stress shielding** — the bone next to a stiff implant carries far less
  strain energy than it did when intact, and bone that is not loaded
  resorbs. This is the mechanism behind late aseptic loosening.

Neither is a property of the implant alone; both are properties of the
implant–bone *system*. This crate makes that explicit by taking the interface
as a first-class input rather than assuming perfect bonding.

## Features

- **Winkler-foundation interface model** — bone as a bed of distributed
  springs, the implant as a rigid punch; interface shear arises from the
  applied load plus frictional resistance at the interface.
- **`InterfaceModel`** — `foundation_stiffness` (cortical ~2–20 N/mm³,
  trabecular ~0.2–2), `contact_area` (mm²), and `friction` (press-fit
  titanium ~0.4–0.6). All three materially change the answer, so all three
  are explicit parameters.
- **`micromotion_analysis`** — per-zone relative displacement, the maximum,
  and an automatic `RiskLevel` classification.
- **`RiskLevel::{Low, Moderate, High}`** — thresholds `< 0.05 mm`,
  `< 0.15 mm`, else high, via `MicromotionResult::classify`.
- **Gruen-style zone analysis** — `stress_shielding_analysis` takes paired
  per-zone strain energy densities (with implant, and intact baseline) and
  returns the per-zone shielding index `1 − SED_implant/SED_intact`, clamped
  to `[0, 1]`.
- **`has_resorption_risk()`** — flags any zone above 0.7 as severely
  shielded, the classic resorption-risk criterion.
- **`micromotion_um`** — reports micromotion in micrometres, which is how the
  thresholds are actually published, using a typed `Length`.
- Zero dependencies beyond `tpt-med-units`.

## Conventions

- `foundation_stiffness` in **N/mm³**; `contact_area` in **mm²**;
  `friction` dimensionless; joint reaction as a typed `Force` in **N**.
- **Micromotion is returned in millimetres** (`max_micromotion`), and

## Usage

```rust
use tpt_med_orthopedics::{
    micromotion_analysis, micromotion_um, stress_shielding_analysis, InterfaceModel,
    MicromotionResult,
};
use tpt_med_units::Force;

fn main() {
    // Press-fit stem in cortical bone.
    let interface = InterfaceModel {
        foundation_stiffness: 8.0, // N/mm^3 (cortical ~2-20)
        contact_area: 400.0,      // mm^2
        friction: 0.5,            // press-fit titanium ~0.4-0.6
    };

    // 3x body weight joint reaction, 30% tangential.
    let reaction = Force::from_n(3.0 * 80.0 * 9.81);
    let zones = [400.0, 400.0, 200.0, 0.0]; // mm^2 per Gruen zone
    let m = micromotion_analysis(&interface, reaction, 0.30, &zones);

    // The zero-area zone is reported as NaN and excluded from the maximum.
    assert!(m.zone_micromotion[3].is_nan());
    assert!(m.max_micromotion > 0.0);
    assert_eq!(m.risk, MicromotionResult::classify(m.max_micromotion));
    println!("max micromotion: {:.1} um ({:?})",
             micromotion_um(&m).to_mm(), m.risk);

    // Stress shielding, per zone, against the intact bone.
    let intact    = [1.0, 1.0, 0.8, 0.5, 0.3, 0.2, 0.1]; // MPa = mJ/mm^3
    let implanted = [0.1, 0.2, 0.6, 0.5, 0.3, 0.2, 0.1];
    let s = stress_shielding_analysis(&intact, &implanted);
    assert_eq!(s.shielding_index.len(), intact.len());
    assert!(s.mean_index() > 0.0);
    assert!(s.has_resorption_risk(), "proximal zones are severely shielded");
}
```


## API Overview

| Item | Purpose |
|---|---|
| `InterfaceModel { foundation_stiffness, contact_area, friction }` | Winkler foundation parameters; cortical ~2–20 N/mm³, trabecular ~0.2–2, friction ~0.4–0.6 |
| `micromotion_analysis(&InterfaceModel, Force, tangential_fraction, zone_areas) -> MicromotionResult` | Per-zone and maximum interface motion |
| `MicromotionResult` | `max_micromotion` (mm), `zone_micromotion` (mm, `NaN` for zero-area zones), `risk` |
| `MicromotionResult::classify(max_micromotion) -> RiskLevel` | `< 0.05` Low, `< 0.15` Moderate, else High (mm) |
| `RiskLevel` | `Low`, `Moderate`, `High`; derives `Ord` for sorting |
| `micromotion_um(&MicromotionResult) -> Length` | Micromotion as a typed `Length` in micrometres |
| `stress_shielding_analysis(&intact_sed, &implanted_sed) -> StressShieldingResult` | Per-zone `1 − SED_implant/SED_intact`, clamped to `[0,1]` |
| `StressShieldingResult` | `with_implant`, `intact`, `shielding_index` |
| `StressShieldingResult::mean_index() -> f64` | Mean index over zones |
| `StressShieldingResult::has_resorption_risk() -> bool` | True if any zone exceeds 0.7 |
| `Force` from `tpt-med-units` | Joint reaction force in newtons |

  thresholds are therefore `0.050 mm` and `0.150 mm`. Use `micromotion_um`
  for the published micrometre form.
- `tangential_fraction` is clamped to `[0, 1]` — it is the fraction of the
  joint reaction acting tangentially to the interface.
- Strain energy density is in **MPa** (= mJ/mm³), so the shielding index is a
  dimensionless ratio.
- Zones with **zero contact area are reported as `NaN` and skipped** when
  taking the maximum. `NaN` here is deliberate: it distinguishes "this zone
  does not exist for this implant" from "this zone has zero micromotion",
  which a `0.0` would conflate.
- `stress_shielding_analysis` asserts equal zone counts rather than
  silently truncating.

## Verification

- **Threshold boundaries are tested at the boundaries.** `classify` is
  asserted at exactly `0.050` and `0.150` mm as well as on either side, since
  an off-by-one at a clinical cut-off is a misclassification.
- **Ordering** — a stiffer foundation and a larger contact area each reduce
  micromotion, asserted as strict inequalities. These are the two
  design-relevant levers, so their sign is a test target, not an assumption.
- **Friction** — higher friction lowers micromotion, asserted strictly.
- **Linearity in load** — doubling the joint reaction doubles the micromotion,
  which pins the unit conversion inside the model.
- **Zero-area zones** are asserted to yield `NaN` and to be excluded from the
  maximum, so a nonexistent zone cannot silently become the governing one.
- **Shielding index** is asserted within `[0,1]` for extreme ratios: an
  implant carrying 10× the intact energy clamps to 0, a zone with no implant
  load clamps to 1.
- **Degenerate baseline** — an intact SED of zero yields index `0.0` rather
  than a division by zero.
- Golden reference dataset:
  `test-data/golden/devices/hip_stem_micromotion.json`.

## Known Limitations

- **Lumped, not continuum.** The Winkler foundation assumes a uniform
  foundation stiffness and a rigid implant. A genuinely flexible stem, or bone
  whose density varies around the implant, is not resolved — pass a
  representative `foundation_stiffness` and treat the result as a screen.
- Zones are supplied by the caller. Defining them is the caller's job, and
  the zone convention (Gruen, Paprosky, or vendor-specific) is not imposed.
- No cyclic or fatigue loading: micromotion is evaluated at a single static
  load, not accumulated over a gait cycle. The loading *magnitude* is the
  caller's choice; a peak-cycle value is not derived here.
- No migration or bone-ingrowth model. The `RiskLevel` classification is a
  threshold rule, not a mechanobiological simulation.

## Related Crates

- [`tpt-med-wear`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-wear) — the other long-term failure mode for a joint replacement.
- [`tpt-med-bone`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-bone) — Wolff's-law remodeling, the time-dependent consequence of stress shielding.
- [`tpt-med-biomechanics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-biomechanics) — supplies the per-zone strain energy densities.
- [`tpt-med-units`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-units) — the `Force` and `Length` types.
- [`tpt-med-implant-sizing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/surgical/tpt-med-implant-sizing) — chooses the component this analysis then evaluates.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Cite the source
of every clinical threshold, and state explicitly whether it is a
literature-typical value or a device-specific requirement.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. Not a
diagnostic device and not a substitute for in-vitro or clinical validation of
any implant.
