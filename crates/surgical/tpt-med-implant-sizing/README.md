# tpt-med-implant-sizing

Automated implant sizing from anatomical landmarks — chart-driven component
selection with interpolation, plus ISB-style alignment proxies for total knee
arthroplasty.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--implant--sizing-orange)](https://crates.io/crates/tpt-med-implant-sizing)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--implant--sizing-blue)](https://docs.rs/tpt-med-implant-sizing)

| | |
|---|---|
| **Layer** | `surgical` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-geometry`](../../core/tpt-med-geometry) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Choosing the right implant size is a *measurement plus lookup* problem that
surgeons do by eye against a printed chart, and a planning system should do
by arithmetic. The structure of the problem is:

1. Extract reproducible **landmarks** from the patient's anatomy.
2. Compute the **measurements** those landmarks define (transepicondylar
   width, AP depth, plateau width).
3. **Interpolate** through a manufacturer's size chart.

The important design choice is that the **size chart is data, not code**. Each
implant family supplies its own `SizeChart`; the selection logic is shared and
generic. That means a new implant system is a data-entry task, not a software
change — and it means this crate never embeds a vendor's proprietary sizing
table.

A second design choice: `select` returns the **interpolated position**
between the two bracketing sizes, not just a label. A measurement that sits
0.05 of the way from size 4 to size 5 is a genuinely different clinical
situation from one 0.45 of the way, and the surgeon should be able to see
which side of the midpoint they are on. Collapsing that to a bare integer
throws away the information the surgeon uses to decide whether to upsize.

## Features

- **`SizeEntry { label, nominal }`** — one chart row: a manufacturer-specific
  size label and the measurement value it represents (mm).
- **`SizeChart { family, entries }`** — an implant family, ascending by
  `nominal`. `select` picks the **smallest size whose nominal ≥ measurement**.
- **`SizeChoice`** — `label`, the fractional `interpolation` between the
  bracketing sizes, and explicit `below_chart` / `above_chart` out-of-range
  flags.
- **`KneeLandmarks`** — the eight points a TKA plan needs: medial and lateral
  epicondyle, trochlea point, posterior condyle, medial and lateral tibial
  plateau, tibial centre, and tibial tubercle.
- **Derived measurements** — `tea_width` (transepicondylar axis length),
  `ap_depth` (projected onto the TEA-normal so it is independent of the
  measurement plane), and `plateau_width`.
- **Alignment proxies** — `femorotibial_angle_deg` (neutral ≈ 0–10°; larger
  indicates varus/valgus) and `tibial_slope_deg` (clamped to `[0, 30]`).

## Usage

```rust
use tpt_med_geometry::Vec3;
use tpt_med_implant_sizing::{size_tka, KneeLandmarks, SizeChart, SizeEntry};

fn chart(family: &str, base: f64, step: f64) -> SizeChart {
    SizeChart {
        family: family.into(),
        // Ascending by nominal, as `select` requires.
        entries: (1..=8)
            .map(|i| SizeEntry { label: i, nominal: base + step * i as f64 })
            .collect(),
    }
}

fn landmarks() -> KneeLandmarks {
    KneeLandmarks {
        medial_epicondyle: Vec3::new(-35.0, 0.0, 0.0),
        lateral_epicondyle: Vec3::new(35.0, 0.0, 0.0),   // TEA = 70 mm
        trochlea_point: Vec3::new(0.0, 10.0, -60.0),
        posterior_condyle: Vec3::new(0.0, 10.0, 0.0),   // AP depth = 60 mm
        tibial_medial: Vec3::new(-24.0, 10.0, -70.0),
        tibial_lateral: Vec3::new(24.0, 10.0, -70.0),    // plateau = 48 mm
        tibial_center: Vec3::new(0.0, 10.0, -70.0),
        tibial_tubercle: Vec3::new(0.0, -8.0, -85.0),
    }
}

fn main() {
    let femoral = chart("tka-femoral-example", 50.0, 4.0);   // 54..82 mm
    let tibial = chart("tka-tibial-example", 40.0, 3.5);    // 43.5..68 mm
    let lm = landmarks();

    // Individual measurements.
    assert!((lm.tea_width() - 70.0).abs() < 1e-9);
    assert!((lm.plateau_width() - 48.0).abs() < 1e-9);
    assert!(lm.ap_depth() > 0.0);

    // The full bundle.
    let s = size_tka(&lm, &femoral, &tibial).unwrap();
    assert_eq!(s.femoral_size, 6);   // 0.6*70 + 0.4*60 = 66 mm -> chart 66..70
    assert_eq!(s.tibial_size, 2);    // 48 mm
    assert!(s.femorotibial_angle_deg >= 0.0);
    assert!((0.0..=30.0).contains(&s.tibial_slope_deg));

    // Out-of-range measurements are flagged, not silently clamped.
    let tiny = chart("tiny", 500.0, 1.0);
    let choice = tiny.select(10.0).unwrap();
    assert!(choice.below_chart);
    let huge = chart("huge", 1.0, 1.0);
    assert!(huge.select(1000.0).unwrap().above_chart);

    // An empty chart yields None rather than a fabricated size.
    let empty = SizeChart { family: "empty".into(), entries: vec![] };
    assert!(empty.select(60.0).is_none());
}
```

## API Overview

| Item | Purpose |
|---|---|

## Verification

- **Chart interpolation** — a measurement exactly on a chart entry selects that
  entry, and one exactly between two entries yields `interpolation ≈ 0.5`.
- **Tie-breaking** — at exactly `t = 0.5` the **larger** size is selected, and
  the test asserts that specifically. This is a deliberate clinical decision,
  and an arbitrary or unspecified tie-break is a real defect.
- **Monotonicity** — the selected label is non-decreasing as the measurement
  increases across the whole chart, asserted entry by entry.
- **Out-of-range flags** — a measurement below the first entry sets
  `below_chart` and clamps to the first label; above the last, `above_chart`
  and the last label. Both are asserted, and neither silently returns a
  mid-chart size.
- **Empty chart** returns `None`, asserted, so a missing data file cannot
  produce a fabricated recommendation.
- **Measurement geometry** — `tea_width` and `plateau_width` are asserted on
  axis-aligned landmarks with exact expected values; `ap_depth` is asserted
  to be invariant when the landmarks are translated, and to equal the
  TEA-perpendicular component (not the raw point distance), which is the
  property that makes it scan-position independent.
- **Angle clamping** — `tibial_slope_deg` is asserted within `[0, 30]` even
  for pathological landmark configurations, and `femorotibial_angle_deg` is
  asserted finite and non-negative.
- **Blend arithmetic** — the femoral measurement is asserted equal to
  `0.6·TEA + 0.4·AP` for a known landmark set, so the weighting cannot drift.

## Known Limitations

- **Landmark extraction is the caller's job.** This crate measures; it does not
  find the epicondyles. Automatic landmark detection from a CT is the hard
  part of the problem and is not implemented here.
- **Charts are supplied by the caller**, and no vendor chart ships with this
  crate. The `family` field is an identifier, not a lookup into bundled data.
- **Knee only.** Hip, shoulder and ankle sizing are not implemented; the
  `KneeLandmarks` structure is knee-specific by design.
- **Alignment proxies are not surgical targets.** `femorotibial_angle_deg` and
  `tibial_slope_deg` are screening quantities for deformity assessment, not a
  mechanical-alignment target. The real constraint — that the reconstructed
  joint line and the flexion gap match — is not checked.
- **No soft-tissue or ligament balance assessment**, and no check that the
  selected size leaves acceptable gap balancing. Sizing is necessary for a
  good plan and not sufficient.
- **Single measurement per chart.** A real femoral decision considers TEA, AP
  depth and posterior condylar offset jointly, and vendor instructions vary on
  which takes precedence.

## Related Crates

- [`tpt-med-surgical-planning`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/surgical/tpt-med-surgical-planning) — plans the osteotomies this crate then sizes against.
- [`tpt-med-orthopedics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-orthopedics) — evaluates micromotion and stress shielding for the selected component.
- [`tpt-med-wear`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-wear) — bearing-area consequences of the size choice.
- [`tpt-med-core`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-core) — the `Landmark` type for recording landmark provenance.
- [`tpt-med-geometry`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-geometry) — `Vec3` landmark positions.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Do **not** commit
proprietary manufacturer sizing charts; ship them as caller-supplied data.
Changing the femoral blend weights or the tie-break rule is a
behaviour-changing change and requires an [RFC](../../../rfcs).

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. **Not for
clinical decision-making.** Component selection must follow the manufacturer's
instructions for use and the surgeon's judgement.

| `SizeEntry { label, nominal }` | One chart row: a size label and the measurement (mm) it represents |
| `SizeChart { family, entries }` | An implant family; `entries` **must** be ascending by `nominal` |
| `SizeChart::select(measurement) -> Option<SizeChoice>` | Smallest size whose `nominal ≥ measurement`; `None` if the chart is empty |
| `SizeChoice { label, below_chart, above_chart, interpolation }` | Selected label, out-of-range flags, and the fractional position between bracketing sizes |
| `KneeLandmarks` | Eight `Vec3` points: medial/lateral epicondyle, trochlea point, posterior condyle, medial/lateral tibial plateau, tibial centre, tibial tubercle |
| `KneeLandmarks::tea_width() -> f64` | Transepicondylar axis length (mm) |
| `KneeLandmarks::ap_depth() -> f64` | AP depth projected perpendicular to the TEA (mm) |
| `KneeLandmarks::plateau_width() -> f64` | Tibial plateau width (mm) |
| `KneeLandmarks::femorotibial_angle_deg() -> f64` | Alignment proxy; neutral ≈ 0–10° |
| `KneeLandmarks::tibial_slope_deg() -> f64` | Posterior tibial slope, clamped to `[0, 30]` |
| `size_tka(&KneeLandmarks, &SizeChart, &SizeChart) -> Option<TkaSizing>` | Femoral (0.6·TEA + 0.4·AP) and tibial (plateau) sizing plus all measurements and alignment proxies |
| `TkaSizing` | `femoral_size`, `tibial_size`, `tea_width_mm`, `ap_depth_mm`, `plateau_width_mm`, `femorotibial_angle_deg`, `tibial_slope_deg` |
| `Vec3` from `tpt-med-geometry` | Landmark positions in patient space (mm) |

- **`size_tka`** — the bundle: both component labels, all three measurements
  and both alignment proxies in one `TkaSizing` value.
- No dependencies beyond `tpt-med-geometry`.

## Conventions

- All measurements in **millimetres**; angles in **degrees** (as returned by
  the accessors — note this differs from the workspace-wide radian convention
  in `tpt-med-geometry`, because these values are read off a surgeon's
  protractor and printed on a chart).
- `SizeChart.entries` **must be ascending by `nominal`**; `select` relies on
  the ordering.
- **Ties round up.** At exactly `t = 0.5` between two sizes, the larger size
  is selected. A chart step is a manufacturing increment and undersizing
  causes more harm than oversizing in this specific decision.
- Femoral selection uses the blend `0.6 · TEA + 0.4 · AP depth`; tibial
  selection uses plateau width directly. These weights are the workspace
  convention and are stated explicitly rather than buried.
- `ap_depth` is measured **perpendicular to the TEA**, not along a fixed
  anatomical axis, so it is unaffected by how the scan was positioned.
- `size_tka` returns `None` if either chart yields no selection — a caller
  must handle a missing or empty chart explicitly.

