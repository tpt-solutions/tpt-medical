# Changelog

All notable changes to `tpt-med-orthopedics` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- **`ElasticHalfSpace` — the continuum step of the interface-compliance
  ladder**: the bone bed as an elastic half-space under a rigid circular
  punch, with the classical Boussinesq settlement
  `δ = F(1−ν²)/(2aE)` (Johnson, *Contact Mechanics*, 1985) and stiffness
  `K = 2aE/(1−ν²)` set by the bone's own modulus instead of a
  caller-invented foundation constant. Uniformises into
  `InterfaceModel`'s N/mm³ units (`K/A`, documented as a screening
  reduction) and composes in series with an implant interface layer —
  where a stiff Ti-alloy layer is nearly transparent (the continuum
  dominates) and a soft cement mantle bites exactly as the series law
  says, both asserted against hand-computed values alongside the exact
  closed form, the settlement's linearity in load, and the physical
  1/E settlement scaling on halved bone modulus.
- `CompliantImplant` (+ `effective_stiffness` / `interface_model`): the
  **screening step from rigid punch toward compliant implant** — the
  Winkler foundation in series with the implant's interface-layer
  stiffness (`k_i = E/t`), so a porous coating or cement mantle softens
  the interface and raises the micromotion a caller screens. The rigid
  punch is the stiff-implant limit and is reproduced unchanged; full
  continuum coupling remains the fem-adapter's job.
- `GruenZone`, `gruen_zone`, `gruen_zones_in_order`: **built-in Gruen
  zone definitions** (Gruen, McNeice & Amstutz 1979) — geometric
  classification of a point into zones 1-7 from the stem axis end points,
  a mediolateral direction and a tip-band fraction, with stable
  `gruen-N` keys for reports.
- Crate README explaining that micromotion and stress shielding are properties
  of the implant–bone *system* rather than of the implant alone, which is why
  the interface is a first-class input instead of an assumed perfect bond.

- `GaitCycle` + `micromotion_over_cycle`: **cyclic loading** — micromotion
  is evaluated across a full load cycle rather than at one static load, and
  both the peak and the per-zone motion *amplitude* (`max − min`) are
  reported, which is the quantity cyclic fibrous-tissue screening keys on.
  `GaitCycle::iso_double_hump` supplies an ISO 14243-style double-hump axial
  load profile (heel-strike and push-off peaks, ≈2.6 × body weight) as a
  screening waveform.
- `MigrationModel`: the **time-dependent consequence of micromotion** — a
  closed-form logarithmic migration law `x(N) = x_bed·ln(1 +
  k(δ_amp − δ_th)N/x_bed)`: per-cycle migration proportional to the motion
  amplitude above a stability threshold, decaying exponentially as the
  implant beds in. `velocity_per_year` / `is_at_risk` implement the
  RSA-style > 0.2 mm/year continued-migration flag.
- Four new tests: double-hump shape (two humps, trough, swing unload,
  periodicity, peak 2.6 × BW), cyclic peak equal to the static result at
  the peak load with positive per-zone amplitudes, closed form vs numerical
  integration of the rate law, and the velocity decay/at-risk boundary/log
  growth shape.

### Planned
- Full continuum coupling (the series-compliance screening and the
  Boussinesq half-space punch are delivered; a genuinely flexible stem in
  a finite bone geometry needs the fem-adapter).

### Notes
- `zone_micromotion` reports **`NaN` for zero-area zones** deliberately, to
  distinguish "this zone does not exist for this implant" from "this zone has
  zero micromotion", which a `0.0` would conflate. NaN zones are excluded from
  the maximum.
- `stress_shielding_analysis` asserts equal zone counts rather than silently
  truncating to the shorter slice.
- The classification thresholds (0.050 mm and 0.150 mm) are
  literature-typical values, not device-specific requirements, and the README
  says so.

## [0.1.0] - 2026-09-22

### Added
- **`InterfaceModel { foundation_stiffness, contact_area, friction }`** — the
  Winkler-foundation parameters, with published ranges documented: cortical
  ~2–20 N/mm³, trabecular ~0.2–2, and press-fit titanium friction ~0.4–0.6.
  All three materially change the answer, so all three are explicit inputs
  rather than hidden constants.
- `micromotion_analysis(&InterfaceModel, Force, tangential_fraction,
  zone_areas) -> MicromotionResult` — per-zone relative interface
  displacement under a joint reaction force, modelling bone as a bed of
  distributed springs and the implant as a rigid punch. Interface shear
  arises from the applied load plus frictional resistance.
- `MicromotionResult` — `max_micromotion` (mm), `zone_micromotion` (mm) and
  `risk`.
- `MicromotionResult::classify(max_micromotion) -> RiskLevel` and the ordered
  `RiskLevel` enum: `Low` below 0.050 mm, `Moderate` below 0.150 mm, `High`
  above. These correspond to the classic osseointegration threshold
  (< 150 µm) and the stricter primary-stability band for cementless press-fit
  components (< 50 µm).
- `micromotion_um(&MicromotionResult) -> Length` — reporting in micrometres,
  which is how the thresholds are actually published, as a typed `Length`.
- **Stress-shielding assessment** — `stress_shielding_analysis(&intact_sed,
  &implanted_sed) -> StressShieldingResult` over Gruen-style zones, returning
  per-zone `with_implant`, `intact` and the shielding index
  `1 − SED_implant/SED_intact`, clamped to `[0, 1]`.
- `StressShieldingResult::mean_index()` and `::has_resorption_risk()`, the
  latter flagging any zone above 0.7 as severely shielded — the classic
  resorption-risk criterion.
- Typed `Force` and `Length` from `tpt-med-units`; no other dependencies.
- `#![forbid(unsafe_code)]`.

### Verification
- **Threshold boundaries are tested at the boundaries** — `classify` is
  asserted at exactly 0.050 and 0.150 mm as well as on either side, since an
  off-by-one at a clinical cut-off is a misclassification.
- **Ordering** — a stiffer foundation and a larger contact area each reduce
  micromotion, asserted as strict inequalities. These are the two
  design-relevant levers, so their sign is a test target, not an assumption.
- **Friction** — higher friction lowers micromotion, asserted strictly.
- **Linearity in load** — doubling the joint reaction doubles the micromotion,
  pinning the unit conversion inside the model.
- **Zero-area zones** yield `NaN` and are excluded from the maximum, so a
  nonexistent zone cannot become the governing one.
- **Shielding index clamping** asserted for extreme ratios (10× implant energy
  clamps to 0; a zone with no implant load clamps to 1), and a zero intact SED
  yields `0.0` rather than a division by zero.
- Golden dataset `test-data/golden/devices/hip_stem_micromotion.json`.

### Known limitations
- Lumped, not continuum: uniform foundation stiffness and a rigid implant. A
  flexible stem or spatially varying bone is not resolved.
- Zones are supplied by the caller; the zone convention is not imposed.
- No cyclic loading; micromotion is evaluated at a single static load, and the
  loading magnitude is the caller's choice.
- No migration or bone-ingrowth model. `RiskLevel` is a threshold rule, not a
  mechanobiological simulation.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
