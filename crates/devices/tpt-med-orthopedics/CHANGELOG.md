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
- `GruenZone`, `gruen_zone`, `gruen_zones_in_order`: **built-in Gruen
  zone definitions** (Gruen, McNeice & Amstutz 1979) — geometric
  classification of a point into zones 1-7 from the stem axis end points,
  a mediolateral direction and a tip-band fraction, with stable
  `gruen-N` keys for reports.
- Crate README explaining that micromotion and stress shielding are properties
  of the implant–bone *system* rather than of the implant alone, which is why
  the interface is a first-class input instead of an assumed perfect bond.

### Planned
- Cyclic loading: micromotion accumulated over a gait cycle rather than
  evaluated at a single static load.
- A migration model, so the time-dependent consequence of micromotion can be
  followed rather than classified at a threshold.
  invent one.
- Continuum coupling, so an implant with realistic compliance can be
  evaluated rather than modelled as a rigid punch.

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
