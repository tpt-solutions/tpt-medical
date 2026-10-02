# Changelog

All notable changes to `tpt-med-implant-sizing` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- **Ligament balance proper: `LigamentModel` + `assess_ligament_balance`**
  — the soft-tissue half the gap-balance check deliberately left out.
  Each collateral is a 1-D spring with a slack range (origin, insertion,
  slack length, stiffness, mediolateral side); the screen pairs every
  ligament with the planned positions of its two attachments (what
  resection depths, component thicknesses and fragment moves do to them)
  and reports anatomical vs planned length, elongation, spring tension,
  per-ligament tautness, and the mediolateral imbalance against a
  tolerance. The documented trap is asserted: an over-released side
  balances numerically while slack, so `is_balanced` must be read
  together with the `is_taut` flags. Verified: exact zero through slack
  and exact linearity beyond, exact planned-length geometry on moved
  attachments, side-sum imbalance arithmetic, and the over-release trap.
- **Shoulder stem sizing: `ShoulderLandmarks` + `size_shoulder`** — the
  humeral landmarks (head centre, anatomical neck, isthmus edges, canal
  entry, elbow epicondyles) and the measurements they define: endosteal
  canal width at the isthmus (drives the stem size), head offset
  (head-centre-to-canal-axis perpendicular distance, the standard
  geometric definition), and a head-height axial projection (a proxy for
  the calcar-to-apex head height — documented as a proxy, since the
  clinical definition references the calcar and the head apex).
  `ShoulderSizing` reports the precedence-resolved stem decision plus
  two alignment outputs, neck–shaft angle (normal ≈ 130–140°) and
  humeral retroversion (the in-plane neck-axis/TEA angle, normal ≈
  20–30°); the alignment quantities are diagnostics, not sizing inputs.
  The caller-supplied stem chart feeds the same precedence-resolved
  multi-measurement machinery as the knee and hip paths; no clinical
  chart data ships with the crate. Verified: exact measurements on
  hand-placed landmarks, a neck axis built at 135° to the shaft reading
  exactly 135°, a 25° in-plane rotation reading exactly 25° retroversion
  without disturbing any other measurement, precedence resolution on
  disagreement, and monotone sizing in canal width.
- **Ankle component sizing: `AnkleLandmarks` + `size_ankle`** — the
  tibial plafond edges (medial/lateral, anterior/posterior), the talar
  dome edges and centre, and a proximal tibial axis point. The tibial
  decision weighs plafond width and depth (width takes precedence), the
  talar decision maps dome width, and the tibiotalar deviation — the
  angle between the tibial anatomical axis and the dome line's
  perpendicular, 0° at neutral — is reported alongside in `AnkleSizing`.
  Verified: exact measurements on axis-aligned landmarks, a 30° axis
  tilt constructed on a `sin 30° = 0.5` triangle reading exactly 30°,
  precedence resolution with the out-of-chart vote recorded, and
  monotone tibial sizing in plafond width.

### Added
- `ResectionPlan` / `GapReport` / `check_gap_balance`: the **geometric
  half of soft-tissue assessment** — extension and flexion gaps from
  resection-vs-component-thickness arithmetic, flagged for overstuffed
  components (negative gap) and flexion/extension imbalance beyond a
  tolerance. Ligament tension itself stays out (needs soft-tissue
  structures). Three-assertion test covering balanced, overstuffed and
  flexion-open plans.
- `size_from_measurements` + `MeasurementInput` / `MultiMeasurementDecision`
  / `Resolution` / `MeasurementVote`: **multi-measurement sizing with an
  explicit vendor precedence rule** — each measurement maps through the
  chart, every vote is recorded (with out-of-chart flags), and
  disagreement resolves by lowest precedence rank, with equal ranks
  up-sizing conservatively. `KneeLandmarks::posterior_condylar_offset`
  (a TEA-perpendicular screening proxy, documented as such) and
  `femoral_measurements()` (TEA → AP depth → PCO precedence) supply the
  femoral triple. Four new tests: the PCO proxy against a hand-computed
  offset, consensus, precedence resolution, and the up-size tie-break
  with out-of-chart recording.
- `SizeChart::validate` (+ `ChartError`): **schema validation for
  caller-supplied charts** — non-empty, finite positive nominals in strict
  ascending order, unique labels — so a mis-transcribed chart is caught at
  load rather than silently producing recommendations.
- Crate README documenting the two design decisions that matter: the size chart
  is **data, not code** (so a new implant system is data entry, and no
  vendor's proprietary table is embedded), and `select` returns the
  **interpolated position** rather than only a label, because a measurement
  0.05 of the way from size 4 to size 5 is a different clinical situation from
  one 0.45 of the way.

### Added
- **Hip stem sizing: `HipLandmarks` + `size_hip`** — the femoral-side
  landmarks (head centre, lesser trochanter, isthmus edges, canal entry)
  and the three geometric measurements they define: endosteal canal
  width at the isthmus (drives the distal size), femoral offset
  (head-centre-to-canal-axis perpendicular distance, the standard
  geometric definition), and a head-to-lesser-trochanter axial
  projection (a femoral-side leg-length proxy — documented as a proxy,
  since the clinical definition references a pelvis landmark). The
  caller-supplied stem chart feeds the same precedence-resolved
  multi-measurement machinery as the knee path; no clinical chart data
  ships with the crate. Verified: exact geometric measurements on
  hand-placed landmarks, precedence resolution on disagreement, and
  monotone sizing in canal width.

### Planned
- Automatic landmark detection from a CT, which is the hard part of the
  problem and is not attempted here.
- Stability assessment beyond static tension (varus/valgus stress
  response, pivot shift) — needs kinematics this crate does not model.

### Notes
- **Ties round up.** At exactly `t = 0.5` between two sizes the larger size is
  selected, and this is asserted specifically. A chart step is a manufacturing
  increment, and undersizing causes more harm than oversizing in this decision.
  Changing the tie-break is a behaviour-changing change requiring an RFC.
- The femoral blend `0.6·TEA + 0.4·AP depth` is a **workspace convention**,
  stated explicitly rather than buried, and is asserted in the tests so it
  cannot drift.
- Do **not** commit proprietary manufacturer sizing charts; they are
  caller-supplied data by design.

## [0.1.0] - 2026-09-22

### Added
- `SizeEntry { label, nominal }` — one chart row: a manufacturer-specific size
  label and the measurement value in millimetres it represents (typically
  transepicondylar width or AP depth for femoral charts, plateau width for
  tibial).
- `SizeChart { family, entries }` — an implant family with entries **ascending
  by `nominal`**, which `select` relies on. `select(measurement)` returns the
  smallest size whose `nominal ≥ measurement`, and `None` for an empty chart.
- `SizeChoice { label, below_chart, above_chart, interpolation }` — the
  selected label, explicit out-of-range flags, and the fractional position
  between the bracketing sizes so the caller can decide to size up or down.
- `KneeLandmarks` — the eight `Vec3` points a TKA plan needs: medial and
  lateral epicondyle, trochlea point, posterior condyle, medial and lateral
  tibial plateau, tibial centre, and tibial tubercle.
- **Derived measurements** — `tea_width` (transepicondylar axis length),
  `ap_depth` (projected onto the TEA-normal so it is independent of how the
  scan was positioned), and `plateau_width`.
- **Alignment proxies** — `femorotibial_angle_deg` (neutral ≈ 0–10°; larger
  indicates varus/valgus deformity) and `tibial_slope_deg` (clamped to
  `[0, 30]`).
- `size_tka(&KneeLandmarks, &SizeChart, &SizeChart) -> Option<TkaSizing>` —
  femoral sizing on the `0.6·TEA + 0.4·AP` blend, tibial on plateau width,
  with both component labels, all three measurements and both alignment
  proxies in one `TkaSizing` bundle. `None` if either chart yields no
  selection.
- No dependencies beyond `tpt-med-geometry`; `#![forbid(unsafe_code)]`.

### Verification
- **Chart interpolation** — a measurement exactly on a chart entry selects that
  entry, and one exactly between two entries yields `interpolation ≈ 0.5`.
- **Tie-breaking** — at exactly `t = 0.5` the larger size is selected, asserted
  specifically, because an arbitrary or unspecified tie-break is a real
  defect.
- **Monotonicity** — the selected label is non-decreasing as the measurement
  increases across the whole chart, asserted entry by entry.
- **Out-of-range flags** — below the first entry sets `below_chart` and clamps
  to the first label; above the last, `above_chart` and the last label. Neither
  silently returns a mid-chart size.
- **Empty chart** returns `None`, so a missing data file cannot produce a
  fabricated recommendation.
- **Measurement geometry** — `tea_width` and `plateau_width` asserted on
  axis-aligned landmarks with exact expected values; `ap_depth` asserted
  invariant under translation and equal to the TEA-perpendicular component
  rather than the raw point distance.
- **Angle clamping** — `tibial_slope_deg` within `[0, 30]` even for
  pathological landmarks, and `femorotibial_angle_deg` finite and
  non-negative.
- **Blend arithmetic** — the femoral measurement asserted equal to
  `0.6·TEA + 0.4·AP`, so the weighting cannot drift.

### Known limitations
- Landmark extraction is the caller's job; automatic detection from a CT is not
  implemented.
- Charts are caller-supplied and no vendor chart ships with the crate.
- Knee only; `KneeLandmarks` is knee-specific by design.
- Alignment proxies are screening quantities, not surgical targets, and the
  reconstructed joint line and flexion gap are not checked.
- No soft-tissue or ligament balance assessment.
- One measurement per chart; real femoral sizing weighs several and vendor
  instructions differ on precedence.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
