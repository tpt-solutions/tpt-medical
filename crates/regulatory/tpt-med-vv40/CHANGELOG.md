# Changelog

All notable changes to `tpt-med-vv40` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- `EvidenceMetric` / `Acceptance` with
  `VerificationActivity::{with_metric, metrics_adequate}`: **structured
  evidence** — numeric metrics with mechanical acceptance bands, checked
  automatically instead of reading `results` as prose.
- `CredibilityAssessment::to_json`: deterministic JSON serialisation of the
  assessment (goals, activities, verdict, unmet goals), so an assessment can
  live inside a submission bundle next to the `tpt-med-fda` package.
- Crate README with the full risk × influence → goal matrix tabulated, and a
  worked example (the femur stance-screening question at medium risk /
  significant influence) showing `evaluate()` returning empty.

### Planned
- Per-component credibility rollup, so a large model assembled from small
  verified parts has a defined composite credibility.
  `VerificationActivity` or `ValidationActivity`, so adequacy can be checked
  mechanically rather than by reading `results` as prose.
- Multi-question assessments with an explicit aggregation rule, for models
  used for several questions of interest.
  submission bundle next to the `tpt-med-fda` package.

### Notes
- **Any change to the matrix bands is behaviour-changing and requires an RFC
  and a V&V re-run.** It invalidates every existing assessment, since an
  assessment that passed yesterday's matrix may fail today's.
- Risk and influence are **self-declared** by the caller. The crate cannot
  judge whether `ModelRisk::Medium` is the honest rating for your question,
  and this is the most likely place for a credibility assessment to go quietly
  wrong.
- The matrix is an **engineering interpolation** of the standard's
  Low/Medium/High structure, not a normative table. For a submission, the
  licensed standard text governs.

## [0.1.0] - 2026-09-22

### Added
- `ModelRisk` — `Low`, `Medium`, `High`; derives `Ord`.
- `ModelInfluence` — `Contributing`, `Significant`, `Direct`; derives `Ord`.
- **The risk-influence matrix** — `CredibilityAssessment::goals()` returning
  `CredibilityGoals { min_verification_types, min_validation_types,
  quantitative_validation_required, independent_code_review_required }`. Bands
  (verification types / validation types):

  | Risk \ Influence | Contributing | Significant | Direct |
  |---|---|---|---|
  | Low | 1 / 0 | 1 / 1 | 2 / 1 |
  | Medium | 1 / 1 | 3 / 2 | 3 / 2 |
  | High | 2 / 1 | 3 / 2 | 4 / 3 |

- **Activity taxonomy** — `VerificationType` (`CodeVerification`,
  `CalculationVerification`, `SensitivityAnalysis`,
  `UncertaintyQuantification`; V&V 40 §5.2) and `ValidationType` (`InVitro`,
  `InVivo`, `ClinicalData`, `BenchmarkModel`; §5.3), each with a stable
  `key()` for reports.
- `AgreementLevel` — `Quantitative`, `Qualitative`, `NotDemonstrated`.
- `VerificationActivity { kind, description, results }` and
  `ValidationActivity { kind, reference, metrics, agreement }`, so evidence is
  recorded with its source and the agreement level actually reached — including
  when that level is `NotDemonstrated`.
- `CredibilityAssessment { question_of_interest, risk, influence,
  verification, validation }` with the question stored **verbatim**; an
  assessment whose question of interest cannot be stated precisely cannot be
  evaluated at all.
- `evaluate() -> Vec<String>` — the specific unmet goals, described rather than
  merely counted, so a gap is actionable. `is_credible()` is the boolean
  summary.
- No dependencies, so an assessment can run in CI, in a WASM build, or inside a
  submission bundle.

### Verification
The matrix itself is the thing under test, since a wrong matrix makes every
downstream credibility claim wrong in a way nobody would notice:
- **Monotonicity** — asserted exhaustively across all 9 cells: raising risk
  with influence fixed never lowers a goal, and raising influence with risk
  fixed never lowers a goal. This is the single most important property of the
  matrix.
- **Goal-band correctness** — every one of the 9 cells asserted against its
  expected counts and flags, so a transcription error in any single cell is
  caught.
- **Coverage counts distinct activity types, not activities** — five
  `CodeVerification` activities satisfy a requirement for one. Asserted
  explicitly, because counting activities would let repeated evidence
  masquerade as broad coverage.
- **Quantitative requirement** — an assessment whose only validation is
  `Qualitative` fails when quantitative agreement is required and passes when
  it is not.
- **`NotDemonstrated` is always a failure**, producing an unmet goal naming
  the reference, regardless of cell. Recorded-but-failed evidence must not read
  as coverage.
- **Independent code review** — required at the `2 / 1` band and above, and
  satisfied by the presence of a `CodeVerification` activity.
- **`key()` stability** — every activity type asserted to produce its own
  distinct stable string, since these are written into reports and downstream
  tooling.
- Golden dataset `test-data/golden/regulatory/vv40_credibility_matrix.json`.

### Known limitations
- Not the whole standard: credibility criteria beyond the matrix, the
  assessment process, and decision logic for when a model is *not* adequate are
  out of scope.
- The matrix is interpolated, not a normative reproduction.
- Risk and influence are self-declared.
- Evidence is free text, so the crate can check that an activity exists but
  not that its evidence is adequate.
- No per-component credibility rollup, and one assessment covers one question
  of interest.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
