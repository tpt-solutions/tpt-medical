# tpt-med-vv40

ASME V&V 40 credibility assessment matrices — machine-checkable
risk × influence → credibility-goal evaluation for computational models.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--vv40-orange)](https://crates.io/crates/tpt-med-vv40)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--vv40-blue)](https://docs.rs/tpt-med-vv40)

| | |
|---|---|
| **Layer** | `regulatory` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | none (`std` only) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

ASME V&V 40-2018 asks a question that is easy to state and easy to answer
badly: *how much verification and validation evidence does this model need?*
The answer depends on two things about the **question of interest**, not about
the model:

- **Model risk** — how bad is a wrong answer? If a wrong answer harms a
  patient or sinks a submission, the risk is high.
- **Model influence** — how much does the model decide the outcome? If the
  output is one input among many, the influence is low.

Risk × influence lands in a cell of a matrix that prescribes **credibility
goals**: how many distinct verification activity types, how many validation
activity types, whether *quantitative* agreement is required, and whether an
independent code review is mandatory.

Every simulation vendor has this matrix somewhere — usually as a slide. The
problem with a slide is that it is not checkable. This crate makes the
assessment **an evaluation, not a paragraph**: declare the question's risk and
influence, attach the activities actually performed, call `evaluate()`, and get
back the specific unmet goals. A reviewer can then argue about the evidence
rather than about whether anyone remembered the matrix.

## Features

- **`ModelRisk`** — `Low`, `Medium`, `High`, ordered.
- **`ModelInfluence`** — `Contributing`, `Significant`, `Direct`, ordered.
- **`VerificationType`** — `CodeVerification`, `CalculationVerification`,
  `SensitivityAnalysis`, `UncertaintyQuantification`, each with a stable
  `key()` for reports (V&V 40 §5.2).
- **`ValidationType`** — `InVitro`, `InVivo`, `ClinicalData`,
  `BenchmarkModel`, each with a stable `key()` (§5.3).
- **`AgreementLevel`** — `Quantitative`, `Qualitative`, `NotDemonstrated`.
- **`VerificationActivity` / `ValidationActivity`** — the evidence, with a
  description, results, reference, metrics and agreement level attached. An
  activity that demonstrated nothing is *recorded as such*, rather than
  quietly omitted.
- **`CredibilityAssessment`** — `question_of_interest` (verbatim), `risk`,
  `influence`, and the activity lists.
- **`goals()`** — the matrix lookup.
- **`evaluate() -> Vec<String>`** — the specific unmet goals, empty when
  credible. `is_credible()` is the boolean summary.
- **`AssessmentRollup`** — an explicit, conservative aggregation of several
  assessments: the per-component credibility rollup for a model assembled
  from verified parts, and the multi-question assessment for one model asked
  several questions. Members must each meet their own goals; the composite
  rating (max of declared and member ratings) demands pooled evidence.
- No dependencies, so the assessment can run in CI, in a WASM build, or in a
  submission bundle.

## Conventions

- The matrix is **monotone along both axes**: raising risk or influence never
  lowers a requirement. This is asserted as an invariant, and it is the
  property that makes the matrix trustworthy.
- Goal bands, as *verification types / validation types*:

  | Risk \ Influence | Contributing | Significant | Direct |
  |---|---|---|---|
  | **Low** | 1 / 0 | 1 / 1 | 2 / 1 |
  | **Medium** | 1 / 1 | 3 / 2 | 3 / 2 |
  | **High** | 2 / 1 | 3 / 2 | 4 / 3 |

  `(Low, Direct)` and `(High, Contributing)` share the `2 / 1` band.
- **Quantitative validation** is required for the three highest cells:
  `(Medium, Direct)`, `(High, Significant)` and `(High, Direct)`.
- **Independent code review** is required for the `2 / 1` band and above:
  `(Low, Direct)`, `(Medium, Significant)`, `(High, Contributing)` and
  everything above them.
- `evaluate` returns **descriptions of unmet goals**, not just a count, so a
  gap is actionable.
- `question_of_interest` is stored verbatim and is expected to be the actual
  question, not a project name. An assessment whose question of interest
  cannot be stated precisely cannot be evaluated at all.

## Usage

The worked example shipped with the workspace: the femur stance-screening
question, at **medium risk / significant influence**.

```rust
use tpt_med_vv40::{
    AgreementLevel, CredibilityAssessment, ModelInfluence, ModelRisk, ValidationActivity,
    ValidationType, VerificationActivity, VerificationType,
};

fn assessment() -> CredibilityAssessment {
    CredibilityAssessment {
        question_of_interest: "Peak von Mises stress in the femur under a 3x \
                              body-weight stance reaction, for pre-operative \
                              fixation risk screening."
            .into(),
        risk: ModelRisk::Medium,
        influence: ModelInfluence::Significant,
        verification: vec![
            VerificationActivity {
                kind: VerificationType::CodeVerification,
                description: "Analytic stress and FEM stress cross-check.".into(),
                results: "agrees to 1e-6 relative".into(),
            },
            VerificationActivity {
                kind: VerificationType::CalculationVerification,
                description: "Mesh and quadrature refinement study.".into(),
                results: "<2% change at 2x refinement".into(),
            },
            VerificationActivity {
                kind: VerificationType::SensitivityAnalysis,
                description: "Threshold and modulus sensitivity.".into(),
                results: "peak stress varies 18% over 130-300 HU".into(),
            },
        ],
        validation: vec![
            ValidationActivity {
                kind: ValidationType::InVitro,
                reference: "Sawbone femur, ISO 7206 loading".into(),
                metrics: vec!["peak strain".into(), "crest displacement".into()],
                agreement: AgreementLevel::Quantitative,
            },
            ValidationActivity {
                kind: ValidationType::BenchmarkModel,
                reference: "Published CT-FEM corpus".into(),
                metrics: vec!["peak von Mises".into()],
                agreement: AgreementLevel::Qualitative,
            },
        ],
    }
}

fn main() {
    let a = assessment();

    // The matrix: 3 verification types, 2 validation types, quantitative
    // agreement and an independent code review are all required here.
    let g = a.goals();
    assert_eq!(g.min_verification_types, 3);
    assert_eq!(g.min_validation_types, 2);
    assert!(g.quantitative_validation_required);
    assert!(g.independent_code_review_required);

    // All goals are met by the activities above.
    let unmet = a.evaluate();
    assert!(unmet.is_empty(), "unmet goals: {unmet:?}");
    assert!(a.is_credible());
}
```

Escalating the risk of the *same* question raises the bar:

```rust
use tpt_med_vv40::{CredibilityAssessment, ModelInfluence, ModelRisk};

fn main() {
    let mut a = CredibilityAssessment {
        question_of_interest: "Surgical plan resection volume.".into(),
        risk: ModelRisk::Low,
        influence: ModelInfluence::Contributing,
        verification: vec![],
        validation: vec![],
    };

    // Low risk, contributing influence: one verification type, no validation.
    let g = a.goals();
    assert_eq!(g.min_verification_types, 1);
    assert_eq!(g.min_validation_types, 0);
    assert!(!a.evaluate().is_empty()); // no evidence collected at all

    // Same question, high risk and direct influence: far more evidence.
    a.risk = ModelRisk::High;
    a.influence = ModelInfluence::Direct;
    let g = a.goals();
    assert_eq!(g.min_verification_types, 4);
    assert_eq!(g.min_validation_types, 3);
    assert!(g.quantitative_validation_required);
    assert!(g.independent_code_review_required);
}
```

## API Overview

| Item | Purpose |
|---|---|
| `ModelRisk` | `Low`, `Medium`, `High`; derives `Ord` |
| `ModelInfluence` | `Contributing`, `Significant`, `Direct`; derives `Ord` |
| `VerificationType` | `CodeVerification`, `CalculationVerification`, `SensitivityAnalysis`, `UncertaintyQuantification`; `key()` for reports (V&V 40 §5.2) |
| `ValidationType` | `InVitro`, `InVivo`, `ClinicalData`, `BenchmarkModel`; `key()` (§5.3) |
| `AgreementLevel` | `Quantitative`, `Qualitative`, `NotDemonstrated` |
| `VerificationActivity { kind, description, results }` | One completed verification activity and its evidence |
| `ValidationActivity { kind, reference, metrics, agreement }` | One completed validation activity, its reference data and the agreement reached |
| `CredibilityAssessment { question_of_interest, risk, influence, verification, validation }` | The assessment for one question of interest |
| `::goals() -> CredibilityGoals` | The risk × influence matrix lookup |
| `::evaluate() -> Vec<String>` | Descriptions of the unmet goals; empty when credible |
| `::is_credible() -> bool` | Boolean summary of `evaluate()` |
| `CredibilityGoals` | `min_verification_types`, `min_validation_types`, `quantitative_validation_required`, `independent_code_review_required` |
| `AssessmentRollup { name, declared, members }` | Explicit aggregation of several assessments (component rollup or multi-question study) |
| `::composite_risk() / ::composite_influence() / ::composite_goals()` | The composite rating — max of the declared and member ratings — and its goals |
| `::evaluate() -> Vec<String>` | Member-prefixed unmet goals plus composite unmet goals; empty when the rollup is credible |
| `RollupMember { label, assessment }` | One named component or question inside the rollup |

## Verification

The matrix itself is the thing under test, since a wrong matrix makes every
downstream credibility claim wrong in a way nobody would notice:

- **Monotonicity** — asserted exhaustively across all 9 risk × influence
  cells: raising risk with influence fixed never lowers any goal, and raising
  influence with risk fixed never lowers any goal. This is the single most
  important property of the matrix.
- **Goal-band correctness** — every one of the 9 cells is asserted against its
  expected verification/validation counts and flags, so a transcription error
  in any one cell is caught.
- **Coverage semantics** — `evaluate` counts **distinct activity types**, not
  activities. Five `CodeVerification` activities satisfy a requirement for one,
  and this is asserted explicitly, because counting activities would let a
  repeated evidence type masquerade as broad coverage.
- **Quantitative requirement** — an assessment whose only validation is
  `Qualitative` is asserted to fail when quantitative agreement is required,
  and to pass when it is not.
- **NotDemonstrated is always a failure** — a validation activity with
  `AgreementLevel::NotDemonstrated` produces an unmet goal naming the
  reference, regardless of the risk/influence cell. Recorded-but-failed
  evidence must not read as coverage.
- **Independent code review** — asserted required at the `2/1` band and above,
  and satisfied by the presence of a `CodeVerification` activity.
- **`key()` stability** — every `VerificationType` and `ValidationType` is
  asserted to produce its own distinct stable string, since these are written
  into reports and downstream tooling.
- Worked example and golden dataset:
  `test-data/golden/regulatory/vv40_credibility_matrix.json`, the femur
  stance-screening question at medium risk / significant influence.

- **Rollup aggregation** — a rollup of individually credible members is
  credible at their maximum rating; a member failure is reported under its
  own label; a declared composite rating above the members' maximum raises
  the composite goals and can fail the pooled evidence; a below-member
  declaration is clamped up; an empty rollup is not credible; and the
  composite goals are monotone as members are added.

## Known Limitations

- **Not the whole standard.** V&V 40 also covers credibility criteria beyond
  the risk–influence matrix, the full credibility-assessment process, and
  decision logic for when a model is *not* adequate. This crate encodes the
  matrix and the activity taxonomy only.
- **The matrix is interpolated, not quoted.** The bands are an engineering
  interpretation of the standard's Low/Medium/High structure. It is not a
  reproduction of a normative table, and for a submission the licensed
  standard text governs.
- **Self-declared risk and influence.** The caller supplies the risk and
  influence ratings; the crate cannot judge whether `ModelRisk::Medium` is the
  honest rating for your question. This is the most likely place for a
  credibility assessment to go quietly wrong.
- **Evidence is a string.** `results` and `reference` are free text, so the
  crate can check that an activity *exists* but not that its evidence is
  adequate. Attaching structured metrics to CI is the next step and is not
  done.
- **The rollup rule is fixed and conservative.**
  [`AssessmentRollup`](#api-overview) aggregates by conjunction — every member
  must meet its own goals, and the pooled evidence must meet the composite
  goals derived from the maximum rating. There is no weighted or
  majority-vote aggregation, by design: credibility is not a vote.
- **Composite ratings can only be declared upward.** A declared composite
  risk/influence below a member's is clamped up to the member's; the rollup
  never lets a declaration undercut a part's own assessment.

## Related Crates

- [`tpt-med-fda`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/regulatory/tpt-med-fda) — the audit trail that records the evidence these activities produce.
- [`tpt-med-biomechanics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-biomechanics) and the rest of the `solid` layer — the code whose verification the assessment scopes.
- [`tpt-med-hemodynamics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/fluid/tpt-med-hemodynamics) — a second assessed question of interest.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Any change to the
matrix bands is a **behaviour-changing** change that invalidates existing
assessments and requires an [RFC](../../../rfcs) and a V&V re-run. Cite the
clause of the standard behind every band.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA. This
crate is a structured way to *reason about* credibility; it is not a
compliance determination, and ASME V&V 40's licensed text governs.
