//! ASME V&V 40 credibility assessment matrices.
//!
//! ASME V&V 40-2018 frames simulation credibility through the **model risk**
//! (how bad is a wrong answer?) and **model influence** (how much does the
//! model decide the outcome?) of a question of interest. Risk × influence
//! dictates the *credibility goals* — how much verification and validation
//! evidence is required. This crate encodes the risk-influence-credibility
//! matrix and the supporting activity taxonomy so assessments are
//! machine-checkable rather than prose.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// Model risk of the question of interest.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ModelRisk {
    /// A wrong answer is tolerable / easily caught downstream.
    Low,
    /// A wrong answer misleads but is recoverable.
    Medium,
    /// A wrong answer can harm a patient or sink the submission.
    High,
}

/// Model influence on the decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ModelInfluence {
    /// Model output is one input among several.
    Contributing,
    /// Model output substantially determines the decision.
    Significant,
    /// Model output alone decides (or a wrong output alone causes harm).
    Direct,
}

/// Verification activity types (V&V 40 §5.2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerificationType {
    /// Code verification (units, regression, analytical benchmarks).
    CodeVerification,
    /// Calculation verification (grid/time convergence, solver quality).
    CalculationVerification,
    /// Sensitivity analysis of inputs to the question of interest.
    SensitivityAnalysis,
    /// Uncertainty quantification on inputs and outputs.
    UncertaintyQuantification,
}

impl VerificationType {
    /// Stable key for reports.
    pub fn key(self) -> &'static str {
        match self {
            VerificationType::CodeVerification => "code_verification",
            VerificationType::CalculationVerification => "calculation_verification",
            VerificationType::SensitivityAnalysis => "sensitivity",
            VerificationType::UncertaintyQuantification => "uq",
        }
    }
}

/// Validation activity types (V&V 40 §5.3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationType {
    /// Bench/in-vitro experiment against the model.
    InVitro,
    /// In-vivo comparison.
    InVivo,
    /// Published clinical/literature data.
    ClinicalData,
    /// Benchmark against an established reference model.
    BenchmarkModel,
}

impl ValidationType {
    /// Stable key for reports.
    pub fn key(self) -> &'static str {
        match self {
            ValidationType::InVitro => "in_vitro",
            ValidationType::InVivo => "in_vivo",
            ValidationType::ClinicalData => "clinical",
            ValidationType::BenchmarkModel => "benchmark",
        }
    }
}

/// Agreement level between model and reference data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AgreementLevel {
    /// Quantitative acceptance criterion met (with documented margin).
    Quantitative,
    /// Qualitative trend agreement only.
    Qualitative,
    /// Agreement not demonstrated.
    NotDemonstrated,
}

/// One completed verification activity with its evidence.
#[derive(Debug, Clone)]
pub struct VerificationActivity {
    /// Activity type.
    pub kind: VerificationType,
    /// Description of what was done.
    pub description: String,
    /// Result summary (metric values, references to test ids).
    pub results: String,
    /// Numeric metrics with acceptance criteria (checked mechanically by
    /// [`VerificationActivity::metrics_adequate`]).
    pub metrics: Vec<EvidenceMetric>,
}

/// A numeric evidence metric with a mechanical acceptance criterion, so
/// adequacy is checked by comparison rather than by reading `results` as
/// prose.
#[derive(Debug, Clone)]
pub struct EvidenceMetric {
    /// Metric name (e.g. `"max_abs_error_mpa"`).
    pub name: String,
    /// Computed value.
    pub value: f64,
    /// Acceptance band: `value` must satisfy every bound present.
    pub acceptance: Acceptance,
}

/// Acceptance bounds for an [`EvidenceMetric`].
#[derive(Debug, Clone, Copy, Default)]
pub struct Acceptance {
    /// `value <= max` (when present).
    pub max: Option<f64>,
    /// `value >= min` (when present).
    pub min: Option<f64>,
}

impl Acceptance {
    /// A band with only an upper bound.
    pub fn at_most(max: f64) -> Self {
        Self {
            max: Some(max),
            min: None,
        }
    }

    /// True when every present bound is satisfied.
    pub fn is_met_by(&self, value: f64) -> bool {
        let max_ok = self.max.is_none_or(|m| value <= m);
        let min_ok = self.min.is_none_or(|m| value >= m);
        max_ok && min_ok
    }
}

impl VerificationActivity {
    /// True when every attached metric meets its acceptance band.
    /// Activities with no metrics are trivially adequate (their evidence
    /// lives in `results` prose).
    pub fn metrics_adequate(&self) -> bool {
        self.metrics.iter().all(|m| m.acceptance.is_met_by(m.value))
    }

    /// Attaches a numeric metric.
    pub fn with_metric(&mut self, metric: EvidenceMetric) -> &mut Self {
        self.metrics.push(metric);
        self
    }
}

/// One completed validation activity with its evidence.
#[derive(Debug, Clone)]
pub struct ValidationActivity {
    /// Activity type.
    pub kind: ValidationType,
    /// Source of the reference data.
    pub reference: String,
    /// Comparison metric(s) used.
    pub metrics: Vec<String>,
    /// Agreement reached.
    pub agreement: AgreementLevel,
}

/// A credibility assessment for one question of interest.
#[derive(Debug, Clone)]
pub struct CredibilityAssessment {
    /// The question the model answers, verbatim.
    pub question_of_interest: String,
    /// Model risk.
    pub risk: ModelRisk,
    /// Model influence.
    pub influence: ModelInfluence,
    /// Completed verification activities.
    pub verification: Vec<VerificationActivity>,
    /// Completed validation activities.
    pub validation: Vec<ValidationActivity>,
}

/// Credibility goals derived from the V&V 40 risk-influence matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CredibilityGoals {
    /// Minimum number of distinct verification activity types required.
    pub min_verification_types: usize,
    /// Minimum number of distinct validation activity types required.
    pub min_validation_types: usize,
    /// Whether quantitative agreement is required (vs qualitative).
    pub quantitative_validation_required: bool,
    /// Whether an independent code-review record is required.
    pub independent_code_review_required: bool,
}

impl CredibilityAssessment {
    /// Risk × influence → credibility goals, per the V&V 40 matrix
    /// (monotone: higher risk/influence ⇒ more evidence). The full matrix
    /// interpolates the standard's Low/Medium/High bands.
    pub fn goals(&self) -> CredibilityGoals {
        use ModelInfluence::{Contributing as C, Direct as D, Significant as S};
        use ModelRisk::{High as RH, Low as RL, Medium as RM};
        match (self.risk, self.influence) {
            (RL, C) => CredibilityGoals {
                min_verification_types: 1,
                min_validation_types: 0,
                quantitative_validation_required: false,
                independent_code_review_required: false,
            },
            (RL, S) | (RM, C) => CredibilityGoals {
                min_verification_types: 1,
                min_validation_types: 1,
                quantitative_validation_required: false,
                independent_code_review_required: false,
            },
            (RL, D) | (RM, S) | (RH, C) => CredibilityGoals {
                min_verification_types: 2,
                min_validation_types: 1,
                quantitative_validation_required: false,
                independent_code_review_required: true,
            },
            (RM, D) | (RH, S) => CredibilityGoals {
                min_verification_types: 3,
                min_validation_types: 2,
                quantitative_validation_required: true,
                independent_code_review_required: true,
            },
            (RH, D) => CredibilityGoals {
                min_verification_types: 4,
                min_validation_types: 3,
                quantitative_validation_required: true,
                independent_code_review_required: true,
            },
        }
    }

    /// Evaluates the assessment: goals vs collected evidence. Returns the
    /// unmet goal descriptions (empty ⇒ credible per the matrix).
    pub fn evaluate(&self) -> Vec<String> {
        let goals = self.goals();
        let mut unmet = Vec::new();

        let mut v_kinds: Vec<VerificationType> = Vec::new();
        for a in &self.verification {
            if !v_kinds.contains(&a.kind) {
                v_kinds.push(a.kind);
            }
        }
        if v_kinds.len() < goals.min_verification_types {
            unmet.push(format!(
                "verification coverage {} < {} required",
                v_kinds.len(),
                goals.min_verification_types
            ));
        }
        if goals.independent_code_review_required
            && !v_kinds.contains(&VerificationType::CodeVerification)
        {
            unmet.push("independent code verification required".into());
        }

        let mut val_kinds: Vec<ValidationType> = Vec::new();
        let mut quantitative = false;
        for a in &self.validation {
            if !val_kinds.contains(&a.kind) {
                val_kinds.push(a.kind);
            }
            if a.agreement == AgreementLevel::Quantitative {
                quantitative = true;
            }
            if a.agreement == AgreementLevel::NotDemonstrated {
                unmet.push(format!(
                    "validation against {} not demonstrated",
                    a.reference
                ));
            }
        }
        if val_kinds.len() < goals.min_validation_types {
            unmet.push(format!(
                "validation coverage {} < {} required",
                val_kinds.len(),
                goals.min_validation_types
            ));
        }
        if goals.quantitative_validation_required && !quantitative {
            unmet.push("quantitative validation agreement required".into());
        }
        unmet
    }

    /// True when all goals are met.
    pub fn is_credible(&self) -> bool {
        self.evaluate().is_empty()
    }
}

impl CredibilityAssessment {
    /// Serialises the assessment (goals, activities, evaluation verdict) to
    /// JSON, so it can live inside a submission bundle next to the
    /// `tpt-med-fda` package. The output is deterministic (no HashMap
    /// iteration order) and contains no PHI by construction — activities
    /// are described in study terms, not patient terms.
    pub fn to_json(&self) -> String {
        fn esc(s: &str) -> String {
            let mut out = String::with_capacity(s.len() + 2);
            for c in s.chars() {
                match c {
                    '"' => out.push_str("\\\""),
                    '\\' => out.push_str("\\\\"),
                    '\n' => out.push_str("\\n"),
                    c => out.push(c),
                }
            }
            out
        }
        let risk = match self.risk {
            ModelRisk::Low => "low",
            ModelRisk::Medium => "medium",
            ModelRisk::High => "high",
        };
        let influence = match self.influence {
            ModelInfluence::Contributing => "contributing",
            ModelInfluence::Significant => "significant",
            ModelInfluence::Direct => "direct",
        };
        let goals = self.goals();
        let unmet = self.evaluate();
        let mut out = String::new();
        out.push_str("{\"question_of_interest\":\"");
        out.push_str(&esc(&self.question_of_interest));
        out.push_str("\",\"risk\":\"");
        out.push_str(risk);
        out.push_str("\",\"influence\":\"");
        out.push_str(influence);
        out.push_str("\",\"goals\":{");
        out.push_str(&format!(
            "\"min_verification_types\":{},\"min_validation_types\":{},\"quantitative_validation_required\":{},\"code_review_required\":{}",
            goals.min_verification_types,
            goals.min_validation_types,
            goals.quantitative_validation_required,
            goals.independent_code_review_required
        ));
        out.push_str("},\"verification\":[");
        for (i, a) in self.verification.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"type\":\"{}\",\"description\":\"{}\",\"results\":\"{}\",\"metrics_adequate\":{}}}",
                a.kind.key(),
                esc(&a.description),
                esc(&a.results),
                a.metrics_adequate()
            ));
        }
        out.push_str("],\"validation\":[");
        for (i, a) in self.validation.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let agreement = match a.agreement {
                AgreementLevel::Quantitative => "quantitative",
                AgreementLevel::Qualitative => "qualitative",
                AgreementLevel::NotDemonstrated => "not_demonstrated",
            };
            out.push_str(&format!(
                "{{\"type\":\"{}\",\"reference\":\"{}\",\"agreement\":\"{}\"}}",
                a.kind.key(),
                esc(&a.reference),
                agreement
            ));
        }
        out.push_str("],\"credible\":");
        out.push_str(if unmet.is_empty() { "true" } else { "false" });
        out.push_str(",\"unmet\":[");
        for (i, u) in unmet.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push('"');
            out.push_str(&esc(u));
            out.push('"');
        }
        out.push_str("]}");
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn activity(kind: VerificationType) -> VerificationActivity {
        VerificationActivity {
            kind,
            description: "performed".into(),
            results: "passed".into(),
            metrics: Vec::new(),
        }
    }

    fn validation(kind: ValidationType, agreement: AgreementLevel) -> ValidationActivity {
        ValidationActivity {
            kind,
            reference: "bench".into(),
            metrics: vec!["rmse".into()],
            agreement,
        }
    }

    #[test]
    fn evidence_metrics_are_checked_mechanically() {
        let mut a = activity(VerificationType::CalculationVerification);
        assert!(a.metrics_adequate(), "no metrics = trivially adequate");
        a.with_metric(EvidenceMetric {
            name: "max_abs_error_mpa".into(),
            value: 0.4,
            acceptance: Acceptance::at_most(0.5),
        });
        assert!(a.metrics_adequate());
        a.with_metric(EvidenceMetric {
            name: "cg_iterations".into(),
            value: 250.0,
            acceptance: Acceptance {
                min: Some(10.0),
                max: Some(500.0),
            },
        });
        assert!(a.metrics_adequate());
        a.metrics[0].value = 0.6;
        assert!(!a.metrics_adequate(), "breached upper bound must fail");
    }

    #[test]
    fn json_serialization_round_trips_the_verdict() {
        let a = CredibilityAssessment {
            question_of_interest: "peak stress below yield?".into(),
            risk: ModelRisk::Medium,
            influence: ModelInfluence::Significant,
            verification: vec![
                activity(VerificationType::CodeVerification),
                activity(VerificationType::CalculationVerification),
            ],
            validation: vec![validation(
                ValidationType::InVitro,
                AgreementLevel::Quantitative,
            )],
        };
        let json = a.to_json();
        assert!(json.contains("\"risk\":\"medium\""));
        assert!(json.contains("\"influence\":\"significant\""));
        assert!(json.contains("\"min_verification_types\":2"));
        assert!(json.contains("\"credible\":true"));
        assert!(json.contains("\"quantitative\""));
        // A failing assessment carries its unmet goals.
        let failing = CredibilityAssessment {
            question_of_interest: "q".into(),
            risk: ModelRisk::High,
            influence: ModelInfluence::Direct,
            verification: vec![],
            validation: vec![],
        };
        let json = failing.to_json();
        assert!(json.contains("\"credible\":false"));
        assert!(json.contains("\"unmet\":["));
    }

    #[test]
    fn matrix_is_monotone() {
        // Goals never decrease as risk increases (fixed influence) or as
        // influence increases (fixed risk) — checked axis-by-axis, since
        // row-major neighbours of the 3×3 matrix are not ordered pairs.
        use ModelInfluence::{Contributing, Direct, Significant};
        use ModelRisk::{High, Low, Medium};
        let goals = |r, i| {
            CredibilityAssessment {
                question_of_interest: "q".into(),
                risk: r,
                influence: i,
                verification: vec![],
                validation: vec![],
            }
            .goals()
        };
        for influence in [Contributing, Significant, Direct] {
            for w in [(Low, influence), (Medium, influence), (High, influence)].windows(2) {
                let g0 = goals(w[0].0, w[0].1);
                let g1 = goals(w[1].0, w[1].1);
                assert!(
                    g1.min_verification_types >= g0.min_verification_types,
                    "{g0:?} -> {g1:?}"
                );
                assert!(g1.min_validation_types >= g0.min_validation_types);
            }
        }
        for risk in [Low, Medium, High] {
            for w in [(risk, Contributing), (risk, Significant), (risk, Direct)].windows(2) {
                let g0 = goals(w[0].0, w[0].1);
                let g1 = goals(w[1].0, w[1].1);
                assert!(
                    g1.min_verification_types >= g0.min_verification_types,
                    "{g0:?} -> {g1:?}"
                );
                assert!(g1.min_validation_types >= g0.min_validation_types);
            }
        }
    }

    #[test]
    fn high_risk_direct_demands_full_evidence() {
        let a = CredibilityAssessment {
            question_of_interest: "stent radial force drives device selection".into(),
            risk: ModelRisk::High,
            influence: ModelInfluence::Direct,
            verification: vec![
                activity(VerificationType::CodeVerification),
                activity(VerificationType::CalculationVerification),
                activity(VerificationType::SensitivityAnalysis),
                activity(VerificationType::UncertaintyQuantification),
            ],
            validation: vec![
                validation(ValidationType::InVitro, AgreementLevel::Quantitative),
                validation(ValidationType::ClinicalData, AgreementLevel::Quantitative),
                validation(ValidationType::BenchmarkModel, AgreementLevel::Qualitative),
            ],
        };
        assert!(a.is_credible(), "{:?}", a.evaluate());
        let goals = a.goals();
        assert_eq!(goals.min_verification_types, 4);
        assert_eq!(goals.min_validation_types, 3);
        assert!(goals.quantitative_validation_required);
    }

    #[test]
    fn gaps_are_reported() {
        let a = CredibilityAssessment {
            question_of_interest: "q".into(),
            risk: ModelRisk::High,
            influence: ModelInfluence::Direct,
            verification: vec![activity(VerificationType::CodeVerification)],
            validation: vec![validation(
                ValidationType::InVitro,
                AgreementLevel::NotDemonstrated,
            )],
        };
        let unmet = a.evaluate();
        assert!(unmet.len() >= 3, "{unmet:?}");
        assert!(!a.is_credible());
    }

    #[test]
    fn low_risk_low_influence_needs_minimal_evidence() {
        let a = CredibilityAssessment {
            question_of_interest: "explorative what-if".into(),
            risk: ModelRisk::Low,
            influence: ModelInfluence::Contributing,
            verification: vec![activity(VerificationType::CodeVerification)],
            validation: vec![],
        };
        assert!(a.is_credible());
    }
}
