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

#[cfg(test)]
mod tests {
    use super::*;

    fn activity(kind: VerificationType) -> VerificationActivity {
        VerificationActivity {
            kind,
            description: "performed".into(),
            results: "passed".into(),
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
