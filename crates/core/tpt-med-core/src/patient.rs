//! Patient model and demographics with HIPAA-safe redaction.

use crate::anatomy::AnatomicalModel;
use crate::audit::AuditSubject;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// Opaque, non-identifying patient token.
///
/// Constructed from local identifiers via a non-cryptographic hash; the raw
/// identifier is discarded and never stored, displayed, or logged. This is
/// a *de-identification* measure (it prevents accidental PHI leakage through
/// logs and error paths), not a security boundary: treat any mapping table
/// between real identities and tokens as PHI under your institution's
/// HIPAA/GDPR controls.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PatientId(u64);

impl PatientId {
    /// Derives an opaque token from a local identifier. The mapping is not
    /// reversible within this crate and the source string is not retained.
    pub fn from_hash(local_id: &str) -> Self {
        let mut hasher = DefaultHasher::new();
        local_id.hash(&mut hasher);
        Self(hasher.finish())
    }

    /// Token for a fully de-identified (anonymized) dataset.
    pub const ANONYMOUS: PatientId = PatientId(0);

    /// Hex rendering for exports (still non-identifying).
    pub fn to_hex(self) -> String {
        format!("{:016x}", self.0)
    }
}

impl core::fmt::Display for PatientId {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "patient:{}", self.to_hex())
    }
}

/// Recorded demographics. Values may be used for population-level scaling
/// (e.g. body-weight-based loads) but render only as coarse ranges in audit
/// output.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Demographics {
    /// Age in whole years.
    pub age_years: u8,
    /// Recorded sex.
    pub sex: Sex,
    /// Body mass in kg.
    pub weight_kg: f64,
    /// Stature in cm.
    pub height_cm: f64,
}

impl Demographics {
    /// Estimated body weight in newtons (local gravity 9.81 m/s²). Used for
    /// body-weight-multiple loading conventions.
    pub fn body_weight_force(&self) -> f64 {
        self.weight_kg * 9.81
    }
}

impl AuditSubject for Demographics {
    fn audit_token(&self) -> String {
        // Decade buckets only: sufficient for population scaling context,
        // insufficient to identify.
        let decade = (self.age_years / 10) * 10;
        let bmi = self.weight_kg / (self.height_cm / 100.0).powi(2);
        let bmi_band = match bmi {
            b if b < 18.5 => "underweight",
            b if b < 25.0 => "normal",
            b if b < 30.0 => "overweight",
            _ => "obese",
        };
        format!("demographics:age~{decade}s:bmi-band={bmi_band}")
    }
}

/// Recorded sex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    /// Male.
    Male,
    /// Female.
    Female,
    /// Other / not recorded.
    Other,
}

/// A patient-specific simulation subject: identity token, demographics, and
/// the assembled anatomical model.
#[derive(Debug, Clone)]
pub struct PatientModel {
    /// Opaque patient token.
    pub id: PatientId,
    /// Demographics (optional for fully anonymized pipelines).
    pub demographics: Option<Demographics>,
    /// Assembled anatomy.
    pub anatomy: AnatomicalModel,
}

impl PatientModel {
    /// Creates a patient model from parts.
    pub fn new(id: PatientId) -> Self {
        Self {
            id,
            demographics: None,
            anatomy: AnatomicalModel::default(),
        }
    }

    /// Attaches demographics (builder style).
    pub fn with_demographics(mut self, demographics: Demographics) -> Self {
        self.demographics = Some(demographics);
        self
    }

    /// Body-weight load scale in newtons, if demographics are present.
    pub fn body_weight_force(&self) -> Option<f64> {
        self.demographics.map(|d| d.body_weight_force())
    }
}

impl AuditSubject for PatientModel {
    fn audit_token(&self) -> String {
        self.id.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn id_is_stable_but_opaque() {
        let a = PatientId::from_hash("mrn-12345");
        let b = PatientId::from_hash("mrn-12345");
        assert_eq!(a, b);
        assert_ne!(a, PatientId::from_hash("mrn-54321"));
        assert!(!format!("{a}").contains("12345"));
    }

    #[test]
    fn body_weight_scaling() {
        let d = Demographics {
            age_years: 70,
            sex: Sex::Male,
            weight_kg: 80.0,
            height_cm: 175.0,
        };
        assert!((d.body_weight_force() - 784.8).abs() < 1e-9);
        let model = PatientModel::new(PatientId::ANONYMOUS).with_demographics(d);
        assert!((model.body_weight_force().unwrap() - 784.8).abs() < 1e-9);
    }
}
