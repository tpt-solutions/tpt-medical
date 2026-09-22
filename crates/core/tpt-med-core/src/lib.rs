//! Core domain types for medical simulation: patient models, anatomical
//! taxonomy, and HIPAA-safe audit traits.
//!
//! # Privacy model
//!
//! Patient identifiers are opaque, non-printing types ([`PatientId`]) and
//! demographic values are never carried into log output: anything that
//! implements [`AuditSubject`] must redact itself to a stable,
//! non-identifying token. Real patient data must never enter this stack —
//! all shipped test data is synthetic.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod anatomy;
pub mod audit;
pub mod patient;

pub use anatomy::{
    AnatomicalModel, AnatomicalRegion, AortaSegment, BoneType, ImplantType, Landmark, OrganType,
    SoftTissueType, VesselType,
};
pub use audit::{AuditAction, AuditEvent, AuditSubject};
pub use patient::{Demographics, PatientId, PatientModel, Sex};

/// Crate-level result alias used across the tpt-medical stack.
pub type Result<T, E> = core::result::Result<T, E>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn patient_id_does_not_leak_identity() {
        let id = PatientId::from_hash("real-name-jane-doe");
        // Debug/Display render the opaque token, never the source identity.
        let rendered = format!("{id:?}");
        assert!(!rendered.contains("jane"), "{rendered}");
        assert!(!rendered.contains("doe"), "{rendered}");
    }

    #[test]
    fn demographics_redact_in_audit_output() {
        let d = Demographics {
            age_years: 67,
            sex: Sex::Female,
            weight_kg: 72.0,
            height_cm: 168.0,
        };
        let token = d.audit_token();
        assert!(!token.contains("67"));
        assert!(token.starts_with("demographics:"));
    }
}
