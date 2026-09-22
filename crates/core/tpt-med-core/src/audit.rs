//! HIPAA-safe audit traits.
//!
//! Anything entering an audit log implements [`AuditSubject`] and must
//! render a non-identifying token. [`AuditEvent`] captures a domain action
//! without payload values that could carry PHI.

/// A type that can render itself as a non-identifying audit token.
///
/// Contract: the token must not contain direct identifiers (names, MRNs,
/// dates of birth, free-text notes) and must be stable across processes so
/// audit trails remain joinable.
pub trait AuditSubject {
    /// Non-identifying, stable token for audit logs.
    fn audit_token(&self) -> String;
}

/// Domain actions recorded in audit trails (21 CFR Part 11 §11.10(e) wants
/// operation-level recording; see `tpt-med-fda` for the signed trail).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AuditAction {
    /// Record created.
    Create,
    /// Record modified.
    Modify,
    /// Record deleted (soft-delete in regulated pipelines).
    Delete,
    /// Record approved by an authorized reviewer.
    Approve,
    /// Record rejected by an authorized reviewer.
    Reject,
    /// Record exported outside the application boundary.
    Export,
    /// Simulation executed.
    Simulate,
}

impl core::fmt::Display for AuditAction {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let name = match self {
            AuditAction::Create => "create",
            AuditAction::Modify => "modify",
            AuditAction::Delete => "delete",
            AuditAction::Approve => "approve",
            AuditAction::Reject => "reject",
            AuditAction::Export => "export",
            AuditAction::Simulate => "simulate",
        };
        f.write_str(name)
    }
}

/// A minimal audit event at the domain level; the regulatory layer
/// (`tpt-med-fda`) timestamps, chains, and signs these.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuditEvent {
    /// Who performed the action, as a non-identifying token.
    pub actor_token: String,
    /// What object class was acted on (e.g. `"patient_model"`, `"mesh"`).
    pub object_class: String,
    /// Non-identifying object token (hash or opaque id).
    pub object_token: String,
    /// The action.
    pub action: AuditAction,
    /// Free-text reason for the action (user-supplied, may be reviewed).
    ///
    /// Callers in regulated deployments are responsible for ensuring the
    /// reason field policy complies with their SOPs (e.g. no PHI in reasons).
    pub reason: String,
}

impl AuditEvent {
    /// Convenience constructor.
    pub fn new(
        actor_token: impl Into<String>,
        object_class: impl Into<String>,
        object_token: impl Into<String>,
        action: AuditAction,
        reason: impl Into<String>,
    ) -> Self {
        Self {
            actor_token: actor_token.into(),
            object_class: object_class.into(),
            object_token: object_token.into(),
            action,
            reason: reason.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn action_display() {
        assert_eq!(AuditAction::Simulate.to_string(), "simulate");
        assert_eq!(AuditAction::Approve.to_string(), "approve");
    }

    #[test]
    fn event_fields() {
        let e = AuditEvent::new(
            "user:op-1",
            "mesh",
            "mesh:deadbeef",
            AuditAction::Export,
            "planning review",
        );
        assert_eq!(e.object_class, "mesh");
        assert_eq!(e.action, AuditAction::Export);
    }
}
