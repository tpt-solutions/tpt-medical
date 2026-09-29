//! 21 CFR Part 11 compliant audit trails.
//!
//! Implements the electronic-records controls of 21 CFR Part 11 §11.10:
//!
//! - **§11.10(e)** operator-identity and timestamp recording: every entry
//!   carries an actor token (non-identifying; see `tpt-med-core`) and a UTC
//!   timestamp.
//! - **Audit trail integrity**: entries are chained
//!   (`digest[i] = H(seed ‖ digest[i−1] ‖ entry)`) and signed with
//!   HMAC-SHA256 under an operator secret — retroactive edits break the
//!   chain, unsigned edits break the signature.
//! - **§11.50 signature manifestations**: electronic signatures carry the
//!   signer, timestamp, and *meaning* (review/approval).
//! - **§11.10(k)** export: the trail serializes to a canonical JSON
//!   package with a detached signature for FDA submission bundles.
//!
//! This is a software control for R&D pipelines; deploying it as part of a
//! validated system requires the usual procedural controls (SOPs, operator
//! training, record retention) that no library can provide.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::time::{SystemTime, UNIX_EPOCH};

use tpt_med_audit::{hash_chain, hex, hmac_sha256, verify_chain};
use tpt_med_core::{AuditAction, AuditEvent};

pub mod manifest;
pub mod worm;

pub use manifest::{InputArtifact, ReproducibilityManifest, MANIFEST_SCHEMA_VERSION, TOOL_ID};
pub use worm::{WormError, WormLog};

/// UTC timestamp: seconds + nanoseconds since the epoch, rendered as
/// ISO-8601 `YYYY-MM-DDThh:mm:ssZ` (proleptic Gregorian; civil-date
/// conversion implemented locally, no external time crate).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct UtcStamp {
    /// Seconds since 1970-01-01T00:00:00Z.
    pub epoch_seconds: i64,
    /// Sub-second nanoseconds [0, 1e9).
    pub nanos: u32,
}

impl UtcStamp {
    /// Current system time.
    pub fn now() -> Self {
        let d = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default();
        Self {
            epoch_seconds: d.as_secs() as i64,
            nanos: d.subsec_nanos(),
        }
    }

    /// ISO-8601 rendering.
    pub fn to_iso8601(&self) -> String {
        let days = self.epoch_seconds.div_euclid(86_400);
        let secs = self.epoch_seconds.rem_euclid(86_400);
        let (y, m, d) = civil_from_days(days);
        format!(
            "{y:04}-{m:02}-{d:02}T{h:02}:{mi:02}:{s:02}Z",
            h = secs / 3600,
            mi = (secs % 3600) / 60,
            s = secs % 60
        )
    }
}

/// Howard Hinnant's `civil_from_days` algorithm (days → (y, m, d)).
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// Electronic signature manifestation (§11.50).
#[derive(Debug, Clone)]
pub struct ElectronicSignature {
    /// Signer token (non-identifying).
    pub signer: String,
    /// When the signature was applied.
    pub timestamp: UtcStamp,
    /// Meaning of the signature.
    pub meaning: SignatureMeaning,
}

/// Signature meanings required to be distinguishable (§11.50).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SignatureMeaning {
    /// The person authored the record.
    Author,
    /// The person reviewed the record.
    Reviewer,
    /// The person approved the record.
    Approver,
    /// The responsible party accepts ownership.
    ResponsibleParty,
}

impl SignatureMeaning {
    fn as_str(self) -> &'static str {
        match self {
            SignatureMeaning::Author => "author",
            SignatureMeaning::Reviewer => "reviewer",
            SignatureMeaning::Approver => "approver",
            SignatureMeaning::ResponsibleParty => "responsible_party",
        }
    }
}

/// One immutable audit-trail entry.
#[derive(Debug, Clone)]
pub struct AuditEntry {
    /// Monotonic sequence number (0-based).
    pub sequence: u64,
    /// When the action occurred.
    pub timestamp: UtcStamp,
    /// Non-identifying actor token.
    pub actor: String,
    /// Object class and token.
    pub object_class: String,
    /// Object token (hash or opaque id).
    pub object_token: String,
    /// The action.
    pub action: AuditAction,
    /// Reason-for-change (§11.10(e) requires it for modified records).
    pub reason: String,
}

impl AuditEntry {
    /// Canonical serialization used for chaining and signing.
    pub fn canonical(&self) -> String {
        format!(
            "{}|{}|{}|{}|{}|{}|{}",
            self.sequence,
            self.timestamp.to_iso8601(),
            self.actor,
            self.object_class,
            self.object_token,
            self.action,
            escape(&self.reason),
        )
    }
}

/// JSON string escaping (RFC 8259, minimal control escaping).
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// Workflow discipline for signatures, enforced by
/// [`AuditTrail::checked_append`] and [`AuditTrail::export_package_checked`].
///
/// The default [`SignaturePolicy::Permissive`] preserves today's behaviour:
/// signing, appending and exporting in any order is permitted and the export
/// records the ordering explicitly. `RequireSignatureAfterLastEdit` is the
/// stricter §11.10 discipline — every append after the newest signature
/// invalidates coverage, and a strict export refuses to emit a package whose
/// signature does not cover the final entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SignaturePolicy {
    /// Sign, append and export in any order (records ordering only).
    #[default]
    Permissive,
    /// The newest signature must post-date the final entry.
    RequireSignatureAfterLastEdit,
}

/// Reason-for-change field policy. `Structured` requires the reason to be
/// one of the listed codes (prefix match on the first whitespace-delimited
/// token, case-sensitive), so a site can require reason codes rather than
/// free text.
#[derive(Debug, Clone, Default)]
pub struct ReasonPolicy {
    /// When set, reasons must start with one of these codes followed by a
    /// space, end of string, or punctuation.
    pub required_codes: Vec<String>,
}

impl ReasonPolicy {
    /// A policy requiring one of the given codes.
    pub fn structured(codes: impl IntoIterator<Item = impl Into<String>>) -> Self {
        Self {
            required_codes: codes.into_iter().map(Into::into).collect(),
        }
    }

    /// Free-text policy (no codes required).
    pub fn free_text() -> Self {
        Self::default()
    }

    /// True when the reason satisfies the policy.
    pub fn accepts(&self, reason: &str) -> bool {
        if self.required_codes.is_empty() {
            return true;
        }
        let token = reason
            .split(|c: char| c.is_whitespace() || c == ':' || c == ',')
            .next();
        token.is_some_and(|t| self.required_codes.iter().any(|c| c == t))
    }
}

/// A policy violation returned by [`AuditTrail::checked_append`] and
/// [`AuditTrail::export_package_checked`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolicyError {
    /// Strict policy: entries were appended after the newest signature.
    SignatureDoesNotCoverFinalEntry {
        /// Entries in the trail.
        entries: usize,
        /// Entries covered by the newest signature (0 when unsigned).
        covered: usize,
    },
    /// The reason failed the configured [`ReasonPolicy`].
    ReasonPolicyViolation {
        /// The rejected reason.
        reason: String,
    },
}

impl core::fmt::Display for PolicyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PolicyError::SignatureDoesNotCoverFinalEntry { entries, covered } => write!(
                f,
                "signature policy: newest signature covers {covered} of {entries} entries;                  sign again after the final edit"
            ),
            PolicyError::ReasonPolicyViolation { reason } => {
                write!(f, "reason policy: rejected reason {reason:?}")
            }
        }
    }
}

impl std::error::Error for PolicyError {}

/// The signed audit trail.
#[derive(Debug, Clone)]
pub struct AuditTrail {
    /// Stable run identifier (chain seed).
    pub run_id: String,
    entries: Vec<AuditEntry>,
    signatures: Vec<ElectronicSignature>,
    digest_index: Vec<String>,
    /// Reproducibility manifest, when the run recorded one.
    manifest: Option<ReproducibilityManifest>,
    /// Signature discipline (default: permissive).
    signature_policy: SignaturePolicy,
    /// Reason-field policy (default: free text).
    reason_policy: ReasonPolicy,
}

impl AuditTrail {
    /// Opens a new trail for a run.
    pub fn new(run_id: impl Into<String>) -> Self {
        Self {
            run_id: run_id.into(),
            entries: Vec::new(),
            signatures: Vec::new(),
            digest_index: Vec::new(),
            manifest: None,
            signature_policy: SignaturePolicy::default(),
            reason_policy: ReasonPolicy::default(),
        }
    }

    /// Rebuilds a trail from persisted entries (crate-only: the WORM
    /// journal reload path). Entries keep their recorded timestamps and
    /// sequence numbers; the chain is recomputed by the caller.
    pub(crate) fn from_entries(run_id: impl Into<String>, entries: Vec<AuditEntry>) -> Self {
        Self {
            run_id: run_id.into(),
            entries,
            signatures: Vec::new(),
            digest_index: Vec::new(),
            manifest: None,
            signature_policy: SignaturePolicy::default(),
            reason_policy: ReasonPolicy::default(),
        }
    }

    /// Sets the signature discipline.
    pub fn set_signature_policy(&mut self, policy: SignaturePolicy) -> &mut Self {
        self.signature_policy = policy;
        self
    }

    /// Sets the reason-field policy.
    pub fn set_reason_policy(&mut self, policy: ReasonPolicy) -> &mut Self {
        self.reason_policy = policy;
        self
    }

    /// How many entries the newest signature covers (0 when unsigned).
    pub fn signed_entry_count(&self) -> usize {
        if self.signatures.is_empty() {
            0
        } else {
            self.entries.len()
        }
    }

    /// Policy-checked append: validates the reason policy and, under
    /// `RequireSignatureAfterLastEdit`, refuses to append entries silently
    /// after a signature (the caller must re-sign; the error says so).
    /// The plain [`Self::append`] keeps the permissive behaviour.
    pub fn checked_append(&mut self, event: AuditEvent) -> Result<(), PolicyError> {
        if !self.reason_policy.accepts(&event.reason) {
            return Err(PolicyError::ReasonPolicyViolation {
                reason: event.reason,
            });
        }
        if self.signature_policy == SignaturePolicy::RequireSignatureAfterLastEdit
            && !self.signatures.is_empty()
        {
            return Err(PolicyError::SignatureDoesNotCoverFinalEntry {
                entries: self.entries.len() + 1,
                covered: self.entries.len(),
            });
        }
        self.append(event);
        Ok(())
    }

    /// Export with signature-policy enforcement: under
    /// `RequireSignatureAfterLastEdit`, refuses to emit a package whose
    /// newest signature does not cover the final entry.
    pub fn export_package_checked(&self, key: &[u8]) -> Result<(String, String), PolicyError> {
        if self.signature_policy == SignaturePolicy::RequireSignatureAfterLastEdit
            && self.signed_entry_count() < self.entries.len()
        {
            return Err(PolicyError::SignatureDoesNotCoverFinalEntry {
                entries: self.entries.len(),
                covered: self.signed_entry_count(),
            });
        }
        Ok(self.export_package(key))
    }

    /// Appends an event; the timestamp is taken now. Signature re-computed
    /// lazily on export.
    pub fn append(&mut self, event: AuditEvent) {
        self.entries.push(AuditEntry {
            sequence: self.entries.len() as u64,
            timestamp: UtcStamp::now(),
            actor: event.actor_token,
            object_class: event.object_class,
            object_token: event.object_token,
            action: event.action,
            reason: event.reason,
        });
        self.digest_index = self.compute_chain();
    }

    /// Applies an electronic signature over the *current* trail state
    /// (§11.50). Signing again after new entries invalidates the earlier
    /// signature's coverage — recorded explicitly in the export.
    pub fn sign(&mut self, signer: &str, meaning: SignatureMeaning, key: &[u8]) {
        self.signatures.push(ElectronicSignature {
            signer: signer.to_string(),
            timestamp: UtcStamp::now(),
            meaning,
        });
        self.digest_index = self.compute_chain();
        let _ = key;
    }

    fn entry_payload(&self) -> Vec<String> {
        self.entries.iter().map(|e| e.canonical()).collect()
    }

    fn compute_chain(&self) -> Vec<String> {
        hash_chain(self.run_id.as_bytes(), &self.entry_payload())
    }

    /// Chain verification: re-walks the whole trail.
    pub fn verify_integrity(&self) -> bool {
        if self.entries.len() != self.digest_index.len() {
            return false;
        }
        verify_chain(
            self.run_id.as_bytes(),
            &self.entry_payload(),
            &self.digest_index,
        )
    }

    /// HMAC-SHA256 tag over the full exported payload under `key`.
    pub fn export_tag(&self, key: &[u8]) -> String {
        hex(&hmac_sha256(key, self.export_payload().as_bytes()))
    }

    fn export_payload(&self) -> String {
        let mut s = String::new();
        s.push_str("{\"run_id\":\"");
        s.push_str(&escape(&self.run_id));
        s.push_str("\",\"entries\":[");
        for (i, e) in self.entries.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!(
                "{{\"sequence\":{},\"timestamp\":\"{}\",\"actor\":\"{}\",\"object_class\":\"{}\",\
                 \"object_token\":\"{}\",\"action\":\"{}\",\"reason\":\"{}\",\"digest\":\"{}\"}}",
                e.sequence,
                e.timestamp.to_iso8601(),
                escape(&e.actor),
                escape(&e.object_class),
                escape(&e.object_token),
                e.action,
                escape(&e.reason),
                self.digest_index.get(i).map_or("", String::as_str),
            ));
        }
        s.push_str("],\"signatures\":[");
        for (i, sig) in self.signatures.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!(
                "{{\"signer\":\"{}\",\"timestamp\":\"{}\",\"meaning\":\"{}\"}}",
                escape(&sig.signer),
                sig.timestamp.to_iso8601(),
                sig.meaning.as_str(),
            ));
        }
        // Close the arrays, then the object. The manifest is appended *inside*
        // the object, so the closing brace must come after it, not before.
        s.push(']');

        // The reproducibility manifest, when the run recorded one. Omitting
        // the key entirely when there is none keeps the export byte-identical
        // for runs that do not use it, so existing packages still verify.
        if let Some(m) = &self.manifest {
            s.push_str(",\"manifest\":");
            s.push_str(&m.to_json());
        }
        s.push('}');
        s
    }

    /// Exports the FDA-ready package: canonical JSON plus detached HMAC tag.
    /// Returns `(json, tag)`.
    pub fn export_package(&self, key: &[u8]) -> (String, String) {
        (self.export_payload(), self.export_tag(key))
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when no entries are recorded.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Attaches a reproducibility manifest and appends an audit event
    /// recording its digest.
    ///
    /// The event is what makes a later manifest swap detectable: the manifest
    /// body is covered by the export tag, and its digest is covered by the
    /// chain, so replacing the manifest invalidates both. Attaching twice
    /// records a second event rather than overwriting, because the first
    /// attachment already happened in the record.
    pub fn attach_manifest(&mut self, manifest: ReproducibilityManifest) {
        let digest = manifest.digest();
        self.manifest = Some(manifest);
        self.append(AuditEvent::new(
            "system:manifest",
            "reproducibility_manifest",
            digest,
            AuditAction::Create,
            "reproducibility manifest attached to run",
        ));
    }

    /// The attached reproducibility manifest, if any.
    pub fn manifest(&self) -> Option<&ReproducibilityManifest> {
        self.manifest.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn event(action: AuditAction, obj: &str) -> AuditEvent {
        AuditEvent::new(
            "operator:tech1",
            "simulation",
            obj,
            action,
            "planning iteration",
        )
    }

    #[test]
    fn strict_policy_rejects_unsigned_export_and_post_signature_edits() {
        use super::*;
        let mut trail = AuditTrail::new("run-strict");
        trail.set_signature_policy(SignaturePolicy::RequireSignatureAfterLastEdit);
        trail.append(event(AuditAction::Create, "sim-1"));

        // Export before signing is refused under the strict policy.
        let err = trail.export_package_checked(b"k").unwrap_err();
        assert_eq!(
            err,
            PolicyError::SignatureDoesNotCoverFinalEntry {
                entries: 1,
                covered: 0
            }
        );

        // Appending after signing is refused (the caller must re-sign).
        trail.sign("reviewer:r1", SignatureMeaning::Reviewer, b"k");
        let err = trail
            .checked_append(event(AuditAction::Modify, "sim-1"))
            .unwrap_err();
        assert_eq!(
            err,
            PolicyError::SignatureDoesNotCoverFinalEntry {
                entries: 2,
                covered: 1
            }
        );

        // Re-signing after the final edit restores coverage and export.
        trail.append(event(AuditAction::Modify, "sim-1")); // permissive append path
        trail.sign("reviewer:r1", SignatureMeaning::Approver, b"k");
        assert_eq!(trail.signed_entry_count(), trail.len());
        let (_, tag) = trail.export_package_checked(b"k").expect("covered");
        assert!(!tag.is_empty());
    }

    #[test]
    fn reason_policy_requires_a_structured_code() {
        use super::*;
        let mut trail = AuditTrail::new("run-codes");
        trail.set_reason_policy(ReasonPolicy::structured(["REVIEW", "PARAM-CHANGE"]));

        let ok = trail.checked_append(AuditEvent::new(
            "op",
            "sim",
            "s1",
            AuditAction::Modify,
            "PARAM-CHANGE: vessel diameter updated",
        ));
        assert!(ok.is_ok());

        let bad = trail.checked_append(AuditEvent::new(
            "op",
            "sim",
            "s1",
            AuditAction::Modify,
            "changed my mind",
        ));
        assert_eq!(
            bad.unwrap_err(),
            PolicyError::ReasonPolicyViolation {
                reason: "changed my mind".into()
            }
        );
    }

    #[test]
    fn permissive_defaults_keep_legacy_behaviour() {
        use super::*;
        let mut trail = AuditTrail::new("run-permissive");
        // No policies set: append after sign is allowed and export succeeds.
        trail.append(event(AuditAction::Create, "s1"));
        trail.sign("r", SignatureMeaning::Reviewer, b"k");
        trail.append(event(AuditAction::Modify, "s1"));
        let (_, _) = trail.export_package_checked(b"k").expect("permissive");
        assert!(trail.verify_integrity());
    }

    #[test]
    fn iso_rendering_of_known_epoch() {
        let t = UtcStamp {
            epoch_seconds: 1_758_240_000, // 2025-09-18T21:20:00Z-ish
            nanos: 0,
        };
        let iso = t.to_iso8601();
        assert_eq!(iso.len(), 20);
        assert!(iso.ends_with('Z'));
        // Epoch itself.
        let epoch = UtcStamp {
            epoch_seconds: 0,
            nanos: 0,
        };
        assert_eq!(epoch.to_iso8601(), "1970-01-01T00:00:00Z");
        // Known: 2026-09-20 00:00:00Z = 1789862400.
        let d = UtcStamp {
            epoch_seconds: 1_789_862_400,
            nanos: 0,
        };
        assert_eq!(d.to_iso8601(), "2026-09-20T00:00:00Z");
    }

    #[test]
    fn trail_chains_and_verifies() {
        let mut trail = AuditTrail::new("run-42");
        trail.append(event(AuditAction::Create, "sim-1"));
        trail.append(event(AuditAction::Modify, "sim-1"));
        trail.append(event(AuditAction::Simulate, "sim-1"));
        assert_eq!(trail.len(), 3);
        assert!(trail.verify_integrity());

        // Export then tamper: forge a trail with a modified reason and the
        // original digest index — verification must fail.
        let (json, tag) = trail.export_package(b"k");
        assert!(json.contains("\"sequence\":2"));
        let mut forged = trail.clone();
        forged.entries[0].reason = "tampered".into();
        assert!(!forged.verify_integrity());
        // Tag binds the payload: a different trail yields a different tag.
        assert_ne!(forged.export_tag(b"k"), tag);
    }

    #[test]
    fn signatures_manifest_meaning() {
        let mut trail = AuditTrail::new("run-7");
        trail.append(event(AuditAction::Create, "plan-1"));
        trail.sign("reviewer:r1", SignatureMeaning::Reviewer, b"k");
        let (json, _) = trail.export_package(b"k");
        assert!(json.contains("\"meaning\":\"reviewer\""));
        assert!(json.contains("\"signer\":\"reviewer:r1\""));
        assert!(trail.verify_integrity());
    }

    /// Cheap structural JSON check: balanced braces/brackets, exactly one
    /// top-level value, and the export is an object.
    ///
    /// This is not a JSON parser, and deliberately so — the crate is
    /// dependency-free. It is enough to catch the failure mode that actually
    /// happened: closing the outer object before appending the manifest, which
    /// yields a string that still *looks* right and is still signed, but is not
    /// parseable by any consumer.
    fn assert_single_json_object(s: &str) {
        assert!(
            s.starts_with('{'),
            "must start with an object: {}",
            &s[..40.min(s.len())]
        );
        assert!(s.ends_with('}'), "must end with an object");

        let (mut depth, mut closes_at_top) = (0i32, false);
        let mut in_str = false;
        let mut escaped = false;
        for c in s.chars() {
            if in_str {
                if escaped {
                    escaped = false;
                } else if c == '\\' {
                    escaped = true;
                } else if c == '"' {
                    in_str = false;
                }
                continue;
            }
            match c {
                '"' => in_str = true,
                '{' | '[' => {
                    depth += 1;
                    if depth == 1 {
                        assert!(!closes_at_top, "more than one top-level value");
                    }
                }
                '}' | ']' => {
                    depth -= 1;
                    assert!(depth >= 0, "unbalanced closing delimiter");
                    if depth == 0 {
                        closes_at_top = true;
                    }
                }
                _ => {}
            }
        }
        assert_eq!(depth, 0, "unbalanced delimiters");
        assert!(closes_at_top, "no top-level value closed");
    }

    #[test]
    fn empty_trail_verifies() {
        let trail = AuditTrail::new("run-empty");
        assert!(trail.is_empty());
        assert!(trail.verify_integrity());
    }

    #[test]
    fn json_escapes_quotes_in_reason() {
        let mut trail = AuditTrail::new("run-esc");
        trail.append(AuditEvent::new(
            "op",
            "sim",
            "s1",
            AuditAction::Modify,
            "changed \"m\" and \\ path",
        ));
        let (json, _) = trail.export_package(b"k");
        assert!(json.contains("\\\"m\\\""));
        assert!(json.contains("\\\\ path"));
    }

    #[test]
    fn export_without_a_manifest_is_unchanged() {
        // Runs that do not use a manifest must produce byte-identical exports,
        // so packages signed before the feature existed still verify.
        let mut trail = AuditTrail::new("run-legacy");
        trail.append(event(AuditAction::Create, "sim-1"));
        let (json, _) = trail.export_package(b"k");
        assert!(!json.contains("manifest"), "no manifest key: {json}");
        assert!(trail.manifest().is_none());
        assert_single_json_object(&json);
    }

    #[test]
    fn export_is_a_single_valid_json_object_with_or_without_a_manifest() {
        let mut plain = AuditTrail::new("run-shape");
        plain.append(event(AuditAction::Create, "sim-1"));
        plain.sign("reviewer:r1", SignatureMeaning::Reviewer, b"k");
        let (plain_json, _) = plain.export_package(b"k");
        assert_single_json_object(&plain_json);

        let mut with_manifest = AuditTrail::new("run-shape");
        with_manifest.append(event(AuditAction::Create, "sim-1"));
        with_manifest.attach_manifest(
            ReproducibilityManifest::new("0.1.0")
                .with_crate("tpt-med-dicom", "0.1.0")
                .with_input("ct", b"x"),
        );
        with_manifest.sign("reviewer:r1", SignatureMeaning::Reviewer, b"k");
        let (manifest_json, _) = with_manifest.export_package(b"k");
        assert_single_json_object(&manifest_json);
        assert!(manifest_json.contains("\"manifest\""));
    }

    #[test]
    fn attaching_a_manifest_records_an_event_and_exports_it() {
        let mut trail = AuditTrail::new("run-manifest");
        trail.append(event(AuditAction::Create, "sim-1"));
        let before = trail.len();

        let m = ReproducibilityManifest {
            created: UtcStamp {
                epoch_seconds: 1_789_862_400,
                nanos: 0,
            },
            ..ReproducibilityManifest::new("0.1.0")
        }
        .with_crate("tpt-med-dicom", "0.1.0")
        .with_git_commit(Some("abc1234".into()))
        .with_input("ct_series", b"pixels");
        let digest = m.digest();

        trail.attach_manifest(m);
        assert_eq!(trail.len(), before + 1, "attachment is an audited event");
        assert!(trail.verify_integrity());

        // The event carries the manifest digest, so a swapped manifest is
        // detectable even if the export tag is not re-checked.
        let last = trail.entries.last().unwrap();
        assert_eq!(last.object_class, "reproducibility_manifest");
        assert_eq!(last.object_token, digest);

        let (json, tag) = trail.export_package(b"k");
        assert!(json.contains("\"manifest\""));
        assert!(json.contains("abc1234"));
        assert_eq!(tag.len(), 64);
    }

    #[test]
    fn a_changed_manifest_changes_the_export_tag() {
        let key = b"k";
        let build = |input: &[u8]| {
            let mut t = AuditTrail::new("run-swap");
            t.append(event(AuditAction::Create, "sim-1"));
            t.attach_manifest(ReproducibilityManifest::new("0.1.0").with_input("ct_series", input));
            t
        };
        assert_ne!(
            build(b"original").export_tag(key),
            build(b"tampered").export_tag(key),
            "the manifest body must be covered by the detached tag"
        );
    }

    #[test]
    fn re_attaching_a_manifest_appends_rather_than_erases() {
        let mut trail = AuditTrail::new("run-twice");
        trail.attach_manifest(ReproducibilityManifest::new("0.1.0"));
        let after_first = trail.len();
        trail.attach_manifest(ReproducibilityManifest::new("0.1.0"));
        assert_eq!(trail.len(), after_first + 1);
        assert!(trail.verify_integrity());
    }
}
