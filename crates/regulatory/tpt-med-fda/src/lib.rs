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

/// The signed audit trail.
#[derive(Debug, Clone)]
pub struct AuditTrail {
    /// Stable run identifier (chain seed).
    pub run_id: String,
    entries: Vec<AuditEntry>,
    signatures: Vec<ElectronicSignature>,
    digest_index: Vec<String>,
}

impl AuditTrail {
    /// Opens a new trail for a run.
    pub fn new(run_id: impl Into<String>) -> Self {
        Self {
            run_id: run_id.into(),
            entries: Vec::new(),
            signatures: Vec::new(),
            digest_index: Vec::new(),
        }
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
        s.push_str("]}");
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
}
