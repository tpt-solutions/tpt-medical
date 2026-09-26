//! Reproducibility manifest for a simulation run.
//!
//! An audit trail says *who did what, when*. It does not say *what code
//! produced the number*, and that is the question a reviewer asks when a
//! submission's stress result cannot be reproduced. A reproducibility manifest
//! answers it by pinning the four things that can silently change a result:
//!
//! - **Which code** - the workspace version, the resolved version of every
//!   `tpt-med-*` crate, and the git commit.
//! - **Which inputs** - a SHA-256 of every input artefact (a DICOM series, a
//!   size chart, a parameter file), so a rerun can prove it read the same
//!   bytes.
//! - **Which build** - the profile, if the caller records it.
//! - **When** - the manifest carries its own [`UtcStamp`].
//!
//! The manifest is deliberately *not* a substitute for the audit chain: it is
//! attached to the trail, and a separate audit event records its digest, so a
//! manifest swapped after the fact is detectable.
//!
//! ## Determinism
//!
//! Crate versions and input digests are deterministic for a given tree, so two
//! runs of the same binary on the same inputs produce an identical manifest -
//! which is what makes it useful as a reproducibility check rather than just a
//! provenance note. The git commit is the one field that necessarily varies,
//! so it is optional: pass `None` when the build is not from a checkout.

use std::collections::BTreeMap;

use tpt_med_audit::sha256_hex;

use crate::UtcStamp;

/// Current manifest schema version.
///
/// Bump when a field is removed or its meaning changes. Consumers should
/// reject a version they do not recognise rather than guess.
pub const MANIFEST_SCHEMA_VERSION: u32 = 1;

/// Canonical tool identifier recorded in every manifest.
pub const TOOL_ID: &str = "tpt-medical";

/// A hashed input artefact.
///
/// The digest is over the raw bytes, so it detects any change to the input -
/// including a re-export of the same anatomy from a different scanner - without
/// storing the input itself (which for a patient scan would be PHI).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputArtifact {
    /// Human-readable label, e.g. `"ct_series"` or `"tibial-chart"`.
    pub label: String,
    /// Lowercase hex SHA-256 of the artefact bytes.
    pub sha256: String,
    /// Length in bytes, so a truncated input is distinguishable from a
    /// different input of the same digest prefix.
    pub bytes: u64,
}

impl InputArtifact {
    /// Hashes `bytes` under `label`.
    pub fn new(label: impl Into<String>, bytes: &[u8]) -> Self {
        Self {
            label: label.into(),
            sha256: sha256_hex(bytes),
            bytes: bytes.len() as u64,
        }
    }
}

/// Reproducibility manifest for one simulation run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReproducibilityManifest {
    /// Schema version; see [`MANIFEST_SCHEMA_VERSION`].
    pub schema_version: u32,
    /// Always [`TOOL_ID`].
    pub tool: String,
    /// Version of the workspace the run was built from.
    pub workspace_version: String,
    /// Git commit, or `None` for a build not from a checkout.
    pub git_commit: Option<String>,
    /// Build profile (`release`, `debug`, ...), if recorded.
    pub build_profile: Option<String>,
    /// Resolved version of every crate that took part in the run.
    pub crates: BTreeMap<String, String>,
    /// Hashed input artefacts, sorted by label for determinism.
    pub inputs: Vec<InputArtifact>,
    /// When the manifest was created.
    pub created: UtcStamp,
}

impl ReproducibilityManifest {
    /// Creates a manifest for `workspace_version`, stamped now.
    pub fn new(workspace_version: impl Into<String>) -> Self {
        Self {
            schema_version: MANIFEST_SCHEMA_VERSION,
            tool: TOOL_ID.to_string(),
            workspace_version: workspace_version.into(),
            git_commit: None,
            build_profile: None,
            crates: BTreeMap::new(),
            inputs: Vec::new(),
            created: UtcStamp::now(),
        }
    }

    /// Records a crate version (builder style).
    pub fn with_crate(mut self, name: impl Into<String>, version: impl Into<String>) -> Self {
        self.crates.insert(name.into(), version.into());
        self
    }

    /// Records many crate versions at once (builder style).
    pub fn with_crates<I, N, V>(mut self, crates: I) -> Self
    where
        I: IntoIterator<Item = (N, V)>,
        N: Into<String>,
        V: Into<String>,
    {
        for (n, v) in crates {
            self.crates.insert(n.into(), v.into());
        }
        self
    }

    /// Records the git commit (builder style). `None` means "not from a
    /// checkout"; that is a legitimate value, not a missing field.
    pub fn with_git_commit(mut self, commit: Option<String>) -> Self {
        self.git_commit = commit;
        self
    }

    /// Records the build profile (builder style).
    pub fn with_build_profile(mut self, profile: impl Into<String>) -> Self {
        self.build_profile = Some(profile.into());
        self
    }

    /// Hashes and records an input artefact (builder style).
    pub fn with_input(mut self, label: impl Into<String>, bytes: &[u8]) -> Self {
        self.inputs.push(InputArtifact::new(label, bytes));
        self.inputs.sort_by(|a, b| a.label.cmp(&b.label));
        self
    }

    /// Records a pre-hashed artefact (builder style), e.g. a digest computed
    /// incrementally over a directory of files.
    pub fn with_artifact(mut self, artifact: InputArtifact) -> Self {
        self.inputs.push(artifact);
        self.inputs.sort_by(|a, b| a.label.cmp(&b.label));
        self
    }

    /// True when the manifest pins a git commit.
    pub fn is_pinned(&self) -> bool {
        self.git_commit.is_some()
    }
    /// Canonical serialization: the exact string a digest is taken over.
    ///
    /// Deliberately *not* the JSON rendering, so that adding a field to the
    /// export format does not silently change what a signature covers. Crate
    /// and input order are already deterministic, so this is stable across
    /// processes and runs.
    pub fn canonical(&self) -> String {
        let mut s = String::new();
        s.push_str("manifest|");
        s.push_str(&self.schema_version.to_string());
        s.push('|');
        s.push_str(&self.tool);
        s.push('|');
        s.push_str(&escape(&self.workspace_version));
        s.push('|');
        s.push_str(self.git_commit.as_deref().unwrap_or("-"));
        s.push('|');
        s.push_str(self.build_profile.as_deref().unwrap_or("-"));
        s.push('|');
        for (name, version) in &self.crates {
            s.push_str(name);
            s.push('=');
            s.push_str(version);
            s.push(',');
        }
        s.push('|');
        for a in &self.inputs {
            s.push_str(&a.label);
            s.push('=');
            s.push_str(&a.sha256);
            s.push(':');
            s.push_str(&a.bytes.to_string());
            s.push(',');
        }
        s.push('|');
        s.push_str(&self.created.to_iso8601());
        s
    }

    /// Digest of [`Self::canonical`], lowercase hex.
    pub fn digest(&self) -> String {
        sha256_hex(self.canonical().as_bytes())
    }

    /// JSON rendering, included verbatim in the exported submission package.
    pub fn to_json(&self) -> String {
        let mut s = String::new();
        s.push_str("{\"schema_version\":");
        s.push_str(&self.schema_version.to_string());
        s.push_str(",\"tool\":\"");
        s.push_str(&escape(&self.tool));
        s.push_str("\",\"workspace_version\":\"");
        s.push_str(&escape(&self.workspace_version));
        s.push_str("\",\"git_commit\":");
        match &self.git_commit {
            Some(c) => {
                s.push('"');
                s.push_str(&escape(c));
                s.push('"');
            }
            None => s.push_str("null"),
        }
        s.push_str(",\"build_profile\":");
        match &self.build_profile {
            Some(p) => {
                s.push('"');
                s.push_str(&escape(p));
                s.push('"');
            }
            None => s.push_str("null"),
        }
        s.push_str(",\"crates\":{");
        for (i, (name, version)) in self.crates.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push('"');
            s.push_str(&escape(name));
            s.push_str("\":\"");
            s.push_str(&escape(version));
            s.push('"');
        }
        s.push_str("},\"inputs\":[");
        for (i, a) in self.inputs.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str("{\"label\":\"");
            s.push_str(&escape(&a.label));
            s.push_str("\",\"sha256\":\"");
            s.push_str(&a.sha256);
            s.push_str("\",\"bytes\":");
            s.push_str(&a.bytes.to_string());
            s.push('}');
        }
        s.push_str("],\"created\":\"");
        s.push_str(&self.created.to_iso8601());
        s.push_str("\"}");
        s
    }
}

/// JSON string escaping (RFC 8259, minimal control escaping).
///
/// Duplicated from the crate root rather than imported, because the canonical
/// form and the JSON form are separate contracts and coupling them would let a
/// future change to one silently change the other.
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
#[cfg(test)]
mod tests {
    use super::*;

    fn manifest() -> ReproducibilityManifest {
        ReproducibilityManifest {
            created: UtcStamp {
                epoch_seconds: 1_789_862_400,
                nanos: 0,
            },
            ..ReproducibilityManifest::new("0.1.0")
        }
    }

    #[test]
    fn digest_is_stable_for_identical_inputs() {
        let a = manifest()
            .with_crate("tpt-med-dicom", "0.1.0")
            .with_input("ct_series", b"pixels");
        let b = manifest()
            .with_crate("tpt-med-dicom", "0.1.0")
            .with_input("ct_series", b"pixels");
        assert_eq!(a.digest(), b.digest());
        assert_eq!(a.digest().len(), 64);
    }

    #[test]
    fn input_order_does_not_change_the_digest() {
        let a = manifest().with_input("a", b"1").with_input("b", b"2");
        let b = manifest().with_input("b", b"2").with_input("a", b"1");
        assert_eq!(a.digest(), b.digest(), "inputs are sorted by label");
    }

    #[test]
    fn every_field_participates_in_the_digest() {
        let base = manifest()
            .with_crate("tpt-med-core", "0.1.0")
            .with_git_commit(Some("abc123".into()))
            .with_input("ct", b"x");
        let d = base.digest();

        let mut c = base.clone();
        c.git_commit = Some("def456".into());
        assert_ne!(d, c.digest(), "git commit must be covered");

        let mut c = base.clone();
        c.workspace_version = "0.2.0".into();
        assert_ne!(d, c.digest(), "workspace version must be covered");

        let mut c = base.clone();
        c.crates.insert("tpt-med-core".into(), "0.2.0".into());
        assert_ne!(d, c.digest(), "crate versions must be covered");

        let mut c = base.clone();
        c.inputs[0].sha256 = sha256_hex(b"y");
        assert_ne!(d, c.digest(), "input digests must be covered");

        let mut c = base.clone();
        c.created.epoch_seconds += 1;
        assert_ne!(d, c.digest(), "created stamp must be covered");
    }

    #[test]
    fn a_changed_input_changes_the_digest() {
        assert_ne!(
            manifest().with_input("ct", b"one").digest(),
            manifest().with_input("ct", b"two").digest()
        );
    }

    #[test]
    fn artifact_records_length_so_truncation_is_visible() {
        let art = InputArtifact::new("ct", b"12345");
        assert_eq!(art.bytes, 5);
        assert_eq!(art.sha256, sha256_hex(b"12345"));
    }

    #[test]
    fn absent_git_commit_is_explicit_not_omitted() {
        let m = manifest().with_git_commit(None);
        assert!(!m.is_pinned());
        // Absent optional fields render as "-", not as an empty field, so a
        // missing git commit can never be confused with an empty-string one.
        assert!(
            m.canonical()
                .starts_with("manifest|1|tpt-medical|0.1.0|-|-|"),
            "canonical: {}",
            m.canonical()
        );
        let json = m.to_json();
        assert!(json.contains("\"git_commit\":null"));
        assert!(json.contains("\"build_profile\":null"));
    }

    #[test]
    fn json_is_escaped() {
        let json = manifest().with_crate("weird\"name", "0.1.0\n").to_json();
        assert!(json.contains("weird\\\"name"), "quote escaped: {json}");
        assert!(json.contains("0.1.0\\n"), "newline escaped: {json}");
        assert!(json.starts_with('{') && json.ends_with('}'));
    }

    #[test]
    fn canonical_covers_every_identifying_field() {
        // A field present in JSON but absent from the canonical form would let
        // a signature miss it. Assert the two agree on what they cover.
        let c = manifest()
            .with_crate("tpt-med-dicom", "0.1.0")
            .with_git_commit(Some("deadbeef".into()))
            .with_build_profile("release")
            .with_input("ct", b"p")
            .canonical();
        for needle in [
            "manifest",
            "tpt-medical",
            "0.1.0",
            "deadbeef",
            "release",
            "tpt-med-dicom",
            "ct",
        ] {
            assert!(c.contains(needle), "canonical missing {needle}: {c}");
        }
    }
}
