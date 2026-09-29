# tpt-med-fda

21 CFR Part 11 compliant audit trails — immutable append-only event log,
electronic signatures with meaning, and a signed export package.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--fda-orange)](https://crates.io/crates/tpt-med-fda)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--fda-blue)](https://docs.rs/tpt-med-fda)

| | |
|---|---|
| **Layer** | `regulatory` |
| **Status** | Alpha, `0.1.0` |
| **Scope** | RFC 0003 — FDA audit trail |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-audit`](../tpt-med-audit), [`tpt-med-core`](../../core/tpt-med-core) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

If you intend to put simulation results in an FDA submission, the record of
how they were produced has to satisfy 21 CFR Part 11 §11.10. Concretely, that
means: operator identity and time recorded on every operation, an audit trail
you cannot quietly rewrite, and electronic signatures that say what they mean.
Those are four specific, testable obligations, and this crate implements them
directly.

What it does **not** do is pretend to be a compliance certificate. Part 11
compliance is a property of a *system* — procedures, training, record
retention, access control — and no library can supply those. What a library
can do is make the software controls auditable and make the failure modes
detectable, and that is the scope here.

## Regulatory Mapping

| Control | Implementation |
|---|---|
| **§11.10(e)** — operator identity and time of each operation | `AuditEntry` carries a non-identifying `actor` token (from the [`tpt-med-core`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-core) privacy model) and a `UtcStamp` taken at `append` time |
| **Audit trail integrity** | Entries are chained, `digest[i] = H(run_id ‖ digest[i−1] ‖ canonical(entry))`; edits break the chain, and the export is bound by a detached HMAC-SHA256 tag |
| **§11.50** — signature manifestations | `ElectronicSignature` records signer, timestamp, **and meaning**: `Author`, `Reviewer`, `Approver`, `ResponsibleParty` — the meanings must be distinguishable, so they are an enum, not free text |
| **§11.10(k)** — export | `export_package` returns canonical JSON plus a detached tag, for a submission bundle |

## Features

- **`AuditTrail::new(run_id)`** — one trail per simulation run; the `run_id` is
  the chain seed, so two runs can never have interchangeable digests.
- **`append(AuditEvent)`** — append-only. There is no remove, no update, and
  no setter for the entry list. Modifying an existing record means appending a
  new `Modify` event with a reason, which is what §11.10(e) requires.
- **`sign(signer, meaning, key)`** — an electronic signature over the current
  trail state. Signing again after new entries explicitly narrows what the
  earlier signature covered, and the export records that.
- **`verify_integrity()`** — re-walks the whole chain; detects edits,
  truncation and seed splicing.
- **`export_tag(key)` / `export_package(key)`** — canonical JSON with a
  detached HMAC-SHA256 tag. The JSON is hand-serialised with RFC 8259
  escaping so the output is byte-stable and the tag is reproducible.
- **`AuditEntry::canonical()`** — the exact string the chain commits to. Being
  explicit about the canonical form is what makes a signature portable
  between implementations.
- **Local UTC implementation** — `UtcStamp` and a proleptic-Gregorian
  civil-date conversion (Hinnant's algorithm), so there is no time-crate
  dependency and no timezone database to get wrong.
- **Reproducibility manifest** — `ReproducibilityManifest` answers the question
  the audit trail does not: *what code and what inputs produced this number?*
  It pins the workspace version, the version of every participating crate, the
  git commit, the build profile, and a SHA-256 of every input artefact.
  `AuditTrail::attach_manifest` attaches it **and records an audit event
  carrying its digest**, so a manifest swapped after the fact is detectable
  from the chain alone.
- **Not a PHI store** — actor tokens are non-identifying by construction, and
  the manifest hashes its inputs rather than embedding them, so a patient scan
  never enters the record.

- **`WormLog`** — an append-only, write-once file journal for a trail:
  every entry is stored with its running hash-chain digest and `fsync`ed
  before the append returns; reopening verifies the header, every chain
  link, and sequence continuity, refusing torn tails (`TornTail`), retroactive
  edits and cross-run splices (`ChainBroken` at the first bad link).
  `into_trail` rebuilds a live, self-verifying `AuditTrail`.

## Conventions

- `sequence` is a **0-based monotonic counter** assigned at `append`.
- Timestamps are **UTC**, rendered ISO-8601 as `YYYY-MM-DDThh:mm:ssZ`.
- `reason` is free text and **is the caller's policy surface**: a regulated
  deployment must ensure it carries no PHI. This is documented rather than
  enforced, because the correct policy is site-specific.
- Signing is **not** mutually exclusive with appending; a trail may carry
  several signatures with different meanings.
- The export tag is computed **over the export payload**, so it binds the
  entries, their digests *and* the signature manifestations together.

## Usage

```rust
use tpt_med_audit::{hex, hmac_sha256};
use tpt_med_core::{AuditAction, AuditEvent, PatientId};
use tpt_med_fda::{AuditTrail, SignatureMeaning};

fn main() {
    let mut trail = AuditTrail::new("run-2026-09-26-femur-001");
    assert!(trail.is_empty());

    // Append-only. Every operation is an event with an actor and a timestamp.
    let patient = PatientId::ANONYMOUS.to_string();
    trail.append(AuditEvent::new(
        "user:op-1", "patient_model", patient,
        AuditAction::Create, "phantom study",
    ));
    trail.append(AuditEvent::new(
        "user:op-1", "mesh", "mesh:deadbeef",
        AuditAction::Simulate, "stance screening",
    ));
    assert_eq!(trail.len(), 2);

    // A modification is a NEW event with a reason, never an in-place edit.
    trail.append(AuditEvent::new(
        "user:op-2", "mesh", "mesh:deadbeef",
        AuditAction::Modify, "re-thresholded at 300 HU",
    ));

    // Signatures with distinguishable meanings (§11.50).
    let key = [42u8; 32];
    trail.sign("user:reviewer-7", SignatureMeaning::Reviewer, &key);
    trail.sign("user:qa-1", SignatureMeaning::Approver, &key);

    // Chain intact?
    assert!(trail.verify_integrity());

    // FDA-submission package: canonical JSON + detached signature.
    let (payload, tag) = trail.export_package(&key);
    assert!(payload.contains("\"run_id\""));
    assert!(payload.contains("\"signatures\""));
    assert_eq!(tag.len(), 64);

    // The tag is a plain HMAC over the payload, so a third party can check it.
    assert_eq!(hex(&hmac_sha256(&key, payload.as_bytes())), tag);
}
```

### Pinning a run for reproducibility

```rust
use tpt_med_fda::{AuditTrail, ReproducibilityManifest, SignatureMeaning};
use tpt_med_core::{AuditAction, AuditEvent};

fn main() {
    let mut trail = AuditTrail::new("run-2026-09-26-001");

    // Which code, which inputs, which build.
    let manifest = ReproducibilityManifest::new("0.1.0")
        .with_crates([("tpt-med-dicom", "0.1.0"), ("tpt-med-fda", "0.1.0")])
        .with_git_commit(Some("57422ff".into()))
        .with_build_profile("release")
        .with_input("ct_series", b"<serialised slice bytes>");

    let digest = manifest.digest();
    trail.attach_manifest(manifest);

    // The manifest event is on the chain, so a swap is detectable without
    // even re-checking the detached tag.
    assert_eq!(trail.entries.len(), 1);

    trail.append(AuditEvent::new(
        "user:op-1", "mesh", "mesh:abc", AuditAction::Simulate, "stance",
    ));
    trail.sign("user:qa-1", SignatureMeaning::Approver, b"k");
    assert!(trail.verify_integrity());

    let (payload, tag) = trail.export_package(b"k");
    assert!(payload.contains("\"manifest\""));
    assert!(payload.contains("57422ff"));
    println!("manifest {digest}, tag {tag}");
}
```

## API Overview

| Item | Purpose |
|---|---|
| `AuditTrail::new(run_id)` | Open a trail; `run_id` is the chain seed |
| `::append(&mut AuditEvent)` | Append-only; assigns the sequence number and timestamps the entry |
| `::sign(&mut signer, meaning, key)` | Electronic signature over the current state (§11.50) |
| `::verify_integrity() -> bool` | Re-walk the chain; detects edits, truncation, seed splicing |
| `::export_tag(&key) -> String` | Detached HMAC-SHA256 tag over the export payload (hex) |
| `::export_package(&key) -> (String, String)` | Canonical JSON payload plus detached tag |
| `::len()`, `::is_empty()` | Entry count |
| `AuditTrail::run_id` | The chain seed |
| `AuditEntry` | `sequence`, `timestamp`, `actor`, `object_class`, `object_token`, `action`, `reason` |
| `AuditEntry::canonical() -> String` | The exact string the chain commits to |
| `AuditTrail::attach_manifest(ReproducibilityManifest)` | Attach a manifest and append an audited event carrying its digest |
| `AuditTrail::manifest() -> Option<&ReproducibilityManifest>` | The attached manifest, if any |
| `ReproducibilityManifest::new(workspace_version)` | Manifest stamped now; builder-style `with_crate`/`with_crates`/`with_git_commit`/`with_build_profile`/`with_input`/`with_artifact` |
| `ReproducibilityManifest::{canonical, digest, to_json, is_pinned}` | The string a digest covers, that digest, the export rendering, and whether a commit is pinned |
| `InputArtifact::new(label, bytes)` | SHA-256 plus byte length of an input |
| `MANIFEST_SCHEMA_VERSION`, `TOOL_ID` | `1`, and `"tpt-medical"` |
| `ElectronicSignature` | `signer`, `timestamp`, `meaning` |
| `SignatureMeaning` | `Author`, `Reviewer`, `Approver`, `ResponsibleParty` |
| `UtcStamp` | `epoch_seconds`, `nanos`; `now()`, `to_iso8601()` |
| `AuditAction`, `AuditEvent` from `tpt-med-core` | The domain-level events appended here |

## Verification

- **Canonical form** — `AuditEntry::canonical()` is asserted stable and
  unambiguous, including the `|`-separated field order, so two entries cannot
  canonicalise to the same string by moving a delimiter.
- **JSON escaping** — quotes, backslashes, newlines and **control characters
  below 0x20** in `reason` and the token fields are escaped per RFC 8259. A
  raw newline in a reason string would otherwise produce a payload no parser
  accepts, and a payload nobody can parse cannot be signed.
- **Time conversion** — the proleptic-Gregorian `civil_from_days` is asserted
  against known boundaries (leap days, century non-leap years, the 1970 epoch,
  and a far-future date), because a wrong date in a regulatory record is
  indefensible.
- **Chain integrity** — `verify_integrity()` is asserted `true` after appends
  and signatures, and `false` after a retroactive edit, an entry removal, and
  a substituted `run_id`.
- **Signatures are distinguishable** — each `SignatureMeaning` is asserted to
  produce its own wire string (`author`, `reviewer`, `approver`,
  `responsible_party`), which is the §11.50 requirement.
- **Export determinism** — two exports of an unchanged trail are asserted
  byte-identical, so the tag is reproducible and the package is diffable.
- **Tag binding** — the tag is asserted to cover the whole payload: changing a
  single character of the export invalidates it.
- Golden reference datasets: `test-data/golden/regulatory/fda_audit_trail.json`
  and `fda_export_example.json`.

## Known Limitations

- **Non-repudiation needs the operational half.** `attach_anchor` records
  that a payload's digest was submitted to an RFC 3161 authority or a
  transparency log (`tpt-med-audit::AnchorRecord`), and the `ed25519`
  feature of `tpt-med-audit` adds asymmetric signatures — but a full
  non-repudiation *claim* still requires a protected signing key and a
  real anchor service (HSM, TSA, log), which are deployment concerns.
- **No key management.** Generation, storage, rotation and revocation are the
  caller's responsibility.
- **No access control or authentication.** The crate does not verify *who* an
  actor token corresponds to; that belongs in the identity layer. An actor
  token is a claim, not an authenticated identity.
- **Signature coverage narrowing is recorded, not enforced.** Signing, then
  appending, then exporting is permitted and the export reflects the ordering.
  A workflow requiring "sign after the final edit" must enforce that itself.
- **Persistence is entry-level, and local.** `WormLog` gives a write-once,
  chain-verified file journal of the *entries* (per-line digests, `fsync`
  per append, torn-tail and splice detection). It is a single local file:
  no replication, no retention schedule, and signatures/policies remain
  trail-side state the caller re-applies after a reload.
- **No retention or archival policy**, no legal hold, no record retention
  schedule.
- **Per-crate versions are the caller's job to populate correctly.**
  `ReproducibilityManifest::with_crates` already records one version per
  crate, not just the workspace version — but the crate has no way to know
  which version is *true* for a given build; it only stores what it is given.
  The `fda-package` example resolves real per-crate versions from
  `Cargo.lock` at build time (`tpt_med_examples::crate_versions`, generated by
  `examples/build.rs`) rather than stamping every crate with the calling
  binary's own `CARGO_PKG_VERSION`, so the manifest stays correct once crates
  stop versioning as a single unit. A caller that skips this and hardcodes one
  version string for every crate — as this example itself did before — gets a
  manifest that is wrong the day versions diverge, with no error to signal it.
- **The git commit is caller-supplied.** There is no build-time embedding, so a
  build not from a checkout records `None` — which is a legitimate, explicit
  value rather than a gap, but it does mean a released binary cannot be tied
  back to a commit without the operator supplying it.
- **Input digests are only as good as what you feed them.** The manifest hashes
  the bytes you pass; it does not know that a DICOM series is the same anatomy
  re-exported, and it cannot detect a hash fed in by a caller that hashed
  something else.
- **Part 11 compliance is a system property.** This crate implements software
  controls only.

## Related Crates

- [`tpt-med-audit`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/regulatory/tpt-med-audit) — the SHA-256 and HMAC-SHA256 primitives underneath.
- [`tpt-med-core`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-core) — `AuditEvent`, `AuditAction`, and the non-identifying privacy model.
- [`tpt-med-vv40`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/regulatory/tpt-med-vv40) — the credibility assessment that scopes how much evidence a run needs.
- [`tpt-med-surgical-planning`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/surgical/tpt-med-surgical-planning) — produces the operation log that a plan records here.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Changes to the
canonical form, the export schema or the signing semantics are
**behaviour-changing** and require an [RFC](../../../rfcs) — a signature
computed over a different canonical string will not verify against an existing
package. **Never commit PHI**, including in `reason` strings.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA. This
crate implements software controls that support a 21 CFR Part 11 process; a
validated system additionally requires procedural controls — SOPs, operator
training, access control, record retention — that no library can provide.
