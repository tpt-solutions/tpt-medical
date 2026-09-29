# Changelog

All notable changes to `tpt-med-fda` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- `SignaturePolicy` (permissive vs `RequireSignatureAfterLastEdit`),
  `ReasonPolicy` (free text vs structured reason codes), and the
  policy-checked `checked_append` / `export_package_checked`: **workflow
  discipline enforcement** — a site can require sign-after-final-edit and
  structured reason codes, with `PolicyError` naming the violation. The
  permissive default preserves the historic append/sign/export ordering.
- Crate README with an explicit **regulatory mapping table** (§11.10(e), audit
  trail integrity, §11.50, §11.10(k)) and a statement that Part 11 compliance
  is a property of a *system* — procedures, training, retention, access
  control — and that no library can supply those.

- `WormLog` / `WormError` (the new `worm` module): an **append-only
  persistence layer with WORM semantics** for an `AuditTrail` — the journal
  is created once (`create_new`, never reopened for a second lifetime),
  each entry is appended with its running hash-chain digest and `fsync`ed
  before the call returns, and `open` verifies the header, per-line chain
  digests, and sequence continuity before handing the log back. A crash
  mid-append leaves a line without its newline and is refused as
  `TornTail` rather than silently truncated; a retroactive edit or a
  record spliced from another run is refused as `ChainBroken` at the first
  bad link. `into_trail` rebuilds a live, self-verifying trail (entries
  keep their recorded timestamps; signatures and policies are trail-side
  state re-applied by the caller). Seven new tests: round-trip into a
  verifying trail, process-restart continuation, retroactive-edit
  detection, torn-tail refusal, `create_new` on an existing path,
  cross-run splice detection, and free-text escaping round-trip.

### Planned
- External anchoring of the detached tag (HSM, transparency log, RFC 3161) to
  close the non-repudiation gap identified in RFC 0003.

### Notes
- **Changes to the canonical form, the export schema or the signing semantics
  are behaviour-changing and require an RFC.** A signature computed over a
  different canonical string will not verify against an existing package, so
  the golden export file in `test-data/golden/regulatory/` is the canary.
- `reason` is free text and is **the caller's policy surface**: a regulated
  deployment must ensure it carries no PHI. This is documented rather than
  enforced, because the correct policy is site-specific.
- **Never commit PHI**, including in `reason` strings.
- In-memory trails are the caller's to persist or not; the `worm` module's
  `WormLog` is the write-once journal when durability is wanted.

## [0.1.0] - 2026-09-22

### Added
- **§11.10(e)** — `AuditEntry` records a non-identifying `actor` token (from the
  `tpt-med-core` privacy model) and a `UtcStamp` taken at `append` time, so
  operator identity and time are captured on every operation.
- **`AuditTrail::new(run_id)`** — one trail per simulation run, with `run_id` as
  the chain seed, so two runs can never have interchangeable digests.
- **`append(AuditEvent)`** — strictly append-only. There is no remove, no
  update and no entry-list setter; modifying a record means appending a new
  `Modify` event with a reason, which is what §11.10(e) requires. A 0-based
  monotonic `sequence` is assigned at append.
- **§11.50** — `ElectronicSignature { signer, timestamp, meaning }` and
  `SignatureMeaning` (`Author`, `Reviewer`, `Approver`, `ResponsibleParty`).
  The meanings are an enum, not free text, because §11.50 requires them to be
  distinguishable.
- **Audit trail integrity** — entries are chained via
  `digest[i] = H(run_id ‖ digest[i-1] ‖ canonical(entry))` over
  `AuditEntry::canonical()`; `verify_integrity()` re-walks the whole chain and
  detects edits, truncation and seed splicing.
- **§11.10(k)** — `export_package(key)` returns canonical JSON plus a detached
  HMAC-SHA256 tag, and `export_tag(key)` returns the tag alone. The JSON is
  hand-serialised with RFC 8259 escaping so the output is byte-stable and the
  tag is reproducible; the tag covers the entries, their digests **and** the
  signature manifestations together.
- `UtcStamp { epoch_seconds, nanos }` with `now()` and `to_iso8601()`
  (`YYYY-MM-DDThh:mm:ssZ`), and a local proleptic-Gregorian
  `civil_from_days` (Hinnant's algorithm) — so there is no time-crate
  dependency and no timezone database to get wrong.
- `len()` and `is_empty()`.
- **Reproducibility manifest** (`manifest` module) —
  `ReproducibilityManifest` pins *which code* (workspace version, per-crate
  versions, git commit), *which inputs* (`InputArtifact`: SHA-256 plus byte
  length, so a truncated input is distinguishable from a different one),
  *which build* (profile), and *when*. `canonical()` is the exact string the
  digest is taken over — deliberately **not** the JSON rendering, so adding an
  export field cannot silently change what a signature covers.
- `AuditTrail::attach_manifest(ReproducibilityManifest)` attaches a manifest
  **and appends an audit event carrying its digest**, so a manifest swapped
  after the fact is detectable from the chain alone. `AuditTrail::manifest()`
  reads it back.
- The manifest appears in `export_package` under a `"manifest"` key **only when
  one is attached**, so runs that do not use it produce byte-identical exports
  and packages signed before this feature existed still verify.
- `InputArtifact`, `MANIFEST_SCHEMA_VERSION` and `TOOL_ID` re-exported from the
  crate root.
- No storage, no network, no external time source. Not a PHI store: actor
  tokens are non-identifying by construction.

### Verification
- **Canonical form** — asserted stable and unambiguous, including the
  `|`-separated field order, so two entries cannot canonicalise to the same
  string by moving a delimiter.
- **JSON escaping** — quotes, backslashes, newlines and **control characters
  below 0x20** in `reason` and the token fields are escaped per RFC 8259. A raw
  newline in a reason string would otherwise produce a payload no parser
  accepts, and a payload nobody can parse cannot be signed.
- **Time conversion** — the civil-date algorithm is asserted against known
  boundaries (leap days, century non-leap years, the 1970 epoch, a far-future
  date), because a wrong date in a regulatory record is indefensible.
- **Chain integrity** — `true` after appends and signatures; `false` after a
  retroactive edit, an entry removal, and a substituted `run_id`.
- **Signature distinctness** — each `SignatureMeaning` produces its own wire
  string (`author`, `reviewer`, `approver`, `responsible_party`).
- **Export determinism** — two exports of an unchanged trail are byte-identical,
  so the tag is reproducible and the package is diffable.
- **Tag binding** — changing a single character of the export invalidates the
  tag.
- **Manifest** — `digest_is_stable_for_identical_inputs`,
  `input_order_does_not_change_the_digest` (inputs are sorted by label),
  `a_changed_input_changes_the_digest`, and `every_field_participates_in_the_digest`
  (git commit, workspace version, crate versions, input digests and the
  timestamp each change the digest independently).
  `absent_git_commit_is_explicit_not_omitted` pins that a missing optional field
  renders as `-` in the canonical form and `null` in JSON, so it can never be
  confused with an empty string.
- **Backward compatibility** — `export_without_a_manifest_is_unchanged` asserts
  a run with no manifest emits no `manifest` key at all, so packages signed
  before this feature existed still verify.
- **Manifest binding** — `a_changed_manifest_changes_the_export_tag` asserts
  the manifest body is covered by the detached HMAC, and
  `attaching_a_manifest_records_an_event_and_exports_it` asserts the recorded
  event's object token equals the manifest digest, so a swap is detectable from
  the chain alone.
- Golden datasets `test-data/golden/regulatory/fda_audit_trail.json` and
  `fda_export_example.json`.

### Known limitations
- Symmetric keys give integrity and attribution, **not** non-repudiation
  (RFC 0003).
- No key management.
- No access control or authentication: an actor token is a claim, not an
  authenticated identity.
- Signature coverage narrowing is recorded, not enforced.
- In-memory only; no retention, archival or legal-hold policy.
- Part 11 compliance is a system property; this crate implements software
  controls only.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
