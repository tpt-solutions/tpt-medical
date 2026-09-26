# Changelog

All notable changes to `tpt-med-core` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README documenting the privacy model explicitly, including the
  distinction between de-identification (accidental PHI leakage) and a
  security boundary.

### Notes
- Adding a variant to any anatomical taxonomy enum (`BoneType`, `VesselType`,
  `OrganType`, `SoftTissueType`, `ImplantType`, `AnatomicalRegion`,
  `AortaSegment`) is a **semver-minor** change that **requires an RFC**,
  because downstream regulatory tooling must enumerate regions exhaustively.
- `PatientId::from_hash` uses `DefaultHasher` and is a de-identification
  measure, **not** a security control. Replacing it with a keyed hash is
  semver-minor but changes every previously issued token, so it also needs an
  RFC and a migration note.

## [0.1.0] - 2026-09-22

### Added
- **Patient models** — `PatientId` (opaque, non-printing token) with
  `from_hash`, `ANONYMOUS` and `to_hex`; `Demographics`
  (`age_years`, `sex`, `weight_kg`, `height_cm`) with `body_weight_force()`;
  `Sex`; and `PatientModel` with `new`, `with_demographics` and
  `body_weight_force()`.
- **Anatomical taxonomy** — closed enums `BoneType`, `VesselType`,
  `AortaSegment`, `OrganType`, `SoftTissueType`, `ImplantType` and
  `AnatomicalRegion`, plus `Landmark` and `AnatomicalModel` with
  `add_region`/`region(key)` and `add_landmark`/`landmark(key)` lookup.
- **HIPAA-safe audit traits** — the `AuditSubject` trait
  (`audit_token() -> String`), the `AuditAction` vocabulary
  (`Create`, `Modify`, `Delete`, `Approve`, `Reject`, `Export`, `Simulate`),
  and `AuditEvent` with a `new` constructor.
- `Result<T, E>` crate-level result alias.

### Privacy
- `PatientId` never stores or renders the source identifier; `Display` emits
  `patient:<16 hex digits>`, which contains no PHI. There is no field on
  `PatientModel` capable of holding a raw identifier.
- `AuditSubject` for `Demographics` renders decade and BMI-band buckets only
  (`demographics:age~70s:bmi-band=normal`) — never a precise age, weight or
  height.
- Actor tokens are stable across processes so audit trails stay joinable
  without being identifying.

### Verification
- `id_is_stable_but_opaque` — the same local ID yields the same token, different
  IDs do not collide, and the `Display` output never contains the source
  string.
- `body_weight_scaling` — pins the 80 kg → 784.8 N conversion exactly.
- `action_display` and `event_fields` — the stable lowercase audit vocabulary
  that the hash chain commits to.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
