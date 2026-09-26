# tpt-med-core

Patient models, anatomical taxonomy, and HIPAA-safe audit traits — the domain
vocabulary every other crate in the stack speaks.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--core-orange)](https://crates.io/crates/tpt-med-core)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--core-blue)](https://docs.rs/tpt-med-core)

| | |
|---|---|
| **Layer** | `core` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-geometry`](../tpt-med-geometry), [`tpt-med-units`](../tpt-med-units) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Two problems sit at the bottom of every medical simulation stack and neither
is numerical:

1. **Vocabulary.** What is a "transepicondylar axis"? Is this region a
   `Femur` or a "femur"? A stringly-typed taxonomy is unmaintainable, and
   worse, it is *unauditable*: regulatory tooling must be able to enumerate
   every region a model could contain, exhaustively.
2. **Privacy.** The moment a patient identifier reaches a log line, a
   `Display` impl, or an error message, the whole stack is in PHI territory.
   That failure mode is easy to introduce and hard to audit.

This crate fixes both by making the taxonomy a **closed enum** and by making
PHI leakage require an explicit, reviewable `AuditSubject` implementation.

## Features

- **Opaque patient tokens** — `PatientId` is a newtype over a non-cryptographic
  hash. The source identifier is discarded, never stored, and never rendered.
  `Display` emits `patient:<16 hex digits>`, which contains no PHI.
- **Demographics with a body-weight load scale** — `Demographics::body_weight_force()`
  and `PatientModel::body_weight_force()` implement the "×N body weight"
  loading convention used across the stack.
- **Closed anatomical taxonomy** — `BoneType`, `VesselType`, `AortaSegment`,
  `OrganType`, `SoftTissueType`, `ImplantType`, `AnatomicalRegion`, plus
  `Landmark` and `AnatomicalModel` with region/landmark lookup by stable key.
- **Audit traits** — the `AuditSubject` trait, the `AuditAction` vocabulary
  (create/modify/delete/approve/reject/export/simulate), and `AuditEvent`,
  which the regulatory layer (`tpt-med-fda`) timestamps, chains and signs.
- No allocation-heavy structures, no `unsafe`, no I/O.

## Conventions

- **Region enums are intentionally closed.** Adding an anatomical region is a
  semver-minor change that requires an RFC, because downstream regulatory
  tooling must enumerate regions exhaustively.
- **De-identification is not a security boundary.** `PatientId::from_hash` uses
  `DefaultHasher` — it prevents *accidental* PHI leakage through logs and error
  paths, it does not resist an adversary. Any mapping table between real
  identities and tokens is itself PHI under your institution's HIPAA/GDPR
  controls.
- **Audit tokens must be coarse.** `Demographics::audit_token()` renders as
  `demographics:age~70s:bmi-band=normal` — decade bucket plus BMI band, never a
  precise age, weight or height.

## Usage

```rust
use tpt_med_core::{AuditAction, AuditEvent, AuditSubject, Demographics, PatientId,
                   PatientModel, Sex};

fn main() {
    // The raw identifier is hashed and discarded.
    let id = PatientId::from_hash("mrn-12345");
    assert!(!id.to_string().contains("12345"));
    assert_eq!(id.to_string(), format!("patient:{}", id.to_hex()));

    let demo = Demographics {
        age_years: 70,
        sex: Sex::Male,
        weight_kg: 80.0,
        height_cm: 175.0,
    };
    // 80 kg * 9.81 m/s^2
    assert!((demo.body_weight_force() - 784.8).abs() < 1e-9);

    let patient = PatientModel::new(id).with_demographics(demo);
    assert!(patient.body_weight_force().is_some());

    // Audit output is de-identified by construction.
    assert_eq!(demo.audit_token(), "demographics:age~70s:bmi-band=normal");
    assert_eq!(patient.audit_token(), id.to_string());

    // Domain-level events; tpt-med-fda adds timestamps, chaining and signing.
    let event = AuditEvent::new(
        "user:op-1",
        "patient_model",
        id.to_string(),
        AuditAction::Simulate,
        "stance screening",
    );
    assert_eq!(event.action.to_string(), "simulate");
}
```



## API Overview

| Item | Purpose |
|---|---|
| `PatientId` | Opaque token; `from_hash`, `ANONYMOUS`, `to_hex`, `Display` |
| `Demographics` | `age_years`, `sex`, `weight_kg`, `height_cm`; `body_weight_force()` |
| `Sex` | `Male` / `Female` / `Other` |
| `PatientModel` | `id`, `demographics`, `anatomy`; `new`, `with_demographics`, `body_weight_force` |
| `AuditSubject` | Trait: `audit_token() -> String`, non-identifying and stable |
| `AuditAction` | `Create`, `Modify`, `Delete`, `Approve`, `Reject`, `Export`, `Simulate` |
| `AuditEvent` | `actor_token`, `object_class`, `object_token`, `action`, `reason`; `new` |
| `BoneType` | `Femur`, `Tibia`, `Fibula`, `Patella`, `Pelvis`, `Humerus`, `Radius`, `Ulna`, `Vertebra{level}`, `Skull`, `Mandible`, `Custom(name)` |
| `VesselType` | `Artery{name}`, `Vein{name}`, `Aorta{segment}`, `Coronary{branch}` |
| `AortaSegment` | `Root`, `Ascending`, `Arch`, `DescendingThoracic`, `Abdominal` |
| `OrganType` | `Heart`, `Lung{left}`, `Liver`, `Kidney{left}`, `Brain` |
| `SoftTissueType` | `ArterialWall`, `Skin`, and further soft-tissue classes |
| `ImplantType` | Implant classes the device crates reason about |
| `AnatomicalRegion` | Region enum; `key()` gives the stable lookup/export string |
| `Landmark` | Named anatomical landmark with position |
| `AnatomicalModel` | `add_region`, `region(key)`, `add_landmark`, `landmark(key)` |
| `Result<T, E>` | Crate-level result alias |

## Verification

- `id_is_stable_but_opaque` — the same local ID always yields the same token,
  different IDs do not collide, and the `Display` output never contains the
  source string.
- `body_weight_scaling` — pins the 80 kg → 784.8 N conversion exactly, so the
  loading convention cannot drift.
- `action_display` / `event_fields` — the audit vocabulary and its stable
  lowercase wire form, which is what the hash chain commits to.

The privacy properties here are structural rather than tested after the fact:
there is no field on `PatientModel` that can hold a raw identifier, and no
`Debug` or `Display` impl that renders one.

## Known Limitations

- **De-identification is not a security boundary.** `PatientId::from_hash` uses
  `DefaultHasher`, which prevents *accidental* PHI leakage through logs and
  error paths but does not resist an adversary. Any mapping table between real
  identities and tokens is itself PHI under your institution's controls.
- **No persistence.** `PatientModel` and `AnatomicalModel` are in-memory
  structures; serialisation, storage and access control are the caller's.
- **No landmark detection.** `Landmark` and `AnatomicalRegion` are containers;
  finding the epicondyles or the transepicondylar axis in a CT is not
  implemented here.
- **No validation of anatomical consistency.** A model can be built with
  overlapping or anatomically impossible regions; nothing cross-checks that a
  `Femur` region does not sit inside a `Liver` region.
- **Demographics are coarse by design.** `Demographics` carries only age, sex,
  mass and stature, which is enough for body-weight scaling and not enough for
  anything requiring a body composition model.
- **Closed taxonomies are a maintenance cost.** Adding a region is an RFC,
  which is the right friction, but it means a user modelling a structure the
  taxonomy does not know must fall back to `Custom(String)` and accept the loss
  of exhaustive enumeration.

## Related Crates

- [`tpt-med-fda`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/regulatory/tpt-med-fda) — timestamps, chains and HMAC-signs the `AuditEvent`s defined here (RFC 0003).
- [`tpt-med-geometry`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-geometry) — `Aabb` region extents and `Vec3` landmark positions.
- [`tpt-med-units`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-units) — `Force` returned by body-weight scaling.
- [`tpt-med-bone`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-bone) — consumes `BoneType` for reference material properties.
- [`tpt-med-surgical-planning`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/surgical/tpt-med-surgical-planning) — landmarks feed osteotomy and sizing.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Adding an
anatomical region variant is a semver-minor change and requires an RFC. New
types that enter an audit log must implement `AuditSubject`.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Privacy & Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. Real patient
data must never enter this stack; all shipped test data is synthetic.
