# Changelog

All notable changes to `tpt-med-surgical-planning` are documented here. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README explaining why the audit trail is the point of this crate: the
  planning primitives are not hard, but producing the operation log *by
  construction* rather than reconstructing it from screenshots afterwards is.

- `execute_with_report` + `SurgeryReport` (`CutMeasurement`,
  `MoveMeasurement`, `StepMeasurement`): **measurement reporting** — each
  cut records its resection volume (discarded voxels × voxel volume) and
  cut depth (deepest discarded voxel centre below the plane); each move
  records the prescribed vs achieved centroid displacement, whose
  difference is the achieved alignment error (sub-voxel plans quantise
  onto the grid, and the report is what surfaces that). The report is
  index-aligned with the audit log and carries the total resection volume;
  `OsteotomyCut::apply_measured` exposes the per-cut numbers standalone.
  Three new tests (volume/depth on a plane cut, alignment error for a
  grid-aligned and a sub-voxel translation, multi-cut totals).

- `rfcs/0011-per-fragment-addressing.md` (Draft): the design for
  per-fragment addressing — `move_fragment_named`, fragment identity as the
  producing cut's name, split/keep semantics for later cuts, and the
  `CutAfterMove` rejection — scoped but not implemented, per the RFC
  process.

### Planned
- Curved and freeform resections, saw-kerf width, and multi-plane wedges.
- Implant component placement with a bone–implant interface, and bone graft or
  defect reconstruction.
- Soft-tissue structures, so a plan can be checked for collateral damage to
  ligaments, capsules and neurovascular bundles.

### Notes
- `NaN` is the sentinel for a discarded voxel, and every operation skips `NaN`
  voxels rather than treating them as zero. This is what keeps `count_above`
  honest: a resected model reports the resected volume.
- Discarded voxels are compacted to the bounding box of the kept region, so
  fragment volumes stay correct rather than carrying a full-size field of dead
  cells.
- Fragment transforms are **rigid** (rotation plus translation). No scaling and
  no shear — a plan that could scale bone would not be a surgical plan.
- Adding a `PlanStep` variant is breaking for downstream exhaustive matches.

## [Unreleased]

### Added
- Crate README explaining why the audit trail is the point of this crate: the
  planning primitives are not hard, but producing the operation log *by
  construction* rather than reconstructing it from screenshots afterwards is.

- `rfcs/0011-per-fragment-addressing.md` (Draft): the design for
  per-fragment addressing — `move_fragment_named`, fragment identity as the
  producing cut's name, split/keep semantics for later cuts, and the
  `CutAfterMove` rejection — scoped but not implemented, per the RFC
  process.

### Planned
- Curved and freeform resections, saw-kerf width, and multi-plane wedges.
- Implant component placement with a bone–implant interface, and bone graft or
  defect reconstruction.
- Soft-tissue structures, so a plan can be checked for collateral damage to
  ligaments, capsules and neurovascular bundles.
- Measurement reporting: resection volumes, cut depths and achieved alignment
  errors, recorded alongside the steps in the audit log.

### Notes
- `NaN` is the sentinel for a discarded voxel, and every operation skips `NaN`
  voxels rather than treating them as zero. This is what keeps `count_above`
  honest: a resected model reports the resected volume.
- Discarded voxels are compacted to the bounding box of the kept region, so
  fragment volumes stay correct rather than carrying a full-size field of dead
  cells.
- Fragment transforms are **rigid** (rotation plus translation). No scaling and
  no shear — a plan that could scale bone would not be a surgical plan.
- Adding a `PlanStep` variant is breaking for downstream exhaustive matches.

## [0.1.0] - 2026-09-22

### Added
- **`VoxelModel { dims, spacing, origin, values }`** — the surgical substrate:
  a scalar voxel field (HU or a mask label) with millimetre pitch and a
  patient-space origin, deliberately the same data structure the imaging crates
  produce. Accessors `index`, `center` and `count_above`; the index is
  `(z·ny+y)·nx+x`, matching `tpt-med-meshing`.
- **`OsteotomyCut { plane, fragment_name, keep_positive }`** — a plane
  resection keeping the `signed_distance ≥ 0` side (`≤ 0` when
  `keep_positive` is false). Discarded voxels are set to `NaN` and the kept
  region is compacted to its bounding box. Every cut carries a
  `fragment_name` so the audit log is readable by a human reviewer, not merely
  parseable.
- **`FragmentTransform { rotation_axis, rotation_angle, translation }`** —
  rigid reposition, with `apply_to_point` and `apply_to_model`. A tibial
  tubercle distalisation or a Le Fort advancement is one call.
- **`PlanStep`** — `Cut(OsteotomyCut)` or `Move(FragmentTransform)`.
- **`VirtualSurgery`** — `new(base)`, `cut(..)` and `move_fragment(..)`
  appending steps in plan order (builder style, returning `&mut Self`), and
  `execute() -> (VoxelModel, Vec<String>)` returning the operated model **and
  the audit log in execution order**.
- `base_model()` returning the pre-operative model for side-by-side planning
  views, and `fragments() -> BTreeMap<usize, String>` returning the recorded
  labels ordered and deduplicated.
- No dependencies beyond `tpt-med-geometry`; `#![forbid(unsafe_code)]`.

### Verification
- **Cut correctness** — a cut at a known coordinate leaves exactly the expected
  voxel count, and every retained voxel is on the correct side of the plane.
- **Compaction** — the returned model's bounding box tightly bounds the kept
  region; a compaction bug leaves dead space that is invisible until a volume
  or a screenshot is wrong.
- **`NaN` discipline** — discarded voxels are `NaN`, and every subsequent
  operation is asserted to skip them rather than treating them as zero.
- **Transform exactness** — a π/2 rotation about `z` maps `X` to `Y`, and a
  zero transform is the identity, pinning the right-handed axis–angle
  convention.
- **Volume conservation** — a rigid `FragmentTransform` does not change the
  non-`NaN` voxel count, as a strict equality.
- **Audit log completeness** — log length equals step count, entries appear in
  plan order, and each is prefixed by its step kind.
- **Idempotence** — `execute` called twice returns identical models, so the
  plan does not consume its own steps.
- **Base model immutability** — `base_model()` is unchanged after `execute`.

### Known limitations
- Plane cuts only; no curved resection, saw-kerf width or multi-plane wedge.
- Rigid fragments only; no implant placement, graft or reconstruction.
- **A `Move` applies to the whole assembled model**, not to a single named
  fragment, so repositioning two fragments independently needs two
  `VirtualSurgery` invocations. A deliberate v0 simplification and a known gap.
- No soft tissue, so no collateral-damage check.
- The audit log records the steps, not the resulting volumes, depths or
  alignment errors.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
