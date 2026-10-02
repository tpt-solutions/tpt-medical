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
- **`SoftTissueStructure` / `StructureRisk` — structures-at-risk
  screening**, the first soft-tissue capability: capsules (a segment
  with an envelope radius) registered via `VirtualSurgery::watch_structures`
  and screened at execution against every cut surface in the plan —
  plane, wedge, cylinder and mesh alike — with a per-(structure, cut)
  clearance (signed distance to the removed region minus the envelope;
  negative = breached) in `SurgeryReport::structures_at_risk`. The
  removed-region SDFs are closed form for planes/wedges/cylinders and
  sampled at ≤ half-radius stride for meshes (documented error bound);
  the mesh closure's kerf convention was written so the removed side's
  sign is negative exactly where the apply path resects. Capsule
  stand-ins, not deformable tissue. Verified: a crossing capsule
  breaches at the exact end-sample depth, a parallel one clears by
  exactly distance-minus-radius, a cylinder's radial margin is exact,
  an envelope poking through a kept contour breaches while the same
  structure clears the cylinder it lives inside, and unregistered
  structures produce no entries.
- **`ImplantPlacement` / `GraftReconstruction` — implant component
  placement and bone graft / defect reconstruction.** A closed contour
  mesh defines the component or graft (pose baked into the vertices);
  placement turns every voxel inside the mesh into the component (marker
  value, replaced bone counted as resection) and grafting fills every
  *empty* voxel inside the mesh (a graft fills a defect, never resects).
  Both report the **bone–implant/graft interface area** — shared faces
  with surviving bone, counted from a finished occupancy mask so the
  number cannot depend on scan order — and both are cut-sequenced
  ([`PlanError::CutAfterMove`]/`ModelAlreadySplit`) and mesh-validated at
  build time. No fixation mechanics: placement and interface
  bookkeeping, with stem/cement stress left to `tpt-med-orthopedics`.
  Verified: exact 8 mm³ resection/occupation and 24 mm² interface for an
  embedded 2×2×2 component; a contoured cavity cut-then-grafted with the
  same surface restores the model's voxel count exactly; bone voxels are
  never grafted; invalid meshes and cut-after-move sequencing are
  rejected at build time.
- **`MeshCut` — freeform (anatomically contoured) resections**, completing
  the curved/freeform item: the cut surface is a caller-supplied closed
  triangle mesh (the patient-matched contour no plane, wedge or cylinder
  expresses), kept on the interior or the exterior with the same kerf,
  `rfcs/0011` retention, measurement and audit conventions as the other
  cuts, and a `PlanStep::Mesh`/`VirtualSurgery::mesh` path validated at
  build time (`PlanError::InvalidMesh`). Side determination: unsigned
  closest-triangle distance signed by +x ray parity — which is only
  meaningful for a watertight, consistently wound surface, so
  `MeshCut::validate` enforces closure (every undirected edge shared by
  exactly two triangles), winding (each directed edge once — a flipped
  triangle is caught before it can silently invert a region),
  non-degeneracy and non-emptiness. Verified: a box contour keeps exactly
  its interior voxels and, on a fine grid, reproduces the analytic box
  volume and the kerf-shrunk volume exactly; exterior+retention conserves
  the voxel count; an 80-face icosphere's voxelised interior matches the
  mesh's own exact signed-tetrahedron volume within 5 % (the ideal
  sphere's volume is deliberately *not* the reference — the polyhedron
  sits apothem-deep below it); closest-point-on-triangle is hand-checked
  in its vertex/edge/face regions; and the alignment limit (a mesh edge
  threading voxel-centre rows) is documented and kept out of the test
  geometry.
- **`CylindricalCut` (+ `PlanStep::Cylinder`, `VirtualSurgery::cylinder`):
  the first curved resection surface.** The cylinder of `radius` about an
  axis (`axis_origin` + `axis_direction`) — the shape a reamer or burr
  actually leaves, and the surface a rotational (derotation) osteotomy
  swings about — keeping the core within `radius` (`keep_inside`, e.g. a
  reaming or core decompression) or the annulus outside it. Everything the
  plane cut has carries over: an optional saw-kerf width removing the
  radial slab within `kerf_width / 2` of the wall on both sides,
  `DiscardedSide::RetainAs` (a cylindrical core kept as a named fragment
  for grafting), the same measurement conventions (resection volume;
  depth measured past the kept face — the wall shifted by half the kerf),
  build-time plan validation, and the `cylinder:` audit label. The
  retained-side compaction is now one shared helper for the plane and
  cylinder cuts, so the two cannot drift. Five new tests: exact core
  count/volume/depth, keep-outside resection, kerf slab, offset + non-unit
  axis, retain-and-named-move integration, and the plan-validation
  rejections.
- **Per-fragment addressing, first slice** (`rfcs/0011` option 3,
  implemented following the maintainer's direction to proceed):
  `DiscardedSide::{Resect, RetainAs { name }}` on `OsteotomyCut` — a cut
  can retain its discarded side as a **named, independently addressable
  fragment** — plus `PlanStep::MoveNamed` /
  `VirtualSurgery::move_fragment_named` and build-time validation
  (`PlanError::{UnknownFragment, CutAfterMove, ModelAlreadySplit,
  EmptyRetainedName}`), so `execute` stays infallible. The executor
  tracks named fragments and composes them into the operated model;
  **plans that never use retention are unaffected** (single-fragment
  composition is the identity — every existing plan, log and
  measurement is byte-identical, asserted by the untouched prior
  suite). A retained side is not counted as resection. Fragment
  moves that collide record last-write-wins in the composition (the
  per-fragment measurement report is the honest statement of the
  quantisation). Five new tests: the retain-and-move workflow,
  retention-vs-resection accounting, and the three build-time
  rejections. The RFC's own rule 2/rule 4 contradiction (found while
  implementing) is recorded in the RFC with option 3 superseding the
  original recommendation.
- `OsteotomyCut::kerf_width` and `WedgeCut` (+ `PlanStep::Wedge`,
  `VirtualSurgery::wedge`): **saw-kerf width and multi-plane closed
  wedges**. The kerf removes a slab of the given width centred on the
  cutting plane (the kept boundary shifts outward by half the kerf, and
  the resection measurement grows by the slab — depth is now measured past
  the kept face, which also fixes `keep_negative` cuts reporting a depth
  of 0). A wedge removes exactly the region on the discarded side of both
  planes — the intersection two sequential single-sided cuts cannot
  express — with its own `wedge:` audit label and measurement. Three new
  tests (kerf slab + depth, exact wedge intersection through the plan,
  kerf-offset wedges).
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

### Planned
- Implicit/spline-contoured surfaces (expressible today by tessellating
  into `MeshCut`; no native implicit cut).
- Independently addressable implant fragments (implants currently live in
  the operated fragment's grid under a marker value, with the pose baked
  into the contour vertices).
- Deformable soft tissue: the structures-at-risk screen is capsule
  stand-ins against the cut field, not tissue that moves.
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
