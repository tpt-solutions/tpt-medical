# tpt-med-surgical-planning

Osteotomy cuts and virtual surgery on voxel anatomical models — plane resections,
rigid fragment reposition, and an auditable operation log.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--surgical--planning-orange)](https://crates.io/crates/tpt-med-surgical-planning)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--surgical--planning-blue)](https://docs.rs/tpt-med-surgical-planning)

| | |
|---|---|
| **Layer** | `surgical` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-geometry`](../../core/tpt-med-geometry) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Surgical planning software is expensive, closed, and — for the tools that do
exist — usually requires uploading the patient's CT to a vendor's cloud. The
planning primitives themselves are not hard: resect on a plane, reposition a
fragment, record what you did. What is hard, and what makes this crate worth
having, is the **audit trail**.

A virtual surgery here is an ordered list of `PlanStep`s. Executing it returns
the operated model *and* a log of every step in execution order. That log is
exactly the artefact a 21 CFR Part 11 pipeline needs, and it is produced by
construction rather than reconstructed afterwards from screenshots.

The substrate is a voxel model, matching the CT-derived geometry the rest of
the stack already produces. No surface mesh, no segmentation watershed, no
BSP tree — a resection is a half-space test, which is exact and fast on a
voxel grid.

## Features

- **`VoxelModel`** — dimensions, pitch, patient-space origin, and a scalar
  field (HU or a mask label). This is the surgical substrate, and it is
  deliberately the same data structure the imaging crates produce.
- **`OsteotomyCut`** — a plane resection keeping the `signed_distance ≥ 0`
  side. Discarded voxels are set to `NaN` and the kept region is compacted to
  its bounding box, so fragment volumes stay correct rather than carrying a
  full-size field of dead cells.
- **`FragmentTransform`** — rigid rotation about an axis plus translation, with
  `apply_to_point` and `apply_to_model`. A tibial tubercle distalisation or a
  Le Fort advancement is one call.
- **`VirtualSurgery`** — `cut` and `move_fragment` in plan order, `execute`
  returning `(operated_model, audit_log)`, plus `base_model` for side-by-side
  pre/post views and `fragments` for the recorded labels.
- **Named fragments** — every cut carries a `fragment_name` so the audit log
  is readable by a human reviewer, not just parseable.
- Voxel index `(z·ny + y)·nx + x`, matching `tpt-med-meshing`.

## Conventions

- `spacing` in **mm** per axis; `origin` is the patient-space position of
  voxel `(0,0,0)`.
- A cut keeps the region where `plane.signed_distance(voxel_centre) ≥ 0` when

## Usage

```rust
use tpt_med_geometry::{Plane, Vec3};
use tpt_med_surgical_planning::{FragmentTransform, OsteotomyCut, VirtualSurgery, VoxelModel};

fn main() {
    // A small 8^3 synthetic bone block, 1 mm voxels, 300 HU.
    let nx = ny = nz = 8;
    let model = VoxelModel {
        dims: (nx, ny, nz),
        spacing: (1.0, 1.0, 1.0),
        origin: Vec3::ZERO,
        values: vec![300.0; nx * ny * nz],
    };
    assert_eq!(model.count_above(200.0), 512);

    // Plan: distal femoral resection at z = 3.5 mm, keeping the proximal part.
    let mut plan = VirtualSurgery::new(model);
    plan.cut(OsteotomyCut {
        plane: Plane::from_point_normal(Vec3::new(0.0, 0.0, 3.5), Vec3::new(0.0, 0.0, -1.0))
            .unwrap(),
        fragment_name: "femoral-resection".into(),
        keep_positive: true,
    });

    // Then reposition the fragment: 3 degrees about the transepicondylar
    // axis, plus 4 mm distally.
    plan.move_fragment(FragmentTransform {
        rotation_axis: Vec3::new(0.0, 1.0, 0.0),
        rotation_angle: 3.0_f64.to_radians(),
        translation: Vec3::new(0.0, 0.0, -4.0),
    });

    // Execute -> operated model plus the audit trail, in execution order.
    let (operated, log) = plan.execute();
    assert_eq!(log.len(), 2);
    assert!(log[0].starts_with("cut:"));
    assert!(log[1].starts_with("move:"));
    assert!(operated.count_above(200.0) < 512); // bone was resected away

    // Recorded fragments, and the untouched base for side-by-side views.
    assert_eq!(plan.fragments().len(), 2);
    assert_eq!(plan.base_model().dims, (nx, ny, nz));
}
```

## API Overview

| Item | Purpose |
|---|---|
| `VoxelModel { dims, spacing, origin, values }` | The surgical substrate; index `(z·ny+y)·nx+x` |

## Verification

- **Cut correctness** — a cut at a known coordinate is asserted to leave
  exactly the expected voxel count, and every retained voxel is asserted to be
  on the correct side of the plane.
- **Compaction** — after a cut, the returned model's bounding box is asserted
  to tightly bound the kept region. A compaction bug would leave dead space
  that is invisible until a volume or a screenshot is wrong.
- **`NaN` discipline** — discarded voxels are `NaN`, and every subsequent
  operation (cut, transform, `count_above`) is asserted to skip them rather
  than treating them as zero.
- **Transform exactness** — a rotation of exactly π/2 about `z` maps `X` to
  `Y`, and a zero transform is the identity, asserted to floating-point
  tolerance. This pins the rotation convention (right-handed, axis–angle).
- **Volume conservation** — a rigid `FragmentTransform` does not change the
  count of non-`NaN` voxels. A transform that loses voxels is a real bug, so
  this is a strict equality test.
- **Audit log completeness** — the log length equals the step count, entries
  appear in plan order, and each is prefixed by its step kind.
- **Idempotence** — `execute` called twice on the same plan returns identical
  models, since the plan must not consume its own steps.
- **Base model immutability** — `base_model()` is unchanged after `execute`.

## Known Limitations

- **Plane cuts only.** No freeform or curved resection, no saw-kerf width, no
  multi-plane wedge. Real surgical planning wants all three.
- **Rigid fragments only.** No implant component placement with a bone-implant
  interface, no bone graft, no defect reconstruction.
- **Fragment transforms act on the whole current model.** A `Move` applies to
  the assembled model rather than to a single named fragment, so a plan that
  repositions two different fragments independently needs two
  `VirtualSurgery` invocations or an extension to `PlanStep`. This is a
  deliberate v0 simplification and a known gap.
- **No soft tissue.** Only the bone voxel model is planned; ligaments,
  capsules and neurovascular structures are absent, so no plan can be checked
  for collateral damage.
- **No measurement or intersection reporting.** The audit log records the
  *steps*, not the resulting resection volumes, cut depths or alignment
  errors — those are computed by the caller.

## Related Crates

- [`tpt-med-implant-sizing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/surgical/tpt-med-implant-sizing) — chooses the component once the cuts are planned.
- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) and [`tpt-med-meshing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-meshing) — produce the base voxel model.
- [`tpt-med-orthopedics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-orthopedics) — evaluates the fixated result.
- [`tpt-med-fda`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/regulatory/tpt-med-fda) — signs the audit trail this crate produces.
- [`tpt-med-geometry`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-geometry) — `Plane` and `Vec3`.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). New plan step
types or a per-fragment addressing model require an [RFC](../../../rfcs).
**Never commit real patient data** — all fixtures are synthetic.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. **Not for
clinical decision-making.** This is an engineering and research tool, not a
surgical planning system.

| `VoxelModel::index(x, y, z) -> Option<usize>` | Bounds-checked flat index |
| `VoxelModel::center(x, y, z) -> Vec3` | Patient-space voxel centre (mm) |
| `VoxelModel::count_above(threshold) -> usize` | Voxel count above a threshold; `NaN` voxels never count |
| `OsteotomyCut { plane, fragment_name, keep_positive }` | Plane resection; discarded voxels become `NaN`, kept region is compacted |
| `OsteotomyCut::apply(&VoxelModel) -> VoxelModel` | The resection |
| `FragmentTransform` | `rotation_axis`, `rotation_angle`, `translation` — rigid only |
| `FragmentTransform::apply_to_point(Vec3) -> Vec3` | Transform a point |
| `FragmentTransform::apply_to_model(&VoxelModel) -> VoxelModel` | Transform a whole model |
| `PlanStep` | `Cut(OsteotomyCut)` or `Move(FragmentTransform)` |
| `VirtualSurgery::new(base)` | Start a plan on the pre-operative model |
| `VirtualSurgery::cut(..)`, `::move_fragment(..)` | Append steps in plan order (builder style, returns `&mut Self`) |
| `VirtualSurgery::execute() -> (VoxelModel, Vec<String>)` | Operated model plus the audit log in execution order |
| `VirtualSurgery::base_model() -> &VoxelModel` | Pre-operative model, for side-by-side views |
| `VirtualSurgery::fragments() -> BTreeMap<usize, String>` | Recorded fragment labels, ordered and deduplicated |
| `Plane`, `Vec3` from `tpt-med-geometry` | The cutting plane and geometry types |

  `keep_positive` is true, and `≤ 0` otherwise.
- `NaN` is the sentinel for a discarded voxel, and every operation skips `NaN`
  voxels rather than treating them as zero. This keeps `count_above` honest:
  a resected model reports the resected volume, not a volume full of
  `NaN`-compared-as-false noise.
- Transforms are **rigid** (rotation + translation). No scaling, no shear —
  a plan that could scale bone would not be a surgical plan.
- `execute` is pure with respect to the plan: it can be called repeatedly and
  always returns the same result from the same base model.
