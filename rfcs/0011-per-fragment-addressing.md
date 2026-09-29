# RFC 0011: Per-Fragment Addressing in Surgical Plans

- **Status:** Draft
- **Started:** 2026-09-30
- **Crates:** `tpt-med-surgical-planning`

## Summary

A `VirtualSurgery` plan today applies every [`PlanStep::Move`] to **the whole
current model** — there is no way to reposition one named fragment of an
assembled model independently, which is the operation real virtual surgery
exists for (a tibial tubercle moved distally while the shaft stays put; two
Le Fort segments advanced by different amounts). This RFC designs the
smallest change that gives plans per-fragment addressing without invalidating
existing plans or the audit-log contract.

## Motivation

The crate's own README names this "the largest known gap": a plan that
repositions two fragments independently currently needs two `VirtualSurgery`
invocations with hand-cropped models, and the resulting audit logs cannot
describe one operation. The `cut` step already produces a *named* fragment
(`OsteotomyCut::fragment_name`) — the name exists for the log but addresses
nothing. The design question is where fragment identity lives, because
voxel models have no concept of connectivity: a "fragment" is whatever a cut
kept, and subsequent cuts can split it further.

## Guide-level explanation

Fragments become first-class, named state on the plan:

```rust
pub struct VirtualSurgery { /* ... */ }

impl VirtualSurgery {
    /// Moves one named fragment. The fragment is the kept side of an
    /// earlier `OsteotomyCut` with `fragment_name == name` — errors with
    /// `PlanError::UnknownFragment` if no such cut has run.
    pub fn move_fragment_named(
        &mut self,
        name: impl Into<String>,
        transform: FragmentTransform,
    ) -> &mut Self;
}
```

Execution semantics, in one sentence per rule:

1. **Fragment identity is the name of the cut that produced the kept side.**
   A cut with `fragment_name = "distal"` and `keep_positive = false` means
   the negative side *is* `distal`. The discarded side is unnamed unless a
   second cut names it.
2. **A later cut splits every fragment it intersects.** If cut `B` runs on
   a model already containing fragment `A`, then what remains of `A` after
   `B` keeps `A`'s name — unless `B` runs with `fragment_name = C`, in which
   case the kept side of `B` is `C` and the rest of `A` stays `A`.
3. **Moves compose in plan order.** A fragment's world position is the
   composition of every `move_fragment_named` that targeted it, applied
   after every cut that touched it. Moving a fragment does not move other
   fragments, even where their voxel grids overlap.
4. **The whole-model `move_fragment` keeps today's behaviour exactly** (it
   is `move_fragment_named` over every fragment, as one step) so existing
   plans, logs and golden examples are untouched.

The audit log gains one step form, `move:distal:…`, carrying the fragment
name, so a submitted plan record states *what was moved*, not merely *that
something moved* — which is the difference a reviewer cares about.

## Reference-level explanation

Internally the executor tracks a `Vec<NamedFragment>` rather than a single
`VoxelModel`:

```rust
struct NamedFragment {
    name: String,
    model: VoxelModel,
}
```

- `PlanStep::Cut` runs against the *union* of fragments (the composition is
  voxel-exact today only because all fragments share the base grid; see
  Drawbacks), producing one renamed kept fragment and leaving the remainder
  under its previous names.
- `PlanStep::Move(FragmentTransform)` applies to all fragments (today's
  behaviour). `PlanStep::MoveNamed { fragment, transform }` applies to the
  one fragment.
- The final operated model is re-composed for the existing
  `execute()` return value; composition is a scatter into a common bounding
  grid, identical to what `FragmentTransform::apply_to_model` already does.

`PlanStep` remains a non-exhaustive-friendly enum: adding a `MoveNamed`
variant is semver-minor for this crate because the crate is
pre-1.0 and the variant list is documented as exhaustively matchable only
within the workspace (the README states the convention).

## Drawbacks

- **Grid unification.** Fragments moved apart no longer share one grid; the
  re-composition scatter is nearest-neighbour by voxel centre, so two
  fragments reassembled with a sub-voxel gap can alias. The measurement
  report's alignment-error field is the honest statement of this
  quantisation; a caller needing exactness plans on the voxel grid.
- **Contact semantics stay undefined.** Nothing in v0 prevents two moved
  fragments from overlapping voxel-wise. A collision check is a separate
  decision (it needs a tolerance policy) and is explicitly out of scope.
- **Cut-on-moved-fragment.** Cutting a fragment that has already moved
  requires transforming the cutting plane into the fragment's moved frame —
  well-defined but easy to get wrong, so v0 *rejects* cuts after a
  per-fragment move on the same fragment rather than risk a silently wrong
  osteotomy (`PlanError::CutAfterMove`).

## Rationale and alternatives

- *Alternative: tag voxels with fragment ids in the value array.* Rejected:
  it overloads the voxel value (which carries HU/labels that downstream
  crates consume) and forces every consumer to know the tag scheme.
- *Alternative: separate models per fragment with explicit intersection
  ops.* More general, but it changes the executor's contract and every
  caller for a capability one enum variant delivers.
- *Alternative: leave it to callers.* That is today's answer, and it is why
  the audit log cannot describe a two-fragment operation.

## Prior art

Surgical-planning systems (commercial and academic) universally model
osteotomies as named rigid bodies with a scene graph; the naming here
mirrors the crate's existing `fragment_name` so the API stays
plan-as-data rather than growing a scene-graph abstraction.

## Implementation review (2026-10-01): a contradiction found

Preparing this RFC for implementation surfaced a genuine conflict between
its own rules. **Rule 2** (a later cut *splits* every fragment it
intersects: the kept side takes the cut's name, the remainder keeps its
prior name) means the remainder of a *named* fragment is **retained**.
But the executor's cut semantics — unchanged since v0 — **resect** the
discarded side outright (`f64::NAN`, gone from the model). For any
existing multi-cut plan, rule 2 would retain bone that today's executor
removes, violating **rule 4** (existing plans unchanged). A concrete
case: `[cut "distal" (keep z≤0), cut "flake" (keep x≥0)]` — today the
second cut resects its negative side of everything; under rule 2 the
portion belonging to named fragment `distal` would survive.

Resolution options, for review:

1. **Discard-for-the-unnamed-base only.** Cuts split named fragments
   (rule 2 as drafted); the remainder is resected only when the
   intersected fragment is the *unnamed base* — which reproduces today's
   semantics exactly for every existing plan (their cuts all run on the
   base, since nothing is named until a cut names its kept side). Cost:
   the meaning of a cut depends on whether it hits named material, which
   the audit log must then record.
2. **Cuts always resect; fragment identity tracks kept sides only.**
   Rule 2 is dropped: each cut's remainder is resected regardless of
   names, and only kept sides become fragments. Preserves rule 4
   trivially, but a plan can never hold two fragments simultaneously
   unless a single cut's *kept side itself* is later split — which
   discards material each time. That makes per-fragment addressing
   nearly useless for real osteotomy planning, where both pieces are
   kept and moved.

Option 1 is the working recommendation: it is the only reading under
which per-fragment addressing is useful *and* existing plans are
untouched. Implementation waits on this being accepted.

## Unresolved questions

- Should fragment names form a hierarchy (a cut splitting `distal` produces
  `distal/prox`)? v0 says no: flat names, and splitting keeps the parent
  name.
- Does the FDA export need a per-fragment final-pose table? Likely yes for
  real submissions; deferred until the measurement report stabilises.
- The rule-2/rule-4 contradiction above needs a maintainer decision
  (option 1 recommended) before implementation starts.
