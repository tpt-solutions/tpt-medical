# RFC 0009: Nonlinear FEM Substrate Adapter

- **Status:** Draft
- **Started:** 2026-09-27
- **Crates:** `tpt-med-tissue`, `tpt-med-biomechanics` (adapter), new
  `tpt-med-fem-adapter` (proposed)

## Summary

Scopes the actual engineering work behind RFC 0002's "Integration lands as
an adapter crate ... behind a cargo feature" line, and RFC 0004 Level 3's
"3D superelastic FEM with contact" — both of which named the destination
without designing the path. Having now read the pinned substrate crates'
real 0.1.0 APIs (`tpt-fem-hyperelastic`, `tpt-fem-mesh`, `tpt-fem-element`,
`tpt-fem-contact`, `tpt-fem-solve`), this RFC records what is and is not
already there, and proposes the adapter architecture for the gap: a 3D
`Hex8` nonlinear hyperelastic assembly, which the substrate does not yet
ship, built from primitives it does ship.

## Motivation

Two already-accepted RFCs point at this integration without designing it:
RFC 0002 says nonlinear hyperelasticity and contact integration "lands as an
adapter crate ... behind a cargo feature," and RFC 0004 Level 3 requires "3D
superelastic FEM with contact" before Nitinol stent deployment can leave
`experimental`. Neither says *how* — reasonably, since at the time neither
had read the substrate's actual 0.1.0 surface closely enough to know what
that "how" would need to cover. This RFC does that reading and reports back.

**Question of interest:** what would a caller of `tpt-med-tissue` or
`tpt-med-biomechanics` actually get from enabling a `substrate-fem` cargo
feature, and what does *not* exist yet that a naive reading of RFC 0002/0004
might suggest already does? **Model risk:** this RFC itself computes and
decides nothing about tissue mechanics — it is an architecture and gap
analysis. Its risk is entirely in what it might get wrong about the
substrate's actual capabilities, which is why every claim below is backed by
a direct read of the pinned crate's source (see citations), not the
substrate's own RFC-level descriptions of itself. **Model influence:** none
directly; it gates whether any future 3D-contact stent or joint-replacement
simulation can honestly claim substrate-grade fidelity, which is exactly RFC
0004 Level 3's own gate.

## Detailed design

### What the substrate already provides (verified by reading the pinned 0.1.0 source)

- **`tpt-fem-element`**: reference-element shape functions and gradients for
  `Hex8` (and `Tet4`, `Quad4`, `Line2`, plus quadratic `Hex20`/`Hex27`/
  `Tet10`). A 3D hex element's geometry and interpolation are fully
  supported.
- **`tpt-fem-mesh`**: `CellType::Hex` meshes, node/element storage, Gmsh
  import. A mesh built from `tpt-med-meshing`'s `VoxelHexMesh` output could
  be re-expressed as a `tpt_fem_mesh::Mesh` (a format-conversion adapter,
  not a numerics one).
- **`tpt-fem-hyperelastic`**: **stress functions only** —
  `neo_hookean_piola`, `mooney_rivlin_piola`, `ogden_piola` compute the first
  Piola-Kirchhoff stress `P = ∂Ψ/∂F` for a single deformation gradient `F`,
  given an externally-supplied incompressibility pressure `p` (a Lagrange
  multiplier the caller must already have solved for, not one the function
  derives). The crate's only assembled solver is `solve_hyperelastic_bar`:
  a **1-D two-node bar** (`Line2`, one axial DOF per node) large-deformation
  Newton solve. There is **no 3D hex (or tet) hyperelastic assembly** in
  this crate at 0.1.0 — a caller wanting one has to build it.
- **`tpt-fem-solve`**: a generic Newton driver (`newton`, taking caller-
  supplied residual and tangent closures) plus arc-length continuation. This
  is genuinely reusable for a 3D assembly — it does not assume 1D.
- **`tpt-fem-sparse`**: `Coo` triplet assembly and a sparse solve, likewise
  generic and reusable.
- **`tpt-fem-contact`**: **DOF-level unilateral constraints** (`x_dof ≥
  lower`) via penalty or augmented-Lagrangian, plus `contact_pairs`
  (Octree-based nearest-node pairing between two point sets). This is a
  *linear* contact-constraint layer operating on a caller-supplied base
  stiffness `Coo` — it has no built-in coupling to a nonlinear (large-
  deformation) residual/tangent, and (from the read above) no explicit
  friction model; only normal non-penetration is implemented.

### The actual gap

**A 3D `Hex8` nonlinear hyperelastic finite-element assembly does not exist
in the substrate and must be written**, using:

1. `tpt-fem-element::Hex8` for shape functions/gradients and
   `tpt-fem-quadrature` for the Gauss rule (2×2×2, matching the existing
   in-house linear hex core's own scheme per RFC 0002's Design section, for
   a consistent locking-behaviour comparison).
2. `tpt-fem-hyperelastic`'s stress functions for the constitutive response
   at each quadrature point — but the incompressibility pressure `p` these
   functions take as a parameter must be **solved for**, not assumed, for
   any 3D state (unlike the 1-D bar case, there is no closed form for `p`
   in general 3D deformation). This needs either a **mixed `u`-`p`
   formulation** (displacement and pressure both primary unknowns, the
   standard approach for incompressible hyperelasticity) or accepting the
   existing in-house crates' **penalty** approach (`(J-1)²/D1`, already
   implemented in `tpt-med-tissue`) instead of true incompressibility. The
   penalty route reuses `tpt-med-tissue::TissueModel::first_piola` directly
   (which already returns a complete `P` including its own volumetric
   term) and needs no Lagrange-multiplier plumbing — **this RFC recommends
   the penalty route for the first 3D adapter increment**, deferring a
   mixed formulation as a later, separate increment once a real use case
   needs true incompressibility at the substrate level.
3. `tpt-fem-solve::newton` for the nonlinear equilibrium iteration, with a
   caller-supplied tangent — meaning **the tangent stiffness (`∂P/∂F`
   contracted through the element's shape-function gradients) must also be
   derived and assembled**, either analytically per model (more code, no
   per-step cost) or via numerical differentiation of the residual (less
   code, slower and noisier — the existing in-house core's own Newton path,
   if any, and `tpt-med-tissue`'s own `first_piola_numerical` precedent
   for Ogden/HGO are the relevant local prior art on this tradeoff).
4. `tpt-fem-sparse::Coo`/`solve` for the linear solve inside each Newton
   iteration.
5. `tpt-fem-contact`'s `ContactConstraint`/`contact_pairs` layered on top of
   the same `Coo` tangent, once the nonlinear assembly above produces one to
   layer it onto — contact is additive to this design, not a separate
   integration effort, but it cannot start before step 1-4 exist.

### Proposed crate boundary

A **new** crate, `tpt-med-fem-adapter` (name pending; alternatives
considered below), rather than extending `tpt-med-tissue` or
`tpt-med-biomechanics` in place:

- It depends on `tpt-med-tissue` (for `TissueModel`/`first_piola`),
  `tpt-med-meshing` (for `VoxelHexMesh` → `tpt_fem_mesh::Mesh` conversion),
  and the substrate crates (`tpt-fem-element`, `tpt-fem-mesh`,
  `tpt-fem-quadrature`, `tpt-fem-solve`, `tpt-fem-sparse`, optionally
  `tpt-fem-contact`) — a dependency direction `tpt-med-tissue` itself should
  not carry, so the default (dependency-free) build of the tissue-model
  crate is unaffected by whether this adapter exists.
- It is gated by a cargo feature only in the sense that *depending on this
  crate at all* is the opt-in — mirroring how `tpt-med-surgical-planning`
  depends on `tpt-med-meshing` directly rather than `tpt-med-meshing`
  gaining a `surgical-planning` feature. This is a **correction** to RFC
  0002's "behind a cargo feature" phrasing: a separate crate is the
  workspace's actual established pattern for optional substrate-dependent
  functionality (see `tpt-med-nifti` needing no feature flag on
  `tpt-med-dicom` to exist), and is preferred here over a feature flag on an
  existing crate for the same reason — it keeps `tpt-med-tissue`'s and
  `tpt-med-biomechanics`' own `Cargo.toml` unchanged and their default
  builds exactly as dependency-light as today.
- Public API sketch:

```rust
/// Builds a tpt-fem mesh from an in-house voxel-hex mesh (format
/// conversion only — no numerics).
pub fn to_fem_mesh(mesh: &tpt_med_meshing::VoxelHexMesh) -> tpt_fem_mesh::Mesh;

/// Solves 3D nonlinear hyperelastic equilibrium under Dirichlet boundary
/// conditions, using tpt-med-tissue's penalty-formulation models and
/// tpt-fem-solve's Newton driver.
pub fn solve_hyperelastic_3d(
    mesh: &tpt_fem_mesh::Mesh,
    material: &tpt_med_tissue::SoftTissueMaterial,
    dirichlet: &[(usize, f64)],
    options: &tpt_fem_solve::NewtonOptions,
) -> Result<Vec<f64>, FemAdapterError>;
```

  Contact (`solve_hyperelastic_3d_with_contact`, taking additional
  `ContactConstraint`s) is deliberately **not** in this first sketch — see
  Unresolved Questions.

### Alternatives considered

- **A cargo feature on `tpt-med-tissue` itself**, as RFC 0002 literally
  said. Rejected on rereading against the workspace's own established
  precedent (above): every other substrate-dependent capability in this
  workspace (NIfTI, meshing, surgical planning) is its own crate, not a
  feature flag on an unrelated crate. RFC 0002's phrasing predates this
  RFC's closer look at the substrate and at the workspace's own pattern;
  this RFC treats that as a correction, not a contradiction, since RFC
  0002's actual intent (opt-in, default build unaffected) is preserved
  either way.
- **A mixed `u`-`p` formulation from the start**, for true incompressibility
  rather than the penalty approximation. Rejected for the first increment:
  it roughly doubles the unknowns per element and needs an inf-sup-stable
  element pairing (e.g. `Hex8`/constant pressure or a Hex20/Hex8 pairing) —
  a real numerical-methods decision this RFC is not resolving pre-emptively
  when the in-house core's existing penalty formulation already has a
  working, verified (RFC 0002) precedent to reuse directly.
- **Numerical (not analytic) tangent stiffness for the first increment.**
  Seriously considered, not rejected outright — see Unresolved Questions.
  It is slower and noisier per Newton step but roughly halves the new code
  this RFC's design requires, since no per-model tangent derivation is
  needed.

### Drawbacks

- This RFC authorizes **no implementation** by itself — it is gap analysis
  and architecture, deliberately, given how large and currently-unscoped
  the actual assembly work turned out to be once the substrate's real 0.1.0
  surface was read closely. A reader expecting "here is the 3D contact
  solver" from RFC 0002/0004's own language will not find one here either.
- The recommended penalty-formulation route (not a mixed `u`-`p` method)
  means the first 3D substrate adapter increment will **not** be more
  physically correct than the existing in-house linear-elastic voxel core
  is for incompressibility — it upgrades the material law (nonlinear
  hyperelastic vs. linear elastic) and the substrate's Hex8 element/solver
  machinery, but does not yet buy exact incompressibility. That gap is
  named, not hidden.
- No friction model exists in `tpt-fem-contact` at 0.1.0 as read; RFC 0004
  Level 3's acceptance criteria include "a friction sensitivity study,"
  which cannot be attempted at all until either the substrate gains one or
  this workspace writes its own — a second, currently entirely unscoped
  gap this RFC surfaces but does not resolve. **Partially resolved (2026-09-27):**
  the in-house frictional layer now exists (`tpt-med-fem-adapter::friction`);
  the sensitivity study itself is still out of scope here.

## Verification strategy

Not applicable in the way the template expects: this RFC proposes no
numerical scheme to verify, only an architecture and a crate boundary for
work a future RFC or the same RFC's own future revision would need to
detail with a real verification strategy (code verification against the
existing in-house core's own uniaxial/patch tests at matching parameters,
since both would claim to solve the same physical problem; calculation
verification via mesh refinement once a real 3D assembly exists). Recorded
here as a placeholder obligation, not skipped: **any implementation PR
against this RFC's design must add that verification strategy before
merge**, matching every other numerical RFC in this project.

## Unresolved questions

- **Analytic vs. numerical tangent stiffness.** Whoever implements this
  should benchmark both on a small patch-test-sized problem before
  committing — this RFC deliberately does not pick, since the tradeoff
  (code volume vs. per-step cost and Newton robustness) is better judged
  against a real, running prototype than in the abstract.
- **Contact coupling design.** `ContactConstraint`/`contact_pairs` operate
  on a linear `Coo` stiffness; wiring them into a nonlinear Newton iteration
  correctly (constraints re-evaluated each iteration as the geometry moves,
  not fixed at the initial configuration) is real, unscoped numerical-
  methods work this RFC explicitly leaves open. Settled by: whoever picks
  up RFC 0004 Level 3's promotion criteria, since that is the actual
  consumer motivating contact at all.
- **Friction.** Needs either a substrate-side feature request (this
  workspace does not own `tpt-fem-contact`) or an in-house frictional layer
  on top of the normal-contact primitives that exist. Not scoped here.
  **Partially resolved (2026-09-27):** the in-house option is now taken, in
  `tpt-med-fem-adapter`'s `friction` module — a regularized Coulomb law rather
  than an exact return map, chosen so it stays stateless under the
  re-evaluated active set this crate already uses. What remains open is
  validation at study parameter ranges (RFC 0004 Level 3's friction
  sensitivity study), not the mechanism.
- **Crate name.** `tpt-med-fem-adapter` is a placeholder; naming follows
  whatever convention feels least awkward once the crate's actual shape is
  clearer from a prototype (compare how `tpt-med-nifti` and
  `tpt-med-meshing` are named for *what* they do, not *how*).
- **Does this RFC need Acceptance before any prototyping starts, or would a
  small non-committal spike (behind no public API, not wired into the
  workspace manifest) be reasonable first?** A spike would inform the two
  Unresolved Questions above with real numbers rather than guesses.
  Settled by: whoever has the bandwidth to run it.
