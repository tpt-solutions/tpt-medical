# tpt-med-fem-adapter

3-D nonlinear hyperelastic FEM assembly and unilateral contact, built from the
TPT substrate's `tpt-fem` primitives: trilinear hexahedra, a `B^T A B` tangent,
a damped Newton solve, and frictionless contact whose active set is
re-evaluated at every iteration.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--fem--adapter-orange)](https://crates.io/crates/tpt-med-fem-adapter)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--fem--adapter-blue)](https://docs.rs/tpt-med-fem-adapter)

| | |
|---|---|
| **Layer** | `solid` (grouped with `tpt-med-tissue`, the constitutive models it assembles) |
| **Status** | Alpha, `0.1.0` |
| **Scope** | RFC 0009 — the 3-D assembly and contact coupling the substrate does not ship |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.84 |
| **Dependencies** | [`tpt-med-tissue`](../../solid/tpt-med-tissue), `tpt-fem-element`, `tpt-fem-quadrature`, `tpt-fem-sparse`, `tpt-fem-solve`, `tpt-fem-contact` |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

`rfcs/0002-hyperelastic-tissue.md` and `rfcs/0004-nitinol-superelasticity.md`
both said nonlinear hyperelasticity and contact would "land as an adapter
crate", without saying how. `rfcs/0009-nonlinear-fem-substrate-adapter.md` then
read the pinned substrate's actual 0.1.0 source and found the real gap: the
substrate ships `Hex8` shape functions, a hex mesh, hyperelastic *stress
functions*, a generic Newton driver and a COO/CSR assembly — but **no 3-D
hyperelastic assembly**, which is exactly what a 3-D stent or
joint-replacement simulation needs first.

This crate is that assembly, plus the contact coupling RFC 0009 scoped but left
as an open question.

## Features

- **Total-Lagrangian `Hex8` internal force** — `f = int B^T P dV` with `P` from
  `tpt-med-tissue`'s `TissueModel` (every variant, including the Ogden and HGO
  models whose own `P` is a finite difference). Written in index form, with no
  Voigt matrix, so there is no engineering-shear weighting to get wrong.
- **Tangent stiffness `B^T A B`** as the true Hessian of the discrete energy
  (`d2E/du2`), assembled as a `tpt-fem-sparse::Coo`. The material tangent
  `A = dP/dF` is a central difference, so every model — including ones with no
  analytic tangent — is supported without special cases.
- **A second, independent tangent** (`tangent_stiffness_numerical`) that
  differentiates the whole residual, so RFC 0009's "analytic vs. numerical
  tangent" open question is answered by a checked comparison rather than an
  assumption.
- **Damped Newton** with Dirichlet condensation, diagonal equilibration before
  the linear solve, and a convergence measure on the **free**-DOF residual
  (`solve_static`, `residual`).
- **Unilateral contact** (`ContactPairing`) layered onto the nonlinear solve:
  `tpt-fem-contact`'s node pairing and penalty, with the active set recomputed
  from the current geometry at every residual and Jacobian evaluation, so a
  body that separates from the obstacle stops being pushed by it.
- **Regularized Coulomb friction** (`friction` module) on top of that contact,
  since the substrate has none: `f = min(mu * f_n, k_t |s|) s/|s|` against
  tangential slip, with the bound from the same normal penalty. Stateless, so it
  stays correct under the re-evaluated active set, and it ships its exact
  tangent so Newton keeps converging quadratically. See the module docs for why
  regularization rather than a return map, and how to size `k_t`.
- **Explicit failure** where a numerical library might return a plausible
  wrong answer: an inverted element is a `MeshError::InvertedDeformation`, not
  a `NaN` that later surfaces as an unrelated "singular matrix".

## Conventions

- **Units**: lengths in mm, stresses in MPa (the workspace convention). A
  stiffness is therefore MPa and a contact penalty is a stiffness in the same
  units.
- **Volumetric penalty direction**: the in-house convention is
  `W_vol = (J - 1)^2 / d1`, so a **small `d1` is a stiff** penalty. This trips
  up anyone arriving from the usual `1/D1` reading; the verification fixtures
  use `d1 = 0.5` and the crate's own docs say why.
- **Connectivity order**: elements follow `tpt-fem-element`'s `Hex8` reference
  node ordering. `hex_box` emits it; a hand-built mesh must too.
- **Quadrature**: `2` is `2x2x2` (the in-house voxel core's order), `3` is the
  default because a finite-difference material tangent is not a low-order
  polynomial in `xi`.
- **Contact normal**: resolved along a single coordinate axis (`0`/`1`/`2`),
  which is the form `tpt-fem-contact`'s `ContactConstraint` takes.

## Usage

```rust
use tpt_med_fem_adapter::{hex_box, internal_force, solve_static, AssemblyOptions, SolveOptions};
use tpt_med_tissue::{NeoHookeanParams, TissueModel};

let l = 10.0;
let lam = 1.3;
let mesh = hex_box(4, 4, 4, l, l, l)?;
let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.49, d1: 0.5 });

// Fixed base, symmetry plane at z = 0, prescribed axial stretch on the top.
let mut dirichlet: Vec<(usize, f64)> = Vec::new();
for n in mesh.face_nodes(1, false) {
    for c in 0..3 {
        dirichlet.push((mesh.dof(n, c), 0.0));
    }
}
for n in mesh.face_nodes(2, false) {
    dirichlet.push((mesh.dof(n, 2), 0.0));
}
for n in mesh.face_nodes(1, true) {
    dirichlet.push((mesh.dof(n, 1), (lam - 1.0) * l));
}

let opts = SolveOptions {
    assembly: AssemblyOptions::with_quadrature_order(3),
    ..SolveOptions::default()
};
let result = solve_static(
    &mesh, &model, &vec![0.0; mesh.dof_count()], &dirichlet, &opts, None,
)?;

// Reaction on the prescribed top face, reduced to a nominal stress.
let internal = internal_force(&mesh, &model, &result.displacement, &opts.assembly)?;
let reaction: f64 = mesh
    .face_nodes(1, true)
    .iter()
    .map(|&n| internal[mesh.dof(n, 1)])
    .sum();
let nominal = reaction / (l * l);
# Ok::<(), Box<dyn std::error::Error>>(())
```

## API Overview

| Item | Purpose |
|---|---|
| `hex_box(nx, ny, nz, lx, ly, lz)` | Structured `Hex8` box, connectivity in `Hex8` reference order |
| `Hex8Mesh` | Node/element storage, isoparametric Jacobian, `physical_gradients`, `inverted_elements`, `face_nodes` |
| `Constitutive` / `FnModel` | `P = dW/dF`; implemented for `TissueModel` and for any `Fn(&Mat3) -> Mat3` |
| `internal_force(mesh, model, u, opts)` | Global internal force `int B^T P dV` |
| `tangent_stiffness(mesh, model, u, opts)` | Global tangent `B^T A B` as a `Coo` |
| `tangent_stiffness_numerical(...)` | The same, by differencing the whole residual (verification reference) |
| `material_tangent(model, f, h)` | `A = dP/dF` by central differences |
| `solve_static(mesh, model, load, dirichlet, opts, contact)` | Newton solve; returns `SolveResult` |
| `residual(mesh, model, load, u, opts, contact)` | The residual at a configuration, for independent checking |
| `ContactPairing::new(axis, slave, master)` | Unilateral pairing of deforming nodes against a rigid obstacle |
| `ContactConfig { pairing, penalty, friction }` | Which contact to enforce in a solve, how stiff, and with what friction |
| `friction_terms(mesh, pairing, u, penalty, cfg)` | Friction force and tangent for the active contacts, plus per-node stick/slip state |
| `FrictionConfig::new(mu, tangential_stiffness)` | Coulomb coefficient and the regularizing tangential stiffness |
| `solve_load_path(mesh, model, load, dirichlet, opts, contact, path)` | Walks a load-controlled path, one converged increment at a time, with bisection cutback |
| `LoadPathOptions { steps, max_cutbacks }` | Increment count and how deep a failed step may be bisected |
| `hex_box(nx, ny, nz, lx, ly, lz)` | A structured trilinear box; the `Hex8` specialisation of `hex_box_of` |
| `hex_box_of::<E>(...)` | The same, for any hexahedral element — `hex_box_of::<Hex20>` for a quadratic box |
| `tet_box_of::<E>(...)` | A structured tetrahedral box, 6 tets per cell — `tet_box_of::<Tet10>` for quadratic tets |
| `ElementFamily::quadrature_rule(order)` | The rule for an element's reference domain; what the assembly integrates with |
| `Mesh::<E>::default_quadrature_order()` | 2 for a linear hex, 3 for a quadratic one — prefer this over a hard-coded order |
| `LoadStep` / `LoadPath` | The converged points, each with its load factor, residual, iteration count and contact summary |
| `SolveOptions { convergence, assembly }` | Tolerances, iteration cap, quadrature order, FD step |

## Verification

ASME V&V 40 **code verification** (the right answer is known independently of
this code) and **calculation verification** (the error is shown to fall under
refinement). The obligations RFC 0009 recorded as a placeholder are discharged
in `src/tests.rs`.

| Check | What it pins |
|---|---|
| `uniaxial_tension_matches_closed_form` | Nominal stress against `mu (lambda - lambda^-2)`, the same closed form `tpt-med-tissue` and `tpt-fem-hyperelastic` both reproduce; within 3% at 4x4x4 |
| `uniaxial_tension_converges_under_mesh_refinement` | Monotone error decrease over 1x1x1 .. 4x4x4 (measured ratios 1.233, 1.058, 1.019, 1.003) |
| `uniform_deformation_satisfies_patch_identity` | The constant-stress patch identity `f[k,I] = V sum_L G[I,L] P[k,L]`, exact to 1e-9 |
| `uniform_dilatation_matches_the_analytic_volumetric_branch` | The in-house analytic volumetric branch `2J(J-1)/d1 J^-1/3`, exact to 1e-9 |
| `zero_deformation_gives_zero_internal_force` | `F = I` gives an exactly zero internal force |
| `analytic_tangent_matches_numerical_tangent` | The cheap tangent against a differenced one, on a deformed configuration (RFC 0009's open tangent question) |
| `tangent_is_minor_symmetric` | The tangent is the Hessian of a scalar energy, hence symmetric |
| `material_tangent_is_step_size_independent` | `dP/dF` index order and step scaling |
| `contact_jacobian_matches_the_finite_differenced_residual` | The contact-augmented Jacobian is the derivative of the contact-augmented residual |
| `contact_active_set_follows_the_moving_geometry` | The active set is a function of the current geometry (RFC 0009's open contact question) |
| `contact_holds_the_body_out_of_the_obstacle` | A punched block stops at the wall and the reaction balances the punch force |
| `contact_changes_the_answer_versus_no_contact` | The contact terms reach the residual, by contrasting with the unconstrained solve |
| `pulling_away_leaves_the_active_set_empty` | Separation releases the constraint within the same solve |
| `an_inverted_element_is_reported_not_silently_assembled` | Inversion is an error, never a silent `NaN` |
| `friction_force_obeys_the_coulomb_bound` | `f = min(mu f_n, k_t \|s\|) s/\|s\|` in both regimes, and the force always opposes slip |
| `zero_friction_is_exactly_the_frictionless_case` | `mu = 0` or `k_t = 0` gives exactly zero force and tangent |
| `friction_is_not_applied_to_a_separating_node` | No friction without a normal reaction to bound it |
| `a_stuck_node_still_resists_a_tangential_perturbation` | Zero slip is zero *force* but not zero tangent, so a held node cannot creep |
| `friction_tangent_matches_the_finite_differenced_force` | The friction tangent is the exact derivative of the friction force |
| `friction_changes_the_converged_answer` | End-to-end: a tangentially loaded face drifts less with friction, but is not rigidly pinned |
| `load_path_ends_at_the_same_answer_as_a_single_solve` | Stepping changes how the path is walked, not where it arrives |
| `load_path_starts_at_zero_and_advances_monotonically` | The path starts undeformed, load factor only increases, and every point is converged |
| `load_path_displacement_grows_with_load` | Each increment is actually applied, not silently skipped |
| `a_single_step_recovers_a_plain_static_solve` | `steps: 1` is the degenerate path, not a special case |
| `load_path_rejects_a_mis_sized_load_and_an_empty_path` | A bad load length is an error; `steps: 0` is a valid empty path |
| `a_quadratic_box_has_the_right_volume` | Element volume is exact for Hex8/Hex20/Hex27 alike — catches a uniformly wrong-sized box |
| `a_quadratic_box_spans_the_requested_box` | Corner-to-corner extent, node count, and that no node is orphaned |
| `hex20_reproduces_the_uniaxial_closed_form` | A quadratic element's error is mesh-independent and equals the known penalty deviation |
| `hex20_tangent_is_the_derivative_of_the_hex20_residual` | The cheap and full tangents agree on a non-affine deformation |
| `a_wrong_node_count_is_rejected_per_element_type` | An 8-node element in a Hex20 mesh is a named error, not a read past the end |
| `a_quadratic_element_gets_a_higher_default_quadrature_order` | The order floor follows the element type |
| `a_tet_box_fills_its_volume` | Six Kuhn tets per cell, each a sixth of it, summing exactly |
| `a_quadratic_tet_box_also_fills_its_volume` | Same for `Tet10`, which pins the mid-edge nodes |
| `a_tet_box_is_conforming_and_has_no_orphans` | Every node belongs to an element — an orphan would make the system singular |
| `a_tet_box_shares_its_interface_nodes` | Corners shared across cells, mid-edges shared across the tets using them |
| `tet_box_rejects_a_hexahedral_element` | A hex element asked of `tet_box_of` is a named error |
| `a_tet_box_reproduces_the_uniaxial_closed_form` | A tet's error is mesh-independent and equals the known penalty deviation |
| `a_tet_tangent_is_the_derivative_of_the_tet_residual` | The two tangent strategies agree on a non-affine tet deformation |
| `a_tet_element_uses_a_simplex_rule_not_a_tensor_product_one` | Every selectable tet rule is valid on the reference simplex |
| `the_substrate_keast4_tet_rule_is_defective` | Regression pin on an upstream bug — fails intentionally once fixed |

Not verified: mixed `u`-`p` incompressibility, meshing a curved surface from
image data, friction at RFC 0004 Level 3 study parameter ranges (the friction
checks above are single-fixture mechanism tests, not a sensitivity study),
dynamic or quasi-static inertia, and any clinical or ex-vivo data.
`tpt-fem-sparse`'s dense backend makes this a small-problem tool: the linear
solve is `O(n^3)` in DOFs.

## Known Limitations

- **`Hex8`, `Hex20`, `Hex27`, `Tet4`, `Tet10`** (`mesh::Mesh<E>`): the assembly
  is generic over the reference element *and* its family. `hex_box_of::<Hex20>`
  builds a quadratic box; `tet_box_of::<Tet10>` builds a tetrahedral one by
  Freudenthal subdivision. `ElementFamily` is why that works — a tet lives on the
  simplex and needs a Keast rule, not the cube rule a hex takes, and swapping
  them gives a plausible but badly wrong answer. `Hex8Mesh` is a type alias for
  `Mesh<Hex8>`, so nothing existing changes.
- **Penalty incompressibility, not exact.** The volumetric term lives in the
  tissue model, so a stiff penalty plus full integration *volumetrically locks*
  on a coarse mesh. This crate does not buy exact incompressibility over the
  in-house linear-elastic core; it upgrades the material law and the element
  machinery. The fixture's `d1 = 0.5` is a measured compromise, documented in
  `src/tests.rs`.
  **Selective reduced integration is available and helps a lot**
  (`AssemblyOptions::volumetric_quadrature_order`): it splits the volumetric
  penalty out of the fused first Piola and integrates it on a coarser rule,
  cutting locking error on a stiff-penalty coarse mesh from 0.5593 to 0.1218
  (4.6x). It is **opt-in** and off by default, and it *reduces* locking rather
  than removing it — exact incompressibility still needs a mixed `u`-`p`
  formulation, which is Planned. Note it is a silent no-op for a bare
  `FnModel`, which reports no volumetric part; see `Constitutive::volumetric_piola`.
- **Friction is regularized, not an exact return map.** A genuinely stuck node
  carries `k_t * s` rather than a saturated `mu * f_n`, so the result is
  regularization-length dependent: `k_t` must be large enough that
  `mu * f_n / k_t` is far below the displacement accuracy you care about, and
  large `k_t` stiffens the tangential block and costs conditioning. The
  verification fixture uses one `(mu, k_t)` pair; a friction *sensitivity study*
  (RFC 0004 Level 3) is not in this crate. Contact is still resolved along a
  single axis against a rigid, fixed obstacle.
- **No arc-length control, so no limit points.** `solve_load_path` steps the
  load proportionally with bisection cutback, which walks a path that has no
  limit point and stops at one that does. Passing a limit point — the snap-back
  and the softening branch beyond it — needs an arc-length or
  dynamic-relaxation formulation with a load-factor sign convention, which this
  crate does not have. For displacement control, prescribe the displacement and
  call `solve_static` directly.
- **The Newton driver is this crate's, not the substrate's.** See the module
  docs in `src/solver.rs`: `tpt-fem-solve::newton` tests the full residual
  against an absolute tolerance, which a displacement-controlled problem with a
  non-zero reaction can never satisfy.
- **`tangent_stiffness_numerical` is quadratic in problem size** and exists for
  verification, not for production solves.
- **Inverted elements abort the solve.** A line search makes that rare from a
  sane initial guess, but there is no element-level return mapping.

## Related Crates

- [`tpt-med-tissue`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-tissue) — the constitutive models this crate assembles, and the source of the closed forms used to verify it.
- [`tpt-med-biomechanics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-biomechanics) — the linear-elastic voxel core this crate's nonlinear assembly is the large-deformation counterpart to.
- [`rfcs/0009-nonlinear-fem-substrate-adapter.md`](https://github.com/tpt-solutions/tpt-medical/blob/master/rfcs/0009-nonlinear-fem-substrate-adapter.md) — the gap analysis and architecture this crate implements.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Open work this
crate deliberately leaves: friction, a mixed `u`-`p` formulation to remove the
locking noted above, load stepping, and curved elements.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. Not a
diagnostic medical device. The verification recorded here is code verification
against closed forms and mesh refinement — not validation against patient,
ex-vivo or clinical data, which this crate has none of.
