# Changelog

All notable changes to `tpt-med-fem-adapter` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- **Selective reduced integration of the volumetric penalty** (option C of the
  u-p discussion) — a cheap, opt-in mitigation for volumetric locking.
  `AssemblyOptions::volumetric_quadrature_order: Option<usize>`; `None` (the
  default) is the current, fully-verified behaviour, so the existing
  verification suite is unaffected. Measured on a coarse mesh with a stiff
  penalty (`d1 = 0.1`, 2x2x2 elements): locking error **0.5593 -> 0.1218**,
  a 4.6x reduction. This *reduces* locking; it is **not** exact
  incompressibility, which only a mixed `u`-`p` formulation buys.
  - **This required a trait change, which is the part worth knowing.**
    `Constitutive::first_piola` returns one fused `Mat3`, and the
    deviatoric/volumetric split cannot be recovered from it in general —
    under-integrating "all of `P`" would damage the deviatoric response too.
    `Constitutive` therefore gains `volumetric_piola`, and
    `TissueModel::volumetric_first_piola` in `tpt-med-tissue` supplies it
    (`2J(J-1)/d1 * F^-T`; all five model variants share the same penalty, so it
    is one closed form rather than five).
  - **The same split is the first third of a mixed `u`-`p` formulation**, which
    cannot exist until the split does — the pressure unknown *is* the
    volumetric part. This is the concrete reason the cheaper option was taken
    first: its architectural work is not throwaway.
  - **The tangent is split to match.** `internal_force` and `tangent_stiffness`
    both use the deviatoric/volumetric split, via a `VolumetricOnly` view that
    lets the existing finite-difference `material_tangent` difference the
    volumetric law alone. A force that split with a tangent that did not would
    silently degrade Newton to a first-order crawl.
  - **Known sharp edge, documented and tested:** the trait's default
    `volumetric_piola` returns zero, so a bare `FnModel` closure gets a
    **silent no-op** from this option. That is why it is opt-in;
    `reduced_volumetric_integration_is_a_no_op_for_a_law_without_one` pins the
    behaviour so it cannot drift into a surprise.
  - Four new tests: the split is exact against the closed form, SRI measurably
    reduces locking, the no-op case above, and the split tangent is still the
    derivative of the split residual.

- **Regularized Coulomb friction on contact** (`friction` module) — the
  in-house answer to RFC 0009's second unresolved question. `tpt-fem-contact` at
  0.1.0 has no friction model and this workspace does not own that crate, so
  the layer sits on the normal-contact primitives rather than in the substrate.
  `FrictionConfig { mu, tangential_stiffness }` and `friction_terms`, wired
  into the Newton loop through a new `ContactConfig::friction` field.
  - The law is `f = min(mu * f_n, k_t * |s|) * s / |s|` against the tangential
    slip, the bound `mu * f_n` taken from the same normal penalty the contact
    layer already uses. Regularized rather than a return map, because the solver
    recomputes its active set from the current geometry at every residual and
    Jacobian evaluation: a stateless law is correct wherever that is, whereas a
    return map needs its own stick/slip state carried between iterations. The
    module docs state the cost (a stuck node carries `k_t * s` rather than a
    saturated `mu * f_n`) and how to size `k_t` against it.
  - The exact tangent is supplied in both regimes: `-k_t I` on stick,
    `-mu * f_n * (I/|s| - s s^T / |s|^3)` on slip. A node at zero slip gets
    zero force but a *non-zero* `+k_t` diagonal, so a held node cannot creep
    tangentially in the linear solve.
  - `mu = 0` (or `k_t = 0`) is the frictionless case exactly, not
    approximately, and short-circuits before any geometry work.
  - Opt-in and non-breaking: `friction: None` reproduces the previous behaviour
    bit for bit. `ContactSummary` gained `slipping_nodes: Option<usize>`.
- `ContactCandidate` gained a public `node` field. The substrate's
  `ContactConstraint` names only a DOF, which does not identify the node, and
  the frictional layer needs the node's tangential displacement.
- **`Hex20` / `Hex27` and curved geometry** — the assembly is now generic over
  the reference element, so quadratic elements are a type parameter rather than
  a second implementation.
  - `mesh::HexMesh<E: ReferenceElement>` replaces the concrete `Hex8Mesh`
    struct; **`Hex8Mesh` is now a type alias** for `HexMesh<Hex8>`, so every
    existing signature, caller and test is unchanged. `assembly`,
    `internal_force`, `tangent_stiffness`, `solve_static`, `residual`,
    `ContactPairing` and `friction_terms` are all generic over `E`. New code
    should name `HexMesh<Hex20>` rather than adding a parallel `Hex20Mesh`
    type, which would be a second implementation to keep in step.
  - `hex_box_of::<E>` builds a box of any element type. For a serendipity
    element it **compacts the node list** to the nodes some element actually
    references: on a twice-refined grid the face-centre and body-centre
    positions are referenced by no `Hex20` element, and leaving them in gives
    orphan nodes with an exactly zero stiffness row, which makes the condensed
    system singular. This was a real bug found by the verification suite, not a
    hypothetical.
  - `HexMesh::<E>::default_quadrature_order()` returns 2 for a linear element
    and 3 for a quadratic one. Under-integrating a curved element's Jacobian is
    a silent accuracy loss — the mesh just behaves too stiff — so the floor is
    attached to the element type rather than left to the caller to remember.
  - `physical_gradients` now returns `Vec<[f64; 3]>` rather than a fixed 8-row
    array. `E::NUM_NODES` is an associated const and cannot be used as an array
    length on stable.
  - **New `MeshError::WrongElementNodeCount`**: an element whose connectivity
    length does not match the reference element is rejected by name, rather than
    being read past its own end.
  - Curved *reference geometry* is supported (a `Hex20` can represent a curved
    face). Meshing curved surfaces from a DICOM/NIfTI boundary is **not** in
    this crate — there is no isosurface extractor here.
- **`Tet4` / `Tet10` tetrahedra** — a tetrahedron is a different reference
  domain, not just another `E`, so this needed real work beyond the generic
  refactor.
  - **New `ElementFamily` trait** (`mesh`), with `HexFamily` / `TetFamily`
    markers. `ReferenceElement` alone cannot pick a quadrature rule: a hex lives
    on `[-1, 1]^3` and takes a tensor-product Gauss rule, a tet lives on the
    simplex and takes a Keast rule, and substituting one for the other evaluates
    shape functions *outside* the element. `ElementFamily::quadrature_rule` and
    `natural_quadrature_order` are the single place a rule is chosen; the
    assembly, element volume and inversion check all go through it.
  - **`Mesh<E>` renamed from `HexMesh<E>`**, since it is no longer hex-only.
    `Hex8Mesh` remains an alias for `Mesh<Hex8>`, so no caller changes.
  - **`tet_box_of::<E>`** builds a tet mesh by Freudenthal (Kuhn) subdivision:
    six tets per cell sharing the main diagonal, which is what makes the
    decomposition conforming. Mid-edge nodes for a quadratic element are keyed on
    the sorted node pair, so tets sharing an edge — including across internal
    faces and the body diagonal — share the node.
  - **New `MeshError::UnsupportedTetBox`**: asking for a hexahedral element from
    `tet_box_of` is a named error, not a mesh of nonsense.
  - Two real bugs found by the verification suite, both silent: the Kuhn apex
    was numbered as corner 7 instead of 6, which made two of the six
    tetrahedra coplanar (one lying entirely in an `x = 0` face) and gave them
    zero volume; and the assembly was still calling `hex_rule` directly, so tets
    were integrated with a cube rule and a `Tet4` uniaxial solve came out **42x**
    the closed form without erroring.
  - **`tpt-fem-quadrature` 0.1.0's `Keast4` tet rule is defective** and is not
    used. Its weights sum correctly to the reference volume 1/6 and all its
    coordinates are positive, but several of its eleven points fall *outside* the
    reference tetrahedron (measured maximum barycentric sum 1.2607). A `Tet10`
    solve using it returned 42x the closed form rather than failing. Keast3 is
    used instead, which is sufficient: a `Tet10` map is quadratic, so `det J` is
    cubic and a degree-3 rule integrates it exactly. This is worth reporting
    upstream, and `the_substrate_keast4_tet_rule_is_defective` is a regression pin
    that will fail — intentionally — if a future substrate release fixes it.
  - A `Tet10` box has interior nodes on the face and body diagonals, not only on
    the grid lines. Inherent to subdividing a hexahedron, but it means a `Tet10`
    box is not a drop-in replacement for a `Hex20` box of the same dimensions.
- **Load stepping / continuation** (`loadpath` module) —
  `solve_load_path` walks a load-controlled path from zero to the full load, one
  converged increment at a time, each seeded from the last. `LoadPathOptions`
  gives the increment count and the cutback depth; `LoadPath`/`LoadStep` return
  the converged points with per-point residual, iteration count and contact
  summary.
  - Proportional stepping with **bisection cutback**: a failed increment is
    halved, repeatedly, up to `max_cutbacks`, bisecting only the *remaining*
    distance so no converged point is revisited. Only convergence failures are
    retried; a singular system or an inverted element is returned immediately,
    since a smaller step cannot fix those.
  - **Proportional, not arc length.** A limit point is still not reachable: the
    load factor only increases, so the driver stops at one rather than passing
    it. Getting past one needs an arc-length or dynamic-relaxation formulation
    with a load-factor sign convention — a different API, and named as such in
    the module docs rather than implied.
  - Every returned point is a converged equilibrium. There is no partial or
    best-effort path: on failure the error carries the displacement reached
    through the wrapped `SolveError`, so no point on a returned path is unsafe
    to use.
  - `solve_static` is now a thin wrapper over the shared `newton_from` loop,
    which `solve_load_path` reuses so the stepping and the single solve cannot
    drift apart. `solve_static`'s own signature and behaviour are unchanged.

### Planned
- A mixed `u`-`p` formulation for exact incompressibility — now designed
  in `rfcs/0012-mixed-up-formulation.md` (Draft): global pressure DOFs
  (element-level condensation is impossible — the Lagrange block diagonal
  is zero), the Q1/P0 pairing recommended for structured voxel hexes, and
  a five-part verification strategy including the spurious-mode check.
  Implementation pending acceptance.
## [0.1.0] - 2026-09-27

Initial release: the 3-D `Hex8` nonlinear hyperelastic assembly and unilateral
contact coupling scoped by `rfcs/0009-nonlinear-fem-substrate-adapter.md`,
which the substrate at 0.1.0 does not provide.

### Added
- **`Hex8Mesh` / `hex_box`** — trilinear hexahedral mesh with `tpt-fem-element`'s
  reference node ordering, isoparametric Jacobian, `J^-T` physical gradients,
  element volumes, inverted-element detection and face selection, plus a
  structured box builder. Connectivity and DOF-length errors are checked, never
  assumed.
- **`Constitutive` / `FnModel`** — the `P = dW/dF` interface, implemented for
  `tpt-med-tissue`'s `TissueModel` (all variants, including those whose own
  first Piola is a finite difference) and for any `Fn(&Mat3) -> Mat3`.
- **`internal_force`** — total-Lagrangian `f = int B^T P dV`, written in index
  form with no Voigt matrix.
- **`tangent_stiffness`** — `B^T A B`, the exact Hessian of the discrete energy,
  with `A = dP/dF` by central differences. **`tangent_stiffness_numerical`**
  differentiates the whole residual instead, as the independent reference
  RFC 0009 asked for.
- **`solve_static` / `residual`** — damped Newton with Dirichlet condensation and
  diagonal equilibration, converging on the *free*-DOF residual; plus the
  residual as public API so a converged answer can be checked independently.
- **`ContactPairing` / `ContactConfig`** — frictionless unilateral contact whose
  active set is recomputed from the current geometry at every residual and
  Jacobian evaluation, built on `tpt-fem-contact`'s pairing and penalty.
  `SolveResult` reports the active set, the maximum penetration and the total
  reaction.
- **`MeshError::InvertedDeformation`** — an inverted element is an explicit
  error rather than a `NaN` that later surfaces as an unrelated "singular
  matrix" from the linear solver.
- 17 tests: the uniaxial closed form, mesh refinement, the constant-stress
  patch identity, the analytic volumetric branch, both tangent strategies and
  minor symmetry, the contact Jacobian against a differenced residual, and four
  contact scenarios (active set follows geometry, the body is held, contact
  changes the answer, separation releases the constraint).

### Notes on the substrate
- `tpt-fem-solve::newton` is not used: it tests the *full* residual against an
  absolute tolerance, which a displacement-controlled problem with a non-zero
  reaction can never satisfy. The loop here keeps its structure (condense the
  essential DOFs, solve, update) and fixes the convergence measure. The
  substrate's `Coo`, sparse solve, element shape functions, quadrature rules and
  contact primitives are all used as-is.
- `tpt-fem-contact`'s `contact_pairs` returns a non-negative distance, which
  cannot express penetration; the pairing selects the contact partner and the
  signed normal gap is computed here.
- A body resting exactly on the obstacle counts as *active*, which is what
  keeps the first Newton step's linear system non-singular.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
