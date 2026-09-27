# Changelog

All notable changes to `tpt-med-fem-adapter` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Planned
- Friction, either as a substrate feature request or an in-house frictional
  layer on the normal-contact primitives (RFC 0009's second unresolved
  question; `tpt-fem-contact` at 0.1.0 has no friction model).
- A mixed `u`-`p` formulation to replace the penalty volumetric term, which is
  what currently causes volumetric locking on coarse meshes.
- Load stepping / continuation for large-deflection load-controlled paths.
- `Hex20`/`Hex27` and `Tet10` elements, and curved geometry.

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
