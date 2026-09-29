# Changelog

All notable changes to `tpt-med-biomechanics` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- **Per-DOF constraints**: `BoundaryConditions::constrained_dofs` /
  `constrain_dofs(nodes, [x, y, z])` — symmetry planes and roller supports
  constrain individual axes instead of all-3-DOFs `fix_nodes`. Verified by
  an exact symmetry-plane uniaxial state (sigma_zz uniform to 1e-9,
  lateral stress zero).
- **Grid-convergence reporting** (`convergence::convergence_study`,
  `ConvergenceReport`): first-class `CalculationVerification` evidence —
  sorts refinement levels, computes relative errors against an analytic or
  finest-mesh reference, and reports the observed convergence order per
  level pair.
- Crate README documenting the deliberate scope boundary (linear, isotropic,
  static), the mm–N–MPa unit set, and the fidelity ladder to `tpt-fem`.

### Planned
- Nonlinear capability behind a cargo feature via `tpt-fem` /
  `tpt-fem-hyperelastic` / `tpt-fem-contact`, per RFC 0002. The linear core
  stays the default so the WASM footprint is unchanged.
  as a convenience over the current all-3-DOFs `fix_nodes`.
- Direct `tpt-med-tissue` material support, so a single model can mix linear
  bone and hyperelastic soft tissue.
  first-class result rather than an off-script exercise.

### Notes
- `StressResult` and `ElementStress` field layout is public API and is consumed
  by `tpt-med-wasm` and the golden datasets; adding a field is semver-minor,
  changing one is semver-major.
- The hex corner ordering is inherited from `tpt-med-meshing` and is not owned
  here.

## [0.1.0] - 2026-09-22

### Added
- **Element formulation** — 8-node trilinear hexahedra (Q1) with 2×2×2 Gauss
  quadrature: `trilinear_hex_stiffness(&[Vec3; 8], e, nu)` returning the 24×24
  matrix, `isotropic_d(e, nu)` for the 6×6 constitutive matrix,
  `centre_strain_displacement` and `check_element` for validation.
- **Model construction** — `BiomechanicsModel::from_voxel_mesh` (carrying
  per-element HU-derived modulus and Poisson's ratio) and `from_parts` for
  synthetic beam and benchmark meshes.
- **Sparse algebra** — `CsrMatrix::from_triplets` (summing duplicate entries),
  `mul_vec`, `diagonal` for the Jacobi preconditioner, and
  `conjugate_gradient` returning `SolveStats` with the iteration count and final
  residual.
- **Boundary conditions** — `BoundaryConditions` with `fix_nodes` (all 3 DOFs),
  `add_force` (nodal, newtons, duplicates summed) and
  `prescribed_displacements` applied after `fix_nodes`.
- `BiomechanicsModel::solve(&BoundaryConditions, tolerance, max_iterations)`
  → `Result<StressResult, SolverError>`, assembling, solving and
  post-processing.
- **Post-processing** — `StressResult` (displacements, `stats`, per-element
  `stresses`) with `max_von_mises`, `mean_von_mises`, `max_displacement` and
  `critical_element`; `ElementStress` with `von_mises`, `principal_stresses`
  (analytic symmetric eigenvalues via `tpt-med-geometry`, sorted descending)
  and `hydrostatic`. Strains and stresses are Voigt `[xx, yy, zz, xy, yz, xz]`.
- **Model diagnostics** — `rigid_mode_residual()` and `coupling()`, which catch
  the two classic FEM bugs (unconstrained rigid-body motion, and a load that
  engages the structure only weakly) before they produce plausible nonsense.
- `SolverError::Invalid` for empty meshes, material-array length mismatches
  and degenerate elements; no panics, no `NaN`s.

### Verification
Following ASME V&V 40, with verification against closed-form solutions rather
than stored snapshots:
- **Uniaxial tension** on a single hex and a block — stress matches `E·ε`
  exactly and transverse contraction matches `ν`.
- **Cantilever beam** — tip deflection against the Euler–Bernoulli solution
  `δ = FL³/(3EI)`.
- **Patch test** — a uniform-strain patch reproduces constant stress.
- **Rigid-mode residual** is asserted negligible on a properly constrained
  model.
- Degenerate (inverted or zero-Jacobian) elements are asserted to return
  `SolverError::Invalid`.
- Golden datasets `test-data/golden/solid/femur_loading.json` and
  `lumbar_spine_compression.json`, with documented analytical basis.

### Known limitations
- Linear small-strain formulation only: no hyperelasticity, no large
  deformation, no contact, no geometric stiffness. Screening-level
  patient-specific stress analysis is the intended use.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
