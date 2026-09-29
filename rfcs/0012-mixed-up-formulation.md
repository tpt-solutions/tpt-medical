# RFC 0012: Mixed u-p Formulation in the FEM Adapter

- **Status:** Draft
- **Started:** 2026-10-01
- **Crates:** `tpt-med-fem-adapter`, `tpt-med-tissue`

## Summary

The fem-adapter's only defence against volumetric locking is selective
reduced integration (SRI), a *penalty* mitigation: the near-incompressible
limit is reached through a large `d1`, and the pressure field stays an
implied quantity. A mixed `u`-`p` formulation carries the pressure as a
**field** and enforces incompressibility exactly. This RFC records what an
implementation requires — the architecture constraints found while
assessing the work, the element-pairing decision the crate's CHANGELOG
deferred, and a verification strategy — so the remaining work is
engineering rather than discovery.

## Motivation

`rfcs/0009` delivered the nonlinear `Hex8` core with SRI and left "a mixed
`u`-`p` formulation to the locking" open, noting (via the deviatoric/
volumetric split now in `tpt-med-tissue`) that it is "add a pressure
unknown rather than redesign the trait". That note under-sells the work:
the pressure unknown cannot be condensed per element, so the solver loop
itself changes. This RFC is the design that falls out of actually trying.

## Architecture constraints (found, not assumed)

1. **Element-level static condensation is unavailable.** For a Lagrange
   multiplier the element block matrix is
   `[K_uu, K_up; K_pu, 0]` — the pressure diagonal is *zero*, so the
   Schur complement `K_uu − K_up K_pp⁻¹ K_pu` does not exist. Pressure
   DOFs must go **global**, which touches:
   - the DOF numbering (`mesh.dof_count()` → node DOFs + one per
     currently-active element),
   - the Dirichlet condensation in `newton_from` (only node DOFs are
     condensable; pressures have no essential BCs),
   - `SolveOptions::convergence` (the residual norm must state whether
     pressure rows are included — they have different units),
   - `solve_load_path` (same loop; seeds and increments gain pressure
     entries).
2. **Both stress evaluations must be pressure-parameterised.** The
   `Constitutive` trait splits `first_piola`/`volumetric_first_piola`; the
   mixed form evaluates the deviatoric part alone and takes the pressure
   from the field (`σ = σ_dev + p·J·F⁻ᵀ`-family, exact form per law).
   The trait needs a `first_piola_pressure_jacobian`-shaped method or a
   new `MixedConstitutive` trait implemented for the tissue models.
3. **The Poisson-like gauge problem does not arise** (unlike the CFD
   pressure): the pressure here has physical units and no gauge freedom,
   but it does have a **null space if any element can expand freely** —
   with at least one Dirichlet-constrained node per load path the mixed
   system is regular; the verification suite must include a
   fully-supported case to prove it.

## Element pairing decision (the deferred call)

Two candidates named in RFC 0009:

- **`Hex8` + element-constant pressure (Q1/P0).** inf-sup is *not*
  satisfied on general meshes, but for the structured voxel hex meshes
  this workspace exclusively generates, the practical locking behaviour
  is well-behaved, and it adds exactly one DOF per element. Batch size:
  the whole adapter's assembly, solver, contact and load-path paths.
- **`Hex20`/`Hex8p` (quadratic displacement, linear pressure).** inf-sup
  stable by construction, but requires curvilateral `Hex20` support,
  which the adapter does not have (tet10/hex20 landing covered other
  needs), roughly doubling the assembly work before the mixed part even
  starts.

**Recommendation: Q1/P0 on the existing `Hex8`.** Justification: the
meshes are structured and uniform (voxel-derived), where Q1/P0's inf-sup
deficiency is mildest; SRI remains available as the fallback; and the
verification strategy below catches the pairing's known failure mode
(spurious pressure modes) directly. If a future unstructured-mesh need
arrives, the pairing decision can be revisited without re-designing the
solver (the pressure space is a parameter of the assembly, not of the
loop).

## Implementation shape

- `Mesh` gains a pressure-DOF map (element index → global pressure slot).
- `AssemblyOptions` gains `mixed_pressure: bool` (off by default; SRI and
  mixed are mutually exclusive — running both is a configuration error).
- Internal force per element: `f_u = ∫ Bᵀ σ_dev J dV + ∫ N_p-mapped...`
  (pressure rows: `r_p = ∫ J N_p dV − V̄_p`, the element-mean incompressibility
  constraint weighted by the pressure shape function).
- Newton: assembled as one matrix; the linear solve is the existing
  `tpt-fem-sparse::Coo` path (size grows by element count — acceptable
  for the voxel-mesh scale this crate targets).
- Contact couples only to displacement DOFs; the contact block keeps its
  current shape with pressure rows/columns zero.

## Verification strategy

1. **Closed-form incompressible uniaxial:** one `Hex8` under imposed
   stretch λ with the pressure unknown active must reproduce
   `σ = μ(λ² − 1/λ)` (incompressible Neo-Hookean) with `P₃₃`-direction
   stress free where traction-free — and, the point of the exercise,
   with the lateral stretch *free* the pressure must equal `−μ/λ`
   (hand-derived) rather than a penalty artefact.
2. **Patch test:** constant-strain state reproduces constant stress and
   constant pressure exactly (Q1/P0 passes on uniform hexes).
3. **Locking benchmark:** a 4×4×4 column at ν = 0.4999 under bending-like
   load: mixed tip displacement must converge monotonically under
   refinement where SRI visibly over-stiffens; the SRI run is recorded
   alongside as the baseline (golden-file the comparison).
4. **Spurious-mode check:** the classic Q1/P0 failure — check the
   pressure field for checkerboard modes on a structured mesh under a
   non-trivial load; on uniform voxel hexes with mean-pressure weighting
   the mode is suppressed, and the test asserts it.
5. **Grand test:** an existing `fem-adapter` contact scenario re-run in
   mixed mode must match the SRI result in the compliant limit (small
   `d1` → penalty ≈ constraint) — cross-validating the two paths.

## Drawbacks

- The solver, load-path, contact and assembly all grow a DOF family —
  this is the "not a redesign but not small" cost, now made concrete.
- Q1/P0's inf-sup failure on *distorted* meshes is inherited; a future
  body-fitted mesh need would force the Hex20 pairing after all.
- Two formulations (SRI and mixed) must stay consistent in their shared
  deviatoric path, or verification comparisons drift.

## Alternatives considered

- **B-bar / mean-dilatation.** Avoids new DOFs, but it is still a penalty
  mitigation (the volumetric energy is evaluated at the element-mean
  dilatation) — it does not deliver exact incompressibility and largely
  duplicates what SRI already buys. Rejected as redundant.
- **Nitsche's method for the constraint.** Avoids the multiplier field
  but adds a penalty-weighted boundary formulation whose consistency and
  conditioning analysis is heavier than the multiplier it replaces.
  Rejected.
