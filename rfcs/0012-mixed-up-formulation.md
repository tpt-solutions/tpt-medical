# RFC 0012: Mixed u-p Formulation in the FEM Adapter

- **Status:** Accepted and implemented (2026-10-01). The mixed Q1/P0
  solve, contact coupling, and all five verification items are
  delivered; cutback load continuation remains as an incremental driver
  improvement — see the implementation note at the end.
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

## Implementation attempt (2026-10-01): findings

A working prototype of the Q1/P0 path was built and debugged to the point
where the remaining blocker is precisely identifiable. Findings, all
verified numerically on a single-element uniaxial problem:

1. **The constraint side is settled; the deviatoric side needed a
   correction.** The constraint is the element-mean volume
   `r_p = ∫(J − 1) dV`, and the multiplier's constraint stress in the
   `u` rows is `p·cof(F)` with the **pointwise** cofactor (it
   differentiates the true-`J` term `p·(J−1)`). Two bugs the prototype
   hit and fixed are worth recording: `J̄` is the volume-weighted **mean**
   (forgetting the division by `V` inflates every gradient by
   `(V/J)^{1/3}` — a 100× Jacobian error caught by a
   finite-difference-vs-assembled-Jacobian consistency check), and the
   multiplier's constraint stress must NOT be evaluated at `F̄` (that
   square-scales `K_up` by the same factor).

   **Correction to the deviatoric treatment as first attempted.** The
   prototype evaluated the law at `F̄ = (J̄/J)^{1/3}F` and still measured
   a ~17% stiff nominal stress. Post-mortem: for the workspace's laws —
   whose deviatoric energy already carries the isochoric split
   internally (`C10·J^{−2/3}I₁`-shaped) — evaluating at `F̄` is a
   mathematical **identity**: substituting `det F̄ = J̄` into the model's
   own `J^{−2/3}` factor reproduces `J^{−2/3}(F)`·pointwise exactly, so
   the "modification" changes nothing and the pointwise-dilatation
   locking remains. The mean-dilatation method proper **substitutes the
   mean into the dilatation factor**: the deviatoric Piola is evaluated
   with `J̄` in place of the pointwise `J` in the model's isochoric
   factors — e.g. for the Neo-Hookean branch `P = 2c₁₀·J̄^{−2/3}·F` (the
   constraint stress `p·cof(F)` supplies what the dropped `F⁻ᵀ` term
   would have contributed). This requires a per-law
   mean-dilatation-stress method (the generic trait cannot express
   "replace J inside the law"), which is the actual remaining
   implementation work on the assembly side, alongside the solver
   findings below.

2. **The pure-Lagrange saddle needs a pressure regularization and a
   continuation strategy — this is the remaining blocker.** With
   `K_pp = 0`, the first Newton step from a zero pressure guess has
   `δp ~ δ_load/ε` (enormous for exact-incompressibility-scale ε), and
   residual-norm line searches reject the consistent coupled step, after
   which the iteration stalls short of convergence or drifts onto a
   higher-energy inhomogeneous branch (observed: an hourglass-like mode
   giving a 9× soft or 18% stiff nominal stress depending on the
   deviatoric treatment). The candidate resolutions, in the order they
   should be tried: (a) the **perturbed Lagrangian** with a compliance
   scaled to the material (`J̄ − 1 = ε·p`, ε ~ 10⁻⁵–10⁻⁶·μ) plus
   fine increments; (b) a **Uzawa/augmented** outer loop on the pressure
   (inner penalty-type Newton at fixed p, then a multiplier update);
   (c) an energy-line-search or trust-region Newton on the condensed
   Schur complement. Each needs its own verification pass before the
   closed-form uniaxial test can be asserted at exact-incompressibility
   tolerances.

3. **The Jacobian-consistency check is mandatory tooling.** The prototype's
   finite-difference-vs-assembled comparison caught both assembly bugs in
   seconds; it should ship as a permanent test of whatever implementation
   lands, alongside the closed-form uniaxial verification.

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

## Implementation (2026-10-01): candidate (a) works; the blocker is resolved

The withdrawn prototype's remaining blocker — the pure-Lagrange saddle's
global convergence — was resolved by candidate (a) as ranked: the
perturbed Lagrangian with a material-scaled compliance (`J̄ − 1 = ε̃·p`,
`ε̃ = 1e-8`) regularizes the pressure diagonal enough for plain Newton
from a zero pressure guess, and fine equal load increments (the pressure
field carried between them) handle the moderate-deformation range the
verification covers. The mean-dilatation deviatoric treatment landed as
the corrected form this RFC's findings prescribe: the substituted stress
is the derivative of the substituted energy — it carries **no** `F^{-T}`
term, because the hydrostatic part the pointwise `J`-dependence produced
is exactly what the constraint stress `p·cof(F)` reinstates through the
pressure field; the closed-form uniaxial test asserts the RFC's
hand-derived `p = −μ/λ` directly, confirming the convention. One further
assembly bug of the findings' own class (the `J̄` integral missing the
reference determinant — off by the reference-hex volume, 8×) was caught
by that same closed-form test before shipping. The tangent is the whole
mixed residual differenced, making Jacobian consistency structural
(finding 3's mandatory tooling, sublimated). Verification items 1-4 pass
(items 1 and 2 at exact-incompressibility tolerances; item 3 shows the
full ladder — full integration locks ~18×, SRI over-stiffens ~25%, and
the mixed constraint holds `J̄ = 1` to 1e-8); item 5's contact coupling
is implemented (u-rows only, active set frozen once per Newton
iteration — re-evaluating it inside the line search flips the set
mid-descent and cycles the iteration), but the free-contact grand
cross-validation still cycles engaged/separated under the saddle-point
Newton: the step from the engaged state lands on the wall surface,
deactivating the set, and the step from the separated state dives back
through. Freezing and small increments were verified necessary but not
sufficient; the fix that landed is the
third-ranked candidate, an **activation-tolerance ratchet**: the
pairing's tolerance keeps a node whose step overshoots the wall by the
equilibrium-penetration scale engaged and pulling back, breaking the
cycle. Item 5 then passes in full — the free-contact punch mixed solve
matches SRI's lateral bulge within a few percent at near-incompressible
d1 with the constraint at `J̄ = 1` to 1e-8 — completing all five
verification items. Cutback load continuation remains as an incremental
driver improvement (the linear equal-increment driver shipped with the
mixed solve covers the verification range).
