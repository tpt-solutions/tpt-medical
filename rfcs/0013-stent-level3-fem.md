# RFC 0013: Stent deployment Level 3 — 3D superelastic FEM with frictional contact

- **Status:** Accepted (2026-10-08, maintainer direction); slices 1 (constitutive model), 2 (per-point state, staged driver), 3 (radial wall contact) 4 (crown-ring fixture, ring studies) and 5 (two-group fixture, axial wall profile) implemented
- **Started:** 2026-10-03
- **Crates:** `tpt-med-fem-adapter` (constitutive model, deployment driver),
  `tpt-med-stents` (Level 1/2 cross-check fixtures and metrics)

## Summary

Completes RFC 0004's fidelity ladder: a 3D superelastic finite-element
deployment of a stent — crimp, balloon expansion, and frictional contact
against a vessel — implemented as a Souza–Auricchio constitutive model
behind a `superelastic` cargo feature in `tpt-med-fem-adapter`, driven
through the adapter's existing Newton loop, contact candidates, and
friction state. `tpt-med-stents` stays the metrics crate: Level 3 exists to
be *cross-checked against* the Level 1/2 ring models under RFC 0004's
acceptance items, not to replace them.

## Motivation

RFC 0004 scoped a three-level fidelity ladder and shipped Levels 1–2 (the
1D superelastic ring and tapered ring groups). Level 3 was named and left
open: "3D superelastic FEM with contact … a multi-month effort". What it
serves: radial-force and dogboning claims that must survive contact with a
real (frictional, non-compliant) vessel rather than a Winkler tube law, and
strain fields at the strut crowns for the fatigue path RFC 0004 names as
its follow-up.

**Question of interest:** does a stent design's radial-force/recoil/dogboning
behaviour, computed with 3D contact against a frictional vessel, agree with
the Level 1/2 ring models within their stated envelope — and where it
disagrees, what does the ring model's tube-law assumption cost?

**Model risk:** medium — a device-sizing-adjacent output (radial force feeds
selection), but the Level 1/2 models it is cross-checked against remain the
screening path.

**Model influence:** "supporting" under ASME V&V 40 until RFC 0004's
acceptance items 1–3 are discharged; promotion out of `experimental` is
gated on exactly those items, so the influence claim cannot outrun the
evidence this RFC's verification section produces.

## Detailed design

### Where the code lives

`tpt-med-fem-adapter`, behind an off-by-default `superelastic` cargo
feature, for the same reason the mixed `u`-`p` formulation landed there:
the adapter owns the `Constitutive` trait, the Newton loop, the contact
candidates, and the friction state, and the substrate
(`tpt-fem-hyperelastic` 0.1.0) does not ship 3D constitutive integration
(the gap RFC 0009 recorded). `tpt-med-stents` gains *no dependency* on the
adapter; the cross-check lives as fixtures in the adapter's test suite
comparing against the Level 1/2 formulas (its dev-dependency on
`tpt-med-stents` is the only edge, mirroring how
`tpt-med-biomechanics` was kept out of the adapter).

### Constitutive model: Souza–Auricchio with return mapping

The `Constitutive` trait is stateless per call (`first_piola(&self, f)`),
so the superelastic model carries its internal variables explicitly and
the deployment driver advances them — the same internal-variable pattern
the viscoelastic crate's Simo-type model uses.

- **State per quadrature point:** martensite volume fraction `ξ ∈ [0, 1]`
  and the transformation strain tensor `E_tr` (deviatoric, direction
  fixed at activation).
- **Energy** (Souza–Auricchio form): `Ψ = ξ·Ψ_m(E_e) + (1−ξ)·Ψ_a(E_e) +
  ξ(1−ξ)·Ψ_mix + (β/2)|E_tr|²` with `E_e = E − E_tr`, `Ψ_a/m` Neo-Hookean
  in the respective modulus (the crate's existing incompressible
  Neo-Hookean at the austenite/martensite moduli), and the mixing term
  the standard quartic `β` form.
- **Kinetics:** the thermodynamic drive `X = −∂Ψ/∂E_tr` reduced to its
  equivalent stress; `ξ` evolves on the forward plateau `σ_ms → σ_mf`
  (loading) and reverse `σ_af → σ_as` (unloading) with the same
  cosine-hardening interface the 1D model uses — the 1D and 3D kinetics
  share the plateau bounds so a uniaxial 3D state reproduces the 1D
  stress–strain loop by construction (that is verification item 1).
- **Return mapping:** a closest-point projection on the transformation
  surface per quadrature point per Newton iteration, with the algorithmic
  tangent `∂P/∂F` central-differenced (`material_tangent` already does
  exactly this for stateless laws; the stateful variant differences the
  *returned* update, which is consistent because the driver re-evaluates
  the residual with updated internal variables).

Units: mm / N / MPa / s, the fem-adapter's own conventions. No
`unsafe`: the model is pure Rust over `Mat3`/`Tensor4`.

### Mesh and geometry expectations

The crate does **not** mesh CAD: the caller supplies the stent's Hex8 (or
Tet10) mesh — a laser-cut pattern, or the shipped verification fixture (a
single CrownRing: N struts of Hex8 elements forming one closed ring, with
the crown regions refined). Two fixtures ship in-crate:

- `crown_ring(n_struts)` — the Level-1-comparable geometry: one ring, one
  crown per strut, sized so the Level 1 crown-stiffness abstraction has a
  defined correspondence.
- `two_group_ring()` — the Level-2-comparable geometry for dogboning.

The vessel is an **analytical rigid cylinder** (the contact candidates'
master side takes surface points; a rigid analytical wall avoids meshing
the vessel and matches the Level 1/2 tube-law comparison exactly).
Frictional contact uses `FrictionConfig`/`FrictionState` unchanged.

### Deployment driver

`Deployment3D::new(mesh, model, vessel, options) -> Result<Self, _>` with
three stages, each a load step through the existing Newton loop:

1. **Crimp:** prescribe the outer surface nodes radially onto the crimp
   diameter (displacement-controlled, `loadpath`'s cutback continuation).
2. **Release/expand:** remove the crimp constraints; the superelastic
   recovery drives expansion until contact with the vessel wall engages.
3. **Equilibrate:** contact + friction converge; report.

Metrics are computed from the converged state in
`tpt-med-stents`-compatible definitions (radial force from the contact
reactions — `total_reaction` already sums them; recoil, dogboning from
per-group ring diameters), so the cross-check numbers are the same
*definitions*, not re-derived ones.

### API surface

- `SouzaAuricchio { params } : Constitutive` (+ `params: NitinolParams3D`,
  the 1D crate's defaults extended with the 3D-only constants `β`, `Y`
  (drive threshold scale), documented against the Souza–Auricchio
  literature).
- `Deployment3D` + `Deployment3DOptions { friction: FrictionConfig,
  loadpath: LoadPathOptions, contact: ContactPairing options }`.
- Feature-gated (`superelastic`); zero effect on the default build.

### Error handling

Typed errors following the adapter's own: `SuperelasticError::{
NonMonotonePlateau, ReturnMappingFailed { drive, xi }, …}` — the 1D
crate's non-monotone-plateau rejection carries over; a failed return
mapping (no admissible `ξ` update) is an error, not a silent clamp, by
the same standard the mixed formulation set.

## Alternatives considered

- **Implement in `tpt-fem-hyperelastic`** (the substrate RFC 0004 named):
  rejected for now — RFC 0009 recorded that the substrate at 0.1.0 has no
  3D constitutive integration, and the adapter is where the equivalent
  gap was already filled once (Newton loop, contact, friction, mixed
  u-p). Promotion into the substrate remains the documented upgrade path
  if `tpt-sci` grows the capability.
- **A refined 1D model instead** (per-crown 3D correction factors):
  rejected — it cannot produce the strut strain fields the fatigue
  follow-up needs, which is Level 3's reason to exist.
- **Body-fitted vessel mesh** instead of the analytical wall: rejected for
  the first slice — it breaks the like-for-like comparison with the Level
  1/2 tube law and doubles the contact cost; noted as the natural second
  slice.

## Drawbacks

- A stateful `Constitutive` model strains the trait's stateless shape:
  the deployment driver must own and advance per-quadrature internal
  variables outside the trait, which is more caller machinery than the
  stateless laws need (the viscoelastic crate carries the same cost).
- Souza–Auricchio's return mapping is sensitive to the drive-threshold
  parameterisation; the 1D-calibrated plateau bounds constrain it but do
  not fully determine the 3D constants (`β`, `Y`), which the verification
  must fit and document rather than assume.
- The feature adds a compile-time surface to the adapter that the
  default build never exercises (the mixed formulation already pays this
  cost; CI gains a feature pass the way the dicom codecs did).

## Verification strategy

- **Code verification:**
  1. *Uniaxial loop identity:* a single-element uniaxial state must
     reproduce the `tpt-med-stents` 1D superelastic loop (loading plateau,
     hysteresis, loop closure) within tolerance — the 3D kinetics sharing
     the 1D plateau bounds makes this an identity, not an agreement.
  2. *Patch/constant-stress identity* through the existing adapter tests
     (already shipped for the stateless laws; re-run for the stateful
     model at fixed internal state).
  3. *Return-mapping spot checks:* drive/ξ pairs against hand-computed
     plateau fractions.
- **Calculation verification:**
  4. *Mesh refinement* on `crown_ring`: radial force convergence with
     reported observed order (no analytic reference exists for the ring —
     the refinement trend is the evidence).
  5. *Level 1/2 cross-check* (RFC 0004 acceptance item 1): radial force
     within 15% of `simulate_deployment` for the straight-tube case,
     documented with the parameter set recorded per the RFC 0004 audit
     requirement.
- **Validation:**
  6. *F2394-style radial-stiffness curve* on a published benchmark stent
     geometry (acceptance item 2) — the validation datum the acceptance
     gate names.
  7. *Friction sensitivity study* as a `tpt-med-vv40` assessment
     (acceptance item 3).
- **Golden data:** the crown-ring radial-force series becomes
  `test-data/golden/stents/level3_crown_ring.json` under
  `scripts/diff-golden.sh`.
- **Unverified, stated plainly:** vessel compliance (the analytical wall
  is rigid), stent-vessel migration/foreshortening dynamics, and fatigue
  (RFC 0004's named follow-up). The `experimental` gate stays on until
  all three acceptance items pass.

## Verification results (2026-10-08, crown ring, 6 crowns, 0.2 x 0.15 mm strut)

Fixture: `fixtures::crown_ring` (swept sinusoidal centreline, `R = 4 mm`,
`A = 1 mm`); wall: `RadialWall` at `k = 1e5 N/mm`; material: default
`SuperelasticParams`. Studies live in `ring_studies.rs` (release-profile,
`#[ignore]`d, command in the module docs).

- **Code verification.** Rigid-mode pins carry no load: pin reactions stay
  below `1e-9 N` against forces of `0.3-6 N` (and `< 2 %` on the tiny CI
  ring). Uniaxial loop identity, return-map spot checks, the wall Jacobian
  and the driver's cutback/release behaviour are covered in the earlier
  slices' suites.
- **Calculation verification (item 4).** *Section mesh is irrelevant:*
  1x2, 2x2, 2x3 elements across the strut move the 0.5 % radial stiffness by
  0.65 % (44.91 / 44.66 / 44.62 N/mm at 8 stations per crown). *The centreline
  resolution is the whole error:* the faceted polyline is a stiff polygon.
  Stiffness at 0.5 % compression against stations per crown `m`:
  `m = 12, 16, 24, 32` -> `19.01, 12.13, 7.745, 6.279 N/mm` (per mm of
  diameter). Observed order `2.13` on `(16, 24, 32)` and `2.34` on
  `(12, 16, 24)` — second order, as chord geometry predicts — with Richardson
  limits `4.54` and `4.98 N/mm`. **The finest mesh (`m = 32`) is still about
  38 % above the limit, and the cross-check mesh (`m = 24`) about 60 %: the
  absolute radial force of this fixture is not converged and must not be
  quoted as a prediction.** Hex8 with a faceted centreline needs either more
  stations or curved elements (Hex20) — the answer to the RFC's Tet10-vs-Hex8
  question is "neither is the issue; the centreline is".
- **Level 1 cross-check (item 5), what it does and does not show.** Level 1's
  `crown_stiffness` is an *input*, so the check calibrates it from the 3-D
  ring's own 0.5 % stiffness and tests Level 1's **linearity** — the ratio
  `F_3D / F_L1` is insensitive to the overall stiffness scale, which is why it
  can be read even though the absolute force is not converged. Ratio at
  1-4 % of diameter: `0.983 / 0.965 / 0.952 / 0.940` (`m = 16`) and
  `0.982 / 0.964 / 0.950 / 0.938` (`m = 24`) — mesh-insensitive and inside
  15 %. Beyond that the agreement depends on when the strut transforms, which
  moves outward as the mesh softens toward convergence: the 15 % window ends
  at **6 %** of diameter for `m = 16` and **8 %** for `m = 24`
  (`xi_max` 0.07 / 0.08 at the last point), so 6-8 % is a *lower bound* on the
  converged window. **Acceptance item 1 is therefore discharged only in this
  sense (linearity within 15 % over >= 6 % of diameter); it is not a
  validation of Level 1's absolute crown stiffness**, which needs an
  independent derivation (beam theory or a converged 3-D value) and is open.
- **A recorded false start.** A first ring with 4 stations per crown showed a
  superelastic *plateau* from 3 % compression and a 1 % agreement window. Both
  were locking artefacts of the stiff polygon (struts over-strained, early
  transformation, ring ~50x too stiff). The plateau is real material
  behaviour but only appears in this fixture when the struts are strained past
  transformation — on a resolved ring that is beyond ~6 % compression.
- **Level 2 cross-check: two ring groups, dogboning (2026-10-08).**
  Fixture: `fixtures::two_group_ring` (two unconnected crown rings, axial
  spacing `2A + 1 mm`, each pinned independently — Level 2's groups have no
  axial coupling, and an unconnected stack is the 3-D object it describes) in
  a `RadialWall` with an axial radius profile (new `with_profile`; the slope
  enters the contact gradient and Jacobian, finite-difference verified).
  Study: 6 crowns, 12 stations per crown (so, per the refinement results
  above, the absolute force is not converged — only ratios are read),
  crown stiffness calibrated at a small uniform compression
  (`3.168 N/mm/crown`), then a stepped lumen compressing ring 0 by 1 % and
  ring 1 by 3 % of diameter (inside the elastic window).
  *Per-group radial force, 3-D / Level 2:* ring 0 `1.532 / 1.559 N` (`0.983`),
  ring 1 `4.462 / 4.676 N` (`0.954`) — inside 15 % and consistent with the
  Level 1 linearity result. *Dogboning:* 3-D `0.02000`, Level 2 `0.01898`,
  lumen step `0.02000`. The 3-D value is the imposed step (a rigid wall
  leaves the stent at the wall diameter); Level 2's `5 %` smaller figure is
  its acute-recoil heuristic (the more-compressed group springs back more),
  which **a rigid wall cannot exercise** — testing that term needs a
  compliant vessel and is open. *Superposition:* ring 1 in the stack carries
  `4.4620 N` against `4.4620 N` alone (`0.00 %`), confirming the stack is two
  independent rings and the wall's profile is transparent inside each ring's
  axial band. Acceptance-wise this extends item 1's linearity result to the
  grouped model; it is not a validation of Level 2's recoil law.
- **Still open:** acceptance item 2 (F2394-style curve on a benchmark
  geometry — needs a published stent geometry, not this sinusoid) and item 3
  (friction sensitivity as a V&V 40 assessment — the wall's friction exists
  and is verified at the unit level, the study does not).

## Unresolved questions

- The 3D-only constants (`β`, `Y`): fitted against what uniaxial+
  torsion dataset, and is the fit documented in `tpt-med-vv40` or the
  parameter-hygiene section of RFC 0004? (Settled at implementation time
  by whoever picks the slice up; the RFC requires the fit be cited, not
  assumed.)
- Tet10 vs Hex8 for the crown regions (the adapter supports both) — a
  numerical-methods call to be settled by the mesh-refinement study in
  verification item 4.
- Whether the fatigue follow-up RFC (Goodman-style criteria on Level 3
  strain fields) should be authored alongside, so the strain output shape
  is chosen once rather than re-shaped later.
