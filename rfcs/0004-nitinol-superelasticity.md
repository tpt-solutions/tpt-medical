# RFC 0004: Nitinol Superelasticity & Stent Deployment

- **Status:** Accepted
- **Started:** 2026-09-20
- **Crate:** `tpt-med-stents`

## Summary

Select the constitutive treatment of Nitinol and the deployment-simulation
fidelity ladder, from the shipped 1D superelastic ring model up to 3D
superelastic FEM with frictional contact on the `tpt-fem` substrate.

## Motivation

Stent design questions — radial force, chronic outward force, acute recoil,
dogboning, conformability — are exactly the metrics bench testing
(ASTM F2394, F2079) measures, and simulation must speak the same language.
Nitinol's superelastic plateau makes linear-elastic approaches wrong by
construction, while full 3D superelastic contact FEM is a multi-month
validation project. A fidelity ladder lets screening happen today with
honest, versioned caveats.

## Design

### Level 1 — shipped: 1D superelastic material + ring deployment

- **Material**: simplified Lagoudas-style superelastic model. Austenite
  elasticity → cosine transformation interface between σ_ms…σ_mf (forward)
  and σ_af…σ_as (reverse) → martensite elasticity. Strain partitioning on
  the plateaus is solved by bisection; the reverse branch interpolates
  between the plateau-entry strain and the loop-closure strain σ_af/E_a so
  the hysteresis loop closes at the origin (superelastic recovery).
  Verification: loop closure, plateau ordering, mid-strain hysteresis
  (unloading stress < loading stress).
- **Deployment**: N crown springs with Winkler-type vessel reaction.
  Outputs: radial force, contact pressure, acute recoil (elastic
  fraction-limited), dogboning (0 for the uniform ring; non-zero once
  tapered ring geometry lands). Metrics named after the bench standards so
  correlation is direct.

### Level 2 — tapered/ring groups

Multiple rings with per-ring stiffness: captures dogboning and
edge-flare metrics. Pure extension of Level 1; no new material machinery.

### Level 3 — 3D superelastic FEM with contact

Crimp → balloon expansion → vessel contact with friction, on
`tpt-fem-hyperelastic` (superelastic constitutive integration) +
`tpt-fem-contact` (frictional contact), pinned in the workspace manifest.
Acceptance before promotion out of `experimental`:

1. Level 1/2 agreement within 15% on radial force for a straight tube.
2. ASTM F2394-style radial-stiffness curve reproduction on a published
   benchmark stent geometry.
3. Friction sensitivity study documented as V&V 40 evidence
   (`tpt-med-vv40` assessment shipped with the crate).

## Parameter hygiene

Default parameters (E_a = 55 GPa, E_m = 28 GPa, ε_L = 5%, plateau stresses
480/560/380/260 MPa) are literature-typical for laser-cut stent tubing at
body temperature, cited in code; they are **starting points**, not
vendor data. Deployment studies must record the parameter set in the audit
trail (`tpt-med-fda`) — the golden dataset includes the reference set.

## Alternatives considered

- **Souza-Auricchio model** (3D, smeared): the standard FEM choice; adopted
  implicitly at Level 3 via `tpt-fem-hyperelastic` rather than re-implemented
  here.
- **Viscoelastic vessel wall** in the reaction model: deferred — the rigid /
  linear tube law covers screening; FSI is a Phase-4-scale effort.

## Unresolved questions

- Fatigue (mean/alternating strain per FDA guidance) — needs Level 3 plus a
  Goodman-style criteria module; candidate follow-up RFC.
- Resorbable scaffolds (degradation-coupled stiffness) — out of scope until
  a partner use case exists.
