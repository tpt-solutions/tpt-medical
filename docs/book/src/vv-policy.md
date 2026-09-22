# Verification & validation policy

Following ASME V&V 40:

1. **Code verification** — every solver ships analytical tests:
   hyperelastic stresses vs closed forms (deviatoric Cauchy comparison) and
   finite-difference references; FEM exact uniaxial/patch tests; Poiseuille
   profile + WSS tube law; SHA-256/HMAC standard vectors; Windkessel
   exponential decay; superelastic loop closure.
2. **Calculation verification** — golden datasets under
   `test-data/golden/` anchor regression behaviour; PRs changing golden
   values must justify the numerics and bump versions.
3. **Validation** — literature-band screens where bench data are unavailable
   (documented in each golden file); curated experimental validation is the
   roadmap, gated by V&V 40 credibility assessments.

Known, documented numerical behaviours: Q1 hexes with full 2x2x2 integration
shear-lock (cantilever band), stair-step voxel walls over-estimate WSS, and
the voxel CFD projection is screening-grade. All are recorded in the golden
notes rather than hidden.
