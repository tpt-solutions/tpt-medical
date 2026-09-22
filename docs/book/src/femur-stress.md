# Femur stress screening

Phase-2 milestone: `femur-stress-analysis` example — segment + mesh the
synthetic CT, fix the distal face, load the proximal head with a 3x
body-weight stance reaction (ISO 7206-style loading convention), solve with
Jacobi-preconditioned CG, and report peak/mean von Mises plus a yield
screen against cortical 110 MPa.

Solver outputs (`StressResult`): nodal displacements, per-element
stress/strain (Voigt), von Mises, principal stresses (analytic eigenvalues),
and CG statistics. Verified against exact uniaxial, patch, and
cantilever-vs-beam references (see `vv-policy.md`).
