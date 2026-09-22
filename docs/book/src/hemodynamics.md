# Hemodynamics and WSS/OSI

`tpt-med-hemodynamics` solves incompressible Navier-Stokes on voxel domains
with a staggered-grid projection method (SOR pressure solve, Dirichlet
pressure anchor at the outlet, Neumann walls). Blood rheology: Newtonian,
Carreau-Yasuda, Casson.

Post-processing: wall shear stress on wall-adjacent faces and the
oscillatory shear index (OSI) accumulated over a pulsatile cycle.

```rust
use tpt_med_hemodynamics::*;
let domain = FluidDomain::cylinder(24, 12, 4.5, 0.5, 0);
let mut solver = HemodynamicsSolver::new(domain, BloodModel::CARREAU_YASUDA_BLOOD, 30.0, SolverConfig::default());
solver.run_steady(500, 1e-4);
let field = extract_wss(&solver);
```

Verification: mass conservation, paraboloid Poiseuille profile, wall shear
vs the tube law `4μQ/πR³`, shear-thinning resistance ordering. The
`carotid-wss-screening` example runs a stenosed tube under a pulsatile
waveform.
