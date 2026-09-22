# Stent deployment

`tpt-med-stents` ships a 1D superelastic Nitinol model (cosine
transformation interface, loop closure at the origin) and a ring deployment
model producing radial force, contact pressure, acute recoil, and dogboning
— the ASTM F2394/F2079 metric family.

Fidelity ladder (rfcs/0004): Level 1 ring model (shipped) → tapered ring
groups → 3D superelastic FEM with frictional contact on `tpt-fem` (pinned).
Default parameters are literature-typical starting points, not vendor data.

The same model ships in the browser via `tpt-med-wasm::wasm_deploy_stent`,
both in the [viewer](./wasm.md) and as the white-label
`<tpt-stent-simulator>` web component under `web/stent-simulator/`.
