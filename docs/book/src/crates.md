# Crates

See the README's crate table for the full list with status. Layering:

- **core**: `tpt-med-units`, `tpt-med-geometry`, `tpt-med-core` (domain
  types, audit traits), `tpt-med-wasm` (browser bindings).
- **imaging**: `tpt-med-dicom` (parsing + HU mapping), `tpt-med-meshing`
  (segmentation, voxel-to-hex meshing, CSV export).
- **solid**: `tpt-med-biomechanics` (hex FEM), `tpt-med-tissue`,
  `tpt-med-bone`, `tpt-med-viscoelastic`, `tpt-med-cartilage`.
- **fluid**: `tpt-med-hemodynamics` (voxel CFD, WSS/OSI),
  `tpt-med-cardiovascular` (Windkessel, FFR).
- **devices**: `tpt-med-stents`, `tpt-med-orthopedics`, `tpt-med-wear`.
- **surgical**: `tpt-med-surgical-planning`, `tpt-med-implant-sizing`.
- **regulatory**: `tpt-med-audit` (SHA-256/HMAC), `tpt-med-fda`,
  `tpt-med-vv40`.

Every crate is `#![forbid(unsafe_code)]` and `std`-only; the TPT substrate
crates (`tpt-fem`, `tpt-science`, `tpt-engineering`, `tpt-math`) are pinned
in the root manifest as the sanctioned nonlinear/high-fidelity integration
points.
