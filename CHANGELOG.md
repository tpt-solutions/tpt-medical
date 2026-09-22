# Changelog

All notable changes to `tpt-medical` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning follows
[SemVer](https://semver.org). Crates version and release independently on a
6-week cadence.

## [Unreleased]

### Added — Phase 0: Scaffolding
- Workspace root manifest, `MIT OR Apache-2.0` dual licensing
  (`LICENSE-MIT`, `LICENSE-APACHE`), `cargo-deny` configuration enforcing a
  permissive-only license chain, `rustfmt.toml`, `clippy.toml`.
- Community files: README, CONTRIBUTING (DCO, CLA-free), SECURITY policy,
  Code of Conduct.
- CI: `ci.yml` (fmt/clippy/test/build), `license.yml` (cargo-deny),
  `benchmark.yml`, `docs.yml`, `release.yml`.
- Issue templates (bug report, feature request, RFC) and PR template.
- RFCs 0001–0005 (DICOM ingestion, hyperelastic tissue, FDA audit trail,
  Nitinol superelasticity, cardiac electrophysiology).

### Added — Phase 1: Foundation
- `tpt-med-units` — type-safe unit system for biomedical simulation.
- `tpt-med-geometry` — vectors, 3×3 matrices, planes, AABBs, anatomical
  coordinate transforms (LPS/RAS/ISB).
- `tpt-med-core` — patient models, anatomical region taxonomy, HIPAA-safe
  audit traits.
- `tpt-med-dicom` — pure-Rust DICOM (uncompressed transfer syntaxes) parser,
  CT series assembly, Hounsfield Unit mapping with published HU→density→E
  correlations.
- `tpt-med-meshing` — threshold segmentation, voxel-to-hexahedral meshing,
  Laplacian surface smoothing, CSV export.
- Milestone CLI `dicom-to-mesh`: DICOM CT series → patient-specific bone
  mesh CSV. Synthetic CT generator included; committed test data is
  synthetic only.

### Added — Phase 2: Solid Biomechanics
- `tpt-med-tissue` — Neo-Hookean, Mooney-Rivlin, Yeoh, Ogden strain energy
  functions with analytic first Piola–Kirchhoff stress and verification
  against closed-form uniaxial solutions.
- `tpt-med-bone` — linear elastic bone materials, HU-based property
  assignment (Morgan–Keaveny correlations), Wolff's-law density remodeling.
- `tpt-med-biomechanics` — trilinear hexahedral FEM core (2×2×2 Gauss
  quadrature, CSR assembly, conjugate-gradient solver), nodal/pressure BCs,
  von Mises post-processing; verified against uniaxial tension and cantilever
  beam analytical solutions.
- Milestone: `femur-stress-analysis` example (physiological stance loading).

### Added — Phase 3: Advanced Tissue Models
- `tpt-med-tissue` — Holzapfel–Gasser–Ogden model with fiber dispersion.
- `tpt-med-viscoelastic` — Prony-series relaxation, generalized Maxwell,
  frequency-domain storage/loss moduli.
- `tpt-med-cartilage` — linear biphasic (Mow) confined-compression model.
- Milestone: `arterial-wall-inflation` example.

### Added — Phase 4: Hemodynamics
- `tpt-med-hemodynamics` — incompressible Navier–Stokes on voxel domains
  (staggered MAC projection), Newtonian / Carreau–Yasuda / Casson blood
  models, wall shear stress and oscillatory shear index; Poiseuille
  verification test.
- `tpt-med-cardiovascular` — 2- and 3-element Windkessel, FFR calculator,
  cardiac flow waveforms.
- Milestone: `carotid-bifurcation-wss` example.

### Added — Phase 5: Implants & Devices
- `tpt-med-stents` — Nitinol superelasticity (simplified Lagoudas),
  crimp/expansion deployment, radial force, recoil and dogboning metrics.
- `tpt-med-orthopedics` — implant–bone micromotion, osseointegration risk,
  stress-shielding index.
- `tpt-med-wear` — Archard and Cross–Land wear laws with cycle
  extrapolation to 10⁷ gait cycles.
- Milestone: `stent-deployment` example.

### Added — Phase 6: Surgical Planning
- `tpt-med-surgical-planning` — plane-based osteotomy on voxel models,
  fragment transforms, virtual surgery plans.
- `tpt-med-implant-sizing` — landmark-driven sizing with interpolation over
  manufacturer size charts.
- Milestone: `knee-replacement-planning` example (virtual TKA).

### Added — Phase 7: Regulatory & Compliance
- `tpt-med-audit` — SHA-256 hash chains, HMAC-SHA256 audit signatures,
  tamper-evident verification.
- `tpt-med-fda` — 21 CFR Part 11 audit trails: immutable append-only event
  log, electronic signatures with meaning, FDA-ready export.
- `tpt-med-vv40` — ASME V&V 40 risk-influence-credibility assessment
  matrices.
- Milestone: `fda-package` example (export simulation package with signed
  audit trail).

### Added — Phase 8: WASM & Ecosystem
- `tpt-med-wasm` — `wasm-bindgen` bindings for the imaging→mesh→solve
  pipeline and stent deployment; zero-cloud in-browser execution.
- `web/viewer` — dependency-free WebGL2 surgical planning viewer.
- Milestone: in-browser patient-specific simulation pipeline.
- End-to-end browser glue: shared `web/pkg/` wasm-bindgen output
  (`scripts/build-web.ps1`/`.sh`), viewer wired for DICOM upload →
  `WasmMeshPipeline` → `wasm_solve_stance_load` + `wasm_deploy_stent`,
  WebGL2 rendering of WASM-built meshes, synthetic-CT demo loader.
- `web/stent-simulator` — white-label `<tpt-stent-simulator>` custom
  element (attribute branding, themes, `tpt-deploy`/`tpt-error` events,
  pressure-sweep chart, npm packaging metadata + README).
- `tpt-med-wasm` API: vessel-compliance parameter and `dogboning` output
  on `wasm_deploy_stent`; `element_nodes()` connectivity export on
  `WasmMeshPipeline`; uniform traction distribution fix in
  `wasm_solve_stance_load`.

### Added — Cross-cutting
- Golden reference datasets (`test-data/golden/{solid,fluid,devices,regulatory}`)
  with documented analytical or literature basis.
- Verification tests: uniaxial tension vs analytical Neo-Hookean, Poiseuille
  flow vs analytical solution, stent radial stiffness benchmark scaffold.
- Benchmark suite (`benches/`) and mdBook user guide (`docs/book/`).
