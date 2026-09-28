# Changelog

All notable changes to `tpt-medical` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning follows
[SemVer](https://semver.org). Crates version and release independently on a
6-week cadence.

**This file records workspace-level and cross-cutting changes only.** Each
crate has its own `CHANGELOG.md` recording what changes for its consumers;
those are the source of truth for per-crate behaviour. When a change is
user-visible in one crate, record it in that crate's changelog, and record it
here only if it affects the workspace as a whole.

## [Unreleased]

### Added
- **`tpt-med-fem-adapter`** (new crate, `solid` layer) — the 3-D nonlinear
  hyperelastic `Hex8` assembly and unilateral contact coupling that
  `rfcs/0009-nonlinear-fem-substrate-adapter.md` scoped and the substrate at
  0.1.0 does not ship. Total-Lagrangian internal force, a `B^T A B` tangent
  (settling RFC 0009's open analytic-vs-numerical question with a checked
  comparison), a damped Newton solve over `tpt-fem-sparse`, and frictionless
  contact whose active set is recomputed from the current geometry at every
  iteration. 17 tests: uniaxial closed form, mesh refinement, the exact
  constant-stress patch identity, the exact analytic volumetric branch, both
  tangent strategies, and four contact scenarios. The RFC's placeholder
  verification obligation is discharged; friction, a mixed `u`-`p`
  formulation, load stepping and curved elements remain open.
- **New crate: `tpt-med-nifti`** (`rfcs/0006-nifti-ingestion.md`) — pure-Rust
  parsing of NIfTI-1 volumes, so a research-space CT/MR export
  (`dcm2niix`, FSL, FreeSurfer, ANTs, most public imaging datasets) can
  enter the pipeline without converting back to DICOM first. Single-file
  `.nii` **and** dual-file `.hdr`/`.img` layouts, sform/qform geometry
  (NIfTI's own precedence), 8 datatypes, `scl_slope`/`scl_inter` value
  scaling, and `.nii.gz` behind an off-by-default `gzip` feature (pure-Rust
  `flate2`; streaming, geometry-capped inflate). NIfTI-2 remains out of
  scope — see the crate's own CHANGELOG and RFC 0006's Unresolved Questions
  (all three now annotated with how they were resolved) for why each
  decision was made rather than guessed at. Workspace member count: 24.
- **NIfTI enters the mesh pipeline.** `tpt-med-meshing` gained
  `SegmentationMask::threshold_nifti(&NiftiVolume, min_hu)` — a direct
  constructor, not a shared volume trait, because `NiftiVolume` is still
  the only non-DICOM source (the shape RFC 0006's Unresolved Question told
  the first implementer to pick). RAS → LPS conversion happens at that
  boundary, so masks and meshes from either source format share one
  patient frame.
- **CI: `tpt-med-dicom` compressed-pixel-data feature coverage.** Its `rle`,
  `jpeg`, `jpeg-ls` and `jpeg2000` cargo features are all off by default, so
  `cargo clippy --workspace` and `cargo test --workspace` never compile them.
  The `clippy` and `test` jobs now run a dedicated
  `-p tpt-med-dicom --features rle,jpeg,jpeg-ls,jpeg2000` pass so those code
  paths are actually checked, not just locally. The same pattern covers
  `tpt-med-nifti`'s off-by-default `gzip` feature (a
  `-p tpt-med-nifti --features gzip` pass in both jobs).
- **Per-crate documentation set.** Every one of the 23 workspace members now
  ships a comprehensive `README.md` and a `CHANGELOG.md`:
  - `README.md` — overview, why the crate exists, features, conventions and
    units, runnable usage examples, an API overview table, the verification
    evidence and what it actually proves, known limitations, related crates,
    license, and the regulatory disclaimer.
  - `CHANGELOG.md` — Keep a Changelog format, with `[Unreleased]` (planned work
    and the semver classification of behaviour-changing areas) and
    `[0.1.0] - 2026-09-22`.
- **Per-crate crates.io metadata.** Each crate now declares its own
  `keywords` and `categories` instead of inheriting workspace-wide generics,
  plus an explicit `readme` and (for publishable crates) a
  `documentation` link to docs.rs. All category slugs are validated against
  the live crates.io category registry.
- **CI guard** (`crate-docs` job) enforcing the per-crate documentation
  convention: every workspace member must have a `README.md` and a
  `CHANGELOG.md`, and its `keywords`/`categories` must satisfy the crates.io
  rules (1–5 keywords of ≤ 20 characters each; ≤ 5 categories; categories drawn
  from the published registry). The list of valid categories is vendored into
  the job so the check does not depend on network access.
- Root `README.md` crate table now links each crate's README and CHANGELOG
  directly.
- [`docs/book/src/crates.md`](docs/book/src/crates.md) links the per-crate
  documentation and states the required README section set.

### Changed
- Root `README.md` crate table gained a **Docs** column and an introductory
  paragraph describing the per-crate documentation contract.
- `[workspace.package] keywords` and `categories` in the root manifest are now
  documented as a fallback for any new crate; all 23 existing members override
  them.

### Fixed
- `deny.toml` used four `[licenses]` keys (`copyleft`, `allow-osi-fsf-free`,
  `default`, `unlicensed`) that current `cargo-deny` (0.16+) removed —
  `cargo deny check` failed outright with a config-validation error rather
  than checking anything. Removed; `version = 2` already denies by default
  anything not in the `allow` list, so this is a syntax fix, not a policy
  change. Found while adding `tpt-med-dicom`'s new JPEG-family optional
  dependencies and wanting to verify their licenses against the gate; CI's
  `cargo-deny-action` would have hit the same failure on the next push
  regardless of those dependencies.

### Phase 9: platform review follow-ups
Tracked in `todo.md`. Summary of the user-visible change here:

#### Fixed
- `README.md` Quick Start used `--bone-threshold`, a flag `dicom-to-mesh` does
  not accept. The adjacent Rust snippet also did not compile (missing
  `Path::new`, and `voxels_to_hex_mesh` is a `&self` method returning
  `Result`). Both corrected — the first command and the first snippet a new
  user runs now work.
- `tpt-med-fda`: attaching a reproducibility manifest produced **invalid
  JSON**. The export closed the outer object before appending the manifest, so
  the package was signed but unparseable by any consumer. Found by
  `scripts/diff-golden.sh`, and now covered by a structural single-JSON-object
  test in the crate.
- Removed `crates/imaging/tpt-med-dicom/src/dbg_test.rs`, an orphaned debug
  test not wired into the crate and ending in an unconditional `panic!`.
- Seven per-crate READMEs had their sections out of order, and one had lost a
  chunk of its Features list. `scripts/check-crate-docs.sh` now asserts
  section *order*, not just presence, and all 23 are repaired.
- `scripts/check-crate-docs.sh` reported all ten required sections of five
  crates as "has no content" on a CRLF working tree. The "section has a body"
  check compared the heading with awk's `$0 == sec`, which never matches when
  the line carries a trailing `\r`; the repo pins only `*.sh` to LF, so any
  `core.autocrlf=true` checkout (the default on Windows) hit this. The awk now
  strips the carriage return before matching. The failure was a false positive
  in the section-body check only — the section-presence, order, and code-fence
  checks were already CRLF-safe — and CI on Linux checkouts was unaffected.

#### Added
- **`tpt-med-fda`: reproducibility manifest** — `ReproducibilityManifest`
  records the workspace version, per-crate versions, git commit, build
  profile, and a SHA-256 plus byte length of every input artefact.
  `AuditTrail::attach_manifest` also appends an audit event carrying the
  manifest digest, so a manifest swapped after the fact is detectable from the
  chain alone, and the manifest body is covered by the export's detached
  HMAC. Runs that do not attach one export byte-identically, so packages
  signed before this change still verify. Wired into the `fda-package`
  milestone.
- `scripts/diff-golden.sh` — renders a before/after numeric-drift table for
  `test-data/golden/` against a PR, with each row judged against the
  `tolerance_percent` the golden file itself declares. CI job
  `golden-drift` fails the build on drift beyond tolerance, on added or
  removed fields, and on unparseable JSON.
- `scripts/bench-baseline.sh` + `benches/baseline.txt` — benchmark suite with
  a stored per-stage baseline; a regression beyond the threshold (25 % by
  default) now fails CI instead of being printed and ignored. Regeneration is
  a deliberate manual `workflow_dispatch` act.
- `templates/` — crate, example-binary and RFC templates. Copy-based rather
  than `cargo generate`, because the thing worth templating is the house
  style, and the contract the templates promise is the one CI already
  enforces on real crates.
- `.github/workflows/pages.yml` — publishes `web/viewer/` as a live
  zero-cloud demo on GitHub Pages, built from the WASM engine against the
  committed **synthetic** CT data. Kept separate from `docs.yml` because a
  WASM + wasm-bindgen build is far more expensive than rustdoc.
- `.gitattributes` pinning `*.sh` to LF, so CI shell scripts do not break on
  a CRLF Windows checkout.

#### Governance
- `CONTRIBUTING.md` rewritten: **this project does not accept external pull
  requests.** Work arrives as GitHub issues; changes land from maintainer
  branches. The RFC process now runs through the issue tracker, with a
  maintainer committing the accepted RFC. `.github/PULL_REQUEST_TEMPLATE.md`
  and the RFC/feature issue templates were updated to match, and the root
  README's Contributing section now says so.

### Post-Phase 9 roadmap
Tracked in `todo.md`. Summary of the user-visible change here:

#### Added
- **Per-crate versions in the reproducibility manifest.** `examples/build.rs`
  line-scans the workspace `Cargo.lock` into a generated
  `tpt_med_examples::crate_versions::CRATE_VERSIONS` table (no new
  dependency); the `fda-package` milestone now looks up each participating
  crate's real resolved version there, instead of stamping every crate with
  the calling binary's own `CARGO_PKG_VERSION`. `ReproducibilityManifest`
  already had a `crates` map for this (`with_crates`); the gap was only in how
  the example populated it. Every crate still resolves to the same `0.1.0`
  workspace version today, so this is invisible in output until crates start
  releasing independently — see `tpt-med-fda`'s README "Known Limitations".

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

### Added — Phase 8 extension: business wedge
- `web/stent-sim/` — white-label stent deployment simulator: live
  deployment metrics (COF, contact pressure, recoil) and a pressure-sweep
  chart over the `tpt-med-wasm` engine, with vessel compliance and a
  `WHITE_LABEL` rebranding config block.
- `web/demo/` — end-to-end in-browser demo: CT upload → WASM
  segment/mesh → FEM stance solve (×2 decimated real-time preview) →
  CT-measured canal/lumen diameter (mid-slice enclosed-void analysis) →
  stent deployment, with WebGL wireframe rendering and viewer hand-off.
- `tpt-med-wasm` — new APIs: `mesh_csv`, `enclosed_void_diameter`,
  `solve_stance_load_decimated` (max-pool decimation for real-time preview).
- `web/build-wasm.sh` — one-command engine build + wasm-bindgen glue
  generation (version-matched CLI); `web/README.md` run instructions.

### Added — Cross-cutting
- Golden reference datasets (`test-data/golden/{solid,fluid,devices,regulatory}`)
  with documented analytical or literature basis.
- Verification tests: uniaxial tension vs analytical Neo-Hookean, Poiseuille
  flow vs analytical solution, stent radial stiffness benchmark scaffold.
- Benchmark suite (`benches/`) and mdBook user guide (`docs/book/`).
