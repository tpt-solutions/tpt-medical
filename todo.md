# tpt-medical — TODO

**Organization:** TPT Solutions (`tpt-solutions`)
**License:** MIT OR Apache-2.0 (dual)
**Description:** Fully open-source, pure-Rust computational biomechanics and medical device simulation engine, compiling to WASM for zero-cloud, browser-based surgical planning.

---

## Phase 0: Repo Scaffolding & Substrate Status

### Workspace & Tooling
- [x] Root `Cargo.toml` with `[workspace]` + `[workspace.package]` (edition 2021, rust-version, license = "MIT OR Apache-2.0", authors = ["TPT Solutions"], repository)
- [x] `LICENSE-MIT`
- [x] `LICENSE-APACHE`
- [x] `deny.toml` (allow MIT, Apache-2.0, BSD-2-Clause, BSD-3-Clause, ISC, Zlib, Unicode-3.0; deny copyleft/unlicensed)
- [x] `rustfmt.toml`
- [x] `clippy.toml`

### Docs & Community Files
- [x] `README.md`
- [x] `CONTRIBUTING.md`
- [x] `SECURITY.md`
- [x] `CODE_OF_CONDUCT.md`
- [x] `CHANGELOG.md`

### CI/CD
- [x] `.github/workflows/ci.yml` (fmt, clippy -D warnings, test dev+release, build x86_64/aarch64/wasm32, MSRV 1.82 — raised from 1.75: const-float in `tpt-med-units` needs 1.82, `wasm-bindgen` 0.2.128 needs 1.77)
- [x] `.github/workflows/license.yml` (cargo deny check licenses)
- [x] `.github/workflows/benchmark.yml`
- [x] `.github/workflows/docs.yml`
- [x] `.github/workflows/release.yml`
- [x] `.github/ISSUE_TEMPLATE/bug_report.md`
- [x] `.github/ISSUE_TEMPLATE/feature_request.md`
- [x] `.github/ISSUE_TEMPLATE/rfc.md`
- [x] `.github/PULL_REQUEST_TEMPLATE.md`

### Directory Layout
- [x] `crates/core/`, `crates/imaging/`, `crates/solid/`, `crates/fluid/`, `crates/devices/`, `crates/surgical/`, `crates/regulatory/`
- [x] `examples/`
- [x] `test-data/dicom/`, `test-data/nifti/`, `test-data/meshes/`, `test-data/golden/`
- [x] `benches/`
- [x] `docs/book/`, `docs/rfc/`, `docs/api/`
- [x] `rfcs/`

### Substrate Dependency Status (external GitHub repos)
- [x] Confirm/pin `tpt-math` version (linalg-fixed, optimize-general, prob-dist, signal-filter) — **published on crates.io at 0.1.0; pinned `=0.1.0` in root `[workspace.dependencies]`**
- [x] Confirm/pin `tpt-engineering` version (material libraries: PEEK, Nitinol, Titanium, UHMWPE) — **`tpt-eng-materials`/`tpt-eng-biomech` 0.1.0; pinned `=0.1.0`**
- [x] Confirm/pin `tpt-science` version — **`tpt-sci-cfd-core`/`tpt-sci-hemodynamics` 0.1.0; pinned `=0.1.0`. No longer blocks Phase 4** (v0 ships an in-house voxel CFD solver; substrate is the high-fidelity upgrade path)
- [x] Confirm/pin `tpt-fem` version — **`tpt-fem`/`tpt-fem-hyperelastic`/`tpt-fem-contact` 0.1.0; pinned `=0.1.0`. No longer blocks Phase 2** (v0 ships the in-house linear voxel-hex core; substrate integration behind a cargo feature per `rfcs/0002`)
- [x] Re-check substrate maturity before starting each blocked phase below — **done 2026-09-20: all four repos published at 0.1.0; bump policy documented in root Cargo.toml (one PR per bump, V&V re-run)**

## Phase 1: Foundation
*Months 1-3 — ✅ complete*
- [x] `tpt-med-core` — patient models, anatomical coordinate systems, HIPAA-safe audit traits
- [x] `tpt-med-units` — unit system
- [x] `tpt-med-geometry` — geometry primitives
- [x] `tpt-med-dicom` — DICOM parsing, Hounsfield Unit (HU) mapping (uncompressed LE syntaxes; RFC 0001 roadmap for compressed/multiframe)
- [x] `tpt-med-meshing` — voxel-to-hex meshing
- [x] **Milestone:** CLI tool that takes DICOM CT scan, outputs bone mesh CSV (`dicom-to-mesh`)

## Phase 2: Solid Biomechanics
*Months 4-6 — ✅ complete (unblocked; linear core)*
- [x] `tpt-med-biomechanics` — FEM solver core (trilinear hex, 2×2×2 Gauss, CSR + Jacobi CG)
- [x] `tpt-med-bone` — linear elastic bone mechanics
- [x] `tpt-med-tissue` — Neo-Hookean, Mooney-Rivlin (+ Yeoh, Ogden)
- [x] **Milestone:** Femur stress analysis under physiological loading (`femur-stress-analysis`: 13k elements, 36 MPa peak vM at 3×BW)

## Phase 3: Advanced Tissue Models
*Months 7-9 — ✅ complete*
- [x] `tpt-med-tissue` — Holzapfel-Gasser-Ogden (HGO) model for arteries
- [x] `tpt-med-viscoelastic` — Prony series relaxation
- [x] `tpt-med-cartilage` — biphasic/poroelastic models (confined compression, Mow)
- [x] **Milestone:** Arterial wall inflation simulation (HGO verification suite + golden `arterial_wall_inflation.json`)

## Phase 4: Hemodynamics
*Months 10-12 — ✅ complete (unblocked; voxel CFD core)*
- [x] `tpt-med-hemodynamics` — Navier-Stokes CFD (staggered MAC + SOR projection), Wall Shear Stress, Oscillatory Shear Index
- [x] `tpt-med-cardiovascular` — Windkessel model, Fractional Flow Reserve (FFR)
- [x] **Milestone:** Carotid bifurcation Wall Shear Stress calculation (`carotid-wss-screening`: stenosed tube, throat WSS 242% of field mean)

## Phase 5: Implants & Devices
*Months 13-15 — ✅ complete (Level-1 fidelity per RFC 0004)*
- [x] `tpt-med-stents` — Nitinol superelasticity, crimping/expansion deployment
- [x] `tpt-med-orthopedics` — micromotion analysis, stress shielding
- [x] `tpt-med-wear` — Archard/Cross-Land wear laws
- [x] **Milestone:** Stent deployment simulation with artery contact (`stent-deployment`)

## Phase 6: Surgical Planning
*Months 16-18 — ✅ complete*
- [x] `tpt-med-surgical-planning` — osteotomy cuts, virtual surgery
- [x] `tpt-med-implant-sizing` — automated sizing from anatomy
- [x] **Milestone:** Virtual total knee replacement planning (`knee-replacement-planning`)

## Phase 7: Regulatory & Compliance
*Months 19-21 — ✅ complete*
- [x] `tpt-med-fda` — 21 CFR Part 11 immutable, cryptographically signed audit trails
- [x] `tpt-med-vv40` — ASME V&V 40 credibility assessment matrices
- [x] `tpt-med-audit` — cryptographic signature verification (SHA-256 + HMAC-SHA256, FIPS/RFC vectors)
- [x] **Milestone:** Export FDA-submission-ready simulation package (`fda-package` → `test-data/golden/regulatory/fda_export_example.json`)

## Phase 8: WASM & Ecosystem
*Months 22-24 — ✅ complete (alpha)*
- [x] `tpt-med-wasm` — browser compilation of the imaging→mesh→solve stack + stent model
- [x] Web-based surgical planning viewer (WebGL2, no VTK/OpenGL) — `web/viewer/`
- [x] **Milestone:** Zero-cloud, in-browser patient-specific simulation (pipeline validated for `wasm32-unknown-unknown`; wasm-bindgen glue step documented in `docs/book/src/wasm.md`)

---

## Cross-Cutting / Ongoing

### Verification & Validation (ASME V&V 40)
- [x] `test-data/golden/solid/` — femur_loading, lumbar_spine_compression, arterial_wall_inflation
- [x] `test-data/golden/fluid/` — carotid_bifurcation_cfd, aortic_aneurysm_flow, coronary_ffr
- [x] `test-data/golden/devices/` — stent_expansion, hip_stem_micromotion, knee_wear_10mcycles
- [x] `test-data/golden/regulatory/` — fda_audit_trail, vv40_credibility_matrix
- [x] Verification test: uniaxial tension vs. analytical Neo-Hookean solution (deviatoric Cauchy comparison; penalty formulations carry model-internal pressure)
- [x] Verification test: Poiseuille flow vs. analytical CFD solution (paraboloid curvature + WSS tube law)
- [x] Verification test: stent radial stiffness vs. ASTM F2394 published data (scaffold + literature band in place; Level-3 FEM correlation pending)

### ASTM / ISO Standards Automation
- [x] ASTM F2028 — dynamic evaluation of total knee replacements (wear screening path in `tpt-med-wear`)
- [x] ASTM F2079 — intrinsic securement of endovascular stents (radial-force metrics in `tpt-med-stents`)
- [x] ASTM F2394 — securement of self-expanding stents (recoil/dogboning/radial force outputs)
- [x] ISO 7206 — hip joint prostheses (loading convention referenced in femur/hip workflows)
- [x] ISO 14879 — wear of total knee-replacement prostheses (mm³/Mc screening limit + flag)

### RFCs
- [x] `rfcs/0001-dicom-ingestion.md` (Accepted)
- [x] `rfcs/0002-hyperelastic-tissue.md` (Accepted)
- [x] `rfcs/0003-fda-audit-trail.md` (Accepted)
- [x] `rfcs/0004-nitinol-superelasticity.md` (Accepted)
- [x] `rfcs/0005-cardiac-electrophysiology.md` (Draft)

### Business Wedge (stretch, ties to Phase 5 + Phase 8 completion)
- [x] White-label web component: web-based stent deployment simulator for MedTech companies — **`web/stent-simulator/` ships `<tpt-stent-simulator>` (attribute branding, `tpt-deploy` events, npm packaging metadata + README with support path); engine glue in shared `web/pkg/` via `scripts/build-web`**
- [x] Real-time WASM FEM solver demo: CT scan upload → in-browser stent expansion in patient-specific artery — **`web/viewer/` wires the full glue: DICOM series (or synthetic CT demo) → `WasmMeshPipeline` → `wasm_solve_stance_load` + `wasm_deploy_stent`, WebGL2 render of WASM-built meshes; node smoke test verified end-to-end**

---

## Phase 9: Platform Review Follow-ups (2026-09-26)
Findings from a full-workspace review (build/clippy/test run + doc/RFC read).

### Bugs
- [x] `README.md` Quick Start used the wrong CLI flag (`--bone-threshold`
      instead of `--threshold`) for `dicom-to-mesh`. **Fixed** — and the same
      review pass found the adjacent Rust snippet also did not compile
      (missing `Path::new`, and `voxels_to_hex_mesh` is a `&self` method
      returning `Result`); both corrected.
- [x] Per-crate `README.md` files for all 23 members, which had
      `readme = "README.md"` in their manifests but no file, breaking
      `cargo package`. **Done** — every member now has a comprehensive
      README and a CHANGELOG, enforced by `scripts/check-crate-docs.sh`
      (section set, order, and crates.io keyword/category rules).
- [x] Removed `crates/imaging/tpt-med-dicom/src/dbg_test.rs` — an orphaned
      debug test, not wired into the crate by any `mod`, ending in an
      unconditional `panic!`.
- [x] `test-data/nifti/` is scoped rather than removed: it is now a tracked
      roadmap item below (see "Post-Phase 9 roadmap") instead of an orphan
      directory that no code or RFC refers to.

### Governance
- [x] `CONTRIBUTING.md` rewritten for an issues-only model: no external pull
      requests, work arrives as issues, and the RFC process runs through the
      issue tracker with a maintainer committing the accepted RFC.
- [x] `.github/PULL_REQUEST_TEMPLATE.md` reframed for maintainer-authored
      PRs that must reference an issue or accepted RFC; RFC and feature
      issue templates updated to match.

### Larger Initiatives
- [x] Live "try it now" demo of `web/viewer/` on GitHub Pages
      (`.github/workflows/pages.yml`), built from the WASM engine and the
      committed **synthetic** CT data, linked from the README. Separate from
      `docs.yml` because a WASM + wasm-bindgen build is far more expensive
      than rustdoc and only needs re-running when the engine changes.
- [x] `templates/` — a crate template, an example-binary template, and an RFC
      template, based on the `examples/src/bin/*.rs` milestones. Deliberately
      copy-based rather than `cargo generate`: the thing worth templating is
      the house style, and the contract the templates promise is the same one
      CI enforces on real crates.
- [x] Golden-dataset CI diffing tool (`scripts/diff-golden.sh`): renders a
      before/after numeric-drift table for `test-data/golden/` against a PR,
      with each row judged against the `tolerance_percent` the golden file
      itself declares. Wired into CI as the `golden-drift` job.
- [x] `benchmark.yml` wired to a stored baseline (`benches/baseline.txt`) via
      `scripts/bench-baseline.sh`, so a regression fails CI instead of just
      being printed. Threshold is deliberately loose (25 % by default,
      `TPT_BENCH_THRESHOLD_PCT`) to catch algorithmic regressions without
      failing on shared-runner jitter. Regeneration is a manual
      `workflow_dispatch` act.
- [x] Reproducibility manifest per simulation run in `tpt-med-fda`:
      `ReproducibilityManifest` pins the workspace version, per-crate
      versions, git commit, build profile, and a SHA-256 of every input
      artefact. `AuditTrail::attach_manifest` records its digest as an audit
      event, so a swapped manifest is detectable from the chain alone, and it
      is covered by the export's detached HMAC. Wired into the `fda-package`
      milestone. **Known gap:** it records the *workspace* version, which is
      only exact while the workspace versions as a unit — see the crate's
      README.

---

## Post-Phase 9 roadmap
Scoped but not started. Nothing here is implemented; each is a real
follow-up rather than an orphan.

- **NIfTI ingestion** (`tpt-med-nifti`) — research-space volumes alongside
  DICOM, so a NIfTI export from 3D Slicer or a public dataset can enter the
  pipeline without an external conversion step. `test-data/nifti/` is
  reserved for its synthetic fixtures. Needs an RFC; the RAS coordinate
  handling already exists in `tpt-med-geometry`.
- **Compressed DICOM transfer syntaxes** — JPEG, JPEG-LS, JPEG 2000 and RLE
  pixel data, currently rejected with `DicomError::CompressedPixelData`. The
  single biggest practical gap: many clinical archives store CT as JPEG 2000.
  Must land behind a named feature with `cargo deny`-approved dependencies.
- **Per-crate versions in the reproducibility manifest** — needs a build
  script or a generated version table, before the independent-release cadence
  makes the workspace version wrong.
- **Quantitative CT calibration** — replace the linear HU→density
  approximation with a phantom-calibrated relation, so absolute density and
  modulus stop being screening estimates.
- **Nonlinear FEM integration** — `tpt-fem` / `tpt-fem-hyperelastic` /
  `tpt-fem-contact` behind a cargo feature, per `rfcs/0002`.
- **Cardiac electrophysiology** — `rfcs/0005-cardiac-electrophysiology.md`
  is still Draft.

