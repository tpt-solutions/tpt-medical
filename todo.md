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
- [x] `rfcs/0005-cardiac-electrophysiology.md` (Accepted 2026-09-27 — Stage 1
      fleshed out to concrete API/numerics/verification; Stages 2/3 remain
      roadmap-depth by design and are not authorized by this acceptance)
- [x] `rfcs/0007-bmd-apparent-density-conversion.md` (Accepted 2026-09-27)
- [x] `rfcs/0008-phantom-rod-sampling.md` (Accepted 2026-09-27)

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
Scoped but not started unless marked done. Each is a real follow-up rather
than an orphan.

- [x] **Per-crate versions in the reproducibility manifest** — `examples/build.rs`
  parses the workspace `Cargo.lock` into `tpt_med_examples::crate_versions`, a
  generated `(crate name, resolved version)` table; `fda-package` looks up
  each participating crate's real version there instead of stamping every
  entry with the calling binary's own `CARGO_PKG_VERSION`. Dependency-free
  (line-scans `Cargo.lock` rather than pulling in a TOML parser or
  `cargo_metadata`). `ReproducibilityManifest` itself already supported
  per-crate versions (`with_crates`); the gap was only in how the example
  populated it. See `tpt-med-fda`'s README "Known Limitations" for the
  remaining caveat: a *caller* still has to use a correct source of versions,
  the crate cannot enforce that.
- [x] **NIfTI ingestion** (`tpt-med-nifti`, `rfcs/0006-nifti-ingestion.md`) —
  pure-Rust parsing of uncompressed single-file NIfTI-1 (`.nii`) volumes:
  sform/qform geometry (RAS, spec precedence), 8 datatypes,
  `scl_slope`/`scl_inter` scaling. `SyntheticNiftiBuilder` fixtures, no real
  dataset needed. Three explicit follow-ups this RFC deliberately left
  open, tracked below.
- [x] **`.nii.gz` support** — gzip is by far the most common NIfTI file
  extension in practice; v0 rejects it with `NiftiError::Gzipped` rather
  than decompress it. Needs a vetted pure-Rust inflate dependency, the same
  kind of decision `tpt-med-dicom`'s RLE/JPEG/JPEG-LS/JPEG 2000 features
  already made three times over. **Done (2026-09-27)** — `gzip` cargo
  feature (off by default) on `tpt-med-nifti`, `flate2` with the pure-Rust
  `rust_backend` (miniz_oxide; MIT OR Apache-2.0; `rust-version` 1.67,
  under the workspace MSRV). Single-file decompression is streaming and
  geometry-capped (the header's declared extent bounds what is ever
  allocated), the typed `Gzipped` error stays for feature-off builds, a
  new `CorruptGzip` names decompression failures, and gzipped `.hdr`/`.img`
  parts work off the same feature. CI clippy/test gained a
  `-p tpt-med-nifti --features gzip` pass.
- [x] **`tpt-med-nifti` → `tpt-med-meshing` integration** —
  `SegmentationMask::threshold_hu` is concretely typed to
  `tpt_med_dicom::DicomSeries` today. Wiring a `NiftiVolume` through it needs
  its own API-design decision (a shared trait? an adapter?), deliberately
  left open by RFC 0006 rather than bundled into the ingestion RFC.
  **Done (2026-09-27)** — `SegmentationMask::threshold_nifti(&NiftiVolume,
  min_hu)`, a direct constructor on the meshing side (no shared trait:
  `NiftiVolume` is still the only non-DICOM source, exactly the narrower
  shape the RFC said to pick). NIfTI's RAS origin/directions are converted
  through `ras_to_lps` into the LPS patient frame `threshold_hu` produces,
  so masks and meshes are interchangeable across sources; `tpt-med-meshing`
  gained the `tpt-med-nifti` dependency; tests cover the threshold result,
  the frame conversion (including the `voxel_position`/`voxel_center`
  invariant), and an end-to-end `voxels_to_hex_mesh` run.
- [x] **Dual-file `.hdr`/`.img` NIfTI-1 support** — lower priority than
  `.nii.gz` (rarer in current tooling), real gap if a workspace member's
  dataset uses it. **Done (2026-09-27)** —
  `NiftiVolume::{parse_dual_file, parse_dual_bytes}` accept the `ni1`
  magic, honour `vox_offset` as an offset into the `.img` (0 is the norm),
  and reject the other layout's magic with a pointer to the entry point
  that can read it; `SyntheticNiftiBuilder::build_dual` writes the pair for
  tests; gzipped parts work when `gzip` is on. Header/value decoding was
  factored into one shared `decode_header`/`decode_values` so the
  single-file, dual-file, and gzip paths cannot drift apart.
- [x] **RLE Lossless pixel data** (`1.2.840.10008.1.2.5`) — decoded behind the
  new `rle` cargo feature in `tpt-med-dicom` (`src/rle.rs`): PackBits
  segments per PS3.5 Annex G, geometry-capped so a corrupt run can't expand
  unbounded, single-frame only (a non-empty Basic Offset Table is rejected).
  Needs no new dependency, so nothing new for `cargo deny` to approve. The
  parser now recognises the encapsulated pixel-data item stream (PS3.5 Annex
  A.4) generally, so this is the template for the remaining codecs below.
- [x] **JPEG-family compressed DICOM transfer syntaxes** — three new opt-in
  cargo features on `tpt-med-dicom`, one per codec, each off by default:
  - `jpeg`: Baseline/Extended DCT (lossy, `.50`/`.51`) and Lossless Process
    14/SV1 (exact, `.57`/`.70`), via [`jpeg-decoder`](https://crates.io/crates/jpeg-decoder)
    (image-rs — the same codec crate covers both coding processes).
  - `jpeg-ls`: Lossless and Near-Lossless (`.80`/`.81`), via
    [`pure_jpegls`](https://crates.io/crates/pure_jpegls).
  - `jpeg2000`: Lossless Only and lossless-or-lossy (`.90`/`.91`), via
    [`pdfluent-jpeg2000`](https://crates.io/crates/pdfluent-jpeg2000)
    (`image`/`simd` extras disabled, so no further dependencies pulled in).

  All three reuse the encapsulated-fragment collection the `rle` feature's
  parser work already generalised — no parser changes needed. See
  `tpt-med-dicom`'s CHANGELOG and README for per-codec detail, and the
  follow-up below that this work surfaced rather than closed out.
- [x] **`jpeg2000` feature: fix the signed-`PixelRepresentation` gap** —
  `pdfluent-jpeg2000` applies the unsigned DC level-shift to every component
  regardless of the codestream's own signed bit (its source says so
  explicitly). Since the shift is a fixed, known offset (`2^(precision-1)`),
  `jpeg2000::decode_frame` now re-reads the SIZ marker's `Ssiz` byte directly
  from the raw codestream (the crate discards it) and undoes the shift itself
  when the component really is signed, instead of rejecting all signed
  objects outright. A file where the codestream's own signed bit and the
  dataset's `PixelRepresentation` disagree is non-conformant and still
  rejected, since there is no safe way to resolve that disagreement.
- **JPEG 2000 Part 2 multi-component and JPIP-referenced pixel data** — still
  rejected with `DicomError::CompressedPixelData`; no decoder exists for
  either. Multi-component is a real but uncommon gap (multi-channel or
  wavelet-transformed color JPEG 2000); JPIP is a network reference to pixel
  data elsewhere, not pixel data itself, and would need its own transport
  story before decoding matters.
- [x] **Quantitative CT calibration (fit)** — `QctCalibration::fit` in
  `tpt-med-dicom` fits a real HU→density line by ordinary least squares from
  a calibration phantom's measured `(HU, known_value)` rod points, replacing
  the fixed two-point screening line for anyone who supplies real
  measurements. `HounsfieldMapper::hu_to_density` is unchanged (still that
  fixed line by default; the new type is additive, not a behaviour change —
  no RFC needed for that reason). Deliberately does **not** solve the two
  harder problems the roadmap item implied and which still need real
  scoping work — see the two follow-ups below.
- [x] **Automatic calibration-phantom rod detection** — `QctCalibration::fit`
  needs `(HU, known_value)` points handed to it; nothing in this crate finds
  a calibration phantom in a series or samples its rod ROIs. Needs an RFC:
  phantom geometry varies by manufacturer (Mindways QCT Pro, CIRS/Image
  Analysis, …) and there's no single detection algorithm across them.
  **RFC accepted and implemented (2026-09-27):**
  `rfcs/0008-phantom-rod-sampling.md` → `tpt-med-dicom/src/phantom.rs` —
  `locate_phantom_centroid` (manufacturer-agnostic, thresholded connected-
  component centroid) plus a caller-supplied, cited `PhantomModel` (rod
  layout + known values) sampled by `sample_phantom_rods`; rotation and
  slice selection stay caller-supplied in v0 rather than blind-detected. 13
  new tests, `cargo test`/`clippy`/`fmt` clean, README/CHANGELOG updated.
- [x] **BMD → apparent-density conversion** — most clinical QCT phantoms report
  rod values as bone mineral density (mg/cm³ K₂HPO₄- or CaHA-equivalent), not
  apparent (whole-tissue) density, which is what the Morgan–Keaveny modulus
  power laws in `hounsfield.rs` expect. Converting one to the other needs a
  documented, protocol-specific relation (the literature has more than one);
  `QctCalibration` deliberately does not pick one silently. **RFC accepted
  and implemented (2026-09-27):** `rfcs/0007-bmd-apparent-density-conversion.md`
  → `tpt-med-dicom/src/bmd.rs` — two-stage `BmdToAshDensity` → `AshFraction`
  conversion, composed as `BmdToApparentDensity`; ships **no** built-in
  preset relations, and both stages refuse to construct without a non-empty
  citation string. `HounsfieldMapper::bmd_to_youngs_modulus` is the
  composing convenience. 8 new tests, `cargo test`/`clippy`/`fmt` clean,
  README/CHANGELOG updated.
- [x] **Nonlinear FEM integration — scoped and given a first safe slice
  (2026-09-27)**. Read the pinned substrate crates' actual 0.1.0 APIs
  (`tpt-fem-hyperelastic`, `tpt-fem-mesh`, `tpt-fem-element`,
  `tpt-fem-contact`, `tpt-fem-solve`) and found the real gap RFC 0002/0004
  didn't detail: no 3D `Hex8` nonlinear hyperelastic assembly exists in the
  substrate yet (only stress functions + a 1-D bar Newton solve); `Hex8`
  element/mesh support exists, contact is DOF-level-constraint-only with no
  friction. Wrote `rfcs/0009-nonlinear-fem-substrate-adapter.md` (Draft) to
  scope the real adapter-crate architecture for that gap — not implemented,
  by design (see the RFC's Drawbacks). Implemented the narrow, safe first
  slice instead: `tpt-med-tissue`'s new `substrate-cross-check` cargo
  feature (off by default) cross-checks the in-house closed-form
  incompressible-Neo-Hookean uniaxial stress against
  `tpt-fem-hyperelastic::solve_hyperelastic_bar`'s independent 1-D bar
  Newton solve, agreeing to `1e-9`. Adds `tpt-fem-hyperelastic`/`tpt-fem-mesh`
  as optional deps (pinned `=0.1.0` in root `[workspace.dependencies]`,
  `cargo deny check licenses` clean); zero effect on the default build.
  `cargo test`/`clippy`/`fmt` clean workspace-wide with the feature on and
  off, `check-crate-docs.sh` passes for all 25 members.
- [x] **Cardiac electrophysiology RFC fleshed out (2026-09-27)** — Stage 1
  (monodomain on voxel geometry, `tpt-med-electrophysiology`) now has a
  concrete API sketch, Mitchell–Schaeffer kinetics with a cited default
  parameter set, explicit numerics (7-point Laplacian, RK2, CFL stability
  bound), and a full verification-strategy section (manufactured-solution
  stencil check, grid-convergence CV, single-cell restitution shape check).
  Stages 2/3 (ECG forward problem, ablation screening) deliberately stay at
  roadmap depth pending Stage 1 actually shipping. `rfcs/0005-cardiac-electrophysiology.md`
  is **Accepted and implemented (2026-09-27)** for Stage 1 only; Stages 2/3
  are not authorized by this acceptance. New crate
  `crates/fluid/tpt-med-electrophysiology`: `MitchellSchaefferParams`,
  `MonodomainTissue` (explicit RK2, diffusion-CFL + reaction-stiffness
  stability bound — the latter added after development caught the
  diffusion-only bound letting the reaction term diverge to `NaN`), and
  `S1S2Protocol` (single-cell restitution, the code-verification fixture).
  New golden dataset `test-data/golden/electrophysiology/monodomain_restitution.json`.
  15 unit tests + 1 doctest, `cargo test`/`clippy`/`fmt` clean workspace-wide,
  README/CHANGELOG written, `check-crate-docs.sh` passes for all 25 members.

## In progress / newly tracked (2026-09-27)

- [x] **Finish verifying the JPEG 2000 Part 2 multi-component change.**
  Implementation is in: `tpt-med-dicom` gained `TransferSyntax::Jpeg2000Part2MultiComponentLossless`/
  `Jpeg2000Part2MultiComponent` (`1.2.840.10008.1.2.4.92`/`.93`), routed
  through the existing `jpeg2000::decode_frame` path (Part 2 extends Part
  1's codestream syntax rather than replacing it; verified by reading
  `pdfluent-jpeg2000`'s marker-parsing loop, which rejects any marker code
  it does not recognise rather than skipping it, so a genuine Part 2
  extended multi-component transform is refused, not mis-decoded). Adds
  tests in `tags.rs` and updates `jpeg2000.rs`/README/CHANGELOG docs.
  **Closed 2026-09-27:** all four commands rerun clean
  (`cargo test -p tpt-med-dicom --features jpeg2000` 56 passed,
  `cargo clippy --all-targets --features jpeg2000 -- -D warnings`,
  `cargo fmt -- --check`, `scripts/check-crate-docs.sh` 25/25), and the
  default-feature build still passes (46 tests) with the new
  `#[cfg(test)]`-gated fixture. The noted asymmetry is closed too:
  `series::jpeg2000_encapsulated_pixel_data_tests` gives JPEG 2000 the
  same Part-10-byte-stream-to-`DicomSlice` coverage `jpeg-ls` already had,
  for `.90`/`.91`/`.92`/`.93` alike — including the signed-bit level-shift
  correction and the non-conformant-mismatch rejection on the real parse
  path, not just `decode_frame`'s unit tests. The 2x2 J2C fixture is now
  one `pub(crate)` constant in `jpeg2000.rs` shared by both test sites
  rather than a near-identical hand-rolled copy per site.
- [x] **Implement RFC 0009's nonlinear FEM adapter, including contact
  coupling** (user-selected scope, 2026-09-27). New crate
  `crates/solid/tpt-med-fem-adapter` (RFC 0009's placeholder name kept; the
  crate's shape now matches the name), 17 tests, `cargo test`/`clippy`/`fmt`
  clean, `scripts/check-crate-docs.sh` 26/26. What shipped, against RFC 0009's
  four items:
  1. **3-D `Hex8` nonlinear hyperelastic assembly** (item 1) —
     `tpt-fem-element` shape functions and `hex_rule` quadrature, with
     `tpt-med-tissue::TissueModel::first_piola` per quadrature point. Total
     Lagrangian, penalty incompressibility as the RFC recommended.
  2. **Tangent stiffness** (item 2) — the RFC's open "analytic vs. numerical"
     question is settled *with numbers*: `tangent_stiffness` differentiates the
     constitutive law only (`A = dP/dF` by central differences) and
     `tangent_stiffness_numerical` differentiates the whole residual; the
     suite checks them against each other on a deformed configuration and checks
     minor symmetry. The cheap one is what the solver uses.
  3. **`tpt-fem-sparse::Coo`/`solve` wiring into a nonlinear equilibrium
     iteration** (item 3) — with one documented deviation: the loop is this
     crate's, not `tpt-fem-solve::newton`'s. Read from the pinned source, that
     driver tests the *full* residual against an *absolute* tolerance, which a
     displacement-controlled problem can never satisfy (the residual at a
     prescribed DOF is the reaction, non-zero by definition), so it always
     reports `MaxIterations` at a perfectly converged solution. Same structure
     (condense the essential DOFs, solve, update), free-DOF convergence measure,
     plus a line search and diagonal equilibration before the solve.
  4. **Contact coupling** (item 4) — `tpt-fem-contact`'s `contact_pairs` selects
     the contact partner and `penalty_contact` supplies the linearised penalty,
     with the active set recomputed from the current geometry at every residual
     and Jacobian evaluation. The RFC's open coupling design is thereby settled:
     a frozen active set is wrong, and separation must release the constraint
     within the same solve (both are tested).
  **Verification, the RFC's other placeholder obligation, is discharged**: code
  verification against the same closed form `tpt-med-tissue` and
  `tpt-fem-hyperelastic` both reproduce (`mu (lambda - lambda^-2)`), the
  exact constant-stress patch identity, the exact analytic volumetric branch,
  and calculation verification by mesh refinement (measured ratios 1.233, 1.058,
  1.019, 1.003 for 1x1x1 .. 4x4x4).
  **Three findings worth keeping** (all now in the crate's README/CHANGELOG):
  the in-house volumetric convention is `W_vol = (J-1)^2/d1`, so a *small* `d1`
  is a *stiff* penalty — the large `d1` used by `substrate-cross-check` would be
  nearly compressible; a stiff penalty plus full integration locks volumetrically,
  which is why the fixture's `d1 = 0.5` is a measured compromise; and an
  inverted element produces a `NaN` from the model's `J^-2/3` that the linear
  solver then reports as an unrelated "singular matrix", so inversion is now an
  explicit error.
  Still open, tracked in the new crate's CHANGELOG: friction, a mixed `u`-`p`
  formulation to remove the locking, load stepping, and curved elements.