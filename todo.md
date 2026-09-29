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
- [x] **JPEG 2000 Part 2 multi-component transfer syntaxes** (`.92`/`.93`) —
  routed through the same `jpeg2000::decode_frame` path as `.90`/`.91`, closed
  below. What remains open is the codestream work, split out into the two
  follow-ups immediately after.
- [ ] **JPEG 2000 Part 2 multi-component codestream decoding** — the transfer
  syntaxes decode, but "multi-component" in their name is not a decode
  capability. A codestream using a real Part 2 extended multi-component
  transform (multi-channel or wavelet-transformed color) is still refused —
  `pdfluent-jpeg2000` rejects markers it does not recognise rather than skipping
  them, and this crate has no `SamplesPerPixel`/`PlanarConfiguration` handling
  anywhere. A genuine gap, uncommon for this crate's CT/MR HU scope.
  **Spiked 2026-09-27** (`rfcs/0010-jpeg2000-part2-spike.md`): the rejection is
  verified at source, no pure-Rust route exists today, and the gating input is
  a purchase of ISO/IEC 15444-2 rather than engineering time. Recommend
  deferring; the spike records the decision order if it is ever revived.
- [ ] **JPIP-referenced pixel data** — still rejected with
  `DicomError::CompressedPixelData`; no decoder exists, and that is the correct
  outcome. JPIP is a network reference to pixel data held elsewhere, not a local
  format, so this needs a transport story (a JPIP client, resolved at the
  archive boundary) before any decoding question is even in scope.
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
  Still open, tracked in the new crate's CHANGELOG: a mixed `u`-`p` formulation to
  the locking. Friction, load stepping and curved (`Hex20`/`Hex27`/`Tet10`)
  elements have since landed there, so that CHANGELOG is the authority on what
  is still missing.

## Aspirational backlog (per-crate `Planned` sections)

The phased roadmap above is complete, but each crate also carries a
`### Planned` section in its own CHANGELOG recording work that is
scoped-but-unstarted. Those lived only in twenty separate files, which
made "what is left?" a question needing twenty file reads to answer.
They are consolidated here as one checkbox each, and the CHANGELOGs stay
authoritative — this section is the tracker, not the source of truth.

None of these block the roadmap. They are listed so they are tracked
rather than forgotten, and so the size of what remains is visible rather
than inferred from a clean-looking list of completed phases.

### `tpt-med-dicom`

The ingestion roadmap's own open items — Part 2 multi-component
codestream decoding and JPIP — are tracked above rather than repeated
here, so there is one checkbox per piece of work. Part 2 was spiked in
`rfcs/0010-jpeg2000-part2-spike.md`, which recommends deferring it: the
gating input is a purchase of ISO/IEC 15444-2, not engineering time.

- [x] **Automatic rotation detection and a small built-in library of named, cited `PhantomModel`s for common commercial phantoms — explicitly out of scope for `rfcs/0008-phantom-rod-sampling.md`'s v0 mechanism; see that RFC's Unresolved Questions.**
  **Detection half delivered (2026-10-01): `detect_phantom_rotation` — least-squares sweep + grid-and-shrink refinement over in-plane rotations about the located centroid, RMS residual as the reject signal. The vendor-phantom library stays open below.**
- [ ] **Built-in library of named, cited `PhantomModel`s for common commercial phantoms (the remaining half of the rotation-detection item) — must come from manufacturer datasheets with citations, not baked-in approximations.**

### `tpt-med-audit`

- [x] **A `verify_chain_detailed` returning per-index status, so a forensic tool can report the first broken link rather than only that the chain fails.**
  **Done (2026-09-29): `verify_chain_detailed`/`LinkReport`/`LinkStatus` in `tpt-med-audit` — first broken link, expected-vs-stored digests, truncation as `Missing`.**
- [x] **Asymmetric signatures (Ed25519) for non-repudiation, behind a clearly named feature, alongside the existing symmetric HMAC path.**
  **Done (2026-09-30): `ed25519` cargo feature (off by default) — `SigningKey`/`VerifyingKey` over `ed25519-dalek`, pinned against the RFC 8032 §7.1 TEST 2 constants (cross-checked against an independent OpenSSL derivation when recorded); CI gained a feature-gated pass.**
- [x] **External anchoring helpers for RFC 3161 timestamping or a transparency log, which RFC 0003 identifies as the way to close the non-repudiation gap.**
  **Done (2026-09-30): `anchor` module — `AnchorRecord::anchor` binds a payload's SHA-256 digest to a TSA/transparency-log submission, `covers` checks the binding mechanically, and the token stays opaque (service-side verification, by design).**

### `tpt-med-biomechanics`

- [x] **Nonlinear capability behind a cargo feature via `tpt-fem` / `tpt-fem-hyperelastic` / `tpt-fem-contact`, per RFC 0002. The linear core stays the default so the WASM footprint is unchanged.**
  **Closed 2026-09-30 by architectural decision: the capability shipped as the `tpt-med-fem-adapter` crate (RFC 0009's accepted implementation). A second Newton assembly inside `biomechanics` was rejected in review — the crate's CHANGELOG Notes carries the rationale; this crate contributes the small-strain inclusion path instead.**
- [x] **Multi-constraint boundary conditions (symmetry planes, roller constraints) as a convenience over the current all-3-DOFs `fix_nodes`.**
  **Done (2026-09-29): `BoundaryConditions::constrain_dofs(nodes, [x,y,z])` — verified by an exact symmetry-plane uniaxial state.**
- [x] **Direct `tpt-med-tissue` material support, so a single model can mix linear bone and hyperelastic soft tissue.**
  **Done (2026-09-30): `ElementMaterial::SoftTissue` + `from_parts_mixed` — a tissue model linearized at `F = I` (`linearized_engineering_constants`) joins linear bone in one solve; `ν ≥ ½` parameter sets are build errors, not clamps.**
- [x] **Grid and time convergence reporting (`CalculationVerification` evidence) as a first-class result rather than an off-script exercise.**
  **Done (2026-09-29): `convergence::convergence_study` + `ConvergenceReport` — sorted levels, relative errors, observed order (locked to 2.0 on an h² series).**

### `tpt-med-bone`

- [x] **Spatial remodeling: drive per-element density from a solved strain energy density field rather than the current single lumped value per call.**
  **Done (2026-09-29): `BoneRemodelingModel::remodel_field(densities, stimuli, dt, viable)` — per-voxel update with viable clamping.**
- [x] **A disuse/resorption-deadline model, and temperature- or load-rate-dependent remodeling.**
  **Done (2026-09-30): `ResorptionDeadline` (+ `update_density_with_deadline`/`remodel_field_with_deadline`) — caller-held per-voxel disuse counter, reset by reloading, multiplier past the deadline; and `RateAugmentation` — log-scaled, saturation-capped stimulus multiplier above a reference load rate. Temperature dependence documented out (immaterial at core temperature for screening cases).**
- [x] **QCT phantom calibration hooks, so a study can supply a calibrated density→modulus relation instead of the default power law.**
  **Done (2026-09-29): `ModulusLaw` trait + `PowerLaw` + `BoneMaterial::from_hu_with_law` — closures and phantom-fitted laws accepted.**

### `tpt-med-cardiovascular`

- [x] **Four-element and non-linear pressure–flow relations for systemic circulation modelling.**
  **Done (2026-09-30): `FourElementWindkessel` (inertance `L` in the pressure-driven two-state form — the formulation where the inertance adds ring-down physics; verified against the theoretical envelope decay and DC-gain superposition) + `WaterfallResistor` (Starling-resistor non-linear pressure–flow). Flow-driven CFD coupling stays 3-element, documented.**
- [x] **Waveform-based instantaneous-hyperbolic FFR, alongside the pressure-ratio definition implemented here.**
  **Done (2026-09-30): `InstantaneousWaveFreeRatio` — Pd/Pa averaged over the wave-free diastolic window (caller-supplied fractions, 0.90 threshold), the resting waveform-based counterpart to the hyperemic pressure ratio.**
- [x] **Patient-specific waveform fitting, rather than the fixed analytic shapes.**
  **Done (2026-09-30): `MeasuredFlowWaveform::fit` — truncated Fourier series (mean + per-harmonic amplitude/phase) fitted by real DFT to sampled flow, with `fit_rms` as the truncation measure; verified by exact recovery of the analytic carotid series.**
- [x] **A direct coupling helper so a `tpt-med-hemodynamics` solve can step the Windkessel state in lockstep with the CFD time step.**
  **Done (2026-09-29): `CoupledWindkessel` — stateful RK4 advance returning the outlet pressure per CFD step; steady-state + decay tests.**

### `tpt-med-cartilage`

- [ ] **Additional boundary conditions: unconfined compression and shear.**
  **Shear half delivered (2026-10-01): `solid_shear_modulus` — first-order biphasic shear is volumetrically silent, so the fluid never pressurises and the response is a closed form at all times (permeability-independent; verified). Unconfined compression remains open above.**
  **Unconfined limits delivered (2026-10-01): `unconfined_equilibrium_modulus` (E_s, equal to H_A at ν_s = 0) plus the documented rigid instantaneous response — both exact and tested. The transient itself remains the classical Bessel-series solution.**
- [ ] **Nonlinear biphasic theory, and a coupling between permeability and strain.**
- [ ] **A lubrication/repulsion term for the contact interface, so the model can be driven by a contact solver rather than a prescribed step load.**
- [ ] **Fibrous-cartilage support (a fibre-reinforced solid matrix).**

### `tpt-med-electrophysiology`

- [ ] **Stage 2 (ECG/EGM forward problem via pseudo-bidomain lead-field projection) and Stage 3 (ablation screening), both kept at roadmap depth in `rfcs/0005-cardiac-electrophysiology.md` pending Stage 1 usage and, for Stage 3, clinical-data validation.**
- [x] **Anisotropic (fiber-direction) conductivity — needs a fiber-field source (atlas or DTI derivation) this crate has no source for yet.**
  **Done (2026-09-30): the *acceptance* half — `FiberConductivity` + `set_anisotropy` take a caller-supplied per-voxel field (axisymmetric projections, whole-field validation, stability bound from the fiber maximum; isotropic-limit equivalence verified). The *source* stays external by design, exactly as the item framed it.**
- [ ] **Promotion path to `tpt-science`'s electrophysiology crate for ionic-model breadth beyond Mitchell-Schaeffer, once a real workflow needs it.**

### `tpt-med-fda`

- [x] **Enforce "sign after the final edit": today signing, appending and exporting is permitted and the export records the ordering, so a workflow requiring the stricter discipline must enforce it itself.**
  **Done (2026-09-29): `SignaturePolicy::RequireSignatureAfterLastEdit` + `checked_append`/`export_package_checked` with `PolicyError`.**
- [x] **An append-only persistence layer with WORM semantics, so the trail survives a process restart without a caller-supplied store.**
  **Done (2026-09-30): `worm::WormLog` — write-once journal (`create_new`), per-entry chain digest + `fsync` per append, reopen verifies header/chain/sequence; torn tails refused as `TornTail`, retroactive edits and cross-run splices as `ChainBroken` at the first bad link; `into_trail` rebuilds a self-verifying trail. Signatures/policies remain trail-side state.**
- [x] **External anchoring of the detached tag (HSM, transparency log, RFC 3161) to close the non-repudiation gap identified in RFC 0003.**
  **Done (2026-09-30): `AuditTrail::attach_anchor` — digest binding checked at attach, the anchor recorded as an audit event (removal breaks the chain), `PolicyError::AnchorDigestMismatch` on swapped payloads.**
- [x] **Reason-field policy enforcement, so a site can require a structured reason code rather than free text.**
  **Done (2026-09-29): `ReasonPolicy::structured([...])` — code-prefix validation in `checked_append`.**

### `tpt-med-fem-adapter`

- [ ] **A mixed `u`-`p` formulation for exact incompressibility. The deviatoric/volumetric split it requires now exists** (added for selective reduced integration, `tpt-med-tissue::TissueModel::volumetric_first_piola`), **so this is "add a pressure unknown" rather than "redesign the trait". The remaining decision is the inf-sup-stable element pairing — `Hex8`/constant pressure, or the `Hex20`/`Hex8` pairing RFC 0009 named — a numerical-methods call, not a mechanical one.**
  **Prototype attempted and withdrawn (2026-10-01): the residual assembly (mean-dilatation F̄, pointwise-cof constraint stress, global pressure DOFs) is settled and recorded in `rfcs/0012`'s new implementation-notes section — including two assembly bugs and their consistency-check tooling. The remaining blocker is the nonlinear saddle solver's global convergence (first-step pressure explosion, hourglass-branch drift); three candidate resolutions are ranked in the RFC. Implementation stays open pending that solver design.**
  **Designed (2026-10-01): `rfcs/0012-mixed-up-formulation.md` (Draft) — records why element-level condensation is impossible (zero Lagrange diagonal ⇒ global pressure DOFs through the shared Newton loop), recommends Q1/P0 on the existing Hex8 for structured voxel meshes, and lays out the five-part verification (closed-form incompressible uniaxial, patch test, locking benchmark vs SRI, spurious-mode check, SRI cross-validation). Implementation pending acceptance.**

### `tpt-med-hemodynamics`

- [x] **Conjugate-gradient or multigrid pressure solve to replace the Jacobi sweeps, which currently dominate the per-step cost.**
  **Done (2026-09-30): `PressureSolver::ConjugateGradient` — Jacobi-preconditioned, matrix-free on the masked grid, verified against a manufactured solution of the discrete operator (which caught a real sign bug: the SOR fixed point is the SPD negative Laplacian, so CG's rhs enters negated) and by the full Poiseuille march under CG. ~40 iterations vs SOR's 400-sweep cap on the verification tube. Default stays SOR so golden datasets and the bench baseline are unchanged; multigrid remains future work.**
- [ ] **Conjugate heat transfer and wall compliance, enabling a coupled fluid–structure boundary.**
- [ ] **Optional local wall refinement, so peak WSS at a geometric corner stops being resolution dependent.**
- [x] **Coupling to `tpt-med-cardiovascular` for a driven, rather than prescribed, boundary condition.**
  **Done (2026-10-01): `step_coupled` — the Windkessel's pressure is imposed as the projection's Dirichlet anchor (both Poisson solvers start from the anchor level, fixing a real spurious-gradient defect at large anchors) and the measured outlet flow advances the 0-D model; verified by exact replay of the flow history through the boundary model's RK4. Explicit staggered coupling; the reported field is gauge-relative to the boundary reference.**

### `tpt-med-implant-sizing`

- [ ] **Automatic landmark detection from a CT, which is the hard part of the problem and is not attempted here.**
- [x] **Multi-measurement charts, so a femoral decision can weigh TEA, AP depth and posterior condylar offset jointly, with a vendor-specific precedence rule.**
  **Done (2026-09-30): `size_from_measurements` — every measurement votes, disagreement resolves by lowest precedence rank (equal ranks up-size), every vote recorded with out-of-chart flags; `posterior_condylar_offset` + `femoral_measurements()` supply the TEA/AP/PCO triple (PCO as a documented TEA-perpendicular proxy).**
- [x] **Soft-tissue and ligament balance assessment, and a check that the selected size leaves acceptable gap balancing. Sizing is necessary for a good plan and not sufficient.**
  **Gap-check half delivered (2026-10-01): `check_gap_balance` — extension/flexion gaps from resection-vs-thickness arithmetic, overstuffed and imbalance flags. Ligament tension proper stays open (needs soft-tissue structures).**
- [ ] **Ligament balance assessment proper (tension and stability), split out of the delivered gap-balance check — needs soft-tissue structures the crate does not model.**
- [ ] **Hip, shoulder and ankle sizing beyond the knee-specific `KneeLandmarks`.**
- [x] **Schema validation and reporting for a caller-supplied chart, so a mis-transcribed chart is caught rather than silently producing a recommendation.**
  **Done (2026-09-29): `SizeChart::validate` + `ChartError` — empty/non-positive/unsorted/duplicate-label detection.**

### `tpt-med-meshing`

- [x] **Per-voxel material overrides, so a caller can supply a QCT-calibrated or region-specific modulus instead of the default HU correlation.**
  **Done (2026-09-29): `MedicalMesher::voxels_to_hex_mesh_with_overrides` — per-voxel modulus map over the HU correlation.**
- [x] **Optional node deduplication across disconnected components.**
  **Done (2026-09-29): `VoxelHexMesh::weld_nodes(tolerance)` — position-hashed welding with connectivity remap.**

### `tpt-med-nifti`

- [ ] **NIfTI-2 (the 2011 540-byte header with 64-bit dimensions) — still out of scope per RFC 0006; revisit if a workspace dataset needs it.**

### `tpt-med-orthopedics`

- [x] **Cyclic loading: micromotion accumulated over a gait cycle rather than evaluated at a single static load.**
  **Done (2026-09-30): `micromotion_over_cycle` — peak plus per-zone motion amplitude (`max − min`) across a full load cycle; `GaitCycle::iso_double_hump` supplies an ISO 14243-style screening waveform (≈2.6 × BW double hump).**
- [x] **A migration model, so the time-dependent consequence of micromotion can be followed rather than classified at a threshold.**
  **Done (2026-09-30): `MigrationModel` — closed-form logarithmic migration `x(N) = x_bed·ln(1 + k(δ−δ_th)N/x_bed)` (bedding-in decay), verified against numerical integration; `velocity_per_year`/`is_at_risk` implement the RSA-style >0.2 mm/year flag.**
- [x] **Built-in zone definitions (Gruen, Paprosky) so callers are not left to invent one.**
  **Done (2026-09-29): `GruenZone`/`gruen_zone`/`gruen_zones_in_order` — geometric zones 1-7 (Paprosky remains future work; the Gruen half of the item is delivered).**
- [ ] **Continuum coupling, so an implant with realistic compliance can be evaluated rather than modelled as a rigid punch.**

### `tpt-med-stents`

- [ ] ****Level 2** — tapered ring groups, which would give `dogboning` a real value instead of the structural `0.0` it reports today.**
- [ ] ****Level 3** — 3D superelastic FEM with frictional contact via `tpt-fem-hyperelastic` / `tpt-fem-contact`.**
- [x] **Foreshortening, and per-crown stiffness variation for a non-uniform ring.**
  **Done (2026-09-30): `StentModel::foreshortening` — diamond-cell geometry with a `link_fraction` calibration landing in the published few-percent band (NaN beyond the developed-length limit); `simulate_deployment_with_crowns`/`NonUniformDeployment` — per-crown stiffness slice with proportional force split and peak-crown share (empty slice reproduces the uniform ring exactly).**
- [x] **Cyclic degradation of `ε_L` over 10⁶ cycles, to support fatigue and accelerated-dilation life claims.**
  **Done (2026-09-30): `StrainLifeLaw::nitinol_screening` — log-log strain-life law (0.4 % at 10⁷, factor-of-two per four decades per the published fatigue band) with a conservative `survives(N, ε)` verdict. Screening interpolation, not a device S–N curve.**
- [ ] **Direct coupling to a `tpt-med-hemodynamics` solution in the same solve, rather than a prescribed vessel law.**

### `tpt-med-surgical-planning`

- [x] **Per-fragment addressing, so a `PlanStep::Move` can target a single named fragment rather than the whole assembled model. This is the largest known gap and needs an RFC.**
  **Implemented (2026-10-01) as RFC 0011's option-3 first slice, following maintainer direction to proceed: `DiscardedSide::{Resect, RetainAs { name }}` + `PlanStep::MoveNamed`/`move_fragment_named` with build-time `PlanError` validation; plans without retention are byte-identical (rule 4 asserted by the untouched prior suite). The implementation also exposed that the RFC's original option 1 was itself flawed — recorded in the RFC. Remaining (tracked in the crate's CHANGELOG): cuts after named moves (grid unification), multi-fragment retaining cuts, collision handling beyond last-write-wins.**
- [ ] **Curved and freeform resections, saw-kerf width, and multi-plane wedges.**
  **Kerf and wedges delivered (2026-10-01): `OsteotomyCut::kerf_width` (symmetric slab removal, boundary shifted by half the kerf, depth measured past the kept face) and `WedgeCut`/`PlanStep::Wedge` (the exact two-plane intersection sequential cuts cannot express). Curved/freeform resections remain open.**
- [ ] **Implant component placement with a bone–implant interface, and bone graft or defect reconstruction.**
- [ ] **Soft-tissue structures, so a plan can be checked for collateral damage to ligaments, capsules and neurovascular bundles.**
- [x] **Measurement reporting: resection volumes, cut depths and achieved alignment errors, recorded alongside the steps in the audit log.**
  **Done (2026-09-30): `execute_with_report` → `SurgeryReport` with `CutMeasurement` (resection volume, cut depth) and `MoveMeasurement` (prescribed vs achieved centroid displacement = the alignment error, surfacing sub-voxel quantisation), index-aligned with the audit log plus total resection.**

### `tpt-med-tissue`

- [x] **Second-order tangent moduli per model, which the nonlinear `tpt-fem` upgrade path needs for Newton convergence.**
  **Done (2026-09-30): `MaterialTangent` (`A[i][j](k,l) = ∂P_ij/∂F_kl`) via `material_tangent` — analytic for Neo-Hookean and Yeoh (the `q′` term needs a second `β`) and the model-independent `volumetric_tangent`; central differences for Mooney–Rivlin/Ogden/HGO (same split as `first_piola`). Verified by FD agreement, major symmetry `A_ij,kl = A_kl,ij` on all five models (the minor symmetry does not hold), and the volumetric FD cross-check.**
- [x] **Plane-stress and reduced-order wrappers over the full 3×3 `F` interface.**
  **Done (2026-10-01): `ReducedPlaneModel` with `PlaneCondition::{PlaneStrain, PlaneStress}` — plane strain pins `F₃₃ = 1`; plane stress solves `P₃₃ = 0` by bracketed bisection on the scalar `F₃₃` (verified against the incompressible closed form `F₃₃ = 1/det F₂ₓ₂` and traction-freeness).**
- [ ] **Fiber-family rotation in HGO (collagen crimp), and the two-family elastin/collagen parameterisation used in some literature.**

### `tpt-med-viscoelastic`

- [x] **Temperature shifting of relaxation times (WLF and Arrhenius relations), so `τᵢ` values can be generated rather than supplied.**
  **Done (2026-09-29): `TemperatureShift::Wlf/Arrhenius` + `shifted_material` — τᵢ(T) generated from cited C1/C2 or Eₐ.**
- [ ] **Non-linear hyperviscoelastic formulations, applying the Prony series to the hyperelastic energy in finite strain rather than in the linear-viscoelastic regime.**
  **Fung-type QLV delivered (2026-10-01): `QuasiLinearViscoelastic` — the Prony kernel convolved with a hyperelastic stress history (step histories exact; linear-regime equivalence with `PronyIntegrator` asserted). The internal-variable finite-strain formulation remains open as the thermodynamically complete alternative.**
- [x] **A time-integration helper, so callers driving a finite-element inner loop do not each reimplement the recurrence.**
  **Done (2026-09-29): `PronyIntegrator` — exact exponential recurrence; matches G(dt) exactly and the analytic ramp response.**

### `tpt-med-vv40`

- [x] **Per-component credibility rollup, so a large model assembled from small verified parts has a defined composite credibility.**
  **Done (2026-09-30): `AssessmentRollup`/`RollupMember` — conjunction over members, composite rating = max of declared and member ratings (a declaration can only raise), pooled evidence vs composite goals, member-prefixed unmet goals + JSON.**
- [x] **Structured evidence: attach numeric metrics with acceptance criteria to a `VerificationActivity` or `ValidationActivity`, so adequacy can be checked mechanically rather than by reading `results` as prose.**
  **Done (2026-09-29): `EvidenceMetric`/`Acceptance` + `metrics_adequate()` — mechanical acceptance-band checking.**
- [x] **Multi-question assessments with an explicit aggregation rule, for models used for several questions of interest.**
  **Done (2026-09-30): the same `AssessmentRollup` serves this shape — one named entry per question, each judged on its own risk/influence, the study credible only when all are and the pooled evidence meets the composite goals. No majority-vote aggregation, by design.**
- [x] **Optional serialisation of an assessment to JSON, so it can live inside a submission bundle next to the `tpt-med-fda` package.**
  **Done (2026-09-29): `CredibilityAssessment::to_json` — deterministic JSON with goals, activities, verdict and unmet goals.**

### `tpt-med-wear`

- [x] **Wear-debris-induced damage feedback, so wear changes the contact geometry and pressures — without which the runaway that ends real implant life is not captured.**
  **Done (2026-09-30): the feedback loop in `simulate_wear_with_contact` — accumulated per-zone wear depths drive each block's pressure solve. The built-in Winkler foundation is self-stabilising (verified: load migration off the faster-worn zone, coupled total below the prescribed run); runaway modes enter through the same `ContactSolver` trait.**
- [x] **A coupling helper to a contact solver, so pressures and sliding distances can be solved rather than supplied.**
  **Done (2026-09-30): `ContactSolver` trait with the built-in `WinklerContact` (penetration solved so zones carry the total load; worn zones shed load; contact loss flagged). A real contact solver implements the trait.**
- [x] **Uncertainty propagation over the wear coefficient, which scatters over orders of magnitude between studies and which a defensible screening study should quantify.**
  **Done (2026-09-29): `WearUncertainty` + `simulate_wear_uncertainty` — low/central/high band and `relative_band`.**
- [x] **A run-in period and activity-level variation, so gait extrapolation is not strictly linear in cycle count.**
  **Done (2026-09-29): `WearSchedule` + `simulate_wear_with_schedule` — cycle-exact decaying run-in and activity bands.**

**67 items across 20 crates.**
