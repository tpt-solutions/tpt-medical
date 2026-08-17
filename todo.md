# tpt-medical — TODO

**Organization:** TPT Solutions (`tpt-solutions`)
**License:** MIT OR Apache-2.0 (dual)
**Description:** Fully open-source, pure-Rust computational biomechanics and medical device simulation engine, compiling to WASM for zero-cloud, browser-based surgical planning.

---

## Phase 0: Repo Scaffolding & Substrate Status

### Workspace & Tooling
- [ ] Root `Cargo.toml` with `[workspace]` + `[workspace.package]` (edition 2021, rust-version, license = "MIT OR Apache-2.0", authors = ["TPT Solutions"], repository)
- [ ] `LICENSE-MIT`
- [ ] `LICENSE-APACHE`
- [ ] `deny.toml` (allow MIT, Apache-2.0, BSD-2-Clause, BSD-3-Clause, ISC, Zlib, Unicode-3.0; deny copyleft/unlicensed)
- [ ] `rustfmt.toml`
- [ ] `clippy.toml`

### Docs & Community Files
- [ ] `README.md`
- [ ] `CONTRIBUTING.md`
- [ ] `SECURITY.md`
- [ ] `CODE_OF_CONDUCT.md`
- [ ] `CHANGELOG.md`

### CI/CD
- [ ] `.github/workflows/ci.yml` (fmt, clippy -D warnings, test, build)
- [ ] `.github/workflows/license.yml` (cargo deny check licenses)
- [ ] `.github/workflows/benchmark.yml`
- [ ] `.github/workflows/docs.yml`
- [ ] `.github/workflows/release.yml`
- [ ] `.github/ISSUE_TEMPLATE/bug_report.md`
- [ ] `.github/ISSUE_TEMPLATE/feature_request.md`
- [ ] `.github/ISSUE_TEMPLATE/rfc.md`
- [ ] `.github/PULL_REQUEST_TEMPLATE.md`

### Directory Layout
- [ ] `crates/core/`, `crates/imaging/`, `crates/solid/`, `crates/fluid/`, `crates/devices/`, `crates/surgical/`, `crates/regulatory/`
- [ ] `examples/`
- [ ] `test-data/dicom/`, `test-data/nifti/`, `test-data/meshes/`, `test-data/golden/`
- [ ] `benches/`
- [ ] `docs/book/`, `docs/rfc/`, `docs/api/`
- [ ] `rfcs/`

### Substrate Dependency Status (external GitHub repos)
- [ ] Confirm/pin `tpt-math` version (linalg-fixed, optimize-general, prob-dist, signal-filter)
- [ ] Confirm/pin `tpt-engineering` version (material libraries: PEEK, Nitinol, Titanium, UHMWPE)
- [ ] Confirm/pin `tpt-science` version — **in progress upstream**; blocks Phase 4 (Hemodynamics/CFD)
- [ ] Confirm/pin `tpt-fem` version — **in progress upstream**; blocks Phase 2+ (FEM solver core)
- [ ] Re-check substrate maturity before starting each blocked phase below

---

## Phase 1: Foundation
*Months 1-3*
- [ ] `tpt-med-core` — patient models, anatomical coordinate systems, HIPAA-safe audit traits
- [ ] `tpt-med-units` — unit system
- [ ] `tpt-med-geometry` — geometry primitives
- [ ] `tpt-med-dicom` — DICOM parsing, Hounsfield Unit (HU) mapping
- [ ] `tpt-med-meshing` — voxel-to-hex meshing
- [ ] **Milestone:** CLI tool that takes DICOM CT scan, outputs bone mesh CSV

## Phase 2: Solid Biomechanics
*Months 4-6 — ⚠ blocked on `tpt-fem` maturity*
- [ ] `tpt-med-biomechanics` — FEM solver core
- [ ] `tpt-med-bone` — linear elastic bone mechanics
- [ ] `tpt-med-tissue` — Neo-Hookean, Mooney-Rivlin models
- [ ] **Milestone:** Femur stress analysis under physiological loading

## Phase 3: Advanced Tissue Models
*Months 7-9*
- [ ] `tpt-med-tissue` — Holzapfel-Gasser-Ogden (HGO) model for arteries
- [ ] `tpt-med-viscoelastic` — Prony series relaxation
- [ ] `tpt-med-cartilage` — biphasic/poroelastic models
- [ ] **Milestone:** Arterial wall inflation simulation

## Phase 4: Hemodynamics
*Months 10-12 — ⚠ blocked on `tpt-science` maturity*
- [ ] `tpt-med-hemodynamics` — Navier-Stokes CFD, Wall Shear Stress, Oscillatory Shear Index
- [ ] `tpt-med-cardiovascular` — Windkessel model, Fractional Flow Reserve (FFR)
- [ ] **Milestone:** Carotid bifurcation Wall Shear Stress calculation

## Phase 5: Implants & Devices
*Months 13-15*
- [ ] `tpt-med-stents` — Nitinol superelasticity, crimping/expansion deployment
- [ ] `tpt-med-orthopedics` — micromotion analysis, stress shielding
- [ ] `tpt-med-wear` — Archard/Cross-Land wear laws
- [ ] **Milestone:** Stent deployment simulation with artery contact

## Phase 6: Surgical Planning
*Months 16-18*
- [ ] `tpt-med-surgical-planning` — osteotomy cuts, virtual surgery
- [ ] `tpt-med-implant-sizing` — automated sizing from anatomy
- [ ] **Milestone:** Virtual total knee replacement planning

## Phase 7: Regulatory & Compliance
*Months 19-21*
- [ ] `tpt-med-fda` — 21 CFR Part 11 immutable, cryptographically signed audit trails
- [ ] `tpt-med-vv40` — ASME V&V 40 credibility assessment matrices
- [ ] `tpt-med-audit` — cryptographic signature verification
- [ ] **Milestone:** Export FDA-submission-ready simulation package

## Phase 8: WASM & Ecosystem
*Months 22-24*
- [ ] `tpt-med-wasm` — browser compilation of full simulation stack
- [ ] Web-based surgical planning viewer (WebGL/WebGPU, no VTK/OpenGL)
- [ ] **Milestone:** Zero-cloud, in-browser patient-specific simulation

---

## Cross-Cutting / Ongoing

### Verification & Validation (ASME V&V 40)
- [ ] `test-data/golden/solid/` — femur_loading, lumbar_spine_compression, arterial_wall_inflation
- [ ] `test-data/golden/fluid/` — carotid_bifurcation_cfd, aortic_aneurysm_flow, coronary_ffr
- [ ] `test-data/golden/devices/` — stent_expansion, hip_stem_micromotion, knee_wear_10mcycles
- [ ] `test-data/golden/regulatory/` — fda_audit_trail, vv40_credibility_matrix
- [ ] Verification test: uniaxial tension vs. analytical Neo-Hookean solution
- [ ] Verification test: Poiseuille flow vs. analytical CFD solution
- [ ] Verification test: stent radial stiffness vs. ASTM F2394 published data

### ASTM / ISO Standards Automation
- [ ] ASTM F2028 — dynamic evaluation of total knee replacements
- [ ] ASTM F2079 — intrinsic securement of endovascular stents (intramedullary rods)
- [ ] ASTM F2394 — securement of self-expanding stents
- [ ] ISO 7206 — hip joint prostheses
- [ ] ISO 14879 — wear of total knee-replacement prostheses

### RFCs
- [ ] `rfcs/0001-dicom-ingestion.md`
- [ ] `rfcs/0002-hyperelastic-tissue.md`
- [ ] `rfcs/0003-fda-audit-trail.md`
- [ ] `rfcs/0004-nitinol-superelasticity.md`
- [ ] `rfcs/0005-cardiac-electrophysiology.md`

### Business Wedge (stretch, ties to Phase 5 + Phase 8 completion)
- [ ] White-label web component: web-based stent deployment simulator for MedTech companies
- [ ] Real-time WASM FEM solver demo: CT scan upload → in-browser stent expansion in patient-specific artery
