# tpt-medical

A fully open-source, pure-Rust computational biomechanics and medical device
simulation engine, compiling to WebAssembly for zero-cloud, browser-based
surgical planning.

**Organization:** TPT Solutions (`tpt-solutions`) ·
**License:** MIT OR Apache-2.0 ·
**Crate prefix:** `tpt-med-*`

---

## Why

Medical device design and surgical planning require FDA-mandated FEA and
fluid dynamics, but legacy tools are proprietary black boxes (Abaqus, Ansys),
carry academic license traps (FEBio, OpenSim), or demand cloud upload of
patient DICOM data (HIPAA/GDPR risk). `tpt-medical` eliminates all three:

- **No license barrier** — permissive dual licensing, no copyleft contamination.
- **Zero-cloud patient privacy** — compiles to WASM; imaging and simulation run
  entirely in the browser or on a local hospital workstation.
- **Reproducible FDA submissions** — open, auditable mathematics with
  21 CFR Part 11 audit trails and ASME V&V 40 credibility tooling.

## Crates

Every crate is documented individually: a comprehensive **README** (overview,
features, conventions, runnable usage, API table, verification evidence,
dependencies, known limitations) and a **Keep a Changelog**-format
**CHANGELOG** recording what changes for its consumers. Crates version and
release independently, so those two files are the per-crate source of truth.

| Crate | Description | Docs | Status |
|---|---|---|---|
| [`tpt-med-core`](crates/core/tpt-med-core) | Patient models, anatomical coordinate systems, HIPAA-safe audit traits | [README](crates/core/tpt-med-core/README.md) · [CHANGELOG](crates/core/tpt-med-core/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-units`](crates/core/tpt-med-units) | Type-safe unit system (mm, MPa, N, g/cm³) | [README](crates/core/tpt-med-units/README.md) · [CHANGELOG](crates/core/tpt-med-units/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-geometry`](crates/core/tpt-med-geometry) | Geometry primitives, anatomical transforms | [README](crates/core/tpt-med-geometry/README.md) · [CHANGELOG](crates/core/tpt-med-geometry/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-wasm`](crates/core/tpt-med-wasm) | WASM bindings for the in-browser stack | [README](crates/core/tpt-med-wasm/README.md) · [CHANGELOG](crates/core/tpt-med-wasm/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-dicom`](crates/imaging/tpt-med-dicom) | DICOM parsing & Hounsfield Unit mapping | [README](crates/imaging/tpt-med-dicom/README.md) · [CHANGELOG](crates/imaging/tpt-med-dicom/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-meshing`](crates/imaging/tpt-med-meshing) | Segmentation masks & voxel-to-hex meshing | [README](crates/imaging/tpt-med-meshing/README.md) · [CHANGELOG](crates/imaging/tpt-med-meshing/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-biomechanics`](crates/solid/tpt-med-biomechanics) | Voxel hexahedral FEM solver core | [README](crates/solid/tpt-med-biomechanics/README.md) · [CHANGELOG](crates/solid/tpt-med-biomechanics/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-tissue`](crates/solid/tpt-med-tissue) | Hyperelastic tissue models (Neo-Hookean, Mooney-Rivlin, Ogden, HGO) | [README](crates/solid/tpt-med-tissue/README.md) · [CHANGELOG](crates/solid/tpt-med-tissue/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-bone`](crates/solid/tpt-med-bone) | Bone mechanics & Wolff's-Law remodeling | [README](crates/solid/tpt-med-bone/README.md) · [CHANGELOG](crates/solid/tpt-med-bone/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-viscoelastic`](crates/solid/tpt-med-viscoelastic) | Prony-series viscoelasticity | [README](crates/solid/tpt-med-viscoelastic/README.md) · [CHANGELOG](crates/solid/tpt-med-viscoelastic/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-cartilage`](crates/solid/tpt-med-cartilage) | Biphasic/poroelastic cartilage models | [README](crates/solid/tpt-med-cartilage/README.md) · [CHANGELOG](crates/solid/tpt-med-cartilage/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-hemodynamics`](crates/fluid/tpt-med-hemodynamics) | Navier–Stokes CFD, WSS & OSI | [README](crates/fluid/tpt-med-hemodynamics/README.md) · [CHANGELOG](crates/fluid/tpt-med-hemodynamics/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-cardiovascular`](crates/fluid/tpt-med-cardiovascular) | Windkessel models, FFR | [README](crates/fluid/tpt-med-cardiovascular/README.md) · [CHANGELOG](crates/fluid/tpt-med-cardiovascular/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-stents`](crates/devices/tpt-med-stents) | Nitinol superelasticity, stent deployment | [README](crates/devices/tpt-med-stents/README.md) · [CHANGELOG](crates/devices/tpt-med-stents/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-orthopedics`](crates/devices/tpt-med-orthopedics) | Micromotion & stress-shielding analysis | [README](crates/devices/tpt-med-orthopedics/README.md) · [CHANGELOG](crates/devices/tpt-med-orthopedics/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-wear`](crates/devices/tpt-med-wear) | Archard / Cross-Land wear laws | [README](crates/devices/tpt-med-wear/README.md) · [CHANGELOG](crates/devices/tpt-med-wear/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-surgical-planning`](crates/surgical/tpt-med-surgical-planning) | Osteotomy cuts, virtual surgery | [README](crates/surgical/tpt-med-surgical-planning/README.md) · [CHANGELOG](crates/surgical/tpt-med-surgical-planning/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-implant-sizing`](crates/surgical/tpt-med-implant-sizing) | Automated implant sizing from anatomy | [README](crates/surgical/tpt-med-implant-sizing/README.md) · [CHANGELOG](crates/surgical/tpt-med-implant-sizing/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-audit`](crates/regulatory/tpt-med-audit) | Cryptographic audit signing & verification | [README](crates/regulatory/tpt-med-audit/README.md) · [CHANGELOG](crates/regulatory/tpt-med-audit/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-fda`](crates/regulatory/tpt-med-fda) | 21 CFR Part 11 audit trails | [README](crates/regulatory/tpt-med-fda/README.md) · [CHANGELOG](crates/regulatory/tpt-med-fda/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-vv40`](crates/regulatory/tpt-med-vv40) | ASME V&V 40 credibility assessment | [README](crates/regulatory/tpt-med-vv40/README.md) · [CHANGELOG](crates/regulatory/tpt-med-vv40/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-examples`](examples) | Runnable end-to-end example binaries | [README](examples/README.md) · [CHANGELOG](examples/CHANGELOG.md) | 🚧 Alpha |
| [`tpt-med-benches`](benches) | Workspace benchmark suite | [README](benches/README.md) · [CHANGELOG](benches/CHANGELOG.md) | 🚧 Alpha |

## Quick Start

Build the workspace and run the Phase 1 milestone CLI, which converts a DICOM
CT series into a patient-specific bone mesh CSV:

```console
cargo run -p tpt-med-examples --bin dicom-to-mesh -- \
    test-data/dicom/synthetic_ct --output femur_mesh.csv --threshold 200
```

Use the solver stack in Rust:

```rust
use std::path::Path;

use tpt_med_dicom::DicomSeries;
use tpt_med_meshing::{MedicalMesher, SegmentationMask};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let series = DicomSeries::load_from_dir(Path::new("test-data/dicom/synthetic_ct"))?;
    let mask = SegmentationMask::threshold_hu(&series, 200.0);
    let mesh = MedicalMesher::default().voxels_to_hex_mesh(&mask)?;
    mesh.write_csv(Path::new("bone_mesh.csv"))?;
    Ok(())
}
```

## Zero-Cloud Patient Privacy

All simulations run entirely locally; the WASM build (`tpt-med-wasm`) runs the
full imaging-to-simulation pipeline inside the browser. Patient DICOM scans
never leave the hospital network. The repository's test data is **synthetic
only** — never commit real patient data.

## Browser Demos

Build the engine glue once, then serve the repository root:

```console
./scripts/build-web.sh          # or scripts\build-web.ps1 on Windows
python -m http.server 8080
```

- **Viewer** — `http://localhost:8080/web/viewer/`: upload a DICOM series
  (or the synthetic CT demo), mesh it, solve a stance load, and inspect
  von Mises stress in WebGL2; plus in-browser stent deployment.
- **White-label stent simulator** — `http://localhost:8080/web/stent-simulator/`:
  the embeddable `<tpt-stent-simulator>` custom element with live branding
  controls (see [`web/stent-simulator/README.md`](web/stent-simulator/README.md)).

> **Try it without installing anything:** a live build of the viewer is
> published to GitHub Pages on every merge to `master` that touches the engine
> or the web assets, against the committed **synthetic** CT data. See the
> repository's *Demo* workflow for the published URL, or run it locally with
> the two commands above.

## Verification & Validation

Following ASME V&V 40, every solver ships with verification tests against
analytical solutions (uniaxial tension, Euler–Bernoulli beam deflection,
Poiseuille flow, Prony relaxation) and curated golden reference datasets under
[`test-data/golden/`](test-data/golden). Standards automation targets ASTM
F2028 / F2079 / F2394 and ISO 7206 / 14879.

## Substrate

Built on the TPT engineering stack (pinned in the root `Cargo.toml`):
`tpt-math` (fixed-size linalg, optimization, probability, signal filtering),
`tpt-engineering` (biocompatible material libraries), `tpt-science`
(bio-fluid dynamics), and `tpt-fem` (nonlinear FEM, contact). Medical crates
are deliberately dependency-light (`std` only) to keep the WASM footprint
small; substrate integration points are documented per crate and in the
RFCs.

## Contributing

Bug reports, feature requests, RFC proposals, and documentation corrections
are all welcome, and all arrive as **GitHub issues** — this project does not
accept external pull requests. See [CONTRIBUTING.md](CONTRIBUTING.md) for the
process and the [RFC process](rfcs) for significant design changes. New
constitutive models, solver algorithms, and regulatory features require an
accepted RFC, and numerical behaviour is treated as API. Everything is
DCO-signed and CLA-free.

## License

Licensed under either of

- MIT license ([LICENSE-MIT](LICENSE-MIT))
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))

at your option.

## Regulatory Disclaimer

This software is for research and development purposes only. It is **not**
cleared or approved by the FDA or any other regulatory body for clinical
diagnostic or treatment use. Regulatory submissions require validation under
ASME V&V 40 and applicable standards.
