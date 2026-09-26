# tpt-med-dicom

Pure-Rust DICOM parsing and Hounsfield Unit mapping for CT-based
patient-specific modelling. No `dicom-rs`, no C bindings, no network.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--dicom-orange)](https://crates.io/crates/tpt-med-dicom)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--dicom-blue)](https://docs.rs/tpt-med-dicom)

| | |
|---|---|
| **Layer** | `imaging` |
| **Status** | Alpha, `0.1.0` |
| **Scope** | RFC 0001 — DICOM ingestion |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-geometry`](../../core/tpt-med-geometry), [`tpt-med-units`](../../core/tpt-med-units) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

The entire patient-specific pipeline starts by reading a CT series. The two
viable Rust options are a heavyweight binding to `dcmtk` (build fragility,
platform matrix, C toolchain) or a pure-Rust parser. We take the pure-Rust
route for three reasons:

1. **Zero-cloud privacy.** A WASM build of a C-binding parser is not a thing.
2. **WASM footprint.** The whole imaging → mesh → solve stack has to fit in a
   ~190 KB `.wasm`; a DICOM library dominates that budget instantly.
3. **Honesty about failure.** Compressed transfer syntaxes are *rejected with a
   typed error*, never silently mis-parsed. A wrong stress field that looks
   plausible is worse than a failed load.

## Features

- **Transfer syntaxes:** implicit VR little endian (`1.2.840.10008.1.2`) and
  explicit VR little endian (`1.2.840.10008.1.2.1`) — the two uncompressed
  syntaxes carried by the overwhelming majority of archive exports.
- **Single-frame CT/MR slices.** Sequences are parsed and skipped safely,
  including undefined-length sequences with depth tracking.
- **Series assembly** — `DicomSeries::load_from_dir` sorts slices by projected
  position along the slice normal, which is the only robust ordering for
  oblique acquisitions; `from_slices` accepts an already-ordered vector.
- **Hounsfield Unit mapping** — `stored × RescaleSlope + RescaleIntercept`, then
  HU → apparent density → Young's modulus through published power laws.
- **Synthetic CT generation** — `SyntheticCtBuilder` and `femur_phantom` write
  real, parseable DICOM files, so the whole test suite runs without PHI.
- **Element encoder** (`encode_element_explicit`) so tests and the synthetic
  generator produce byte-accurate output.

## Explicit Non-Features

Compressed/encapsulated pixel data (JPEG, JPEG-LS, JPEG 2000, RLE) is **not
decoded**. You get `DicomError::CompressedPixelData`, not garbage. Decompress
at the archive boundary — that is a deliberate architectural line, not an
oversight. Multi-frame objects, private tags with odd VRs, and DICOM
networking (C-STORE, DICOMweb) are likewise out of scope for v0.

## Conventions

- Pixels are stored as `i32` **raw stored values**; HUs are derived per slice.
- Geometry follows the DICOM patient coordinate system (**LPS**), exposed via
  `ImageFrame` from tags (0020,0032) and (0020,0037).
- `HU = stored × RescaleSlope + RescaleIntercept`.
- Apparent density uses the linear CT approximation: `ρ = (HU + 1000) / 1000`
  g/cm³, so water is 1.0 and air is 0.0. Quantitative-CT phantom calibration
  replaces this in regulated pipelines.
- Modulus power laws (Morgan–Keaveny style):
  - cortical: `E = 10500 · ρ^2.0`
  - trabecular: `E = 6850 · ρ^1.49`
  - non-positive densities clamp to zero modulus (air voxels).
- `DEFAULT_BONE_THRESHOLD_HU = 200.0` — the literature band is 130–300 HU;
  200 is a common middle choice for appendicular CT. Tune per protocol.



## Usage

### Loading a CT series

```rust
use std::path::Path;
use tpt_med_dicom::{BoneRegion, DicomSeries, HounsfieldMapper};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let series = DicomSeries::load_from_dir(Path::new("test-data/dicom/synthetic_ct"))?;

    let (nx, ny, nz) = series.dims();
    println!("{nx} x {ny} x {nz}, spacing {:?}", series.pixel_spacing);

    // Hounsfield Units via rescale slope/intercept.
    if let Some(hu) = series.hu_at(nx / 2, ny / 2, nz / 2) {
        println!("centre HU = {hu}");
    }

    // HU -> density -> Young's modulus.
    let density = HounsfieldMapper::hu_to_density(1000.0);
    assert!((density.to_gcm3() - 2.0).abs() < 1e-12);

    let e = HounsfieldMapper::hu_to_youngs_modulus(1000.0, BoneRegion::Cortical);
    assert!((e.to_mpa() - 10_500.0 * 4.0).abs() < 1.0);

    // Air must never produce a negative modulus.
    let air = HounsfieldMapper::hu_to_youngs_modulus(-1000.0, BoneRegion::Trabecular);
    assert_eq!(air.to_mpa(), 0.0);

    Ok(())
}
```

### Generating synthetic CT (no PHI, ever)

```rust
use tpt_med_dicom::SyntheticCtBuilder;

fn main() -> std::io::Result<()> {

## API Overview

| Item | Purpose |
|---|---|
| `DicomSeries::load_from_dir(&Path)` | Load and order a series from a directory (slices are sorted by projection along the slice normal) |
| `DicomSeries::from_slices(Vec<DicomSlice>)` | Assemble from pre-ordered slices; validates consistency |
| `DicomSeries::{dims, hu_volume, hu_at}` | `(nx, ny, nz)`; flat HU volume; bounds-checked HU access |
| `DicomSeries::{pixel_spacing, slices, modality}` | Geometry and per-slice metadata |
| `DicomSlice::{hu_at, hu_plane, normal, frame}` | Per-slice HU, plane normal, and DICOM `ImageFrame` |
| `DicomParser::{parse_file, parse_bytes}` | Parse one Part-10 file into a `DicomSlice` |
| `DicomParser::encode_element_explicit(tag, vr, value)` | Byte-accurate explicit-VR element encoder |
| `DicomElement::{as_us, as_is, as_ds_first, as_ds_vec, as_text}` | Typed accessors; each returns a typed `DicomError` on mismatch |
| `TransferSyntax::from_uid`, `::implicit_vr(tag)` | Transfer syntax detection and implicit-VR tag lookup |
| `Vr::from_ascii`, `::ascii`, `::uses_long_length` | VR handling (long-form value lengths) |
| `HounsfieldMapper::hu_to_density` | HU → apparent density (g/cm³) |
| `HounsfieldMapper::density_to_youngs_modulus` | Density → modulus via the cortical or trabecular law |
| `HounsfieldMapper::hu_to_youngs_modulus` | Combined convenience |
| `HounsfieldMapper::DEFAULT_BONE_THRESHOLD_HU` | 200.0 HU |
| `BoneRegion::{Cortical, Trabecular}` | Which correlation applies |
| `SyntheticCtBuilder` | Builder for synthetic series (dims, spacing, thickness, origin, patient id, arbitrary HU function) |
| `SyntheticCtSeries::{parse, write_to_dir}` | Round-trip: write real DICOM files, or parse them back |
| `SyntheticCtBuilder::femur_phantom` | Ready-made femoral shaft phantom |
| `DicomError` | `Io`, `NotDicom`, `UnexpectedEof`, `CompressedPixelData`, `UnknownTransferSyntax`, `UnsupportedVr`, `BadValue`, `InconsistentSeries` |
| `Result<T>` | Crate result alias |

## Verification

- Both transfer syntaxes parse the committed synthetic series
  (`test-data/dicom/synthetic_ct/`, 24 slices) to identical HU volumes — a
  direct test that the implicit and explicit VR paths agree.
- Rescale slope/intercept handling is pinned, including negative intercepts.
- Density and modulus correlations are locked to the exact power laws above,
  with a dedicated air-voxel clamping test.
- The synthetic generator round-trips: `write_to_dir` then `load_from_dir`
  reproduces the HU volume exactly, which is what makes the rest of the
  workspace testable without real patient data.
- `CompressedPixelData` is asserted to be returned, never swallowed.

## Related Crates

- [`tpt-med-meshing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-meshing) — thresholds the HU volume and builds the hex mesh.
- [`tpt-med-geometry`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-geometry) — `ImageFrame` and `Vec3` for the slice normal.
- [`tpt-med-units`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-units) — the `Density` and `Modulus` results.
- [`tpt-med-bone`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-bone) — builds on `BoneRegion` and the modulus correlations.
- [`tpt-med-wasm`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-wasm) — runs this parser in the browser on raw DICOM bytes.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). New tags,
transfer syntaxes or sequences follow the roadmap in
`rfcs/0001-dicom-ingestion.md`. **Never commit real patient data** — all
fixtures are synthetic.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. Not a
diagnostic medical device; Hounsfield Unit output is for research modelling.

    SyntheticCtBuilder::femur_phantom(64, 64, 40)
        .patient_id("SYNTHETIC-0001")
        .build("1.2.826.0.1.3680043.9.7484.1.1")
        .write_to_dir(std::path::Path::new("test-data/dicom/synthetic_ct"))
}
```
