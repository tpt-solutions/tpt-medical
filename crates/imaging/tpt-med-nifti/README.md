# tpt-med-nifti

Pure-Rust NIfTI-1 volume parsing for research-space imaging. Zero external
crates by default, no C bindings, no network; the optional `gzip` feature
adds one pure-Rust dependency (`flate2`).

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--nifti-orange)](https://crates.io/crates/tpt-med-nifti)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--nifti-blue)](https://docs.rs/tpt-med-nifti)

| | |
|---|---|
| **Layer** | `imaging` |
| **Status** | Alpha, `0.1.0` |
| **Scope** | RFC 0006 — NIfTI ingestion |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.84 |
| **Dependencies** | [`tpt-med-geometry`](../../core/tpt-med-geometry); optional [`flate2`](https://crates.io/crates/flate2) (`gzip` feature, pure-Rust backend) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

`tpt-med-dicom` covers DICOM, which is what hospital archives export — but
it is not what research tooling uses once a scan leaves the archive.
`dcm2niix`, FSL, FreeSurfer, ANTs, and most public imaging datasets work in
NIfTI. Without this crate, using one of those volumes in this pipeline means
converting it back to DICOM with an external tool first. Three reasons this
is its own crate rather than a feature of `tpt-med-dicom`:

1. **Different format, different failure modes.** NIfTI's fixed 348-byte
   header, sform/qform affine, and `scl_slope`/`scl_inter` scaling have
   nothing in common with DICOM's tag/VR system — sharing a crate would mean
   sharing nothing but the name.
2. **No PHI, less at stake for privacy, but the honesty-about-failure
   principle still applies.** A NIfTI header that this crate cannot decode
   confidently is rejected with a typed error, never silently mis-parsed.
3. **Same dependency-light reasoning as RFC 0001.** A general-purpose `nifti`
   crate exists on crates.io; the same weight/WASM-footprint/audit-surface
   argument that kept `tpt-med-dicom` dependency-free applies here, and the
   format is if anything simpler to hand-roll than DICOM (one fixed-size
   header, no tag/VR system).

## Features

- **NIfTI-1, both file layouts.** Single-file `.nii` (magic `n+1`) and
  dual-file `.hdr`/`.img` (magic `ni1`, via `parse_dual_file` /
  `parse_dual_bytes`, with `vox_offset` honoured as an offset into the
  `.img`). The 348-byte header (`nifti1.h`, unchanged since 2004) is
  decoded directly; endianness is detected from `sizeof_hdr` (NIfTI has no
  separate byte-order flag).
- **`.nii.gz` behind the off-by-default `gzip` feature** — gzip is by far
  the most common form in practice. Decompression is a *streaming,
  geometry-capped* inflate: the header's own declared extent bounds what
  the reader will ever allocate (the same posture `tpt-med-dicom`'s RLE
  decoder takes against an expanding run), via `flate2`'s pure-Rust
  backend (`miniz_oxide`, no C in the build). Without the feature a
  gzipped input — single-file or either dual-file part — is named with
  `NiftiError::Gzipped`, never misparsed.
- **sform and qform geometry**, with the spec's own precedence: sform when
  `sform_code > 0`, else qform (quaternion → rotation, per `nifti1.h`'s
  documented formula) when `qform_code > 0`, else an axis-aligned
  spacing-only fallback ("Analyze-compatible" mode).
- **8 datatypes**: `uint8`, `int8`, `int16`, `uint16`, `int32`, `uint32`,
  `float32`, `float64` — every datatype a CT or MR NIfTI export realistically
  uses. Others are rejected by name, not guessed at.
- **Value scaling** — `value = raw * scl_slope + scl_inter` when
  `scl_slope != 0`, mirroring DICOM's `RescaleSlope`/`RescaleIntercept`
  convention (including "zero slope means no scaling").
- **RAS geometry**, interoperable with `tpt-med-dicom`'s LPS convention via
  `tpt_med_geometry::ras_to_lps`.
- **Synthetic volume generation** (`SyntheticNiftiBuilder`) — writes real,
  parseable `.nii` bytes for both sform and qform geometry, so the test
  suite needs no real dataset (there is no PHI risk for NIfTI the way there
  is for DICOM, but a hand-verified fixture is still how format-decode
  correctness gets checked here, same as `tpt-med-dicom`).

## Explicit Non-Features

- **No acquisition metadata** (`descrip`, `aux_file`, `intent_*`,
  slice-timing) — geometry and scaled values only, per RFC 0006's scope.

## Conventions

- Values are stored as `f64`, already scaled by `scl_slope`/`scl_inter`. The
  original on-disk type is kept on `NiftiVolume::datatype` for reference.
- Geometry is **RAS** (NIfTI's native convention): `+x` Right, `+y` Anterior,
  `+z` Superior — the opposite sign convention from DICOM's LPS.
- Voxel spacing and rotation are kept separate (`voxel_spacing`, `rotation`
  with unit columns), matching `tpt_med_dicom::ImageFrame`'s
  direction-cosines-plus-spacing shape rather than a single baked-in affine
  — even though NIfTI's own sform *is* a single baked-in affine, decomposed
  back into the two on the way in.
- This crate assumes no shear in the sform affine (each axis decomposed as a
  spacing magnitude times a unit direction). A sheared sform is a real gap,
  not something silently handled — see Known Limitations.

## Usage

### Loading a volume

```rust
use tpt_med_nifti::NiftiVolume;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let vol = NiftiVolume::parse_file(std::path::Path::new("scan.nii"))?;

    let (nx, ny, nz) = vol.dims;
    println!("{nx} x {ny} x {nz}, spacing {:?}", vol.voxel_spacing);

    if let Some(v) = vol.value_at(nx / 2, ny / 2, nz / 2) {
        println!("centre value = {v}");
    }

    // RAS mm position of that same voxel.
    let p = vol.voxel_position(nx / 2, ny / 2, nz / 2);
    println!("centre position (RAS mm) = {p:?}");

    Ok(())
}
```

### Gzipped and dual-file inputs

With the `gzip` feature enabled, `.nii.gz` reads through the same entry
points; without it the input is named with `NiftiError::Gzipped`. The
dual-file layout is two files that must be handed over together:

```rust
use tpt_med_nifti::NiftiVolume;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // `gzip` feature on: compressed single-file, same call as plain .nii.
    let gz = NiftiVolume::parse_file(std::path::Path::new("scan.nii.gz"))?;

    // Dual-file pair: header file + voxel-data file (either may be
    // gzipped when the `gzip` feature is on).
    let pair = NiftiVolume::parse_dual_file(
        std::path::Path::new("scan.hdr"),
        std::path::Path::new("scan.img"),
    )?;
    println!("{} x {} x {} (gz), {} voxels (dual)",
             gz.dims.0, gz.dims.1, gz.dims.2,
             pair.dims.0 * pair.dims.1 * pair.dims.2);
    Ok(())
}
```

### A CT-derived NIfTI into the existing HU pipeline

If the volume is a CT export whose writer baked DICOM's
`RescaleSlope`/`RescaleIntercept` into `scl_slope`/`scl_inter` (as `dcm2niix`
does), `values` are already HU and can go straight into
`tpt_med_dicom::HounsfieldMapper` — which takes a bare `f64`, with no
DICOM-specific coupling:

```rust
use tpt_med_dicom::{BoneRegion, HounsfieldMapper};
use tpt_med_nifti::NiftiVolume;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let vol = NiftiVolume::parse_file(std::path::Path::new("ct_export.nii"))?;
    if let Some(hu) = vol.value_at(0, 0, 0) {
        let e = HounsfieldMapper::hu_to_youngs_modulus(hu, BoneRegion::Cortical);
        println!("modulus at voxel (0,0,0): {} MPa", e.to_mpa());
    }
    Ok(())
}
```

This crate cannot verify that assumption for you — NIfTI carries no
`Modality`-equivalent tag saying a value *is* HU.

## API Overview

| Item | Purpose |
|---|---|
| `NiftiVolume::{parse_file, parse_bytes}` | Parse a single-file NIfTI volume, version 1 or 2 (`.nii`, or `.nii.gz` with the `gzip` feature) |
| `NiftiVolume::{parse_dual_file, parse_dual_bytes}` | Parse a dual-file `.hdr`/`.img` pair, either version (either part gzippable with `gzip`) |
| `NiftiVolume::version` -> `NiftiVersion` | Which format version (`V1`/`V2`) the volume was parsed from |
| `NiftiVolume::{dims, voxel_spacing, origin, rotation, datatype, values}` | Parsed geometry and data |
| `NiftiVolume::value_at(i, j, k)` | Bounds-checked value access |
| `NiftiVolume::voxel_position(i, j, k)` | RAS-mm position of a voxel centre |
| `Datatype` | The 8 supported on-disk voxel types |
| `SyntheticNiftiBuilder` | Builder for synthetic fixtures (`build` → `.nii`, `build_dual` → `.hdr`/`.img`; sform or qform, u8 or f32) |
| `NiftiError` | `Io`, `NotNifti`, `Gzipped`, `CorruptGzip`, `UnsupportedDatatype`, `BitpixMismatch`, `UnexpectedEof`, `BadValue` |
| `Result<T>` | Crate result alias |

## Verification

- The qform quaternion → rotation formula is checked against a hand-derived
  closed-form case (180° and 90° rotations about the k-axis, worked by hand
  for RFC 0006's review, not copied from an external tool) and a property
  test that the resulting matrix is orthonormal with determinant `+1` across
  a spread of unit quaternions — catching an algebra transcription error
  independent of any single worked example.
- The header/voxel-data decode is checked by round-tripping
  `SyntheticNiftiBuilder` through `parse_bytes` for both sform and qform
  geometry, `u8` and `f32` datatypes, scaled and unscaled values, and the
  axis-aligned fallback (both form codes absent).
- The dual-file path is round-tripped against the single-file one on
  identical geometry and data, plus a non-zero `vox_offset` into the
  `.img`; with `gzip` on, a compressed volume (and a gzipped `.hdr`/`.img`
  pair) is asserted to parse to the same result as the plain bytes it
  wraps.
- Malformed input is asserted to be rejected, not padded or guessed at:
  truncated header, truncated voxel data, an unsupported datatype code, a
  `bitpix`/`datatype` mismatch, a short dual-file `.img`, the wrong magic
  at either entry point (with a pointer to the right one), gzip magic
  bytes without the feature, a corrupt gzip stream and a gzip of a
  truncated volume with it — each has a dedicated test.
- **Not verified**: behaviour against a real-world writer's output
  (`dcm2niix`, FSL, SPM, ITK-SNAP). Only hand-built and self-round-tripped
  fixtures are exercised — see `rfcs/0006-nifti-ingestion.md`'s "What remains
  unverified".

## Known Limitations

- **No NIfTI-2, and no gzip without the `gzip` feature.** See Explicit
  Non-Features.
- **No shear in sform geometry.** Spacing is recovered as each affine
  column's norm and rotation as that column normalised; a genuinely sheared
  sform (non-orthogonal columns) is decomposed incorrectly rather than
  detected and rejected. No real-world writer this crate has been checked
  against produces one, but that is not the same as a guarantee.
- **No qform-quaternion validation.** A non-unit `(quatern_b, quatern_c,
  quatern_d)` (i.e. `b²+c²+d² > 1`) is clamped rather than rejected, per
  `nifti1.h`'s own formula — see `rfcs/0006-nifti-ingestion.md` Drawbacks.
- **No acquisition metadata.** `descrip`, `aux_file`, `cal_min`/`cal_max`,
  `intent_*`, and slice timing are not decoded — this crate answers "what
  are the values and where are they in space," nothing about acquisition.

## Related Crates

- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) — the DICOM equivalent; `HounsfieldMapper` and `QctCalibration` take a bare `f64`, usable with a NIfTI-derived value with no coupling to this crate.
- [`tpt-med-meshing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-meshing) — `SegmentationMask::threshold_nifti` turns a `NiftiVolume` into a bone mask in the same LPS patient frame the DICOM path produces.
- [`tpt-med-geometry`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-geometry) — `Vec3`, `Mat3`, and `ras_to_lps` for interop with DICOM's patient frame.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Remaining
format support (NIfTI-2) follows the roadmap in
`rfcs/0006-nifti-ingestion.md`.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. Not a
diagnostic medical device.
