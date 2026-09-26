# Changelog

All notable changes to `tpt-med-nifti` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently**
on a six-week cadence. Workspace-level and cross-cutting changes are
recorded in the [workspace CHANGELOG](../../../CHANGELOG.md); this file
records only what changes for consumers of this crate.

## [Unreleased]

### Planned
Follow the roadmap in `rfcs/0006-nifti-ingestion.md`:
- `.nii.gz` support, behind a named feature with a `cargo deny`-approved
  pure-Rust inflate dependency — the same pattern `tpt-med-dicom`'s
  RLE/JPEG/JPEG-LS/JPEG 2000 features already established.
- Dual-file `.hdr`/`.img` NIfTI-1 support.
- `tpt-med-meshing` integration — deliberately left open by RFC 0006, since
  it means either generalising `SegmentationMask::threshold_hu`'s
  `DicomSeries`-concrete signature or an explicit adapter, and that decision
  deserves its own review.

## [0.1.0] - 2026-09-27

### Added
- **`NiftiVolume::{parse_file, parse_bytes}`** — parses an uncompressed,
  single-file NIfTI-1 (`.nii`) volume: dimensions, voxel spacing, RAS
  orientation, and per-voxel values already scaled by
  `scl_slope`/`scl_inter` (mirroring DICOM's `RescaleSlope`/`RescaleIntercept`
  convention, including "zero slope means no scaling").
- **Geometry resolution follows the NIfTI-1 spec's own precedence**: sform
  affine (`sform_code > 0`) preferred, else qform quaternion
  (`qform_code > 0`), else an axis-aligned spacing-only fallback
  ("Analyze-compatible" mode, both codes 0). The qform quaternion → rotation
  formula is `nifti1.h`'s own documented one.
- **8 datatypes**: `uint8`, `int8`, `int16`, `uint16`, `int32`, `uint32`,
  `float32`, `float64`. Others are rejected as `NiftiError::UnsupportedDatatype`,
  and a `bitpix`/`datatype` mismatch is rejected as `NiftiError::BitpixMismatch`
  — neither is silently coerced.
- **Endianness auto-detection** from `sizeof_hdr` (NIfTI has no separate
  byte-order flag; a reader tries both interpretations and takes whichever
  gives 348).
- **`SyntheticNiftiBuilder`** — writes real, parseable `.nii` bytes (sform or
  qform geometry, `u8` or `f32` data) so the test suite needs no real
  dataset, mirroring `tpt_med_dicom::SyntheticCtBuilder`'s role for DICOM.
- **`NiftiError`**: `Io`, `NotNifti`, `Gzipped` (a `.nii.gz` is named, not
  misparsed), `UnsupportedDatatype`, `BitpixMismatch`, `UnexpectedEof`,
  `BadValue`.

### Notes
- Values are `f64` regardless of on-disk type; the original type is kept on
  `NiftiVolume::datatype` for reference.
- Geometry is RAS (NIfTI's native convention); `tpt_med_geometry::ras_to_lps`
  converts into the LPS space `tpt_med_dicom::ImageFrame` uses.
- No shear is assumed in a sform affine — see the README's Known Limitations.
