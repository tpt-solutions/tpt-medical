# Changelog

All notable changes to `tpt-med-nifti` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently**
on a six-week cadence. Workspace-level and cross-cutting changes are
recorded in the [workspace CHANGELOG](../../../CHANGELOG.md); this file
records only what changes for consumers of this crate.

## [Unreleased]

### Added
- **`gzip` cargo feature** (off by default) — decompresses `.nii.gz` and
  gzipped `.hdr`/`.img` parts via `flate2`'s pure-Rust backend
  (`rust_backend` = miniz_oxide; `cargo deny`-approved MIT OR Apache-2.0,
  no C in the build), following the named-feature pattern
  `tpt-med-dicom`'s RLE/JPEG/JPEG-LS/JPEG 2000 features established.
  Single-file decompression is **streaming and geometry-capped**: the
  header's own declared extent bounds what the reader will ever allocate,
  so a corrupt or hostile stream cannot expand unbounded (the same posture
  as the RLE decoder), and a stream shorter than its header claims fails as
  `UnexpectedEof` rather than being padded. Without the feature the typed
  `Gzipped` error remains, now pointing at the feature to enable.
- **New error variant `NiftiError::CorruptGzip(std::io::Error)`** — a gzip
  stream that could not be decompressed (bad gzip header, corrupt deflate
  stream, checksum mismatch), named rather than misparsed. Adding a variant
  is semver-minor under 0.x and exhaustive matches on `NiftiError` gain an
  arm.
- **Dual-file `.hdr`/`.img` NIfTI-1 support** —
  `NiftiVolume::{parse_dual_file, parse_dual_bytes}` accept the `ni1`
  magic, honour `vox_offset` as an offset into the `.img` (0 is the norm),
  and give an actionable error when handed the wrong layout's magic at
  either entry point (single-file entry point ↔ dual-file pointer, and
  vice versa). `SyntheticNiftiBuilder::build_dual` writes the pair for
  tests; with `gzip` on, either part may be gzipped (`.hdr.gz`/`.img.gz`,
  which the spec allows).
- **Meshing integration** — `tpt-med-meshing` gained
  `SegmentationMask::threshold_nifti(&NiftiVolume, min_hu)`: a narrow
  direct constructor rather than a shared trait, the shape RFC 0006's
  Unresolved Question said to pick while `NiftiVolume` stays the only
  non-DICOM source. NIfTI's RAS geometry is converted through
  `ras_to_lps` into the same LPS patient frame `threshold_hu` produces, so
  a mask is interchangeable regardless of which format it came from. This
  crate gained no new dependency of its own.

### Planned
- NIfTI-2 (the 2011 540-byte header with 64-bit dimensions) — still out
  of scope per RFC 0006; revisit if a workspace dataset needs it.

### Notes
- The dimension product and the voxel-data extent are now computed with
  checked arithmetic on every path (they were unchecked in 0.1.0's
  single-file path), so an adversarial header on a 32-bit target
  (`wasm32`) yields a typed `BadValue` instead of a wrapped extent.

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
