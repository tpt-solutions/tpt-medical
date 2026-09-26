# Changelog

All notable changes to `tpt-med-dicom` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Planned
- Follow the ingestion roadmap in `rfcs/0001-dicom-ingestion.md`:
  - Encapsulated/compressed pixel data (JPEG, JPEG-LS, JPEG 2000, RLE) —
    currently rejected with `DicomError::CompressedPixelData`. Decoding must
    land behind a clearly named feature with the decompression dependencies
    `cargo deny` has approved, and must never silently fall back.
  - Multi-frame objects and DICOM networking (C-STORE, DICOMweb).
  - Quantitative-CT phantom calibration to replace the linear HU→density
    approximation in regulated pipelines.
- A NIfTI reader, so research-space volumes can enter the pipeline alongside
  DICOM without an external conversion step.

### Notes
- The HU→density relation is the linear CT approximation
  (`ρ = (HU + 1000) / 1000`). **Replacing it with a calibrated relation is
  semver-minor but changes every mesh modulus downstream**, so it requires an
  RFC and a V&V re-run.
- `DEFAULT_BONE_THRESHOLD_HU` is a screening default (200 HU, inside the
  130–300 HU literature band), not a clinical default.

## [0.1.0] - 2026-09-22

### Added
- **Transfer syntaxes** — implicit VR little endian (`1.2.840.10008.1.2`) and
  explicit VR little endian (`1.2.840.10008.1.2.1`), with
  `TransferSyntax::from_uid` and `TransferSyntax::implicit_vr(tag)`.
- **Part-10 parser** — `DicomParser::parse_file` and `parse_bytes`, handling
  the 128-byte preamble and the `DICM` magic. Sequences are parsed and skipped
  safely, including undefined-length sequences with depth tracking.
- `DicomElement` with typed accessors `as_us`, `as_is`, `as_ds_first`,
  `as_ds_vec` and `as_text`, each returning a typed `DicomError` on VR
  mismatch.
- `DicomParser::encode_element_explicit(tag, vr, value)` — a byte-accurate
  explicit-VR encoder, used by the synthetic generator and by tests.
- **Tag definitions** — patient name and ID, study/series/instance UIDs,
  `Modality`, `InstanceNumber`, `ImagePositionPatient`,
  `ImageOrientationPatient`, `SliceThickness`, `Rows`, `Columns`,
  `PixelSpacing`, `BitsAllocated`, `PixelRepresentation`,
  `RescaleIntercept`, `RescaleSlope` and `PixelData`.
- `Vr` with `from_ascii`, `ascii` and `uses_long_length`.
- **Series assembly** — `DicomSeries::load_from_dir` sorts slices by projection
  along the slice normal, the only robust ordering for oblique acquisitions;
  `from_slices` accepts a pre-ordered vector and validates consistency.
  `dims`, `hu_volume`, `hu_at` and `pixel_spacing` are the accessors.
- **Hounsfield Unit mapping** — `HounsfieldMapper::hu_to_density`,
  `density_to_youngs_modulus` and `hu_to_youngs_modulus`, with `BoneRegion`
  selecting the correlation. Power laws: cortical `E = 10500·ρ^2.0`,
  trabecular `E = 6850·ρ^1.49`; non-positive densities clamp to zero modulus.
  `DEFAULT_BONE_THRESHOLD_HU = 200.0`.
- **Synthetic CT generation** — `SyntheticCtBuilder` (dimensions, spacing,
  thickness, origin, patient id, arbitrary per-voxel HU function),
  `femur_phantom`, and `SyntheticCtSeries::{parse, write_to_dir}`, which
  writes real, parseable DICOM files so the workspace is testable without PHI.
- `DicomError` covering `Io`, `NotDicom`, `UnexpectedEof`,
  `CompressedPixelData`, `UnknownTransferSyntax`, `UnsupportedVr`, `BadValue`
  and `InconsistentSeries`.

### Known limitations
- Encapsulated/compressed pixel data is **rejected, never mis-parsed** — a
  wrong stress field that looks plausible is worse than a failed load.
  Decompression belongs at the archive boundary.
- Single-frame CT/MR slices only.

### Verification
- Both transfer syntaxes parse the committed synthetic series to identical HU
  volumes, which directly tests that the implicit and explicit VR paths agree.
- Rescale slope/intercept handling is pinned, including negative intercepts.
- Density and modulus correlations are locked to the exact power laws, with a
  dedicated air-voxel clamping test.
- The synthetic generator round-trips: `write_to_dir` then `load_from_dir`
  reproduces the HU volume exactly.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
