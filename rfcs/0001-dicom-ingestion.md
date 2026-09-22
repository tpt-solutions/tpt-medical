# RFC 0001: DICOM Ingestion

- **Status:** Accepted
- **Started:** 2026-09-20
- **Crate:** `tpt-med-dicom` (implementation), `tpt-med-meshing` (downstream)

## Summary

Define the scope, transfer-syntax policy, and evolution path for pure-Rust
DICOM ingestion, from the v0 uncompressed parser to a validated, QCT-capable
ingestion layer.

## Motivation

Patient-specific modelling starts from DICOM CT. The ecosystem options are:

- **dicom-rs** (MIT): complete, actively maintained — but brings a large
  dependency tree for pipelines that only need a handful of tags and pixel
  data, and its API churns across 0.x releases.
- **DCMTK / GDCM bindings**: C++ FFI breaks the WASM/zero-cloud story and
  complicates the permissive-only license audit.

A dependency-light parser keeps the WASM footprint small, the license chain
trivially clean, and the audit surface (SECURITY-relevant parsing code)
in-repo and readable — which matters when the parser is a submission asset.

## Design

### v0 (implemented)

- Transfer syntaxes: implicit VR LE (`1.2.840.10008.1.2`) and explicit VR LE
  (`1.2.840.10008.1.2.1`). Encapsulated/compressed syntaxes are **rejected**
  with `DicomError::CompressedPixelData` — never silently mis-parsed.
- Single-frame CT slices; sequences parsed and skipped safely (undefined
  length handled with depth tracking).
- Elements required for geometry and HU mapping: ImagePosition/Orientation,
  PixelSpacing, SliceThickness, RescaleSlope/Intercept, Rows/Columns/
  BitsAllocated/PixelRepresentation, PixelData, Modality, SeriesInstanceUID,
  InstanceNumber, StudyDate, PatientID (carried but **never logged** — see
  `tpt-med-core`'s privacy model).
- Series assembly: sort along the slice normal (row × column cosines),
  verify consistent rows/columns.
- HU mapping: `stored × slope + intercept`; density via the linear CT
  approximation; moduli via Morgan–Keaveny-style power laws split at
  ρ = 1.3 g/cm³ (cortical/trabecular).

### v1 roadmap

1. **Decompression gateway**: accept JPEG/JPEG-LS/JPEG 2000 by delegating to
   an opt-in codec crate behind a feature flag, keeping the default build
   dependency-free. The gateway is the *only* place pixel decoding happens.
2. **Enhanced (multiframe) CT** — PS3.3 C.7.6.6: per-frame functional groups
   map onto the existing `DicomSlice` list.
3. **QCT calibration**: phantom-based HU→mgHA/cm³ calibration replacing the
   linear screening approximation when a calibration object is present in
   the series.
4. **Regulatory evidence pack**: the parser's conformance statement
   (handled tags, VR table) generated from the tag table so submitted
   documentation cannot drift from the code.

### Explicitly out of scope

Structured reporting, presentation states, MPPS/network roles (C-GET/
C-MOVE), and modality worklist — ingestion only.

## Alternatives considered

- **Adopt dicom-rs wholesale**: rejected for v0 on dependency-weight and
  API-stability grounds; revisit for the decompression gateway (its codec
  crates are candidates for the feature-gated integration).
- **Vendor-neutral pixel pipeline (NIfTI-first)**: keeps DICOM at the edge,
  but every hospital export starts as DICOM; deferring DICOM pushes
  de-identification complexity onto users.

## Unresolved questions

- Should de-identification (tag scrubbing on export) live in this crate or
  in `tpt-med-fda`'s export path? Current lean: ingestion scrubbing with an
  auditable "scrubbed by" event.
