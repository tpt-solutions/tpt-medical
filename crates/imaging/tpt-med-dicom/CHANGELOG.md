# Changelog

All notable changes to `tpt-med-dicom` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- **Multi-frame (Enhanced CT/MR) ingestion for the uncompressed syntaxes**
  — RFC 0001's v1 roadmap item 2. `DicomParser::parse_bytes_all` /
  `parse_file_all` return one `DicomSlice` **per frame**, mapping the
  object's per-frame functional groups (PS3.3 C.7.6.6/C.7.6.16) onto the
  slice list: per-frame Plane Position (required — it is what
  distinguishes the frames), with Plane Orientation and Pixel Value
  Transformation rescale honoured per frame when present and the Shared
  Functional Groups (orientation, pixel measures, rescale) re-absorbed
  through the very same match arms as top-level tags, so shared and
  top-level encodings cannot drift. `DicomSeries::load_from_dir`
  assembles multi-frame and single-frame files into one series
  positionally. Sequence handling is now structural: items are resolved
  into element lists (defined- and undefined-length alike, so the
  dcm4che-style and GDCM-style encodings both parse) rather than skipped.
  Guardrails: `parse_bytes`/`parse_file` refuse a multi-frame object with
  an error naming the `_all` entry points rather than silently
  truncating to frame 0; a declared `NumberOfFrames` that disagrees with
  the per-frame item count, and a frame without a Plane Position, are
  named `BadValue`s; multi-frame **encapsulated** pixel data remains a
  named rejection (one decode per payload cannot represent N frames).
  `DicomElement` gained `items` — the resolved sequence items. Eight new
  tests: functional-group mapping (positions, shared geometry, per-frame
  rescale, per-frame pixel payloads), undefined-length equivalence,
  implicit-VR multi-frame, the `_all` refusal, the count mismatch, the
  missing plane position, series assembly, and a single-frame
  `NumberOfFrames = 1` file through the unchanged entry point.
- **Multi-frame over the compressed syntaxes too** — completing the
  roadmap item for every transfer syntax this crate decodes. The Basic
  Offset Table (PS3.5 Annex A.4) is resolved instead of rejected: one u32
  byte offset per frame into the concatenated fragment stream, validated
  (strictly increasing, within the stream, count == `NumberOfFrames`);
  with an empty table, one fragment per frame is accepted when the counts
  agree. Each frame decodes through the same feature-gated codec path,
  with per-frame functional-group geometry exactly as in the native case.
  An empty table whose fragment count disagrees with the frame count is a
  named "boundaries cannot be recovered" error rather than a guess.
  `DicomElement::pixel_fragments` (`PixelFragments { fragments,
  basic_offset_table }`) carries the resolved structure. Five new tests
  (RLE frames through the BOT, the empty-table fragment-per-frame form,
  and the three named rejections).
- **Fixed: Part-10 implicit-VR files failed to parse at all.** The file
  meta-group loop read the first dataset element with explicit-VR
  parsing, so an implicit dataset's 32-bit length header was
  misinterpreted as a VR and rejected (`UnsupportedVr` on
  `(0008,0005)`-shaped headers). The meta/dataset boundary is now found
  by peeking the next tag's group (the existing per-element implicit
  reader was already correct — the bug was only in the hand-over).
  Found by the new implicit-VR multi-frame test, which real archive
  exports would have hit.
- `detect_phantom_rotation` / `RotationDetection`: **phantom rotation
  detection** — the mechanical half of RFC 0008's deferred follow-up. A
  coarse sweep (2π/`coarse_steps`, ≥ 8) over in-plane rotations scored by
  the sampled-vs-known least squares, refined by grid-and-shrink rounds
  (the score's valley is flat-bottomed under pixel quantisation, which
  defeats ternary search). `rms_residual_hu` is the caller's reject
  signal for a layout that matches nothing. The *named vendor-phantom
  library* half of that item stays open deliberately: manufacturer
  rod layouts must come from datasheets, not baked-in guesses. Three new
  tests (known-rotation recovery at the fixture's resolution limit,
  honest residual on a non-matching image, coarse-step validation).
- **`jpeg2000` feature now covers JPEG 2000 Part 2 Multi-component**
  (`1.2.840.10008.1.2.4.92` Lossless Only, `.93` lossless-or-lossy) — two new
  `TransferSyntax` variants routed through the exact same `jpeg2000::decode_frame`
  path as `.90`/`.91`, since Part 2 extends Part 1's codestream syntax rather
  than replacing it. "Multi-component" names the transfer syntax, not a new
  decode capability: `decode_frame`'s existing single-component restriction
  is unchanged, so a `.92`/`.93` file with `SamplesPerPixel = 1` (the common
  CT/MR case) decodes normally, and a genuinely multi-sample-per-pixel
  (color) one is rejected exactly as a color `.90`/`.91` file already would
  be — this crate has no `SamplesPerPixel`/`PlanarConfiguration` handling
  anywhere. A codestream using a real Part 2 extended multi-component
  transform is rejected by the underlying crate's own unrecognised-marker
  check (verified by reading `pdfluent-jpeg2000`'s marker-parsing loop),
  never silently mis-decoded.
- **Series-level integration tests for every JPEG 2000 transfer syntax**
  (`series::jpeg2000_encapsulated_pixel_data_tests`), closing the parity gap
  where only `jpeg-ls` had a Part-10-byte-stream-to-`DicomSlice` test. A real
  Part-10 stream with PS3.5 Annex A.4-shaped encapsulated fragments is built
  for `.90`/`.91`/`.92`/`.93` and parsed by `DicomParser::parse_bytes`, proving
  the transfer-syntax routing and the parser's own fragment collection rather
  than a pre-assembled frame. The signed-bit level-shift correction and its
  non-conformant-mismatch rejection are pinned on that path too, not only in
  `decode_frame`'s unit tests. The 2x2 J2C fixture is now a single
  `pub(crate)` constant in `jpeg2000.rs` shared by both test sites.
- **`rle` cargo feature: RLE Lossless (`1.2.840.10008.1.2.5`) pixel data.**
  Off by default — without it an RLE object still yields
  `DicomError::CompressedPixelData`, unchanged from before. With it,
  `src/rle.rs` decodes the PackBits segments per PS3.5 Annex G, single-frame
  only (a non-empty Basic Offset Table is rejected as unsupported
  multi-frame), and every run is capped to `rows * columns` so a corrupt or
  hostile fragment cannot expand without bound. Needs no new dependency, so
  there is nothing new for `cargo deny` to approve.
- The parser now recognises the general encapsulated pixel-data item stream
  (PS3.5 Annex A.4) for any transfer syntax whose dataset is explicit VR
  little endian — the RLE decoder was the first consumer, and the `jpeg`/
  `jpeg-ls`/`jpeg2000` decoders below reuse the same fragment collection
  without any parser changes.
- **`jpeg` cargo feature: the classic ITU-T T.81 JPEG family**, via
  [`jpeg-decoder`](https://crates.io/crates/jpeg-decoder) (image-rs, MIT OR
  Apache-2.0). Off by default. One dependency covers two DICOM transfer
  syntaxes because the crate implements both T.81 coding processes:
  - Baseline / Extended DCT (`1.2.840.10008.1.2.4.50` / `.51`) — lossy;
    decoded values are an approximation of the originals.
  - Lossless, Process 14 and Process 14 SV1 (`.57` / `.70`, the latter being
    the "default lossless JPEG" transfer syntax) — exact, DPCM + Huffman.

  `src/jpeg.rs`. Single-component (grayscale) frames only.
- **`jpeg-ls` cargo feature: JPEG-LS**, via
  [`pure_jpegls`](https://crates.io/crates/pure_jpegls) (MIT OR Apache-2.0).
  Off by default. Covers JPEG-LS Lossless (`.80`, exact) and Near-Lossless
  (`.81`, bounded per-sample error — not exact). `src/jpeg_ls.rs`.
  Single-component only.
- **`jpeg2000` cargo feature: JPEG 2000**, via
  [`pdfluent-jpeg2000`](https://crates.io/crates/pdfluent-jpeg2000)
  (`hayro-jpeg2000`, Apache-2.0 OR MIT), built with its `image`/`simd` extras
  disabled so it adds no further dependency. Off by default. Covers JPEG
  2000 Lossless Only (`.90`) and JPEG 2000 (`.91`, lossless or lossy — the
  UID alone does not say which). `src/jpeg2000.rs`. **Works around a real bug
  in the underlying crate**: it applies JPEG 2000's unsigned DC level shift
  to every component unconditionally, regardless of the codestream's own
  signed bit (which the crate reads and discards, by its own admission not
  knowing how to handle the signed case). Since the shift is a fixed, known
  offset, `decode_frame` re-reads that bit directly from the raw SIZ marker
  bytes and undoes the shift itself when the component really is signed. A
  file where the codestream's signed bit and the dataset's
  `PixelRepresentation` disagree is non-conformant and has no safe
  resolution, so it is rejected rather than guessed at. See the module doc
  comment.
- `TransferSyntax` gained one variant per newly-recognised UID
  (`JpegBaseline`, `JpegExtended`, `JpegLossless`, `JpegLosslessSv1`,
  `JpegLsLossless`, `JpegLsNearLossless`, `Jpeg2000Lossless`, `Jpeg2000`).
  `from_uid` now recognises all of them unconditionally — the dataset itself
  parses whether or not the matching decode feature is enabled, since only
  the `PixelData` element's own decode is feature-gated.
- **`QctCalibration`: a fitted HU → density calibration.** `QctCalibration::fit`
  takes a calibration phantom's measured `(HU, known_value)` rod points and
  fits a line by ordinary least squares, replacing the fixed two-point
  screening line (water 0 HU → 1.0 g/cm³, air −1000 HU → 0.0) with one
  actually measured against that scan. `HounsfieldMapper::hu_to_density` is
  now defined as `QctCalibration::screening_default().hu_to_apparent_density`,
  so the fixed line is expressible in the same type rather than a separate
  hardcoded formula — this is a refactor, not a behaviour change, and is
  pinned by a dedicated equivalence test.
  `HounsfieldMapper::hu_to_youngs_modulus_calibrated` is the calibrated
  counterpart to the existing `hu_to_youngs_modulus`. Does **not** solve
  automatic phantom-rod detection (points must be supplied, not extracted
  from an image) or BMD→apparent-density conversion (most clinical QCT
  phantoms report bone mineral density, mg/cm³ K₂HPO₄/CaHA-equivalent, a
  different physical quantity from the apparent density the modulus power
  laws expect — see the module docs). New `DicomError::Calibration` variant
  for a fit that can't be made (too few points, degenerate HU spread,
  non-finite input).
- **`BmdToAshDensity` / `AshFraction` / `BmdToApparentDensity`** (`src/bmd.rs`,
  `rfcs/0007-bmd-apparent-density-conversion.md`): the documented BMD→
  apparent-density conversion `QctCalibration`'s own docs named as missing.
  Ships **no** built-in preset relation for any manufacturer or paper — both
  stages (BMD→ash density, ash→apparent density) refuse to construct without
  a non-empty citation string, making provenance structurally mandatory
  rather than an optional doc comment. `HounsfieldMapper::bmd_to_youngs_modulus`
  is the convenience composing a `BmdToApparentDensity` with the existing
  power laws. New `BmdConvention` enum (`K2Hpo4Equivalent`,
  `HydroxyapatiteEquivalent`) is metadata only — it does not convert between
  the two conventions.
- **`locate_phantom_centroid` / `sample_phantom_rods` / `PhantomModel`**
  (`src/phantom.rs`, `rfcs/0008-phantom-rod-sampling.md`): turns a CT scan of
  a calibration phantom into the `(HU, known_value)` pairs `QctCalibration::fit`
  and `BmdToAshDensity` consume. `locate_phantom_centroid` is manufacturer-
  agnostic (thresholded connected-component centroid); `PhantomModel` (rod
  layout + known values) is caller-supplied and cited, same discipline as
  `BmdToAshDensity`. Rotation and slice selection stay caller-supplied in v0
  — see the RFC for why blind detection of either was rejected. New
  `DicomError::Phantom` variant.

### Planned
- Follow the ingestion roadmap in `rfcs/0001-dicom-ingestion.md`:
  - JPEG 2000 Part 2 multi-component *codestream* decoding — the `.92`/`.93`
    transfer syntaxes decode, but "multi-component" is not a decode
    capability (see Added above): a real Part 2 extended multi-component
    transform is still refused, and multi-sample-per-pixel data is still out
    of scope for this crate. (RLE/JPEG/JPEG-LS/JPEG 2000 Part 1 are done —
    see Added above, including correct signed-component handling.)
    **Spiked 2026-09-27** (`rfcs/0010-jpeg2000-part2-spike.md`): the
    rejection is verified at source, no pure-Rust route exists today, and the
    gating input is a purchase of ISO/IEC 15444-2 rather than engineering
    time. Recommend deferring.
  - JPIP-referenced pixel data — no decoder, and none should be written here:
    it is a network reference to pixel data held elsewhere, so it needs a
    transport story (resolved at the archive boundary) before decoding is even
    in scope. Still `DicomError::CompressedPixelData`.
  - DICOM networking (C-STORE, DICOMweb). Multi-frame ingestion is done
    for the uncompressed syntaxes and, behind each codec feature, for the
    compressed ones (see Added above).
- A built-in library of named, cited `PhantomModel`s for common commercial
  phantoms — the remaining half after the rotation detection delivered
  above (`detect_phantom_rotation`); the layouts and known values must come
  from manufacturer datasheets with citations, not baked-in approximations.

### Notes
- The default HU→density relation is still the linear CT approximation
  (`ρ = (HU + 1000) / 1000`, i.e. `QctCalibration::screening_default()`).
  `QctCalibration` is additive — nothing is forced to use it. **Changing what
  `hu_to_density` itself returns by default would be semver-minor but change
  every mesh modulus downstream**, so that (as opposed to offering the opt-in
  calibrated path added here) still requires an RFC and a V&V re-run.
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
