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
| **Dependencies** | [`tpt-med-core`](../../core/tpt-med-core), [`tpt-med-geometry`](../../core/tpt-med-geometry), [`tpt-med-units`](../../core/tpt-med-units) |
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
3. **Honesty about failure.** A compressed transfer syntax this crate cannot
   decode with confidence is *rejected with a typed error*, never silently
   mis-parsed. A wrong stress field that looks plausible is worse than a
   failed load — see the `jpeg2000` feature's signed-component handling below
   for what that looks like in practice: a real bug in the underlying codec
   crate is worked around where it safely can be, and the one case that
   genuinely cannot be resolved (a non-conformant file) is refused rather
   than guessed at.

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
- **`QctCalibration`: phantom-fitted HU → density**, replacing the fixed
  screening line with an ordinary-least-squares fit through measured
  calibration-phantom rod points. `HounsfieldMapper::hu_to_density` is now
  defined in terms of `QctCalibration::screening_default()`, so the fixed
  line and a real fit go through the exact same code path — see
  `hounsfield.rs` for the apparent-density-vs-BMD units caveat before
  feeding a clinical phantom's raw rod values into it.
- **`BmdToAshDensity`/`AshFraction`/`BmdToApparentDensity`: BMD → apparent
  density** (`bmd.rs`), for the common case where a phantom reports bone
  mineral density (mg/cm³, K₂HPO₄- or CaHA-equivalent) rather than apparent
  density directly. Ships no built-in preset relation for any manufacturer —
  both conversion stages refuse to construct without a citation string.
  `HounsfieldMapper::bmd_to_youngs_modulus` composes the result with the
  existing power laws.
- **`locate_phantom_centroid`/`sample_phantom_rods`/`PhantomModel`: phantom
  rod sampling** (`phantom.rs`), turning a CT scan of a calibration phantom
  into the `(HU, known_value)` pairs `QctCalibration::fit`/`BmdToAshDensity`
  consume. The centroid locator is manufacturer-agnostic; the rod layout
  (`PhantomModel`) is caller-supplied and cited, same discipline as above.
- **Synthetic CT generation** — `SyntheticCtBuilder` and `femur_phantom` write
  real, parseable DICOM files, so the whole test suite runs without PHI.
- **Element encoder** (`encode_element_explicit`) so tests and the synthetic
  generator produce byte-accurate output.
- **`rle` feature: RLE Lossless (`1.2.840.10008.1.2.5`) pixel data.** Off by
  default. With it enabled, encapsulated PackBits-compressed pixel data (PS3.5
  Annex G) is decoded directly, single-frame only. Needs no new dependency.
- **`jpeg` feature: classic JPEG family, via [`jpeg-decoder`](https://crates.io/crates/jpeg-decoder)
  (image-rs).** Off by default. Covers Baseline/Extended DCT
  (`1.2.840.10008.1.2.4.50`/`.51`, lossy — decoded pixel values are only an
  approximation of the originals, DCT quantization is lossy by design) and
  JPEG Lossless, Process 14 and Process 14 SV1 (`.57`/`.70`, exact — DPCM +
  Huffman, the same codec's other coding process). Single-component
  (grayscale) frames only.
- **`jpeg-ls` feature: JPEG-LS, via [`pure_jpegls`](https://crates.io/crates/pure_jpegls).**
  Off by default. Covers JPEG-LS Lossless (`.80`, exact) and Near-Lossless
  (`.81`, bounded per-sample error, not exact). Single-component only.
- **`jpeg2000` feature: JPEG 2000, via [`pdfluent-jpeg2000`](https://crates.io/crates/pdfluent-jpeg2000)
  (`hayro-jpeg2000`).** Off by default. Covers JPEG 2000 Lossless Only (`.90`)
  and JPEG 2000 (`.91`, lossless *or* lossy — the UID alone does not say
  which). Built with its `image`/`simd` extras disabled, so it pulls in no
  further dependencies. **Works around a real bug in the underlying crate:**
  it applies JPEG 2000's unsigned DC level shift to every component
  unconditionally, regardless of whether the codestream declares it signed.
  `decode_frame` re-reads that bit directly from the SIZ marker bytes (the
  crate discards it) and undoes the shift itself when the component really
  is signed. A file where the codestream's signed bit and the dataset's
  `PixelRepresentation` disagree is non-conformant and is rejected outright,
  since there is no safe way to resolve that disagreement — see `jpeg2000.rs`.

## Explicit Non-Features

The transfer syntaxes above are opt-in and behind their own cargo feature; the
default build still decodes only the two uncompressed syntaxes, and every
compressed syntax this crate does *not* implement (JPEG 2000 Part 2
multi-component, JPIP, and any retired Process not listed above) still
returns `DicomError::CompressedPixelData`, not garbage — decompress those at
the archive boundary. That is a deliberate architectural line, not an
oversight. Multi-frame objects, private tags with odd VRs, and DICOM
networking (C-STORE, DICOMweb) are likewise out of scope for v0.

## Conventions

- Pixels are stored as `i32` **raw stored values**; HUs are derived per slice.
- Geometry follows the DICOM patient coordinate system (**LPS**), exposed via
  `ImageFrame` from tags (0020,0032) and (0020,0037).
- `HU = stored × RescaleSlope + RescaleIntercept`.
- By default, apparent density uses the linear CT approximation:
  `ρ = (HU + 1000) / 1000` g/cm³, so water is 1.0 and air is 0.0
  (`HounsfieldMapper::hu_to_density`, equivalently
  `QctCalibration::screening_default()`). A `QctCalibration` fitted from a
  real calibration phantom's rods replaces this in regulated pipelines —
  see `hounsfield.rs` for the apparent-density-vs-BMD units caveat.
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
    SyntheticCtBuilder::femur_phantom(64, 64, 40)
        .patient_id("SYNTHETIC-0001")
        .build("1.2.826.0.1.3680043.9.7484.1.1")
        .write_to_dir(std::path::Path::new("test-data/dicom/synthetic_ct"))
}
```

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
| `HounsfieldMapper::hu_to_density` | HU → apparent density (g/cm³) via the fixed screening line |
| `HounsfieldMapper::density_to_youngs_modulus` | Density → modulus via the cortical or trabecular law |
| `HounsfieldMapper::hu_to_youngs_modulus` | Combined convenience (screening line) |
| `HounsfieldMapper::hu_to_youngs_modulus_calibrated` | Combined convenience, via a fitted `QctCalibration` |
| `HounsfieldMapper::DEFAULT_BONE_THRESHOLD_HU` | 200.0 HU |
| `QctCalibration::fit(&[(hu, value)])` | Least-squares line through calibration-phantom rod points |
| `QctCalibration::{screening_default, evaluate, hu_to_apparent_density, slope, intercept}` | The fixed screening line as a `QctCalibration`; evaluate the fit |
| `HounsfieldMapper::bmd_to_youngs_modulus` | Combined convenience, via a `BmdToApparentDensity` |
| `BmdToAshDensity::new(convention, slope, intercept, source)` | Cited BMD → ash-density line |
| `AshFraction::new(fraction, source)` | Cited ash → apparent-density fraction |
| `BmdToApparentDensity::{new, apparent_density, provenance}` | Composed two-stage conversion; `provenance()` returns both citations |
| `BmdConvention::{K2Hpo4Equivalent, HydroxyapatiteEquivalent}` | Which mineral-equivalent convention a BMD value uses |
| `locate_phantom_centroid(&slice, background_max_hu, min_area_px)` | Manufacturer-agnostic phantom centroid, by intensity thresholding |
| `sample_phantom_rods(&slices, &model, centroid, rotation_rad, roi_fraction)` | Per-rod `(mean_hu, known_value)` pairs from a `PhantomModel` |
| `PhantomModel::new(rods, source)`, `PhantomRod` | Cited rod layout (offset, angle, radius, known value) |
| `BoneRegion::{Cortical, Trabecular}` | Which correlation applies |
| `SyntheticCtBuilder` | Builder for synthetic series (dims, spacing, thickness, origin, patient id, arbitrary HU function) |
| `SyntheticCtSeries::{parse, write_to_dir}` | Round-trip: write real DICOM files, or parse them back |
| `SyntheticCtBuilder::femur_phantom` | Ready-made femoral shaft phantom |
| `DicomError` | `Io`, `NotDicom`, `UnexpectedEof`, `CompressedPixelData`, `UnknownTransferSyntax`, `UnsupportedVr`, `BadValue`, `InconsistentSeries`, `Calibration`, `Phantom` |
| `Result<T>` | Crate result alias |

## Verification

- Both transfer syntaxes parse the committed synthetic series
  (`test-data/dicom/synthetic_ct/`, 24 slices) to identical HU volumes — a
  direct test that the implicit and explicit VR paths agree.
- Rescale slope/intercept handling is pinned, including negative intercepts.
- Density and modulus correlations are locked to the exact power laws above,
  with a dedicated air-voxel clamping test.
- `QctCalibration::fit` is checked against a noiseless exact line (recovers
  the true slope/intercept), a noisy set of points with symmetric error
  (recovers the true line within tolerance, proving it minimises squared
  error rather than just interpolating), and is asserted to reject too few
  points, a degenerate (collinear-HU) rod set, and non-finite input. A
  dedicated test pins that `screening_default()` reproduces
  `hu_to_density`'s output exactly, since the latter is now defined in terms
  of the former.
- `BmdToAshDensity`/`AshFraction`/`BmdToApparentDensity` are checked against
  hand-computed round numbers for both stages, and both constructors are
  asserted to reject an empty/whitespace-only citation, non-finite
  coefficients, and an out-of-range ash fraction.
- `locate_phantom_centroid` is checked against a synthetic slice with a
  hand-placed circular high-HU region at a known, non-trivial (row, col) and
  radius, including an off-centre case, and is asserted to return `None`
  both when nothing reaches `min_area_px` and when the only above-threshold
  pixels are small disconnected noise. `sample_phantom_rods` is checked
  against a synthetic slice with rod regions at known `PhantomModel` offsets
  — the returned `(mean_hu, known_value)` pairs reproduce the exact known HU
  values, including averaging correctly across multiple slices — and is
  asserted to reject an out-of-bounds ROI, mismatched slice dimensions, an
  invalid `roi_fraction`, and an empty slice list.
- The synthetic generator round-trips: `write_to_dir` then `load_from_dir`
  reproduces the HU volume exactly, which is what makes the rest of the
  workspace testable without real patient data.
- `CompressedPixelData` is asserted to be returned, never swallowed, for the
  syntaxes still rejected.
- The `rle` feature has dedicated unit tests: literal and repeat PackBits
  runs, 8-bit and 16-bit (signed and unsigned), a run capped to the declared
  geometry, and truncated/malformed frames rejected rather than padded.
- The `jpeg` feature's Lossless (Process 14 SV1) path is checked against a
  hand-built, hand-verified bitstream (not a round-trip through the same
  encoder the decoder is tested against), covering a zero DPCM difference and
  both signs of a nonzero one.
- The `jpeg-ls` and `jpeg2000` features are checked with round-trip tests
  (encode via the same crate, decode via `decode_frame`) across 8-bit, 16-bit,
  and signed pixel representations where applicable, plus dimension-mismatch
  and malformed-input rejection.
- `series::encapsulated_pixel_data_tests::jpeg_ls_lossless_end_to_end` proves
  the whole path, not just the codec module in isolation: a real Part-10 byte
  stream with a `TransferSyntaxUID` of `1.2.840.10008.1.2.4.80` and PS3.5
  Annex A.4-shaped encapsulated fragments, parsed by `DicomParser::parse_bytes`
  and decoded into `DicomSlice::pixel_data`.

## Known Limitations

- **Two transfer syntaxes.** Implicit and explicit VR little endian only. Big
  endian (`1.2.840.10008.1.2.2`) is a real syntax still emitted by some
  archive exports and is not implemented.
- **Compressed pixel data decoders are all opt-in cargo features**, off by
  default: `rle`, `jpeg` (Baseline/Extended lossy, Lossless Process 14/SV1
  exact), `jpeg-ls` (Lossless exact, Near-Lossless bounded-error), `jpeg2000`
  (Lossless Only exact, `.91` either). Without the matching feature, an object
  using that transfer syntax still yields `DicomError::CompressedPixelData`.
  JPEG 2000 Part 2 multi-component and JPIP-referenced pixel data have no
  decoder at all yet.
- **The `jpeg2000` feature only trusts a *conformant* signed
  `PixelRepresentation`.** `pdfluent-jpeg2000` applies the unsigned DC
  level-shift to every component regardless of whether the codestream
  declared it signed (its own source says so); `jpeg2000::decode_frame`
  compensates for this correctly by re-reading the codestream's own SIZ
  signed bit and undoing the shift when needed. What it cannot do is resolve
  a file where that bit and the dataset's `PixelRepresentation` disagree —
  such a file is non-conformant, and `decode_frame` rejects it rather than
  guess which one to believe.
- **Lossy transfer syntaxes decode approximate pixel values, not exact stored
  values**, by construction (`jpeg`'s Baseline/Extended DCT, `jpeg-ls`'s
  Near-Lossless, and whichever encoder wrote a `.91` JPEG 2000 stream lossily).
  A HU value derived from one of these is not the exact number the scanner
  produced. Treat it the same way you would treat any other lossy source.
- **Single-frame only.** Enhanced multi-frame CT (a common Siemens/GE
  representation) and MR object hierarchies are not supported.
- **Incomplete tag coverage.** Only the tags needed for geometry and HU
  mapping are decoded. Window/level, pixel padding, slice position sorting by
  `INSTANCE_NUMBER`, private tags and structured reports are not.
- **`HounsfieldMapper::hu_to_density` is still a screening estimate by
  default**, not quantitative CT — it ignores the scanner's actual
  calibration, so absolute density (and therefore absolute modulus) from it
  alone is not a measurement. `QctCalibration::fit` lets a caller supply a
  real phantom's rod measurements instead; `sample_phantom_rods` +
  `locate_phantom_centroid` (`rfcs/0008-phantom-rod-sampling.md`) turn a
  phantom scan into those points, and `BmdToApparentDensity`
  (`rfcs/0007-bmd-apparent-density-conversion.md`) converts a BMD-convention
  rod value to apparent density. Neither mechanism ships a built-in
  manufacturer preset or literature coefficient — a caller still supplies
  (and cites) the phantom's rod layout and the published conversion
  coefficients; only the pixel-hunting and manual averaging are automated.
  Rotation and slice selection also stay caller-supplied — see the RFC for
  why blind detection of either was rejected.
- **Slice ordering assumes a single acquisition.** Sorting by projection along
  the slice normal handles oblique acquisitions well, but a series containing
  multiple stacks (localiser plus scan, or overlapping repeats) is not
  de-duplicated or clustered.
- **No DICOM networking and no DICOMweb**, so there is no C-STORE receiver and
  no query/retrieve.
- **The synthetic generator is a phantom, not an anatomy.** It produces
  geometric primitives with realistic HU values and DICOM structure, which is
  what the test suite needs, and is not a substitute for a real dataset when
  judging whether a thresholding choice is clinically sensible.

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
