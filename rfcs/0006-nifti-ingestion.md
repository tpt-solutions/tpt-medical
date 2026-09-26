# RFC 0006: NIfTI Ingestion

- **Status:** Accepted
- **Started:** 2026-09-27
- **Crates:** `tpt-med-nifti` (implementation)

## Summary

Add `tpt-med-nifti`, a pure-Rust, zero-dependency parser for uncompressed
single-file NIfTI-1 (`.nii`) volumes, so a research-space CT/MR volume (a
`dcm2niix` export, a public dataset, a 3D Slicer save) can enter the imaging
pipeline without an external conversion step back to DICOM. It does **not**
wire NIfTI volumes into `tpt-med-meshing`'s existing mesher — that mesher is
concretely typed to `tpt_med_dicom::DicomSeries` today, and generalising it is
a separate, explicitly out-of-scope decision (see Unresolved Questions).

## Motivation

`tpt-med-dicom` covers DICOM (RFC 0001). DICOM is what hospital archives
export, but it is not what research tooling uses once a scan is out of the
archive: `dcm2niix`, FSL, FreeSurfer, ANTs, and most public imaging datasets
(e.g. the Visible Human Project derivatives, many QIN/TCIA reprocessed sets)
work in NIfTI. Today, using one of those volumes in this pipeline means
converting it back to DICOM first with an external tool — an awkward step for
something the pipeline should be able to read directly, and mentioned as a
gap in RFC 0001's "Alternatives considered" (rejecting a NIfTI-first design
does not mean rejecting NIfTI ingestion entirely).

**Question of interest:** given a NIfTI-1 volume, what are its voxel
dimensions, physical spacing, and orientation relative to patient anatomy,
and what is the physical value ("intensity", or HU-equivalent for a CT
export) at each voxel? **Model risk:** low for a screening/research tool —
this crate parses a file format, computes no biomechanics itself. **Model
influence:** whatever consumes the parsed volume (e.g. a future NIfTI path
into meshing) inherits the same influence its DICOM equivalent has today, so
correctness here is exactly as consequential as `tpt-med-dicom`'s parser
correctness, even though this RFC's own scope is "just" ingestion.

## Detailed design

### Format scope (v0)

- **NIfTI-1 only.** NIfTI-2 (a 2011 revision with a 540-byte header, for
  volumes needing 64-bit dimension fields) is not implemented — every
  practical CT/MR volume from the tools named above fits in NIfTI-1's 32-bit
  dimension limits by an enormous margin.
- **Single-file `.nii` only** (magic `n+1`, header and voxel data in one
  file). The dual-file `.hdr`/`.img` form (magic `ni1`) needs two files
  coordinated by the caller and is out of scope for v0 — see Unresolved
  Questions.
- **Uncompressed only.** A gzip-compressed `.nii.gz` (by far the most common
  form in practice, since NIfTI files are large and gzip-transparent readers
  are ubiquitous) is detected by its magic bytes (`1F 8B`) and rejected with
  a typed, actionable error — exactly the `tpt-med-dicom` precedent for a
  transfer syntax it does not (yet) decode. Decompression needs a vetted
  dependency decision (a pure-Rust inflate crate) the same way RLE/JPEG/
  JPEG-LS/JPEG 2000 needed one in `tpt-med-dicom`; that is future work behind
  a named feature (`gzip`), not silently bundled into v0's zero-dependency
  build.

### Header parsing

The NIfTI-1 header is a fixed 348-byte struct (unchanged since 2004; this is
`nifti1.h`, not something that has multiple competing revisions to track).
Parsed fields, and only these — matching `tpt-med-dicom`'s "only what
geometry and value mapping need" policy:

| Field | Byte offset | Type | Use |
|---|---|---|---|
| `sizeof_hdr` | 0 | `i32` | Must be 348; also used for endianness detection (see below) |
| `dim[0..8]` | 40 | `i16[8]` | `dim[0]` = number of dimensions; `dim[1..4]` = nx, ny, nz |
| `datatype` | 70 | `i16` | Voxel storage type (see below) |
| `bitpix` | 72 | `i16` | Cross-checked against `datatype`, not trusted alone |
| `pixdim[0..8]` | 76 | `f32[8]` | `pixdim[0]` = `qfac` (±1, qform handedness); `pixdim[1..4]` = voxel spacing (mm) |
| `vox_offset` | 108 | `f32` | Byte offset to voxel data (single-file form) |
| `scl_slope` | 112 | `f32` | Value scaling: `value = raw * scl_slope + scl_inter` if `scl_slope != 0`, else `raw` unchanged (mirrors `RescaleSlope`/`RescaleIntercept` with the same "0 slope means no scaling" convention) |
| `scl_inter` | 116 | `f32` | — |
| `qform_code` | 252 | `i16` | 0 = qform absent |
| `sform_code` | 254 | `i16` | 0 = sform absent |
| `quatern_b/c/d` | 256, 260, 264 | `f32` | Quaternion components (`quatern_a` is derived, not stored) |
| `qoffset_x/y/z` | 268, 272, 276 | `f32` | qform translation (mm, RAS) |
| `srow_x/y/z[0..4]` | 280, 296, 312 | `f32[4]` each | sform affine rows (direct voxel→RAS mm affine; the implicit 4th row is `[0,0,0,1]`) |
| `magic` | 344 | `4` bytes | `"n+1"` (single-file) required; first 3 bytes checked, not the trailing byte, since some real-world writers pad it inconsistently |

Every other field (`descrip`, `aux_file`, `cal_min/max`, `intent_*`, slice
timing, …) is skipped, not decoded — this crate answers "what are the values
and where are they in space," nothing about acquisition metadata.

**Endianness.** NIfTI has no separate byte-order flag; a reader detects it by
trying `sizeof_hdr` both ways. If little-endian interpretation gives 348, the
file is little-endian; if big-endian interpretation gives 348, it is
big-endian; if neither, the file is rejected as not-NIfTI. Every subsequent
multi-byte field is read with whichever endianness matched.

**Datatype.** Supported: `DT_UINT8` (2), `DT_INT16` (4), `DT_INT32` (8),
`DT_FLOAT32` (16), `DT_FLOAT64` (64), `DT_INT8` (256), `DT_UINT16` (512),
`DT_UINT32` (768) — covering every datatype a CT or MR NIfTI export
realistically uses. Others (RGB, complex, 64-bit integer, binary) are
rejected with a typed error naming the unsupported code, not guessed at.
`bitpix` is cross-checked against the expected size for `datatype` and
rejected on mismatch (a corrupt or adversarial header claiming `datatype` and
`bitpix` that disagree must not be trusted for the read-size computation).

### Geometry

Orientation resolution follows the NIfTI-1 spec's own precedence exactly:

1. `sform_code > 0`: the affine is `srow_x/y/z` directly — three rows of
   `[a, b, c, d]` meaning `ras = [a,b,c]·(i,j,k) + d`. This is what
   `dcm2niix` and most modern writers populate.
2. Else, `qform_code > 0`: build the rotation from the quaternion
   (`quatern_a = sqrt(max(0, 1 − b² − c² − d²))`, standard quaternion → 3×3
   rotation, per `nifti1.h`'s documented formula), scale its columns by
   `pixdim[1]`, `pixdim[2]`, `qfac·pixdim[3]` (`qfac` from `pixdim[0]`,
   treated as `1.0` if it is exactly `0` per spec), and translate by
   `qoffset_x/y/z`.
3. Else (`Analyze-compatible` mode, both codes 0): axis-aligned, spacing-only
   affine — `diag(pixdim[1], pixdim[2], pixdim[3])`, origin at the volume
   corner. No rotation information exists in this mode; a caller relying on
   patient-relative orientation from an Analyze-compatible file is trusting
   the file's dimension order alone, same as any Analyze reader.

The resulting affine's rotation part is exposed as `tpt_med_geometry::Mat3`
and its translation as `tpt_med_geometry::Vec3`, both already in **RAS**
(NIfTI's native convention) — `tpt_med_geometry::ras_to_lps` converts into
the same **LPS** space `tpt_med_dicom::ImageFrame` uses, for a caller mixing
NIfTI and DICOM inputs in one patient frame.

### Public API (sketch)

```rust
pub struct NiftiVolume {
    pub dims: (usize, usize, usize),        // nx, ny, nz
    pub voxel_spacing: (f64, f64, f64),     // mm
    pub origin: Vec3,                        // RAS mm, first voxel (0,0,0)
    pub rotation: Mat3,                      // unit columns: voxel axes in RAS
    pub values: Vec<f64>,                    // scaled; index = i + nx*(j + ny*k)
}

impl NiftiVolume {
    pub fn parse_file(path: &Path) -> Result<Self>;
    pub fn parse_bytes(bytes: &[u8]) -> Result<Self>;
    pub fn value_at(&self, i: usize, j: usize, k: usize) -> Option<f64>;
    pub fn voxel_position(&self, i: usize, j: usize, k: usize) -> Vec3; // RAS mm
}
```

No `HounsfieldMapper`-equivalent lives in this crate: `scl_slope`/`scl_inter`
give a scaled physical value, but NIfTI carries no `Modality`-equivalent tag
saying that value *is* HU. A caller who knows their volume is a CT export
(and that the exporting tool baked `RescaleSlope`/`RescaleIntercept` into
`scl_slope`/`scl_inter`, which `dcm2niix` does) can feed `values` straight
into `tpt_med_dicom::HounsfieldMapper` or `QctCalibration` — both take a bare
`f64`, with no DICOM-specific coupling — without this crate depending on
`tpt-med-dicom` at all.

### Error handling

A dedicated `NiftiError` (mirroring `DicomError`'s shape): `Io`, `NotNifti`
(bad magic/endianness), `Gzipped` (detected `.nii.gz`, named as unsupported
rather than misparsed), `UnsupportedDatatype(i16)`, `BitpixMismatch{ datatype,
bitpix }`, `UnexpectedEof`, `BadValue { reason }` (header field out of range —
zero/negative dimension, `vox_offset` before the header, etc.).
`#![forbid(unsafe_code)]`, matching every other crate in the workspace.

### Alternatives considered

- **Depend on an existing `nifti`-family crate.** A mature `nifti` crate
  exists on crates.io, but RFC 0001 already made this exact call for DICOM
  (dependency weight, WASM footprint, API-churn risk, in-repo auditability
  for a crate that is a submission asset) and nothing about NIfTI changes
  that reasoning — the format is, if anything, simpler than DICOM to
  hand-roll (one fixed-size header, no VR/tag system).
- **Wire directly into `tpt-med-meshing` now.** Rejected for this RFC: it
  would mean either changing `SegmentationMask::threshold_hu`'s signature
  (a `DicomSeries`-typed API today) or duplicating the mesher for NIfTI
  input, and either choice deserves its own review rather than riding in on
  an ingestion RFC. See Unresolved Questions.
- **NIfTI-2 support now.** Deferred: no tool in the motivation's list emits
  it by default, and its only reason to exist (dimensions needing more than
  32 bits) does not apply to any CT/MR volume size in this pipeline's
  reach.

### Drawbacks

- Two research-format-adjacent crates (`tpt-med-dicom`, `tpt-med-nifti`) with
  independent, not-yet-unified geometry/value plumbing, rather than one
  shared abstraction — accepted deliberately rather than forcing a premature
  shared trait before a second real consumer (the meshing integration)
  exists to inform its shape.
- `.nii.gz`, by far the most common NIfTI file extension in the wild, is
  explicitly not readable by v0. Anyone handed a typical public dataset
  download will hit `NiftiError::Gzipped` immediately and need to
  decompress externally until the `gzip` feature lands.
- No qform *validation* beyond the quaternion formula itself — a file with
  `qform_code > 0` but a non-unit `(b,c,d)` (i.e. `b²+c²+d² > 1`, which the
  formula clamps rather than rejects) is accepted with a degraded rotation
  rather than refused. Revisit if this proves to matter in practice; it is
  a narrower version of the same "trust the header, verify what's cheap to
  verify" posture `tpt-med-dicom` already takes with DICOM's own headers.

## Verification strategy

- **Code verification:** the quaternion → rotation formula is checked against
  a hand-derived closed-form case (`quatern_d = 1`, `b = c = 0` ⇒ a 180°
  rotation about the k-axis, `diag(-1, -1, 1)` — worked by hand in this RFC's
  review, not copied from an external tool), plus a property test that the
  rotation matrix is orthonormal (`RᵀR = I`) for a spread of unit quaternions,
  which catches an algebra transcription error independent of any single
  worked example.
- **Calculation verification:** not applicable — this crate performs no
  discretisation or iterative solve; it is a direct decode.
- **Format-decode verification**, the same posture as the DICOM codecs
  (RLE/JPEG/JPEG-LS/JPEG 2000 all built a hand-verified minimal fixture
  rather than only round-tripping through their own encoder): a synthetic
  writer (`NiftiVolume::write_to_file`, mirroring `SyntheticCtBuilder`) is
  built alongside the parser so `write_to_file` → `parse_file` round-trips
  exercise the real byte layout, and one hand-built minimal header (magic,
  dims, sform, one voxel) is parsed and checked field-by-field against
  values computed by hand from the bytes, not against the crate's own
  writer.
- **Validation:** no external reference dataset — there is no
  `test-data/nifti/` real-world file (by design; PHI-free synthetic fixtures
  only, matching `tpt-med-dicom`'s policy) and no golden dataset entry, since
  this crate produces no simulation output. `test-data/golden/` is
  unaffected.
- **What remains unverified:** behaviour against real-world files from actual
  writers (`dcm2niix`, FSL, SPM, ITK-SNAP) — only hand-built and
  self-round-tripped fixtures are exercised. Real-writer quirks (the
  known FSL magic-byte padding inconsistency this RFC already accounts for
  is one example of the kind of thing that only shows up against real
  files) are a real gap until a workspace member tries a real dataset
  against this parser and reports back.

## Unresolved questions

- **Meshing integration.** `SegmentationMask::threshold_hu` is concretely
  typed to `tpt_med_dicom::DicomSeries`. Generalising it (a small trait both
  `DicomSeries` and `NiftiVolume` implement? an adapter that builds a
  `DicomSeries`-shaped view over a `NiftiVolume`?) is real API-design work
  this RFC deliberately does not resolve. Settled by: whoever picks up the
  meshing integration, informed by whether a second non-DICOM volume source
  ever materialises (if `tpt-med-nifti` stays the only one, a narrower
  `DicomSeries`-specific adapter may be simpler than a general trait).
- **`.nii.gz` support.** Needs a vetted pure-Rust inflate dependency, the
  same kind of decision RFC 0001's compressed-transfer-syntax work already
  made three times over for `tpt-med-dicom`. Settled by: whoever picks this
  up, following that precedent (named feature, `cargo deny`-approved
  dependency, off by default).
- **Dual-file `.hdr`/`.img` support.** Lower priority than `.nii.gz` — rarer
  in current tooling — but a real gap if a workspace member's dataset uses
  it. Settled by: whoever hits it first.
