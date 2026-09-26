# RFC 0008: Calibration-Phantom Rod Sampling

- **Status:** Accepted
- **Started:** 2026-09-27
- **Crates:** `tpt-med-dicom` (implementation)

## Summary

Add a `phantom` module to `tpt-med-dicom` that turns a CT scan of a
calibration phantom into the `(HU, known_value)` pairs `QctCalibration::fit`
(and, via [RFC 0007](0007-bmd-apparent-density-conversion.md),
`BmdToAshDensity`) already consume, closing the gap those types' own docs
name: nothing in this crate today finds a phantom in a series or samples its
rod ROIs. It deliberately does **not** attempt fully-automatic, manufacturer-
blind phantom detection — the todo item this RFC resolves already names why
that is not one algorithm ("phantom geometry varies by manufacturer... no
single detection algorithm across them"). Instead it splits the problem into
a genuinely universal part (locating a solid phantom's cross-section
centroid by intensity thresholding, which does not depend on rod layout) and
a caller-supplied, per-manufacturer part (rod layout and known values), and
ships v0 as semi-automatic: centroid-located, but rotation- and
slice-selected by the caller.

## Motivation

`QctCalibration::fit` and `BmdToAshDensity::new` (RFC 0007) both take
`(HU, known_value)` points and say nothing about where those points come
from. In practice they come from scanning a calibration phantom — a solid
insert with several parallel rods of known density/BMD — alongside or
instead of the patient, then reading the mean HU inside each rod's
cross-section on the reconstructed image. Today a caller has to do that
sampling by hand, outside the crate, for every scan: find the phantom slice,
find each rod's pixel coordinates, average the HU inside a circular ROI, and
pair it with the rod's datasheet value. That manual step is exactly the kind
of thing this crate exists to make repeatable and auditable instead of ad
hoc.

**Question of interest:** given a CT slice (or slices) containing a
calibration phantom and a description of that phantom's rod layout, what HU
value does each rod read? **Model risk:** low for the sampling step itself
(it is a geometric ROI average, not a biomechanical computation), but it sits
directly upstream of `QctCalibration`/`BmdToAshDensity`, so a wrong rod
correspondence (sampling the wrong rod, or getting two rods' known values
swapped) silently produces a wrong calibration line with the same "looks
fine, is wrong" character RFC 0007 already worried about for BMD conversion.
**Model influence:** the same as `QctCalibration` itself — whatever modulus
or density downstream analysis eventually consumes inherits this step's
correctness.

## Detailed design

### Splitting the universal part from the manufacturer-specific part

A calibration phantom's **outer cross-section** (the water-equivalent
cylinder or slab the rods are embedded in) is, across every manufacturer this
RFC is aware of (Mindways QCT Pro, CIRS/Image Analysis, the European Forearm
Phantom), simply denser than the air and scanner-table padding surrounding it
in the same slice, and close enough to circular or rectangular in cross-
section to have a well-defined centroid. **Locating that centroid does not
need to know the rod layout at all** — it is a thresholded connected-
component centroid, the same primitive `tpt-med-meshing`'s
`SegmentationMask::threshold_hu` already builds on for a completely different
purpose. This is the universal 20% of the problem.

The remaining 80% — how many rods, at what radial offset and angle from the
centroid, and what each one's known value is — is genuinely different per
manufacturer and per phantom model, and sometimes per lot (a datasheet
revision can change stated rod values). There is no way to infer this from
pixel data alone, and guessing at it would be exactly the kind of silent,
unverifiable default RFC 0007 already rejected for BMD coefficients. So this
RFC requires it as caller-supplied, cited data — the same enforcement
pattern.

### Public API (sketch)

```rust
/// One phantom rod: its offset from the phantom's cross-section centroid (in
/// the phantom's own frame, at a stated reference angle), physical radius,
/// and known calibration value.
#[derive(Debug, Clone, Copy)]
pub struct PhantomRod {
    /// Radial distance from centroid, mm.
    pub radial_offset_mm: f64,
    /// Angle from the phantom's reference axis, radians.
    pub angle_rad: f64,
    /// Rod radius, mm (the ROI is sampled at a fraction of this — see below).
    pub radius_mm: f64,
    /// The value `QctCalibration::fit`/`BmdToAshDensity` should pair with
    /// this rod's measured HU (density, BMD, whatever the caller's
    /// downstream calibration expects).
    pub known_value: f64,
}

/// A named, cited rod layout for one phantom model.
///
/// Cannot be constructed without a `source` citation (datasheet, revision,
/// lot if known) — same enforcement as `BmdToAshDensity`/`AshFraction` in
/// RFC 0007, and for the same reason: this is exactly the kind of
/// manufacturer-specific data this crate refuses to bake in silently.
pub struct PhantomModel {
    rods: Vec<PhantomRod>,
    source: String,
}

impl PhantomModel {
    pub fn new(rods: Vec<PhantomRod>, source: impl Into<String>) -> Result<Self>;
    pub fn rods(&self) -> &[PhantomRod];
    pub fn source(&self) -> &str;
}

/// Locates a solid phantom's cross-section centroid in one slice by
/// intensity thresholding — the manufacturer-agnostic half of this RFC.
///
/// `background_max_hu` separates phantom material from surrounding air/
/// table (a phantom's water-equivalent shell reads well above air's
/// −1000 HU; a caller who knows their scan can tighten this past the
/// default). Returns `None` if no connected region above threshold reaches
/// `min_area_px` — i.e. no phantom-sized object was found in this slice.
pub fn locate_phantom_centroid(
    slice: &DicomSlice,
    background_max_hu: f64,
    min_area_px: usize,
) -> Option<(f64, f64)>; // (row, col), sub-pixel centroid

/// Samples one slice's HU under `model`'s rod layout, given the located
/// centroid and the phantom's rotation (0 if the phantom's reference axis
/// is aligned with the image row axis — true by protocol for phantoms
/// scanned with a fiducial notch aligned to the scanner laser, which is how
/// every manufacturer in this RFC's awareness documents positioning).
///
/// Each rod's ROI is a disk of radius `rod.radius_mm * roi_fraction`
/// (`roi_fraction < 1.0`, default 0.8) centred on its computed image-space
/// position — deliberately smaller than the rod's true radius, so partial-
/// volume pixels at the rod's own edge (which read a blend of rod and
/// surrounding phantom material, not the rod's true value) are excluded
/// rather than dragging the mean toward the boundary.
///
/// Returns one `(mean_hu, known_value)` pair per rod, in `model.rods()`
/// order — feed directly into `QctCalibration::fit` or
/// `BmdToAshDensity`-adjacent calibration.
pub fn sample_phantom_rods(
    slice: &DicomSlice,
    model: &PhantomModel,
    centroid: (f64, f64),
    rotation_rad: f64,
    roi_fraction: f64,
) -> Result<Vec<(f64, f64)>>;
```

`sample_phantom_rods` averages several axial slices when the caller passes
more than one (a real phantom scan is usually several slices thick; the rod
is the same physical object across them), reducing per-slice noise — the
function sketch above takes one `slice` for clarity; the actual signature
takes `&[DicomSlice]` and averages per-rod means across them, rejecting the
call if slices disagree on `rows`/`columns`/`pixel_spacing`.

### Rotation is caller-supplied, not detected, in v0

Every phantom this RFC is aware of is positioned in the scanner with a
documented reference orientation (a flat edge or fiducial notch aligned to
the scanner's laser crosshair), which is a scanning-protocol fact, not
something recoverable from the reconstructed image alone without either a
visible fiducial-detection step (itself manufacturer-specific — some
phantoms' fiducials are notches, others are an asymmetric rod pattern) or an
assumption about how the technologist actually positioned it (not
guaranteed, human-operated). v0 asks the caller for `rotation_rad` — `0.0`
covers the documented-default-orientation case, and a caller who scanned
with a known offset (or wants to verify orientation by rotating until sampled
rod HUs best match `model`'s expected ordering — a brute-force search over
`rotation_rad`, cheap since sampling one slice is cheap) can supply their own.
Automatic fiducial-based rotation detection is named as future work (see
Unresolved Questions), not solved here.

### Slice selection is caller-supplied, not detected, in v0

Identifying *which* series or slice range in a study is the phantom scan
(as opposed to patient anatomy) is a study/protocol-level fact — some
protocols scan the phantom as a wholly separate series, others scan it in the
same series as an off-axis insert. Distinguishing these blind, from pixel
data, risks the same "looks confident, is wrong" failure mode as everything
else this RFC and RFC 0007 refuse to guess at (mistaking a real patient
anatomy slice for the phantom would silently corrupt the calibration with
patient tissue HU). v0 requires the caller to pass the already-identified
slice(s) — typically by `SeriesInstanceUID`, already exposed on
`DicomSlice`/`DicomSeries` today.

### Error handling

New `DicomError::Phantom(String)` variant (parallel to the existing
`Calibration` variant) for:

- `PhantomModel::new` with an empty `source` or an empty `rods` list.
- A `PhantomRod` with non-finite or non-positive `radius_mm`, or non-finite
  `radial_offset_mm`/`angle_rad`/`known_value`.
- `sample_phantom_rods` with mismatched slice dimensions/spacing, or a
  computed rod ROI that falls partially or fully outside the slice bounds
  (rather than silently sampling fewer pixels than intended, which would
  quietly bias the mean toward whichever side stayed in-bounds).
- `roi_fraction` outside `(0.0, 1.0]`.

`locate_phantom_centroid` returns `Option`, not `Result` — "no phantom-sized
object in this slice" is an expected, non-exceptional outcome (the caller
handed in a slice, and may hand in slices without a phantom while probing a
series), not a malformed-input error.

### `#![forbid(unsafe_code)]`

Unaffected. Thresholding and centroid computation are plain iteration over
`Vec<f64>`/`Vec<i32>`, matching every other crate in the workspace.

### Alternatives considered

- **Ship a fully-automatic detector that also finds rotation and slice,
  using a small built-in library of known phantom models (by product name)
  and a generic circle/fiducial-fitting algorithm.** Rejected for v0: this is
  the todo item's own stated obstacle ("no single detection algorithm across
  them") — a genuinely different fiducial-recognition problem per
  manufacturer, and building even one of them (say, Mindways) well enough to
  trust without a real scanned phantom to validate against would be
  asserting confidence this RFC cannot back up. Left as explicit future work.
- **Require the caller to supply rod image-space pixel coordinates directly,
  skipping centroid location and the phantom-frame rod layout entirely.**
  Rejected: pushes all the geometry work back onto the caller for every scan
  (recomputing pixel coordinates whenever the phantom is shifted, angled, or
  scanned at different in-plane resolution), when the centroid-plus-rotation
  parameterisation only needs those two numbers per scan and reuses the same
  `PhantomModel` across scans of the same physical phantom.
- **Detect the phantom centroid via a fixed image-space region (e.g. "assume
  centred in the field of view").** Rejected: true for some protocols, false
  for others (an off-axis insert scanned alongside patient anatomy is
  explicitly one of the cases named above), and silently wrong when false is
  worse than requiring one cheap, real thresholding step that works either
  way.

### Drawbacks

- v0 is not "automatic phantom calibration" end to end — a caller still
  supplies the phantom model (cited), the slice(s), and the rotation. This
  RFC removes the per-rod pixel-hunting and manual averaging, not every
  manual step in the workflow.
- `locate_phantom_centroid`'s threshold-and-connected-component approach can
  be fooled by a phantom scanned touching or overlapping another dense
  object in the same slice (a table rail, a second insert) — `min_area_px`
  and `background_max_hu` are blunt controls, not a robust segmentation.
  Acceptable for v0 given the caller already identifies the correct slice
  (a phantom-only or clearly-isolated-insert slice, by construction of how
  phantom scans are actually acquired); revisit if a real workflow hits a
  case where isolation cannot be assumed.
- No visual/diagnostic output (e.g. an overlay image showing where each rod
  was sampled) — a caller who wants to sanity-check the geometry before
  trusting the fit has to do so by comparing `sample_phantom_rods`' returned
  HU values against the phantom's expected HU ordering themselves.

## Verification strategy

- **Code verification:** `locate_phantom_centroid` is checked against a
  synthetic slice (`SyntheticCtBuilder`, already in this crate for other
  tests) with a hand-placed circular high-HU region at a known, non-trivial
  (row, col) and radius — the centroid returned must match the hand-computed
  geometric centroid to sub-pixel tolerance. `sample_phantom_rods` is checked
  against a synthetic slice with several disk regions placed at known
  `PhantomModel` offsets/angles and known HU values — the returned
  `(mean_hu, known_value)` pairs must reproduce those exact HU values (a
  noiseless synthetic case), and a variant with per-pixel noise added must
  still recover the known values within the noise's own bound.
- **Calculation verification:** not applicable — no discretisation or
  iterative solve.
- **Validation:** none against a real phantom scan for the same reason RFC
  0007 validates no specific coefficients — this crate holds no real patient
  or phantom imaging data by policy (synthetic fixtures only). A real-phantom
  validation (scan a real Mindways or CIRS phantom, run
  `sample_phantom_rods`, compare against the datasheet) is explicitly a gap,
  named below.
- **Golden datasets:** none of `test-data/golden/` changes — this is an
  imaging-input utility, not a simulation whose output is golden-checked.
- **What remains unverified:** behaviour against a real scanned phantom
  (partial-volume effects at real rod boundaries, real image noise
  statistics, a real technologist's actual positioning tolerance around the
  documented reference orientation) — only synthetic, geometrically-exact
  fixtures are exercised until someone runs this against a real scan and
  reports back, the same posture RFC 0006 (NIfTI) already took for real-
  writer files.

## Unresolved questions

- **Automatic rotation detection.** Sketched above as a brute-force search
  over `rotation_rad` matching sampled HU ordering against `model`'s expected
  ordering, but not implemented in v0 (ambiguous when a phantom's rods have
  near-identical known values, and untested against real data). Settled by:
  whoever hits a real workflow where the documented reference orientation
  cannot be assumed.
- **A small built-in library of named `PhantomModel`s for common commercial
  phantoms.** Explicitly out of scope for this RFC's mechanism (see
  Alternatives), but once the mechanism exists, adding e.g.
  `PhantomModel::mindways_qct_pro_solid()` behind the same "must cite a
  specific datasheet revision" discipline as a *value*, not a *default*, may
  be worth a small follow-up RFC of its own if several adopters end up
  hand-transcribing the same manufacturer's numbers independently. Settled
  by: whoever has (and can cite) a specific datasheet in hand.
- **Multi-slice averaging weighting.** v0 averages per-rod means uniformly
  across the caller-supplied slices; a caller with slices of varying
  image-noise (e.g. mixed reconstruction kernels in one scan) gets no say in
  weighting. Settled by: whoever hits a real case where uniform averaging is
  measurably worse than a noise-weighted alternative.
