# RFC 0007: BMD → Apparent-Density Conversion

- **Status:** Accepted
- **Started:** 2026-09-27
- **Crates:** `tpt-med-dicom` (implementation), `tpt-med-fda` (provenance capture)

## Summary

Add `BmdToApparentDensity`, a small `tpt-med-dicom` type that converts a
clinical QCT phantom's bone-mineral-density (BMD) rod value — the quantity
`QctCalibration::fit` is fit against when a caller uses a real phantom rather
than the screening default — into the apparent (whole-tissue) density that
`HounsfieldMapper::density_to_youngs_modulus` actually expects. It does not
ship a built-in default relation: BMD-to-apparent-density is a two-stage,
literature-published, protocol- and site-dependent conversion, and this RFC's
central design decision is to make supplying that provenance *structurally
required* rather than optional — the caller must name the published relation
they used, and that name is captured wherever the resulting density is
captured, so a regulatory reviewer can see exactly which conversion produced
a given modulus.

## Motivation

`QctCalibration` already lets a caller fit a real HU-to-value line from a
calibration phantom's rods, and its own docs (`hounsfield.rs`) flag the gap
this RFC closes: most clinical QCT phantoms (Mindways QCT Pro, CIRS/Image
Analysis, the European Forearm Phantom) report rod values as BMD — mineral
concentration in mg/cm³, K₂HPO₄- or hydroxyapatite(CaHA)-equivalent — not as
apparent density in g/cm³. Feeding a BMD number into
`density_to_youngs_modulus` as if it were apparent density silently produces
a modulus that is wrong by roughly the ash fraction of bone (on the order of
40–60%, tissue- and site-dependent) — exactly the "confidently wrong number"
the crate's docs already say it refuses to produce by omission. Without this
conversion, every caller who owns a real clinical phantom is stuck on the
screening default they were trying to move past.

**Question of interest:** given a BMD value from a named phantom convention,
what apparent (whole-tissue) density does it correspond to, for use in a
density-to-modulus power law? **Model risk:** moderate — this is a material
property input to a subject-specific FEM model (femur/lumbar-spine loading,
implant sizing), not the FEM solve itself, but a wrong density propagates
directly into a wrong modulus via `density_to_youngs_modulus`'s power law,
and errors compound (the power-law exponents are 1.49–2.0). **Model
influence:** high for any workflow that adopts a real QCT phantom instead of
the screening default specifically *because* it wants a more trustworthy
modulus — a silently wrong conversion at this step defeats the reason the
caller reached for a calibrated phantom in the first place.

## Detailed design

### Why two stages, not one

The literature does not publish a single BMD→apparent-density number; it
publishes two separate, independently-varying relations, each from a
different body of work:

1. **BMD (phantom-convention mg/cm³) → ash density (g/cm³).** A linear
   relation specific to the phantom's mineral-equivalent convention
   (K₂HPO₄ vs. CaHA read differently for the same physical rod) and to the
   scanner/protocol the calibration was performed under.
2. **Ash density (g/cm³) → apparent (wet, whole-tissue) density (g/cm³).**
   A separate relation via the tissue's ash fraction — the mass fraction of
   dry bone that is mineral, roughly 0.6 for human cortical/trabecular bone
   but reported across a real range (site, species-model, and demineralised
   vs. fresh preparation all move it) rather than as one universal constant.

Collapsing these into one composite slope/intercept would hide which half of
the conversion a caller is trusting, and would make it impossible to reuse
stage 1's output (ash density is itself the input several published
density-to-modulus laws use directly, distinct from this crate's
apparent-density-based Morgan–Keaveny power laws) without redoing the whole
fit. Keeping them separate also matches how the source literature itself is
structured: a caller citing "BMD→ash" and "ash→apparent" is citing two
different papers, not paraphrasing one.

### No built-in default relation

This is the RFC's load-bearing decision, and follows directly from
`QctCalibration`'s own precedent of refusing an assumed default once a real
phantom is in play. `BmdToApparentDensity` (see API below) ships **zero**
named presets (no `BmdToApparentDensity::mindways()`,
no `BmdToApparentDensity::schileo_2008()`) in v0. Every relation is
constructed from caller-supplied coefficients plus a **mandatory, non-empty
citation string** — the type cannot be constructed without one. Reasons:

- The exact published coefficients for any given relation need to be
  transcribed from the primary source (paper, table, and often a specific
  scanner/protocol column within that table) with enough care that this RFC
  is not the place to assert them from memory. Baking in a specific number
  under a name like `schileo_2008()` without that transcription being
  checked against the actual paper would be exactly the silent, unverifiable
  default this design exists to avoid — worse than the gap it closes, because
  it would *look* authoritative.
- Real QCT protocols already vary per scanner and per phantom lot in ways a
  single hard-coded constant cannot track; a site running its own V&V is
  expected to have (or obtain) its own validated coefficients regardless of
  what this crate ships.
- The pattern generalises the same way `QctCalibration::fit` already does:
  ship the *mechanism*, require the *data*, and make the omission a compile-
  or run-time error rather than a wrong number that compiles silently.

### Public API (sketch)

```rust
/// Which mineral-equivalent convention a phantom's BMD value uses. Two rods
/// of the same physical composition read different numbers under the two
/// conventions, so this is not cosmetic — it selects which published
/// BMD→ash relation is valid to apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BmdConvention {
    /// K₂HPO₄-equivalent (e.g. Mindways QCT Pro solid phantom).
    K2Hpo4Equivalent,
    /// Hydroxyapatite/CaHA-equivalent (e.g. CIRS, European Forearm Phantom).
    HydroxyapatiteEquivalent,
}

/// A published (or site-fitted) linear relation from a phantom's BMD
/// (mg/cm³, in `convention`'s units) to ash density (g/cm³).
///
/// Cannot be constructed without a `source` citation — see this RFC's
/// "No built-in default relation" for why that is enforced, not advisory.
pub struct BmdToAshDensity {
    convention: BmdConvention,
    slope: f64,
    intercept: f64,
    source: String,
}

impl BmdToAshDensity {
    /// `slope`/`intercept` in `ash_density_g_cm3 = intercept + slope * bmd_mg_cm3`.
    /// `source` must be non-empty (paper/table/protocol citation) — enforced,
    /// see `DicomError::Calibration`.
    pub fn new(convention: BmdConvention, slope: f64, intercept: f64, source: impl Into<String>) -> Result<Self>;

    pub fn convention(&self) -> BmdConvention;
    pub fn source(&self) -> &str;
    pub fn ash_density(&self, bmd_mg_cm3: f64) -> Density; // g/cm3, an "ash density" newtype role reusing tpt_med_units::Density
}

/// A published (or site-fitted) ash fraction: the dry-mass fraction of bone
/// that is mineral, converting ash density to apparent (whole-tissue)
/// density: `apparent = ash / ash_fraction`.
///
/// Also requires a `source` citation, for the same reason.
pub struct AshFraction {
    fraction: f64, // (0, 1]
    source: String,
}

impl AshFraction {
    pub fn new(fraction: f64, source: impl Into<String>) -> Result<Self>;
    pub fn source(&self) -> &str;
    pub fn apparent_density(&self, ash_density: Density) -> Density;
}

/// The composed two-stage conversion, carrying both citations together so a
/// caller (and, downstream, an audit trail) has one object naming the full
/// provenance of a BMD-derived apparent density.
pub struct BmdToApparentDensity {
    bmd_to_ash: BmdToAshDensity,
    ash_fraction: AshFraction,
}

impl BmdToApparentDensity {
    pub fn new(bmd_to_ash: BmdToAshDensity, ash_fraction: AshFraction) -> Self;
    pub fn apparent_density(&self, bmd_mg_cm3: f64) -> Density;
    /// Both citations, concatenated for a report/audit line.
    pub fn provenance(&self) -> String;
}
```

`HounsfieldMapper::density_to_youngs_modulus` is unchanged — it keeps taking
a `Density` regardless of how the caller obtained it. The new convenience
would be:

```rust
impl HounsfieldMapper {
    pub fn bmd_to_youngs_modulus(
        bmd_mg_cm3: f64,
        conversion: &BmdToApparentDensity,
        region: BoneRegion,
    ) -> Modulus {
        Self::density_to_youngs_modulus(conversion.apparent_density(bmd_mg_cm3), region)
    }
}
```

### Units and conventions

- BMD input is mg/cm³ (the near-universal clinical QCT reporting unit for
  both conventions); ash and apparent density are g/cm³, matching
  `tpt_med_units::Density`'s existing convention (`Density::from_gcm3`).
- `BmdConvention` is metadata carried alongside the fit, not itself a unit
  conversion — `BmdToAshDensity` does not attempt to convert between the two
  conventions (K₂HPO₄-equivalent ↔ CaHA-equivalent), since that conversion is
  its own published, phantom-batch-specific relation and out of scope here
  (a caller who needs it supplies their own `BmdToAshDensity` fitted for
  whichever convention their phantom actually reports).
- No sign or coordinate conventions are involved; this is a scalar value
  conversion, not a geometric one.

### Error handling

Extends `DicomError::Calibration` (already used by `QctCalibration::fit`)
rather than adding a new variant — both are "a calibration input failed a
sanity check," and the existing variant already carries a free-form message:

- `BmdToAshDensity::new` / `AshFraction::new` reject an empty or
  whitespace-only `source`.
- `AshFraction::new` rejects `fraction <= 0.0` or `fraction > 1.0` (an ash
  fraction is a mass fraction; it cannot be zero, negative, or exceed 1) and
  non-finite input.
- `BmdToAshDensity::new` rejects non-finite `slope`/`intercept`, mirroring
  `QctCalibration::fit`'s existing finiteness check.

### `#![forbid(unsafe_code)]`

Unaffected — this is pure arithmetic on `f64`/`String`, no new unsafe
surface, matching every other crate in the workspace.

### Alternatives considered

- **Ship named presets for the well-known phantoms/papers anyway, clearly
  labelled "unverified — check against source before regulatory use."**
  Rejected: a label reduces but does not remove the risk that a caller copies
  the constant and never does check it, and the whole point of this design is
  to make the citation an unavoidable, structural part of using the type, not
  an asterisk on a number that already looks final.
- **Fold the two stages into one caller-supplied linear relation
  (BMD→apparent directly), skipping the ash-density intermediate.** Rejected:
  loses the ability to reuse the ash-density intermediate for a future
  ash-density-based modulus law, and obscures that a caller is actually
  trusting two independently-sourced relations, not one.
- **Make `source` an optional field with a default of `"unspecified"`.**
  Rejected: optional-with-a-default is exactly the shape that lets the
  citation silently not happen; making construction fail without one is the
  entire mechanism this RFC is proposing.

### Drawbacks

- Two extra small types (`BmdToAshDensity`, `AshFraction`) plus a composing
  wrapper (`BmdToApparentDensity`) is more ceremony than a single conversion
  function, for what is arithmetically two multiply-adds. Accepted as the
  direct cost of making the citation mandatory rather than a doc-comment
  suggestion.
- Shipping no presets means every adopter's first use of this type requires
  them to already have (or go find) a specific published relation for their
  phantom and protocol — this RFC closes the "no mechanism exists" gap but
  deliberately does not close the "which paper do I cite" research gap for
  any given site.
- `provenance()` returning a concatenated free-text string is a weak
  machine-readable contract (fine for a human-read audit line, not for
  automated cross-checking that a given citation was, say, actually about the
  right phantom model). Revisit if `tpt-med-fda` ever wants to validate
  citations structurally rather than just record them.

## Verification strategy

- **Code verification:** the conversion is two linear maps composed; verified
  by hand-worked test cases (e.g. a slope/intercept pair chosen so the
  expected ash density and apparent density are round numbers) rather than
  against any external tool, the same posture `QctCalibration`'s own tests
  already take for its OLS fit.
- **Calculation verification:** not applicable — no discretisation or
  iterative solve, a direct evaluation.
- **Validation:** deliberately **none** in this crate for the coefficients
  themselves — validating a specific `BmdToAshDensity`/`AshFraction` pair
  against real phantom data is the adopting site's V&V responsibility, not
  this type's. What this crate validates is that the *mechanism* behaves
  correctly for the arithmetic it performs (see code verification above) and
  that it refuses missing/invalid inputs (empty citation, out-of-range
  fraction, non-finite coefficients).
- **Golden datasets:** none of `test-data/golden/` changes — this is an input
  calibration utility, not a simulation whose output is golden-checked.
- **What remains unverified:** every real published BMD→ash and ash→apparent
  relation a caller might supply. This RFC verifies the conversion mechanism
  works and enforces provenance; it explicitly does not, and cannot, verify
  that any given site's cited numbers are correctly transcribed from their
  source or correctly applicable to their scanner/protocol — that is exactly
  the site's own V&V, which `provenance()` exists to make auditable rather
  than to replace.

## Unresolved questions

- **Should `provenance()` (or the two source strings) be wired into
  `tpt-med-fda`'s `ReproducibilityManifest`/audit trail automatically, or
  left for the caller to record explicitly?** Leaning toward explicit for
  v0 (mirroring how `QctCalibration` itself is not auto-captured today), but
  worth revisiting once a real workflow using this type exists to show
  whether the manual step gets skipped in practice. Settled by: whoever
  wires the first real BMD-based workflow through `tpt-med-fda`.
- **Does a single ash fraction (`AshFraction`) suffice, or does site
  (cortical vs. trabecular) need two separate fractions the way
  `BoneRegion` already splits the modulus power law?** The literature
  reports somewhat different ash fractions for cortical vs. trabecular bone.
  v0's sketch takes one `AshFraction` per `BmdToApparentDensity`, requiring a
  caller who wants per-region fractions to build two `BmdToApparentDensity`
  instances. Settled by: whoever implements this, informed by whether a real
  use case needs both regions from one phantom fit simultaneously.
