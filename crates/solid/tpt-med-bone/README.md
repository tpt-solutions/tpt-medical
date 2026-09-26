# tpt-med-bone

Linear elastic bone mechanics, HU-based property assignment, and Wolff's-law
density remodeling.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--bone-orange)](https://crates.io/crates/tpt-med-bone)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--bone-blue)](https://docs.rs/tpt-med-bone)

| | |
|---|---|
| **Layer** | `solid` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-core`](../../core/tpt-med-core), [`tpt-med-dicom`](../../imaging/tpt-med-dicom), [`tpt-med-geometry`](../../core/tpt-med-geometry), [`tpt-med-units`](../../core/tpt-med-units) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Bone is not a homogeneous isotropic material. Cortical and trabecular bone
differ by an order of magnitude in modulus, and both differ from each other
*within the same bone* — and along different directions. A CT scan resolves
density, and density maps to modulus through published correlations, but the
correlations are region-specific and the resulting material is
anisotropic along the principal loading axis.

This crate is the layer that turns a `BoneType` and a Hounsfield value into a
material a solver can use, and it provides the time-domain remodeling law that
stress-shielding assessment needs.

The **HU → modulus correlations themselves live in
[`tpt-med-dicom`](../tpt-med-dicom/README.md)** (`HounsfieldMapper`); this
crate adds the structural description, the anisotropy classification, and the
remodeling dynamics on top.

## Features

- **`BoneMaterial`** — Young's modulus, Poisson's ratio, yield and ultimate
  stress, tissue class, and an anisotropy model.
- **Reference materials** — `cortical_reference(bone_type)` and
  `trabecular_reference(bone_type)` for literature-typical values per bone.
- **HU-based construction** — `from_hu(hu, bone_type, poissons_ratio)`
  classifies cortical vs. trabecular at 1.3 g/cm³ apparent density (≈300 HU)
  and applies the matching power law automatically.
- **Anisotropy** — `Isotropic`, `TransverselyIsotropic { axis, e_long,
  e_trans }` (e.g. along the femoral shaft) and `Orthotropic { axes, e }`,
  with `effective_modulus(direction)` doing a proper directional projection
  rather than returning a scalar.
- **Wolff's-law remodeling** — `BoneRemodelingModel` with an explicit
  **lazy zone** (the mechanostat's dead band) around the reference stimulus,
  separate apposition and resorption rates, and a viability clamp.
- **Three stimulus measures** — strain energy density, principal strain, and
  damage accumulation rate, so the remodeler can be driven by whatever the
  solver actually produces.
- Typed `Density`/`Modulus` from `tpt-med-units`, so the working unit set
  cannot drift.

## Conventions

- Moduli and stresses in **MPa**; density in **g/cm³**.
- Reference modulus: `reference_modulus(TissueClass)`; cortical is the stiff
  branch, trabecular the compliant one.
- Remodeling response is linear in the *normalized* over/under-stimulus,
  `Δρ = rate · (S/S_ref − 1)`, where `S_ref` is the nearer lazy-zone bound,
  capped at a 3× overshoot and clamped to the viable range.

## Usage

```rust
use tpt_med_bone::{BoneMaterial, BoneRemodelingModel, RemodelingStimulus, TissueClass};
use tpt_med_core::BoneType;
use tpt_med_geometry::Vec3;
use tpt_med_units::Density;

fn main() {
    // HU -> material, with automatic cortical/trabecular classification.
    let cortical = BoneMaterial::from_hu(1000.0, BoneType::Femur, 0.30);
    assert_eq!(cortical.tissue_class, TissueClass::Cortical);

    let trabecular = BoneMaterial::from_hu(300.0, BoneType::Femur, 0.30);
    assert_eq!(trabecular.tissue_class, TissueClass::Trabecular);
    assert!(cortical.youngs_modulus > trabecular.youngs_modulus);

    // Directional stiffness along the femoral shaft axis.
    let anisotropic = BoneMaterial {
        anisotropy: tpt_med_bone::Anisotropy::TransverselyIsotropic {
            axis: Vec3::new(0.0, 0.0, 1.0),
            e_long: 17_500.0,
            e_trans: 11_500.0,
        },
        ..BoneMaterial::cortical_reference(BoneType::Femur)
    };
    let along = anisotropic.effective_modulus(Vec3::new(0.0, 0.0, 1.0));
    let across = anisotropic.effective_modulus(Vec3::new(1.0, 0.0, 0.0));
    assert!((along - 17_500.0).abs() < 1e-6);
    assert!((across - 11_500.0).abs() < 1e-6);

    // Wolff's Law: a year of over-stressed bone apposes density.
    let model = BoneRemodelingModel {
        stimulus: RemodelingStimulus::StrainEnergyDensity,
        ..Default::default()
    };
    let start = Density::from_gcm3(1.0);
    let after = model.simulate_days(start, 0.008 /* above lazy zone */, 365);
    assert!(after.to_gcm3() > start.to_gcm3());

    // ... and an under-stressed (stress-shielded) region loses bone.
    let shielded = model.simulate_days(start, 0.001, 365);
    assert!(shielded.to_gcm3() < start.to_gcm3());
}
```

## API Overview

| Item | Purpose |
|---|---|
| `TissueClass` | `Cortical`, `Trabecular` |
| `Anisotropy` | `Isotropic`, `TransverselyIsotropic { axis, e_long, e_trans }`, `Orthotropic { axes, e }` |
| `BoneMaterial` | `bone_type`, `tissue_class`, `youngs_modulus`, `poissons_ratio`, `yield_stress`, `ultimate_stress`, `anisotropy` |
| `BoneMaterial::cortical_reference(BoneType)` | Literature-typical cortical material for a bone |
| `BoneMaterial::trabecular_reference(BoneType)` | Literature-typical trabecular material |
| `BoneMaterial::from_hu(hu, BoneType, poissons_ratio)` | HU → classified material via the `tpt-med-dicom` correlations (split at 1.3 g/cm³) |
| `BoneMaterial::effective_modulus(direction) -> f64` | Directional modulus; falls back to `youngs_modulus` when isotropic |
| `RemodelingStimulus` | `StrainEnergyDensity`, `PrincipalStrain`, `DamageAccumulation` |
| `BoneRemodelingModel` | `stimulus`, `apposition_rate`, `resorption_rate`, `lazy_zone`, `viable_density`; `Default` = Frost screening values |
| `BoneRemodelingModel::update_density(current, stimulus, dt_days)` | One time step, clamped to the viable range |
| `BoneRemodelingModel::simulate_days(start, stimulus, days)` | Multi-day forward simulation |
| `clamp_viable(density, (low, high))` | Clamp a density into the physiologically viable window |
| `reference_modulus(TissueClass) -> Modulus` | Reference modulus for a tissue class |

## Verification

- **HU → modulus** is locked to the `tpt-med-dicom` power laws, and the
  cortical/trabecular split is asserted on both sides of 1.3 g/cm³ (≈300 HU) —
  including exactly at the boundary.
- **Reference materials** are asserted to be ordered cortical > trabecular for
  every `BoneType`, which catches a transposed constant immediately.
- **`effective_modulus`** is pinned at the symmetry axis, perpendicular to it,
  and for the isotropic fallback, plus a rotational invariance check (the
  modulus must depend on the orientation of the direction vector, not its
  magnitude).
- **Remodeling** is verified against the three regimes the model encodes:
  inside the lazy zone → no change; above → apposition; below → resorption.
  Over/under-shoot asymmetry (`apposition_rate` ≠ `resorption_rate`) is
  asserted, and the viability clamp is tested to bound runaway apposition.
- Monotonicity: `update_density` is a non-decreasing function of the stimulus.
- Golden reference datasets: `test-data/golden/solid/` and
  `test-data/golden/devices/hip_stem_micromotion.json`.

## Known Limitations

- Reference constants are literature-typical population values, not
  patient-calibrated. Quantitative CT with a phantom, or mechanical
  indentation data, is required for a patient-specific claim.
- Remodeling is a lumped, spatially uniform law per call — it does not
  distribute density across a mesh. Driving per-element remodeling from a
  solved strain energy density field is a caller-side loop.
- No explicit disuse or resorption-deadline model; the lazy zone is the only
  memory.
- Default parameters are **Frost-style mechanostat screening values**:
  reference SED stimulus ≈0.004 mJ/mm³ with a ±35 % lazy zone
  `(0.0026, 0.0054)`, apposition 0.003 and resorption 0.002 g/cm³/day.
  These are *starting points for screening, not patient-calibrated values*.

## Related Crates

- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) — owns the HU→density→modulus correlations.
- [`tpt-med-biomechanics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-biomechanics) — consumes `youngs_modulus` and `poissons_ratio` per element.
- [`tpt-med-orthopedics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-orthopedics) — stress-shielding assessment, the main consumer of remodeling.
- [`tpt-med-core`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-core) — the `BoneType` taxonomy.
- [`tpt-med-units`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-units) — the `Density` and `Modulus` types.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Every published
constant must cite its source in a doc comment. New material models require an
[RFC](../../../rfcs).

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.
