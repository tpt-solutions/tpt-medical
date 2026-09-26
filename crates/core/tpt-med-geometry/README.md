# tpt-med-geometry

Geometry primitives and anatomical coordinate transforms for medical
simulation — `Vec3`, `Mat3`, `Plane`, `Aabb`, analytic eigenvalues, and the
DICOM LPS ↔ RAS conversions.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--geometry-orange)](https://crates.io/crates/tpt-med-geometry)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--geometry-blue)](https://docs.rs/tpt-med-geometry)

| | |
|---|---|
| **Layer** | `core` (leaf — no workspace dependencies) |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | none (`std` only) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Every downstream crate needs vectors, matrices and the LPS/RAS conversion, and
all of them must compile to `wasm32-unknown-unknown` without dragging in a
linear-algebra stack. `nalgebra` would work but is heavy for what this domain
actually requires: **small fixed-size 3×3 maths with an analytic symmetric
eigensolver** (principal stresses) and **anatomical coordinate awareness**
(DICOM patient space vs. research space).

This crate provides exactly that, with `#![forbid(unsafe_code)]`, no external
dependencies, and unit tests on every predicate.

## Features

- `Vec3` — `const fn` constructors (`new`, `splat`, `X`/`Y`/`Z`/`ZERO`), dot and
  cross products, `normalize`, `length`, componentwise `min`/`max`/`abs`.
- `Mat3` — `IDENTITY`, `ZERO`, `from_rows`, `from_cols`, `from_array`,
  `diagonal`, indexing, `mul_vec`, `mul_mat`, `transpose`, `det`, `inverse`,
  `rotation_axis_angle`, `scaling`, `trace`, and
  **`symmetric_eigenvalues` → `[f64; 3]`** (principal stresses).
- `Plane` — from point/normal or three points, signed distance, ray intersection.
- `Aabb` — `point`, `from_points`, `expand`, `union`, `center`, `extents`,
  `contains`, `intersects`.
- `coords` — `CoordinateSystem`, `ImageFrame` (DICOM `ImagePositionPatient` +
  `ImageOrientationPatient`), `lps_to_ras`, `ras_to_lps`, `flip_xy`.
- `EPS` — the shared tolerance used by geometric predicates.

## Conventions

- **Angles are radians.**

## Usage

### Coordinates

```rust
use tpt_med_geometry::{lps_to_ras, ras_to_lps, ImageFrame, Vec3};

fn main() {
    // LPS <-> RAS flips x and y; the map is an involution.
    let p = Vec3::new(10.0, 20.0, 30.0);
    assert_eq!(ras_to_lps(lps_to_ras(p)), p);

    // Rebuild the DICOM image frame.
    let frame = ImageFrame {
        origin: Vec3::new(-100.0, -100.0, -250.0),
        row_dir: Vec3::new(1.0, 0.0, 0.0),
        col_dir: Vec3::new(0.0, 1.0, 0.0),
    };
    assert_eq!(frame.normal(), Vec3::new(0.0, 0.0, 1.0));
    let pos = frame.voxel_position(10, 20, (0.5, 0.5));
    assert!((pos.x - -95.0).abs() < 1e-12);
}
```

### Linear algebra and eigenvalues

```rust
use tpt_med_geometry::{Mat3, Vec3};

fn main() {
    let r = Mat3::rotation_axis_angle(Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2);
    let v = r.mul_vec(Vec3::X);
    assert!(v.x.abs() < 1e-12);
    assert!((v.y - 1.0).abs() < 1e-12);

    // Symmetric eigenvalues are the principal stresses.
    let s = Mat3::from_rows(
        Vec3::new(10.0, 2.0, 0.0),
        Vec3::new(2.0, 6.0, 0.0),
        Vec3::new(0.0, 0.0, 2.0),
    );
    let mut e = s.symmetric_eigenvalues();
    e.sort_by(|a, b| b.partial_cmp(a).unwrap());
    assert!(e[0] > 10.0 && e[2] < 2.0);
}
```

### Primitives

```rust
use tpt_med_geometry::{Aabb, Plane, Vec3};

fn main() {
    let plane = Plane::from_point_normal(Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)).unwrap();
    assert!((plane.signed_distance(Vec3::new(0.0, 0.0, 5.0)) - 5.0).abs() < 1e-12);
    assert!(plane
        .ray_intersection(Vec3::new(0.0, 0.0, 10.0), Vec3::new(0.0, 0.0, -1.0))
        .is_some());

    let b = Aabb::point(Vec3::ZERO).expand(Vec3::new(10.0, 10.0, 10.0));
    assert!(b.contains(Vec3::new(5.0, 5.0, 5.0)));
    assert!(b.intersects(&Aabb::point(Vec3::new(20.0, 20.0, 20.0)).expand(Vec3::X)));
    assert_eq!(b.extents(), Vec3::new(10.0, 10.0, 10.0));
}
```

## API Overview

| Item | Purpose |
|---|---|
| `Vec3` | 3-vector; `new`/`splat` are `const fn`, plus `X`/`Y`/`Z`/`ZERO` |
| `Mat3` | 3×3 matrix; `IDENTITY`, `ZERO`, `from_rows`/`from_cols`/`from_array`/`diagonal`, `at`/`set`/`row`/`col`/`to_array` |
| `Mat3::mul_vec`, `::mul_mat` | Right-handed matrix application |
| `Mat3::det`, `::inverse`, `::transpose`, `::trace` | Standard operations |
| `Mat3::rotation_axis_angle`, `::scaling` | Rigid and diagonal constructions |
| `Mat3::symmetric_eigenvalues -> [f64; 3]` | Analytic principal values (principal stresses) |
| `Plane` | `from_point_normal`, `from_three_points`, `signed_distance`, `ray_intersection` |
| `Aabb` | `point`, `from_points`, `expand`, `union`, `center`, `extents`, `contains`, `intersects` |
| `CoordinateSystem::{Lps, Ras}` | Named anatomical coordinate systems |
| `ImageFrame` | `origin`, `row_dir`, `col_dir`; `normal`, `voxel_position`, `rotation` |
| `lps_to_ras`, `ras_to_lps`, `flip_xy` | Axis flips between patient and research space |
| `EPS: f64` | Shared predicate tolerance (`1.0e-12`) |

## Verification

- Unit tests for every constructor, predicate and transform, including the
  LPS↔RAS round trip and frame normalisation.
- `Mat3::inverse` returns `Option` and is tested on singular input: a
  rank-deficient matrix yields `None` rather than a silent `NaN`.
- `symmetric_eigenvalues` is exercised on the stress tensors produced by
  `tpt-med-biomechanics`, so principal-stress regressions surface here.
- `Aabb::contains`/`intersects` are boundary-inclusive; tests pin that
  behaviour explicitly, because an off-by-one in a viewport fit is invisible
  until a mesh clips.

## Related Crates

- [`tpt-med-units`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-units) — typed quantities; feed millimetres into these unitless types.
- [`tpt-med-core`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-core) — uses `Aabb` for anatomical region extents.
- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) — builds an `ImageFrame` from DICOM position/orientation tags.
- [`tpt-med-tissue`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-tissue) — deformation gradients are `Mat3`.
- [`tpt-med-surgical-planning`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/surgical/tpt-med-surgical-planning) — osteotomy planes and fragment transforms.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Numerical
constants must cite their source in a doc comment.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.

- **Lengths are caller-defined.** The workspace convention is millimetres (see
  [`tpt-med-units`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-units));
  these types are deliberately unitless so they compose with any unit type.
- **Right-handed, column-vector convention:** `y = M · x`, rotations are
  right-handed about `axis`.
- **DICOM patient coordinates are LPS** (+x Left, +y Posterior, +z Superior).
  **Research coordinates (NIfTI, 3D Slicer) are RAS** (+x Right, +y Anterior,
  +z Superior). Anatomical directions follow the ISB recommendations
  (Wu & Cavanagh 1995; Wu et al. 2002).
