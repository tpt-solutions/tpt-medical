# tpt-med-meshing

Threshold segmentation and voxel-to-hexahedral meshing — the bridge between a
CT Hounsfield volume and a finite-element mesh.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--meshing-orange)](https://crates.io/crates/tpt-med-meshing)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--meshing-blue)](https://docs.rs/tpt-med-meshing)

| | |
|---|---|
| **Layer** | `imaging` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-core`](../../core/tpt-med-core), [`tpt-med-dicom`](../tpt-med-dicom), [`tpt-med-geometry`](../../core/tpt-med-geometry), [`tpt-med-units`](../../core/tpt-med-units) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

To turn a CT scan into something you can solve, you need a mesh. The
industry-standard tools for that — Gmsh, CGAL, Netgen — are GPL or carry
C++ build chains, and neither is a good fit for a WASM build that must run in
a browser.

The **voxel-to-hex** strategy sidesteps meshing entirely: a thresholded
segmentation becomes a structured hexahedral mesh with **one element per solid
voxel**. This is the standard approach for CT-based biomechanical models, it
produces a mesh that is *automatically conforming* (no boolean operations, no
sliver elements, no remeshing), and it is a few hundred lines of Rust.

The honest cost: a stair-step surface. `smooth_mesh` relaxes it for
visualization, but the topology is whatever the CT resolution says it is. That
is a feature for biomechanics (voxel density is meaningful) and a limitation
for aesthetics.

## Features

- **Threshold segmentation** — `SegmentationMask::threshold_hu` with HU
  statistics (`solid_count`, `mean_solid_hu`) for sanity-checking the result.
- **Voxel → hex meshing** — one hex per solid voxel, **shared corner nodes
  compacted** (a 512³ mask does not allocate 8 nodes per voxel).
- **HU-derived per-element materials** — mean HU, apparent density, Young's
  modulus, `BoneRegion`, and Poisson's ratio per element, using the
  `tpt-med-dicom` correlations.
- **Configurable region split** — `region_split_density` (default 1.3 g/cm³,
  ≈300 HU) decides which power law each element uses; `poissons_ratio`
  defaults to 0.30.
- **Laplacian surface smoothing** — `smooth_mesh(mesh, iterations, relaxation)`
  for the stair-step boundary, with the iteration count as an explicit
  parameter so it is never applied accidentally.
- **Versioned CSV export/import** — the human-readable interchange format, with
  a header that states the units.
- Typed `MeshError` (`Io`, `EmptyMask`, `Inconsistent`) — no panics on bad
  input.

## Conventions

- **Coordinates:** patient space, millimetres, taken from the DICOM
  `ImageFrame` — no re-centring, no reorientation.
- **Voxel indexing:** `index = (z·ny + y)·nx + x`.
- **Hex node ordering:** the 8 corners in `(±x, ±y, ±z)` bit order
  `[000, 100, 110, 010, 001, 101, 111, 011]` relative to the voxel's minimum
  corner. This is fixed and shared with the CSV format and the FEM solver.
- **CSV header:** `tpt-medical voxel hex mesh v1`, units mm and MPa. Node
  order and element rows are documented in
  [`docs/book/src/mesh-csv.md`](../../../docs/book/src/mesh-csv.md).

## Usage

```rust
use std::path::Path;
use tpt_med_dicom::DicomSeries;
use tpt_med_meshing::{smooth_mesh, MedicalMesher, SegmentationMask};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. CT series -> HU volume.
    let series = DicomSeries::load_from_dir(Path::new("test-data/dicom/synthetic_ct"))?;

    // 2. Threshold-segment bone. 200 HU is the workspace default; the
    //    literature band is 130-300 HU. Tune per acquisition protocol.
    let mask = SegmentationMask::threshold_hu(&series, 200.0);
    println!("{} solid voxels, mean HU {:?}", mask.solid_count(), mask.mean_solid_hu());

    // 3. Voxels -> hexahedral mesh with HU-derived moduli.
    let mut mesh = MedicalMesher::default().voxels_to_hex_mesh(&mask)?;
    println!("{} nodes, {} elements", mesh.nodes.len(), mesh.elements.len());
    println!("E: mean {:.0} MPa, max {:.0} MPa", mesh.mean_modulus(), mesh.max_modulus());

    // 4. Relax the stair-step surface (optional, for display).
    smooth_mesh(&mut mesh, 3, 0.5);

    // 5. Export the versioned CSV interchange format.
    mesh.write_csv(Path::new("femur_mesh.csv"))?;

    Ok(())
}
```

From a plain voxel model (no DICOM involved) you can still mesh any boolean
mask, which is what the fluid and surgical-planning crates do.

## API Overview

| Item | Purpose |
|---|---|
| `SegmentationMask::threshold_hu(&DicomSeries, min_hu)` | Threshold an HU volume into a solid/void mask |
| `SegmentationMask::{solid_count, mean_solid_hu}` | Segmentation quality indicators — check these before meshing |
| `SegmentationMask::{index, is_solid}` | Bounds-checked voxel access |
| `SegmentationMask::{voxel_center, node_position}` | Voxel centre and node-grid positions in patient space (mm) |
| `MedicalMesher { region_split_density, poissons_ratio }` | Tunable material assignment; defaults 1.3 g/cm³ and 0.30 |
| `MedicalMesher::voxels_to_hex_mesh(&SegmentationMask)` | The conversion; returns `MeshError::EmptyMask` on a void mask |
| `VoxelHexMesh` | `nodes: Vec<Vec3>`, `elements: Vec<[u32; 8]>`, `materials: Vec<ElementMaterial>` |
| `VoxelHexMesh::{mean_modulus, max_modulus, bounds}` | Mesh-level material statistics and extents |
| `VoxelHexMesh::{write_csv, write_csv_to, parse_csv}` | Versioned CSV export/import |
| `ElementMaterial` | `mean_hu`, `density`, `youngs_modulus`, `region`, `poissons_ratio` |
| `smooth_mesh(&mut VoxelHexMesh, iterations, relaxation)` | Laplacian surface relaxation |
| `MeshError` | `Io`, `EmptyMask`, `Inconsistent` |
| `Result<T>` | Crate result alias |

## Verification

- Element/node counts are pinned on a known mask, including the **node-sharing**
  invariant: a solid block of `n³` voxels must produce strictly fewer than
  `8·n³` nodes. This is the test that catches an accidental per-voxel
  allocation.
- Hex corner ordering is asserted element-by-element, because a permuted corner
  order produces a plausible-looking but wrong stiffness matrix — a silent,
  expensive failure.
- CSV round-trip (`write_csv` → `parse_csv`) preserves nodes, connectivity and
  moduli. Mean HU and Poisson's ratio are intentionally *not* exported.
- `mean_solid_hu` is tested so a mis-set threshold (e.g. 0 HU capturing the
  whole image) is obvious from the output.
- `MeshError::EmptyMask` is covered: meshing a void mask is an error, not a
  panic and not an empty mesh.

## Known Limitations

- **The stair-step surface is the topology.** A voxel mesh is exactly as
  detailed as the CT, so a 0.5 mm scan of a femur gives a visibly faceted
  surface. `smooth_mesh` relaxes the nodes for display, but it does not
  recover curvature the scan never measured, and it moves surface nodes away
  from the data the moduli were computed from.
- **Thresholding is binary.** `threshold_hu` produces solid/void with no
  partial-volume weighting, so a voxel that is half bone and half soft tissue
  is fully one or the other. This is a known source of small systematic errors
  in absolute stress, and is usually acceptable because the resulting modulus
  error is small next to the biological variability.
- **No morphological post-processing.** No hole filling, no island removal, no
  connected-component labelling, no closing or opening. A noisy threshold can
  leave specks and voids that then become elements or holes in the FEM mesh.
- **No adaptive or multi-resolution meshing.** Uniform voxel pitch, one element
  per voxel. A structure that needs refinement gets it by acquiring a finer
  scan, not by refining locally.
- **HU → modulus is inherited from `tpt-med-dicom`** and carries that crate's
  limitations: the linear HU→density approximation, and no per-patient QCT
  calibration. There is no override hook for a caller-supplied modulus field.
- **CSV is lossy by design.** Mean HU and Poisson's ratio are not exported, so
  a round-trip through the file is not lossless. That is deliberate — the
  format is a human-readable interchange, not a model store — but it means the
  CSV cannot be the system of record.
- **No surface mesh export.** There is no STL/OBJ/VTK output, so a
  visualisation or CAD workflow outside the web viewer needs the CSV plus
  external tooling.
- **`smooth_mesh` has no validity guard.** A large relaxation factor on a thin
  structure can invert elements, producing a mesh that the FEM solver will
  correctly reject. The failure surfaces in `tpt-med-biomechanics`, not here.

## Related Crates

- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) — supplies the HU volume and the modulus correlations.
- [`tpt-med-biomechanics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-biomechanics) — consumes `VoxelHexMesh` directly; the hex ordering is the contract between them.
- [`tpt-med-wasm`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-wasm) — exposes mesh buffers to WebGL2.
- [`tpt-med-hemodynamics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/fluid/tpt-med-hemodynamics) — `FluidDomain::from_mask` uses the same voxel-mask idea for lumen geometry.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Changing the hex
corner ordering is a breaking change to the CSV format and the FEM contract —
open an RFC first.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.
