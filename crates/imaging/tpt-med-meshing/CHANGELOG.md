# Changelog

All notable changes to `tpt-med-meshing` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README documenting the voxel-to-hex strategy, its trade-offs (a
  stair-step surface that `smooth_mesh` relaxes for display but that is
  meaningful for biomechanics), and the shared hex corner-ordering contract
  with `tpt-med-biomechanics`.
- **`SegmentationMask::threshold_nifti(&NiftiVolume, min_hu)`** — the
  `tpt-med-nifti` integration RFC 0006 deliberately left open. Same
  threshold rule and same default as `threshold_hu`; NIfTI's native RAS
  origin and direction columns are converted through `ras_to_lps` into the
  LPS patient frame `threshold_hu` produces, so a mask (and the mesh built
  from it) is interchangeable regardless of which format the volume came
  from. Settled as a narrow direct constructor, not a shared volume trait,
  because `NiftiVolume` is currently the only non-DICOM source — a trait
  would be generalisation ahead of a second consumer, which RFC 0006's
  Unresolved Question explicitly warned against. New dependency:
  `tpt-med-nifti` (zero-dependency itself; the `gzip` feature stays off).

### Planned
- Per-voxel material overrides, so a caller can supply a QCT-calibrated or
  region-specific modulus instead of the default HU correlation.
- Optional node deduplication across disconnected components.

### Notes
- The **hex corner ordering** `[000, 100, 110, 010, 001, 101, 111, 011]` is a
  contract shared with `tpt-med-biomechanics` and the CSV format. Changing it
  is semver-major and requires an RFC: a permuted corner order produces a
  plausible-looking but wrong stiffness matrix, and the failure is silent.
- `region_split_density` defaults to 1.3 g/cm³ (≈300 HU). Changing the default
  changes every mesh modulus, so it is treated as a semver-minor change with a
  V&V re-run.

## [0.1.0] - 2026-09-22

### Added
- `SegmentationMask::threshold_hu(&DicomSeries, min_hu)` — threshold
  segmentation into a solid/void mask, with `solid_count` and `mean_solid_hu`
  as quality indicators, plus bounds-checked `index`/`is_solid` and
  `voxel_center`/`node_position` in patient space (mm).
- `MedicalMesher::voxels_to_hex_mesh(&SegmentationMask)` — one hexahedral
  element per solid voxel with **shared corner nodes compacted**.
- `ElementMaterial` per element: `mean_hu`, `density`, `youngs_modulus`,
  `region` and `poissons_ratio`, derived from the `tpt-med-dicom`
  correlations.
- `MedicalMesher` configuration: `region_split_density` (default 1.3 g/cm³,
  the common CT-FEM cortical/trabecular heuristic) and `poissons_ratio`
  (default 0.30).
- `VoxelHexMesh` with `nodes`, `elements` and `materials`, plus `mean_modulus`,
  `max_modulus` and `bounds`.
- **Versioned CSV interchange format** — `write_csv`, `write_csv_to` and
  `parse_csv`, headed `tpt-medical voxel hex mesh v1` with units stated in the
  header. Mean HU and Poisson's ratio are intentionally not exported.
- `smooth_mesh(&mut VoxelHexMesh, iterations, relaxation)` — Laplacian surface
  relaxation, with the iteration count as an explicit parameter so it is never
  applied accidentally.
- `MeshError` (`Io`, `EmptyMask`, `Inconsistent`) and a `Result<T>` alias; no
  panics on bad input.

### Verification
- Element and node counts are pinned on a known mask, including the
  **node-sharing** invariant: a solid block of `n³` voxels must produce
  strictly fewer than `8·n³` nodes.
- Hex corner ordering is asserted element by element, because a permuted
  ordering yields a wrong stiffness matrix silently.
- CSV round-trip preserves nodes, connectivity and moduli.
- Meshing a void mask returns `MeshError::EmptyMask` rather than panicking or
  returning an empty mesh.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
