# Changelog

All notable changes to `tpt-med-benches` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README documenting why the suite is **harness-free**: these are
  pipeline stage timings ("how long does thresholding 48×48×24 take?"), not
  statistical micro-benchmarks of a single function. libtest's regression model
  would measure the wrong thing, and a printed table is directly comparable
  across runs and readable in a CI log.

### Planned
- Automatic regression detection with a per-stage threshold, rather than
  relying on a human reading a diff.
- A large-mesh benchmark sized to production CT volumes, since the current
  problem sizes are chosen to run in seconds in CI.
- A parallel-scaling benchmark, once the solvers are multi-threaded.

### Notes
- `publish = false`. This package exists only to measure the library crates.
- Absolute timings are machine dependent: compare within a machine, not across
  shared CI runners. `cargo bench` uses the `bench` profile, so do not compare
  against a `cargo run` timing.
- When adding a benchmark, keep the harness-free style — pipeline stage
  timings, a printed table, a consumed `sink`. A mixed suite is harder to read
  than a consistent one.

## [0.1.0] - 2026-09-22

### Added
- `dicom-parsing` — on a 48×48×24 femur phantom held **in memory** (no disk
  I/O in the measurement path): parsing 24 slices, threshold segmentation,
  voxel-to-hex meshing, and Laplacian smoothing ×5, each timed separately.
- `hyperelastic-fem` — constitutive stress evaluation (the
  `first_piola_numerical` central-difference path, 2000 iterations) plus hex
  element stiffness assembly and a full CSR/CG solve on a synthetic grid mesh
  built in-benchmark.
- `hemodynamics-cfd` — per-step projection cost on a 24×12×12 cylinder with
  250 SOR sweeps per step and Newtonian blood, after a 10-step warm-up and over
  100 measured steps. Reports inlet and outlet flow alongside the timing, so a
  mass-conservation divergence is visible as a solver regression rather than
  only as a timing change.
- All three use `harness = false` with a plain `main` and `std::time::Instant`,
  printing a labelled table. Each accumulates into a `sink` that is printed, so
  the optimiser cannot delete the work.
- No library API; bench targets only.

### Verification
- All three benchmarks are compiled and run in CI, so a change that breaks the
  build or panics in a hot path fails the pipeline.
- Timing output is archived as a CI artifact, so performance history is
  reviewable across commits rather than only visible in the run log.
- Synthetic data only; nothing here touches real patient data.

### Known limitations
- Small problem sizes, chosen to run in seconds in CI; a large-mesh performance
  claim needs its own benchmark.
- No statistical regression detection — the output is a printed table for
  human comparison, with no automatic failing threshold.
- Single-threaded measurement, so the numbers say nothing about parallel
  scaling.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
