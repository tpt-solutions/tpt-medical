# tpt-med-benches

Workspace benchmark suite for the `tpt-medical` stack — harness-free timings
for the three hot paths: DICOM parsing and meshing, hex FEM, and CFD.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--benches-orange)](https://crates.io/crates/tpt-med-benches)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--benches-blue)](https://docs.rs/tpt-med-benches)

| | |
|---|---|
| **Layer** | `applications` |
| **Status** | Alpha, `0.1.0` — `publish = false` |
| **Harness** | `false` (custom `main`, prints timings) |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | `tpt-med-dicom`, `-meshing`, `-biomechanics`, `-tissue`, `-hemodynamics`, `-geometry`, `-units` |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

A performance regression that nobody notices until a user reports "it got
slow" is a defect that costs more to fix than one CI caught. Three benchmarks
covering the three hot paths — parsing and meshing, the FEM inner loop, and the
CFD projection — make the cost of a change visible at review time.

The design choice that matters is **`harness = false`**. These are *pipeline
stage timings* ("how long does thresholding 48×48×24 take?"), not statistical
micro-benchmarks of a single function. libtest's regression model would be
measuring the wrong thing, and a printed table is directly comparable between
commits and readable in a CI log without tooling.

## What this is

Three benchmark programs covering the three computational hot paths in the
workspace, so a performance regression is visible in CI rather than
discovered by a user.

They use `harness = false` and a plain `main` with `std::time::Instant` rather
than `#[bench]` and libtest. That is a deliberate choice: these are
**pipeline stage timings** ("how long does thresholding 48×48×24 take?"), not
statistical micro-benchmarks of a single function. libtest's regression model
would be measuring the wrong thing, and the printed table is directly
comparable across runs and readable in a CI log.

CI runs these on every push and archives the output, so a change in
performance is a visible diff rather than a memory.

## The benchmarks

| Benchmark | Measures |
|---|---|
| `dicom-parsing` | On a 48×48×24 femur phantom held **in memory** (no disk I/O in the measurement path): parsing 24 slices, threshold segmentation, voxel-to-hex meshing, and Laplacian smoothing ×5 |
| `hyperelastic-fem` | Constitutive stress evaluation (the `first_piola_numerical` central-difference path) plus hex-element stiffness assembly and a full CSR/CG solve on a synthetic grid mesh |
| `hemodynamics-cfd` | Per-step projection cost on a 24×12×12 cylinder with 250 SOR sweeps per step, Newtonian blood; reports inlet/outlet flow as a mass-conservation check alongside the timing |

## Features

- Three benchmarks, one per hot path, each printing a labelled table of
  per-stage timings.
- Results are consumed, not optimised away: every benchmark accumulates into a
  `sink` that is printed.
- No disk I/O inside a measurement path — the phantom is built in memory — so
  a timing reflects computation rather than the filesystem.
- Warm-up iterations before measurement where it matters.
- Output archived as a CI artifact, so history is reviewable across commits.
- Synthetic data only.

## Conventions

- Timings are printed with a `{:>10.2?}` human duration, never a raw
  nanosecond count.
- Problem sizes are chosen to finish in seconds on a shared CI runner. They are
  **not** production sizes, and a large-mesh claim needs its own benchmark.
- Mass conservation is printed alongside the CFD timing, so a solver regression
  appears as a flow divergence rather than only as a timing change.
- `cargo bench` uses the `bench` profile; do not compare these numbers against
  a `cargo run` timing.
- When adding a benchmark, keep the harness-free style: pipeline stage timings,
  a printed table, a consumed `sink`. A mixed suite is harder to read than a
  consistent one.

## Usage

```console
# All three.
cargo bench -p tpt-med-benches

# One at a time.
cargo bench -p tpt-med-benches --bench dicom-parsing
cargo bench -p tpt-med-benches --bench hyperelastic-fem
cargo bench -p tpt-med-benches --bench hemodynamics-cfd
```

Each prints a small labelled table, for example:

```text
dicom-parsing benchmark (48x48x24 phantom):
  parse 24 slices                    ...
  threshold segmentation              ...
  voxel-to-hex meshing               ...
  laplacian smoothing x5             ...
```

## Design Notes

- **No disk I/O in the measurement path.** `dicom-parsing` builds the phantom
  in memory and parses the bytes directly, so the number reflects parsing
  rather than the filesystem.
- **Warm-up before measurement.** `hemodynamics-cfd` runs 10 steps before
  timing 100, so branch prediction and cache state are representative.
- **Results are consumed, not optimised away.** Each benchmark accumulates
  into a `sink` that is printed, so the optimiser cannot delete the work.
- **Synthetic data only.** Nothing here touches real patient data.

## API Overview

This package is `publish = false` and exposes **no library API** — only the
three bench targets above. It exists to measure the library crates, not to be
measured itself.

## Verification

- All three benchmarks are compiled and run in CI, so a change that breaks the
  build or panics in a hot path fails the pipeline.
- `hemodynamics-cfd` asserts mass conservation implicitly by printing inlet
  and outlet flow; a large divergence indicates a solver regression, not just
  a timing change.
- Timing output is archived as a CI artifact, so performance history is
  reviewable across commits.

## Known Limitations

- **Absolute timings are machine-dependent.** Compare within a machine, not
  across CI runners, which are shared and noisy.
- **Small problem sizes.** The meshes and domains here are sized to run in
  seconds in CI, not to represent production CT volumes. A large-mesh
  performance claim needs its own benchmark.
- **No statistical regression detection.** The output is a printed table for
  human comparison; there is no automatic threshold that fails the build.
- **Single-threaded measurement.** The solvers are single-threaded today, so
  the numbers do not say anything about parallel scaling.
- `cargo bench` uses the `bench` profile; make sure you compare like with like
  when comparing against a `cargo run` timing.

## Related Crates

- [`tpt-med-examples`](https://github.com/tpt-solutions/tpt-medical/tree/master/examples) — the correctness counterpart to these performance binaries.
- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) and [`tpt-med-meshing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-meshing) — measured by `dicom-parsing`.
- [`tpt-med-biomechanics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-biomechanics) and [`tpt-med-tissue`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-tissue) — measured by `hyperelastic-fem`.
- [`tpt-med-hemodynamics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/fluid/tpt-med-hemodynamics) — measured by `hemodynamics-cfd`.

## Contributing

See the workspace [CONTRIBUTING.md](../../CONTRIBUTING.md). When adding a
benchmark, keep the harness-free style (pipeline stage timings, a printed
table, a consumed `sink`) rather than switching to `#[bench]`; a mixed suite
is harder to read than a consistent one. **Never commit real patient data.**

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.
