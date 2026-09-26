# tpt-med-examples

Runnable end-to-end example binaries for the `tpt-medical` stack — the fastest
way to see what the whole pipeline does, and the reference for how the crates
are meant to be composed.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--examples-orange)](https://crates.io/crates/tpt-med-examples)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--examples-blue)](https://docs.rs/tpt-med-examples)

| | |
|---|---|
| **Layer** | `applications` |
| **Status** | Alpha, `0.1.0` — `publish = false` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | all 20 library crates |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

Three reasons, in order of importance:

1. **They are the integration tests.** Depending on all 20 library crates at
   once means a breaking change to any public API fails the workspace build
   immediately, in one place, instead of surfacing as a downstream compile
   error in a user's project.
2. **They are the reference composition.** The crate layering is the hardest
   thing about this workspace to get right from documentation alone. Seven
   working programs that thread `dicom → mask → mesh → model → BC → result`
   are worth more than a diagram.
3. **They make the pipeline demonstrable offline.** Each binary runs against
   the committed synthetic CT in seconds, with no licence, no cloud and no
   scanner. "Does this project actually work?" is answerable in one command.

## What this is

Seven command-line programs, one per project phase, each running a real
pipeline over synthetic data. They are not toy snippets: they are the
integration tests of the workspace. If a phase is broken, the corresponding
binary fails.

Every example runs against the **committed synthetic CT series** in
`test-data/dicom/synthetic_ct/`. No example requires real patient data, and no
example will accept a directory of PHI without you deliberately pointing it
there.

## The binaries

| Binary | Phase | What it does |
|---|---|---|
| `gen-synthetic-ct` | 1 | Writes a synthetic DICOM CT series (femur phantom) to a directory, so the rest of the workspace is testable without a scanner |
| `dicom-to-mesh` | 1 | DICOM CT series → threshold-segmented voxel-hex bone mesh CSV |
| `femur-stress-analysis` | 2 | Segments, meshes, applies an ISO 7206-style 3× body-weight stance load, solves, and reports peak/mean von Mises with a yield screen against cortical bone (~110 MPa) |
| `carotid-wss-screening` | 4 | Solves pulsatile flow in a stenosed carotid tube and reports wall shear stress and OSI |
| `stent-deployment` | 5 | Deploys a Nitinol ring into a compliant vessel; reports radial force, contact pressure, recoil and dogboning |
| `knee-replacement-planning` | 6 | Plane osteotomy on a voxel femur, rigid fragment reposition, and component sizing from landmarks |
| `fda-package` | 7 | Runs a full screening analysis under a 21 CFR Part 11 audit trail and exports the signed package |

The Phase-3 arterial-wall inflation milestone is exercised by the
`tpt-med-tissue` and `tpt-med-viscoelastic` test suites and the golden dataset
`test-data/golden/solid/arterial_wall_inflation.json`, rather than by a
dedicated binary.

## Features

- One binary per project phase, from synthetic data generation through
  regulatory export, so the full arc is walkable in order.
- All seven accept `--help`, read a real input path, and print a small labelled
  table of results.
- All seven exit non-zero on failure. There is no path that prints an error and
  returns success, which is what makes them usable in CI.
- No interactive prompts, no network access, no licence check.
- The `fda-package` binary writes into `test-data/golden/regulatory/`, so the
  regulatory export schema is diffed by the test suite like any other artefact.

## Conventions

- All seven take the input directory as the first positional argument and
  flags after it, following the `getopt`-style long-option form
  (`--output femur.csv --threshold 200`).
- Units follow the workspace convention throughout: **mm**, **N**, **MPa**,
  seconds, and HU for imaging values. Each binary prints units with its
  results rather than assuming the reader knows.
- **Synthetic data only** in every committed code path. The default input is
  `test-data/dicom/synthetic_ct/`.
- Output is a plain table on stdout, so a binary can be piped into any
  downstream tool; the only file outputs are the mesh CSV, the synthetic DICOM
  directory and the signed regulatory package.
- Errors go to stderr with a non-zero exit; results go to stdout.

## Usage

```console
# 1. (Re)generate the synthetic CT series.
cargo run -p tpt-med-examples --bin gen-synthetic-ct -- test-data/dicom/synthetic_ct

# 2. DICOM -> patient-specific bone mesh CSV.
cargo run -p tpt-med-examples --bin dicom-to-mesh -- \
    test-data/dicom/synthetic_ct --output femur_mesh.csv --threshold 200 --smooth 5

# 3. FEMur stress screening under a stance load.
cargo run -p tpt-med-examples --bin femur-stress-analysis -- \
    test-data/dicom/synthetic_ct --body-weight 80

# 4. Hemodynamic screening with WSS and OSI.
cargo run -p tpt-med-examples --bin carotid-wss-screening

# 5. Stent deployment metrics.
cargo run -p tpt-med-examples --bin stent-deployment

# 6. Virtual TKA: osteotomy + reposition + sizing.
cargo run -p tpt-med-examples --bin knee-replacement-planning

# 7. Signed FDA submission package.
cargo run -p tpt-med-examples --bin fda-package
```

Every binary accepts `--help`. All of them run offline, on synthetic data, in
a few seconds.

## API Overview

This package is `publish = false` and exposes **no library API** — only the
binaries listed above. Its manifest exists to depend on every library crate at
once, so:

- a breaking change to any public API breaks the workspace build, and
- the examples double as the reference integration for the crate layering.

## Verification

- Each binary is executed in the workspace test run and in the release
  artifact job, across Linux, macOS and Windows.
- The committed synthetic CT series is regenerated by `gen-synthetic-ct`, and
  the two are asserted to be equivalent by the meshing path — a regeneration
  that changes the data would fail downstream.
- `fda-package` **overwrites**
  `test-data/golden/regulatory/fda_export_example.json` on each run. Its
  timestamps and chain digests differ per run by design, so it is a worked
  artifact rather than a regression anchor, and no test diffs it. A change to
  the export schema is caught by review. The *reference* golden files are what
  `scripts/diff-golden.sh` guards mechanically.
- Binaries exit non-zero on failure; there is no path that prints an error and
  returns success.

## Known Limitations

- **Synthetic data only** in every committed path. Running an example against
  real patient data is the caller's decision and their compliance burden.
- Output paths are command-line arguments with no sandboxing; the binaries will
  overwrite files they are pointed at.
- No interactive prompts or TUI. If a parameter is out of a defensible range,
  the binary says so and exits rather than asking.
- Windows, macOS and Linux are covered by CI; other platforms are untested.

## Related Crates

- [`tpt-med-benches`](https://github.com/tpt-solutions/tpt-medical/tree/master/benches) — the performance counterpart to these correctness binaries.
- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) — `gen-synthetic-ct` and `dicom-to-mesh` start here.
- [`tpt-med-fda`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/regulatory/tpt-med-fda) — `fda-package` exercises the audit trail end to end.
- [`tpt-med-wasm`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-wasm) — the browser path, which these binaries mirror natively.

## Contributing

See the workspace [CONTRIBUTING.md](../../CONTRIBUTING.md). A new example
binary is a welcome addition — keep it dependency-light in *concept* (show one
idea clearly), and **never commit real patient data**.

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.
