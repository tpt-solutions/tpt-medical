# Introduction

`tpt-medical` is a fully open-source, pure-Rust computational biomechanics
and medical device simulation engine. It compiles to WebAssembly, so the
full imaging-to-simulation pipeline runs **in the browser or on a local
workstation** — patient DICOM data never needs to touch a cloud server.

The stack covers the complete device-analysis lifecycle:

```text
DICOM ─▶ segmentation ─▶ voxel-hex mesh ─▶ FEM / CFD ─▶ metrics ─▶ signed, audited export
```

## Design pillars

1. **Zero-cloud privacy** — WASM-safe, `std`-only crates; no network code in
   the simulation path.
2. **Permissive-only chain** — `MIT OR Apache-2.0` end to end, enforced by
   `cargo deny`; no GPL/academic-license contamination.
3. **Verifiable mathematics** — every solver ships with analytical
   verification tests (ASME V&V 40 code verification) and curated golden
   reference datasets.
4. **Regulatory readiness** — 21 CFR Part 11 audit trails with hash-chain
   integrity and HMAC signatures; ASME V&V 40 credibility matrices.

## Status

All crates are 🚧 **alpha**: APIs are settling, numerical fidelity is
screening-grade (see each RFC for the fidelity ladder), and the software is
**not** cleared or approved for clinical use.
