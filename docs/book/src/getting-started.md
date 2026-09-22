# Getting started

Requirements: Rust 1.82+ (stable). No other toolchain dependencies.

```console
git clone https://github.com/tpt-solutions/tpt-medical
cd tpt-medical
cargo test --workspace          # full verification suite
cargo run -p tpt-med-examples --bin gen-synthetic-ct -- test-data/dicom/synthetic_ct
cargo run -p tpt-med-examples --bin dicom-to-mesh -- \
    test-data/dicom/synthetic_ct --output femur.csv --threshold 200 --smooth 5
```

The committed test data is synthetic only; regenerate it any time. The
repository's CI runs formatting, Clippy (`-D warnings`), the full test
suite on dev and release profiles, native builds for three OSes plus
`wasm32-unknown-unknown`, and `cargo deny` license enforcement.
