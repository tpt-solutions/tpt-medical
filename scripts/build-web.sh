#!/usr/bin/env sh
# Builds the tpt-med-wasm engine glue into web/pkg (shared by the viewer
# and the white-label stent simulator). Run from the repository root.
set -eu

rustup target add wasm32-unknown-unknown
cargo build -p tpt-med-wasm --target wasm32-unknown-unknown --release
wasm-bindgen --out-dir web/pkg --target web \
  target/wasm32-unknown-unknown/release/tpt_med_wasm.wasm
echo "wasm-bindgen glue written to web/pkg/"
