#!/usr/bin/env bash
# Builds the in-browser WASM engine and generates the JS glue.
# wasm-bindgen-cli version must match the wasm-bindgen crate in Cargo.lock.
set -euo pipefail
cd "$(dirname "$0")/.."
ver=$(grep -A1 '^name = "wasm-bindgen"$' Cargo.lock | sed -n 's/^version = "\(.*\)"/\1/p')
rustup target add wasm32-unknown-unknown
cargo build -p tpt-med-wasm --target wasm32-unknown-unknown --release
if ! command -v wasm-bindgen >/dev/null || [ "$(wasm-bindgen --version | awk '{print $2}')" != "$ver" ]; then
  echo "installing wasm-bindgen-cli $ver ..."
  cargo install wasm-bindgen-cli --version "$ver" --force
fi
wasm-bindgen --out-dir web/pkg --target web \
  target/wasm32-unknown-unknown/release/tpt_med_wasm.wasm
echo "web/pkg/ ready"
