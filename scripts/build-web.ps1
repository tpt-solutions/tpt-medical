# Builds the tpt-med-wasm engine glue into web/pkg (shared by the viewer
# and the white-label stent simulator). Run from the repository root.
$ErrorActionPreference = "Stop"

rustup target add wasm32-unknown-unknown | Out-Null
cargo build -p tpt-med-wasm --target wasm32-unknown-unknown --release
wasm-bindgen --out-dir web/pkg --target web `
  target/wasm32-unknown-unknown/release/tpt_med_wasm.wasm
Write-Host "wasm-bindgen glue written to web/pkg/"
