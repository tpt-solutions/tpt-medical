# Web applications

Browser apps over the `tpt-med-wasm` engine — zero-cloud: uploaded patient
data is processed inside the tab and never transmitted.

| App | Purpose |
|---|---|
| [`viewer/`](viewer/index.html) | WebGL2 surgical planning viewer for mesh CSVs (orbit, zoom, modulus colour ramp). |
| [`stent-sim/`](stent-sim/index.html) | **White-label stent deployment simulator** — live radial-force/COF/contact-pressure/recoil metrics over the WASM engine; rebrand via the `WHITE_LABEL` config block at the top of `index.html`. |
| [`demo/`](demo/index.html) | End-to-end demo: CT upload → WASM segment/mesh → in-browser FEM stance solve → measured lumen → stent deployment, rendered wireframe. |

## Build the engine (once)

```console
rustup target add wasm32-unknown-unknown
cargo build -p tpt-med-wasm --target wasm32-unknown-unknown --release
wasm-bindgen --out-dir web/pkg --target web \
    target/wasm32-unknown-unknown/release/tpt_med_wasm.wasm
```

(`web/pkg/` is generated and git-ignored; wasm-bindgen-cli version must match
the `wasm-bindgen` crate in `Cargo.lock` exactly.)

## Run

Serve the repository root over HTTP (ES modules do not load from `file://`):

```console
python -m http.server 8000
# viewer:    http://localhost:8000/web/viewer/
# stent-sim: http://localhost:8000/web/stent-sim/
# demo:      http://localhost:8000/web/demo/
```

The demo's "load the committed synthetic series" button fetches
`test-data/dicom/synthetic_ct/` (synthetic data only — never load real
patient scans on a shared machine). Both apps bust the engine HTTP cache
per load so rebuilt glue is picked up without a hard refresh.
