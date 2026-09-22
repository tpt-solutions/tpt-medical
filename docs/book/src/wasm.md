# WebAssembly and the browser apps

The entire imaging→mesh→solve stack and the stent deployment model compile
to `wasm32-unknown-unknown`; patient bytes enter the module and only derived
scalars/buffers (meshes, stress statistics, deployment metrics) cross back
to JavaScript.

## Building the engine glue

```console
rustup target add wasm32-unknown-unknown
./scripts/build-web.sh        # macOS / Linux
# or
.\scripts\build-web.ps1       # Windows
```

The script runs the two underlying commands:

```console
cargo build -p tpt-med-wasm --target wasm32-unknown-unknown --release
wasm-bindgen --out-dir web/pkg --target web \
  target/wasm32-unknown-unknown/release/tpt_med_wasm.wasm
```

`web/pkg/` (git-ignored, ~190 KB `.wasm`) is shared by both browser apps
below. Requires [`wasm-bindgen-cli`](https://crates.io/crates/wasm-bindgen-cli)
**matching the pinned version in `crates/core/tpt-med-wasm/Cargo.toml`**
(currently `=0.2.128`):

```console
cargo install wasm-bindgen-cli --version 0.2.128 --locked
```

## Serving

ES modules and `fetch` need HTTP — serve from the repository root:

```console
python -m http.server 8080
```

| URL | App |
|---|---|
| `http://localhost:8080/web/viewer/` | End-to-end planning viewer |
| `http://localhost:8080/web/stent-simulator/` | White-label stent component demo |

## Viewer: CT → mesh → FEM, fully in-browser

`web/viewer/index.html` is a dependency-free WebGL2 viewer (orbit, zoom,
modulus colour ramp). Without the engine it still loads mesh CSVs from
`dicom-to-mesh`. With `web/pkg/` built it adds:

1. **Simulation (WASM)** — select a DICOM series (or fetch the synthetic CT
   demo), threshold-segment bone, voxel-hex mesh, then solve a stance load
   (pinned distal band, uniform proximal traction) and report peak/mean
   von Mises stress, displacement, and CG iterations.
2. **Stent deployment (WASM)** — parameterised ring deployment with
   pressure–diameter vessel compliance; reports equilibrium diameter,
   radial force, contact pressure, recoil, and dogboning.

If `web/pkg/` is missing the page still works from CSV uploads and shows
`Engine: not built`.

## White-label stent component

`web/stent-simulator/` packages the deployment model as the
`<tpt-stent-simulator>` custom element — brand name, logo, accent colour,
theme, and disclaimer are attributes; parameters and results are
attributes, methods, and `tpt-deploy`/`tpt-error` custom events. See
[`web/stent-simulator/README.md`](https://github.com/tpt-solutions/tpt-medical/blob/master/web/stent-simulator/README.md)
for the full API and the npm packaging notes.

## API surface (`tpt-med-wasm`)

| Binding | Purpose |
|---|---|
| `WasmMeshPipeline` | Framed DICOM bytes + HU threshold → mesh buffers and statistics |
| `wasm_solve_stance_load` | Linear hex FEM solve of a pipeline mesh under a stance load |
| `wasm_deploy_stent` | Ring deployment into a compliant vessel (`D(p) = lumen + c·p`) |

DICOM framing: concatenated files, each prefixed with a `u32`
little-endian length (see `frameDicom` in the viewer).
