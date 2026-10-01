# tpt-med-wasm

WebAssembly bindings for the `tpt-medical` zero-cloud simulation stack — the
imaging → mesh → solve pipeline and the Nitinol stent deployment model,
exposed to JavaScript.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--wasm-orange)](https://crates.io/crates/tpt-med-wasm)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--wasm-blue)](https://docs.rs/tpt-med-wasm)

| | |
|---|---|
| **Layer** | `core` (bindings) |
| **Target** | `wasm32-unknown-unknown` (also builds natively for tests) |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.84 (`wasm-bindgen 0.2.128` needs 1.77) |
| **Dependencies** | `wasm-bindgen` (exact-pinned), plus `tpt-med-dicom`, `-meshing`, `-biomechanics`, `-stents`, `-geometry`, `-units` |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

The privacy claim of this project — *patient DICOM data never leaves the
hospital* — is only credible if the whole imaging-to-simulation pipeline can
run **inside the browser sandbox**. That requires shipping the parser, the
mesher and the FEM core to `wasm32-unknown-unknown`.

This crate is the thin, deliberately small binding layer that makes those
crates callable from JavaScript. It adds no algorithms of its own: it owns
serialization (flat typed arrays, so the JS side never parses a CSV) and error
translation (Rust errors become JS exceptions).

**Data flow:** a patient's DICOM bytes enter the module; only derived results
leave it — mesh buffers, stress scalars, deployment metrics. No pixel data, no
identifiers, no network calls.

## Features

- `WasmMeshPipeline` — framed DICOM bytes + HU threshold → threshold-segmented
  voxel-hex mesh with HU-derived per-element moduli.
- `wasm_solve_stance_load` — linear FEM solve of a pipeline mesh under a
  physiological stance load.
- `wasm_deploy_stent` — Nitinol ring deployment into a pressure–diameter
  compliant vessel, reporting the ASTM F2394-style metric family.
- Flat `Float32Array`/`Uint32Array` exports (`node_positions`,
  `element_nodes`, `element_moduli`) for direct WebGL2 upload.
- Errors surface as `JsValue` so `try`/`catch` works naturally in JS.
- `crate-type = ["cdylib", "rlib"]` — the same code runs under `cargo test`
  natively, so the bindings are covered by CI, not just by a browser.

## Conventions

| Parameter | Unit |
|---|---|
| DICOM payload | bytes (framed, see below) |
| Diameters (`expanded`, `crimped`, `lumen`) | mm |
| `crown_stiffness` | N/mm per crown |
| `pressure_mpa` | MPa |
| `vessel_compliance_mm_per_mpa` | mm/MPa — slope of `D(p) = lumen + c·p` (0 ⇒ rigid lumen) |
| `load_newtons` | N |
| Stresses / moduli | MPa |

**DICOM framing.** The constructor takes a concatenation of DICOM files, each
prefixed with a **4-byte little-endian length**. JavaScript builds this with a
`DataView` — see `frameDicom` in `web/viewer/`. This avoids a length-delimited
container format in the hot path and keeps the payload a single `Uint8Array`.



## Building the WebAssembly

```console
rustup target add wasm32-unknown-unknown
./scripts/build-web.sh          # macOS / Linux
# or
.\scripts\build-web.ps1         # Windows
```

The script runs the two underlying commands:

```console
cargo build -p tpt-med-wasm --target wasm32-unknown-unknown --release
wasm-bindgen --out-dir web/pkg --target web \
  target/wasm32-unknown-unknown/release/tpt_med_wasm.wasm
```

The `wasm-bindgen` **CLI version must match the exact pin in `Cargo.toml`**
(currently `=0.2.128`):

```console
cargo install wasm-bindgen-cli --version 0.2.128 --locked
```

The generated `web/pkg/` directory (git-ignored, ~190 KB `.wasm`) is shared by
both browser applications: the end-to-end planning viewer (`web/viewer/`) and
the white-label `<tpt-stent-simulator>` component (`web/stent-simulator/`).

## Usage

### JavaScript

```javascript
import init, {
  WasmMeshPipeline,
  wasm_solve_stance_load,
  wasm_deploy_stent,
} from "./pkg/tpt_med_wasm.js";

// Concatenate DICOM files, each prefixed with a u32 little-endian length.
function frameDicom(files) {
  const total = files.reduce((n, f) => n + 4 + f.byteLength, 0);
  const out = new Uint8Array(total);
  const dv = new DataView(out.buffer);
  let off = 0;
  for (const f of files) {
    dv.setUint32(off, f.byteLength, true);
    out.set(f, off + 4);
    off += 4 + f.byteLength;
  }
  return out;
}

await init();

try {
  // 1. DICOM bytes + HU threshold -> meshed pipeline.
  const pipeline = new WasmMeshPipeline(frameDicom(ctFiles), 200.0);
  console.log(pipeline.element_count, pipeline.mean_modulus);

  // 2. Direct WebGL2 upload buffers.
  gl.bindBuffer(gl.ARRAY_BUFFER, ...pipeline.node_positions());   // mm
  gl.bindBuffer(gl.ARRAY_BUFFER, ...pipeline.element_nodes());     // connectivity
  gl.bindBuffer(gl.ARRAY_BUFFER, ...pipeline.element_moduli());    // MPa

  // 3. Physiological stance load solve.
  const stress = wasm_solve_stance_load(pipeline, 2400.0); // N
  console.log(stress.max_von_mises, stress.mean_von_mises,
              stress.max_displacement, stress.iterations);

  // 4. Stent deployment into a compliant vessel.
  const stent = wasm_deploy_stent(
    4.0,    // expanded diameter (mm)
    1.2,    // crimped diameter (mm)
    8,      // crowns
    0.9,    // crown stiffness (N/mm)
    3.6,    // lumen diameter (mm)
    0.013,  // intraluminal pressure (MPa)
    0.35,   // vessel compliance (mm/MPa)
  );
  console.log(stent.diameter, stent.radial_force,
              stent.contact_pressure, stent.recoil, stent.dogboning);
} catch (e) {
  console.error("pipeline failed", e); // DICOM/solve errors arrive as JS values
}
```

## API Overview

| Binding | Purpose |
|---|---|
| `WasmMeshPipeline::new(payload, threshold_hu)` | Framed DICOM bytes + HU threshold → meshed pipeline (throws `JsValue` on failure) |
| `WasmMeshPipeline::node_count` | Number of shared corner nodes |
| `WasmMeshPipeline::element_count` | Number of hexahedral elements |
| `WasmMeshPipeline::mean_modulus`, `::max_modulus` | Element modulus statistics (MPa) |
| `WasmMeshPipeline::node_positions -> Vec<f32>` | Flat `[x, y, z, ...]` buffer (mm) |
| `WasmMeshPipeline::element_nodes -> Vec<u32>` | Flat connectivity, 8 node indices per hex, element-major |
| `WasmMeshPipeline::element_moduli -> Vec<f32>` | Per-element Young's modulus (MPa) for colour mapping |
| `wasm_solve_stance_load(pipeline, load_newtons)` | Linear hex FEM solve; bottom band fixed, top band loaded uniformly |
| `WasmStressResult::max_von_mises`, `::mean_von_mises`, `::max_displacement`, `::iterations` | Solve results (MPa, mm, iteration count) |
| `wasm_deploy_stent(expanded, crimped, crowns, stiffness, lumen, pressure, compliance)` | Ring deployment into a compliant vessel |
| `WasmStentResult::diameter`, `::radial_force`, `::contact_pressure`, `::recoil`, `::dogboning` | ASTM F2394-style deployment metrics |

## Verification

- Native unit tests run the bindings under `cargo test` (the crate is an `rlib`
  as well as a `cdylib`), including `framed_payload_splits_into_files`, which
  pins the length-prefix framing.
- CI type-checks the crate for `wasm32-unknown-unknown` on every push, so a
  host-only construct cannot land.
- The uniform-traction distribution in `wasm_solve_stance_load` and the
  `element_nodes()` connectivity export are both covered by regression tests —
  both were historical bugs, which is why they are pinned explicitly.

## Known Limitations

- **Bindings only.** This crate adds no algorithms of its own. Everything it
  exposes lives in `tpt-med-dicom`, `tpt-med-meshing`,
  `tpt-med-biomechanics` and `tpt-med-stents`, so its API is only as good as
  those crates' and it inherits their limitations.
- **Three bindings, not a general API surface.** There is no binding for
  `tpt-med-tissue`, `tpt-med-hemodynamics`, `tpt-med-cardiovascular`,
  `tpt-med-orthopedics`, `tpt-med-wear`, `tpt-med-surgical-planning`,
  `tpt-med-implant-sizing`, `tpt-med-audit` or `tpt-med-fda`. A browser
  application that wants a CFD or regulatory result cannot get it from here
  today, and adding one is a deliberate scope decision, not an oversight.
- **WASM is the only real target.** A `cdylib` for a native host is built for
  tests, not as a supported FFI surface; do not depend on the symbol layout.
- **Synchronous and single-threaded.** Every binding blocks. A large solve
  will freeze the browser main thread, and the caller must move it to a Web
  Worker — which means the objects are not `Send` and must be constructed
  inside the worker.
- **Flat arrays only.** Nodes, elements and moduli are returned as
  `Float32Array`/`Uint32Array`, which loses `f64` precision on positions and
  moduli. Fine for rendering; not sufficient to reconstruct an exact model
  from the browser.
- **No incremental or partial results.** `wasm_solve_stance_load` returns only
  after the full solve, so a viewer cannot stream progress or show a
  convergence trace.
- **`~190 KB` `.wasm` is a hard budget.** Adding a binding that pulls in a
  heavier dependency will need a fresh size review, and the medical crates are
  deliberately dependency-light to keep this number small.
- **`wasm-bindgen` version lockstep.** The CLI must match the exact pin; this
  is documented but is still a foot-gun for new contributors.

## Related Crates

- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) — parses the incoming payload.
- [`tpt-med-meshing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-meshing) — builds the hex mesh.
- [`tpt-med-biomechanics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-biomechanics) — the FEM solver behind `wasm_solve_stance_load`.
- [`tpt-med-stents`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-stents) — the deployment model behind `wasm_deploy_stent`.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). If you change the
`wasm-bindgen` pin, update `scripts/build-web.*` and
`docs/book/src/wasm.md` in the same PR — the CLI version is a hard
requirement, not a suggestion.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.
