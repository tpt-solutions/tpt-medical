# Changelog

All notable changes to `tpt-med-wasm` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- `WasmMeshPipeline::new` now ingests **Enhanced (multi-frame) CT files**:
  each framed payload is parsed with the dicom crate's frame-list entry
  point, so a multi-frame upload contributes all of its frames to the
  meshed volume instead of being refused. Single-frame uploads are
  unchanged.
- Crate README with a complete JavaScript usage example, the unit table for
  every binding parameter, and the DICOM framing specification.

### Fixed
- Documented that the `wasm-bindgen` CLI version must match the exact pin in
  `Cargo.toml`; a mismatch produces a glue module that fails at load time with
  an opaque error. Any change to the pin must update `scripts/build-web.*` and
  `docs/book/src/wasm.md` in the same PR.

### Notes
- The **DICOM payload framing is public API**: concatenated files, each
  prefixed with a 4-byte little-endian length. Changing it breaks every
  existing JavaScript caller, so it is semver-major.
- The MSRV floor is 1.82, set by `wasm-bindgen 0.2.128` (needs 1.77) and by
  `const` float arithmetic in `tpt-med-units`.

## [0.1.0] - 2026-09-22

### Added
- `WasmMeshPipeline` — constructed from framed DICOM bytes and an HU threshold;
  threshold-segments bone, builds the voxel-hex mesh, and applies the HU-derived
  per-element moduli. Exposes `node_count`, `element_count`, `mean_modulus`
  and `max_modulus` as getters.
- **Flat buffer exports for direct WebGL2 upload** — `node_positions`
  (`Vec<f32>`, mm), `element_nodes` (`Vec<u32>`, 8 indices per hex,
  element-major) and `element_moduli` (`Vec<f32>`, MPa).
- `wasm_solve_stance_load(pipeline, load_newtons)` — linear hex FEM solve
  pinning the bottom node band and sharing the load uniformly over the top
  band. Returns `WasmStressResult` with `max_von_mises`, `mean_von_mises`,
  `max_displacement` and `iterations`.
- `wasm_deploy_stent(expanded, crimped, crowns, stiffness, lumen, pressure,
  compliance)` — Nitinol ring deployment into a pressure–diameter compliant
  vessel, `D(p) = lumen + compliance·p`. Returns `WasmStentResult` with
  `diameter`, `radial_force`, `contact_pressure`, `recoil` and `dogboning`.
- Errors are surfaced as `JsValue`, so failures arrive as JS exceptions and
  are catchable with `try`/`catch`.
- `crate-type = ["cdylib", "rlib"]`, so the bindings are covered by native
  `cargo test` rather than only by a browser.
- Documented build step: `cargo build --target wasm32-unknown-unknown` followed
  by `wasm-bindgen --out-dir web/pkg --target web`, wrapped by
  `scripts/build-web.sh` / `.ps1`.

### Fixed
- **Uniform traction distribution** in `wasm_solve_stance_load` — the stance
  load is shared uniformly over the loaded node band rather than applied to a
  single node.
- **`element_nodes()` connectivity export** — previously the element-to-node
  mapping was not reachable from JavaScript, so a mesh could not be rendered.

### Verification
- Native unit tests run the bindings under `cargo test`, including
  `framed_payload_splits_into_files`, which pins the length-prefix framing.
- CI type-checks the crate for `wasm32-unknown-unknown` on every push, so a
  host-only construct cannot land.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
