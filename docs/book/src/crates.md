# Crates

See the [README crate table](https://github.com/tpt-solutions/tpt-medical/blob/master/README.md#crates)
for the full list with status. Layering:

- **core**: `tpt-med-units`, `tpt-med-geometry`, `tpt-med-core` (domain
  types, audit traits), `tpt-med-wasm` (browser bindings).
- **imaging**: `tpt-med-dicom` (parsing + HU mapping), `tpt-med-meshing`
  (segmentation, voxel-to-hex meshing, CSV export).
- **solid**: `tpt-med-biomechanics` (hex FEM), `tpt-med-tissue`,
  `tpt-med-bone`, `tpt-med-viscoelastic`, `tpt-med-cartilage`.
- **fluid**: `tpt-med-hemodynamics` (voxel CFD, WSS/OSI),
  `tpt-med-cardiovascular` (Windkessel, FFR).
- **devices**: `tpt-med-stents`, `tpt-med-orthopedics`, `tpt-med-wear`.
- **surgical**: `tpt-med-surgical-planning`, `tpt-med-implant-sizing`.
- **regulatory**: `tpt-med-audit` (SHA-256/HMAC), `tpt-med-fda`,
  `tpt-med-vv40`.
- **applications**: `tpt-med-examples` (binaries, `publish = false`),
  `tpt-med-benches` (benchmark suite, `publish = false`).

Every crate is `#![forbid(unsafe_code)]` and `std`-only; the TPT substrate
crates (`tpt-fem`, `tpt-science`, `tpt-engineering`, `tpt-math`) are pinned
in the root manifest as the sanctioned nonlinear/high-fidelity integration
points.

## Per-crate documentation

Each crate directory carries two files, and they are the per-crate source of
truth. Crates version and release independently, so neither the root
`CHANGELOG.md` nor this guide describes a single crate's behaviour in enough
detail to use it.

### `README.md`

Required sections, in this order:

1. **Title, one-line description, and badges** (crates.io, docs.rs).
2. **Metadata table** — layer, status, license, MSRV, dependencies, and a link
   to the crate's own `CHANGELOG.md`.
3. **Why** — the problem the crate solves, and the design decision that
   resolves it. State the alternatives and why they were rejected.
4. **Features** — a bulleted list of what actually ships.
5. **Conventions** — units, sign conventions, indexing, coordinate systems, and
   anything where a caller could reasonably be wrong.
6. **Usage** — runnable examples using the real API.
7. **API Overview** — a table of every public item and its purpose.
8. **Verification** — what is verified, against what reference, and what that
   check *proves*. Not "we have tests".
9. **Known Limitations** — what this crate cannot do. Required; a README that
   only lists strengths is not documentation.
10. **Related Crates**, **Contributing**, **License**, and the regulatory
    disclaimer.

### `CHANGELOG.md`

[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format, with:

- A header stating that the crate versions independently, and linking the
  workspace changelog.
- `[Unreleased]` — planned work, and an explicit note of which changes are
  semver-breaking (this is where an "adding a variant to this enum is breaking"
  warning belongs, before someone discovers it).
- A version section per release, with `### Added` / `### Changed` /
  `### Fixed` / `### Removed` subsections, plus `### Verification` and
  `### Known limitations` where relevant.
- Compare and release links at the bottom.

Both files are enforced by the `crate-docs` CI job, which also checks each
manifest for crates.io-compliant `keywords` and `categories`. Categories must
be slugs from the published registry — note that **`medical-science` is not a
valid crates.io category**, which is why the crates here use `science` plus a
domain keyword instead.

