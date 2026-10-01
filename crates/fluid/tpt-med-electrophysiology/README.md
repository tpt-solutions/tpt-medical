# tpt-med-electrophysiology

Monodomain cardiac electrophysiology on voxel geometry: Mitchell-Schaeffer
ionic kinetics, explicit RK2 reaction-diffusion on a `tpt-med-meshing`
segmentation mask, and an S1-S2 single-cell restitution protocol. Zero
external crates.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--electrophysiology-orange)](https://crates.io/crates/tpt-med-electrophysiology)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--electrophysiology-blue)](https://docs.rs/tpt-med-electrophysiology)

| | |
|---|---|
| **Layer** | `fluid` (grouped with the cardiac domain, alongside `tpt-med-cardiovascular`) |
| **Status** | Alpha, `0.1.0` — Stage 1 only |
| **Scope** | RFC 0005 Stage 1 — monodomain on voxel geometry |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.84 |
| **Dependencies** | [`tpt-med-meshing`](../../imaging/tpt-med-meshing) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

`rfcs/0005-cardiac-electrophysiology.md` scopes cardiac EP as a staged
ladder: Stage 0 (0D timing, already in `tpt-med-cardiovascular`), Stage 1
(monodomain tissue conduction, this crate), Stage 2 (ECG forward problem),
Stage 3 (ablation screening). This crate is Stage 1 only, deliberately: it
answers "given a tissue mask and a pacing protocol, what is the local
activation time at each point?" — a screening-fidelity conduction question,
not a diagnostic or treatment-planning one. Splitting it out from
`tpt-med-cardiovascular` rather than extending that crate mirrors why
`tpt-med-nifti` is its own crate from `tpt-med-dicom`: genuinely different
numerics (an explicit PDE time-stepper, not a lumped-parameter ODE), and a
different, much lower validation maturity that should not be hidden inside
an already-shipped crate's version history.

## Features

- **Mitchell-Schaeffer ionic kinetics** (`MitchellSchaefferParams`) — a
  two-variable model (Mitchell & Schaeffer, 2003) chosen over ten-plus-
  variable ionic detail (Luo-Rudy, ten Tusscher) for screening fidelity with
  tractable parameter hygiene. `human_ventricular_default()` is the
  commonly-reproduced generic-ventricular fit from the original paper's
  Table 1 — a screening default, not a validated patient- or species-
  specific one; `new()` validates any custom parameter set (finite,
  positive time constants, `v_gate` in `(0, 1)`).
- **`MonodomainTissue`** — isotropic, homogeneous-diffusivity monodomain
  reaction-diffusion, built directly from a `tpt_med_meshing::SegmentationMask`
  so the same CT/MRI-derived geometry that feeds structural meshing feeds
  electrophysiology. No-flux (Neumann) boundary at the mask's own edge and
  at any tissue/non-tissue voxel interface.
- **Explicit RK2 (Heun) time stepping**, with `max_stable_dt()` computing
  the smaller of two independent bounds: the diffusion CFL condition
  (accounting for anisotropic voxel spacing) and a reaction-stiffness bound
  (`0.2 * tau_in`, the upstroke's own fast timescale) — `step(dt)` errs
  rather than silently diverging if `dt` exceeds either.
- **`activation_map()`** — first `V >= v_gate` crossing time per voxel, the
  output a Stage 2 lead-field projection or a Stage 3 lesion-line evaluation
  would consume.
- **`S1S2Protocol`** — a single-cell (0D, no diffusion) S1-S2 restitution
  protocol: `s1_beats` paced stimuli at `s1_cycle_length_ms`, then one test
  stimulus per `s2_coupling_intervals_ms` entry, each measured independently
  from the same paced state. Used here as this crate's own code-verification
  fixture (see Verification) — not a tissue-level restitution behaviour,
  which additionally depends on conduction and the diffusion coupling.

## Explicit Non-Features

- **Anisotropy is accepted, not sourced.** `set_anisotropy` takes a
  caller-supplied per-voxel fiber field (`FiberConductivity`, axisymmetric
  about the fiber direction) with whole-field validation; the crate still
  has no source for the field itself (atlas/DTI derivation stays external).
  Isotropic by default: diffusivity remains a single
  scalar; real atrial/ventricular conduction is anisotropic along fiber
  direction. Needs a fiber-field source (an atlas or DTI derivation) this
  crate has no source for — see Known Limitations and the RFC's Unresolved
  Questions.
- **No ECG/EGM forward problem (Stage 2) and no ablation screening
  (Stage 3).** Both are named, staged, and deliberately kept at roadmap
  depth in the RFC pending Stage 1 shipping and, for Stage 3, clinical-data
  validation this project does not have.
- **No operator splitting.** Reaction and diffusion are integrated together
  by one unsplit RK2 step. Mitchell-Schaeffer's kinetics are mild enough for
  this at screening fidelity; a stiffer future ionic model might need
  Strang splitting, which is a solver-internal change, not an API change.

## Conventions

- Time is in **ms**, matching the kinetics' own time constants; length is in
  **mm**, matching `SegmentationMask::spacing`; diffusivity is therefore in
  **mm²/ms**.
- Voltage `V` and gate `h` are **dimensionless**, in `[0, 1]` at rest and
  during a normal action potential (not real mV) — this is the model's own
  convention, not a units simplification this crate introduces.
- Grid indexing matches `SegmentationMask` exactly: flat index
  `(z * ny + y) * nx + x`, so a mask and a tissue built from it always agree
  on which voxel is which.
- `stimulate()` sets the active stimulus for exactly the next `step()` call
  and is cleared automatically afterward — call it again every step a pulse
  should remain active, and stop calling it once the pulse ends.

## Usage

### Stimulating a slab and reading the activation map

```rust
use tpt_med_electrophysiology::{MitchellSchaefferParams, MonodomainTissue};
use tpt_med_geometry::Vec3;
use tpt_med_meshing::SegmentationMask;

fn main() -> Result<(), tpt_med_electrophysiology::EpError> {
    let (nx, ny, nz) = (20, 1, 1);
    let mask = SegmentationMask {
        dims: (nx, ny, nz),
        origin: Vec3::ZERO,
        row_dir: Vec3::new(1.0, 0.0, 0.0),
        col_dir: Vec3::new(0.0, 1.0, 0.0),
        slice_dir: Vec3::new(0.0, 0.0, 1.0),
        spacing: (0.25, 1.0, 1.0),
        voxels: vec![true; nx * ny * nz],
        hu: vec![0.0; nx * ny * nz],
    };

    let params = MitchellSchaefferParams::human_ventricular_default();
    let mut tissue = MonodomainTissue::from_mask(&mask, 0.05, params)?;
    let dt = tissue.max_stable_dt() * 0.9;

    // Stimulate one end for a 1ms burst.
    let stim_steps = (1.0 / dt).ceil() as u64;
    for _ in 0..stim_steps {
        tissue.stimulate(|x, _, _| x == 0, 3.0);
        tissue.step(dt)?;
    }
    // Let the wave propagate.
    for _ in 0..(100.0 / dt).ceil() as u64 {
        tissue.step(dt)?;
    }

    let map = tissue.activation_map();
    println!("activation at x=10: {:?} ms", map[10]);
    Ok(())
}
```

### An S1-S2 restitution curve

```rust
use tpt_med_electrophysiology::{MitchellSchaefferParams, S1S2Protocol};

fn main() -> Result<(), tpt_med_electrophysiology::EpError> {
    let params = MitchellSchaefferParams::human_ventricular_default();
    let protocol = S1S2Protocol {
        s1_cycle_length_ms: 300.0,
        s1_beats: 3,
        // Must clear the steady-state APD (~263ms for the default
        // parameters) or every S2 lands in the refractory period.
        s2_coupling_intervals_ms: vec![300.0, 320.0, 350.0, 400.0, 500.0],
        stimulus_amplitude: 1.0,
        stimulus_duration_ms: 1.0,
    };
    for (di, apd) in protocol.run(&params)? {
        println!("DI={di:.1}ms -> APD={apd:.1}ms");
    }
    Ok(())
}
```

## API Overview

| Item | Purpose |
|---|---|
| `MitchellSchaefferParams::{new, human_ventricular_default}` | Validated kinetic parameters, or the cited screening default |
| `MonodomainTissue::from_mask(&mask, diffusivity, params)` | Builds a resting tissue from a `SegmentationMask` |
| `MonodomainTissue::{dims, time, v_at, h_at}` | Grid dimensions, elapsed time, per-voxel state |
| `MonodomainTissue::stimulate(region, amplitude)` | Sets the active stimulus for the next `step()` only |
| `MonodomainTissue::max_stable_dt()` | The larger of the diffusion CFL and reaction-stiffness bounds, whichever is smaller |
| `MonodomainTissue::step(dt)` | One explicit RK2 step; errs over the stability bound |
| `MonodomainTissue::activation_map()` | First `V >= v_gate` crossing time per voxel |
| `S1S2Protocol::run(&params)` | `(diastolic_interval, apd)` pairs from the S1-S2 protocol |
| `EpError` | `EmptyMask`, `NonFiniteParameter`, `UnstableTimeStep`, `InvalidDiffusivity`, `RestitutionFailed` |
| `Result<T>` | Crate result alias |

## Verification

- **Single-cell kinetics**: `S1S2Protocol` is checked to produce a
  restitution curve where APD increases monotonically with the preceding
  diastolic interval — the qualitative shape Mitchell & Schaeffer's own
  paper reports for this parameter set (no closed-form restitution curve
  exists for this model to check against instead). A too-weak stimulus and
  a cycle length inside the refractory period are both asserted to behave
  sensibly (an explicit error for the former, a `Result` rather than a
  panic for the latter).
- **Diffusion stencil**: `max_stable_dt()` is checked to reject a `dt` above
  its own computed bound, and a `dt` just under it is checked to keep an
  unstimulated tissue exactly at rest and a stimulated 1D line's activation
  times strictly increasing with distance from the stimulus site — a
  conduction wave, not simultaneous or backward activation.
- **Stability bound derivation**: the reaction-stiffness term
  (`0.2 * tau_in`) is not a guess kept only because tests pass — it was
  added after an initial diffusion-only bound produced `NaN` state during
  development (an explicit step exceeding the upstroke's own timescale
  overshoots the excitable manifold instead of saturating). See
  `rfcs/0005-cardiac-electrophysiology.md`'s Verification strategy.
- **Golden dataset**: `test-data/golden/electrophysiology/monodomain_restitution.json`
  pins the exact `S1S2Protocol` output for the parameters in Usage above,
  checked by `scripts/diff-golden.sh`. It is a code-verification fixture
  (self-consistent numerical reproduction, no external reference exists),
  not a validation against patient or ex-vivo data.
- **Not verified**: grid-convergence of tissue-level conduction velocity
  (named in the RFC as a Richardson-extrapolation check, not yet reduced to
  a checked-in number), anisotropic conduction (not modelled), and
  everything Stage 2/3 would need.

## Known Limitations

- **Isotropic conductivity only.** A real fidelity limit, not a
  simplification with no cost — see Explicit Non-Features.
- **No tissue-level restitution study.** `S1S2Protocol` is single-cell; a
  tissue's own restitution additionally depends on conduction and wavefront
  curvature, which this protocol does not simulate.
- **`stimulate()`'s one-step contract is easy to get wrong.** Forgetting to
  call it again on a step where a pulse should still be active silently
  ends the pulse early, rather than erroring — there is no way to detect
  "the caller meant to keep stimulating" from inside the tissue.
- **No visualisation or activation-map export helpers.** `activation_map()`
  returns a flat `Vec<Option<f64>>`; turning that into an image or a WebGL
  overlay (the Stage 3 product vision) is out of scope here.

## Related Crates

- [`tpt-med-meshing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-meshing) — `SegmentationMask` is this crate's only input; the same mask that feeds structural meshing feeds electrophysiology.
- [`tpt-med-cardiovascular`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/fluid/tpt-med-cardiovascular) — Stage 0 (0D timing/Windkessel), the boundary-condition scaffolding this crate's future stages will consume.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Stages 2/3
follow the roadmap in `rfcs/0005-cardiac-electrophysiology.md`; anisotropic
conductivity and the promotion path to `tpt-science`'s electrophysiology
crate are tracked in that RFC's Unresolved Questions.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. Not a
diagnostic medical device. This crate's Stage 1 scope is explicitly
screening-fidelity with population-level parameters — see
`rfcs/0005-cardiac-electrophysiology.md`'s ASME V&V 40 discussion before
using it for anything beyond illustration.
