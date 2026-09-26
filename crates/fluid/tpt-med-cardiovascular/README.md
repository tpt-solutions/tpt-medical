# tpt-med-cardiovascular

Lumped cardiovascular models — 2- and 3-element Windkessel boundary
conditions, fractional flow reserve (FFR), and cardiac flow waveforms.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--cardiovascular-orange)](https://crates.io/crates/tpt-med-cardiovascular)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--cardiovascular-blue)](https://docs.rs/tpt-med-cardiovascular)

| | |
|---|---|
| **Layer** | `fluid` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | none (leaf within `fluid`) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

A CFD solve of the coronary or carotid vasculature needs *boundary
conditions*, and a physiological one is not "pressure = 0 at the outlet". The
standard choice is a Windkessel model: a lumped two- or three-element
representation of the entire cardiovascular tree distal to the outlet,
integrated in time alongside the CFD solve.

This crate is that lumped model, plus the two scalar clinical outputs that
depend on it — **fractional flow reserve** and **cardiac flow waveforms** —
and nothing else. It is small on purpose: it has no dependencies, it is
closed-form where theory allows, and it can be stepped at the CFD time step
without a performance argument.

## Features

- **3-element Windkessel (Westerhof)** — characteristic impedance `Rc`,
  peripheral resistance `Rp`, compliance `C`, outflow pressure `p_out`,
  obeying `C dp/dt = (1 + Rc/Rp)·Q − (p − p_out)/Rp`.
- **2-element model** — `two_element(r_p, c, p_out)` is `Rc = 0`, not a
  separate code path, so the two models cannot drift apart.
- **Two integrators** — `step_rk4` (explicit RK4, the reference) and
  `simulate` (the semi-implicit scheme used online by the CFD boundary
  coupling). RK4 is the reference because the semi-implicit scheme is
  unconditionally stable but not the most accurate.
- **Analytic quantities** — `steady_state_pressure(flow)` and `time_constant()`,
  so a boundary condition can be sanity-checked without integrating anything.
- **`FractionalFlowReserve`** — `calculate(p_distal, p_aortic)` and
  `is_ischemic(ffr)` against the ≤ 0.80 threshold.
- **`FlowWaveform`** — `flow(t)` for the carotid and coronary beds
  (`carotid_default`, `coronary_default`), which drive pulsatile CFD runs.

## Conventions

- Resistances in **MPa·s/mm³**, compliance in **mm³/MPa**, pressures in
  **MPa**, flow in **mm³/s**, time in **seconds**.

## Usage

### Windkessel boundary condition for a CFD solve

```rust
use tpt_med_cardiovascular::{FlowWaveform, WindkesselModel};

fn main() {
    // 3-element (Westerhof) model.
    let wk = WindkesselModel {
        r_c: 0.05,     // characteristic impedance, MPa*s/mm^3
        r_p: 1.20,     // peripheral resistance
        c: 0.002,      // compliance, mm^3/MPa
        p_out: 0.010,  // venous pressure, MPa (~10 mmHg)
    };

    // Analytic checks, before integrating anything.
    let q_steady = 1.0;                       // mm^3/s
    let p_expected = wk.p_out + q_steady * (wk.r_c + wk.r_p);
    assert!((wk.steady_state_pressure(q_steady) - p_expected).abs() < 1e-12);
    assert!(wk.time_constant() > 0.0);

    // RK4 is the reference integrator.
    let p = wk.step_rk4(p_expected, q_steady, 1.0e-4);
    assert!((p - p_expected).abs() < 1e-3);

    // Full cycle against a coronary waveform.
    let wave = FlowWaveform::coronary_default();
    let dt = 1.0e-4;
    let trace = wk.simulate(p_expected, dt, 10_000, |t| wave.flow(t));
    let peak = trace.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!(peak > p_expected); // systolic overshoot
    assert!(trace.iter().all(|p| p.is_finite()));
}
```

### Fractional flow reserve

```rust
use tpt_med_cardiovascular::FractionalFlowReserve;

fn main() {
    // 80% stenosis: 100 -> 20 mmHg across the lesion.
    let ffr = FractionalFlowReserve::calculate(0.0266, 0.1333); // 20 / 75 mmHg
    assert!((ffr - 0.8).abs() < 0.01);
    assert!(!FractionalFlowReserve::is_ischemic(ffr));

    // Severe lesion.
    let severe = FractionalFlowReserve::calculate(0.0400, 0.1333); // 30 / 100 mmHg
    assert!(FractionalFlowReserve::is_ischemic(severe));
}
```

## API Overview

| Item | Purpose |
|---|---|
| `WindkesselModel { r_c, r_p, c, p_out }` | 3-element (Westerhof) parameters |
| `WindkesselModel::two_element(r_p, c, p_out)` | 2-element model (`Rc = 0`) |
| `::dp_dt(p, flow) -> f64` | `dp/dt` at a state point (MPa/s) |
| `::steady_state_pressure(flow) -> f64` | `p_out + Q(Rc + Rp)` — the analytic fixed point |
| `::time_constant() -> f64` | Characteristic time constant (s) |
| `::step_rk4(p, flow, dt) -> f64` | Explicit RK4 step — the reference integrator |
| `::simulate(p0, dt, steps, flow_fn) -> Vec<f64>` | Semi-implicit integration; `flow_fn` supplies `Q(t)` |
| `FractionalFlowReserve::calculate(p_distal, p_aortic) -> f64` | The pressure ratio |
| `FractionalFlowReserve::is_ischemic(ffr) -> bool` | `ffr <= 0.80` |
| `FlowWaveform::flow(t) -> f64` | Instantaneous flow (mm³/s) at time `t` |
| `FlowWaveform::carotid_default()` | Carotid-bed waveform |
| `FlowWaveform::coronary_default()` | Coronary-bed waveform |

## Verification

Lumped circulation models have closed-form solutions, and all of them are
asserted:

- **Fixed point:** `dp_dt` is zero at `steady_state_pressure(flow)`, for a
  sweep of `flow` and across 2- and 3-element configurations.
- **Time constant:** the analytic `time_constant()` is compared against the
  time for the explicit RK4 solution to decay from an initial perturbation
  back to within 1/e of equilibrium — the integrator and the closed form must
  agree, not merely be self-consistent.
- **Two-element limit:** setting `Rc = 0` reproduces `two_element(...)`
  exactly, element by element, which is what keeps the two models from
  drifting.
- **FFR boundary:** the ischaemia classifier is tested *at* 0.80 (asserted
  ischaemic — the threshold is inclusive), just above, and just below.
  An off-by-one here is a clinical misclassification.
- **FFR degeneracies:** zero aortic pressure does not produce `NaN` or a
  panic, and FFR > 1 (possible with measurement noise) is reported rather
  than silently clamped.
- **Waveform sanity:** both waveforms are strictly positive, bounded, and
  return to a diastolic baseline within a cycle; `flow(0)` equals
  `flow(cycle)` so the waveform is continuous.
- Golden reference dataset: `test-data/golden/fluid/coronary_ffr.json`.

## Known Limitations

- **Lumped, not distributed.** A Windkessel boundary condition reproduces
  the *global* impedance of the distal tree; it does not reproduce wave
  reflection shape or the timing of the reflected wave, so it is the wrong
  tool for pulse-wave-velocity or augmentation-index work.
- Two- and three-element only. The four-element (Windkessel– Westerhof with
  characteristic impedance split) and the non-linear
  pressure–flow relationships used in systemic circulation modelling are not
  implemented.
- FFR here is the *pressure-ratio* definition. The waveform-based
  instantaneous-hyperbolic FFR used in some catheter workflows is not
  implemented.
- Waveforms are fixed analytic shapes, not measured from a specific patient.
  Per-patient waveform fitting is a caller-side concern.

## Related Crates

- [`tpt-med-hemodynamics`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/fluid/tpt-med-hemodynamics) — the CFD solver these boundary conditions drive.
- [`tpt-med-stents`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/devices/tpt-med-stents) — stent deployment, the main post-CFD use case.
- [`tpt-med-units`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-units) — the `Pressure`/`FlowRate` conventions these units follow.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). Cite the source
of every physiological constant and the reference for any waveform shape. New
lumped-parameter models require an [RFC](../../../rfcs).

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use. **FFR is a
clinical measurement**; the value computed here is a research estimate and
must not be used for diagnosis or treatment decisions.

- `dp_dt` is in MPa/s.
- FFR is the **dimensionless pressure ratio** `p_distal / p_aortic`, with
  distal pressure taken *after* the stenosis and aortic pressure proximal to
  it.
- Ischaemia threshold: **FFR ≤ 0.80** is the accepted clinical cut-off for
  physiologically significant stenosis.
- Waveforms are per cardiac cycle; `flow(t)` expects `t` in the same units as
  the cycle length used to build the waveform.

