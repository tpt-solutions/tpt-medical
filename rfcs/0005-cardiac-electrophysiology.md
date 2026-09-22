# RFC 0005: Cardiac Electrophysiology

- **Status:** Draft
- **Started:** 2026-09-20
- **Crates:** `tpt-med-cardiovascular` (0D today), future `tpt-med-electrophysiology`

## Summary

Scope and staged design for cardiac electrophysiology (EP): from the 0D
timing/waveform models already shipped, through monodomain tissue
simulation, toward patient-specific arrhythmia screening.

## Motivation

The design document lists electrophysiology as the substrate capability
(`tpt-science` electrophys + `tpt-math-signal-filter` for ECG/PPG boundary
conditions). EP is the largest untouched domain in the stack and the one
where the zero-cloud story is most compelling (atrial fibrillation ablation
planning runs today on huge cloud clusters). It is also the furthest from
the current validation base, so the entry point must be deliberately small.

## Design

### Stage 0 — shipped (in `tpt-med-cardiovascular`)

Cardiac timing as data: `FlowWaveform` (harmonic synthesis with mean +
pulsatility), Windkessel 2/3-element models with RK4 integration, FFR
calculator. These provide the boundary conditions hemodynamics needs and the
timing scaffolding EP will consume.

### Stage 1 — monodomain on voxel geometry (new crate)

- Ionic model: **Mitchell–Schaeffer** (two-variable, captures the action
  potential and restitution with 6 parameters) rather than ten Tushode-tier
  ionic detail; screening fidelity, tractable parameter hygiene.
- Tissue: monodomain reaction-diffusion, explicit/explicit- RK2 time
  stepping on the voxel grid (reuses the mask infra from `tpt-med-meshing`
  for atrial geometry from CT/MRI).
- Stimulation protocol: S1–S2 pacing per the restitutions-screening
  convention.
- Verification: plane-wave conduction velocity vs analytic grid corrections
  (conductivity tensor vs CV calibration), restitution curve stability.
- Dependencies: std-only, consistent with the workspace philosophy;
  `tpt-science`'s electrophys crate is the pinned integration point for
  ionic-model breadth once promoted.

### Stage 2 — ECG/EGM forward problem

Pseudo-bidomain lead-field projection to synthesize body-surface ECGs from
tissue state; validated against published 12-lead reference
recordings (QT/conduction intervals as quantitative metrics).

### Stage 3 — ablation screening workflow

Lesion-line placement on patient geometry → re-induction of the arrhythmia
circuit → success metric. This stage is where ASME V&V 40 risk/influence
assessment becomes stringent (High/Direct) and is explicitly gated on
clinical-data validation availability.

## Standards context

- Verification follows the monodomain convergence literature; conduction
  velocity calibration is a code-verification metric with analytic answers.
- Validation realism: EP tissue parameters are population-level unless
  patient-specific data (latency maps) exist; the V&V 40 assessment must
  reflect that honestly (model influence "contributing" until validated
  against patient data).

## Alternatives considered

- **Bidomain immediately**: the extracellular potential matters only for
  the ECG stage; monodomain + lead fields is the standard pragmatic ladder.
- **Coupled excitation–contraction**: valuable for cardiac resynchronization
  questions but doubles validation burden; deferred past Stage 2.

## Unresolved questions

- Anisotropy handling on voxel grids (fiber fields need an atlas or DTMRI
  derivation — geometry pipeline work in `tpt-med-segmentation`).
- Real-time constraint: Stage 3 in-browser (WASM + WebGL visualisation of
  activation maps) is the product vision; feasibility study needed at
  0.25 mm atrial resolution.
