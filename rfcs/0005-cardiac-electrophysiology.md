# RFC 0005: Cardiac Electrophysiology

- **Status:** Accepted — Stage 1 implemented (2026-09-27); Stage 2's
  unbounded-bath lead-field slice implemented (2026-10-03, see the Stage 2
  section's status note); the torso-model remainder and Stage 3 stay at
  roadmap depth
- **Started:** 2026-09-20
- **Crates:** `tpt-med-cardiovascular` (0D today), new `tpt-med-electrophysiology`
  (Stage 1)

## Summary

Scope and staged design for cardiac electrophysiology (EP): from the 0D
timing/waveform models already shipped, through a monodomain tissue solver
on voxel geometry (Stage 1, the concrete near-term scope this revision makes
buildable), toward the ECG forward problem and ablation screening (Stages 2–3,
deliberately kept at roadmap depth, gated on validation data this project
does not yet have).

## Motivation

The design document lists electrophysiology as the substrate capability
(`tpt-science` electrophys + `tpt-math-signal-filter` for ECG/PPG boundary
conditions). EP is the largest untouched domain in the stack and the one
where the zero-cloud story is most compelling (atrial fibrillation ablation
planning runs today on huge cloud clusters). It is also the furthest from
the current validation base, so the entry point must be deliberately small.

**Question of interest (Stage 1):** given a tissue mask, a conductivity, and
a pacing protocol, what is the local activation time and action-potential
morphology at each point in the tissue? **Model risk:** low at Stage 1 — no
Stage 1 output feeds a device-sizing or surgical-planning decision; it is a
screening-fidelity conduction model whose own docs must say so. **Model
influence:** "supporting" in ASME V&V 40 terms for Stage 1 (illustrative
conduction behaviour, population-level parameters, no patient-specific
calibration), rising to "contributing" only if a later stage calibrates
against patient-specific latency-map data, and to "controlling" only if
Stage 3's ablation screening is ever used to select a lesion set for a real
procedure — a step this RFC does not authorize and explicitly gates on
validation availability (see Stage 3 below).

## Detailed design

### Stage 0 — shipped (in `tpt-med-cardiovascular`)

Cardiac timing as data: `FlowWaveform` (harmonic synthesis with mean +
pulsatility), Windkessel 2/3-element models with RK4 integration, FFR
calculator. These provide the boundary conditions hemodynamics needs and the
timing scaffolding EP will consume. Unchanged by this revision.

### Stage 1 — monodomain on voxel geometry (new crate `tpt-med-electrophysiology`)

#### Ionic model

**Mitchell–Schaeffer** (Mitchell & Schaeffer, *Bull. Math. Biol.* 65:767–793,
2003) — a two-variable model chosen over ten-plus-variable ionic detail
(Luo–Rudy, ten Tusscher) for screening fidelity with tractable parameter
hygiene: six parameters, each independently meaningful (an upstroke time
constant, a repolarization time constant, gate-opening/closing time
constants, and a gating threshold), rather than a dozen-plus conductances
requiring a full voltage-clamp dataset to fit responsibly.

Per-cell (0D) kinetics, in the model's own dimensionless voltage `V ∈ [0, 1]`
and gate `h ∈ [0, 1]`:

```
dV/dt = J_in(V, h) + J_out(V) + J_stim
dh/dt = (1 - h) / τ_open   if V <  v_gate
      = -h / τ_close       if V >= v_gate

J_in(V, h)  =  h · V² · (1 - V) / τ_in
J_out(V)    = -V / τ_out
```

Parameters: `τ_in`, `τ_out`, `τ_open`, `τ_close`, `v_gate`. A sixth
"parameter" some presentations count separately, the stimulus current
`J_stim`, is a per-call input, not a tissue property. Default parameters are
the original paper's own fitted values for a generic ventricular myocyte —
cited in code alongside the type, and documented as **screening defaults**,
not patient- or even validated-species-specific, exactly the posture RFC
0004 already takes for its Nitinol parameter set.

#### Tissue: monodomain reaction-diffusion on the voxel grid

```
∂V/∂t = D · ∇²V + J_in(V, h) + J_out(V) + J_stim
```

- **Isotropic, homogeneous `D` for v0.** Anisotropic (fiber-direction)
  conductivity needs a fiber-direction field this project has no source for
  yet (an atlas or DTI derivation) — carried over from the previous revision
  as an explicit Unresolved Question, not solved by this one.
- **Discretisation:** 7-point (3D, face-neighbour) finite-difference
  Laplacian on the same voxel grid `tpt-med-meshing::SegmentationMask`
  already produces from CT/MRI, with no-flux (Neumann) boundary at the
  mask's own boundary — a voxel outside the mask contributes no diffusive
  flux, so tissue current cannot leak through a boundary that has no
  physical tissue beyond it.
- **Time stepping:** explicit RK2 (Heun's method), applied to the whole
  coupled `(V, h)` system per voxel — no operator splitting for v0 (Strang
  splitting is the standard choice for stiffer ionic models; Mitchell–
  Schaeffer's own kinetics are mild enough that unsplit RK2 is adequate at
  screening fidelity, and splitting is easy to add later without an API
  change since it is purely a solver-internal detail).
- **Stability bound:** explicit diffusion on a 7-point 3D stencil requires
  `Δt ≤ Δx² / (6D)`; `step` computes this from the tissue's own `Δx`/`D` and
  returns a typed error rather than silently stepping into an unstable,
  diverging state if the caller requests a larger `Δt`.

#### Public API (sketch)

```rust
/// Mitchell-Schaeffer kinetic parameters. `human_ventricular_default()`
/// cites Mitchell & Schaeffer (2003) directly in its doc comment — a
/// screening default, not a validated patient- or species-specific fit.
pub struct MitchellSchaefferParams {
    pub tau_in: f64,
    pub tau_out: f64,
    pub tau_open: f64,
    pub tau_close: f64,
    pub v_gate: f64,
}

impl MitchellSchaefferParams {
    pub fn human_ventricular_default() -> Self;
}

/// A monodomain tissue on a voxel grid, built directly from a
/// `tpt-med-meshing` segmentation mask so the same CT/MRI-derived geometry
/// that feeds structural meshing feeds electrophysiology.
pub struct MonodomainTissue { /* dims, voxel size, mask, D, V[], h[], params */ }

impl MonodomainTissue {
    pub fn from_mask(
        mask: &tpt_med_meshing::SegmentationMask,
        diffusivity: f64,
        params: MitchellSchaefferParams,
    ) -> Result<Self>;

    /// Injects `J_stim` at voxels where `region` returns true, for use as an
    /// S1 or S2 stimulus.
    pub fn stimulate(&mut self, region: impl Fn(usize, usize, usize) -> bool, amplitude: f64);

    /// One RK2 step. Errs (`EpError::UnstableTimeStep`) rather than
    /// stepping if `dt` violates the diffusion stability bound.
    pub fn step(&mut self, dt: f64) -> Result<()>;

    /// First V-upstroke-crossing time per voxel, `None` if never activated —
    /// the activation map a Stage-2 lead-field projection or a Stage-3
    /// lesion-line evaluation would consume.
    pub fn activation_map(&self) -> Vec<Option<f64>>;
}

/// S1-S2 restitution protocol on a single 0D cell (no diffusion) — the
/// direct analogue of how Mitchell & Schaeffer's own paper characterises
/// the model, used here as this crate's code-verification fixture rather
/// than a tissue-level feature.
pub struct S1S2Protocol {
    pub s1_cycle_length_ms: f64,
    pub s1_beats: u32,
    pub s2_coupling_intervals_ms: Vec<f64>,
}

impl S1S2Protocol {
    /// Returns `(diastolic_interval, action_potential_duration)` pairs, one
    /// per `s2_coupling_intervals_ms` entry — a restitution curve.
    pub fn run(&self, params: &MitchellSchaefferParams) -> Vec<(f64, f64)>;
}
```

#### Error handling

New `EpError` (mirroring `DicomError`'s shape): `UnstableTimeStep { dt,
max_stable_dt }`, `EmptyMask` (a `SegmentationMask` with no tissue voxels —
nothing to build a tissue on), `NonFiniteParameter { name }` (a
`MitchellSchaefferParams` field that is NaN/infinite — same finiteness
discipline `QctCalibration`/RFC 0007 already apply to their own numeric
inputs).

#### `#![forbid(unsafe_code)]`

Applies, matching every other crate in the workspace; the explicit stencil
and RK2 update are plain array iteration, no unsafe surface needed.

#### Dependencies

std-only for Stage 1, consistent with the workspace philosophy of shipping a
tractable in-house core before an optional substrate upgrade path (the same
posture `tpt-med-hemodynamics`/`tpt-fem` already established).
`tpt-science`'s electrophysiology crate remains the pinned integration point
for ionic-model breadth (a full ten Tusscher-tier model, say) once promoted,
behind a cargo feature — not a Stage 1 dependency.

### Stage 2 — ECG/EGM forward problem

Pseudo-bidomain lead-field projection to synthesize body-surface ECGs from
tissue state (an `activation_map` plus the underlying `V(t)` per voxel);
validated against published 12-lead reference recordings (QT/conduction
intervals as quantitative metrics). Kept at roadmap depth: the lead-field
formulation, torso-conductivity assumptions, and electrode-placement
convention are real design decisions this revision does not make. A
follow-up RFC picks this up once Stage 1 has shipped and been exercised.

  **First slice implemented (2026-10-03, maintainer direction): the
  unbounded-bath source integral** — `LeadFieldProjection` in
  `tpt-med-electrophysiology` evaluates the Geselowitz/pseudo-bidomain
  integral face-wise over the tissue with no torso model. The closure
  question this section flagged as a real design decision is settled for
  this slice by the integration-by-parts identity: no bath faces (a
  fictitious `u = 0` jump at the tissue boundary double-counts the surface
  term the interior faces already carry), and the tests pin the identity's
  observable consequences — closed fronts cancel, boundary-open fronts
  carry the solid-angle signal with a dipolar far field, and an open
  uniform front converges to an independent solid-angle quadrature. What
  remains at roadmap depth is exactly the torso half: geometry,
  conductivity heterogeneity, and electrode transfer impedances, which a
  quantitative 12-lead comparison needs and a follow-up RFC must own.

### Stage 3 — ablation screening workflow

Lesion-line placement on patient geometry → re-induction of the arrhythmia
circuit → success metric. This stage is where ASME V&V 40 risk/influence
assessment becomes stringent (High/Direct) and is explicitly gated on
clinical-data validation availability — this revision does not authorize
Stage 3 for anything beyond research/screening use, and a future RFC
proposing Stage 3 must itself state what validation data closes that gate,
not merely that Stage 1/2 shipped successfully.

### Alternatives considered

- **Bidomain immediately.** The extracellular potential matters only for the
  ECG stage; monodomain + lead fields is the standard pragmatic ladder, and
  bidomain's added cost (a second PDE, a harder linear solve per step) buys
  nothing Stage 1's question of interest needs.
- **A richer ionic model at Stage 1** (Luo–Rudy, ten Tusscher). Rejected for
  the same reason RFC 0004 chose a simplified superelastic law over full
  Souza–Auricchio at Level 1: a dozen-plus conductances need a voltage-clamp
  dataset to fit responsibly, and Stage 1's question of interest (conduction
  timing and morphology at screening fidelity) does not need that detail.
  `tpt-science`'s electrophysiology crate is the named upgrade path once a
  real use case needs it.
- **Coupled excitation–contraction.** Valuable for cardiac resynchronization
  questions but doubles validation burden; deferred past Stage 2, same as
  the previous revision of this RFC already concluded.
- **Operator-split (Strang) time integration at Stage 1.** The more common
  choice in production monodomain codes, and cheaper per step for stiffer
  ionic models. Deferred, not rejected: Mitchell–Schaeffer's mild kinetics
  make unsplit RK2 adequate at screening fidelity, and switching later is a
  solver-internal change, not an API change — revisit only if profiling or a
  stiffer future ionic model actually needs it.

### Drawbacks

- Isotropic conductivity is a real fidelity limit, not a simplification with
  no cost: real atrial/ventricular conduction is anisotropic along fiber
  direction, and Stage 1's conduction-velocity outputs will read wrong by a
  direction-dependent factor for any tissue where that matters. Named,
  tracked as the top Unresolved Question, not hidden behind "screening
  fidelity" language.
- No operator splitting means Stage 1's time step is bounded by the
  diffusion CFL condition regardless of the ionic kinetics' own stiffness
  margin — for a fine voxel grid (sub-millimeter atrial wall thickness) this
  can force a smaller `Δt`, and thus more steps, than a split scheme would
  need. Accepted for v0's simplicity; revisit if a real workflow's runtime
  is dominated by this.
- The `S1S2Protocol` restitution curve is a single-cell (0D) characterisation
  used for verification; it is not itself a tissue-level restitution
  behaviour (which additionally depends on conduction, memory effects at the
  wavefront, and the diffusion coupling) — a caller must not treat it as
  standing in for a tissue-level restitution study.

## Verification strategy

- **Code verification:**
  - Single-cell kinetics: the `S1S2Protocol` restitution curve (APD vs. DI)
    is compared, in shape (monotonically increasing APD with DI, saturating
    at long DI), against Mitchell & Schaeffer's own published restitution
    curve for their fitted parameter set — a qualitative-shape check backed
    by a specific cited figure, not an invented closed-form target (the
    model has no simple closed-form restitution curve to check against
    instead).
  - Diffusion stencil: the 7-point Laplacian is checked against a
    manufactured solution (a sinusoidal `V` field with a known analytic
    `∇²V`, no ionic kinetics) on a small grid — a direct discretisation-error
    check independent of the ionic model.
  - Stability bound: a deliberately-oversized `Δt` (above the computed
    `Δx²/(6D)` bound) must be rejected by `step`, and a test confirms that a
    `Δt` just under the bound remains numerically stable (bounded `V`) over
    many steps while one just over it is rejected outright rather than
    silently diverging.
- **Calculation verification:** planar-wave conduction velocity, measured
  from `activation_map` along a 1D-like corridor of a 3D grid, is checked
  for grid convergence via Richardson extrapolation (halving `Δx` and
  confirming the measured CV converges toward a fine-grid extrapolated
  value at the expected first-order rate for this discretisation) — not
  against a closed-form analytic CV, which the Mitchell–Schaeffer kinetics
  do not admit in general; the extrapolated fine-grid value is the reference,
  named as such rather than mislabelled "analytic."
- **Validation:** none at Stage 1 — no patient-specific or ex-vivo
  conduction-velocity dataset is used, consistent with the "supporting"
  V&V 40 influence level claimed above. Validation against real tissue data
  is explicitly deferred to whichever later stage first calibrates against
  patient-specific measurements.
- **Golden datasets:** a new `test-data/golden/electrophysiology/` entry
  (planar-wave conduction velocity + single-cell restitution curve on a
  small synthetic slab) is added once Stage 1 ships, exercised by
  `scripts/diff-golden.sh` the same as every other domain's golden dataset.
- **What remains unverified:** anisotropic conduction (not modelled, so
  nothing to verify yet), tissue-level restitution (only single-cell is
  characterised), and everything Stage 2/3 would need (lead-field accuracy,
  re-induction success against real ablation outcomes) — named here rather
  than left implicit, matching this RFC's own V&V 40 "supporting, not
  controlling" claim for Stage 1.

## Unresolved questions

- **Anisotropy handling on voxel grids.** Fiber fields need an atlas or
  DTI derivation — geometry-pipeline work, not solver work, and a
  prerequisite for Stage 1's isotropic limitation to be lifted. Settled by:
  whoever picks up a fiber-field source, informed by which imaging modality
  (DTI-MRI vs. a generic rule-based atlas) a real dataset actually has
  available.
- **Real-time constraint.** Stage 3 in-browser (WASM + WebGL visualisation
  of activation maps) is the product vision; a feasibility study at
  clinically-relevant atrial resolution (sub-millimeter) is needed before
  committing to it as a Stage 3 requirement rather than an aspiration.
  Settled by: whoever runs that feasibility study, once Stage 1's actual
  per-step cost is measured (this revision does not benchmark Stage 1
  itself, since it has not shipped).
- **Promotion path to `tpt-science`'s electrophysiology crate.** Named as
  the upgrade path for ionic-model breadth, but the concrete trigger
  (a specific ionic model a real workflow needs that Mitchell–Schaeffer
  cannot represent) is not yet identified. Settled by: whoever hits that
  workflow.
