//! Lumped cardiovascular models: Windkessel boundary conditions, FFR, and
//! cardiac flow waveforms.
//!
//! The 3-element Windkessel (Westerhof) — characteristic impedance `Rc`,
//! peripheral resistance `Rp`, arterial compliance `C` — obeys
//!
//! ```text
//! C dp/dt = (1 + Rc/Rp)·Q − (p − p_out)/Rp
//! ```
//!
//! with `Rc = 0` reducing to the 2-element model. Integrators: explicit
//! RK4 (reference) and the semi-implicit scheme used online by the CFD
//! boundary coupling.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// A Windkessel model.
#[derive(Debug, Clone, Copy)]
pub struct WindkesselModel {
    /// Characteristic (proximal) resistance, MPa·s/mm³.
    pub r_c: f64,
    /// Peripheral (distal) resistance, MPa·s/mm³.
    pub r_p: f64,
    /// Arterial compliance, mm³/MPa.
    pub c: f64,
    /// Outflow (venous) pressure, MPa.
    pub p_out: f64,
}

impl WindkesselModel {
    /// Two-element Windkessel (`Rc = 0`).
    pub fn two_element(r_p: f64, c: f64, p_out: f64) -> Self {
        Self {
            r_c: 0.0,
            r_p,
            c,
            p_out,
        }
    }

    /// dp/dt at a state point (MPa/s).
    pub fn dp_dt(&self, p: f64, flow: f64) -> f64 {
        ((1.0 + self.r_c / self.r_p) * flow - (p - self.p_out) / self.r_p) / self.c
    }

    /// Steady-state pressure for constant flow (MPa).
    pub fn steady_state_pressure(&self, flow: f64) -> f64 {
        self.p_out + flow * (self.r_c + self.r_p)
    }

    /// Diastolic decay time constant `τ = Rp·C` (s).
    pub fn time_constant(&self) -> f64 {
        self.r_p * self.c
    }

    /// Advances one step with RK4 (flow held constant over the step).
    pub fn step_rk4(&self, p: f64, flow: f64, dt: f64) -> f64 {
        let k1 = self.dp_dt(p, flow);
        let k2 = self.dp_dt(p + 0.5 * dt * k1, flow);
        let k3 = self.dp_dt(p + 0.5 * dt * k2, flow);
        let k4 = self.dp_dt(p + dt * k3, flow);
        p + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    }

    /// Simulates `steps` of duration `dt` with a per-step flow callback,
    /// returning the pressure history.
    pub fn simulate(&self, p0: f64, dt: f64, steps: usize, flow: impl Fn(f64) -> f64) -> Vec<f64> {
        let mut p = p0;
        let mut history = Vec::with_capacity(steps);
        for s in 0..steps {
            p = self.step_rk4(p, flow(s as f64 * dt), dt);
            history.push(p);
        }
        history
    }
}

/// Four-element Windkessel: the 3-element model plus an **inertance** `L`
/// (blood and wall inertia) in the series branch — the classic four-element
/// form used for systemic circulation modelling, e.g. Stergiopoulos,
/// Young & Westerhof (1999). Unlike the 3-element model, it is genuinely
/// second-order, so it is **pressure-driven**: a prescribed inlet (aortic)
/// pressure produces the flow, as state `(p, Q)`
///
/// ```text
/// C dp/dt = Q − (p − p_out)/Rp          (parallel Rp‖C bank)
/// L dQ/dt = p_in − p − Rc·Q             (series Rc + L branch)
/// ```
///
/// `L` must be positive: at `L = 0` the flow becomes algebraic and this
/// formulation is singular (use [`WindkesselModel`] for the flow-driven
/// series branch). Steady state is exact and independent of `L`:
/// `Q = (p_in − p_out)/(Rc + Rp)`, `p = p_in − Rc·Q`.
#[derive(Debug, Clone, Copy)]
pub struct FourElementWindkessel {
    /// Characteristic (proximal) resistance, MPa·s/mm³.
    pub r_c: f64,
    /// Inertance of the series branch, MPa·s²/mm³.
    pub l: f64,
    /// Peripheral (distal) resistance, MPa·s/mm³.
    pub r_p: f64,
    /// Arterial compliance, mm³/MPa.
    pub c: f64,
    /// Outflow (venous) pressure, MPa.
    pub p_out: f64,
}

impl FourElementWindkessel {
    /// The four-element model as a 3-element [`WindkesselModel`] plus an
    /// inertance.
    pub fn new(base: WindkesselModel, inertance: f64) -> Self {
        Self {
            r_c: base.r_c,
            l: inertance,
            r_p: base.r_p,
            c: base.c,
            p_out: base.p_out,
        }
    }

    /// State derivatives `(dp/dt, dQ/dt)` at `(p, Q)` for inlet pressure
    /// `p_in`.
    pub fn derivs(&self, p: f64, q: f64, p_in: f64) -> (f64, f64) {
        let dp = (q - (p - self.p_out) / self.r_p) / self.c;
        let dq = (p_in - p - self.r_c * q) / self.l;
        (dp, dq)
    }

    /// Advances one RK4 step with `p_in` held constant over the step,
    /// returning the new `(p, Q)`.
    pub fn step_rk4(&self, p: f64, q: f64, p_in: f64, dt: f64) -> (f64, f64) {
        let (k1p, k1q) = self.derivs(p, q, p_in);
        let (k2p, k2q) = self.derivs(p + 0.5 * dt * k1p, q + 0.5 * dt * k1q, p_in);
        let (k3p, k3q) = self.derivs(p + 0.5 * dt * k2p, q + 0.5 * dt * k2q, p_in);
        let (k4p, k4q) = self.derivs(p + dt * k3p, q + dt * k3q, p_in);
        (
            p + dt / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p),
            q + dt / 6.0 * (k1q + 2.0 * k2q + 2.0 * k3q + k4q),
        )
    }

    /// Steady state `(p, Q)` for a constant inlet pressure `p_in`.
    /// Independent of `L` and `C` — inertia and compliance only shape the
    /// transient.
    pub fn steady_state(&self, p_in: f64) -> (f64, f64) {
        let q = (p_in - self.p_out) / (self.r_c + self.r_p);
        (p_in - self.r_c * q, q)
    }

    /// Simulates `steps` of duration `dt` with a per-step inlet-pressure
    /// callback, returning the `(p, Q)` state history.
    pub fn simulate(
        &self,
        p0: f64,
        q0: f64,
        dt: f64,
        steps: usize,
        p_in: impl Fn(f64) -> f64,
    ) -> Vec<(f64, f64)> {
        let mut state = (p0, q0);
        let mut history = Vec::with_capacity(steps);
        for s in 0..steps {
            state = self.step_rk4(state.0, state.1, p_in(s as f64 * dt), dt);
            history.push(state);
        }
        history
    }
}

/// Vascular waterfall (Starling-resistor) pressure–flow relation: a vessel
/// that collapses once its transmural pressure falls below a critical
/// closing pressure, after which flow is **independent of downstream
/// pressure** (Permutt & Bromberger-Barnea; the "vascular waterfall" used in
/// systemic and cerebral circulation modelling).
///
/// ```text
/// Q = max(0, (p_upstream − p_collapse) / R)
/// ```
#[derive(Debug, Clone, Copy)]
pub struct WaterfallResistor {
    /// Resistance above the collapse threshold, MPa·s/mm³.
    pub r: f64,
    /// Critical closing (collapse) pressure, MPa.
    pub p_collapse: f64,
}

impl WaterfallResistor {
    /// Flow through the collapsible segment (mm³/s). `p_downstream` is
    /// accepted for call-site clarity and deliberately unused: downstream
    /// independence is the defining property of the waterfall.
    pub fn flow(&self, p_upstream: f64, p_downstream: f64) -> f64 {
        let _ = p_downstream;
        if p_upstream <= self.p_collapse {
            0.0
        } else {
            (p_upstream - self.p_collapse) / self.r
        }
    }
}

/// Fractional flow reserve.
#[derive(Debug, Clone, Copy)]
pub struct FractionalFlowReserve;

impl FractionalFlowReserve {
    /// `FFR = Pd / Pa` measured under maximal hyperemia.
    pub fn calculate(p_distal: f64, p_aortic: f64) -> f64 {
        p_distal / p_aortic
    }

    /// Clinical classification (0.80 cutoff, standard CARES-like
    /// convention): ≥ 0.80 non-ischemic, < 0.80 ischemia-inducing.
    pub fn is_ischemic(ffr: f64) -> bool {
        ffr < 0.80
    }
}

/// The instantaneous wave-free ratio (iFR): the distal/proximal pressure
/// ratio restricted to the **wave-free period** of diastole, when
/// microvascular resistance is naturally low and stable — a resting
/// (non-hyperemic) alternative to FFR built on pressure *waveforms*
/// (Davies et al. 2012; the commercial convention averages over a window
/// beginning 5 ms after end-systole and ending 5 ms before the next
/// systole, detected from the dP/dt minima).
///
/// The waveform form here takes sampled aortic and distal pressures plus
/// the window as cycle fractions, so the caller (or a detector upstream)
/// owns the window detection convention.
#[derive(Debug, Clone, Copy)]
pub struct InstantaneousWaveFreeRatio {
    /// Wave-free window start as a fraction of the cardiac cycle
    /// (clinical convention ≈ 0.45).
    pub window_start: f64,
    /// Wave-free window end as a fraction of the cycle (≈ 0.95).
    pub window_end: f64,
    /// Clinical threshold: iFR < 0.90 flags ischemia.
    pub ischemic_threshold: f64,
}

impl Default for InstantaneousWaveFreeRatio {
    fn default() -> Self {
        Self {
            window_start: 0.45,
            window_end: 0.95,
            ischemic_threshold: 0.90,
        }
    }
}

impl InstantaneousWaveFreeRatio {
    /// Mean `Pd/Pa` over the wave-free window, from evenly sampled
    /// waveforms (exactly one cycle, `pa[i]`/`pd[i]` at `t = i/n · T`).
    /// Samples with non-finite or non-positive aortic pressure are skipped
    /// rather than poisoning the mean; an empty window returns `NaN`.
    pub fn calculate(&self, pa: &[f64], pd: &[f64]) -> f64 {
        assert_eq!(pa.len(), pd.len(), "pressure waveforms must pair");
        let n = pa.len();
        let mut sum = 0.0;
        let mut count = 0usize;
        for i in 0..n {
            let t = (n as f64).recip() * i as f64;
            if t < self.window_start || t >= self.window_end {
                continue;
            }
            if pa[i].is_finite() && pd[i].is_finite() && pa[i] > 0.0 {
                sum += pd[i] / pa[i];
                count += 1;
            }
        }
        if count == 0 {
            f64::NAN
        } else {
            sum / count as f64
        }
    }

    /// The ischemia classification at the configured threshold.
    pub fn is_ischemic(&self, ifr: f64) -> bool {
        ifr < self.ischemic_threshold
    }
}

/// Pulsatile flow waveform `Q(t)` over one cardiac cycle (mm³/s).
#[derive(Debug, Clone)]
pub struct FlowWaveform {
    /// Cycle length (s).
    pub cycle: f64,
    /// Mean flow over the cycle (mm³/s).
    pub mean: f64,
    /// Pulsatile (zero-mean) component amplitude fraction of mean.
    pub pulsatility: f64,
    /// Number of harmonics used to synthesize the pulsatile component.
    pub harmonics: usize,
}

impl FlowWaveform {
    /// Evaluates `Q(t)` (periodic).
    pub fn flow(&self, t: f64) -> f64 {
        if self.mean == 0.0 {
            return 0.0;
        }
        let omega = 2.0 * core::f64::consts::PI / self.cycle;
        let mut q = self.mean;
        for n in 1..=self.harmonics {
            q += self.mean * self.pulsatility * (n as f64).recip() * (omega * (n as f64) * t).sin();
        }
        q
    }

    /// A canonical carotid-like waveform: mean 6 mL/s, pulsatility 0.6,
    /// 1 Hz, 6 harmonics.
    pub fn carotid_default() -> Self {
        Self {
            cycle: 1.0,
            mean: 6000.0 / 60.0, // 6 mL/s = 100 mm³/s
            pulsatility: 0.6,
            harmonics: 6,
        }
    }

    /// A canonical coronary-like waveform: mean 60 mL/min/1.4e6 mm³...
    /// expressed as mean 1.0 mm³/ms = 1000 mm³/s scale-free testing value.
    pub fn coronary_default() -> Self {
        Self {
            cycle: 0.8,
            mean: 80.0,
            pulsatility: 0.5,
            harmonics: 4,
        }
    }
}

/// A flow waveform **fitted to measurement**: mean plus per-harmonic
/// amplitude and phase from a discrete Fourier transform of one cycle of
/// evenly sampled data — the patient-specific counterpart to the fixed
/// analytic [`FlowWaveform`] shapes, which remain the screening default.
///
/// ```text
/// Q(t) = mean + Σₙ Aₙ·cos(2π n t/T − φₙ)
/// ```
#[derive(Debug, Clone)]
pub struct MeasuredFlowWaveform {
    /// Cycle length (s).
    pub cycle: f64,
    /// Cycle-mean flow (mm³/s).
    pub mean: f64,
    /// `(amplitude, phase)` per harmonic `n = 1..=k`, phase in radians.
    pub harmonics: Vec<(f64, f64)>,
}

impl MeasuredFlowWaveform {
    /// Least-squares fit of the truncated Fourier series above to one
    /// cycle of evenly sampled flow data. The first `keep_harmonics`
    /// harmonics are retained; sample counts beyond `2·k` make the
    /// truncation meaningful rather than exact interpolation.
    pub fn fit(cycle: f64, samples: &[f64], keep_harmonics: usize) -> Self {
        let n = samples.len();
        assert!(n > 0, "cannot fit an empty waveform");
        assert!(keep_harmonics > 0, "at least one harmonic is required");
        let mean = samples.iter().sum::<f64>() / n as f64;
        let mut harmonics = Vec::with_capacity(keep_harmonics);
        for h in 1..=keep_harmonics {
            let mut a = 0.0f64;
            let mut b = 0.0f64;
            for (j, &s) in samples.iter().enumerate() {
                let theta = 2.0 * core::f64::consts::PI * h as f64 * j as f64 / n as f64;
                a += s * theta.cos();
                b += s * theta.sin();
            }
            a *= 2.0 / n as f64;
            b *= 2.0 / n as f64;
            harmonics.push((a.hypot(b), b.atan2(a)));
        }
        Self {
            cycle,
            mean,
            harmonics,
        }
    }

    /// Evaluates the fitted waveform (periodic).
    pub fn flow(&self, t: f64) -> f64 {
        let omega = 2.0 * core::f64::consts::PI / self.cycle;
        let mut q = self.mean;
        for (n, &(amplitude, phase)) in self.harmonics.iter().enumerate() {
            q += amplitude * (omega * (n + 1) as f64 * t - phase).cos();
        }
        q
    }

    /// RMS residual of the fit against the samples it was (or would be)
    /// fitted from — the caller's truncation-error measure.
    pub fn fit_rms(&self, samples: &[f64]) -> f64 {
        let n = samples.len();
        let se: f64 = samples
            .iter()
            .enumerate()
            .map(|(j, &s)| {
                let t = self.cycle * j as f64 / n as f64;
                let e = s - self.flow(t);
                e * e
            })
            .sum();
        (se / n as f64).sqrt()
    }
}

/// A stateful Windkessel boundary condition for lockstep coupling with a
/// `tpt-med-hemodynamics` time step: the CFD outlet `flow` advances the
/// 0-D model, and the returned pressure feeds back as the outlet boundary
/// value for the next step.
///
/// The coupling is explicit (pressure lags flow by one step), which is the
/// standard stable choice when the 0-D time constant is large relative to
/// the CFD step; `WindkesselModel::time_constant` is the number to check.
#[derive(Debug, Clone)]
pub struct CoupledWindkessel {
    model: WindkesselModel,
    pressure: f64,
    steps: u64,
}

impl CoupledWindkessel {
    /// Initialises the boundary at a steady pressure for `mean_flow`.
    pub fn new(model: WindkesselModel, mean_flow: f64) -> Self {
        Self {
            model,
            pressure: model.steady_state_pressure(mean_flow),
            steps: 0,
        }
    }

    /// Advances one CFD step of duration `dt` with the instantaneous outlet
    /// flow, returning the new outlet pressure (MPa) to prescribe.
    pub fn advance(&mut self, flow: f64, dt: f64) -> f64 {
        self.pressure = self.model.step_rk4(self.pressure, flow, dt);
        self.steps += 1;
        self.pressure
    }

    /// Current outlet pressure (MPa).
    pub fn pressure(&self) -> f64 {
        self.pressure
    }

    /// Steps taken since initialisation.
    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// The wrapped model (for time-constant checks at the call site).
    pub fn model(&self) -> &WindkesselModel {
        &self.model
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wk() -> WindkesselModel {
        WindkesselModel {
            r_c: 0.0002,
            r_p: 0.001,
            c: 1000.0,
            p_out: 0.001,
        }
    }

    #[test]
    fn coupled_boundary_reaches_the_same_steady_state() {
        let model = WindkesselModel {
            r_c: 0.0002,
            r_p: 0.001,
            c: 1000.0,
            p_out: 0.001,
        };
        let q = 100.0;
        let mut bc = CoupledWindkessel::new(model, q);
        let dt = 0.01;
        let mut last = bc.pressure();
        for _ in 0..2000 {
            last = bc.advance(q, dt);
        }
        assert!(
            (last - model.steady_state_pressure(q)).abs() < 1e-4,
            "settled {last} vs steady {}",
            model.steady_state_pressure(q)
        );
        assert_eq!(bc.steps(), 2000);
        // Withheld flow: pressure decays toward p_out monotonically.
        let mut prev = f64::INFINITY;
        for _ in 0..500 {
            let p = bc.advance(0.0, dt);
            assert!(p < prev, "not decaying: {p} after {prev}");
            prev = p;
        }
    }

    #[test]
    fn steady_state_and_time_constant() {
        let w = wk();
        let q = 100.0;
        assert!((w.time_constant() - 1.0).abs() < 1e-12);
        let p_inf = w.steady_state_pressure(q);
        let history = w.simulate(0.0, 0.01, 2000, |_| q);
        let final_p = history[history.len() - 1];
        assert!(
            (final_p - p_inf).abs() < 1e-4,
            "final {final_p} vs steady {p_inf}"
        );
    }

    #[test]
    fn diastolic_decay_is_exponential() {
        // With zero flow the pressure decays as p − p_out ∝ e^{−t/τ}.
        let w = WindkesselModel::two_element(0.001, 1000.0, 0.0);
        let p0 = 1.0;
        let history = w.simulate(p0, 0.005, 400, |_| 0.0);
        for (i, &p) in history.iter().enumerate() {
            let t = (i + 1) as f64 * 0.005;
            let expected = p0 * (-t / w.time_constant()).exp();
            assert!((p - expected).abs() < 1e-4, "t={t}: {p} vs {expected}");
        }
    }

    #[test]
    fn two_element_is_three_element_without_rc() {
        let w3 = wk();
        let w2 = WindkesselModel::two_element(w3.r_p, w3.c, w3.p_out);
        // Steady state differs by the Rc drop.
        let q = 50.0;
        assert!(
            (w3.steady_state_pressure(q) - (w2.steady_state_pressure(q) + q * w3.r_c)).abs()
                < 1e-12
        );
    }

    #[test]
    fn ffr_classification() {
        assert!((FractionalFlowReserve::calculate(0.78, 1.0) - 0.78).abs() < 1e-12);
        assert!(FractionalFlowReserve::is_ischemic(0.75));
        assert!(!FractionalFlowReserve::is_ischemic(0.85));
    }

    #[test]
    fn waveform_is_periodic_and_positive() {
        let wf = FlowWaveform::carotid_default();
        assert!((wf.flow(0.0) - wf.flow(wf.cycle)).abs() < 1e-6);
        for i in 0..100 {
            let t = wf.cycle * i as f64 / 100.0;
            assert!(wf.flow(t) > 0.0, "Q(t={t}) = {}", wf.flow(t));
        }
    }

    fn wk4() -> FourElementWindkessel {
        FourElementWindkessel::new(
            WindkesselModel {
                r_c: 0.0002,
                r_p: 0.001,
                c: 1000.0,
                p_out: 0.001,
            },
            1.0e-5, // inertance: underdamped at ω ≈ 10 rad/s for these Rc/Rp/C
        )
    }

    #[test]
    fn four_element_settles_on_the_exact_steady_state() {
        let w = wk4();
        let p_in = 0.013; // ~13 kPa
        let (p_ss, q_ss) = w.steady_state(p_in);
        assert!((q_ss - (p_in - w.p_out) / (w.r_c + w.r_p)).abs() < 1e-15);
        assert!((p_ss - (p_in - w.r_c * q_ss)).abs() < 1e-15);
        let (p, q) = w
            .simulate(0.0, 0.0, 0.005, 4000, |_| p_in)
            .pop()
            .expect("non-empty");
        assert!(
            (p - p_ss).abs() < 1e-6 && (q - q_ss).abs() < 1e-4,
            "settled ({p}, {q}) vs steady ({p_ss}, {q_ss})"
        );
        // Steady state is independent of the inertance.
        let no_inertia = FourElementWindkessel::new(
            WindkesselModel {
                r_c: 0.0002,
                r_p: 0.001,
                c: 1000.0,
                p_out: 0.001,
            },
            1.0e-3,
        );
        let (p2, q2) = no_inertia.steady_state(p_in);
        assert!((p2 - p_ss).abs() < 1e-15 && (q2 - q_ss).abs() < 1e-15);
    }

    #[test]
    fn four_element_inertia_rings_down() {
        // Underdamped: with Rc = 0 and L < 4·Rp²·C the flow overshoots the
        // steady state and oscillates while decaying — the physics the
        // inertance element exists for. The 3-element (first-order) model
        // cannot do this.
        let w = FourElementWindkessel::new(
            WindkesselModel {
                r_c: 0.0,
                r_p: 0.001,
                c: 1000.0,
                p_out: 0.0,
            },
            1.0e-5,
        );
        let p_in = 0.01;
        let q_ss = (p_in - w.p_out) / w.r_p;
        let history = w.simulate(0.0, 0.0, 0.005, 300, |_| p_in);
        let excess: Vec<f64> = history.iter().map(|&(_, q)| q - q_ss).collect();
        // A sign change of q − q_ss is the overshoot; a second one is the ring.
        let sign_changes = excess.windows(2).filter(|w| w[0] * w[1] < 0.0).count();
        assert!(sign_changes >= 2, "no ring-down: {sign_changes} crossings");
        // And the oscillation decays at the theoretical envelope rate:
        // for Rc = 0 the envelope is e^{−t/(2·Rp·C)}, τ_env = 2 s, so peaks
        // one ring apart (~1.26 s here) must fall to ≈ e^{−0.63} ≈ 0.53.
        let peak_early = excess.iter().take(40).fold(0.0, |a, &x| x.abs().max(a));
        let peak_late = excess
            .iter()
            .skip(260)
            .take(40)
            .fold(0.0, |a, &x| x.abs().max(a));
        let ratio = peak_late / peak_early;
        assert!(
            (ratio - 0.53).abs() < 0.1,
            "envelope ratio {ratio} (peaks {peak_early} -> {peak_late})"
        );
    }

    #[test]
    fn four_element_dc_gain_matches_the_resistance_divider() {
        // For a linear system, the time-average of Q over full periods of a
        // periodic p_in equals the response to the mean p_in — the DC gain
        // (p̄_in − p_out)/(Rc + Rp). This checks the whole transient against
        // superposition, not just the endpoint.
        let w = wk4();
        let p_mean = 0.013;
        let amp = 0.002;
        let period = 1.0;
        let dt = 0.005;
        let steps_per_cycle = (period / dt) as usize;
        let settle = 20 * steps_per_cycle;
        let average_over = 20 * steps_per_cycle;
        let p_in = |t: f64| p_mean + amp * (2.0 * core::f64::consts::PI * t / period).sin();
        let history = w.simulate(0.0, 0.0, dt, settle + average_over, p_in);
        let mean_q = history[settle..].iter().map(|&(_, q)| q).sum::<f64>() / average_over as f64;
        let expected = (p_mean - w.p_out) / (w.r_c + w.r_p);
        assert!(
            (mean_q - expected).abs() < 1e-5 * expected.abs(),
            "mean flow {mean_q} vs DC gain {expected}"
        );
    }

    #[test]
    fn ifr_averages_only_over_the_wave_free_window() {
        // Distal pressure tracks 90% of aortic inside the window but only
        // 70% outside: the iFR must read 0.90, not the whole-cycle mean.
        let ifr = InstantaneousWaveFreeRatio::default();
        let n = 100;
        let mut pa = Vec::with_capacity(n);
        let mut pd = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / n as f64;
            let pulse = 1.0 + 0.2 * (2.0 * core::f64::consts::PI * t).sin();
            pa.push(10.0 * pulse);
            pd.push(if (0.45..0.95).contains(&t) {
                9.0 * pulse
            } else {
                7.0 * pulse
            });
        }
        let value = ifr.calculate(&pa, &pd);
        assert!((value - 0.90).abs() < 1e-9, "iFR {value}");
        // The whole-cycle mean, for contrast, is dragged down by the
        // extrasystolic samples.
        let whole: f64 = (0..n).map(|i| pd[i] / pa[i]).sum::<f64>() / n as f64;
        assert!(
            whole < value - 0.05,
            "windowing must matter: {whole} vs {value}"
        );
        // Classification at the 0.90 clinical threshold (strict).
        assert!(!ifr.is_ischemic(0.90));
        assert!(ifr.is_ischemic(0.89));
        // Degenerate window: NaN, not a division by zero.
        let empty = InstantaneousWaveFreeRatio {
            window_start: 0.9,
            window_end: 0.9,
            ..ifr
        };
        assert!(empty.calculate(&pa, &pd).is_nan());
    }

    #[test]
    fn dft_fit_recovers_an_exact_fourier_series() {
        // Sample the analytic carotid waveform (a pure 6-harmonic sine
        // series) and fit with 6 harmonics: the fit is exact and the
        // reconstructed waveform agrees pointwise.
        let wf = FlowWaveform::carotid_default();
        let n = 200;
        let samples: Vec<f64> = (0..n)
            .map(|i| wf.flow(wf.cycle * i as f64 / n as f64))
            .collect();
        let fitted = MeasuredFlowWaveform::fit(wf.cycle, &samples, 6);
        assert!((fitted.mean - wf.mean).abs() < 1e-9);
        assert!(
            fitted.fit_rms(&samples) < 1e-9,
            "rms {}",
            fitted.fit_rms(&samples)
        );
        for i in [0, 17, 55, 130, 199] {
            let t = wf.cycle * i as f64 / n as f64;
            assert!(
                (fitted.flow(t) - wf.flow(t)).abs() < 1e-6,
                "t={t}: {} vs {}",
                fitted.flow(t),
                wf.flow(t)
            );
        }
        // Periodicity of the fitted form.
        assert!((fitted.flow(0.0) - fitted.flow(fitted.cycle)).abs() < 1e-9);
    }

    #[test]
    fn dft_truncation_error_shrinks_with_harmonics() {
        let wf = FlowWaveform::carotid_default();
        let n = 128;
        let samples: Vec<f64> = (0..n)
            .map(|i| wf.flow(wf.cycle * i as f64 / n as f64))
            .collect();
        let rms2 = MeasuredFlowWaveform::fit(wf.cycle, &samples, 2).fit_rms(&samples);
        let rms4 = MeasuredFlowWaveform::fit(wf.cycle, &samples, 4).fit_rms(&samples);
        let rms6 = MeasuredFlowWaveform::fit(wf.cycle, &samples, 6).fit_rms(&samples);
        assert!(rms2 > rms4 && rms4 > rms6, "{rms2} {rms4} {rms6}");
        assert!(rms6 < 1e-9);
    }

    #[test]
    fn waterfall_flow_is_independent_of_downstream_pressure() {
        let wf = WaterfallResistor {
            r: 0.001,
            p_collapse: 0.002,
        };
        // Above the collapse threshold: linear in upstream pressure…
        assert!((wf.flow(0.012, 0.001) - 10.0).abs() < 1e-12);
        assert!((wf.flow(0.013, 0.001) - 11.0).abs() < 1e-12);
        // …and blind to the downstream pressure, even when it rises above
        // the collapse pressure (the waterfall decouples the segments).
        assert_eq!(wf.flow(0.012, 0.001), wf.flow(0.012, 0.005));
        // Below the threshold: collapsed, zero flow regardless of suction.
        assert_eq!(wf.flow(0.002, -1.0), 0.0);
        assert_eq!(wf.flow(0.001, 0.0), 0.0);
    }
}
