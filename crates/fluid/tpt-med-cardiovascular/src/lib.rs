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
}
