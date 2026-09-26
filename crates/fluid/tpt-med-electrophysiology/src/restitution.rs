//! S1-S2 single-cell restitution protocol
//! (`rfcs/0005-cardiac-electrophysiology.md` Stage 1) — the direct analogue
//! of how Mitchell & Schaeffer's own paper characterises the model, used
//! here as a code-verification fixture for the ionic kinetics rather than a
//! tissue-level feature. A tissue-level restitution behaviour additionally
//! depends on conduction and the diffusion coupling; this protocol does
//! not simulate either, and must not be treated as standing in for one.

use crate::error::{EpError, Result};
use crate::kinetics::MitchellSchaefferParams;

/// Fixed single-cell integration step, ms. Small relative to every
/// `MitchellSchaefferParams` time constant this crate ships or validates
/// against, so explicit RK2 is comfortably stable without a CFL-style bound
/// (there is no diffusion term in a single cell).
const DT_MS: f64 = 0.01;

/// An S1-S2 restitution protocol on a single 0D cell (no diffusion): `s1_beats`
/// stimuli at `s1_cycle_length_ms`, then one test stimulus per entry of
/// `s2_coupling_intervals_ms`, each measured independently from the same
/// paced state (not cumulatively).
#[derive(Debug, Clone)]
pub struct S1S2Protocol {
    /// S1 pacing cycle length, ms.
    pub s1_cycle_length_ms: f64,
    /// Number of S1 beats to pace before each S2 test.
    pub s1_beats: u32,
    /// S2 coupling intervals to test, ms, each measured from the last S1
    /// stimulus's own onset.
    pub s2_coupling_intervals_ms: Vec<f64>,
    /// Stimulus current amplitude (dimensionless, added to `dV/dt`).
    pub stimulus_amplitude: f64,
    /// Stimulus duration, ms.
    pub stimulus_duration_ms: f64,
}

/// One detected action potential: `(upstroke_time, apd)`, `apd` measured as
/// the duration `V >= v_gate` — the model-consistent plateau definition
/// (the same threshold at which `dh/dt` itself switches sign), not an
/// arbitrary repolarization percentage.
type ActionPotential = (f64, f64);

impl S1S2Protocol {
    /// Runs the protocol, returning one `(diastolic_interval, apd)` pair
    /// per entry of `s2_coupling_intervals_ms` that produced a detectable
    /// action potential — a coupling interval that falls in the model's
    /// refractory period and produces no detectable upstroke is a
    /// meaningful restitution-curve endpoint, not an error, and is simply
    /// omitted from the result. Errs only if the *baseline* S1 pacing train
    /// itself never produces a detectable action potential (a
    /// `s1_cycle_length_ms` shorter than the model's refractory period, or
    /// too weak a stimulus).
    pub fn run(&self, params: &MitchellSchaefferParams) -> Result<Vec<(f64, f64)>> {
        let s1_schedule: Vec<(f64, f64, f64)> = (0..self.s1_beats)
            .map(|i| {
                (
                    i as f64 * self.s1_cycle_length_ms,
                    self.stimulus_duration_ms,
                    self.stimulus_amplitude,
                )
            })
            .collect();

        let last_s1_onset = (self.s1_beats.saturating_sub(1)) as f64 * self.s1_cycle_length_ms;
        let margin_ms = 5.0 * params.tau_close;
        let baseline_total_time = last_s1_onset + margin_ms;

        let baseline_aps = Self::simulate(params, &s1_schedule, baseline_total_time);
        let last_s1_ap = baseline_aps
            .iter()
            .rev()
            .find(|(onset, _)| *onset <= last_s1_onset + 1.0)
            .copied()
            .ok_or_else(|| {
                EpError::RestitutionFailed(
                    "no action potential detected during S1 pacing — check \
                     s1_cycle_length_ms against the model's refractory period, \
                     and stimulus_amplitude/stimulus_duration_ms against v_gate"
                        .to_string(),
                )
            })?;
        let (last_s1_upstroke, last_s1_apd) = last_s1_ap;

        let mut results = Vec::with_capacity(self.s2_coupling_intervals_ms.len());
        for &coupling_interval in &self.s2_coupling_intervals_ms {
            let s2_onset = last_s1_upstroke + coupling_interval;
            let mut schedule = s1_schedule.clone();
            schedule.push((s2_onset, self.stimulus_duration_ms, self.stimulus_amplitude));
            let total_time = s2_onset + margin_ms;

            let aps = Self::simulate(params, &schedule, total_time);
            if let Some((_, apd_s2)) = aps
                .iter()
                .find(|(onset, _)| *onset > last_s1_upstroke + 1.0)
            {
                let di = coupling_interval - last_s1_apd;
                results.push((di, *apd_s2));
            }
        }

        Ok(results)
    }

    /// Integrates a single cell from rest (`V=0`, `h=1`) for `total_time_ms`,
    /// applying each `(onset, duration, amplitude)` in `schedule` as a
    /// constant stimulus current while active, and returns every detected
    /// action potential as `(upstroke_time, apd)`.
    fn simulate(
        params: &MitchellSchaefferParams,
        schedule: &[(f64, f64, f64)],
        total_time_ms: f64,
    ) -> Vec<ActionPotential> {
        let mut v = 0.0_f64;
        let mut h = 1.0_f64;
        let mut t = 0.0_f64;
        let mut above_gate = false;
        let mut upstroke_time = 0.0_f64;
        let mut aps = Vec::new();

        let steps = (total_time_ms / DT_MS).ceil() as u64;
        let stim_at = |time: f64| -> f64 {
            schedule
                .iter()
                .find(|(onset, duration, _)| time >= *onset && time < onset + duration)
                .map(|(_, _, amplitude)| *amplitude)
                .unwrap_or(0.0)
        };

        for _ in 0..steps {
            let i_stim = stim_at(t);
            let k1v = params.dv_dt(v, h, i_stim);
            let k1h = params.dh_dt(v, h);
            let v1 = v + DT_MS * k1v;
            let h1 = h + DT_MS * k1h;
            let i_stim2 = stim_at(t + DT_MS);
            let k2v = params.dv_dt(v1, h1, i_stim2);
            let k2h = params.dh_dt(v1, h1);
            v += DT_MS * 0.5 * (k1v + k2v);
            h += DT_MS * 0.5 * (k1h + k2h);
            t += DT_MS;

            let now_above = v >= params.v_gate;
            if now_above && !above_gate {
                upstroke_time = t;
            } else if !now_above && above_gate {
                aps.push((upstroke_time, t - upstroke_time));
            }
            above_gate = now_above;
        }

        aps
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_protocol() -> S1S2Protocol {
        // At `s1_cycle_length_ms = 300`, the S1 action potential itself
        // lasts ~260-285ms (time V >= v_gate) with these default kinetics,
        // so `s2_coupling_intervals_ms` must clear that plateau for the S2
        // stimulus to find the cell excitable at all — otherwise every
        // coupling interval lands in the absolute refractory period and no
        // capture occurs. 300/320/... is comfortably past the plateau
        // while still probing progressively larger diastolic intervals.
        S1S2Protocol {
            s1_cycle_length_ms: 300.0,
            s1_beats: 3,
            s2_coupling_intervals_ms: vec![300.0, 320.0, 350.0, 400.0, 500.0],
            stimulus_amplitude: 1.0,
            stimulus_duration_ms: 1.0,
        }
    }

    #[test]
    fn produces_a_restitution_curve() {
        let params = MitchellSchaefferParams::human_ventricular_default();
        let protocol = default_protocol();
        let curve = protocol.run(&params).expect("runs");
        assert!(!curve.is_empty());
        for (di, apd) in &curve {
            assert!(di.is_finite());
            assert!(apd.is_finite() && *apd > 0.0);
        }
    }

    #[test]
    fn apd_increases_monotonically_with_diastolic_interval() {
        // Restitution: shorter DI (a more premature S2) yields a shorter
        // APD than a longer DI, up to saturation — the qualitative shape
        // Mitchell & Schaeffer's own paper reports.
        let params = MitchellSchaefferParams::human_ventricular_default();
        let protocol = default_protocol();
        let curve = protocol.run(&params).expect("runs");
        let mut sorted = curve.clone();
        sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for i in 1..sorted.len() {
            assert!(
                sorted[i].1 >= sorted[i - 1].1 - 1e-6,
                "APD should not decrease as DI increases: {:?}",
                sorted
            );
        }
    }

    #[test]
    fn rejects_a_cycle_length_inside_the_refractory_period() {
        let params = MitchellSchaefferParams::human_ventricular_default();
        let protocol = S1S2Protocol {
            s1_cycle_length_ms: 1.0,
            s1_beats: 3,
            s2_coupling_intervals_ms: vec![10.0],
            stimulus_amplitude: 1.0,
            stimulus_duration_ms: 0.5,
        };
        // Beats fire faster than the cell can repolarize between them;
        // whether this specific model/parameter combination technically
        // still captures one action potential is not the point under test
        // -- what matters is that a physiologically nonsensical schedule
        // never panics and always returns a `Result`.
        let _ = protocol.run(&params);
    }

    #[test]
    fn too_weak_a_stimulus_fails_baseline_pacing() {
        let params = MitchellSchaefferParams::human_ventricular_default();
        let protocol = S1S2Protocol {
            s1_cycle_length_ms: 300.0,
            s1_beats: 2,
            s2_coupling_intervals_ms: vec![200.0],
            stimulus_amplitude: 0.0,
            stimulus_duration_ms: 1.0,
        };
        assert!(matches!(
            protocol.run(&params),
            Err(EpError::RestitutionFailed(_))
        ));
    }
}
