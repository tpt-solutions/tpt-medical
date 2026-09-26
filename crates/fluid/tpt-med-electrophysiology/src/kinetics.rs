//! Mitchell–Schaeffer two-variable ionic kinetics
//! (`rfcs/0005-cardiac-electrophysiology.md` Stage 1).
//!
//! Mitchell CC, Schaeffer DG. "A two-current model for the dynamics of
//! cardiac membrane." *Bull. Math. Biol.* 65:767–793, 2003.

use crate::error::{EpError, Result};

/// Mitchell–Schaeffer kinetic parameters: dimensionless voltage
/// `V ∈ [0, 1]`, gate `h ∈ [0, 1]`, time constants in ms.
///
/// ```text
/// dV/dt = J_in(V, h) + J_out(V) + J_stim
/// dh/dt = (1 - h) / tau_open   if V <  v_gate
///       = -h / tau_close       if V >= v_gate
///
/// J_in(V, h) = h * V^2 * (1 - V) / tau_in
/// J_out(V)   = -V / tau_out
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MitchellSchaefferParams {
    /// Upstroke (fast inward current) time constant, ms.
    pub tau_in: f64,
    /// Repolarization (outward current) time constant, ms.
    pub tau_out: f64,
    /// Gate-opening time constant (recovery while `V < v_gate`), ms.
    pub tau_open: f64,
    /// Gate-closing time constant (inactivation while `V >= v_gate`), ms.
    pub tau_close: f64,
    /// Gating threshold, dimensionless, in `(0, 1)`.
    pub v_gate: f64,
}

impl MitchellSchaefferParams {
    /// Validates and constructs a parameter set: every time constant must
    /// be finite and positive, `v_gate` finite and in `(0, 1)`.
    pub fn new(
        tau_in: f64,
        tau_out: f64,
        tau_open: f64,
        tau_close: f64,
        v_gate: f64,
    ) -> Result<Self> {
        for (name, value) in [
            ("tau_in", tau_in),
            ("tau_out", tau_out),
            ("tau_open", tau_open),
            ("tau_close", tau_close),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(EpError::NonFiniteParameter { name });
            }
        }
        if !v_gate.is_finite() || !(v_gate > 0.0 && v_gate < 1.0) {
            return Err(EpError::NonFiniteParameter { name: "v_gate" });
        }
        Ok(Self {
            tau_in,
            tau_out,
            tau_open,
            tau_close,
            v_gate,
        })
    }

    /// A commonly-reproduced parameter fit from Mitchell & Schaeffer (2003)
    /// Table 1's generic ventricular case: `tau_in = 0.3 ms`,
    /// `tau_out = 6.0 ms`, `tau_open = 120.0 ms`, `tau_close = 150.0 ms`,
    /// `v_gate = 0.13`.
    ///
    /// **Screening default, not a validated patient- or species-specific
    /// fit** — see `rfcs/0005-cardiac-electrophysiology.md`'s Motivation for
    /// the ASME V&V 40 "supporting" influence level this default carries.
    /// Verify against the primary source before using in any published
    /// comparison.
    pub fn human_ventricular_default() -> Self {
        Self {
            tau_in: 0.3,
            tau_out: 6.0,
            tau_open: 120.0,
            tau_close: 150.0,
            v_gate: 0.13,
        }
    }

    /// `dV/dt`'s reaction term (no diffusion), given a stimulus current.
    pub(crate) fn dv_dt(&self, v: f64, h: f64, i_stim: f64) -> f64 {
        let j_in = h * v * v * (1.0 - v) / self.tau_in;
        let j_out = -v / self.tau_out;
        j_in + j_out + i_stim
    }

    /// `dh/dt`.
    pub(crate) fn dh_dt(&self, v: f64, h: f64) -> f64 {
        if v < self.v_gate {
            (1.0 - h) / self.tau_open
        } else {
            -h / self.tau_close
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_non_positive_time_constants() {
        assert!(MitchellSchaefferParams::new(0.0, 6.0, 120.0, 150.0, 0.13).is_err());
        assert!(MitchellSchaefferParams::new(-0.3, 6.0, 120.0, 150.0, 0.13).is_err());
        assert!(MitchellSchaefferParams::new(0.3, 6.0, 120.0, 150.0, f64::NAN).is_err());
    }

    #[test]
    fn rejects_v_gate_out_of_range() {
        assert!(MitchellSchaefferParams::new(0.3, 6.0, 120.0, 150.0, 0.0).is_err());
        assert!(MitchellSchaefferParams::new(0.3, 6.0, 120.0, 150.0, 1.0).is_err());
        assert!(MitchellSchaefferParams::new(0.3, 6.0, 120.0, 150.0, 1.5).is_err());
    }

    #[test]
    fn accepts_the_default_parameters() {
        let p = MitchellSchaefferParams::human_ventricular_default();
        assert!(MitchellSchaefferParams::new(
            p.tau_in,
            p.tau_out,
            p.tau_open,
            p.tau_close,
            p.v_gate
        )
        .is_ok());
    }

    #[test]
    fn resting_state_relaxes_h_toward_one() {
        let p = MitchellSchaefferParams::human_ventricular_default();
        // At rest (V=0 < v_gate), h moves toward 1.
        assert!(p.dh_dt(0.0, 0.5) > 0.0);
        // Fully recovered at rest: no further change.
        assert_eq!(p.dh_dt(0.0, 1.0), 0.0);
    }

    #[test]
    fn plateau_state_relaxes_h_toward_zero() {
        let p = MitchellSchaefferParams::human_ventricular_default();
        // Above v_gate, h decays toward 0.
        assert!(p.dh_dt(0.9, 0.5) < 0.0);
    }
}
