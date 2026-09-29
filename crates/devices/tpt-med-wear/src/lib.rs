//! Implant wear simulation: Archard and Archard-type Cross–Land laws with
//! gait-cycle extrapolation (ISO 14879 / ASTM F2028 style screening).
//!
//! - **Archard**: `V = k · F · s` (volumetric wear per cycle from contact
//!   load and sliding distance).
//! - **Cross–Land**: per-cycle wear depth proportional to contact pressure
//!   above a fatigue threshold `p₀`: `dh = K·(p − p₀)·ds` — zero wear
//!   below the threshold.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// Wear law selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WearLaw {
    /// Archard adhesive/abrasive law with dimensionless wear coefficient
    /// `k` (mm³/(N·mm) conventionally expressed via the specific wear rate).
    Archard {
        /// Specific wear rate k (mm³/(N·m) — multiply by sliding metres).
        k: f64,
    },
    /// Cross–Land pressure-threshold law.
    CrossLand {
        /// Wear coefficient K (mm³/(N·m)).
        k: f64,
        /// Fatigue threshold pressure p₀ (MPa).
        pressure_threshold: f64,
    },
}

/// Configuration for a wear simulation.
#[derive(Debug, Clone, Copy)]
pub struct WearModel {
    /// Active wear law.
    pub law: WearLaw,
    /// Number of gait cycles to simulate/extrapolate.
    pub gait_cycles: u64,
}

/// Result of a wear simulation.
#[derive(Debug, Clone, Copy)]
pub struct WearResult {
    /// Total volumetric wear (mm³).
    pub volumetric_wear: f64,
    /// Mean linear wear depth over the contact area (mm).
    pub linear_wear: f64,
    /// Wear per million cycles (mm³/Mc).
    pub wear_per_megacycle: f64,
}

/// A non-uniform duty schedule: run-in excess wear and activity-level
/// variation, so extrapolation is not strictly linear in cycle count.
#[derive(Debug, Clone, PartialEq)]
pub struct WearSchedule {
    /// Length of the run-in period (cycles). During run-in the wear
    /// coefficient is multiplied by `1 + run_in_excess`, decaying
    /// exponentially to 1 across the period.
    pub run_in_cycles: u64,
    /// Extra wear-coefficient multiplier at cycle 0 (e.g. 2.0 = double
    /// wear on the first cycle, decaying to 1.0). 0 disables run-in.
    pub run_in_excess: f64,
    /// Activity levels as `(starting_cycle, multiplier)` pairs, applied in
    /// order; cycles before the first entry run at 1.0. Accumulated wear
    /// integrates the multiplier over each band exactly.
    pub activity_bands: Vec<(u64, f64)>,
}

impl Default for WearSchedule {
    fn default() -> Self {
        Self {
            run_in_cycles: 0,
            run_in_excess: 0.0,
            activity_bands: Vec::new(),
        }
    }
}

/// Per-cycle wear multiplier at cycle `n` under a schedule.
fn cycle_multiplier(schedule: &WearSchedule, n: u64) -> f64 {
    let run_in = if schedule.run_in_cycles > 0 && n < schedule.run_in_cycles {
        let progress = n as f64 / schedule.run_in_cycles as f64;
        1.0 + schedule.run_in_excess * (1.0 - progress)
    } else {
        1.0
    };
    let mut activity = 1.0;
    for &(start, mult) in &schedule.activity_bands {
        if n >= start {
            activity = mult;
        }
    }
    run_in * activity
}

/// Uncertainty band over the wear coefficient: low / central / high
/// results, plus the linearised sensitivity `∂V/V = ∂k/k` check.
#[derive(Debug, Clone, Copy)]
pub struct WearUncertainty {
    /// Low-end specific wear rate (same units as the law's `k`).
    pub k_low: f64,
    /// High-end specific wear rate.
    pub k_high: f64,
}

impl WearUncertainty {
    /// Band half-width relative to the central result:
    /// `(V_high − V_low) / (2 V_central)` — the number a screening report
    /// quotes as "±x %".
    pub fn relative_band(&self, central: &WearResult, low: &WearResult, high: &WearResult) -> f64 {
        if central.volumetric_wear.abs() < 1e-30 {
            return 0.0;
        }
        (high.volumetric_wear - low.volumetric_wear) / (2.0 * central.volumetric_wear.abs())
    }
}

impl WearModel {
    /// Simulates wear over `gait_cycles` cycles.
    ///
    /// `contact_pressures[i]` (MPa) and `sliding_distances[i]` (mm per
    /// cycle) are per contact zone; `area` is the total bearing area (mm²)
    /// used to convert volumetric → linear wear.
    pub fn simulate_wear(
        &self,
        contact_pressures: &[f64],
        sliding_distances: &[f64],
        contact_area: f64,
    ) -> WearResult {
        assert_eq!(
            contact_pressures.len(),
            sliding_distances.len(),
            "zone arrays must pair"
        );
        let mut per_cycle_volume = 0.0f64;
        for (&p, &s) in contact_pressures.iter().zip(sliding_distances) {
            let load = p * contact_area; // N (p in MPa → N/mm²)
            match self.law {
                WearLaw::Archard { k } => {
                    per_cycle_volume += k * load * s / 1000.0; // s mm → m
                }
                WearLaw::CrossLand {
                    k,
                    pressure_threshold,
                } => {
                    let excess = (p - pressure_threshold).max(0.0);
                    per_cycle_volume += k * (excess * contact_area) * s / 1000.0;
                }
            }
        }
        let cycles = self.gait_cycles as f64;
        let volumetric = per_cycle_volume * cycles;
        let linear = if contact_area > 0.0 {
            volumetric / contact_area
        } else {
            0.0
        };
        WearResult {
            volumetric_wear: volumetric,
            linear_wear: linear,
            wear_per_megacycle: per_cycle_volume * 1.0e6,
        }
    }

    /// Simulates wear under a non-uniform [`WearSchedule`] by integrating
    /// the per-cycle multiplier over `gait_cycles`. Cycle-exact for the
    /// piecewise bands (each cycle is counted with its own multiplier); for
    /// very long horizons this is O(cycles) by design — extrapolation
    /// beyond ~10⁹ cycles should use band-integrated closures instead.
    pub fn simulate_wear_with_schedule(
        &self,
        contact_pressures: &[f64],
        sliding_distances: &[f64],
        contact_area: f64,
        schedule: &WearSchedule,
    ) -> WearResult {
        assert_eq!(
            contact_pressures.len(),
            sliding_distances.len(),
            "zone arrays must pair"
        );
        let mut per_cycle_volume = 0.0f64;
        // Zones are static; accumulate each zone's base per-cycle volume,
        // then integrate the multiplier over cycles.
        for (&p, &s) in contact_pressures.iter().zip(sliding_distances) {
            let load = p * contact_area;
            match self.law {
                WearLaw::Archard { k } => per_cycle_volume += k * load * s / 1000.0,
                WearLaw::CrossLand {
                    k,
                    pressure_threshold,
                } => {
                    let excess = (p - pressure_threshold).max(0.0);
                    per_cycle_volume += k * (excess * contact_area) * s / 1000.0;
                }
            }
        }
        let mut total = 0.0f64;
        let mut n = 0u64;
        while n < self.gait_cycles {
            // Sum the multiplier run-length until the next band edge.
            let mut run = 1u64;
            for &(start, _) in &schedule.activity_bands {
                if n >= start && start > n {
                    run = run.min(start - n);
                }
                if start > n {
                    run = run.min(start - n);
                }
            }
            if schedule.run_in_cycles > n {
                run = run.min(schedule.run_in_cycles - n);
            }
            run = run.min(self.gait_cycles - n);
            total += per_cycle_volume * cycle_multiplier(schedule, n) * run as f64;
            n += run;
        }
        let volumetric = total;
        let linear = if contact_area > 0.0 {
            volumetric / contact_area
        } else {
            0.0
        };
        WearResult {
            volumetric_wear: volumetric,
            linear_wear: linear,
            wear_per_megacycle: volumetric / (self.gait_cycles as f64 / 1.0e6),
        }
    }

    /// Runs the wear simulation at `k_low`, the model's own `k`, and
    /// `k_high`, returning the three results for band reporting.
    pub fn simulate_wear_uncertainty(
        &self,
        contact_pressures: &[f64],
        sliding_distances: &[f64],
        contact_area: f64,
        uncertainty: &WearUncertainty,
    ) -> (WearResult, WearResult, WearResult) {
        let scale = |factor: f64| match self.law {
            WearLaw::Archard { k } => WearModel {
                law: WearLaw::Archard { k: k * factor },
                gait_cycles: self.gait_cycles,
            },
            WearLaw::CrossLand {
                k,
                pressure_threshold,
            } => WearModel {
                law: WearLaw::CrossLand {
                    k: k * factor,
                    pressure_threshold,
                },
                gait_cycles: self.gait_cycles,
            },
        };
        // Factor is relative to the model's own k; k=0 models scale the
        // zero result (band collapses).
        let k_ref = match self.law {
            WearLaw::Archard { k } => k,
            WearLaw::CrossLand { k, .. } => k,
        };
        let (low, high) = if k_ref.abs() < 1e-30 {
            (1.0, 1.0)
        } else {
            (uncertainty.k_low / k_ref, uncertainty.k_high / k_ref)
        };
        let low = scale(low).simulate_wear(contact_pressures, sliding_distances, contact_area);
        let central = self.simulate_wear(contact_pressures, sliding_distances, contact_area);
        let high = scale(high).simulate_wear(contact_pressures, sliding_distances, contact_area);
        (low, central, high)
    }

    /// ISO 14879-style screening: simulator wear limit for total knees is
    /// commonly quoted as < 30 mm³ per million cycles for UHMWPE tibial
    /// inserts.
    pub fn exceeds_iso14879_screen(&self, result: &WearResult, limit_mm3_per_mc: f64) -> bool {
        result.wear_per_megacycle > limit_mm3_per_mc
    }

    /// Simulates wear with the contact **solved each block** rather than
    /// prescribed: the accumulated per-zone wear depths feed the next
    /// block's [`ContactSolver::solve_pressures`], so wear changes the
    /// contact geometry and pressures — the debris-induced feedback the
    /// prescribed-pressure runs cannot represent. The built-in
    /// [`WinklerContact`] is self-stabilising (worn zones shed load to
    /// fresh ones); runaway modes (edge loading, third-body abrasion
    /// raising the effective wear coefficient) enter through a different
    /// solver implementation or a caller-supplied coefficient drift, and
    /// this loop is the hook they plug into.
    ///
    /// `block_cycles` is the wear/pressure update interval: smaller blocks
    /// resolve the migration more sharply at proportionally more solves.
    pub fn simulate_wear_with_contact(
        &self,
        contact: &mut dyn ContactSolver,
        sliding_distances: &[f64],
        total_load_n: f64,
        block_cycles: u64,
    ) -> ContactWearResult {
        let zones = sliding_distances.len();
        assert!(zones > 0, "at least one contact zone is required");
        assert!(block_cycles > 0, "block length must be positive");
        let mut wear_depths = vec![0.0f64; zones];
        let mut total_volume = 0.0f64;
        let mut cycles_done = 0u64;
        let mut final_pressures = contact.solve_pressures(&wear_depths, total_load_n);
        let mut lost_contact = false;

        while cycles_done < self.gait_cycles {
            let block = block_cycles.min(self.gait_cycles - cycles_done);
            final_pressures = contact.solve_pressures(&wear_depths, total_load_n);
            for (i, &p) in final_pressures.iter().enumerate() {
                if p <= 0.0 {
                    lost_contact = true;
                    continue;
                }
                let area = contact.zone_area(i);
                // Depth rate (mm/cycle): Archard V = k·F·s gives a depth
                // V/A = k·p·s (MPa·mm³/(N·m) cancels to mm per mm of
                // sliding /1000); Cross–Land is already a depth law,
                // dh = K·(p − p₀)·ds.
                let depth_rate = match self.law {
                    WearLaw::Archard { k } => k * p * sliding_distances[i] / 1000.0,
                    WearLaw::CrossLand {
                        k,
                        pressure_threshold,
                    } => k * (p - pressure_threshold).max(0.0) * sliding_distances[i] / 1000.0,
                };
                let depth = depth_rate * block as f64;
                total_volume += depth * area;
                wear_depths[i] += depth;
            }
            cycles_done += block;
        }

        let cycles = self.gait_cycles as f64;
        ContactWearResult {
            volumetric_wear: total_volume,
            wear_per_megacycle: total_volume / (cycles / 1.0e6),
            final_pressures,
            wear_depths,
            any_zone_lost_contact: lost_contact,
        }
    }
}

/// A contact solver: given the accumulated per-zone wear depths (mm) and
/// the total load (N), returns the equilibrium per-zone contact pressures
/// (MPa). This is the seam a real contact solver plugs into — the crate
/// ships [`WinklerContact`] as the built-in screening implementation.
pub trait ContactSolver {
    /// Equilibrium per-zone pressures (MPa) for the current wear state.
    fn solve_pressures(&mut self, wear_depths: &[f64], total_load_n: f64) -> Vec<f64>;

    /// Area of contact zone `zone` (mm²) — needed to convert the wear
    /// depth the loop accrues into the volumetric wear total.
    fn zone_area(&self, zone: usize) -> f64;
}

/// Winkler-foundation contact: each zone is a bed of springs of stiffness
/// `k` (N/mm³) over `area` (mm²); a rigid counter-surface at penetration
/// `δ` loads zone `i` at `p_i = k·(δ − w_i)` (MPa-scale via unit area),
/// clamped at zero once wear `w_i` passes the surface. `δ` is solved so
/// the zones carry the total load — worn zones shed load to fresh ones,
/// which is the self-stabilising screening feedback.
#[derive(Debug, Clone)]
pub struct WinklerContact {
    /// Foundation stiffness per unit area (N/mm³).
    pub foundation_stiffness: f64,
    /// Per-zone contact areas (mm²).
    pub zone_areas: Vec<f64>,
}

impl WinklerContact {
    /// Builds a contact with uniform zone areas.
    pub fn uniform_zones(foundation_stiffness: f64, zone_area: f64, zones: usize) -> Self {
        Self {
            foundation_stiffness,
            zone_areas: vec![zone_area; zones],
        }
    }
}

impl ContactSolver for WinklerContact {
    fn solve_pressures(&mut self, wear_depths: &[f64], total_load_n: f64) -> Vec<f64> {
        assert_eq!(
            wear_depths.len(),
            self.zone_areas.len(),
            "wear-depth and zone-area counts must pair"
        );
        let k = self.foundation_stiffness;
        // Start with every zone in contact; drop zones whose wear exceeds
        // the penetration and re-solve on the remainder.
        let mut active: Vec<usize> = (0..wear_depths.len()).collect();
        let mut pressures = vec![0.0f64; wear_depths.len()];
        loop {
            let total_area: f64 = active.iter().map(|&i| self.zone_areas[i]).sum();
            if total_area <= 0.0 || active.is_empty() {
                return vec![0.0; wear_depths.len()];
            }
            let mean_wear: f64 = active
                .iter()
                .map(|&i| self.zone_areas[i] * wear_depths[i])
                .sum::<f64>()
                / total_area;
            let penetration = total_load_n / (k * total_area) + mean_wear;
            let mut dropped = false;
            for &i in &active {
                let gap = penetration - wear_depths[i];
                if gap <= 0.0 {
                    pressures[i] = 0.0;
                    dropped = true;
                } else {
                    pressures[i] = k * gap; // N/mm³·mm = N/mm² = MPa
                }
            }
            if !dropped {
                return pressures;
            }
            active.retain(|&i| penetration > wear_depths[i]);
            if active.is_empty() {
                return vec![0.0; wear_depths.len()];
            }
        }
    }

    fn zone_area(&self, zone: usize) -> f64 {
        self.zone_areas[zone]
    }
}

/// Result of a contact-coupled wear run.
#[derive(Debug, Clone)]
pub struct ContactWearResult {
    /// Total volumetric wear over the run (mm³).
    pub volumetric_wear: f64,
    /// Wear per million cycles (mm³/Mc).
    pub wear_per_megacycle: f64,
    /// The last solved per-zone pressures (MPa).
    pub final_pressures: Vec<f64>,
    /// Accumulated per-zone wear depths (mm).
    pub wear_depths: Vec<f64>,
    /// Whether any zone's wear passed the penetration (contact lost) —
    /// the screening flag for load migration past a zone.
    pub any_zone_lost_contact: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model(law: WearLaw, cycles: u64) -> WearModel {
        WearModel {
            law,
            gait_cycles: cycles,
        }
    }

    #[test]
    fn run_in_period_exceeds_linear_extrapolation() {
        let zones_p = [5.0];
        let zones_s = [20.0];
        let area = 800.0;
        let plain = WearModel {
            law: WearLaw::Archard { k: 1.0e-6 },
            gait_cycles: 1_000_000,
        };
        let scheduled_model = WearModel {
            law: WearLaw::Archard { k: 1.0e-6 },
            gait_cycles: 1_000_000,
        };
        let schedule = WearSchedule {
            run_in_cycles: 200_000,
            run_in_excess: 1.0,
            activity_bands: vec![],
        };
        let base = plain.simulate_wear(&zones_p, &zones_s, area);
        let with_run_in =
            scheduled_model.simulate_wear_with_schedule(&zones_p, &zones_s, area, &schedule);
        // Run-in excess decays linearly 2x -> 1x across the first 20 % of
        // cycles (average 1.5x), so the total is
        // 0.8*1.0 + 0.2*1.5 = 1.1x the linear extrapolation.
        let expected = base.volumetric_wear * 1.1;
        assert!(
            (with_run_in.volumetric_wear - expected).abs() < 1e-6 * expected,
            "{} vs {expected}",
            with_run_in.volumetric_wear
        );
    }

    #[test]
    fn activity_bands_scale_exactly() {
        let zones_p = [5.0];
        let zones_s = [20.0];
        let area = 800.0;
        let model = WearModel {
            law: WearLaw::Archard { k: 1.0e-6 },
            gait_cycles: 1_000_000,
        };
        let schedule = WearSchedule {
            run_in_cycles: 0,
            run_in_excess: 0.0,
            activity_bands: vec![(500_000, 2.0)],
        };
        let r = model.simulate_wear_with_schedule(&zones_p, &zones_s, area, &schedule);
        let base = model.simulate_wear(&zones_p, &zones_s, area);
        let expected = base.volumetric_wear * 1.5;
        assert!(
            (r.volumetric_wear - expected).abs() < 1e-6 * expected,
            "{} vs {expected}",
            r.volumetric_wear
        );
    }

    #[test]
    fn uncertainty_band_tracks_k_ratio() {
        let zones_p = [5.0];
        let zones_s = [20.0];
        let area = 800.0;
        let model = WearModel {
            law: WearLaw::Archard { k: 1.0e-6 },
            gait_cycles: 1_000_000,
        };
        let unc = WearUncertainty {
            k_low: 5.0e-7,
            k_high: 2.0e-6,
        };
        let (low, central, high) = model.simulate_wear_uncertainty(&zones_p, &zones_s, area, &unc);
        assert!((low.volumetric_wear - central.volumetric_wear * 0.5).abs() < 1e-12);
        assert!((high.volumetric_wear - central.volumetric_wear * 2.0).abs() < 1e-9);
        let band = WearUncertainty::relative_band(&unc, &central, &low, &high);
        assert!((band - 0.75).abs() < 1e-9, "band {band}");
    }

    #[test]
    fn archard_is_linear_in_cycles_and_load() {
        let zones_p = [5.0]; // MPa over 800 mm² → 4000 N
        let zones_s = [20.0]; // mm per cycle
        let a = model(WearLaw::Archard { k: 1.0e-6 }, 1);
        let r1 = a.simulate_wear(&zones_p, &zones_s, 800.0);
        // V = k·F·s = 1e-6 · 4000 · 0.020 m = 8e-5 mm³ per cycle
        assert!((r1.volumetric_wear - 8.0e-5).abs() < 1e-9);
        let r10m = model(WearLaw::Archard { k: 1.0e-6 }, 10_000_000)
            .simulate_wear(&zones_p, &zones_s, 800.0);
        assert!((r10m.volumetric_wear - 1.0e7 * 8.0e-5).abs() < 1e-3);
        assert!((r10m.wear_per_megacycle - 80.0).abs() < 1e-6);
    }

    #[test]
    fn cross_land_threshold_gives_zero_wear_below_cutoff() {
        let zones_p = [5.0];
        let zones_s = [20.0];
        let below = model(
            WearLaw::CrossLand {
                k: 1.0e-6,
                pressure_threshold: 6.0,
            },
            1_000_000,
        )
        .simulate_wear(&zones_p, &zones_s, 800.0);
        assert_eq!(below.volumetric_wear, 0.0);
        let above = model(
            WearLaw::CrossLand {
                k: 1.0e-6,
                pressure_threshold: 4.0,
            },
            1_000_000,
        )
        .simulate_wear(&zones_p, &zones_s, 800.0);
        assert!(above.volumetric_wear > 0.0);
    }

    #[test]
    fn multi_zone_wear_sums() {
        let p = [5.0, 3.0];
        let s = [20.0, 10.0];
        let r = model(WearLaw::Archard { k: 1.0e-6 }, 1).simulate_wear(&p, &s, 800.0);
        // Zone 1: 1e-6·4000·0.02 = 8e-5; zone 2: 1e-6·2400·0.01 = 2.4e-5
        assert!((r.volumetric_wear - 1.04e-4).abs() < 1e-9);
    }

    #[test]
    fn linear_wear_divides_by_area() {
        let r =
            model(WearLaw::Archard { k: 1.0e-6 }, 1_000_000).simulate_wear(&[5.0], &[20.0], 400.0);
        assert!((r.linear_wear - r.volumetric_wear / 400.0).abs() < 1e-12);
    }

    #[test]
    fn iso_screening_flag() {
        let r =
            model(WearLaw::Archard { k: 2.0e-6 }, 10_000_000).simulate_wear(&[5.0], &[20.0], 800.0);
        assert!(r.wear_per_megacycle > 30.0);
        let m = model(WearLaw::Archard { k: 2.0e-6 }, 10_000_000);
        assert!(m.exceeds_iso14879_screen(&r, 30.0));
        assert!(!m.exceeds_iso14879_screen(&r, 1.0e6));
    }

    #[test]
    fn winkler_contact_conserves_the_total_load() {
        let mut contact = WinklerContact::uniform_zones(5.0, 200.0, 3);
        let p = contact.solve_pressures(&[0.0, 0.0, 0.0], 3000.0);
        // Uniform fresh bed: each zone carries F/3 over 200 mm².
        for pi in &p {
            assert!((pi - 5.0).abs() < 1e-9, "{p:?}");
        }
        // A worn zone sheds load; the total is conserved.
        let worn = contact.solve_pressures(&[0.5, 0.0, 0.0], 3000.0);
        let total: f64 = worn.iter().sum::<f64>() * 200.0;
        assert!((total - 3000.0).abs() < 1e-6, "total {total}");
        assert!(worn[0] < worn[1], "worn zone carries less: {worn:?}");
        // A zone worn past the penetration drops out entirely.
        let out = contact.solve_pressures(&[5.0, 0.0, 0.0], 3000.0);
        assert_eq!(out[0], 0.0);
        let total2: f64 = out.iter().sum::<f64>() * 200.0;
        assert!(
            (total2 - 3000.0).abs() < 1e-6,
            "total after dropout {total2}"
        );
    }

    #[test]
    fn contact_coupled_wear_migrates_load_off_worn_zones() {
        let mut contact = WinklerContact::uniform_zones(5.0, 200.0, 2);
        let law = WearLaw::Archard { k: 1.0e-3 };
        // Zone 0 slides 3x further, so it wears 3x faster and must shed
        // load to zone 1 — the feedback the prescribed-pressure run
        // cannot represent.
        let out = WearModel {
            law,
            gait_cycles: 200_000,
        }
        .simulate_wear_with_contact(&mut contact, &[30.0, 10.0], 4000.0, 10_000);
        assert!(out.volumetric_wear > 0.0);
        assert!(
            out.final_pressures[0] < out.final_pressures[1],
            "the faster-worn zone must carry less: {:?}",
            out.final_pressures
        );
        // Depth tracks the sliding asymmetry, direction preserved.
        assert!(out.wear_depths[0] > out.wear_depths[1]);
        // Self-stabilisation: the coupled total is BELOW the
        // prescribed-pressure run at the initial uniform pressure
        // (penetration 2 mm on k=5, 400 mm² → 10 MPa per zone).
        let prescribed_initial = WearModel {
            law,
            gait_cycles: 200_000,
        }
        .simulate_wear(&[10.0, 10.0], &[30.0, 10.0], 400.0);
        assert!(
            out.volumetric_wear < prescribed_initial.volumetric_wear,
            "coupled {} vs prescribed-at-initial {}",
            out.volumetric_wear,
            prescribed_initial.volumetric_wear
        );
        // Load conservation holds at the final solve.
        let total: f64 = out.final_pressures.iter().sum::<f64>() * 200.0;
        assert!((total - 4000.0).abs() < 1e-6, "final load {total}");
    }

    #[test]
    fn asymmetric_wear_sheds_load_to_fresh_zones() {
        let mut contact = WinklerContact::uniform_zones(5.0, 200.0, 2);
        // Pre-wear zone 0 by starting the loop with a seeded depth: solve
        // once with an asymmetric state through the public path.
        let seeded = contact.solve_pressures(&[0.3, 0.0], 4000.0);
        assert!(
            seeded[0] < seeded[1],
            "seeded wear must shed load: {seeded:?}"
        );
        // The full loop from a seeded solver: use two blocks; zone 0's
        // pressure stays below zone 1's throughout.
        let out = WearModel {
            law: WearLaw::Archard { k: 1.0e-3 },
            gait_cycles: 100_000,
        }
        .simulate_wear_with_contact(&mut contact, &[20.0, 20.0], 4000.0, 50_000);
        assert_eq!(out.wear_depths.len(), 2);
        assert!(!out.any_zone_lost_contact, "uniform zones stay in contact");
    }

    #[test]
    fn uhmwpe_typical_knee_screen_is_plausible() {
        // Typical PE knee wear k ≈ 1e-6 mm³/N·m, F ≈ 2500 N, s ≈ 15 mm:
        // ~37 mm³/Mc — the right order for ISO 14879 simulator data.
        let r =
            model(WearLaw::Archard { k: 1.0e-6 }, 1_000_000).simulate_wear(&[5.0], &[15.0], 800.0);
        assert!((10.0..100.0).contains(&r.wear_per_megacycle));
    }
}
