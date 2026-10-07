//! Crown-ring verification studies (RFC 0013 verification items 4 and 5).
//!
//! The 3-D ring is compressed by the radial wall in stages and its total
//! radial force (the sum of normal reactions, the same scalar Level 1 calls
//! `radial_force`) is recorded against the compression of the diameter.
//!
//! Two tiers, because the dense linear backend makes the realistic meshes a
//! release-profile job:
//!
//! - a **tiny ring** in the ordinary test suite (qualitative behaviour and
//!   the "pins carry no load" claim);
//! - **`#[ignore]`d studies** for the cross-check and the refinement trend,
//!   run with
//!   `cargo test --release -p tpt-med-fem-adapter --features superelastic
//!   --lib ring_studies -- --ignored --nocapture`.
//!
//! The recorded numbers live in the RFC and tracker, not in assertions that
//! would fail on a harmless re-meshing. The assertions pin what the studies
//! established: the section mesh is irrelevant (< 2 %), the centreline
//! refinement converges at roughly second order, and Level 1's linear law
//! agrees within 15 % over an elastic window (>= 6 % of diameter).
//!
//! **A lesson recorded where it will be seen:** a first, coarse ring
//! (4 stations per crown) showed a superelastic plateau from 3 % compression
//! and a 1 % agreement window. Both were artefacts — the faceted polyline is
//! a stiff polygon, so the struts were over-strained and transformed early.
//! Refined, the struts stay elastic to ~6 % and the ring is ~50x softer.
//! Never read material behaviour off a ring whose centreline is unresolved.

use crate::deployment::{Deployment3D, StageOptions};
use crate::fixtures::CrownRingSpec;
use crate::solver::{Convergence, SolveOptions};
use crate::superelastic::{SouzaAuricchio, SuperelasticParams};
use crate::wall::{RadialWall, WallSide};
use tpt_med_stents::{simulate_deployment, NitinolParams, StentModel};
use tpt_med_units::Pressure;

/// Normal penalty stiffness (N/mm): ≫ the ring's structural stiffness so the
/// penetration is a negligible fraction of the compression.
const KAPPA: f64 = 1.0e5;

/// One point of a radial-force curve.
#[derive(Debug, Clone, Copy)]
pub struct CurvePoint {
    /// Diameter compression as a fraction of the nominal outer diameter.
    pub fraction: f64,
    /// Compression in mm of diameter.
    pub compression: f64,
    /// Total radial force (N).
    pub force: f64,
    /// Largest martensite fraction in the ring.
    pub xi_max: f64,
    /// Axial/twist pin reaction magnitude (N) — should be ~0.
    pub pin_reaction: f64,
}

/// Compress the ring through `fractions` (ascending) and record the curve.
pub fn radial_curve(spec: &CrownRingSpec, fractions: &[f64]) -> Vec<CurvePoint> {
    let ring = spec.build().expect("ring");
    let r0 = spec.outer_radius();
    let model = SouzaAuricchio::new(SuperelasticParams::default()).expect("model");
    let opts = SolveOptions {
        convergence: Convergence {
            abs_tol: 1e-7,
            rel_tol: 1e-10,
            max_iter: 60,
        },
        ..SolveOptions::default()
    };
    let mut dep = Deployment3D::new(&ring.mesh, model, opts, &ring.pins).expect("dep");
    let wall = RadialWall::new(
        2,
        [0.0; 2],
        r0 + 1e-3,
        WallSide::Inside,
        ring.outer_nodes.clone(),
    )
    .expect("wall");
    dep.set_wall(wall, KAPPA);
    let pin_dofs: Vec<usize> = ring.pins.iter().map(|&(d, _)| d).collect();
    let stage = StageOptions {
        steps: 1,
        max_cutbacks: 10,
    };
    let mut out = Vec::new();
    for &frac in fractions {
        dep.move_wall(r0 * (1.0 - frac)).expect("move");
        dep.prescribe(&ring.pins, stage, None).expect("stage");
        let s = dep.wall_summary().expect("summary").expect("wall");
        let pin_reaction = pin_dofs
            .iter()
            .map(|&d| dep.reaction(&[d]).expect("reaction").abs())
            .fold(0.0, f64::max);
        out.push(CurvePoint {
            fraction: frac,
            compression: 2.0 * r0 * frac,
            force: s.total_reaction,
            xi_max: dep.field().martensite_summary().1,
            pin_reaction,
        });
    }
    out
}

/// Level 1's radial force for the same nominal ring, with the crown
/// stiffness calibrated from the 3-D curve's first point.
pub fn level_one_force(spec: &CrownRingSpec, calibration: &CurvePoint, compression: f64) -> f64 {
    let nominal = 2.0 * spec.outer_radius();
    let k = calibration.force / (spec.crowns as f64 * calibration.compression);
    let stent = StentModel {
        expanded_diameter: nominal,
        crimped_diameter: 0.5 * nominal,
        n_crowns: spec.crowns as u32,
        crown_stiffness: k,
    };
    let lumen = nominal - compression;
    simulate_deployment(
        &stent,
        &NitinolParams::default(),
        |_| lumen,
        Pressure::from_pa(0.0),
    )
    .radial_force
}

fn tiny() -> CrownRingSpec {
    let mut s = CrownRingSpec::new(4);
    s.stations_per_crown = 3;
    s.width_elements = 1;
    s.thickness_elements = 1;
    s
}

#[test]
fn tiny_ring_responds_radially_and_the_pins_carry_no_load() {
    let curve = radial_curve(&tiny(), &[0.005, 0.01, 0.02]);
    assert!(curve[0].force > 0.0);
    assert!(
        curve.windows(2).all(|w| w[1].force > w[0].force),
        "force must rise with compression in the elastic window: {curve:?}"
    );
    // The rigid-mode pins remove modes without carrying the radial load.
    for p in &curve {
        assert!(
            p.pin_reaction < 0.02 * p.force,
            "pin reaction {} vs force {}",
            p.pin_reaction,
            p.force
        );
    }
}

/// Observed convergence order and Richardson limit from three refinements.
///
/// `h` are the (decreasing) mesh sizes and `k` the corresponding values of a
/// monotonically converging quantity. Solves
/// `(k0 - k1) / (k1 - k2) = (h0^p - h1^p) / (h1^p - h2^p)` for `p` by
/// bisection and returns `(p, k_limit)`, or `None` when the three values are
/// not monotone (no honest order can be quoted — the same fallback the
/// peak-WSS resolution study uses).
pub fn observed_order(h: [f64; 3], k: [f64; 3]) -> Option<(f64, f64)> {
    let (d01, d12) = (k[0] - k[1], k[1] - k[2]);
    let prod = d01 * d12;
    if !prod.is_finite() || prod <= 0.0 || !h.windows(2).all(|w| w[0] > w[1]) {
        return None;
    }
    let target = d01 / d12;
    let ratio = |p: f64| (h[0].powf(p) - h[1].powf(p)) / (h[1].powf(p) - h[2].powf(p));
    // ratio(p) is increasing in p over the useful range.
    let (mut lo, mut hi) = (0.05f64, 12.0f64);
    if target < ratio(lo) || target > ratio(hi) {
        return None;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if ratio(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let p = 0.5 * (lo + hi);
    let c = d12 / (h[1].powf(p) - h[2].powf(p));
    Some((p, k[2] - c * h[2].powf(p)))
}

#[test]
fn observed_order_recovers_a_synthetic_power_law() {
    let h: [f64; 3] = [1.0 / 12.0, 1.0 / 16.0, 1.0 / 24.0];
    let k: Vec<f64> = h.iter().map(|x| 5.0 + 300.0 * x.powf(2.5)).collect();
    let (p, limit) = observed_order(h, [k[0], k[1], k[2]]).expect("monotone");
    assert!((p - 2.5).abs() < 1e-6, "{p}");
    assert!((limit - 5.0).abs() < 1e-6, "{limit}");
    // Non-monotone data yields no order.
    assert!(observed_order(h, [3.0, 4.0, 3.5]).is_none());
}

#[test]
#[ignore = "release-profile study; minutes; see module docs"]
fn section_sensitivity_is_small() {
    // Centreline fixed; vary the section mesh.
    let mut ks = Vec::new();
    for (nw, nt) in [(1usize, 2usize), (2, 2), (2, 3)] {
        let mut s = CrownRingSpec::new(6);
        s.stations_per_crown = 8;
        s.width_elements = nw;
        s.thickness_elements = nt;
        let c = radial_curve(&s, &[0.005]);
        let k = c[0].force / c[0].compression;
        println!("m=8 nw={nw} nt={nt}: k(0.5%) = {k:.3} N/mm");
        ks.push(k);
    }
    let spread = (ks.iter().cloned().fold(f64::MIN, f64::max)
        - ks.iter().cloned().fold(f64::MAX, f64::min))
        / ks[1];
    println!("section spread {:.3} %", 100.0 * spread);
    assert!(spread < 0.02, "section mesh moves the answer by {spread}");
}

#[test]
#[ignore = "release-profile study; many minutes; see module docs"]
fn centreline_refinement_order() {
    // The centreline resolution is the dominant discretisation error (the
    // faceted ring is stiffer than the smooth one), so it is the axis the
    // order is quoted on. Section fixed at 1 x 2 (section-insensitive).
    let stations = [12usize, 16, 24, 32];
    let mut ks = Vec::new();
    for m in stations {
        let mut s = CrownRingSpec::new(6);
        s.stations_per_crown = m;
        s.width_elements = 1;
        s.thickness_elements = 2;
        let t0 = std::time::Instant::now();
        let c = radial_curve(&s, &[0.005]);
        let k = c[0].force / c[0].compression;
        println!("m={m:2}: k(0.5%) = {k:9.4} N/mm   [{:?}]", t0.elapsed());
        ks.push(k);
    }
    assert!(
        ks.windows(2).all(|w| w[1] < w[0]),
        "stiffness must fall with refinement: {ks:?}"
    );
    let h = |m: usize| 1.0 / m as f64;
    let (p_c, lim_c) = observed_order([h(16), h(24), h(32)], [ks[1], ks[2], ks[3]])
        .expect("monotone finest triple");
    let (p_f, lim_f) = observed_order([h(12), h(16), h(24)], [ks[0], ks[1], ks[2]])
        .expect("monotone coarser triple");
    println!("order (16,24,32) = {p_c:.2}, limit {lim_c:.3}; order (12,16,24) = {p_f:.2}, limit {lim_f:.3}");
    println!(
        "finest-mesh error vs limit: {:.1} %",
        100.0 * (ks[3] - lim_c) / lim_c
    );
    assert!(p_c > 1.0, "order {p_c}");
}

#[test]
#[ignore = "release-profile study; many minutes; see module docs"]
fn level_one_cross_check_on_the_refined_ring() {
    // Two centreline resolutions, so the *ratio* to Level 1 (which cancels
    // the overall stiffness scale through calibration) can be seen to be, or
    // not be, mesh-sensitive even though the absolute force is not converged.
    for m in [16usize, 24] {
        let mut spec = CrownRingSpec::new(6);
        spec.stations_per_crown = m;
        spec.width_elements = 1;
        spec.thickness_elements = 2;
        let fractions = [0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08];
        let curve = radial_curve(&spec, &fractions);
        let cal = curve[0];
        println!(
            "m={m}: crown stiffness N/mm/crown = {:.4}",
            cal.force / (spec.crowns as f64 * cal.compression)
        );
        println!("frac  d(mm)   F3D(N)   F_L1(N)  ratio  xi_max  pin(N)");
        let mut window = 0.0;
        let mut in_window = true;
        for p in &curve {
            let l1 = level_one_force(&spec, &cal, p.compression);
            let ratio = p.force / l1;
            println!(
                "{:.3} {:.4} {:8.4} {:8.4} {:6.3} {:6.3} {:.2e}",
                p.fraction, p.compression, p.force, l1, ratio, p.xi_max, p.pin_reaction
            );
            if in_window && (ratio - 1.0).abs() <= 0.15 {
                window = p.fraction;
            } else {
                in_window = false;
            }
        }
        println!("m={m}: 15 % agreement window ends at {window:.3} of diameter");
        assert!(window >= 0.06, "m={m}: window {window}");
    }
}
