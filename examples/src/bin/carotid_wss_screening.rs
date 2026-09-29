//! Phase 4 milestone: carotid-like waveform WSS/OSI screening.
//!
//! Runs the voxel CFD solver on a stenosed tube under a pulsatile
//! carotid-like flow waveform, accumulating wall shear stress over the
//! cycle to produce mean WSS and OSI (oscillatory shear index). Post-
//! stenotic low-WSS / high-OSI regions are the classic atheroprone
//! signature this screening targets.
//!
//! ```console
//! cargo run -p tpt-med-examples --bin carotid-wss-screening
//! ```

use tpt_med_hemodynamics::{
    BloodModel, FluidDomain, HemodynamicsSolver, OsiAccumulator, PressureSolver, SolverConfig,
    WssField,
};

fn main() {
    println!("=== tpt-medical carotid-like WSS/OSI screening (Phase 4 milestone) ===");

    // Domain: tube with a 40% diameter stenosis built from a radius
    // profile along x: r(x) = R0·(1 − 0.4·exp(−((x−xm)/3)²)).
    let nx = 24usize;
    let nr = 12usize;
    let spacing = 0.5f64; // mm
    let r0 = 3.5f64; // cells
    let xm = nx as f64 / 2.0;

    let mut mask = vec![false; nx * nr * nr];
    let c1 = (nr as f64 - 1.0) / 2.0;
    for i in 0..nx {
        let stenosis = 1.0 - 0.4 * (-(((i as f64 + 0.5) - xm) / 3.0).powi(2)).exp();
        let r_local = r0 * stenosis;
        for j in 0..nr {
            for k in 0..nr {
                let radial = (((j as f64 - c1).powi(2)) + ((k as f64 - c1).powi(2))).sqrt();
                mask[(i * nr + j) * nr + k] = radial <= r_local;
            }
        }
    }
    let domain = FluidDomain::from_mask((nx, nr, nr), (spacing, spacing, spacing), mask, 0, true);
    let (_, vol) = domain.fluid_stats();
    println!(
        "domain: {nx}x{nr}x{nr} voxels, pitch {spacing} mm, fluid volume {vol:.1} mm³, 40% stenosis"
    );

    let config = SolverConfig {
        dt: 5.0e-4,
        poisson_iterations: 400,
        include_convection: false,
        viscosity_relaxation: 0.2,
        density: 1.06e-3,
        pressure_solver: PressureSolver::default(),
    };
    let mut solver =
        HemodynamicsSolver::new(domain, BloodModel::CARREAU_YASUDA_BLOOD, 30.0, config);

    // Carotid-like pulsatile inflow: mean 30 mm/s with ±60% sinusoidal
    // pulsation at 1 Hz over a 0.5 s march (half cycle sampled).
    let cycle = 1.0f64;
    let mut osi = OsiAccumulator::default();
    let mut wss_history: Vec<WssField> = Vec::new();
    println!("marching pulsatile cycle (shear-thinning blood) ...");
    for step in 0..600 {
        let t = step as f64 * solver.config.dt;
        let phase = 2.0 * core::f64::consts::PI * t / cycle;
        solver.inlet_velocity = 30.0 * (1.0 + 0.6 * phase.sin());
        solver.step();
        if step % 20 == 0 {
            let field = tpt_med_hemodynamics::extract_wss(&solver);
            for tr in &field.tractions {
                osi.sample(*tr, 20.0 * solver.config.dt);
            }
            wss_history.push(field);
        }
    }

    let last = wss_history.last().expect("samples");
    println!("--- results ---");
    println!("  instantaneous mean WSS: {:.3} Pa", last.mean_magnitude());
    println!("  instantaneous max  WSS: {:.3} Pa", last.max_magnitude());
    println!(
        "  cycle OSI (wall-averaged): {:.3} [0=steady, 0.5=reversed]",
        osi.osi()
    );

    // Stenosis screening: WSS at the throat plane should exceed the mean
    // (accelerated flow), the classic critical-plaque threshold reference
    // is ~15 Pa.
    let throat = wss_history
        .last()
        .expect("samples")
        .positions
        .iter()
        .enumerate()
        .filter(|(_, p)| (p.x - xm * spacing).abs() < spacing)
        .map(|(i, _)| last.tractions[i].norm())
        .collect::<Vec<_>>();
    let throat_mean = if throat.is_empty() {
        f64::NAN
    } else {
        throat.iter().sum::<f64>() / throat.len() as f64
    };
    println!(
        "  throat-plane mean WSS: {throat_mean:.3} Pa ({:.0}% of field mean)",
        100.0 * throat_mean / last.mean_magnitude()
    );
}
