//! Hemodynamics CFD benchmark (steady Poiseuille march).
//!
//! Run: `cargo bench -p tpt-med-benches --bench hemodynamics-cfd`

use std::time::Instant;

use tpt_med_hemodynamics::{BloodModel, FluidDomain, HemodynamicsSolver, SolverConfig};

fn main() {
    let domain = FluidDomain::cylinder(24, 12, 4.5, 0.5, 0);
    let config = SolverConfig {
        dt: 5.0e-4,
        poisson_iterations: 250,
        include_convection: false,
        viscosity_relaxation: 0.2,
        density: 1.06e-3,
    };
    let mut solver = HemodynamicsSolver::new(
        domain,
        BloodModel::Newtonian { viscosity: 0.0035 },
        30.0,
        config,
    );

    // Warm up with a few steps, then measure 100.
    for _ in 0..10 {
        solver.step();
    }
    let t = Instant::now();
    let mut sink = 0.0f64;
    for _ in 0..100 {
        sink += solver.step();
    }
    let per_step = t.elapsed() / 100;
    let stats = solver.stats(110);

    println!("hemodynamics-cfd benchmark:");
    println!("  24x12x12 cylinder, 250 SOR sweeps/step");
    println!(
        "  {:<28} {:>10.2?} per step (sink {sink:.4})",
        "projection step", per_step
    );
    println!(
        "  inlet flow {:.1} mm3/s, outlet {:.1} mm3/s",
        stats.inlet_flow, stats.outlet_flow
    );
    println!(
        "  ~{:.1} steps/s single-threaded",
        1.0 / per_step.as_secs_f64()
    );
}
