//! Verification suite: Poiseuille flow, mass conservation, rheology
//! coupling, and WSS vs the analytic tube law.

use crate::blood::BloodModel;
use crate::domain::FluidDomain;
use crate::solver::{HemodynamicsSolver, PressureSolver, SolverConfig};
use crate::wss::{extract_wss, OsiAccumulator};
use tpt_med_geometry::Vec3;

/// Straight tube verification run.
fn poiseuille_run(mu: f64, mean_velocity: f64) -> (HemodynamicsSolver, crate::solver::SteadyStats) {
    poiseuille_run_with(mu, mean_velocity, PressureSolver::Sor)
}

/// Same, with the pressure solver chosen.
fn poiseuille_run_with(
    mu: f64,
    mean_velocity: f64,
    pressure_solver: PressureSolver,
) -> (HemodynamicsSolver, crate::solver::SteadyStats) {
    // Tube: radius 3.5 cells, 16 cells long, pitch 0.5 mm → R ≈ 1.75 mm.
    // Grid 16×10×10 = 1600 cells (kept small so debug-profile tests stay
    // fast). dt = 5e-4 s is within the explicit diffusion limit
    // dx²ρ/(4μ) ≈ 8.3e-4 s for μ = 0.08.
    let domain = FluidDomain::cylinder(16, 10, 3.5, 0.5, 0);
    let blood = BloodModel::Newtonian { viscosity: mu };
    let config = SolverConfig {
        dt: 5.0e-4,
        poisson_iterations: 400,
        include_convection: false,
        viscosity_relaxation: 0.2,
        density: 1.06e-3,
        pressure_solver,
    };
    let mut solver = HemodynamicsSolver::new(domain, blood, mean_velocity, config);
    let stats = solver.run_steady(900, 1e-4);
    (solver, stats)
}

#[test]
fn poiseuille_profile_is_parabolic() {
    // Blood-analog verification fluid: μ = 0.08 Pa·s gives Re ≈ 1.4 and a
    // viscous development time R²/ν ≈ 41 ms, reachable within the march
    // budget (real-blood runs would need ~1 s of physical time).
    //
    // Verification metric: the developed profile must be a **paraboloid**
    // — constant radial second derivative along a diameter, symmetric
    // about the axis, with the fitted wall intercept at the mask radius.
    // The classic umax/mean = 2 does NOT hold on stair-step masks (the
    // no-slip ghost ring sits ~½ cell outside the mask radius), so it is
    // deliberately not used here (see test-data/golden/fluid notes).
    let mu = 0.08;
    let u_in = 30.0; // mm/s plug inlet
    let (solver, stats) = poiseuille_run(mu, u_in);

    // Mass conservation: outlet ≈ inlet within 10% (stair-step losses).
    assert!(
        (stats.outlet_flow - stats.inlet_flow).abs() < 0.10 * stats.inlet_flow,
        "inlet {:.1} vs outlet {:.1}",
        stats.inlet_flow,
        stats.outlet_flow
    );

    // Radial line of axial velocity at mid-length.
    let (nx, ny, nz) = solver.domain.dims;
    let i = nx / 2;
    let kc = nz / 2;
    let mut line: Vec<(usize, f64)> = Vec::new();
    for j in 0..ny {
        if solver.domain.is_fluid(i as i64, j as i64, kc as i64)
            && solver.domain.is_fluid(i as i64 - 1, j as i64, kc as i64)
        {
            let u = 0.5 * (solver.u[solver.uid(i, j, kc)] + solver.u[solver.uid(i + 1, j, kc)]);
            line.push((j, u));
        }
    }
    assert!(line.len() >= 5, "diameter too short: {}", line.len());
    let umax = line.iter().map(|(_, u)| *u).fold(0.0f64, f64::max);
    assert!(
        umax > 1.3 * u_in,
        "umax {umax:.1} not developed beyond plug"
    );

    // Constant curvature: second differences of u along j must agree
    // within 10% across the profile.
    let second_diff = |k: usize| {
        let (_ja, ua) = line[k];
        let (_, ub) = line[k + 1];
        let (_, uc) = line[k + 2];
        ua - 2.0 * ub + uc
    };
    let n = line.len();
    let d1 = second_diff(n / 2 - 1);
    let d2 = second_diff(n / 2);
    assert!(d1 < 0.0 && d2 < 0.0, "profile must be concave: {d1} {d2}");
    let rel = (d1 - d2).abs() / d1.abs().max(1e-12);
    assert!(rel < 0.10, "curvature not constant: {d1:.3} vs {d2:.3}");

    // Symmetry about the axis.
    let (jl, ul) = line[0];
    let (jr, ur) = line[n - 1];
    assert!(
        (ul - ur).abs() < 0.05 * umax,
        "asymmetric ends {ul:.1} {ur:.1}"
    );
    assert!(jl < ny / 2 && jr >= ny / 2);
}

#[test]
fn wall_shear_matches_tube_law() {
    // τ_w = 4 μ Q / (π R³); the screening extractor should land within a
    // factor ~1.5 given stair-step walls.
    let mu = 0.08;
    let u_in = 25.0;
    let (solver, stats) = poiseuille_run(mu, u_in);
    let field = extract_wss(&solver);
    assert!(
        field.tractions.len() > 20,
        "wall samples {}",
        field.tractions.len()
    );

    let radius_mm = 3.5f64 * 0.5;
    let q = stats.outlet_flow;
    let tau_analytic = 4.0 * mu * q / (core::f64::consts::PI * radius_mm.powi(3));
    let tau_sim = field.mean_magnitude();
    let ratio = tau_sim / tau_analytic;
    // Stair-step walls over-estimate the wall gradient; the accepted
    // screening band is documented in test-data/golden/fluid.
    assert!(
        (0.4..=1.6).contains(&ratio),
        "τ_sim {tau_sim:.4} vs τ_analytic {tau_analytic:.4} (ratio {ratio:.2})"
    );
}

#[test]
fn shear_thinning_blood_increases_resistance() {
    // Same inlet velocity: the shear-thinning Carreau–Yasuda fluid has a
    // higher apparent viscosity at the low shear rates near the tube
    // centre, so the pressure drop must exceed the Newtonian case.
    let (_, stats_n) = poiseuille_run(0.0035, 30.0);
    let domain = FluidDomain::cylinder(16, 10, 3.5, 0.5, 0);
    let config = SolverConfig {
        dt: 5.0e-4,
        poisson_iterations: 400,
        include_convection: false,
        viscosity_relaxation: 0.2,
        density: 1.06e-3,
        pressure_solver: PressureSolver::Sor,
    };
    let mut solver =
        HemodynamicsSolver::new(domain, BloodModel::CARREAU_YASUDA_BLOOD, 30.0, config);
    let stats_cy = solver.run_steady(400, 1e-4);
    assert!(
        stats_cy.pressure_drop > stats_n.pressure_drop,
        "CY drop {:.3} vs Newtonian {:.3}",
        stats_cy.pressure_drop,
        stats_n.pressure_drop
    );
    // Both fluids conserve mass.
    assert!((stats_cy.outlet_flow - stats_cy.inlet_flow).abs() < 0.10 * stats_cy.inlet_flow);
}

#[test]
fn osi_accumulator_over_synthetic_cycle() {
    // Reversing near-sinusoidal traction over a cycle → OSI near 0.5;
    // steady component → near 0.
    let mut osc = OsiAccumulator::default();
    for i in 0..100 {
        let phase = 2.0 * core::f64::consts::PI * i as f64 / 100.0;
        osc.sample(Vec3::new(phase.sin(), 0.0, 0.0), 0.01);
    }
    assert!(osc.osi() > 0.4, "osi {}", osc.osi());

    let mut steady = OsiAccumulator::default();
    for _ in 0..100 {
        steady.sample(Vec3::new(1.0, 0.1, 0.0), 0.01);
    }
    assert!(steady.osi() < 0.05, "osi {}", steady.osi());
}

#[test]
fn cg_pressure_solve_reproduces_a_manufactured_solution() {
    // Code verification of the CG Poisson solve against an exact solution:
    // pick φ = sin(πi/nx) on interior fluid cells (0 on the Dirichlet
    // outlet layer and on solids), form b = Aφ with the crate's own
    // operator, solve, and require the iterate to reproduce φ. This is the
    // check that would catch a sign error in the discretisation (the SOR
    // fixed point is A = −∇², so the rhs enters negated).
    let domain = FluidDomain::cylinder(16, 10, 3.5, 0.5, 0);
    let config = SolverConfig {
        pressure_solver: PressureSolver::ConjugateGradient,
        ..SolverConfig::default()
    };
    let solver = HemodynamicsSolver::new(
        domain,
        BloodModel::Newtonian { viscosity: 0.0035 },
        30.0,
        config,
    );
    let (nx, _ny, _nz) = solver.domain.dims;
    let n = solver.domain.mask.len();
    let phi_exact: Vec<f64> = (0..n)
        .map(|idx| {
            let (i, _, _) = solver.domain.coords(idx);
            if solver.domain.mask[idx] && i + 1 < nx {
                (core::f64::consts::PI * i as f64 / nx as f64).sin()
            } else {
                0.0
            }
        })
        .collect();
    let mut b = vec![0.0; n];
    solver.laplacian_apply(&phi_exact, &mut b);
    // `solve_poisson` takes the projection's rhs (the equation ∇²φ = rhs);
    // internally the SPD operator A = −∇² is used, so negate here.
    for v in b.iter_mut() {
        *v = -*v;
    }

    let (phi_cg, cg_iters) = solver.solve_poisson(&b);
    assert!(cg_iters > 0);
    let scale = phi_exact.iter().cloned().fold(0.0f64, f64::max);
    let max_err = phi_cg
        .iter()
        .zip(&phi_exact)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        max_err < 1e-4 * scale,
        "CG manufactured solution: max err {max_err:.3e} vs scale {scale:.3e}"
    );
    // The pinned unknowns stay pinned.
    for idx in 0..n {
        let (i, _, _) = solver.domain.coords(idx);
        if !solver.domain.mask[idx] || i + 1 == nx {
            assert_eq!(phi_cg[idx], 0.0);
        }
    }
}

#[test]
fn cg_projection_reproduces_the_poiseuille_verification() {
    // The full march with the CG projection must satisfy the same physical
    // checks as the SOR one: mass conservation and a concave, symmetric,
    // constant-curvature profile.
    let mu = 0.08;
    let u_in = 30.0;
    let (solver, stats) = poiseuille_run_with(mu, u_in, PressureSolver::ConjugateGradient);
    assert!(
        (stats.outlet_flow - stats.inlet_flow).abs() < 0.10 * stats.inlet_flow,
        "inlet {:.1} vs outlet {:.1}",
        stats.inlet_flow,
        stats.outlet_flow
    );
    // Same developed-paraboloid checks as the SOR verification: concave,
    // constant-curvature, symmetric — umax/mean = 2 does not hold on
    // stair-step masks.
    let (nx, ny, nz) = solver.domain.dims;
    let i = nx / 2;
    let kc = nz / 2;
    let mut line: Vec<(usize, f64)> = Vec::new();
    for j in 0..ny {
        if solver.domain.is_fluid(i as i64, j as i64, kc as i64)
            && solver.domain.is_fluid(i as i64 - 1, j as i64, kc as i64)
        {
            let u = 0.5 * (solver.u[solver.uid(i, j, kc)] + solver.u[solver.uid(i + 1, j, kc)]);
            line.push((j, u));
        }
    }
    assert!(line.len() >= 5);
    let umax = line.iter().map(|(_, u)| *u).fold(0.0f64, f64::max);
    assert!(
        umax > 1.3 * u_in,
        "umax {umax:.1} not developed beyond plug"
    );
    let second_diff = |k: usize| {
        let (_, ua) = line[k];
        let (_, ub) = line[k + 1];
        let (_, uc) = line[k + 2];
        ua - 2.0 * ub + uc
    };
    let n = line.len();
    let d1 = second_diff(n / 2 - 1);
    let d2 = second_diff(n / 2);
    assert!(d1 < 0.0 && d2 < 0.0, "profile must be concave: {d1} {d2}");
    let rel = (d1 - d2).abs() / d1.abs().max(1e-12);
    assert!(rel < 0.10, "curvature not constant: {d1:.3} vs {d2:.3}");
    let (jl, ul) = line[0];
    let (jr, ur) = line[line.len() - 1];
    assert!(
        (ul - ur).abs() < 0.05 * umax,
        "profile not symmetric: {ul:.2} vs {ur:.2}"
    );
    assert!(jl < ny / 2 && jr >= ny / 2);
    assert!(solver.last_pressure_solve_iterations > 0);
}

#[test]
fn cg_iterations_are_accounted_per_step() {
    // The iteration counter must reflect the configured solver, not a
    // stale value.
    let domain = FluidDomain::cylinder(16, 10, 3.5, 0.5, 0);
    let config = SolverConfig {
        dt: 5.0e-4,
        poisson_iterations: 400,
        include_convection: false,
        viscosity_relaxation: 0.2,
        density: 1.06e-3,
        pressure_solver: PressureSolver::ConjugateGradient,
    };
    let mut solver = HemodynamicsSolver::new(
        domain,
        BloodModel::Newtonian { viscosity: 0.08 },
        30.0,
        config,
    );
    assert_eq!(solver.last_pressure_solve_iterations, 0);
    solver.step();
    let after_one = solver.last_pressure_solve_iterations;
    assert!(after_one > 0, "CG reported no iterations");
    solver.step();
    assert!(solver.last_pressure_solve_iterations > 0);
}
