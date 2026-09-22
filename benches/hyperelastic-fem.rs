//! Hyperelastic material + hex-FEM benchmark.
//!
//! Run: `cargo bench -p tpt-med-benches --bench hyperelastic-fem`

use std::time::Instant;

use tpt_med_geometry::{Mat3, Vec3};
use tpt_med_tissue::models::{NeoHookeanParams, TissueModel};
use tpt_med_units::Force;

fn grid_mesh(nx: usize, ny: usize, nz: usize) -> (Vec<Vec3>, Vec<[u32; 8]>) {
    let id = |i: usize, j: usize, k: usize| (k * (ny + 1) + j) * (nx + 1) + i;
    let mut nodes = Vec::new();
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push(Vec3::new(i as f64, j as f64, k as f64));
            }
        }
    }
    let mut elements = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                elements.push([
                    id(i, j, k) as u32,
                    id(i + 1, j, k) as u32,
                    id(i + 1, j + 1, k) as u32,
                    id(i, j + 1, k) as u32,
                    id(i, j, k + 1) as u32,
                    id(i + 1, j, k + 1) as u32,
                    id(i + 1, j + 1, k + 1) as u32,
                    id(i, j + 1, k + 1) as u32,
                ]);
            }
        }
    }
    (nodes, elements)
}

fn main() {
    let mut timings = Vec::new();

    // Constitutive stress evaluation (central-difference heavy path).
    let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 10.0 });
    let lam = 1.4f64;
    let f = Mat3::from_rows(
        Vec3::new(lam, 0.0, 0.0),
        Vec3::new(0.0, 1.0 / lam.sqrt(), 0.0),
        Vec3::new(0.0, 0.0, 1.0 / lam.sqrt()),
    );
    let t0 = Instant::now();
    let mut sink = 0.0f64;
    for _ in 0..2000 {
        sink += model.first_piola_numerical(&f).at(0, 0);
    }
    timings.push(("2000x FD first_piola", t0.elapsed(), sink));

    // FEM: 16x4x4 cantilever under tip load.
    let (nodes, elements) = grid_mesh(16, 4, 4);
    let model = tpt_med_biomechanics::BiomechanicsModel::from_parts(nodes, elements, 1000.0, 0.3);
    let mut bc = tpt_med_biomechanics::BoundaryConditions::default();
    let mut tip = Vec::new();
    for (i, n) in model.nodes.iter().enumerate() {
        if n.x < 1e-9 {
            bc.fix_nodes([i as u32]);
        }
        if n.x > 15.5 {
            tip.push(i as u32);
        }
    }
    let per = Force::from_n(1.0).to_n() / tip.len() as f64;
    for &n in &tip {
        bc.add_force(n, Vec3::new(0.0, 0.0, -per));
    }
    let t1 = Instant::now();
    let result = model.solve(&bc, 1e-10, 50_000).expect("solves");
    let fem_time = t1.elapsed();
    timings.push(("hex FEM 16x4x4 solve", fem_time, result.max_von_mises()));

    println!("hyperelastic-fem benchmark:");
    for (name, d, sink) in &timings {
        println!("  {name:<28} {:>10.2?}  (sink {sink:.4})", d);
    }
}
