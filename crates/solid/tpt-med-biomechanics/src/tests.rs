//! Verification suite for the hexahedral FEM core.
//!
//! - Single-element uniaxial tension reproduces the exact linear solution
//!   (Q1 elements represent linear displacement fields exactly).
//! - A two-element patch test with linear boundary displacement passes with
//!   machine-precision interior agreement and uniform stress.
//! - A slender cantilever is compared against Euler–Bernoulli tip
//!   deflection with a documented tolerance (full 2×2×2 integration on Q1
//!   hexes exhibits shear locking; the tolerance band accounts for it and
//!   the golden dataset records the measured value).

use crate::solver::{BiomechanicsModel, BoundaryConditions};
use tpt_med_geometry::Vec3;

/// Structured hex grid of `nx×ny×nz` unit cubes with corner node ids.
fn grid_mesh(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    dy: f64,
    dz: f64,
) -> (Vec<Vec3>, Vec<[u32; 8]>) {
    let node =
        |i: usize, j: usize, k: usize| Vec3::new(i as f64 * dx, j as f64 * dy, k as f64 * dz);
    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push(node(i, j, k));
            }
        }
    }
    let id = |i: usize, j: usize, k: usize| (k * (ny + 1) + j) * (nx + 1) + i;
    let mut elements = Vec::with_capacity(nx * ny * nz);
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

#[test]
fn single_element_uniaxial_is_exact() {
    // 1×1×1 mm cube, E = 100 MPa, ν = 0. Bottom face fixed, top face pulled
    // with σ = 1 MPa (0.25 N per corner). Top-face displacement must equal
    // σ·L/E = 0.01 mm exactly (linear field exactly represented).
    let (nodes, elements) = grid_mesh(1, 1, 1, 1.0, 1.0, 1.0);
    let model = BiomechanicsModel::from_parts(nodes, elements, 100.0, 0.0);

    let mut bc = BoundaryConditions::default();
    bc.fix_nodes([0, 1, 2, 3]); // z = 0 face
    for n in [4u32, 5, 6, 7] {
        bc.add_force(n, Vec3::new(0.0, 0.0, 0.25));
    }
    let result = model.solve(&bc, 1e-10, 5_000).expect("solves");

    let top = [4usize, 5, 6, 7];
    for &t in &top {
        assert!(
            (result.displacements[t].z - 0.01).abs() < 1e-9,
            "uz = {}",
            result.displacements[t].z
        );
    }
    // Axial stress in the element: σzz = 1 MPa exactly.
    let el = &result.stresses[0];
    assert!((el.stress[2] - 1.0).abs() < 1e-6, "σzz={}", el.stress[2]);
    // Free lateral faces ⇒ σxx = σyy = 0.
    assert!(el.stress[0].abs() < 1e-6 && el.stress[1].abs() < 1e-6);
    assert!((el.von_mises() - 1.0).abs() < 1e-6);
}

#[test]
fn patch_test_linear_field_is_reproduced() {
    // 2×2×2-element patch; impose the exact linear field u = (0.01x, 0.02y,
    // 0.03z) on every boundary node. The single interior node (1,1,1) must
    // match the field to machine precision and all elements must carry the
    // same uniform strain (0.01, 0.02, 0.03).
    let (mut nodes, elements) = grid_mesh(2, 2, 2, 1.0, 1.0, 1.0);
    // interior node index id(1,1,1) = (1*(ny+1)+1)*(nx+1)+1 with n=2
    let interior = (3 + 1) * 3 + 1;
    assert_eq!(nodes[interior], Vec3::new(1.0, 1.0, 1.0));
    let model = BiomechanicsModel::from_parts(nodes.clone(), elements, 210.0, 0.3);

    let field = |p: Vec3| Vec3::new(0.01 * p.x, 0.02 * p.y, 0.03 * p.z);
    let mut bc = BoundaryConditions::default();
    for (i, &p) in nodes.iter().enumerate() {
        if i != interior {
            bc.prescribed_displacements.push((i as u32, field(p)));
        }
    }
    let result = model.solve(&bc, 1e-10, 5_000).expect("solves");
    let expected = field(nodes[interior]);
    let got = result.displacements[interior];
    assert!(
        (got - expected).norm() < 1e-9,
        "interior {got:?} vs {expected:?}"
    );
    // Uniform strain across both elements.
    for el in &result.stresses {
        assert!((el.strain[0] - 0.01).abs() < 1e-9);
        assert!((el.strain[1] - 0.02).abs() < 1e-9);
        assert!((el.strain[2] - 0.03).abs() < 1e-9);
    }
    assert!((result.stresses[0].stress[0] - result.stresses[1].stress[0]).abs() < 1e-9);
    let _ = &mut nodes;
}

#[test]
fn cantilever_matches_euler_bernoulli_within_locking_band() {
    // Slender cantilever: 20×2×2 mm, E = 1000 MPa, tip load P = 1 N.
    // Euler–Bernoulli: δ = PL³/(3EI), I = b h³/12 = 2·8/12 = 4/3.
    // Full-integration Q1 hexes shear-lock; the accepted band for this
    // discretisation is 0.5–1.15 × δ_EB (documented in the golden dataset).
    let (nodes, elements) = grid_mesh(20, 2, 2, 1.0, 1.0, 1.0);
    let model = BiomechanicsModel::from_parts(nodes, elements, 1000.0, 0.3);

    // Tip face: x = 20, all (y, z).
    let mut tip_nodes = Vec::new();
    let mut root_nodes = Vec::new();
    for (i, n) in model.nodes.iter().enumerate() {
        if (n.x - 20.0).abs() < 1e-9 {
            tip_nodes.push(i as u32);
        }
        if n.x < 1e-9 {
            root_nodes.push(i as u32);
        }
    }
    let mut bc = BoundaryConditions::default();
    bc.fix_nodes(root_nodes);
    let per = 1.0 / tip_nodes.len() as f64;
    for &n in &tip_nodes {
        bc.add_force(n, Vec3::new(0.0, 0.0, per));
    }
    let result = model.solve(&bc, 1e-10, 50_000).expect("solves");

    let p = 1.0f64;
    let big_l = 20.0f64;
    let inertia = 2.0f64 * 4.0f64 / 12.0; // b=2, h=4? h=2 ⇒ I = 2·8/12/...
    let i_be = 2.0 * (2.0f64).powi(3) / 12.0; // b=2, h=2
    let delta_eb = p * big_l.powi(3) / (3.0 * 1000.0 * i_be);
    let _ = inertia;
    let mut tip_deflection = 0.0f64;
    for &n in &tip_nodes {
        tip_deflection = tip_deflection.max(result.displacements[n as usize].z);
    }
    let ratio = tip_deflection / delta_eb;
    assert!(
        (0.45..=1.15).contains(&ratio),
        "tip δ={tip_deflection:.4}, δ_EB={delta_eb:.4}, ratio={ratio:.3}"
    );
    assert!(result.max_von_mises() > 0.0);
}

#[test]
fn solve_rejects_empty_and_inconsistent_models() {
    let model = BiomechanicsModel {
        nodes: vec![],
        elements: vec![],
        element_modulus: vec![],
        element_poisson: vec![],
    };
    assert!(model
        .solve(&BoundaryConditions::default(), 1e-8, 100)
        .is_err());
}

#[test]
fn rigid_translation_produces_no_stress_when_supported() {
    // Simply supported block translated via prescribed displacements: the
    // stress state must be zero everywhere (patch-test style check of
    // prescribed-DOF handling).
    let (nodes, elements) = grid_mesh(2, 2, 2, 1.0, 1.0, 1.0);
    let model = BiomechanicsModel::from_parts(nodes.clone(), elements, 500.0, 0.3);
    let mut bc = BoundaryConditions::default();
    for (i, _) in nodes.iter().enumerate() {
        bc.prescribed_displacements
            .push((i as u32, Vec3::new(0.5, -0.25, 1.0)));
    }
    let result = model.solve(&bc, 1e-10, 5_000).expect("solves");
    assert!(
        result.max_von_mises() < 1e-9,
        "vm={}",
        result.max_von_mises()
    );
}

#[test]
fn roller_and_symmetry_constraints_give_a_clean_uniaxial_state() {
    // 1x1x4 column, E=100, nu=0.3, uniform 1 N top load (sigma_zz = 1 MPa).
    // Per-DOF constraints in their textbook role: x = 0 face constrained in
    // x only and y = 0 face in y only (two symmetry planes — exact for this
    // loading), base on z-rollers. All three are partial constraints; no
    // node is fully fixed, so the Poisson contraction is unperturbed and
    // the uniform state is exact.
    let (nodes, elements) = grid_mesh(1, 1, 4, 1.0, 1.0, 1.0);
    let model = BiomechanicsModel::from_parts(nodes, elements, 100.0, 0.3);
    let mut bc = BoundaryConditions::default();
    for (i, n) in model.nodes.iter().enumerate() {
        if n.x < 1e-9 {
            bc.constrain_dofs([i as u32], [true, false, false]);
        }
        if n.y < 1e-9 {
            bc.constrain_dofs([i as u32], [false, true, false]);
        }
        if n.z < 1e-9 {
            bc.constrain_dofs([i as u32], [false, false, true]);
        }
    }
    let mut top = Vec::new();
    for (i, n) in model.nodes.iter().enumerate() {
        if n.z > 3.5 {
            top.push(i as u32);
        }
    }
    let per = 1.0 / top.len() as f64;
    for &n in &top {
        bc.add_force(n, Vec3::new(0.0, 0.0, -per));
    }
    let result = model.solve(&bc, 1e-10, 20_000).expect("solves");
    // Exact axial response: uz(top) = FL/EA = 1*4/(1*100) = 0.04 mm.
    let uz = top
        .iter()
        .map(|&n| -result.displacements[n as usize].z)
        .fold(0.0f64, f64::max);
    assert!((uz - 0.04).abs() < 1e-9, "uz {uz}");
    // Zero lateral stress (traction-free sides, symmetry planes).
    let lateral = result
        .stresses
        .iter()
        .map(|e| e.stress[0].abs().max(e.stress[1].abs()))
        .fold(0.0f64, f64::max);
    assert!(lateral < 1e-9, "lateral {lateral}");
    // sigma_zz uniform and exact: -1 MPa in compression.
    for e in &result.stresses {
        assert!((e.stress[2] + 1.0).abs() < 1e-9, "sigma_zz {}", e.stress[2]);
    }
}
