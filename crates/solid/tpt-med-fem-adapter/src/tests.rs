//! Verification suite (ASME V&V 40 code verification).
//!
//! The obligations RFC 0009 recorded as a placeholder are discharged here:
//!
//! 1. **Code verification** against the closed form the in-house core already
//!    verifies against — nominal incompressible Neo-Hookean uniaxial stress
//!    `mu (lambda - lambda^-2)`, the same target `tpt-med-tissue`'s
//!    `substrate-cross-check` compares to `tpt-fem-hyperelastic` — and against
//!    the in-house analytic volumetric branch under pure dilatation.
//! 2. A **patch-test identity** that is exact for a uniform deformation
//!    (`f[k, I] = V sum_L G[I, L] P[k, L]`), independent of any closed form.
//! 3. **Tangent verification** — the cheap `B^T A B` tangent against the
//!    residual-level finite-difference tangent, and minor symmetry.
//! 4. **Calculation verification** by mesh refinement, the check RFC 0009 said
//!    becomes possible only once a real 3-D assembly exists.
//! 5. **Contact**: the active set follows the moving geometry, a pressed body
//!    does not pass through the obstacle, and the contact reaction balances
//!    the punch force.

use crate::assembly::FnModel;
use crate::assembly::{
    coo_max_abs_diff, internal_force, material_tangent, tangent_stiffness,
    tangent_stiffness_numerical, AssemblyOptions, Constitutive,
};
use crate::contact::{ContactError, ContactPairing};
use crate::friction::{friction_terms, FrictionConfig, FrictionError};
use crate::loadpath::{solve_load_path, LoadPathError, LoadPathOptions};
use crate::mesh::{hex_box, hex_box_of, tet_box_of, MeshError};
use crate::solver::{solve_static, ContactConfig, SolveError, SolveOptions};
use tpt_fem_element::{Hex20, Hex27, Hex8, ReferenceElement, Tet10, Tet4};
use tpt_med_geometry::{Mat3, Vec3};
use tpt_med_tissue::{NeoHookeanParams, TissueModel};

/// Shear modulus used by the uniaxial checks, and the `c10` that produces it.
const MU: f64 = 0.98;
const C10: f64 = 0.49;
/// Volumetric penalty parameter for the compression-sensitive fixtures.
///
/// The in-house convention is `W_vol = (J - 1)^2 / d1`, so a *small* `d1` is a
/// *stiff* penalty. Two failure modes bracket the usable range, both measured
/// rather than guessed:
///
/// - Too stiff (`d1 <~ 0.2`) and trilinear hexahedra with full integration
///   lock: the volumetric locking error dominates.
/// - Too soft (`d1 >~ 1`) and the material is genuinely compressible, so the
///   answer is legitimately not the incompressible closed form at all.
///
/// `d1 = 0.5` sits where the discretisation error falls monotonically with
/// refinement toward the closed form (measured ratios against
/// `mu (lambda - lambda^-2)` at `lambda = 1.3`: 1.233, 1.058, 1.019, 1.003 for
/// 1x1x1 .. 4x4x4).
const D1: f64 = 0.5;

/// The material used by every fixture below.
fn nh() -> TissueModel {
    TissueModel::NeoHookean(NeoHookeanParams { c10: C10, d1: D1 })
}

/// Closed-form nominal (first Piola) uniaxial stress for incompressible
/// Neo-Hookean: `P11 = mu (lambda - lambda^-2)`.
fn closed_form_nominal(lam: f64) -> f64 {
    MU * (lam - lam.powi(-2))
}

/// Solves a homogeneous uniaxial tension block and returns the nominal axial
/// stress (reaction on the prescribed top face / initial area).
///
/// Half-symmetry model: the `y = 0` face is fully fixed, `u_z = 0` on the
/// `z = 0` face (the symmetry plane of a uniaxial test) and `u_y` prescribed on
/// the `y = l` face, with the lateral (`x`) faces free so the block contracts
/// by its Poisson response.
///
/// The symmetry plane is what makes this a *uniaxial strain* test. Without it,
/// the top face's axial displacement is satisfied just as well by a shear
/// state, which for this material is the *lower*-energy minimiser — an earlier
/// version of this fixture measured that shear response, and would have passed
/// against the wrong target.
fn uniaxial_nominal(n: usize, lam: f64, opts: &SolveOptions) -> f64 {
    let l = 10.0;
    let mesh = hex_box(n, n, n, l, l, l).expect("box");
    let mut dirichlet: Vec<(usize, f64)> = Vec::new();
    for &node in &mesh.face_nodes(1, false) {
        for c in 0..3 {
            dirichlet.push((mesh.dof(node, c), 0.0));
        }
    }
    for &node in &mesh.face_nodes(2, false) {
        dirichlet.push((mesh.dof(node, 2), 0.0));
    }
    for &node in &mesh.face_nodes(1, true) {
        dirichlet.push((mesh.dof(node, 1), (lam - 1.0) * l));
    }
    let result = solve_static(
        &mesh,
        &nh(),
        &vec![0.0; mesh.dof_count()],
        &dirichlet,
        opts,
        None,
    )
    .expect("uniaxial solve converges");
    let internal =
        internal_force(&mesh, &nh(), &result.displacement, &opts.assembly).expect("assembly");
    let reaction: f64 = mesh
        .face_nodes(1, true)
        .iter()
        .map(|&node| internal[mesh.dof(node, 1)])
        .sum();
    reaction / (l * l)
}

#[test]
fn uniaxial_tension_matches_closed_form() {
    // Code verification: the 3-D assembly reproduces the same nominal stress the
    // in-house core and the substrate's 1-D bar solve both reproduce. The
    // tolerance covers the discretisation error plus the penalty formulation's
    // own deviation from exact incompressibility, and is largest at the
    // smallest stretch.
    let opts = SolveOptions::default();
    for lam in [1.1, 1.3, 1.6] {
        let nominal = uniaxial_nominal(4, lam, &opts);
        let expected = closed_form_nominal(lam);
        assert!(
            (nominal - expected).abs() / expected < 0.03,
            "lam={lam}: FEM nominal {nominal}, closed form {expected}"
        );
    }
}

#[test]
fn uniaxial_tension_converges_under_mesh_refinement() {
    // Calculation verification: the error against the closed form must fall
    // monotonically as the mesh is refined.
    let opts = SolveOptions::default();
    let lam = 1.3;
    let expected = closed_form_nominal(lam);
    let errors: Vec<f64> = [1usize, 2, 3, 4]
        .iter()
        .map(|&n| ((uniaxial_nominal(n, lam, &opts) - expected) / expected).abs())
        .collect();
    assert!(
        errors.windows(2).all(|w| w[1] < w[0]),
        "error must decrease monotonically with refinement, got {errors:?}"
    );
    assert!(errors[3] < 1.0e-2, "4x4x4 error {} too large", errors[3]);
}

#[test]
fn zero_deformation_gives_zero_internal_force() {
    // At F = I every Neo-Hookean branch of the in-house model vanishes, so the
    // assembled internal force must be exactly zero — a cheap check that the
    // Jacobian, the quadrature weights and the DOF bookkeeping all cancel.
    let mesh = hex_box(2, 2, 2, 4.0, 4.0, 4.0).expect("box");
    let u = vec![0.0; mesh.dof_count()];
    let f = internal_force(&mesh, &nh(), &u, &AssemblyOptions::default()).expect("assembly");
    assert!(f.iter().all(|v| v.abs() < 1.0e-12), "{f:?}");
}

#[test]
fn uniform_deformation_satisfies_patch_identity() {
    // Constant-stress patch test. For a homogeneous deformation
    // u = (F0 - I) X on a single (affine) hexahedron the discrete stress is
    // exactly P(F0) and the gradients are exactly constant, so the internal
    // force is exactly
    //
    //     f[k, I] = V * sum_L G[I, L] P(F0)[k, L]
    //
    // No closed form, no material symmetry and no reference solution: if the
    // Jacobian, the quadrature weights or the DOF bookkeeping are wrong, this
    // fails. It is the identity a hexahedron must reproduce to pass a patch
    // test at all.
    let mesh = hex_box(1, 1, 1, 2.0, 2.0, 2.0).expect("box");
    let f0 = Mat3::from_rows(
        Vec3::new(1.2, 0.1, 0.0),
        Vec3::new(0.0, 0.9, 0.05),
        Vec3::new(0.0, 0.0, 1.1),
    );
    let mut u = vec![0.0; mesh.dof_count()];
    for n in 0..mesh.node_count() {
        let x = mesh.nodes()[n];
        let ux = (f0 * x).to_array();
        for k in 0..3 {
            u[3 * n + k] = ux[k] - x.to_array()[k];
        }
    }
    let opts = AssemblyOptions::default();
    let volume = mesh.element_volume(0).expect("upright element");
    let p = nh().first_piola(&f0);
    let grad = mesh
        .physical_gradients(0, &[0.0, 0.0, 0.0])
        .expect("upright element");
    let f = internal_force(&mesh, &nh(), &u, &opts).expect("assembly");
    for local in 0..8 {
        let node = mesh.elements()[0][local];
        for k in 0..3 {
            let expected: f64 = volume * (0..3).map(|l| grad[local][l] * p.at(k, l)).sum::<f64>();
            let got = f[3 * node + k];
            assert!(
                (got - expected).abs() < 1.0e-9 * expected.abs().max(1.0),
                "node {node} comp {k}: {got} vs {expected}"
            );
        }
    }
}

#[test]
fn uniform_dilatation_matches_the_analytic_volumetric_branch() {
    // Second exact check, against the in-house analytic volumetric branch that
    // `tpt-med-tissue` itself verifies: under pure dilatation the isochoric part
    // vanishes identically and P11 collapses to `2 J (J - 1) / d1 * J^-1/3`. A
    // uniform dilatation is affine, so a single hexahedron represents it
    // exactly and the comparison is limited only by round-off.
    let mesh = hex_box(1, 1, 1, 2.0, 2.0, 2.0).expect("box");
    let volume = mesh.element_volume(0).expect("upright element");
    let grad = mesh
        .physical_gradients(0, &[0.0, 0.0, 0.0])
        .expect("upright element");
    for &j in &[1.02f64, 1.1, 0.95] {
        let scale = j.powf(1.0 / 3.0);
        let u: Vec<f64> = (0..mesh.node_count())
            .flat_map(|n| {
                let x = mesh.nodes()[n].to_array();
                (0..3).map(move |k| (scale - 1.0) * x[k])
            })
            .collect();
        let f = internal_force(&mesh, &nh(), &u, &AssemblyOptions::default()).expect("assembly");
        let expected_p = 2.0 * j * (j - 1.0) / D1 * j.powf(-1.0 / 3.0);
        for local in 0..8 {
            let node = mesh.elements()[0][local];
            for k in 0..3 {
                // The patch identity again, with P diagonal at `expected_p`.
                let expected = volume * expected_p * grad[local][k];
                let got = f[3 * node + k];
                assert!(
                    (got - expected).abs() < 1.0e-9 * expected.abs().max(1.0),
                    "J={j} node {node} comp {k}: {got} vs {expected}"
                );
            }
        }
    }
}

#[test]
fn material_tangent_is_step_size_independent() {
    // The assembly's material tangent is itself a finite difference of `P`;
    // this checks the step scaling and index order by comparing it to a much
    // finer, independent difference of the same quantity.
    let model = nh();
    let f = Mat3::from_rows(
        Vec3::new(1.3, 0.2, 0.0),
        Vec3::new(0.0, 0.85, 0.1),
        Vec3::new(0.0, 0.0, 0.9),
    );
    let coarse = material_tangent(&model, &f, 1.0e-6);
    let fine = material_tangent(&model, &f, 1.0e-8);
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let a = coarse[i][j][k][l];
                    let b = fine[i][j][k][l];
                    assert!(
                        (a - b).abs() < 1.0e-4 * b.abs().max(1.0),
                        "A[{i}][{j}][{k}][{l}]: {a} vs {b}"
                    );
                }
            }
        }
    }
}

#[test]
fn analytic_tangent_matches_numerical_tangent() {
    // The RFC's open "analytic vs. numerical tangent" question, answered with
    // numbers on a deformed configuration: B^T A B agrees with d f_int / du to
    // finite-difference accuracy.
    let mesh = hex_box(1, 1, 1, 3.0, 3.0, 3.0).expect("box");
    let mut u = vec![0.0; mesh.dof_count()];
    for n in 0..mesh.node_count() {
        let x = mesh.nodes()[n].to_array();
        u[3 * n] = 0.4 * x[0];
        u[3 * n + 1] = 0.3 * x[1];
        u[3 * n + 2] = 0.2 * x[2];
    }
    let opts = AssemblyOptions::default();
    let k_analytic = tangent_stiffness(&mesh, &nh(), &u, &opts).expect("tangent");
    let k_numerical = tangent_stiffness_numerical(&mesh, &nh(), &u, &opts).expect("tangent");
    let scale = k_analytic.vals.iter().fold(1e-12f64, |m, v| m.max(v.abs()));
    let diff = coo_max_abs_diff(&k_analytic, &k_numerical);
    assert!(
        diff < 1.0e-4 * scale,
        "tangent disagreement {diff} (scale {scale})"
    );
}

#[test]
fn tangent_is_minor_symmetric() {
    // A consistent hyperelastic tangent is the Hessian of a scalar energy and
    // is therefore symmetric; a sign or index slip in the assembly breaks that
    // immediately.
    let mesh = hex_box(2, 1, 1, 2.0, 2.0, 2.0).expect("box");
    let mut u = vec![0.0; mesh.dof_count()];
    for n in 0..mesh.node_count() {
        u[3 * n + 1] = 0.25 * mesh.nodes()[n].to_array()[1];
    }
    let opts = AssemblyOptions::default();
    let k = tangent_stiffness(&mesh, &nh(), &u, &opts).expect("tangent");
    let csr = k.to_csr();
    let n = mesh.dof_count();
    let get = |row: usize, col: usize| -> f64 {
        csr.row_ptrs[row..row + 1]
            .iter()
            .flat_map(|&start| start..csr.row_ptrs[row + 1])
            .filter(|c| csr.col_ind[*c] == col)
            .map(|c| csr.values[c])
            .sum()
    };
    let mut worst = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            worst = worst.max((get(i, j) - get(j, i)).abs());
        }
    }
    assert!(worst < 1.0e-9, "tangent asymmetry {worst}");
}

#[test]
fn fn_model_adapter_delegates() {
    // The closure adapter must be interchangeable with the tissue models.
    let model = crate::FnModel(|f: &Mat3| *f * 7.0);
    let f = Mat3::from_rows(
        Vec3::new(1.1, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );
    let p: Mat3 = model.first_piola(&f);
    assert!((p.at(0, 0) - 7.7).abs() < 1e-12);
    let wrapped: &dyn Constitutive = &model;
    assert!((wrapped.first_piola(&f).at(0, 0) - 7.7).abs() < 1e-12);
}

#[test]
fn an_inverted_element_is_reported_not_silently_assembled() {
    // The in-house models evaluate `J^-2/3`, which is NaN for `J <= 0`. An
    // inverted element must therefore be an explicit error, not a NaN that
    // surfaces later as an unrelated "singular matrix" from the linear solver.
    let mesh = hex_box(1, 1, 1, 1.0, 1.0, 1.0).expect("box");
    let mut u = vec![0.0; mesh.dof_count()];
    for &node in &mesh.face_nodes(1, true) {
        u[mesh.dof(node, 1)] = -2.0;
    }
    let err = internal_force(&mesh, &nh(), &u, &AssemblyOptions::default())
        .expect_err("inversion must be reported");
    assert!(
        matches!(err, MeshError::InvertedDeformation { .. }),
        "got {err:?}"
    );
}

/// A block pressed down onto a rigid wall at `y = 0`: `hex_box(1, 1, 2, ...)`,
/// top face clamped in `x` and `z`, a prescribed downward `y` displacement on
/// the top face, and the bottom face paired with a wall at its own rest
/// positions.
///
/// Displacement control, not load control: the punch force is then a
/// *reaction* the solver reports rather than a load it applies. A load on a
/// prescribed DOF is condensed out of the system and does nothing at all, and a
/// load-controlled block standing on a wall has no `y` constraint without the
/// very contact it is meant to be testing.
struct WallSetup {
    mesh: crate::Hex8Mesh,
    dirichlet: Vec<(usize, f64)>,
    load: Vec<f64>,
    pairing: ContactPairing,
    slave: Vec<usize>,
    /// The prescribed downward displacement applied to each top node.
    #[allow(dead_code)]
    punch_displacement: f64,
}

fn wall_setup(punch_displacement: f64) -> WallSetup {
    let l = 10.0;
    let mesh = hex_box(1, 1, 2, l, l, l).expect("box");
    let top = mesh.face_nodes(1, true);
    let slave = mesh.face_nodes(1, false);
    let mut dirichlet: Vec<(usize, f64)> = Vec::new();
    for &node in &top {
        dirichlet.push((mesh.dof(node, 0), 0.0));
        dirichlet.push((mesh.dof(node, 2), 0.0));
        dirichlet.push((mesh.dof(node, 1), punch_displacement / top.len() as f64));
    }
    let load = vec![0.0; mesh.dof_count()];
    let master: Vec<Vec3> = slave.iter().map(|&n| mesh.nodes()[n]).collect();
    let pairing = ContactPairing::new(1, slave.iter().copied(), master).expect("axis 1 is valid");
    WallSetup {
        mesh,
        dirichlet,
        load,
        pairing,
        slave,
        punch_displacement,
    }
}

impl WallSetup {
    /// The punch force implied by the prescribed displacement: the sum of the
    /// internal forces at the prescribed top DOFs.
    fn punch_force(&self, displacement: &[f64]) -> f64 {
        let internal = internal_force(
            &self.mesh,
            &nh(),
            displacement,
            &SolveOptions::default().assembly,
        )
        .expect("assembly");
        -self
            .mesh
            .face_nodes(1, true)
            .iter()
            .map(|&node| internal[self.mesh.dof(node, 1)])
            .sum::<f64>()
    }

    /// Lowest `y` reached by the slave (bottom) face.
    fn lowest_slave(&self, displacement: &[f64]) -> f64 {
        self.slave
            .iter()
            .map(|&n| self.mesh.nodes()[n].y + displacement[self.mesh.dof(n, 1)])
            .fold(f64::INFINITY, f64::min)
    }
}

#[test]
fn contact_pairing_rejects_bad_definitions() {
    assert_eq!(
        ContactPairing::new(3, [0], [Vec3::ZERO]),
        Err(ContactError::AxisOutOfRange(3))
    );
    assert_eq!(
        ContactPairing::new(1, [0], Vec::<Vec3>::new()),
        Err(ContactError::EmptyMaster)
    );
}

#[test]

fn contact_active_set_follows_the_moving_geometry() {
    // The specific thing RFC 0009 left unresolved: the active set must be a
    // function of the current geometry, not of the undeformed configuration.
    // Touching the obstacle counts as active (see
    // `ContactPairing::with_activation_tolerance`), so the three states below
    // are: resting on the wall, lifted clear of it, pressed into it.
    let s = wall_setup(-0.3);
    let at_rest = vec![0.0; s.mesh.dof_count()];
    let lifted = {
        let mut u = at_rest.clone();
        for &node in &s.slave {
            u[s.mesh.dof(node, 1)] = 0.1;
        }
        u
    };
    let pressed = {
        let mut u = at_rest.clone();
        for &node in &s.slave {
            u[s.mesh.dof(node, 1)] = -0.1;
        }
        u
    };

    assert_eq!(
        s.pairing
            .active_constraints(&s.mesh, &at_rest)
            .expect("ok")
            .len(),
        s.slave.len(),
        "a face resting on the obstacle must be constrained from the first iteration"
    );
    assert!(s
        .pairing
        .active_constraints(&s.mesh, &lifted)
        .expect("ok")
        .is_empty());
    assert_eq!(
        s.pairing
            .active_constraints(&s.mesh, &pressed)
            .expect("ok")
            .len(),
        s.slave.len()
    );
    assert_eq!(
        s.pairing.max_penetration(&s.mesh, &lifted).expect("ok"),
        0.0
    );
    assert!((s.pairing.max_penetration(&s.mesh, &pressed).expect("ok") - 0.1).abs() < 1e-12);
    for c in s.pairing.active_constraints(&s.mesh, &pressed).expect("ok") {
        assert_eq!(c.lower, 0.0);
    }
}

#[test]
fn contact_holds_the_body_out_of_the_obstacle() {
    // The block is punched down onto a rigid wall with a penalty far above the
    // structural stiffness: it must stop, with a penetration well below the
    // geometric scale, and the wall reaction must carry the punch force.
    let s = wall_setup(-0.3);
    let penalty = 1.0e4;
    let result = solve_static(
        &s.mesh,
        &nh(),
        &s.load,
        &s.dirichlet,
        &SolveOptions::default(),
        Some(ContactConfig {
            radial: None,
            pairing: &s.pairing,
            penalty,
            friction: None,
        }),
    )
    .expect("contact solve converges");

    let summary = result.contact.expect("contact summary present");
    assert_eq!(summary.active_constraints.len(), s.slave.len());
    assert!(
        summary.max_penetration < 1.0e-3,
        "penetrated by {}",
        summary.max_penetration
    );
    let punch = s.punch_force(&result.displacement);
    assert!(
        punch > 0.0,
        "the prescribed downward displacement must compress the block, got {punch}"
    );
    assert!(
        (summary.total_reaction - punch).abs() < 0.01 * punch,
        "wall reaction {} does not balance the punch force {}",
        summary.total_reaction,
        punch
    );
    assert!(
        result.residual_norm < 1.0e-8,
        "residual {}",
        result.residual_norm
    );
}

#[test]
fn contact_changes_the_answer_versus_no_contact() {
    // Without the obstacle the same prescribed displacement drives the bottom
    // face through the wall; with it, the answer is bounded. This is the
    // contrast that shows the contact terms reach the Newton residual at all.
    let s = wall_setup(-0.3);
    let free = solve_static(
        &s.mesh,
        &nh(),
        &s.load,
        &s.dirichlet,
        &SolveOptions::default(),
        None,
    )
    .expect("unconstrained solve converges");
    let min_free = s.lowest_slave(&free.displacement);
    assert!(
        min_free < -1.0e-2,
        "unconstrained body barely moved: {min_free}"
    );

    let constrained = solve_static(
        &s.mesh,
        &nh(),
        &s.load,
        &s.dirichlet,
        &SolveOptions::default(),
        Some(ContactConfig {
            radial: None,
            pairing: &s.pairing,
            penalty: 1.0e4,
            friction: None,
        }),
    )
    .expect("contact solve converges");
    let min_constrained = s.lowest_slave(&constrained.displacement);
    assert!(
        min_constrained > min_free,
        "contact did not hold the body: {min_constrained} vs {min_free}"
    );
    assert!(min_constrained > -1.0e-3, "penetrated to {min_constrained}");
}

#[test]
fn pulling_away_leaves_the_active_set_empty() {
    // The other half of "re-evaluated each iteration": a load that lifts the
    // body off the wall must not leave a stale constraint pushing on it.
    let s = wall_setup(0.3);
    let result = solve_static(
        &s.mesh,
        &nh(),
        &s.load,
        &s.dirichlet,
        &SolveOptions::default(),
        Some(ContactConfig {
            radial: None,
            pairing: &s.pairing,
            penalty: 1.0e4,
            friction: None,
        }),
    )
    .expect("solve converges");
    let summary = result.contact.expect("contact summary present");
    assert!(summary.active_constraints.is_empty());
    assert_eq!(summary.total_reaction, 0.0);
    assert_eq!(summary.max_penetration, 0.0);
}

#[test]
fn solve_rejects_a_mis_sized_load() {
    let mesh = hex_box(1, 1, 1, 1.0, 1.0, 1.0).expect("box");
    let err = solve_static(&mesh, &nh(), &[0.0; 5], &[], &SolveOptions::default(), None)
        .expect_err("length mismatch is an error, not a panic");
    assert!(matches!(
        err,
        SolveError::LoadSizeMismatch {
            expected: 24,
            found: 5
        }
    ));
}

#[test]
fn contact_jacobian_matches_the_finite_differenced_residual() {
    // The contact-augmented Jacobian the solver uses must be the derivative of
    // the contact-augmented residual, penalty term included. Getting this wrong
    // is silent: the solve still creeps toward a smaller residual, just never
    // quadratically, and the failure surfaces as a mysterious
    // `NotConverged` far from the cause.
    let s = wall_setup(-0.3);
    let penalty = 1.0e3;
    let opts = SolveOptions::default();
    let cfg = Some(ContactConfig {
        radial: None,
        pairing: &s.pairing,
        penalty,
        friction: None,
    });
    let n = s.mesh.dof_count();
    let mut u = vec![0.0; n];
    for (d, v) in &s.dirichlet {
        u[*d] = *v;
    }
    for &node in &s.slave {
        u[s.mesh.dof(node, 1)] = -0.01;
    }
    let active = s.pairing.active_constraints(&s.mesh, &u).expect("active");
    assert!(!active.is_empty(), "fixture must be in contact");
    let mut k = tangent_stiffness(&s.mesh, &nh(), &u, &opts.assembly).expect("tangent");
    let (k_c, _) = tpt_fem_contact::penalty_contact(
        &tpt_fem_sparse::Coo::new(),
        &vec![0.0; n],
        &active,
        penalty,
    );
    for i in 0..k_c.len() {
        k.push(k_c.rows[i], k_c.cols[i], k_c.vals[i]);
    }
    let get = |row: usize, col: usize| -> f64 {
        k.rows
            .iter()
            .zip(&k.cols)
            .zip(&k.vals)
            .filter(|((&r, &c), _)| r == row && c == col)
            .map(|(_, &v)| v)
            .sum()
    };
    let step = 1.0e-6;
    let mut probe = u.clone();
    for b in 0..n {
        probe.copy_from_slice(&u);
        probe[b] += step;
        let plus = crate::residual(&s.mesh, &nh(), &s.load, &probe, &opts, cfg).expect("r");
        probe[b] = u[b] - step;
        let minus = crate::residual(&s.mesh, &nh(), &s.load, &probe, &opts, cfg).expect("r");
        for a in 0..n {
            let fd = (plus[a] - minus[a]) / (2.0 * step);
            let an = get(a, b);
            assert!(
                (fd - an).abs() < 1.0e-4 * fd.abs().max(an.abs()).max(1.0),
                "J[{a}][{b}]: finite difference {fd}, assembled {an}"
            );
        }
    }
}

// --- Friction -------------------------------------------------------------
//
// Verification for the regularized Coulomb layer, in the same shape as the rest
// of this suite: the law against its closed form, the friction tangent against a
// differenced force, and a solve-level check that friction changes the
// converged answer in the physically right direction.

/// A `WallSetup` pressed into the wall with the slave face's tangential `x` DOFs
/// *released* and a tangential load applied to them.
///
/// The releasing matters: a prescribed tangential DOF is condensed out of every
/// linear solve, so a fixture that pins the shear by Dirichlet conditions
/// leaves friction nothing to resist and the answer is unchanged by definition.
/// The tangential drive has to be a load for the test to mean anything.
fn tangential_wall_setup(punch: f64, tangential_load: f64) -> WallSetup {
    let s = wall_setup(punch);
    let dirichlet = s
        .dirichlet
        .iter()
        .copied()
        .filter(|(dof, _)| !s.slave.iter().any(|&n| s.mesh.dof(n, 0) == *dof))
        .collect::<Vec<_>>();
    let mut load = vec![0.0; s.mesh.dof_count()];
    for &node in &s.slave {
        load[s.mesh.dof(node, 0)] = tangential_load / s.slave.len() as f64;
    }
    WallSetup {
        dirichlet,
        load,
        ..s
    }
}

/// A configuration with the slave face penetrating by `penetration` and
/// displaced tangentially by `slip`.
fn sheared_configuration(punch: f64, penetration: f64, slip: f64) -> (WallSetup, Vec<f64>) {
    let s = tangential_wall_setup(punch, 1.0);
    let n = s.mesh.dof_count();
    let mut u = vec![0.0; n];
    for (d, v) in &s.dirichlet {
        u[*d] = *v;
    }
    for &node in &s.slave {
        u[s.mesh.dof(node, 0)] = slip;
        u[s.mesh.dof(node, 1)] = -penetration;
    }
    (s, u)
}

#[test]
fn friction_config_rejects_unphysical_parameters() {
    // `matches!` rather than `assert_eq!`: a NaN never compares equal to itself,
    // so the derived `PartialEq` on `FrictionError` cannot assert the NaN case.
    assert!(matches!(
        FrictionConfig::new(-0.1, 1.0),
        Err(FrictionError::InvalidMu(m)) if m == -0.1
    ));
    assert!(matches!(
        FrictionConfig::new(f64::NAN, 1.0),
        Err(FrictionError::InvalidMu(_))
    ));
    assert!(matches!(
        FrictionConfig::new(0.3, -1.0),
        Err(FrictionError::InvalidTangentialStiffness(k)) if k == -1.0
    ));
    assert!(matches!(
        FrictionConfig::new(0.3, f64::INFINITY),
        Err(FrictionError::InvalidTangentialStiffness(_))
    ));
    assert!(FrictionConfig::new(0.0, 0.0).is_ok(), "zero is allowed");
}

#[test]
fn zero_friction_is_exactly_the_frictionless_case() {
    let (s, u) = sheared_configuration(-0.3, 0.01, 0.05);
    for cfg in [
        FrictionConfig::new(0.0, 1.0e4).expect("valid"),
        FrictionConfig::new(0.5, 0.0).expect("valid"),
    ] {
        let terms = friction_terms(&s.mesh, &s.pairing, &u, 1.0e4, cfg).expect("terms");
        assert_eq!(terms.slipping_nodes, 0);
        assert!(
            terms.force.iter().all(|&f| f == 0.0),
            "a zero coefficient or stiffness must give exactly zero force"
        );
        assert_eq!(terms.tangent.len(), 0, "and a zero tangent");
    }
}

#[test]
fn friction_force_obeys_the_coulomb_bound() {
    // Both regimes against the closed form `f = min(mu * f_n, k_t |s|) s/|s|`:
    // slipping saturates at `mu * f_n` and stops growing with slip; sticking is
    // the linear `k_t |s|`.
    let penalty = 1.0e4;
    let mu = 0.5;
    let kt = 1.0e3;
    let cfg = FrictionConfig::new(mu, kt).expect("valid");
    let penetration = 0.01;
    let bound = mu * penalty * penetration;

    for (slip, expect_slipping) in [(0.001, false), (100.0, true)] {
        let (s, u) = sheared_configuration(-0.3, penetration, slip);
        let terms = friction_terms(&s.mesh, &s.pairing, &u, penalty, cfg).expect("terms");
        assert_eq!(terms.slipping_nodes > 0, expect_slipping, "slip {slip}");
        assert!(!terms.states.is_empty(), "fixture must be in contact");

        for st in &terms.states {
            assert!(
                (st.bound - bound).abs() < 1.0e-6 * bound.max(1.0),
                "bound {} vs {bound}",
                st.bound
            );
            let magnitude = st.force[0].hypot(st.force[1]);
            let expected = if expect_slipping {
                bound
            } else {
                (kt * st.slip_norm).min(bound)
            };
            assert!(
                (magnitude - expected).abs() < 1.0e-6 * expected.max(1.0),
                "slip {}: |f| {magnitude} vs {expected}",
                st.slip_norm
            );
            // The force opposes the slip, never reinforces it.
            let opposing = st.force[0] * st.slip[0] + st.force[1] * st.slip[1];
            assert!(
                opposing <= 0.0,
                "friction must oppose slip, got dot {opposing}"
            );
        }
    }
}

#[test]
fn friction_is_not_applied_to_a_separating_node() {
    // The Coulomb bound scales with the normal reaction, so a node that has
    // lifted off must produce no friction at all — otherwise a body retracting
    // from a wall would be dragged back by a force it is not in contact with.
    let (s, mut u) = sheared_configuration(0.3, 0.01, 0.05);
    // Pull the whole slave face clear of the wall.
    for &node in &s.slave {
        u[s.mesh.dof(node, 1)] = 0.5;
    }
    let terms = friction_terms(
        &s.mesh,
        &s.pairing,
        &u,
        1.0e4,
        FrictionConfig::new(0.5, 1.0e3).expect("valid"),
    )
    .expect("terms");
    assert!(
        terms.states.is_empty(),
        "no active node should carry friction, got {}",
        terms.states.len()
    );
    assert!(terms.force.iter().all(|&f| f == 0.0));
}

#[test]
fn a_stuck_node_still_resists_a_tangential_perturbation() {
    // Zero slip means zero *force*, but the tangent must not be zero: the node is
    // held, and a linear solve that let it slide freely would creep. This is the
    // easiest thing to get wrong in a regularized formulation, and it is
    // invisible in the converged force — only in the displacement.
    let (s, mut u) = sheared_configuration(-0.3, 0.01, 0.05);
    for &node in &s.slave {
        u[s.mesh.dof(node, 0)] = 0.0;
    }
    let kt = 1.0e3;
    let terms = friction_terms(
        &s.mesh,
        &s.pairing,
        &u,
        1.0e4,
        FrictionConfig::new(0.5, kt).expect("valid"),
    )
    .expect("terms");
    assert!(!terms.states.is_empty(), "fixture must be in contact");
    for st in &terms.states {
        assert_eq!(st.slip_norm, 0.0);
        assert_eq!(st.force, [0.0, 0.0], "no slip, no force");
    }
    // Each stuck node contributes a +k_t diagonal at each tangential DOF.
    assert_eq!(
        terms.tangent.len(),
        2 * terms.states.len(),
        "expected a 2x2 diagonal per stuck node"
    );
    let sum: f64 = terms.tangent.vals.iter().sum();
    assert!(
        (sum - kt * 2.0 * terms.states.len() as f64).abs() < 1.0e-9 * kt,
        "tangent diagonal should be k_t per tangential DOF, got {sum}"
    );
}

#[test]
fn friction_tangent_matches_the_finite_differenced_force() {
    // The friction tangent is checked against a differenced *friction force*,
    // not against the full solver residual.
    //
    // Differencing the full residual would also be valid, but it cannot isolate
    // this layer: `contact_pairs` re-runs a nearest-point search at every
    // evaluation, so a tangential perturbation of a slave node can re-pair it to
    // a *different* master point, changing the normal gap and hence the Coulomb
    // bound discontinuously. That is a pre-existing property of the substrate's
    // pairing, not of the friction law, and it makes a full-residual difference
    // disagree at the 1e-2 level for reasons unrelated to the derivative under
    // test. Differencing the friction force alone holds the pairing fixed and
    // tests exactly what is claimed.
    //
    // Run on the *stick* branch, where the derivative is the unambiguous
    // `-k_t I`. The slip branch is not differenced here because the `min` kink
    // makes a central difference straddling the switch meaningless; its
    // magnitude is covered by the Coulomb bound test.
    let (s, u) = sheared_configuration(-0.3, 0.01, 0.001);
    let penalty = 1.0e4;
    let cfg = FrictionConfig::new(0.5, 1.0e3).expect("valid");
    let terms = friction_terms(&s.mesh, &s.pairing, &u, penalty, cfg).expect("terms");
    assert_eq!(
        terms.slipping_nodes, 0,
        "fixture must be sticking for this test to mean anything"
    );

    let n = s.mesh.dof_count();
    let get = |row: usize, col: usize| -> f64 {
        terms
            .tangent
            .rows
            .iter()
            .zip(&terms.tangent.cols)
            .zip(&terms.tangent.vals)
            .filter(|((&r, &c), _)| r == row && c == col)
            .map(|(_, &v)| v)
            .sum()
    };
    let step = 1.0e-7;
    let mut probe = u.clone();
    for b in 0..n {
        probe.copy_from_slice(&u);
        probe[b] += step;
        let plus = friction_terms(&s.mesh, &s.pairing, &probe, penalty, cfg).expect("terms");
        probe[b] = u[b] - step;
        let minus = friction_terms(&s.mesh, &s.pairing, &probe, penalty, cfg).expect("terms");
        for a in 0..n {
            let fd = (plus.force[a] - minus.force[a]) / (2.0 * step);
            let an = get(a, b);
            // Only entries the friction force actually reaches are interesting;
            // everywhere else both sides are exactly zero.
            if fd == 0.0 && an == 0.0 {
                continue;
            }
            assert!(
                (fd - an).abs() < 1.0e-4 * fd.abs().max(an.abs()).max(1.0),
                "df[{a}][{b}]: finite difference {fd}, analytic {an}"
            );
        }
    }
}

#[test]
fn friction_changes_the_converged_answer() {
    // End-to-end: a block pressed into the wall and driven tangentially by a
    // load must drift *less* with friction than without. This is the physical
    // statement the layer exists to make, checked on the converged solve rather
    // than the residual, so it also exercises the wiring in `solver`.
    let s = tangential_wall_setup(-0.3, 5.0);
    let penalty = 1.0e4;
    let opts = SolveOptions::default();
    let solve_with = |mu: Option<f64>| {
        solve_static(
            &s.mesh,
            &nh(),
            &s.load,
            &s.dirichlet,
            &opts,
            Some(ContactConfig {
                radial: None,
                pairing: &s.pairing,
                penalty,
                friction: mu.map(|m| FrictionConfig::new(m, 1.0e3).expect("valid")),
            }),
        )
        .expect("converges")
    };
    let free = solve_with(None);
    let held = solve_with(Some(0.5));

    // Mean tangential drift of the slave face.
    let drift = |r: &crate::SolveResult| -> f64 {
        s.slave
            .iter()
            .map(|&node| r.displacement[s.mesh.dof(node, 0)])
            .sum::<f64>()
            / s.slave.len() as f64
    };
    let (d_free, d_held) = (drift(&free), drift(&held));
    assert!(
        d_free > 0.0,
        "the tangential load must actually drive the face, got {d_free}"
    );
    assert!(
        d_held < d_free,
        "friction should hold the face back: {d_held} vs frictionless {d_free}"
    );
    // Friction is not a rigid clamp: it should reduce the drift without pinning
    // the face completely, which is what distinguishes a Coulomb law from a
    // stuck boundary condition.
    assert!(
        d_held > 0.0,
        "friction should not rigidly pin the face, got {d_held}"
    );
    let summary = held.contact.expect("contact ran");
    assert!(summary.slipping_nodes.is_some(), "friction was configured");
    // Friction acts tangentially, so the normal contact is unaffected by it.
    assert!(
        !summary.active_constraints.is_empty(),
        "the wall should still be holding the block"
    );
}

// --- Load path -------------------------------------------------------------
//
// Verification for proportional load stepping, on a genuinely load-controlled
// fixture, since that is the only case the driver exists for.

/// A block clamped at the bottom and pulled upward by a traction load, with
/// every DOF above the clamp free.
fn tensile_setup(traction: f64) -> (crate::Mesh<Hex8>, Vec<f64>, Vec<(usize, f64)>) {
    let l = 10.0;
    let mesh = hex_box(2, 2, 2, l, l, l).expect("box");
    let mut dirichlet = Vec::new();
    for n in mesh.face_nodes(1, false) {
        for c in 0..3 {
            dirichlet.push((mesh.dof(n, c), 0.0));
        }
    }
    let top = mesh.face_nodes(1, true);
    let area = l * l;
    let mut load = vec![0.0; mesh.dof_count()];
    for &n in &top {
        load[mesh.dof(n, 1)] = traction / area;
    }
    (mesh, load, dirichlet)
}

#[test]
fn load_path_ends_at_the_same_answer_as_a_single_solve() {
    // The whole point of stepping is that it must not change where the path
    // arrives. A converged equilibrium is a converged equilibrium regardless of
    // how it was reached, so the full-load point of a finely stepped path and a
    // single full-load solve must agree.
    let (mesh, load, dirichlet) = tensile_setup(0.05);
    let opts = SolveOptions::default();

    let one_shot =
        solve_static(&mesh, &nh(), &load, &dirichlet, &opts, None).expect("single solve converges");
    let path = solve_load_path(
        &mesh,
        &nh(),
        &load,
        &dirichlet,
        &opts,
        None,
        LoadPathOptions {
            steps: 8,
            max_cutbacks: 4,
        },
    )
    .expect("path converges");

    assert_eq!(path.steps.len(), 9, "8 increments plus the zero point");
    assert_eq!(path.reached(), 1.0, "the path must reach full load");
    let last = path.last().expect("a final point");
    let diff: f64 = last
        .displacement
        .iter()
        .zip(&one_shot.displacement)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        diff < 1.0e-6,
        "stepped and single-shot answers differ by {diff}"
    );
}

#[test]
fn load_path_starts_at_zero_and_advances_monotonically() {
    let (mesh, load, dirichlet) = tensile_setup(0.05);
    let path = solve_load_path(
        &mesh,
        &nh(),
        &load,
        &dirichlet,
        &SolveOptions::default(),
        None,
        LoadPathOptions::default(),
    )
    .expect("path converges");

    assert_eq!(path.steps[0].load_factor, 0.0);
    assert_eq!(path.steps[0].residual_norm, 0.0);
    assert!(
        path.steps[0].displacement.iter().all(|&d| d == 0.0),
        "the path must start undeformed"
    );
    for w in path.steps.windows(2) {
        assert!(
            w[1].load_factor > w[0].load_factor,
            "load factor must increase: {} then {}",
            w[0].load_factor,
            w[1].load_factor
        );
    }
    // Every returned point is a converged equilibrium, which is the property
    // that makes a point on the path safe to use.
    for s in &path.steps {
        assert!(
            s.residual_norm < 1.0e-6,
            "point at factor {} is not converged: {}",
            s.load_factor,
            s.residual_norm
        );
    }
}

#[test]
fn load_path_displacement_grows_with_load() {
    // A path that reached full load but moved *less* than a single solve would
    // mean the increments were not actually being applied.
    let (mesh, load, dirichlet) = tensile_setup(0.05);
    let opts = SolveOptions::default();
    let path = solve_load_path(
        &mesh,
        &nh(),
        &load,
        &dirichlet,
        &opts,
        None,
        LoadPathOptions::default(),
    )
    .expect("path converges");
    let top = mesh.face_nodes(1, true);
    let top_y = |u: &[f64]| -> f64 {
        top.iter().map(|&n| u[mesh.dof(n, 1)]).sum::<f64>() / top.len() as f64
    };
    for w in path.steps.windows(2) {
        assert!(
            top_y(&w[1].displacement) > top_y(&w[0].displacement),
            "the block should stretch further at each step"
        );
    }
    assert!(top_y(&path.last().expect("last").displacement) > 0.0);
}

#[test]
fn a_single_step_recovers_a_plain_static_solve() {
    // `steps: 1` is the degenerate path, and it must not be a special case: it
    // has to be the same solve `solve_static` performs.
    let (mesh, load, dirichlet) = tensile_setup(0.02);
    let opts = SolveOptions::default();
    let path = solve_load_path(
        &mesh,
        &nh(),
        &load,
        &dirichlet,
        &opts,
        None,
        LoadPathOptions {
            steps: 1,
            max_cutbacks: 0,
        },
    )
    .expect("converges");
    assert_eq!(path.steps.len(), 2, "the zero point plus one step");
    let direct = solve_static(&mesh, &nh(), &load, &dirichlet, &opts, None).expect("converges");
    let diff: f64 = path.steps[1]
        .displacement
        .iter()
        .zip(&direct.displacement)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        diff < 1.0e-9,
        "steps: 1 diverged from solve_static by {diff}"
    );
}

#[test]
fn load_path_rejects_a_mis_sized_load_and_an_empty_path() {
    let (mesh, _, dirichlet) = tensile_setup(0.01);
    let opts = SolveOptions::default();
    let err = solve_load_path(
        &mesh,
        &nh(),
        &[0.0; 5],
        &dirichlet,
        &opts,
        None,
        LoadPathOptions::default(),
    )
    .expect_err("a bad load length is an error, not a panic");
    assert!(matches!(err, LoadPathError::LoadSizeMismatch { .. }));
    // Zero steps is a legitimate request for an empty path, not an error.
    let empty = solve_load_path(
        &mesh,
        &nh(),
        &vec![0.0; mesh.dof_count()],
        &dirichlet,
        &opts,
        None,
        LoadPathOptions {
            steps: 0,
            max_cutbacks: 0,
        },
    )
    .expect("zero steps is valid");
    assert!(empty.steps.is_empty());
    assert_eq!(empty.reached(), 0.0);
    assert!(empty.last().is_none());
}

// --- Higher-order elements -------------------------------------------------
//
// Hex20/Hex27 support rests on one claim: the assembly was never Hex8-specific,
// it was only *typed* as Hex8. These tests make that claim checkable — a
// quadratic element must produce the same box geometry, reproduce the same
// closed form, and pass the same tangent checks a Hex8 element does. A generic
// refactor that quietly mis-assembled would show up here and nowhere else.

/// A structured box of element type `E`, with the same argument meaning as
/// [`hex_box`]: element counts, not grid cells.
fn box_of<E: ReferenceElement + crate::mesh::ElementFamily>(n: usize, l: f64) -> crate::Mesh<E> {
    hex_box_of::<E>(n, n, n, l, l, l).expect("box")
}

#[test]
fn a_quadratic_box_has_the_right_volume() {
    // The regression this guards is subtle and was real during this work: the
    // node spacing has to be `l / (n * s)` so the last node lands exactly on `l`.
    // Getting it wrong yields a box that is still rectangular, still positively
    // oriented, and uniformly the wrong size. Volume is the invariant that
    // catches that, so it is checked directly rather than through a solve.
    let l = 10.0;
    for (name, v) in [
        (
            "Hex8",
            box_of::<Hex8>(2, l).element_volume(0).expect("Hex8"),
        ),
        (
            "Hex20",
            box_of::<Hex20>(2, l).element_volume(0).expect("Hex20"),
        ),
        (
            "Hex27",
            box_of::<Hex27>(2, l).element_volume(0).expect("Hex27"),
        ),
    ] {
        let expected = (l / 2.0f64).powi(3);
        assert!(
            (v - expected).abs() < 1.0e-9 * expected,
            "{name} element volume {v} vs {expected}"
        );
    }
}

#[test]
fn a_quadratic_box_spans_the_requested_box() {
    // Corner-to-corner extent, which is what a caller actually asked for, plus
    // the node count a twice-refined grid implies.
    let l = 10.0;
    let m = box_of::<Hex20>(2, l);
    let (lo, hi) = m
        .nodes()
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), n| {
            (a.min(n.x), b.max(n.x))
        });
    assert!(lo.abs() < 1.0e-12, "min x {lo}");
    assert!((hi - l).abs() < 1.0e-12, "max x {hi}");
    // A 2x2x2 box of serendipity elements: 8 elements x 20 nodes = 160 slots,
    // 79 of which are shared between neighbouring elements, leaving 81 distinct
    // nodes. A raw twice-refined grid would hold 125 positions, but the
    // face-centre and body-centre ones are referenced by no `Hex20` element, and
    // an orphan node has an exactly zero stiffness row — which would make the
    // condensed system singular. The builder compacts them away, and this count
    // is what proves it did.
    assert_eq!(m.node_count(), 81, "Hex20 2x2x2 node count");
    assert_eq!(m.element_count(), 8);
    assert_eq!(m.dof_count(), 243);
    // Every node must belong to at least one element.
    let referenced: std::collections::HashSet<usize> =
        m.elements().iter().flatten().copied().collect();
    assert_eq!(
        referenced.len(),
        m.node_count(),
        "every node must be referenced by some element"
    );
}

#[test]
fn hex20_reproduces_the_uniaxial_closed_form() {
    // The load-bearing check. If the assembly were wrong for a quadratic element
    // — a dropped mid-edge term, a mis-ordered connectivity — the answer would be
    // wrong in a mesh-dependent way.
    //
    // What the measurement actually shows, and why the assertion is
    // mesh-*independence* rather than convergence: a uniform bar in uniaxial
    // tension has an affine deformation and a uniform stress, so full
    // integration is exact at any mesh. The error against the incompressible
    // closed form is therefore constant, and equal to this crate's known
    // penalty-incompressibility deviation (the `d1 = 0.5` compromise the Hex8
    // suite also documents). Measured 0.0920 at 2x2x2, 3x3x3 and 4x4x4 alike.
    //
    // A *changing* error is the failure signal: it would mean locking (which a
    // quadratic element does not fix — that needs the mixed u-p formulation) or
    // a mis-assembly. Constant-and-known is the passing case.
    let l: f64 = 10.0;
    let lam: f64 = 1.3;
    let model = nh();
    let expected = 2.0 * C10 * (lam - lam.powi(-2));
    let mut errors: Vec<(usize, f64)> = Vec::new();
    for n in [2usize, 3] {
        let mesh = box_of::<Hex20>(n, l);
        // Constraints are placed by *coordinate*, not by "the nodes of this
        // face". A Hex20 box has mid-edge nodes on the side faces that lie on
        // neither the top nor the bottom, and leaving them free is an
        // unrestrained rigid-body mode.
        let mut dirichlet: Vec<(usize, f64)> = Vec::new();
        let tol = 1.0e-9;
        for node in 0..mesh.node_count() {
            let p = mesh.nodes()[node];
            for c in 0..3 {
                if p.to_array()[c].abs() < tol {
                    dirichlet.push((mesh.dof(node, c), 0.0));
                }
            }
            if (p.y - l).abs() < tol {
                dirichlet.push((mesh.dof(node, 1), (lam - 1.0) * l));
            }
        }
        let opts = SolveOptions {
            assembly: AssemblyOptions::with_quadrature_order(
                crate::Mesh::<Hex20>::default_quadrature_order(),
            ),
            ..SolveOptions::default()
        };
        let result = solve_static(
            &mesh,
            &model,
            &vec![0.0; mesh.dof_count()],
            &dirichlet,
            &opts,
            None,
        )
        .expect("Hex20 solve converges");
        let internal =
            internal_force(&mesh, &model, &result.displacement, &opts.assembly).expect("assembly");
        let top = mesh.face_nodes(1, true);
        let reaction: f64 = top.iter().map(|&node| internal[mesh.dof(node, 1)]).sum();
        let nominal = reaction / (l * l);
        errors.push((n, (nominal - expected).abs() / expected));
    }
    for (n, e) in &errors {
        println!("Hex20 {n}x{n}x{n} relative error {e:.4}");
    }
    for w in errors.windows(2) {
        assert!(
            (w[1].1 - w[0].1).abs() < 1.0e-6,
            "error must be mesh-independent for a uniform state (a changing \
             error means locking or a mis-assembly): {} then {}",
            w[0].1,
            w[1].1
        );
    }
    // And it must be the known penalty-incompressibility deviation, not an
    // arbitrary number. The Hex8 suite carries the same `d1 = 0.5` compromise.
    assert!(
        errors.first().expect("a mesh").1 < 0.12,
        "Hex20 deviates from the incompressible closed form by {}",
        errors.first().expect("a mesh").1
    );
}

#[test]
fn hex20_tangent_is_the_derivative_of_the_hex20_residual() {
    // The strongest available check that the quadratic assembly is right, and the
    // one that needs no energy identity: the cheap `B^T A B` tangent against a
    // tangent obtained by differencing the *whole residual* w.r.t. every DOF. This
    // additionally pins a non-affine deformation, which a straight-box uniaxial
    // test never exercises.
    //
    // Deliberately *not* a patch-test energy check. The tempting identity
    // `E = int P : (F - I)` is exact only for a material whose energy is linear
    // in `F`; for a general hyperelastic `W` with `P = dW/dF` the missing term is
    // `int W`, so differencing it disagrees with the assembled force by a
    // constant. Differencing the residual avoids inventing a subtly wrong energy.
    let l = 4.0;
    let mesh = box_of::<Hex20>(1, l);
    let opts =
        AssemblyOptions::with_quadrature_order(crate::Mesh::<Hex20>::default_quadrature_order());
    let n = mesh.dof_count();
    let mut u = vec![0.0; n];
    for node in 0..mesh.node_count() {
        let p = mesh.nodes()[node];
        u[mesh.dof(node, 0)] = 0.1 * p.x + 0.02 * p.y * p.z;
        u[mesh.dof(node, 1)] = 0.05 * p.y - 0.01 * p.x * p.z;
        u[mesh.dof(node, 2)] = 0.03 * p.z + 0.015 * p.x * p.y;
    }
    let cheap = tangent_stiffness(&mesh, &nh(), &u, &opts).expect("tangent");
    let full = tangent_stiffness_numerical(&mesh, &nh(), &u, &opts).expect("numerical tangent");
    let diff = coo_max_abs_diff(&cheap, &full);
    // Scaled by the magnitude of the two matrices rather than absolute: a
    // quadratic element's entries span a wider range than Hex8's, so a fixed
    // threshold would be either vacuous or flaky.
    let scale = cheap
        .vals
        .iter()
        .chain(full.vals.iter())
        .map(|v| v.abs())
        .fold(0.0f64, f64::max)
        .max(1.0);
    assert!(
        diff < 1.0e-5 * scale,
        "Hex20 tangent strategies differ by {diff} (scale {scale})"
    );
}

#[test]
fn a_wrong_node_count_is_rejected_per_element_type() {
    let m = box_of::<Hex20>(1, 1.0);
    let mut elements = m.elements().to_vec();
    elements[0].truncate(8);
    let err = crate::Mesh::<Hex20>::from_parts(m.nodes().to_vec(), elements)
        .expect_err("an 8-node element in a Hex20 mesh must be rejected");
    assert!(matches!(
        err,
        MeshError::WrongElementNodeCount {
            element: 0,
            found: 8,
            expected: 20
        }
    ));
}

#[test]
fn a_quadratic_element_gets_a_higher_default_quadrature_order() {
    // Under-integrating a curved element's Jacobian is a silent accuracy loss,
    // so the floor is attached to the element type rather than left to the
    // caller to remember.
    assert_eq!(crate::Mesh::<Hex8>::default_quadrature_order(), 2);
    assert_eq!(crate::Mesh::<Hex20>::default_quadrature_order(), 3);
    assert_eq!(crate::Mesh::<Hex27>::default_quadrature_order(), 3);
}

// --- Tetrahedra -----------------------------------------------------------
//
// A tetrahedron is not a hexahedron with fewer nodes: it lives on a different
// reference domain and is integrated by a different rule. These check the mesh
// (positive volumes, conforming, no orphans) and the assembly, which are the two
// ways a tet can look fine and still be wrong.

#[test]
fn a_tet_box_fills_its_volume() {
    // Six tets per cell, each a sixth of the cell, summing to the cell exactly.
    // A wrong Kuhn ordering would still give positive volumes but not sum right.
    let l = 10.0;
    for n in [1usize, 2, 3] {
        let m = tet_box_of::<Tet4>(n, n, n, l, l, l).expect("tet box");
        assert_eq!(m.element_count(), 6 * n * n * n, "six tets per cell");
        let total: f64 = (0..m.element_count())
            .map(|e| m.element_volume(e).expect("positively oriented"))
            .sum();
        let expected = l * l * l;
        assert!(
            (total - expected).abs() < 1.0e-9 * expected,
            "Tet4 {n}^3 total volume {total} vs {expected}"
        );
    }
}

#[test]
fn a_quadratic_tet_box_also_fills_its_volume() {
    // Same for `Tet10`. Its mid-edge nodes sit at true edge midpoints, so a
    // quadratic element integrates its own geometry exactly; misplaced mid-nodes
    // would make the total drift.
    let l = 6.0;
    for n in [1usize, 2] {
        let m = tet_box_of::<Tet10>(n, n, n, l, l, l).expect("tet box");
        let total: f64 = (0..m.element_count())
            .map(|e| m.element_volume(e).expect("positively oriented"))
            .sum();
        let expected = l * l * l;
        assert!(
            (total - expected).abs() < 1.0e-9 * expected,
            "Tet10 {n}^3 total volume {total} vs {expected}"
        );
    }
}

#[test]
fn a_tet_box_is_conforming_and_has_no_orphans() {
    // Every node must belong to some element. An orphan has an exactly zero
    // stiffness row and makes the condensed system singular — the same failure
    // the serendipity hex box had, so it is worth checking here too.
    // Every node must belong to some element. An orphan has an exactly zero
    // stiffness row and makes the condensed system singular — the same failure
    // the serendipity hex box had, so it is worth checking here too.
    fn check<E: ReferenceElement + crate::mesh::ElementFamily>(name: &str, m: &crate::Mesh<E>) {
        let referenced: std::collections::HashSet<usize> =
            m.elements().iter().flatten().copied().collect();
        assert_eq!(
            referenced.len(),
            m.node_count(),
            "{name}: every node must be referenced by some element"
        );
        for e in 0..m.element_count() {
            assert_eq!(
                m.elements()[e].len(),
                E::NUM_NODES,
                "{name}: element {e} node count"
            );
        }
    }
    check(
        "Tet4",
        &tet_box_of::<Tet4>(2, 2, 2, 1.0, 1.0, 1.0).expect("Tet4"),
    );
    check(
        "Tet10",
        &tet_box_of::<Tet10>(2, 2, 2, 1.0, 1.0, 1.0).expect("Tet10"),
    );
}

#[test]
fn a_tet_box_shares_its_interface_nodes() {
    // The Kuhn subdivision exists to make internal faces conforming. If each of
    // the six tets created its own nodes, adjacent cells would share nothing and
    // the mesh would be a pile of disconnected cells.
    let m = tet_box_of::<Tet4>(2, 1, 1, 1.0, 1.0, 1.0).expect("tet box");
    // 2x1x1 cells have 3x2x2 = 12 grid corners, shared. Unshared it would be 16.
    assert_eq!(m.node_count(), 12, "2x1x1 cells share their corners");
    // One cell: 8 corners, and 6 tets x 6 edges = 36 edge slots over 19 distinct
    // edges (12 grid, 6 face diagonals, 1 body diagonal), all mid-nodes shared.
    let m10 = tet_box_of::<Tet10>(1, 1, 1, 1.0, 1.0, 1.0).expect("tet box");
    assert_eq!(m10.node_count(), 27, "8 corners + 19 distinct edges");
}

#[test]
fn tet_box_rejects_a_hexahedral_element() {
    // Asked the wrong question, this must be a named error rather than a mesh of
    // nonsense or a panic.
    assert!(matches!(
        tet_box_of::<Hex8>(1, 1, 1, 1.0, 1.0, 1.0),
        Err(MeshError::UnsupportedTetBox {
            nodes_per_element: 8
        })
    ));
    assert!(matches!(
        tet_box_of::<Hex20>(1, 1, 1, 1.0, 1.0, 1.0),
        Err(MeshError::UnsupportedTetBox {
            nodes_per_element: 20
        })
    ));
    assert!(matches!(
        tet_box_of::<Tet4>(0, 1, 1, 1.0, 1.0, 1.0),
        Err(MeshError::EmptyBox { axis: 0 })
    ));
}

#[test]
fn a_tet_box_reproduces_the_uniaxial_closed_form() {
    // End-to-end, same target as the hex suite. A uniform bar has an affine
    // deformation, so the error is again the known penalty-incompressibility
    // deviation and is mesh-independent; what is being verified is that the
    // simplex reference domain and simplex quadrature are applied correctly,
    // which a wrong rule would break badly.
    let l: f64 = 10.0;
    let lam: f64 = 1.3;
    let model = nh();
    let expected = 2.0 * C10 * (lam - lam.powi(-2));
    let opts = SolveOptions {
        assembly: AssemblyOptions::with_quadrature_order(
            crate::Mesh::<Tet4>::default_quadrature_order(),
        ),
        ..SolveOptions::default()
    };
    let mut errors = Vec::new();
    for n in [1usize, 2] {
        let m = tet_box_of::<Tet4>(n, n, n, l, l, l).expect("tet box");
        let mut dirichlet: Vec<(usize, f64)> = Vec::new();
        let tol = 1.0e-9;
        for node in 0..m.node_count() {
            let p = m.nodes()[node];
            for c in 0..3 {
                if p.to_array()[c].abs() < tol {
                    dirichlet.push((m.dof(node, c), 0.0));
                }
            }
            if (p.y - l).abs() < tol {
                dirichlet.push((m.dof(node, 1), (lam - 1.0) * l));
            }
        }
        let result = solve_static(
            &m,
            &model,
            &vec![0.0; m.dof_count()],
            &dirichlet,
            &opts,
            None,
        )
        .expect("tet solve converges");
        let internal =
            internal_force(&m, &model, &result.displacement, &opts.assembly).expect("assembly");
        let top = m.face_nodes(1, true);
        let reaction: f64 = top.iter().map(|&node| internal[m.dof(node, 1)]).sum();
        errors.push(((reaction / (l * l) - expected).abs() / expected, n));
    }
    for (e, n) in &errors {
        println!("Tet4 {n}^3 relative error {e:.4}");
    }
    assert!(
        errors.first().expect("a mesh").0 < 0.12,
        "Tet4 deviates from the closed form by {}",
        errors.first().expect("a mesh").0
    );
    for w in errors.windows(2) {
        assert!(
            (w[1].0 - w[0].0).abs() < 1.0e-6,
            "error must be mesh-independent for a uniform state: {} then {}",
            w[0].0,
            w[1].0
        );
    }
}

#[test]
fn a_tet_tangent_is_the_derivative_of_the_tet_residual() {
    // The same shape as the Hex20 tangent check, on a non-affine deformation, so
    // the simplex gradients and `J^-T` are exercised away from their linear
    // regime.
    let m = tet_box_of::<Tet10>(1, 1, 1, 4.0, 4.0, 4.0).expect("tet box");
    let opts =
        AssemblyOptions::with_quadrature_order(crate::Mesh::<Tet10>::default_quadrature_order());
    let n = m.dof_count();
    let mut u = vec![0.0; n];
    for node in 0..m.node_count() {
        let p = m.nodes()[node];
        u[m.dof(node, 0)] = 0.1 * p.x + 0.02 * p.y * p.z;
        u[m.dof(node, 1)] = 0.05 * p.y - 0.01 * p.x * p.z;
        u[m.dof(node, 2)] = 0.03 * p.z + 0.015 * p.x * p.y;
    }
    let cheap = tangent_stiffness(&m, &nh(), &u, &opts).expect("tangent");
    let full = tangent_stiffness_numerical(&m, &nh(), &u, &opts).expect("numerical tangent");
    let diff = coo_max_abs_diff(&cheap, &full);
    let scale = cheap
        .vals
        .iter()
        .chain(full.vals.iter())
        .map(|v| v.abs())
        .fold(0.0f64, f64::max)
        .max(1.0);
    assert!(
        diff < 1.0e-5 * scale,
        "Tet10 tangent strategies differ by {diff} (scale {scale})"
    );
}

#[test]
fn a_tet_element_uses_a_simplex_rule_not_a_tensor_product_one() {
    // The failure this guards is silent: integrating a simplex against a
    // tensor-product rule on `[-1, 1]^3` returns plausible numbers from points
    // that are not even inside the element. Checked on the rule's own domain —
    // a point of the reference tetrahedron has non-negative barycentric
    // coordinates summing to one.
    let rule = <Tet10 as crate::ElementFamily>::quadrature_rule(
        crate::Mesh::<Tet10>::default_quadrature_order(),
    );
    assert!(!rule.points.is_empty());
    for p in &rule.points {
        assert!(
            p[0] + p[1] + p[2] <= 1.0 + 1.0e-12,
            "tet quadrature point {p:?} lies outside the reference simplex"
        );
        assert!(p[0] >= -1.0e-12 && p[1] >= -1.0e-12 && p[2] >= -1.0e-12);
    }
    // Weights on the reference simplex sum to its volume, 1/6.
    let total: f64 = rule.weights.iter().sum();
    assert!(
        (total - 1.0 / 6.0).abs() < 1.0e-12,
        "tet rule weights sum to {total}, expected 1/6"
    );
    // The hex rule, for contrast, is on the cube and sums to 8.
    let hex = <Hex8 as crate::ElementFamily>::quadrature_rule(2);
    let hex_total: f64 = hex.weights.iter().sum();
    assert!((hex_total - 8.0).abs() < 1.0e-12);
}

// --- Volumetric locking: selective reduced integration ---------------------
//
// Option C from the u-p discussion: split the volumetric first Piola out of the
// fused `P` and integrate it on its own coarser rule. This reduces locking; it
// does not buy exact incompressibility, which only a mixed u-p formulation does.

#[test]
fn the_volumetric_split_is_exact() {
    // The method rests on `P = P_dev + P_vol` being a true decomposition, not an
    // approximation. For a body of known energy that is checkable: the
    // volumetric first Piola is the exact derivative of `(J-1)^2/d1`, and it
    // vanishes at `J = 1` whatever `d1` is.
    let model = nh();
    for lam in [0.9f64, 1.0, 1.1, 1.3] {
        let f = Mat3::from_array([lam, 0.0, 0.0, 0.0, lam, 0.0, 0.0, 0.0, lam]);
        let got = model.volumetric_piola(&f).at(0, 0);
        let expected = 2.0 * lam * lam * (lam.powi(3) - 1.0) / D1;
        assert!(
            (got - expected).abs() < 1.0e-9 * expected.abs().max(1.0),
            "lam={lam}: volumetric P {got} vs {expected}"
        );
    }
    assert_eq!(model.volumetric_piola(&Mat3::IDENTITY), Mat3::ZERO);
}

#[test]
fn reduced_volumetric_integration_reduces_locking() {
    // The point of the exercise, measured rather than asserted: a *stiff* penalty
    // (`d1` small) on a coarse mesh is the regime where locking dominates, and
    // that is where SRI should show its benefit. Target is the closed form
    // `mu (lambda - lambda^-2)`.
    let l: f64 = 10.0;
    let lam: f64 = 1.3;
    let expected = 2.0 * C10 * (lam - lam.powi(-2));
    // Deliberately stiffer than the fixture's `d1 = 0.5`: this is the
    // locking-dominated corner, which is where the method is supposed to help.
    let stiff = TissueModel::NeoHookean(NeoHookeanParams { c10: C10, d1: 0.1 });
    let error_for = |sri: Option<usize>| -> f64 {
        let mesh = hex_box(2, 2, 2, l, l, l).expect("box");
        let mut dirichlet: Vec<(usize, f64)> = Vec::new();
        for n in mesh.face_nodes(1, false) {
            for c in 0..3 {
                dirichlet.push((mesh.dof(n, c), 0.0));
            }
        }
        for n in mesh.face_nodes(1, true) {
            dirichlet.push((mesh.dof(n, 2), 0.0));
            dirichlet.push((mesh.dof(n, 1), (lam - 1.0) * l));
        }
        let mut assembly = AssemblyOptions::with_quadrature_order(3);
        assembly.volumetric_quadrature_order = sri;
        let opts = SolveOptions {
            assembly,
            ..SolveOptions::default()
        };
        let r = solve_static(
            &mesh,
            &stiff,
            &vec![0.0; mesh.dof_count()],
            &dirichlet,
            &opts,
            None,
        )
        .expect("converges");
        let internal =
            internal_force(&mesh, &stiff, &r.displacement, &opts.assembly).expect("assembly");
        let reaction: f64 = mesh
            .face_nodes(1, true)
            .iter()
            .map(|&n| internal[mesh.dof(n, 1)])
            .sum();
        ((reaction / (l * l) - expected) / expected).abs()
    };
    let full = error_for(None);
    let reduced = error_for(Some(1));
    println!("SRI locking error: full={full:.4} reduced={reduced:.4}");
    assert!(
        reduced < full,
        "reduced volumetric integration should reduce locking: {reduced} vs {full}"
    );
    // Substantially, not by a rounding error — a change too small to see would
    // mean the option is not actually wired to anything.
    assert!(
        reduced < 0.5 * full,
        "expected a substantial reduction, got {reduced} vs {full}"
    );
}

#[test]
fn reduced_volumetric_integration_is_a_no_op_for_a_law_without_one() {
    // The documented sharp edge, pinned so it stays documented: a law reporting
    // no volumetric part makes the option do nothing, silently. This is why the
    // option is opt-in.
    let l: f64 = 10.0;
    let model = nh();
    let closure = FnModel(|f: &Mat3| TissueModel::first_piola(&model, f));
    let mesh = hex_box(1, 1, 1, l, l, l).expect("box");
    let mut u = vec![0.0; mesh.dof_count()];
    for n in mesh.face_nodes(1, true) {
        u[mesh.dof(n, 1)] = 0.3;
    }
    let assembly = AssemblyOptions {
        volumetric_quadrature_order: Some(1),
        ..AssemblyOptions::default()
    };
    let f_with = internal_force(&mesh, &closure, &u, &assembly).expect("asm");
    let f_without = internal_force(&mesh, &closure, &u, &AssemblyOptions::default()).expect("asm");
    let diff: f64 = f_with
        .iter()
        .zip(&f_without)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        diff == 0.0,
        "an `FnModel` has no volumetric part, so the option must be a no-op; \
         it changed the force by {diff}"
    );
}
