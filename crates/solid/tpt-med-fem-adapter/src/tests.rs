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

use crate::assembly::{
    coo_max_abs_diff, internal_force, material_tangent, tangent_stiffness,
    tangent_stiffness_numerical, AssemblyOptions, Constitutive,
};
use crate::contact::{ContactError, ContactPairing};
use crate::mesh::{hex_box, MeshError};
use crate::solver::{solve_static, ContactConfig, SolveError, SolveOptions};
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
            pairing: &s.pairing,
            penalty,
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
            pairing: &s.pairing,
            penalty: 1.0e4,
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
            pairing: &s.pairing,
            penalty: 1.0e4,
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
        pairing: &s.pairing,
        penalty,
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
