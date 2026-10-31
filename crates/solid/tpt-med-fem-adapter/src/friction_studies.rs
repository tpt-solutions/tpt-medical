//! Friction sensitivity study (RFC 0004 acceptance item 3; RFC 0013
//! verification item 7).
//!
//! The crown ring is closed radially onto a wall whose Coulomb coefficient
//! `mu` is swept, and then translated axially as a rigid body against the
//! same wall. Two things are measured:
//!
//! - **radial force** at the closed state: the sensitivity of the
//!   deployment's headline number to friction. Radial closing produces no
//!   tangential slip by symmetry, so this should be insensitive to `mu`;
//!   the study records it rather than assuming it.
//! - **pull-out force**: the friction force the wall exerts on the ring as it
//!   slides. Every node that slips sits at the Coulomb bound `mu * f_n`, so
//!   the total is bounded by `mu * N` (N the sum of normal reactions) and
//!   reaches that bound once all contact nodes slide. A stuck node is held by
//!   the regularising spring, so the stick branch is below the bound.
//!
//! Translating every axial DOF by the same amount is a rigid motion of the
//! ring, so it carries no elastic stress: with `mu = 0` the pull-out force is
//! exactly zero, and any pull-out force is friction alone.
//!
//! Run the release-profile sweep with
//! `cargo test --release -p tpt-med-fem-adapter --features superelastic
//! --lib friction_studies -- --ignored --nocapture`.

use crate::deployment::{Deployment3D, StageOptions};
use crate::fixtures::CrownRingSpec;
use crate::friction::FrictionConfig;
use crate::solver::{Convergence, SolveOptions};
use crate::superelastic::{SouzaAuricchio, SuperelasticParams};
use crate::wall::{RadialWall, WallSide};

/// Normal penalty (N/mm), as in the ring studies.
const KAPPA: f64 = 1.0e5;
/// Regularising tangential stiffness per contact node (N/mm).
const K_T: f64 = 1.0e6;

/// One point of a pull-out curve.
#[derive(Debug, Clone, Copy)]
pub struct PullPoint {
    /// Rigid axial translation of the ring (mm).
    pub delta: f64,
    /// Total friction force magnitude opposing the translation (N).
    pub pull_force: f64,
    /// Sum of normal reactions at the wall (N).
    pub normal_total: f64,
    /// Contact nodes saturated at the Coulomb bound.
    pub slipping: usize,
    /// Contact nodes in contact.
    pub active: usize,
}

/// The outcome of one friction coefficient.
#[derive(Debug, Clone)]
pub struct FrictionRun {
    /// Radial force at the closed state (N).
    pub radial_force: f64,
    /// Pull-out curve, in the order the translations were applied.
    pub pull: Vec<PullPoint>,
}

/// Close the ring to `fraction` of its nominal diameter against a wall with
/// friction `mu`, then translate it axially through each of `deltas`
/// (ascending) and record the friction.
pub fn friction_run(spec: &CrownRingSpec, mu: f64, fraction: f64, deltas: &[f64]) -> FrictionRun {
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
    let friction = FrictionConfig::new(mu, K_T).expect("friction");
    let wall = RadialWall::new(
        2,
        [0.0; 2],
        r0 + 1e-3,
        WallSide::Inside,
        ring.outer_nodes.clone(),
    )
    .expect("wall")
    .with_friction(friction);
    dep.set_wall(wall, KAPPA);

    // Closing in staged radius steps: a single jump onto the refined ring at
    // mu = 0.3 exhausts the cutback budget (see the sweep's recorded output).
    let stage = StageOptions {
        steps: 2,
        max_cutbacks: 10,
    };
    let (start, target) = (r0 + 1e-3, r0 * (1.0 - fraction));
    for k in 1..=4 {
        dep.move_wall(start + (target - start) * k as f64 / 4.0)
            .expect("move");
        dep.prescribe(&ring.pins, stage, None).expect("close");
    }
    let radial_force = dep
        .wall_summary()
        .expect("summary")
        .expect("wall")
        .total_reaction;

    // Rigid axial translation: every node's axial DOF moves by the same
    // amount. The radial pins are kept, except their axial members, which the
    // translation supersedes.
    let axial: Vec<usize> = (0..ring.mesh.node_count())
        .map(|n| ring.mesh.dof(n, 2))
        .collect();
    let tangential: Vec<(usize, f64)> = ring
        .pins
        .iter()
        .copied()
        .filter(|(d, _)| !axial.contains(d))
        .collect();
    let slide = StageOptions {
        steps: 4,
        max_cutbacks: 10,
    };
    let mut pull = Vec::new();
    for &delta in deltas {
        let mut targets: Vec<(usize, f64)> = axial.iter().map(|&d| (d, delta)).collect();
        targets.extend(tangential.iter().copied());
        dep.prescribe(&targets, slide, None).expect("slide");
        let friction_forces = dep.wall_friction().expect("friction").expect("wall");
        let s = dep.wall_summary().expect("summary").expect("wall");
        // Friction opposes the motion, so the axial components share a sign.
        let pull_force = friction_forces.iter().map(|f| f[2]).sum::<f64>().abs();
        pull.push(PullPoint {
            delta,
            pull_force,
            normal_total: s.total_reaction,
            slipping: s.slipping,
            active: s.active,
        });
    }
    FrictionRun { radial_force, pull }
}

fn tiny() -> CrownRingSpec {
    let mut s = CrownRingSpec::new(4);
    s.stations_per_crown = 3;
    s.width_elements = 1;
    s.thickness_elements = 1;
    s
}

#[test]
fn frictionless_pull_out_is_exactly_zero() {
    let run = friction_run(&tiny(), 0.0, 0.02, &[0.01, 0.03]);
    for p in &run.pull {
        assert!(p.active > 0, "fixture must be in contact");
        assert!(p.pull_force.abs() < 1e-12, "{p:?}");
    }
}

#[test]
fn pull_out_obeys_the_coulomb_bound_and_friction_stiffens_radial_closing() {
    let mu = 0.3;
    let frictionless = friction_run(&tiny(), 0.0, 0.02, &[]);
    let run = friction_run(&tiny(), mu, 0.02, &[0.005, 0.02]);
    // Crown tips that are held tangentially cannot slide into the free
    // configuration, so friction stiffens the closing ring. The sign is the
    // physics; the magnitude is what the sweep reports.
    println!(
        "radial: mu=0 {:.4} N, mu={mu} {:.4} N",
        frictionless.radial_force, run.radial_force
    );
    assert!(
        run.radial_force > frictionless.radial_force,
        "radial {} vs frictionless {}",
        run.radial_force,
        frictionless.radial_force
    );
    for p in &run.pull {
        assert!(p.pull_force > 0.0, "{p:?}");
        // Each slipping node contributes at most mu * f_n, so the total is
        // bounded by mu times the sum of normal reactions.
        assert!(
            p.pull_force <= mu * p.normal_total * (1.0 + 1e-6),
            "pull {} exceeds the Coulomb bound {}",
            p.pull_force,
            mu * p.normal_total
        );
    }
    let last = run.pull.last().expect("points");
    assert!(
        last.slipping > 0,
        "the larger translation should slide some nodes: {last:?}"
    );
}

/// The friction sensitivity sweep reported in the RFC.
#[test]
#[ignore = "release-profile study; many minutes; see module docs"]
fn friction_sensitivity_sweep() {
    let mut spec = CrownRingSpec::new(6);
    spec.stations_per_crown = 12;
    spec.width_elements = 1;
    spec.thickness_elements = 2;
    let fraction = 0.02;
    let deltas = [0.005, 0.01, 0.02, 0.05];
    let base = friction_run(&spec, 0.0, fraction, &deltas);
    println!("mu   F_radial(N)  dF_radial  delta(mm)  pull(N)  mu*N(N)  pull/(mu*N)  slip/active");
    // mu = 0.5 is not included: the staged closing does not converge there on
    // this mesh (cutback budget exhausted), so no number is recorded for it.
    for mu in [0.0, 0.1, 0.2, 0.3] {
        let run = friction_run(&spec, mu, fraction, &deltas);
        let rel = (run.radial_force - base.radial_force) / base.radial_force;
        for p in &run.pull {
            let bound = mu * p.normal_total;
            let ratio = if bound > 0.0 {
                p.pull_force / bound
            } else {
                0.0
            };
            println!(
                "{mu:.1} {:11.5} {:+9.2e} {:9.3} {:8.4} {:8.4} {:11.4} {}/{}",
                run.radial_force, rel, p.delta, p.pull_force, bound, ratio, p.slipping, p.active
            );
        }
    }
}

/// ASME V&V 40 assessment of the question "how much does wall friction move
/// the radial force of a crown-ring deployment?" (RFC 0004 acceptance item 3;
/// RFC 0013 verification item 7).
///
/// The risk and influence are proposed (Medium risk, Significant influence)
/// and need confirming by the maintainer. The acceptance bounds are proposed
/// too: a 10 % radial-force sensitivity and a 0.999 Coulomb ratio on slipping
/// contact. Validation is left empty on purpose: no bench friction data
/// exists in the workspace, so the assessment must report that gap rather
/// than close it.
#[test]
#[ignore = "release-profile study; minutes; see module docs"]
fn vv40_assessment_of_the_friction_question() {
    use tpt_med_vv40::{
        Acceptance, CredibilityAssessment, EvidenceMetric, ModelInfluence, ModelRisk,
        VerificationActivity, VerificationType,
    };

    let mut spec = CrownRingSpec::new(6);
    spec.stations_per_crown = 12;
    spec.width_elements = 1;
    spec.thickness_elements = 2;
    let deltas = [0.01, 0.02, 0.05];
    let mu = 0.3;
    let base = friction_run(&spec, 0.0, 0.02, &deltas);
    let rough = friction_run(&spec, mu, 0.02, &deltas);

    let frictionless_pull = base
        .pull
        .iter()
        .map(|p| p.pull_force.abs())
        .fold(0.0, f64::max);
    // Taken at the largest translation, where the sliding nodes have
    // saturated; the first step is still in the stick-to-slip transition.
    let last = rough.pull.last().expect("points");
    let coulomb_ratio = last.pull_force / (mu * last.normal_total);
    let sensitivity = (rough.radial_force - base.radial_force).abs() / base.radial_force;

    let mut code = VerificationActivity {
        kind: VerificationType::CodeVerification,
        description: "Zero-friction limit and Coulomb bound on rigid axial translation".into(),
        results: "friction_studies::frictionless_pull_out_is_exactly_zero and \
                  pull_out_obeys_the_coulomb_bound (release run of the refined ring)"
            .into(),
        metrics: Vec::new(),
    };
    code.with_metric(EvidenceMetric {
        name: "frictionless_pull_out_n".into(),
        value: frictionless_pull,
        acceptance: Acceptance::at_most(1e-12),
    })
    .with_metric(EvidenceMetric {
        name: "coulomb_ratio_when_sliding".into(),
        value: coulomb_ratio,
        acceptance: Acceptance {
            max: None,
            min: Some(0.999),
        },
    });
    let mut sens = VerificationActivity {
        kind: VerificationType::SensitivityAnalysis,
        description: "Radial force against wall Coulomb coefficient, mu = 0 to 0.3".into(),
        results: "friction_sensitivity_sweep (2 % closing, 6 crowns x 12 stations x 1 x 2 section)"
            .into(),
        metrics: Vec::new(),
    };
    sens.with_metric(EvidenceMetric {
        name: "radial_force_sensitivity_mu_0_3".into(),
        value: sensitivity,
        acceptance: Acceptance::at_most(0.10),
    });

    let assessment = CredibilityAssessment {
        question_of_interest: "How much does wall friction change the radial force of a \
                               crown-ring stent deployment?"
            .into(),
        risk: ModelRisk::Medium,
        influence: ModelInfluence::Significant,
        verification: vec![code.clone(), sens.clone()],
        validation: Vec::new(),
    };

    assert!(code.metrics_adequate(), "{code:?}");
    assert!(sens.metrics_adequate(), "{sens:?}");
    // Verification goals are met; the validation gap is the finding.
    assert_eq!(
        assessment.evaluate(),
        vec!["validation coverage 0 < 1 required".to_string()],
        "{}",
        assessment.to_json()
    );
    assert!(!assessment.is_credible());
    println!("{}", assessment.to_json());
}
