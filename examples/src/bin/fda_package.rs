//! Phase 7 milestone: FDA-submission-ready simulation package export.
//!
//! Runs an end-to-end screening analysis (synthetic CT → mesh → FEM under
//! stance loading) while recording every operator action into a 21 CFR
//! Part 11 audit trail, applies an electronic review signature, performs
//! an ASME V&V 40 credibility assessment, and exports the signed package.
//!
//! ```console
//! cargo run -p tpt-med-examples --bin fda-package
//! ```

use tpt_med_audit::{hmac_sha256, sha256_hex};
use tpt_med_biomechanics::{BiomechanicsModel, BoundaryConditions};
use tpt_med_core::{AuditAction, AuditEvent};
use tpt_med_dicom::DicomSeries;
use tpt_med_fda::{AuditTrail, ReproducibilityManifest, SignatureMeaning};
use tpt_med_meshing::{MedicalMesher, SegmentationMask};
use tpt_med_vv40::{
    AgreementLevel, CredibilityAssessment, ModelInfluence, ModelRisk, ValidationActivity,
    ValidationType, VerificationActivity, VerificationType,
};

/// Crates whose code took part in this run, recorded in the reproducibility
/// manifest.
///
/// Each crate's *resolved* version comes from `Cargo.lock` via
/// [`tpt_med_examples::crate_versions::version_of`], not from this binary's
/// own `CARGO_PKG_VERSION`. Every crate here pins to the workspace version
/// today (`version.workspace = true`), so the two happen to agree — but the
/// lookup is already correct for the day a crate is released independently
/// and its version diverges from the rest.
const PARTICIPATING_CRATES: &[&str] = &[
    "tpt-med-audit",
    "tpt-med-biomechanics",
    "tpt-med-core",
    "tpt-med-dicom",
    "tpt-med-fda",
    "tpt-med-meshing",
    "tpt-med-units",
    "tpt-med-vv40",
];

/// SHA-256 over the whole input series, in filename order.
///
/// Hashing the concatenation rather than each file separately means the
/// manifest records "this exact series" without a per-file list, and any
/// change to any slice changes the digest.
fn hash_series(dir: &std::path::Path) -> (String, u64) {
    let mut files: Vec<_> = std::fs::read_dir(dir)
        .expect("test data present")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "dcm"))
        .collect();
    files.sort();

    let mut buf = Vec::new();
    for f in &files {
        buf.extend_from_slice(&std::fs::read(f).expect("slice readable"));
    }
    (sha256_hex(&buf), buf.len() as u64)
}

fn main() {
    println!("=== tpt-medical FDA package export (Phase 7 milestone) ===");

    // The submission secret would come from the operator's HSM/key store;
    // the demo uses a fixed synthetic key.
    let signing_key: [u8; 32] = hmac_sha256(b"tpt-medical-demo", b"submission-key");

    let mut trail = AuditTrail::new("sim-run-2026-09-20-001");
    trail.append(AuditEvent::new(
        "operator:engineer-1",
        "patient_model",
        "synthetic-femur",
        AuditAction::Create,
        "loaded synthetic CT series (no PHI)",
    ));

    // Reproducibility manifest: which code, which inputs, which build.
    // Attached before any computation so the manifest is on the record even if
    // the run later fails.
    let series_dir = std::path::Path::new("test-data/dicom/synthetic_ct");
    let (series_sha, series_bytes) = hash_series(series_dir);
    let manifest = ReproducibilityManifest::new(env!("CARGO_PKG_VERSION"))
        .with_crates(
            PARTICIPATING_CRATES
                .iter()
                .map(|c| (*c, tpt_med_examples::crate_versions::version_of(c))),
        )
        .with_git_commit(std::env::var("TPT_GIT_COMMIT").ok())
        .with_build_profile(std::env::var("PROFILE").unwrap_or_else(|_| "unknown".into()))
        .with_artifact(tpt_med_fda::InputArtifact {
            label: "ct_series".into(),
            sha256: series_sha.clone(),
            bytes: series_bytes,
        });
    let manifest_digest = manifest.digest();
    trail.attach_manifest(manifest);
    println!("[0/4] reproducibility manifest {manifest_digest}");
    println!("  input ct_series: {series_sha} ({series_bytes} bytes)");

    trail.append(AuditEvent::new(
        "operator:engineer-1",
        "mesh",
        "synthetic-femur",
        AuditAction::Modify,
        "threshold 200 HU, hex meshing, 3 smoothing iterations",
    ));

    println!("[1/4] building the screening model ...");
    let series = DicomSeries::load_from_dir(std::path::Path::new("test-data/dicom/synthetic_ct"))
        .expect("test data present");
    let mask = SegmentationMask::threshold_hu(&series, 200.0);
    let mesh = MedicalMesher::default()
        .voxels_to_hex_mesh(&mask)
        .expect("bone present");
    let model = BiomechanicsModel::from_voxel_mesh(&mesh);

    let mut bc = BoundaryConditions::default();
    let mut loaded = 0usize;
    let (mut min_z, mut max_z) = (f64::INFINITY, f64::NEG_INFINITY);
    for n in &model.nodes {
        min_z = min_z.min(n.z);
        max_z = max_z.max(n.z);
    }
    for (i, n) in model.nodes.iter().enumerate() {
        if n.z <= min_z + 1e-6 {
            bc.fix_nodes([i as u32]);
        }
        if n.z >= min_z + 0.85 * (max_z - min_z) {
            bc.add_force(
                i as u32,
                tpt_med_geometry::Vec3::new(0.0, 0.0, -2210.0 / 2100.0),
            );
            loaded += 1;
        }
    }

    trail.append(AuditEvent::new(
        "operator:engineer-1",
        "simulation",
        "stance-3bw",
        AuditAction::Simulate,
        "3x body-weight stance loading per ISO 7206 loading convention",
    ));
    println!("[2/4] solving ({loaded} loaded nodes) ...");
    let result = model.solve(&bc, 1e-8, 30_000).expect("solver converges");
    let peak = result.max_von_mises();
    println!("  peak von Mises: {peak:.2} MPa");

    // ASME V&V 40 credibility for this screening question.
    let assessment = CredibilityAssessment {
        question_of_interest: String::from(
            "Does peak cortical stress under 3x body-weight stance loading remain below yield?",
        ),
        risk: ModelRisk::Medium,
        influence: ModelInfluence::Significant,
        verification: vec![
            VerificationActivity {
                kind: VerificationType::CodeVerification,
                description: "analytical uniaxial + patch tests (crate test suite)".into(),
                results: "exact within 1e-9 mm".into(),
                metrics: vec![tpt_med_vv40::EvidenceMetric {
                    name: "patch_test_max_error_mm".into(),
                    value: 1e-9,
                    acceptance: tpt_med_vv40::Acceptance::at_most(1e-6),
                }],
            },
            VerificationActivity {
                kind: VerificationType::CalculationVerification,
                description: "cantilever vs Euler-Bernoulli with locking band".into(),
                results: "ratio 0.85 of analytical, inside documented band".into(),
                metrics: vec![tpt_med_vv40::EvidenceMetric {
                    name: "deflection_ratio_vs_analytic".into(),
                    value: 0.85,
                    acceptance: tpt_med_vv40::Acceptance {
                        min: Some(0.5),
                        max: Some(1.15),
                    },
                }],
            },
        ],
        validation: vec![ValidationActivity {
            kind: ValidationType::BenchmarkModel,
            reference: "synthetic femur phantom reference (test-data/golden/solid)".into(),
            metrics: vec!["peak von Mises vs golden".into()],
            agreement: AgreementLevel::Quantitative,
        }],
    };
    trail.append(AuditEvent::new(
        "operator:engineer-1",
        "credibility",
        "vv40-assessment",
        AuditAction::Approve,
        "risk=medium influence=significant credibility evaluation",
    ));
    println!(
        "[3/4] V&V 40 assessment: credible = {}",
        assessment.is_credible()
    );
    for gap in assessment.evaluate() {
        println!("  unmet: {gap}");
    }

    // Reviewer signature (§11.50).
    trail.sign("reviewer:qa-2", SignatureMeaning::Reviewer, &signing_key);
    assert!(trail.verify_integrity(), "trail integrity");

    // Export.
    println!("[4/4] exporting signed package ...");
    let (json, tag) = trail.export_package(&signing_key);
    let out_dir = std::path::Path::new("test-data/golden/regulatory");
    std::fs::create_dir_all(out_dir).ok();
    let json_path = out_dir.join("fda_export_example.json");
    std::fs::write(&json_path, &json).ok();
    println!("  entries: {}", trail.len());
    println!("  package json: {}", json_path.display());
    println!("  detached HMAC-SHA256 tag: {tag}");

    // Summary block appended to the golden record.
    println!("--- submission summary ---");
    println!("  question: {}", assessment.question_of_interest);
    println!("  peak stress {peak:.2} MPa vs cortical yield 110 MPa");
    println!("  audit entries signed: {}", trail.len());
}
