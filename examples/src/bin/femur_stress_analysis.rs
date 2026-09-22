//! Phase 2 milestone: femur stress analysis under physiological loading.
//!
//! Builds the synthetic proximal-femur model from the DICOM test data,
//! meshes it, assigns HU-based materials, applies a 3× body-weight stance
//! load, and reports peak von Mises stress and displacement.
//!
//! ```console
//! cargo run -p tpt-med-examples --bin femur-stress-analysis
//! ```

use tpt_med_biomechanics::{BiomechanicsModel, BoundaryConditions};
use tpt_med_dicom::DicomSeries;
use tpt_med_meshing::{smooth_mesh, MedicalMesher, SegmentationMask};
use tpt_med_units::Force;

fn main() {
    let ct_dir = std::path::Path::new("test-data/dicom/synthetic_ct");
    let body_weight = Force::from_n(75.0 * 9.81);

    println!("=== tpt-medical femur stress analysis (Phase 2 milestone) ===");
    println!("loading synthetic CT ...");
    let series =
        DicomSeries::load_from_dir(ct_dir).expect("test data present (run gen-synthetic-ct)");

    println!("segmenting + meshing ...");
    let mask = SegmentationMask::threshold_hu(&series, 200.0);
    let mut mesh = MedicalMesher::default()
        .voxels_to_hex_mesh(&mask)
        .expect("phantom contains bone");
    smooth_mesh(&mut mesh, 3, 0.5);
    println!(
        "  {} nodes, {} elements, E ∈ [{:.0}, {:.0}] MPa",
        mesh.nodes.len(),
        mesh.elements.len(),
        mesh.materials
            .iter()
            .map(|m| m.youngs_modulus)
            .fold(f64::INFINITY, f64::min),
        mesh.max_modulus()
    );

    let model = BiomechanicsModel::from_voxel_mesh(&mesh);

    // Boundary conditions: distal end fully fixed; stance joint reaction
    // applied over the proximal head region (top 15% of the model) as a
    // uniform traction equivalent.
    let (min_z, max_z) = {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for n in &model.nodes {
            lo = lo.min(n.z);
            hi = hi.max(n.z);
        }
        (lo, hi)
    };
    let load_band = min_z + 0.85 * (max_z - min_z);

    let mut bc = BoundaryConditions::default();
    let mut fixed = 0usize;
    let mut loaded = Vec::new();
    for (i, n) in model.nodes.iter().enumerate() {
        if n.z <= min_z + 1e-6 {
            bc.fix_nodes([i as u32]);
            fixed += 1;
        }
        if n.z >= load_band {
            loaded.push(i as u32);
        }
    }
    let joint_reaction = Force::body_weights(3.0, body_weight);
    let per_node = joint_reaction.to_n() / loaded.len().max(1) as f64;
    for &n in &loaded {
        bc.add_force(n, tpt_med_geometry::Vec3::new(0.0, 0.0, -per_node));
    }
    println!(
        "loading: 3.0×BW stance reaction = {:.2} kN over {} proximal nodes ({} fixed distal)",
        joint_reaction.to_kn(),
        loaded.len(),
        fixed
    );

    println!("solving (CG) ...");
    let result = model.solve(&bc, 1e-8, 30_000).expect("solver converges");
    println!(
        "  solver: {} iterations, residual {:.2e}",
        result.stats.iterations, result.stats.relative_residual
    );
    println!("  max displacement: {:.4} mm", result.max_displacement());
    println!("  max von Mises: {:.2} MPa", result.max_von_mises());
    println!("  mean von Mises: {:.3} MPa", result.mean_von_mises());
    if let Some(el) = result.critical_element() {
        println!("  critical element: {el}");
    }

    // Simple screening check: cortical yield ~110 MPa (see tpt-med-bone).
    let vm = result.max_von_mises();
    println!(
        "screening vs cortical yield (110 MPa): {:.0}% of yield",
        100.0 * vm / 110.0
    );
}
