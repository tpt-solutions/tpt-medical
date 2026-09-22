//! DICOM parsing + meshing pipeline benchmark.
//!
//! Run: `cargo bench -p tpt-med-benches --bench dicom-parsing`
//! (harness-free: prints timings; CI archives the output).

use std::time::Instant;

use tpt_med_dicom::{DicomParser, DicomSeries};
use tpt_med_meshing::{smooth_mesh, MedicalMesher, SegmentationMask};

fn main() {
    let mut timings = Vec::new();

    // Synthetic series in memory (no disk I/O in the measurement path).
    let series = tpt_med_dicom::synthetic::femur_phantom(48, 48, 24);

    let t0 = Instant::now();
    let mut parsed = Vec::new();
    for s in &series.slices {
        parsed.push(DicomParser::parse_bytes(&s.bytes).expect("parse"));
    }
    timings.push(("parse 24 slices", t0.elapsed()));

    let ct = DicomSeries::from_slices(parsed).expect("series");
    let t1 = Instant::now();
    let mask = SegmentationMask::threshold_hu(&ct, 200.0);
    timings.push(("threshold segmentation", t1.elapsed()));

    let t2 = Instant::now();
    let mesher = MedicalMesher::default();
    let mut mesh = mesher.voxels_to_hex_mesh(&mask).expect("meshes");
    timings.push(("voxel-to-hex meshing", t2.elapsed()));

    let t3 = Instant::now();
    smooth_mesh(&mut mesh, 5, 0.5);
    timings.push(("laplacian smoothing x5", t3.elapsed()));

    println!("dicom-parsing benchmark (48x48x24 phantom):");
    for (name, d) in &timings {
        println!(
            "  {name:<28} {:>10.2?}  ({:.2} ms)",
            d,
            d.as_secs_f64() * 1e3
        );
    }
    let sum: std::time::Duration = timings.iter().map(|d| d.1).sum();
    println!("  {:<28} {:>10.2?}", "total", sum);
}
