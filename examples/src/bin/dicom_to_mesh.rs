//! Phase 1 milestone CLI: DICOM CT series → patient-specific bone mesh CSV.
//!
//! ```console
//! cargo run -p tpt-med-examples --bin dicom-to-mesh -- \
//!     test-data/dicom/synthetic_ct --output femur_mesh.csv --threshold 200
//! ```

use std::path::PathBuf;

use tpt_med_dicom::{DicomSeries, HounsfieldMapper};
use tpt_med_meshing::{smooth_mesh, MedicalMesher, SegmentationMask};

struct Args {
    input: PathBuf,
    output: PathBuf,
    threshold: f64,
    smooth_iterations: u32,
}

fn parse_args() -> Option<Args> {
    let mut input = None;
    let mut output = PathBuf::from("bone_mesh.csv");
    let mut threshold = HounsfieldMapper::DEFAULT_BONE_THRESHOLD_HU;
    let mut smooth_iterations = 0;
    let mut it = std::env::args().skip(1);
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--output" => output = PathBuf::from(it.next()?),
            "--threshold" => threshold = it.next()?.parse().ok()?,
            "--smooth" => smooth_iterations = it.next()?.parse().ok()?,
            "--help" | "-h" => return None,
            other if input.is_none() => input = Some(PathBuf::from(other)),
            _ => return None,
        }
    }
    Some(Args {
        input: input?,
        output,
        threshold,
        smooth_iterations,
    })
}

fn main() {
    let Some(args) = parse_args().or_else(|| {
        eprintln!(
            "usage: dicom-to-mesh <dicom-dir> [--output mesh.csv] \
             [--threshold 200] [--smooth 5]"
        );
        std::process::exit(2);
    }) else {
        unreachable!()
    };

    println!("loading DICOM series from {} ...", args.input.display());
    let series = match DicomSeries::load_from_dir(&args.input) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    };
    let (nx, ny, nz) = series.dims();
    println!(
        "  series uid: {}",
        if series.series_uid.is_empty() {
            "<none>"
        } else {
            &series.series_uid
        }
    );
    println!(
        "  modality: {}, volume: {nx}x{ny}x{nz} voxels",
        series.modality
    );
    println!(
        "  spacing: {:.3} x {:.3} x {:.3} mm",
        series.pixel_spacing.1, series.pixel_spacing.0, series.slice_thickness
    );

    println!("segmenting bone at HU >= {} ...", args.threshold);
    let mask = SegmentationMask::threshold_hu(&series, args.threshold);
    println!("  {} bone voxels", mask.solid_count());

    println!("meshing voxels to hexahedra ...");
    let mesher = MedicalMesher::default();
    let mut mesh = match mesher.voxels_to_hex_mesh(&mask) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    };
    if args.smooth_iterations > 0 {
        println!(
            "smoothing surface ({} iterations) ...",
            args.smooth_iterations
        );
        smooth_mesh(&mut mesh, args.smooth_iterations, 0.5);
    }
    println!(
        "  {} nodes, {} hex elements, E in [{:.0}, {:.0}] MPa (mean {:.0})",
        mesh.nodes.len(),
        mesh.elements.len(),
        mesh.materials
            .iter()
            .map(|m| m.youngs_modulus)
            .fold(f64::INFINITY, f64::min),
        mesh.max_modulus(),
        mesh.mean_modulus()
    );

    if let Err(e) = mesh.write_csv(&args.output) {
        eprintln!("error writing {}: {e}", args.output.display());
        std::process::exit(1);
    }
    println!("wrote {}", args.output.display());
}
