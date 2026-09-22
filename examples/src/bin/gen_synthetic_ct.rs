//! Generates the committed synthetic CT test data (never real patient data).
//!
//! ```console
//! cargo run -p tpt-med-examples --bin gen-synthetic-ct -- \
//!     test-data/dicom/synthetic_ct --cols 48 --rows 48 --slices 24
//! ```

use tpt_med_dicom::synthetic;

fn main() {
    let mut args = std::env::args().skip(1);
    let Some(output) = args.next().map(std::path::PathBuf::from) else {
        eprintln!("usage: gen-synthetic-ct <output-dir> [--cols N] [--rows N] [--slices N]");
        std::process::exit(2);
    };
    let mut get_opt = || -> Option<(String, String)> {
        let k = args.next()?;
        Some((k, args.next()?))
    };
    let (mut cols, mut rows, mut slices) = (48usize, 48usize, 24usize);
    while let Some((k, v)) = get_opt() {
        match k.as_str() {
            "--cols" => cols = v.parse().expect("cols"),
            "--rows" => rows = v.parse().expect("rows"),
            "--slices" => slices = v.parse().expect("slices"),
            other => {
                eprintln!("unknown option {other}");
                std::process::exit(2);
            }
        }
    }

    println!(
        "generating synthetic femur phantom {cols}x{rows}x{slices} into {} ...",
        output.display()
    );
    let series = synthetic::femur_phantom(cols, rows, slices);
    if let Err(e) = series.write_to_dir(&output) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
    println!("wrote {} slice files", series.slices.len());
}
