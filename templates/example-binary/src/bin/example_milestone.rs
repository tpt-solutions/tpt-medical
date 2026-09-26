//! TODO: one line describing what this milestone demonstrates.
//!
//! This binary is the integration test for a whole pipeline, not a snippet.
//! It runs offline against the committed synthetic CT in
//! `test-data/dicom/synthetic_ct/` in a few seconds, and it exits non-zero on
//! failure — that is what makes it usable in CI.
//!
//! ```console
//! cargo run -p tpt-med-examples --bin <name> -- \
//!     test-data/dicom/synthetic_ct --output out.csv
//! ```
//!
//! Run with `--help` for the flags.

use std::path::{Path, PathBuf};

/// Long-option parsing, matching the existing milestone binaries.
///
/// Hand-rolled rather than pulled from a crate: the workspace keeps its
/// dependency footprint near zero so the WASM build stays small, and seven
/// binaries do not justify an argument-parsing dependency.
fn parse_args() -> Result<Args, String> {
    let mut args = Args::default();
    let mut it = std::env::args().skip(1);
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--output" => {
                args.output = Some(PathBuf::from(it.next().ok_or("--output needs a path")?))
            }
            "--help" | "-h" => return Err(usage()),
            other => return Err(format!("unknown argument: {other}\n{}", usage())),
        }
    }
    if args.output.is_none() {
        return Err(format!("--output is required\n{}", usage()));
    }
    Ok(args)
}

struct Args {
    output: Option<PathBuf>,
}

impl Default for Args {
    fn default() -> Self {
        Self { output: None }
    }
}

fn usage() -> String {
    "usage: <name> <dicom-dir> --output <file>".to_string()
}

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };

    // The input directory is the first positional argument, if given.
    let input: Option<&Path> = std::env::args().nth(1).map(Path::new);

    println!("=== tpt-medical: TODO ===");

    // TODO: the pipeline. Follow the shape of the existing milestones:
    //   parse -> threshold -> mesh -> solve -> post-process -> report.
    //
    // Print results as a labelled table on stdout so the binary can be piped
    // into another tool, and print units with every value.
    if let Some(dir) = input {
        println!("  input: {}", dir.display());
    }
    let out = args.output.unwrap();
    println!("  output: {}", out.display());

    // A failure must be a non-zero exit, never a printed error and a success
    // status. That is the only reason this binary is worth having in CI.
    if let Err(e) = std::fs::write(&out, b"") {
        eprintln!("error: cannot write {}: {e}", out.display());
        std::process::exit(1);
    }
}
