// Add the binary to the `examples` package manifest, then run:
//
//   cargo build -p tpt-med-examples --bin <name>
//
// Every milestone binary follows the same shape, so copy the manifest stanza
// from a neighbouring `[[bin]]` entry and change the name and path. Note the
// path uses the *file* name, which is snake_case, while the binary name is
// kebab-case:
//
//   [[bin]]
//   name = "example-milestone"
//   path = "src/bin/example_milestone.rs"
//
// A milestone binary is the integration test for a pipeline. Three rules make
// it worth having:
//
//   1. It runs offline, against `test-data/dicom/synthetic_ct/`, in seconds.
//   2. It exits non-zero on failure. There is no path that prints an error and
//      returns success, because that is the only thing CI can act on.
//   3. It prints a labelled table on stdout, with units on every value.
//
// If the milestone needs a new golden dataset under `test-data/golden/`, add
// the JSON in the documented `golden/v1` shape and run
// `scripts/diff-golden.sh` to show the numerical drift.
