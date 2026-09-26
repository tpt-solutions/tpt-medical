# NIfTI test data

`tpt-med-nifti` (RFC 0006) does not commit fixtures here — its test suite
builds synthetic `.nii` byte buffers in-process via `SyntheticNiftiBuilder`
(`crates/imaging/tpt-med-nifti/src/synthetic.rs`), the same way
`tpt-med-dicom`'s tests use `SyntheticCtBuilder` rather than a committed
file. This directory is reserved should a hand-built binary fixture (e.g. one
exercising a real-world writer's quirk) ever be worth committing instead of
generating in-process.
