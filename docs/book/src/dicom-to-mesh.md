# DICOM to bone mesh

The Phase-1 milestone pipeline:

```rust
use tpt_med_dicom::DicomSeries;
use tpt_med_meshing::{MedicalMesher, SegmentationMask};

let series = DicomSeries::load_from_dir(std::path::Path::new("ct/"))?;
let mask = SegmentationMask::threshold_hu(&series, 200.0); // HU threshold
let mut mesh = MedicalMesher::default().voxels_to_hex_mesh(&mask)?;
smooth_mesh(&mut mesh, 3, 0.5);
mesh.write_csv(std::path::Path::new("femur.csv"))?;
```

- Threshold: 200 HU default (literature band 130-300; tune per protocol).
- One hex element per solid voxel; shared corner nodes are compacted.
- Per-element modulus from HU → density → Morgan-Keaveny-style power law
  (cortical above 1.3 g/cm³ apparent density, trabecular below).
- CLI: `cargo run -p tpt-med-examples --bin dicom-to-mesh -- <dir> --output m.csv --threshold 200 --smooth 5`

Compressed transfer syntaxes are rejected, not mis-parsed — decompress at
the archive boundary (see `rfcs/0001-dicom-ingestion.md`).
