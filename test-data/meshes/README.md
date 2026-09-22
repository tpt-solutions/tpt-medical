# Mesh exports

`synthetic_femur.csv` is the Phase-1 milestone output: the synthetic CT
series segmented at 200 HU, voxel-to-hex meshed, 3 Laplacian smoothing
iterations:

```console
cargo run -p tpt-med-examples --bin dicom-to-mesh -- \
    test-data/dicom/synthetic_ct --output test-data/meshes/synthetic_femur.csv \
    --threshold 200 --smooth 5
```

Format is documented in `crates/imaging/tpt-med-meshing/src/csv.rs`.
