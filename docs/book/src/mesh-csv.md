# Mesh CSV format

Versioned in the header (`tpt-medical voxel hex mesh v1`):

```text
# tpt-medical voxel hex mesh v1
# units: mm, MPa
nodes,<count>
node,<id>,<x>,<y>,<z>
elements,<count>
hex,<id>,<n0>,...,<n7>,<youngs_modulus_mpa>,<density_gcm3>
```

Node ordering in a hex follows `(±x, ±y, ±z)` corner bit order
(`[000,100,110,010,001,101,111,011]`). The format round-trips through
`VoxelHexMesh::parse_csv` (mean HU and Poisson ratio are not exported).
