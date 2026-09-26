# tpt-med-biomechanics

Voxel hexahedral FEM solver core for patient-specific biomechanics — trilinear
Q1 hexes, 2×2×2 Gauss quadrature, CSR assembly, Jacobi-preconditioned
conjugate gradient, and von Mises post-processing.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--biomechanics-orange)](https://crates.io/crates/tpt-med-biomechanics)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--biomechanics-blue)](https://docs.rs/tpt-med-biomechanics)

| | |
|---|---|
| **Layer** | `solid` |
| **Status** | Alpha, `0.1.0` |
| **Scope** | Linear (small-strain) isotropic elasticity |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | [`tpt-med-meshing`](../../imaging/tpt-med-meshing), [`tpt-med-geometry`](../../core/tpt-med-geometry), [`tpt-med-units`](../../core/tpt-med-units) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

The question this crate answers is narrow and deliberately so: *given a
patient-specific voxel mesh with HU-derived moduli and a set of loads, where is
the stress?*

Abaqus answers that question, but it is proprietary, script-hostile, and
cannot run in a browser. This crate answers it in ~600 lines of dependency-free
Rust that compiles to WASM in under a second, is verified against closed-form
analytical solutions, and whose every intermediate number is inspectable.

**Scope is the point.** This is a *linear* core: 8-node trilinear hexes,
isotropic materials, static loading, conjugate gradient. It is not a
general-purpose nonlinear FEM package, and it is not trying to be. Nonlinear
hyperelasticity, large deformation and contact are the documented upgrade path
via `tpt-fem` / `tpt-fem-hyperelastic` / `tpt-fem-contact` (RFC 0002), pinned in
the workspace manifest as the sanctioned integration points.

## Features

- **8-node trilinear hexahedra (Q1)** with 2×2×2 Gauss quadrature — the
  natural match for the voxel meshes `tpt-med-meshing` produces, so the mesh
  and the element are the same object and there is no surface-fitting step.
- **Analytic element stiffness** — `trilinear_hex_stiffness` returns the full
  24×24 matrix; `isotropic_d` builds the 6×6 constitutive matrix.
- **CSR sparse assembly** — `CsrMatrix::from_triplets` with duplicate
  summation, plus `mul_vec` and `diagonal` for preconditioning.
- **Jacobi-preconditioned conjugate gradient** with iteration count and final
  residual reported back (`SolveStats`).
- **Boundary conditions** — `fix_nodes` (all 3 DOFs), `add_force` (nodal, in
  N), and prescribed nodal displacements.
- **Post-processing** — per-element strain/stress in Voigt notation, von Mises,
  principal stresses (analytic symmetric eigenvalues via `tpt-med-geometry`),
  hydrostatic stress, and the critical element index.
- **Model diagnostics** — `rigid_mode_residual()` catches unconstrained
  rigid-body motion. `coupling()` is a placeholder for a second diagnostic
  (detecting a mesh that loads in a way that decouples into a floppy mode)
  and is not yet implemented — see Known Limitations.
- **Input validation** — empty meshes, material-array length mismatches and
  degenerate (zero/inverted-Jacobian) elements are `SolverError::Invalid`,
  never a panic or a `NaN`.

## Conventions

- Displacement-based, small-strain formulation; nodal unknowns, 3 DOFs each.
- Units: **mm**, **N**, **MPa**, **tonnes** (a consistent mm–N–MPa set).
- Stiffness is 24×24 (8 nodes × 3 DOFs) per element; stress/strain are Voigt
  `[xx, yy, zz, xy, yz, xz]` with engineering shear strains.
- Element node ordering is the `tpt-med-meshing` hex convention:
  `[000, 100, 110, 010, 001, 101, 111, 011]`.

## Usage

This crate takes a mesh as plain nodes/elements (typically built by
`tpt-med-meshing` from a DICOM-derived `VoxelHexMesh` via
`BiomechanicsModel::from_voxel_mesh`); the example below builds one directly
with `BiomechanicsModel::from_parts` so it depends on nothing beyond
`tpt-med-geometry` and `tpt-med-units`:

```rust
use tpt_med_biomechanics::{BiomechanicsModel, BoundaryConditions};
use tpt_med_geometry::Vec3;
use tpt_med_units::Force;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // A structured 4x4x10 mm grid of 1 mm cubes, uniform E = 15,000 MPa,
    // nu = 0.3 (representative cortical bone), in place of a DICOM-derived mesh.
    let (nx, ny, nz) = (4usize, 4usize, 10usize);
    let mut nodes = Vec::new();
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push(Vec3::new(i as f64, j as f64, k as f64));
            }
        }
    }
    let id = |i: usize, j: usize, k: usize| (k * (ny + 1) + j) * (nx + 1) + i;
    let mut elements = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                elements.push([
                    id(i, j, k) as u32, id(i + 1, j, k) as u32,
                    id(i + 1, j + 1, k) as u32, id(i, j + 1, k) as u32,
                    id(i, j, k + 1) as u32, id(i + 1, j, k + 1) as u32,
                    id(i + 1, j + 1, k + 1) as u32, id(i, j + 1, k + 1) as u32,
                ]);
            }
        }
    }
    let model = BiomechanicsModel::from_parts(nodes.clone(), elements, 15_000.0, 0.3);

    // Distal face fixed, proximal face loaded (ISO 7206-style stance).
    let mut bc = BoundaryConditions::default();
    let mut fixed = Vec::new();
    let mut loaded = Vec::new();
    let min_z = nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
    let max_z = nodes.iter().map(|n| n.z).fold(f64::NEG_INFINITY, f64::max);
    for (i, n) in nodes.iter().enumerate() {
        if n.z <= min_z + 1e-6 { fixed.push(i as u32); }
        if n.z >= max_z - 1e-6 { loaded.push(i as u32); }
    }
    bc.fix_nodes(fixed);

    // 3x body weight (80 kg patient) distributed over the loaded face.
    let stance = Force::body_weights(3.0, Force::from_n(80.0 * 9.81));
    let per_node = -stance.to_n() / loaded.len().max(1) as f64;
    for i in loaded { bc.add_force(i, Vec3::new(0.0, 0.0, per_node)); }

    let result = model.solve(&bc, 1e-8, 5_000)?;

    println!("max von Mises   {:.2} MPa", result.max_von_mises());
    println!("mean von Mises  {:.2} MPa", result.mean_von_mises());
    println!("max displacement {:.3} mm", result.max_displacement());
    println!("CG iterations: {}", result.stats.iterations);

    // Yield screen against cortical bone (~110 MPa).
    if let Some(id) = result.critical_element() {
        let e = &result.stresses[id as usize];
        println!("critical element {id}: von Mises {:.1}, principal {:?}",
                 e.von_mises(), e.principal_stresses());
    }
    Ok(())
}
```

## API Overview

| Item | Purpose |
|---|---|
| `BiomechanicsModel::from_voxel_mesh(&VoxelHexMesh)` | Adapt a voxel mesh, keeping per-element E and ν |
| `BiomechanicsModel::from_parts(nodes, elements, modulus, poisson)` | Build from raw parts (synthetic beams, benchmark meshes) |
| `BiomechanicsModel::solve(&BoundaryConditions, tolerance, max_iterations)` | Assemble, solve, post-process → `Result<StressResult, SolverError>` |
| `BiomechanicsModel::element_centroid(id)` | Element centre position (mm) |
| `BiomechanicsModel::rigid_mode_residual()` | Diagnostic: residual rigid-body motion; ~0 for a well-constrained model |
| `BiomechanicsModel::coupling()` | Placeholder diagnostic; currently a stub returning `Mat3::IDENTITY` (see Known Limitations) |
| `BoundaryConditions::fix_nodes(iter)` | Fully constrain nodes (all 3 DOFs) |
| `BoundaryConditions::add_force(node, Vec3)` | Nodal force in N; duplicates are summed |
| `BoundaryConditions::prescribed_displacements` | Prescribed nodal displacements (mm), applied after `fix_nodes` |
| `StressResult` | `displacements: Vec<Vec3>`, `stats: SolveStats`, `stresses: Vec<ElementStress>` |
| `StressResult::{max_von_mises, mean_von_mises, max_displacement, critical_element}` | Headline results |
| `ElementStress` | `element`, `strain: [f64; 6]`, `stress: [f64; 6]` in Voigt notation (MPa) |
| `ElementStress::{von_mises, principal_stresses, hydrostatic}` | Derived measures; principal stresses sorted descending |
| `hex::trilinear_hex_stiffness(&[Vec3; 8], e, nu)` (re-exported at crate root) | 24×24 element stiffness matrix |
| `hex::isotropic_d(e, nu)` | 6×6 isotropic constitutive matrix |
| `hex::centre_strain_displacement(&[Vec3; 8])` | Centroidal strain/displacement pair, used for patch tests |
| `hex::check_element(&[Vec3; 8], id, error_sink)` | Geometry validation hook used during solve |
| `sparse::CsrMatrix::from_triplets(n, triplets)` | Build from `(row, col, value)` triplets, summing duplicates |
| `sparse::CsrMatrix::{mul_vec, diagonal}` | Spmv and Jacobi preconditioner extraction |
| `sparse::conjugate_gradient(...)` | Jacobi-preconditioned CG returning `SolveStats` |
| `SolverError` | `Invalid(String)` — empty mesh, length mismatch, degenerate element, singular system |

Only `trilinear_hex_stiffness`, `ElementStress`, `StressResult`,
`BiomechanicsModel`, `BoundaryConditions` and `SolverError` are re-exported at
the crate root (`tpt_med_biomechanics::`); `isotropic_d`,
`centre_strain_displacement`, `check_element` and `CsrMatrix` /
`conjugate_gradient` live under the `hex::` and `sparse::` module paths shown
above.

## Verification

Following ASME V&V 40, this crate ships **code and calculation verification**
against closed-form analytical solutions, not regression snapshots:

- **Uniaxial tension** on a single hex and on a block — stress matches `E·ε`
  exactly, and transverse contraction matches `ν`.
- **Cantilever beam** — tip deflection compared against the Euler–Bernoulli
  analytic solution `δ = FL³/(3EI)`.
- **Patch test** — a uniform-strain patch reproduces constant stress, the
  standard check that the shape functions and quadrature are correct.
- **Rigid-mode residual** is asserted negligible on a properly constrained
  model, catching under-constrained BCs before they produce plausible
  nonsense.
- Degenerate elements (inverted or zero Jacobian) are asserted to produce
  `SolverError::Invalid`, never `NaN`s.
- Golden reference datasets `test-data/golden/solid/femur_loading.json` and
  `lumbar_spine_compression.json`, with documented analytical basis.

The workspace [V&V policy](https://github.com/tpt-solutions/tpt-medical/blob/master/docs/book/src/vv-policy.md)
records which claims are verified, which are validated, and which are neither.

## Performance

On the benchmark meshes in `benches/` (`hyperelastic-fem.rs`), the linear
solve is dominated by CSR assembly and the CG iterations. The matrix is SPD
and well-conditioned for a well-constrained bone model, so Jacobi
preconditioning plus CG converges in tens to low hundreds of iterations on
realistic CT-sized meshes. Single-precision node data would halve the memory
footprint but is not used: stress accuracy in MPa is the product requirement,
and it is worth more than the bytes.

## Known Limitations

- **Linear small-strain only.** This is the central limitation and the reason
  the crate is small. No hyperelasticity, no large deformation, no geometric
  stiffness, no plasticity, no creep. A femoral head under physiological load
  deforms roughly elastically, which is why this is defensible for bone
  screening — but it is not a general FEM package and will be wrong for
  anything soft.
- **Q1 hexes are stiff.** Trilinear hexahedra lock volumetrically and
  over-predict bending stiffness relative to a quadratic element. For a
  bending-dominated structure this biases stress downward near the neutral
  axis. A Q8 or Q20 element would fix it, at real cost in DOFs and assembly
  time.
- **No contact.** Two bodies in a TKA planning workflow cannot touch through
  this solver; contact is a `tpt-fem-contact` responsibility.
- **Isotropic materials per element.** `element_modulus` and `element_poisson`
  are scalars, so the transversely isotropic and orthotropic descriptions in
  `tpt-med-bone` cannot be expressed to the solver. A directional modulus
  would need a full tensor per element and a correspondingly richer element
  formulation.
- **Single-threaded, and the CG is the bottleneck.** Jacobi-preconditioned CG
  is fine for a few hundred thousand DOFs but will not scale to a whole-body
  mesh without a multigrid or a direct solver. No sparse direct factorisation
  is available.
- **No dynamic or transient analysis.** `solve` is static. There is no
  transient, frequency-domain or modal capability.
- **No multi-constraint boundary conditions.** `fix_nodes` removes all 3 DOFs;
  symmetry planes and roller constraints must be emulated by the caller
  applying reaction forces, which is more work and more error-prone than it
  should be.
- **Conjugate gradient can fail to converge within `max_iterations`** and the
  `SolveStats` residual must be inspected. There is no adaptive restart and no
  fallback to a direct solve, so a badly conditioned model returns a
  non-converged result that looks otherwise valid.
- **Verification is against linear analytical solutions**, which is the correct
  standard for a linear code but says nothing about a nonlinear one that will
  eventually be built on top of this.
- **`coupling()` is a stub.** It currently returns `Mat3::IDENTITY`
  unconditionally and does not perform the floppy-mode diagnostic described in
  Features; only `rigid_mode_residual()` is a real diagnostic today.

## Related Crates

- [`tpt-med-meshing`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-meshing) — produces the `VoxelHexMesh`; the hex corner order is the contract between the two crates.
- [`tpt-med-dicom`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/imaging/tpt-med-dicom) — the HU→E correlations behind the per-element materials.
- [`tpt-med-tissue`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-tissue) — hyperelastic models for the nonlinear upgrade path (RFC 0002).
- [`tpt-med-bone`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/solid/tpt-med-bone) — anisotropic bone materials and remodeling.
- [`tpt-med-wasm`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-wasm) — runs this solver in the browser.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). New solver
algorithms require an [RFC](../../../rfcs). Numerical code must ship
verification against an analytical or published reference.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.
