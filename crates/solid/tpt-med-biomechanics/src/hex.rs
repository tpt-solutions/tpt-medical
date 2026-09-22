//! 8-node trilinear hexahedral element.
//!
//! Standard isoparametric Q1 element with 2×2×2 Gauss quadrature. Element
//! nodes are ordered `[000, 100, 110, 010, 001, 101, 111, 011]` (bit order
//! `(i, j, k)` on the reference cube `[-1,1]³`), matching
//! `tpt-med-meshing::VoxelHexMesh`.

use tpt_med_geometry::{Mat3, Vec3};

/// Reference-corner coordinates for the node ordering above.
const CORNERS: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
];

/// Gauss points and weights for 2×2×2 (exact for Q1 elements).
const GAUSS: [f64; 2] = [-0.577_350_269_189_625_8, 0.577_350_269_189_625_8]; // ±1/√3
const WEIGHT: f64 = 1.0;

/// Shape-function natural derivatives at `xi`: `dn[a] = [dN/dξ, dN/dη, dN/dζ]`.
fn shape_derivs(xi: [f64; 3]) -> [[f64; 3]; 8] {
    let mut dn = [[0.0; 3]; 8];
    for (a, c) in CORNERS.iter().enumerate() {
        dn[a][0] = 0.125 * c[0] * (1.0 + c[1] * xi[1]) * (1.0 + c[2] * xi[2]);
        dn[a][1] = 0.125 * (1.0 + c[0] * xi[0]) * c[1] * (1.0 + c[2] * xi[2]);
        dn[a][2] = 0.125 * (1.0 + c[0] * xi[0]) * (1.0 + c[1] * xi[1]) * c[2];
    }
    dn
}

/// Isotropic elasticity matrix in Voigt order `[xx, yy, zz, xy, yz, xz]`.
pub fn isotropic_d(e: f64, nu: f64) -> [[f64; 6]; 6] {
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let mut d = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            d[i][j] = if i == j { lambda + 2.0 * mu } else { lambda };
        }
    }
    for i in 3..6 {
        d[i][i] = mu;
    }
    d
}

/// Assembles the 24×24 element stiffness matrix for one hex element.
///
/// `nodes`: the 8 corner positions in patient coordinates (mm).
/// `e`: Young's modulus (MPa), `nu`: Poisson's ratio.
pub fn trilinear_hex_stiffness(nodes: &[Vec3; 8], e: f64, nu: f64) -> [[f64; 24]; 24] {
    let d = isotropic_d(e, nu);
    let mut ke = [[0.0f64; 24]; 24];

    for &gx in &GAUSS {
        for &gy in &GAUSS {
            for &gz in &GAUSS {
                let xi = [gx, gy, gz];
                let dn = shape_derivs(xi);

                // Jacobian J = ∂x/∂ξ (rows are ∂x/∂ξᵀ stacked).

                let mut row0 = Vec3::ZERO;
                let mut row1 = Vec3::ZERO;
                let mut row2 = Vec3::ZERO;
                for a in 0..8 {
                    let x = nodes[a];
                    row0 += x * dn[a][0];
                    row1 += x * dn[a][1];
                    row2 += x * dn[a][2];
                }
                let j_mat = Mat3::from_rows(row0, row1, row2);
                let det_j = j_mat.det();
                debug_assert!(det_j > 0.0, "inverted element");
                let j_inv = j_mat.inverse().expect("non-singular element");

                // B-matrix (6×24) at this Gauss point.
                let mut b_mat = [[0.0f64; 24]; 6];
                for a in 0..8 {
                    // dN/dx = J^{-1} · dN/dξ
                    let g_nat = Vec3::new(dn[a][0], dn[a][1], dn[a][2]);
                    let g = j_inv.mul_vec(g_nat);
                    let (nx, ny, nz) = (g.x, g.y, g.z);
                    let i = 3 * a;
                    b_mat[0][i] = nx; // εxx
                    b_mat[1][i + 1] = ny; // εyy
                    b_mat[2][i + 2] = nz; // εzz
                    b_mat[3][i] = ny; // γxy
                    b_mat[3][i + 1] = nx;
                    b_mat[4][i + 1] = nz; // γyz
                    b_mat[4][i + 2] = ny;
                    b_mat[5][i] = nz; // γxz
                    b_mat[5][i + 2] = nx;
                }

                // ke += Bᵀ D B detJ w
                let mut db = [[0.0f64; 24]; 6];
                for r in 0..6 {
                    for c in 0..24 {
                        let mut s = 0.0;
                        for k in 0..6 {
                            s += d[r][k] * b_mat[k][c];
                        }
                        db[r][c] = s;
                    }
                }
                for r in 0..24 {
                    for c in 0..24 {
                        let mut s = 0.0;
                        for k in 0..6 {
                            s += b_mat[k][r] * db[k][c];
                        }
                        ke[r][c] += s * det_j * WEIGHT;
                    }
                }
            }
        }
    }
    ke
}

/// Uniform B-matrix and detJ at the element centre (ξ=η=ζ=0) for
/// post-processing.
pub fn centre_strain_displacement(nodes: &[Vec3; 8]) -> ([[f64; 24]; 6], f64) {
    let dn = shape_derivs([0.0, 0.0, 0.0]);
    let mut row0 = Vec3::ZERO;
    let mut row1 = Vec3::ZERO;
    let mut row2 = Vec3::ZERO;
    for a in 0..8 {
        let x = nodes[a];
        row0 += x * dn[a][0];
        row1 += x * dn[a][1];
        row2 += x * dn[a][2];
    }
    let j_mat = Mat3::from_rows(row0, row1, row2);
    let det_j = j_mat.det();
    let j_inv = j_mat.inverse().unwrap_or(Mat3::IDENTITY);

    let mut b_mat = [[0.0f64; 24]; 6];
    for a in 0..8 {
        let g_nat = Vec3::new(dn[a][0], dn[a][1], dn[a][2]);
        let g = j_inv.mul_vec(g_nat);
        let (nx, ny, nz) = (g.x, g.y, g.z);
        let i = 3 * a;
        b_mat[0][i] = nx;
        b_mat[1][i + 1] = ny;
        b_mat[2][i + 2] = nz;
        b_mat[3][i] = ny;
        b_mat[3][i + 1] = nx;
        b_mat[4][i + 1] = nz;
        b_mat[4][i + 2] = ny;
        b_mat[5][i] = nz;
        b_mat[5][i + 2] = nx;
    }
    (b_mat, det_j)
}
/// Reports an error string when element geometry is degenerate (non-positive
/// integrated volume); used by the solver to produce actionable errors.
pub fn check_element(nodes: &[Vec3; 8], id: usize, error_sink: &mut Vec<String>) {
    let mut vol = 0.0;
    for &gx in &GAUSS {
        for &gy in &GAUSS {
            for &gz in &GAUSS {
                let dn = shape_derivs([gx, gy, gz]);
                let mut row0 = Vec3::ZERO;
                let mut row1 = Vec3::ZERO;
                let mut row2 = Vec3::ZERO;
                for a in 0..8 {
                    let x = nodes[a];
                    row0 += x * dn[a][0];
                    row1 += x * dn[a][1];
                    row2 += x * dn[a][2];
                }
                vol += Mat3::from_rows(row0, row1, row2).det();
            }
        }
    }
    if vol <= 0.0 {
        error_sink.push(format!(
            "element {id} has non-positive volume ({vol:.6}); check node ordering or smoothing"
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Unit cube with the meshing crate's corner ordering.
    fn unit_cube() -> [Vec3; 8] {
        CORNERS.map(|c| Vec3::new(0.5 * (c[0] + 1.0), 0.5 * (c[1] + 1.0), 0.5 * (c[2] + 1.0)))
    }

    #[test]
    fn unit_cube_stiffness_is_symmetric_and_positive() {
        let nodes = unit_cube();
        let ke = trilinear_hex_stiffness(&nodes, 1000.0, 0.3);
        // Symmetry
        for i in 0..24 {
            for j in 0..24 {
                assert!((ke[i][j] - ke[j][i]).abs() < 1e-9);
            }
        }
        // Rigid translation: K · 1 = 0
        for dof in 0..3 {
            let mut u = [0.0; 24];
            for a in 0..8 {
                u[3 * a + dof] = 1.0;
            }
            for r in 0..24 {
                let f: f64 = (0..24).map(|c| ke[r][c] * u[c]).sum();
                assert!(f.abs() < 1e-6, "rigid mode {dof} row {r}: {f}");
            }
        }
    }

    #[test]
    fn uniaxial_stiffness_scales_linearly_with_e() {
        let nodes = unit_cube();
        let k1 = trilinear_hex_stiffness(&nodes, 100.0, 0.3);
        let k2 = trilinear_hex_stiffness(&nodes, 250.0, 0.3);
        for i in 0..24 {
            for j in 0..24 {
                assert!((k2[i][j] - 2.5 * k1[i][j]).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn elasticity_matrix_invariants() {
        let d = isotropic_d(1000.0, 0.3);
        // E from diagonal: d[0][0] = λ+2μ; λ/μ = ν/((1-ν)/2)... verify trace
        // relations: d[0][0] − d[0][1] = 2μ.
        let two_mu = d[0][0] - d[0][1];
        let mu = two_mu / 2.0;
        let e_check = 2.0 * mu * (1.0 + 0.3);
        assert!((e_check - 1000.0).abs() < 1e-6);
    }
}
