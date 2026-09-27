//! Nonlinear hyperelastic internal-force and tangent-stiffness assembly for
//! [`Hex8Mesh`] elements.
//!
//! This is the 3-D `Hex8` assembly the TPT substrate does not ship at 0.1.0
//! (`rfcs/0009-nonlinear-fem-substrate-adapter.md`, "actual gap" item 1). It is
//! built entirely from substrate primitives that do exist — `tpt_fem_element`
//! shape functions and gradients, `tpt_fem_quadrature` tensor-product rules via
//! `tpt_fem_element::hex_rule`, and `tpt_fem_sparse::Coo` for the global
//! assembly — and takes its per-quadrature-point stress from the in-house
//! [`TissueModel::first_piola`].
//!
//! # Formulation
//!
//! Total Lagrangian, penalty incompressibility (the volumetric term lives
//! inside the tissue model, so no pressure unknown is needed — the "penalty
//! formulation recommended for this first increment" of RFC 0009). With
//! `G = dN/dX` the gradients of the *reference* configuration,
//! `F = I + u_{i,I} G_{j,I}` and `P = P(F)`:
//!
//! ```text
//! f_int(u)[k, I] = int  sum_L  P(F)[k, L] G[I, L] dV
//! K[(k,I),(m,J)]  = int  sum_{L,M} A(F)[k,L,m,M] G[I,L] G[J,M] dV  (k == m)
//! ```
//!
//! where `A = dP/dF` is the material tangent. The first line is
//! `dE/du` for `E = int W(F) dV`; the second is `d²E/du²`, which is exactly
//! `B^T A B` written in index form — no Voigt matrix, so no engineering-shear
//! weighting can be got wrong, and no initial-stress ("geometric") term,
//! which belongs to the updated-Lagrangian formulation, not this one. Being
//! the true Hessian of the discrete energy is also why the assembled tangent
//! is symmetric, which the test suite checks.
//!
//! # Tangent strategy
//!
//! RFC 0009 deliberately left "analytic vs. numerical tangent" unresolved.
//! Both are implemented here so the choice can be made on numbers:
//!
//! - [`tangent_stiffness`] differentiates the *constitutive law* only
//!   (`A = dP/dF` by central differences, 18 evaluations per quadrature
//!   point), then contracts it with the shape-function gradients. One
//!   assembly per Newton step.
//! - [`tangent_stiffness_numerical`] differentiates the *whole residual* with
//!   respect to every DOF, so it is consistent by construction but costs two
//!   full assemblies per DOF per Newton step.
//!
//! The crate's verification suite checks the two against each other (and
//! against minor symmetry) rather than asserting one in the abstract.

use crate::mesh::{HexMesh, MeshError};
use tpt_fem_element::hex_rule;
use tpt_fem_element::ReferenceElement;
use tpt_fem_sparse::Coo;
use tpt_med_geometry::Mat3;
use tpt_med_tissue::TissueModel;

/// A fourth-order tensor `T[i][j][k][l]`.
pub type Tensor4 = [[[[f64; 3]; 3]; 3]; 3];

/// Anything that can return a first Piola-Kirchhoff stress for a deformation
/// gradient.
///
/// Implemented for [`TissueModel`] (the in-house constitutive models) and for
/// [`FnModel`], which wraps any `Fn(&Mat3) -> Mat3` so callers are not forced
/// to introduce a named type per material.
pub trait Constitutive {
    /// First Piola-Kirchhoff stress `P = dW/dF`.
    fn first_piola(&self, f: &Mat3) -> Mat3;
}

impl Constitutive for TissueModel {
    fn first_piola(&self, f: &Mat3) -> Mat3 {
        TissueModel::first_piola(self, f)
    }
}

/// Adapter letting a bare closure act as a [`Constitutive`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FnModel<F>(pub F);

impl<F: Fn(&Mat3) -> Mat3> Constitutive for FnModel<F> {
    fn first_piola(&self, f: &Mat3) -> Mat3 {
        (self.0)(f)
    }
}

/// Quadrature order and finite-difference step for the assembly routines.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AssemblyOptions {
    /// Tensor-product Gauss order per axis. `2` is `2x2x2` (the order the
    /// in-house voxel core uses); `3` is the default because the
    /// finite-difference material tangent is not a low-order polynomial in
    /// `xi` and order 2 under-integrates it.
    pub quadrature_order: usize,
    /// Central-difference step for the material tangent.
    pub fd_step: f64,
}

impl Default for AssemblyOptions {
    fn default() -> Self {
        Self {
            quadrature_order: 3,
            fd_step: 1.0e-7,
        }
    }
}

impl AssemblyOptions {
    /// Options with the given quadrature order and the default `fd_step`.
    pub fn with_quadrature_order(quadrature_order: usize) -> Self {
        Self {
            quadrature_order,
            ..Self::default()
        }
    }
}

/// Deformation gradient `F` of element `element` at reference point `xi`.
///
/// Returns `None` when the element is inverted or degenerate there (the same
/// condition under which the physical gradients are undefined).
pub fn element_deformation_gradient<E: ReferenceElement>(
    mesh: &HexMesh<E>,
    element: usize,
    u: &[f64],
    xi: &[f64; 3],
) -> Option<Mat3> {
    let grad = mesh.physical_gradients(element, xi)?;
    let mut f = Mat3::IDENTITY;
    for (local, &node) in mesh.elements()[element].iter().enumerate() {
        for k in 0..3 {
            for l in 0..3 {
                let v = f.at(k, l) + u[3 * node + k] * grad[local][l];
                f.set(k, l, v);
            }
        }
    }
    Some(f)
}

/// Material tangent `A[i][j][k][l] = dP_ij/dF_kl` by central differences.
///
/// Eighteen perturbed evaluations of the constitutive law are enough regardless
/// of which model is supplied, which is why this adapter can serve every
/// [`TissueModel`] variant — including the Ogden and HGO variants whose
/// analytic first Piola is itself a finite-difference approximation.
///
/// The step is taken as `h * max(1, |F_kl|)` so the derivative is equally
/// accurate on a unit-scale configuration and on a strongly stretched one; the
/// bare `h` would lose half its digits in the second case.
pub fn material_tangent(model: &dyn Constitutive, f: &Mat3, h: f64) -> Tensor4 {
    let mut out = [[[[0.0f64; 3]; 3]; 3]; 3];
    for k in 0..3 {
        for l in 0..3 {
            let step = h * f.at(k, l).abs().max(1.0);
            let mut fp = *f;
            let plus = fp.at(k, l) + step;
            fp.set(k, l, plus);
            let mut fm = *f;
            let minus = fm.at(k, l) - step;
            fm.set(k, l, minus);
            let pp = model.first_piola(&fp);
            let pm = model.first_piola(&fm);
            for i in 0..3 {
                for j in 0..3 {
                    out[i][j][k][l] = (pp.at(i, j) - pm.at(i, j)) / (plus - minus);
                }
            }
        }
    }
    out
}

/// Assembles the internal force vector `f_int(u) = int B^T P dV`.
///
/// Elements that are inverted or degenerate at an integration point contribute
/// nothing there (their `J^-T` is undefined), which is why
/// [`Hex8Mesh::inverted_elements`] exists: a caller driving a large
/// deformation should check it rather than silently assemble a partial load.
///
/// # Errors
///
/// [`MeshError::DofCountMismatch`] if `u` is not `3 * node_count` long.
pub fn internal_force<E: ReferenceElement>(
    mesh: &HexMesh<E>,
    model: &dyn Constitutive,
    u: &[f64],
    opts: &AssemblyOptions,
) -> Result<Vec<f64>, MeshError> {
    if u.len() != mesh.dof_count() {
        return Err(MeshError::DofCountMismatch {
            expected: mesh.dof_count(),
            found: u.len(),
        });
    }
    let mut force = vec![0.0; mesh.dof_count()];
    let rule = hex_rule(opts.quadrature_order);
    for e in 0..mesh.element_count() {
        for (xi, w) in rule.points.iter().zip(&rule.weights) {
            let Some(grad) = mesh.physical_gradients(e, xi) else {
                return Err(MeshError::DegenerateElement {
                    element: e,
                    jacobian_determinant: mesh.jacobian(e, xi).det(),
                });
            };
            let dv = w * mesh.jacobian(e, xi).det();
            let f = element_deformation_gradient(mesh, e, u, xi).ok_or(
                MeshError::DegenerateElement {
                    element: e,
                    jacobian_determinant: mesh.jacobian(e, xi).det(),
                },
            )?;
            let j = f.det();
            if !j.is_finite() || j <= 0.0 {
                // An inverted element has no valid deformation gradient for the
                // constitutive law (the in-house models evaluate J^-2/3, which is
                // NaN for J <= 0). Failing loudly here beats assembling a NaN
                // matrix and having the linear solver report it as singular.
                return Err(MeshError::InvertedDeformation {
                    element: e,
                    deformation_determinant: j,
                });
            }
            let p = model.first_piola(&f);
            for (local, &node) in mesh.elements()[e].iter().enumerate() {
                for k in 0..3 {
                    let mut acc = 0.0;
                    for l in 0..3 {
                        acc += p.at(k, l) * grad[local][l];
                    }
                    force[3 * node + k] += dv * acc;
                }
            }
        }
    }
    Ok(force)
}

/// Assembles the tangent stiffness `B^T A B = d²E/du²` as a global `Coo`.
///
/// Costs one constitutive evaluation per quadrature point plus the material
/// tangent's 18. See [`tangent_stiffness_numerical`] for the residual-level
/// alternative.
///
/// # Errors
///
/// [`MeshError::DofCountMismatch`] if `u` is not `3 * node_count` long.
pub fn tangent_stiffness<E: ReferenceElement>(
    mesh: &HexMesh<E>,
    model: &dyn Constitutive,
    u: &[f64],
    opts: &AssemblyOptions,
) -> Result<Coo, MeshError> {
    if u.len() != mesh.dof_count() {
        return Err(MeshError::DofCountMismatch {
            expected: mesh.dof_count(),
            found: u.len(),
        });
    }
    let mut stiffness = Coo::with_capacity(mesh.element_count() * 24 * 24);
    let rule = hex_rule(opts.quadrature_order);
    for e in 0..mesh.element_count() {
        for (xi, w) in rule.points.iter().zip(&rule.weights) {
            let Some(grad) = mesh.physical_gradients(e, xi) else {
                return Err(MeshError::DegenerateElement {
                    element: e,
                    jacobian_determinant: mesh.jacobian(e, xi).det(),
                });
            };
            let dv = w * mesh.jacobian(e, xi).det();
            let f = element_deformation_gradient(mesh, e, u, xi).ok_or(
                MeshError::DegenerateElement {
                    element: e,
                    jacobian_determinant: mesh.jacobian(e, xi).det(),
                },
            )?;
            let a = material_tangent(model, &f, opts.fd_step);
            for (i, &ni) in mesh.elements()[e].iter().enumerate() {
                for (j, &nj) in mesh.elements()[e].iter().enumerate() {
                    for k in 0..3 {
                        for m in 0..3 {
                            // No Kronecker delta here: the material tangent of a
                            // general hyperelastic law couples the two material
                            // index pairs, so cross-component stiffness terms are
                            // real and required (they vanish only for a
                            // diagonal-`A` material).
                            let mut acc = 0.0;
                            for l in 0..3 {
                                for n in 0..3 {
                                    acc += a[k][l][m][n] * grad[i][l] * grad[j][n];
                                }
                            }
                            stiffness.push(3 * ni + k, 3 * nj + m, dv * acc);
                        }
                    }
                }
            }
        }
    }
    Ok(stiffness)
}

/// Assembles the tangent stiffness by central-differencing the *whole residual*
/// with respect to every degree of freedom.
///
/// This is the definitionally consistent tangent (`df_int/du`, not a
/// hand-derived contraction of the material law), so it is the reference the
/// cheaper [`tangent_stiffness`] is checked against. It costs two full
/// internal-force assemblies per DOF per Newton step — quadratic in problem
/// size, and only appropriate for the small verification problems this crate
/// ships, not for a production solve.
///
/// # Errors
///
/// [`MeshError::DofCountMismatch`] if `u` is not `3 * node_count` long.
pub fn tangent_stiffness_numerical<E: ReferenceElement>(
    mesh: &HexMesh<E>,
    model: &dyn Constitutive,
    u: &[f64],
    opts: &AssemblyOptions,
) -> Result<Coo, MeshError> {
    let n = mesh.dof_count();
    let mut stiffness = Coo::with_capacity(n * n);
    let mut probe = u.to_vec();
    let step = opts.fd_step;
    for b in 0..n {
        probe.copy_from_slice(u);
        probe[b] = u[b] + step;
        let plus = internal_force(mesh, model, &probe, opts)?;
        probe[b] = u[b] - step;
        let minus = internal_force(mesh, model, &probe, opts)?;
        probe[b] = u[b];
        for a in 0..n {
            stiffness.push(a, b, (plus[a] - minus[a]) / (2.0 * step));
        }
    }
    Ok(stiffness)
}

/// Maximum absolute difference between two `Coo` matrices over the union of
/// their `(row, col)` supports, with entries missing from one side treated as
/// zero.
///
/// Used by this crate's own verification of the two tangent strategies against
/// each other.
pub fn coo_max_abs_diff(a: &Coo, b: &Coo) -> f64 {
    let lookup = |m: &Coo, row: usize, col: usize| -> f64 {
        m.rows
            .iter()
            .zip(&m.cols)
            .zip(&m.vals)
            .filter(|((&r, &c), _)| r == row && c == col)
            .map(|(_, &v)| v)
            .sum()
    };
    let mut worst: f64 = 0.0;
    let mut seen = std::collections::HashSet::new();
    for ((&r, &c), _) in a.rows.iter().zip(&a.cols).zip(&a.vals) {
        if seen.insert((r, c)) {
            worst = worst.max((lookup(a, r, c) - lookup(b, r, c)).abs());
        }
    }
    for ((&r, &c), _) in b.rows.iter().zip(&b.cols).zip(&b.vals) {
        if seen.insert((r, c)) {
            worst = worst.max((lookup(a, r, c) - lookup(b, r, c)).abs());
        }
    }
    worst
}
