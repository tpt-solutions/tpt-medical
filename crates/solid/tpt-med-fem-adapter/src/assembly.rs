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

use crate::mesh::{Mesh, MeshError};
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
    /// The **volumetric part** of `P`, so a caller can integrate it separately
    /// from the deviatoric part.
    ///
    /// # Why this exists
    ///
    /// Selective reduced integration — the cheap mitigation for volumetric
    /// locking — works by integrating the volumetric response on a lower-order
    /// rule than the deviatoric one. That is only possible if the two can be
    /// told apart. `first_piola` returns them fused, and the split cannot be
    /// recovered from a single `Mat3` in general: under-integrating "all of `P`"
    /// would damage the deviatoric response, not just the volumetric one.
    ///
    /// It is also the first third of a mixed `u`-`p` formulation, which cannot
    /// exist until the split does — the pressure unknown *is* the volumetric
    /// part.
    ///
    /// # The default is a sharp edge, and it is left visible
    ///
    /// The default is zero, meaning "this law has no separable volumetric part".
    /// That is right for a genuinely incompressible-by-construction law and
    /// silently wrong for a bare closure that happens to include a penalty, so
    /// `FnModel` keeps this default: a user who wraps a penalty-incompressible
    /// closure in `FnModel` and switches on reduced volumetric integration gets
    /// **no effect and no warning**, because the trait cannot distinguish that
    /// closure from one with genuinely no volumetric term.
    ///
    /// That is why [`AssemblyOptions::volumetric_quadrature_order`] is opt-in and
    /// says so in its own docs. A `FnModel` user who wants the behaviour
    /// implements `Constitutive` directly and supplies the split.
    fn volumetric_piola(&self, _f: &Mat3) -> Mat3 {
        Mat3::ZERO
    }
    /// First Piola-Kirchhoff stress at quadrature point `point` of element
    /// `element`.
    ///
    /// The assembly routines call this rather than [`first_piola`] so a law
    /// with **per-point internal state** (the superelastic model's martensite
    /// fraction) can look its own state up. The default ignores the indices
    /// and delegates, so every stateless law is unaffected. `point` indexes
    /// the element family's rule at `AssemblyOptions::quadrature_order`, in
    /// rule order.
    ///
    /// [`first_piola`]: Constitutive::first_piola
    fn first_piola_at(&self, _element: usize, _point: usize, f: &Mat3) -> Mat3 {
        self.first_piola(f)
    }
}

impl Constitutive for TissueModel {
    fn first_piola(&self, f: &Mat3) -> Mat3 {
        TissueModel::first_piola(self, f)
    }

    fn volumetric_piola(&self, f: &Mat3) -> Mat3 {
        TissueModel::volumetric_first_piola(self, f)
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
    /// Quadrature order for the **volumetric part alone**, or `None` to
    /// integrate it on the same rule as everything else.
    ///
    /// `Some(order)` selects *selective reduced integration* (SRI), the cheap
    /// mitigation for volumetric locking: the deviatoric response keeps full
    /// order while the volumetric penalty is integrated on a coarser rule, which
    /// removes most of the spurious compressive stiffness a stiff penalty
    /// produces under full integration. It reduces locking; it does not eliminate
    /// it, and it is **not** exact incompressibility. A mixed `u`-`p`
    /// formulation is the only thing that is.
    ///
    /// # It only does something if the law reports a volumetric part
    ///
    /// See [`Constitutive::volumetric_piola`]. A law that returns the default
    /// zero here — notably a bare [`FnModel`] closure — makes this option a
    /// silent no-op. `None` (the default) is the current, fully-verified
    /// behaviour, so reduced integration is opt-in and the existing verification
    /// suite is unaffected by it.
    ///
    /// # The patch-test interaction
    ///
    /// The constant-stress patch identity is exact for the *deviatoric*
    /// response. Under SRI the volumetric part is deliberately under-integrated,
    /// so any test asserting an exact patch identity on a state with non-uniform
    /// volume change will see a small deviation. That is the point of the
    /// method, not a defect in it — the volumetric response is no longer
    /// integrated consistently, which is the entire mechanism.
    pub volumetric_quadrature_order: Option<usize>,
}

impl Default for AssemblyOptions {
    fn default() -> Self {
        Self {
            quadrature_order: 3,
            fd_step: 1.0e-7,
            volumetric_quadrature_order: None,
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

    /// Options with selective reduced integration of the volumetric part.
    ///
    /// See [`AssemblyOptions::volumetric_quadrature_order`] — including the
    /// caveat that it is a no-op for a law reporting no volumetric part.
    pub fn with_volumetric_quadrature_order(order: usize) -> Self {
        Self {
            volumetric_quadrature_order: Some(order),
            ..Self::default()
        }
    }
}

/// Deformation gradient `F` of element `element` at reference point `xi`.
///
/// Returns `None` when the element is inverted or degenerate there (the same
/// condition under which the physical gradients are undefined).
pub fn element_deformation_gradient<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
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
    material_tangent_of(&|g: &Mat3| model.first_piola(g), f, h)
}

/// [`material_tangent`] of an arbitrary stress function, so the assembler can
/// difference the point-indexed [`Constitutive::first_piola_at`].
fn material_tangent_of(stress: &dyn Fn(&Mat3) -> Mat3, f: &Mat3, h: f64) -> Tensor4 {
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
            let pp = stress(&fp);
            let pm = stress(&fm);
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
pub fn internal_force<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
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
    // The rule must come from the element family, not from a hexahedral
    // default: a tetrahedron lives on a different reference domain and a
    // tensor-product cube rule would evaluate its shape functions outside the
    // element, giving a plausible but badly wrong answer rather than an error.
    let rule = E::quadrature_rule(opts.quadrature_order);
    for e in 0..mesh.element_count() {
        for (q, (xi, w)) in rule.points.iter().zip(&rule.weights).enumerate() {
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
            // Selective reduced integration: the deviatoric part keeps the full
            // rule, the volumetric part is integrated on its own coarser rule.
            // With `None` both halves share one rule and this reduces to the
            // original single `P` integration, exactly.
            let p_full = model.first_piola_at(e, q, &f);
            let p_vol = match opts.volumetric_quadrature_order {
                Some(_) => model.volumetric_piola(&f),
                None => Mat3::ZERO,
            };
            let p = p_full - p_vol;
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
        if let Some(order) = opts.volumetric_quadrature_order {
            let vrule = E::quadrature_rule(order);
            for (xi, w) in vrule.points.iter().zip(&vrule.weights) {
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
                let p_vol = model.volumetric_piola(&f);
                for (local, &node) in mesh.elements()[e].iter().enumerate() {
                    for k in 0..3 {
                        let mut acc = 0.0;
                        for l in 0..3 {
                            acc += p_vol.at(k, l) * grad[local][l];
                        }
                        force[3 * node + k] += dv * acc;
                    }
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
pub fn tangent_stiffness<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
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
    // The rule must come from the element family, not from a hexahedral
    // default: a tetrahedron lives on a different reference domain and a
    // tensor-product cube rule would evaluate its shape functions outside the
    // element, giving a plausible but badly wrong answer rather than an error.
    let rule = E::quadrature_rule(opts.quadrature_order);
    for e in 0..mesh.element_count() {
        for (q, (xi, w)) in rule.points.iter().zip(&rule.weights).enumerate() {
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
            // The tangent must be split the same way the force is, or the
            // Jacobian stops being the derivative of the residual and Newton
            // degrades to a first-order crawl with no error. `a_vol` is
            // differenced from the volumetric law alone for the same reason the
            // force splits.
            let a = match opts.volumetric_quadrature_order {
                Some(_) => {
                    let full = material_tangent_of(
                        &|g: &Mat3| model.first_piola_at(e, q, g),
                        &f,
                        opts.fd_step,
                    );
                    let vol = material_tangent(&VolumetricOnly(model), &f, opts.fd_step);
                    sub4(&full, &vol)
                }
                None => {
                    material_tangent_of(&|g: &Mat3| model.first_piola_at(e, q, g), &f, opts.fd_step)
                }
            };
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
        if let Some(order) = opts.volumetric_quadrature_order {
            let vrule = E::quadrature_rule(order);
            for (xi, w) in vrule.points.iter().zip(&vrule.weights) {
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
                let a_vol = material_tangent(&VolumetricOnly(model), &f, opts.fd_step);
                for (i, &ni) in mesh.elements()[e].iter().enumerate() {
                    for (j, &nj) in mesh.elements()[e].iter().enumerate() {
                        for k in 0..3 {
                            for m in 0..3 {
                                let mut acc = 0.0;
                                for l in 0..3 {
                                    for n in 0..3 {
                                        acc += a_vol[k][l][m][n] * grad[i][l] * grad[j][n];
                                    }
                                }
                                stiffness.push(3 * ni + k, 3 * nj + m, dv * acc);
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(stiffness)
}

/// A [`Constitutive`] view that reports *only* the volumetric first Piola, so
/// [`material_tangent`] can difference the volumetric law in isolation.
struct VolumetricOnly<'a>(&'a dyn Constitutive);

impl Constitutive for VolumetricOnly<'_> {
    fn first_piola(&self, f: &Mat3) -> Mat3 {
        self.0.volumetric_piola(f)
    }
}

/// Componentwise difference of two fourth-order tensors.
fn sub4(a: &Tensor4, b: &Tensor4) -> Tensor4 {
    let mut out = [[[[0.0f64; 3]; 3]; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    out[i][j][k][l] = a[i][j][k][l] - b[i][j][k][l];
                }
            }
        }
    }
    out
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
pub fn tangent_stiffness_numerical<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
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
