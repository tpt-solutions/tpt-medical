//! Mixed `u`-`p` formulation for near/exact incompressibility — Q1/P0 on
//! the existing `Hex8`, per `rfcs/0012-mixed-up-formulation.md`.
//!
//! # The formulation
//!
//! Each element carries one element-constant pressure unknown `p_e` (the
//! Q1/P0 pairing RFC 0009 deferred and this RFC recommends for the
//! structured voxel meshes this workspace generates). The element's
//! deviatoric response is the tissue model's
//! [`TissueModel::mean_dilatation_first_piola`] evaluated with the
//! element-mean dilatation `J̄_e` substituted for the pointwise `J` (the
//! mean-dilatation treatment proper — see that method's doc comment for
//! why the substituted stress carries no `F^{-T}` term), and the
//! incompressibility constraint is enforced through the pressure:
//!
//! ```text
//! R_u = ∫ B^T P̄(F, J̄) dV  +  p_e · ∫ B^T cof(F) dV  −  f_ext
//! R_p = (J̄_e − 1)  −  ε̃ · p_e
//! ```
//!
//! with `cof(F) = J F^{-T}` evaluated **pointwise** (it differentiates the
//! true-`J` term `p(J̄−1)` — evaluating it at a mean-deformed `F` was one
//! of the assembly bugs the prototype hit), and `ε̃` the
//! **perturbed-Lagrangian compliance** — RFC 0012's candidate (a), the
//! first of the three ranked resolutions for the saddle solver. A pure
//! Lagrange multiplier has a zero pressure diagonal, and the first Newton
//! step from a zero pressure guess then produces an enormous `δp` that
//! line searches reject; the compliance `K_pp = −ε̃` regularizes exactly
//! that. It is a *perturbation*, not a penalty: `ε̃` enters nowhere in the
//! deviatoric path, and the constraint residual is satisfied to
//! `ε̃ · |p|`, so verification at exact-incompressibility tolerances runs
//! with `ε̃ = 1e-8` and asserts the closed forms directly.
//!
//! # The tangent is the whole residual, differenced
//!
//! `dR_p/du` involves the derivative of the element-mean `J̄` (a global
//! per-element quantity), and RFC 0012's prototype found two distinct
//! assembly bugs — a missing `1/V` in the `J̄` gradient and a
//! mean-vs-pointwise mixup in the constraint stress — through
//! finite-difference-vs-assembled checks. This implementation therefore
//! differenciates the **whole mixed residual** with respect to every DOF
//! for the Newton matrix: consistency with the residual is a structural
//! property, not something to re-derive per law. The cost (two residual
//! assemblies per DOF per iteration) is the same trade
//! [`crate::assembly::tangent_stiffness_numerical`] already makes for the
//! small verification problems this crate ships.
//!
//! # Scope
//!
//! Contact does not couple to the mixed path yet (RFC 0012's grand
//! cross-validation test against SRI-with-contact is the remaining piece);
//! SRI and mixed are mutually exclusive by construction — the mixed path
//! never reads the law's volumetric penalty at all.

use crate::mesh::{Mesh, MeshError};
use crate::solver::{Convergence, SolveError};
use tpt_fem_element::ReferenceElement;
use tpt_med_geometry::Mat3;
use tpt_med_tissue::TissueModel;

/// Options for [`solve_mixed_static`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MixedOptions {
    /// Perturbed-Lagrangian compliance `ε̃` (dimensionless: the constraint
    /// reads `J̄ − 1 = ε̃·p`). 1e-8 keeps the pressure error orders of
    /// magnitude below any engineering tolerance while keeping `K_pp`
    /// regular; RFC 0012 ranked this regularization first of the three
    /// saddle-solver candidates.
    pub compliance: f64,
    /// Tensor-product Gauss order per axis for the volume integrals.
    pub quadrature_order: usize,
    /// Central-difference step for the numerical tangent.
    pub fd_step: f64,
    /// Newton convergence settings. The convergence norm spans the free
    /// displacement rows and the (dimensionless) pressure rows.
    pub convergence: Convergence,
    /// Equal load increments to walk from zero to the full load, each
    /// warm-started from the previous converged state (the pressure field
    /// carried along). `1` is a single full-load solve — enough for small
    /// stretches; RFC 0012 prescribed fine increments alongside the
    /// compliance for larger ones.
    pub increments: usize,
}

impl Default for MixedOptions {
    fn default() -> Self {
        Self {
            compliance: 1.0e-8,
            quadrature_order: 3,
            fd_step: 1.0e-7,
            convergence: Convergence::default(),
            increments: 1,
        }
    }
}

impl MixedOptions {
    /// Options with the given quadrature order and the other defaults.
    pub fn with_quadrature_order(order: usize) -> Self {
        Self {
            quadrature_order: order,
            ..Self::default()
        }
    }
}

/// The converged result of a mixed `u`-`p` solve.
#[derive(Debug, Clone)]
pub struct MixedSolveResult {
    /// Converged nodal displacement vector (`3 * node_count`).
    pub displacement: Vec<f64>,
    /// Converged element pressure, in element order.
    pub pressure: Vec<f64>,
    /// Element-mean dilatation at the converged configuration (≈ 1 to
    /// `ε̃ · |p|` — this is the incompressibility error, directly).
    pub mean_dilatation: Vec<f64>,
    /// Free-DOF residual norm at the returned state.
    pub residual_norm: f64,
    /// Newton iterations taken.
    pub newton_iterations: usize,
}

/// Total DOF count of the mixed system: displacement DOFs plus one
/// element-constant pressure per element.
pub fn mixed_dof_count<E: ReferenceElement + crate::mesh::ElementFamily>(mesh: &Mesh<E>) -> usize {
    mesh.dof_count() + mesh.element_count()
}

/// Element-mean dilatation `J̄ = (1/V)∫J dV` on the options' quadrature rule.
fn mean_dilatation<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    element: usize,
    u: &[f64],
    opts: &MixedOptions,
) -> Result<f64, MeshError> {
    let rule = E::quadrature_rule(opts.quadrature_order);
    let mut volume = 0.0f64;
    let mut j_volume = 0.0f64;
    for (xi, w) in rule.points.iter().zip(&rule.weights) {
        let det = mesh.jacobian(element, xi).det();
        volume += w * det;
        let f = crate::assembly::element_deformation_gradient(mesh, element, u, xi).ok_or(
            MeshError::DegenerateElement {
                element,
                jacobian_determinant: det,
            },
        )?;
        // The numerator carries the same reference determinant as the
        // denominator: dV = det(J_ref) dxi on the [-1,1]^3 reference hex.
        j_volume += w * f.det() * det;
    }
    Ok(j_volume / volume)
}

/// The mixed residual: displacement rows first, then one dimensionless
/// pressure row per element.
#[allow(clippy::too_many_arguments)]
fn mixed_residual<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &TissueModel,
    load: &[f64],
    u: &[f64],
    pressure: &[f64],
    opts: &MixedOptions,
) -> Result<Vec<f64>, MeshError> {
    let n_u = mesh.dof_count();
    if u.len() != n_u || load.len() != n_u || pressure.len() != mesh.element_count() {
        return Err(MeshError::DofCountMismatch {
            expected: n_u,
            found: u.len(),
        });
    }
    let mut r = vec![0.0f64; mixed_dof_count(mesh)];
    r[..n_u].copy_from_slice(load);
    let rule = E::quadrature_rule(opts.quadrature_order);
    for e in 0..mesh.element_count() {
        // Element-mean dilatation on the same rule the deviatoric part
        // integrates on.
        let j_bar = {
            let mut vol = 0.0;
            let mut jv = 0.0;
            for (xi, w) in rule.points.iter().zip(&rule.weights) {
                let det = mesh.jacobian(e, xi).det();
                vol += w * det;
                let f = crate::assembly::element_deformation_gradient(mesh, e, u, xi).ok_or(
                    MeshError::DegenerateElement {
                        element: e,
                        jacobian_determinant: det,
                    },
                )?;
                // Reference determinant on both integrals: dV = det dxi.
                jv += w * f.det() * det;
            }
            jv / vol
        };
        // Displacement rows: deviatoric internal force + constraint stress.
        for (xi, w) in rule.points.iter().zip(&rule.weights) {
            let det = mesh.jacobian(e, xi).det();
            let Some(grad) = mesh.physical_gradients(e, xi) else {
                return Err(MeshError::DegenerateElement {
                    element: e,
                    jacobian_determinant: det,
                });
            };
            let f = crate::assembly::element_deformation_gradient(mesh, e, u, xi).ok_or(
                MeshError::DegenerateElement {
                    element: e,
                    jacobian_determinant: det,
                },
            )?;
            let j = f.det();
            if !j.is_finite() || j <= 0.0 {
                return Err(MeshError::InvertedDeformation {
                    element: e,
                    deformation_determinant: j,
                });
            }
            let p_bar = model.mean_dilatation_first_piola(&f, j_bar);
            let f_inv_t = f.inverse().map(|i| i.transpose()).unwrap_or(Mat3::ZERO);
            let dv = w * det;
            for (local, &node) in mesh.elements()[e].iter().enumerate() {
                for k in 0..3 {
                    let mut acc = 0.0;
                    for l in 0..3 {
                        acc +=
                            (p_bar.at(k, l) + pressure[e] * j * f_inv_t.at(k, l)) * grad[local][l];
                    }
                    r[3 * node + k] -= dv * acc;
                }
            }
        }
        // Pressure row, dimensionless: J̄ − 1 − ε̃·p.
        r[n_u + e] = j_bar - 1.0 - opts.compliance * pressure[e];
    }
    Ok(r)
}

/// The mixed Newton tangent: central differences of the whole residual
/// with respect to every DOF (displacement and pressure), so consistency
/// with [`mixed_residual`] is structural. See the module docs for why
/// RFC 0012's findings mandate this over a hand-contracted Jacobian.
fn mixed_tangent<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &TissueModel,
    load: &[f64],
    u: &[f64],
    pressure: &[f64],
    opts: &MixedOptions,
) -> Result<Vec<Vec<f64>>, MeshError> {
    let n = mixed_dof_count(mesh);
    let mut k = vec![vec![0.0f64; n]; n];
    for b in 0..n {
        let (up, um);
        if b < mesh.dof_count() {
            let mut probe = u.to_vec();
            probe[b] = u[b] + opts.fd_step;
            up = mixed_residual(mesh, model, load, &probe, pressure, opts)?;
            probe[b] = u[b] - opts.fd_step;
            um = mixed_residual(mesh, model, load, &probe, pressure, opts)?;
        } else {
            let e = b - mesh.dof_count();
            let mut probe = pressure.to_vec();
            probe[e] = pressure[e] + opts.fd_step;
            up = mixed_residual(mesh, model, load, u, &probe, opts)?;
            probe[e] = pressure[e] - opts.fd_step;
            um = mixed_residual(mesh, model, load, u, &probe, opts)?;
        }
        for a in 0..n {
            k[a][b] = (up[a] - um[a]) / (2.0 * opts.fd_step);
        }
    }
    Ok(k)
}

/// Solves the mixed `u`-`p` equilibrium `R(u, p) = 0` by Newton, with the
/// pressure regularized per the perturbed-Lagrangian compliance.
///
/// `dirichlet` condenses displacement DOFs only; pressure DOFs have no
/// essential conditions (RFC 0012's architecture constraint 1). The
/// external load applies to displacement DOFs; `load` is
/// `3 * node_count` long.
///
/// # Errors
///
/// As [`crate::solve_static`], plus [`SolveError::LoadSizeMismatch`] when
/// `load` is not `3 * node_count` long.
pub fn solve_mixed_static<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &TissueModel,
    load: &[f64],
    dirichlet: &[(usize, f64)],
    opts: &MixedOptions,
) -> Result<MixedSolveResult, SolveError> {
    let n_u = mesh.dof_count();
    if load.len() != n_u {
        return Err(SolveError::LoadSizeMismatch {
            expected: n_u,
            found: load.len(),
        });
    }
    let increments = opts.increments.max(1);
    let mut u = vec![0.0f64; n_u];
    for (dof, value) in dirichlet {
        u[*dof] = *value;
    }
    let mut pressure = vec![0.0f64; mesh.element_count()];
    let mut residual_norm = f64::INFINITY;
    let mut iterations = 0usize;
    for step in 1..=increments {
        let factor = step as f64 / increments as f64;
        let scaled: Vec<f64> = load.iter().map(|v| v * factor).collect();
        let r = newton_mixed_from(mesh, model, &scaled, dirichlet, opts, &u, &pressure)?;
        u = r.displacement;
        pressure = r.pressure;
        residual_norm = r.residual_norm;
        iterations += r.newton_iterations;
    }
    let mean_dilatation = (0..mesh.element_count())
        .map(|e| mean_dilatation(mesh, e, &u, opts))
        .collect::<Result<Vec<_>, _>>()
        .map_err(SolveError::Mesh)?;
    Ok(MixedSolveResult {
        displacement: u,
        pressure,
        mean_dilatation,
        residual_norm,
        newton_iterations: iterations,
    })
}

/// One full-load Newton solve from the supplied state (the continuation
/// increment of [`solve_mixed_static`]).
#[allow(clippy::too_many_arguments)]
fn newton_mixed_from<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &TissueModel,
    load: &[f64],
    dirichlet: &[(usize, f64)],
    opts: &MixedOptions,
    u_seed: &[f64],
    p_seed: &[f64],
) -> Result<MixedSolveResult, SolveError> {
    let n_u = mesh.dof_count();
    let n = mixed_dof_count(mesh);
    let fixed: std::collections::HashSet<usize> = dirichlet.iter().map(|(i, _)| *i).collect();
    let free: Vec<usize> = (0..n).filter(|i| !fixed.contains(i)).collect();
    let load_norm = load.iter().map(|x| x * x).sum::<f64>().sqrt();
    let tolerance = opts.convergence.abs_tol + opts.convergence.rel_tol * load_norm;

    let mut u = u_seed.to_vec();
    for (dof, value) in dirichlet {
        u[*dof] = *value;
    }
    let mut pressure = p_seed.to_vec();

    let mut iterations = 0usize;
    let mut residual_norm = f64::INFINITY;
    for _ in 0..opts.convergence.max_iter {
        iterations += 1;
        let r = mixed_residual(mesh, model, load, &u, &pressure, opts).map_err(SolveError::Mesh)?;
        residual_norm = free.iter().map(|&i| r[i] * r[i]).sum::<f64>().sqrt();
        if residual_norm <= tolerance {
            break;
        }
        let k = mixed_tangent(mesh, model, load, &u, &pressure, opts).map_err(SolveError::Mesh)?;
        // Condense fixed displacement DOFs; pressures are always free.
        let mut condensed = vec![vec![0.0f64; free.len()]; free.len()];
        for (a, &ra) in free.iter().enumerate() {
            for (b, &cb) in free.iter().enumerate() {
                condensed[a][b] = k[ra][cb];
            }
        }
        // Diagonal equilibration, as in the penalty solver: the mixed
        // diagonal spans displacement-stiffness and compliance entries.
        let mut scale = vec![1.0f64; free.len()];
        for (i, row) in condensed.iter().enumerate() {
            let diag = row[i];
            if diag.abs() > 0.0 && diag.is_finite() {
                scale[i] = 1.0 / diag.abs().sqrt();
            }
        }
        for (i, row) in condensed.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v *= scale[i] * scale[j];
            }
        }
        // Row-scaled right-hand side (once — the matrix carries s_i·s_j,
        // the vector only s_i).
        let rhs: Vec<f64> = free.iter().zip(&scale).map(|(&i, s)| r[i] * s).collect();
        // Dense solve with partial pivoting — the verification problems
        // this module ships are small, and a dense LU keeps the module
        // self-contained (the penalty solver's sparse path is unchanged).
        let delta = solve_dense(&condensed, &rhs).ok_or(SolveError::Singular)?;
        let step: Vec<f64> = delta.iter().zip(&scale).map(|(v, s)| v * s).collect();

        // Damped update on the full free-DOF norm.
        let mut alpha = 1.0f64;
        let mut accepted = false;
        for _ in 0..crate::solver::MAX_HALVINGS_MIXED {
            let mut trial_u = u.clone();
            let mut trial_p = pressure.clone();
            for (i, &dof) in free.iter().enumerate() {
                if dof < n_u {
                    trial_u[dof] -= alpha * step[i];
                } else {
                    trial_p[dof - n_u] -= alpha * step[i];
                }
            }
            if let Ok(r_trial) = mixed_residual(mesh, model, load, &trial_u, &trial_p, opts) {
                if r_trial.iter().all(|v| v.is_finite()) {
                    let norm = free
                        .iter()
                        .map(|&i| r_trial[i] * r_trial[i])
                        .sum::<f64>()
                        .sqrt();
                    if norm < residual_norm {
                        u = trial_u;
                        pressure = trial_p;
                        accepted = true;
                        break;
                    }
                }
            }
            alpha *= 0.5;
        }
        if !accepted {
            for (i, &dof) in free.iter().enumerate() {
                if dof < n_u {
                    u[dof] -= step[i];
                } else {
                    pressure[dof - n_u] -= step[i];
                }
            }
        }
    }
    if residual_norm > tolerance {
        return Err(SolveError::NotConverged {
            iterations,
            residual_norm,
            displacement: u,
        });
    }
    let mean_dilatation = (0..mesh.element_count())
        .map(|e| mean_dilatation(mesh, e, &u, opts))
        .collect::<Result<Vec<_>, _>>()
        .map_err(SolveError::Mesh)?;
    Ok(MixedSolveResult {
        displacement: u,
        pressure,
        mean_dilatation,
        residual_norm,
        newton_iterations: iterations,
    })
}

/// Dense LU with partial pivoting; `Ok(x)` for `A x = b` when `A` is
/// non-singular at the pivoting threshold.
fn solve_dense(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1.0e-13 {
            return None;
        }
        a.swap(pivot, col);
        b.swap(pivot, col);
        for row in (col + 1)..n {
            let factor = a[row][col] / a[col][col];
            if factor != 0.0 {
                for k in col..n {
                    a[row][k] -= factor * a[col][k];
                }
                b[row] -= factor * b[col];
            }
        }
    }
    let mut x = vec![0.0f64; n];
    for row in (0..n).rev() {
        let mut acc = b[row];
        for k in (row + 1)..n {
            acc -= a[row][k] * x[k];
        }
        x[row] = acc / a[row][row];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::{hex_box, Hex8Mesh};
    use tpt_med_tissue::NeoHookeanParams;

    /// Displacement-controlled uniaxial stretch of a single element with
    /// traction-free lateral faces - RFC 0012's verification item 1. The
    /// closed forms for incompressible Neo-Hookean (mu = 2 c10):
    ///
    /// - lateral stretch `lambda_t = lambda^(-1/2)` (the constraint
    ///   determines it - the point of the exercise),
    /// - Cauchy `sigma11 = mu (lambda^2 - 1/lambda)`, `sigma22 = 0`,
    /// - pressure `p = -mu/lambda` (the RFC's hand-derived value, which
    ///   the substituted deviatoric stress `2 c10 F` reproduces exactly:
    ///   `sigma22 = 2 c10 lambda_t^2 + p = 0` gives it directly).

    #[test]
    fn mixed_uniaxial_matches_the_incompressible_closed_form() {
        let mesh: Hex8Mesh = hex_box(1, 1, 1, 1.0, 1.0, 1.0).expect("mesh");
        let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 1.0 });
        let lambda = 1.2;
        let mut dirichlet = Vec::new();
        for n in mesh.face_nodes(0, false) {
            dirichlet.push((3 * n, 0.0)); // x = 0 face: no axial motion
        }
        for n in mesh.face_nodes(1, false) {
            dirichlet.push((3 * n + 1, 0.0)); // y = 0 plane pinned
        }
        for n in mesh.face_nodes(2, false) {
            dirichlet.push((3 * n + 2, 0.0)); // z = 0 plane pinned
        }
        for n in mesh.face_nodes(0, true) {
            dirichlet.push((3 * n, lambda - 1.0)); // stretched face
        }
        let opts = MixedOptions::default();
        let result = solve_mixed_static(
            &mesh,
            &model,
            &vec![0.0; mesh.dof_count()],
            &dirichlet,
            &opts,
        )
        .expect("converges");

        // Incompressibility, to the compliance's own scale.
        assert!(
            (result.mean_dilatation[0] - 1.0).abs() < 1.0e-6,
            "J_bar = {}",
            result.mean_dilatation[0]
        );
        // Lateral stretch free and determined by the constraint.
        let mu = 1.0;
        let lambda_t_expected = lambda.powf(-0.5);
        for n in mesh.face_nodes(1, true) {
            let uy = result.displacement[3 * n + 1];
            assert!(
                (uy - (lambda_t_expected - 1.0)).abs() < 1.0e-4,
                "node {n} u_y = {uy}, expected {}",
                lambda_t_expected - 1.0
            );
        }
        // Pressure: p = -mu/lambda, the RFC's hand-derived value.
        let p_expected = -mu / lambda;
        assert!(
            (result.pressure[0] - p_expected).abs() < 1.0e-4,
            "p = {} vs {p_expected}",
            result.pressure[0]
        );
        // Cauchy stress at the (uniform) element: sigma = (P_bar + p cof) F^T / J.
        let f = crate::assembly::element_deformation_gradient(
            &mesh,
            0,
            &result.displacement,
            &[0.5, 0.5, 0.5],
        )
        .expect("valid element");
        let f_inv_t = f.inverse().expect("non-singular");
        let j = f.det();
        let p_bar = model.mean_dilatation_first_piola(&f, result.mean_dilatation[0]);
        let mut p_total = Mat3::ZERO;
        for i in 0..3 {
            for jj in 0..3 {
                let v = p_bar.at(i, jj) + result.pressure[0] * j * f_inv_t.at(i, jj);
                p_total.set(i, jj, v);
            }
        }
        let cauchy = |k: usize| (0..3).map(|m| p_total.at(k, m) * f.at(k, m)).sum::<f64>() / j;
        let sigma11 = cauchy(0);
        let sigma22 = cauchy(1);
        assert!(
            (sigma11 - mu * (lambda * lambda - 1.0 / lambda)).abs() < 1.0e-4,
            "sigma11 = {sigma11}"
        );
        assert!(sigma22.abs() < 1.0e-4, "traction-free lateral: {sigma22}");
    }

    /// RFC 0012's verification item 2: a uniform-strain state reproduces
    /// constant stress and constant pressure - Q1/P0 passes the patch on
    /// uniform hexes.
    #[test]
    fn mixed_patch_state_gives_constant_pressure() {
        let mesh: Hex8Mesh = hex_box(2, 1, 1, 2.0, 1.0, 1.0).expect("mesh");
        let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 1.0 });
        let lambda = 1.1;
        let mut dirichlet = Vec::new();
        for n in mesh.face_nodes(0, false) {
            dirichlet.push((3 * n, 0.0));
        }
        for n in mesh.face_nodes(1, false) {
            dirichlet.push((3 * n + 1, 0.0));
        }
        for n in mesh.face_nodes(2, false) {
            dirichlet.push((3 * n + 2, 0.0));
        }
        for n in mesh.face_nodes(0, true) {
            dirichlet.push((3 * n, lambda - 1.0));
        }
        let opts = MixedOptions::default();
        let result = solve_mixed_static(
            &mesh,
            &model,
            &vec![0.0; mesh.dof_count()],
            &dirichlet,
            &opts,
        )
        .expect("converges");
        assert_eq!(result.pressure.len(), 2);
        assert!(
            (result.pressure[0] - result.pressure[1]).abs() < 1.0e-9,
            "uniform state must give constant pressure: {:?}",
            result.pressure
        );
        for jd in &result.mean_dilatation {
            assert!((jd - 1.0).abs() < 1.0e-6, "J_bar {jd}");
        }
    }

    /// RFC 0012's verification item 3: a near-incompressible slender column
    /// under transverse load. Mixed displacement matches SRI at the same
    /// near-incompressible penalty (both approach the exact constraint),
    /// while full integration visibly locks.
    #[test]
    fn mixed_column_matches_sri_where_full_integration_locks() {
        use crate::assembly::AssemblyOptions;
        use crate::solver::SolveOptions;
        // Two elements through the thickness: Q1/P0's element-constant
        // pressure needs resolution across the bending gradient (a single
        // element cannot represent the linear bending pressure, which is
        // the pairing's known coarse-mesh limitation, not a solver bug).
        let mesh: Hex8Mesh = hex_box(1, 2, 4, 1.0, 1.0, 4.0).expect("mesh");
        // Near-incompressible penalty (volumetric/shear ratio 1000, i.e.
        // nu ~ 0.4995). Stiffer penalties were tried and are the reason
        // the load-path driver exists — but at 1e-6 the penalty solve
        // itself inverts elements during cutback; 1e-3 carries a ~0.1%
        // compliance error, far below the comparison tolerance below.
        let d1 = 1.0e-3;
        let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1 });

        // Transverse tip load on the z-max face, x direction. Small enough
        // for a single-increment solve (the mixed driver has no
        // continuation yet; RFC 0012's load-path integration is future
        // work).
        let mut load = vec![0.0; mesh.dof_count()];
        for n in mesh.face_nodes(2, true) {
            load[3 * n] = 0.0025; // 4 nodes, 0.01 N total
        }
        let mut dirichlet = Vec::new();
        for n in mesh.face_nodes(2, false) {
            dirichlet.push((3 * n, 0.0));
            dirichlet.push((3 * n + 1, 0.0));
            dirichlet.push((3 * n + 2, 0.0));
        }

        let mixed_opts = MixedOptions {
            quadrature_order: 2,
            increments: 8,
            ..MixedOptions::default()
        };
        let mixed = solve_mixed_static(&mesh, &model, &load, &dirichlet, &mixed_opts)
            .expect("mixed converges");
        // The penalty runs go through the load-path driver: a stiff
        // penalty in one increment is exactly what the continuation
        // machinery exists for. Their convergence carries a looser
        // relative tolerance deliberately: the penalty formulation's own
        // compliance error (the d1 = 1e-6 constraint violation it exists
        // to approximate) is orders of magnitude above 1e-6 relative, so
        // tightening the Newton tolerance further refines a number whose
        // dominant error is the penalty's, not the iteration's.
        let penalty_convergence = crate::solver::Convergence {
            abs_tol: 1.0e-8,
            rel_tol: 1.0e-6,
            max_iter: 50,
        };
        let sri_opts = SolveOptions {
            assembly: AssemblyOptions::with_volumetric_quadrature_order(1),
            convergence: penalty_convergence,
        };
        let path = crate::LoadPathOptions {
            steps: 8,
            ..crate::LoadPathOptions::default()
        };
        let sri = crate::solve_load_path(&mesh, &model, &load, &dirichlet, &sri_opts, None, path)
            .expect("sri path converges");
        let full = crate::solve_load_path(
            &mesh,
            &model,
            &load,
            &dirichlet,
            &SolveOptions {
                convergence: penalty_convergence,
                ..SolveOptions::default()
            },
            None,
            crate::LoadPathOptions {
                steps: 8,
                ..crate::LoadPathOptions::default()
            },
        )
        .expect("full path converges");

        let tip = |u: &[f64]| {
            mesh.face_nodes(2, true)
                .iter()
                .map(|&n| u[3 * n])
                .sum::<f64>()
                / 4.0
        };
        let (u_mixed, u_sri, u_full) = (
            tip(&mixed.displacement),
            tip(&sri.last().expect("path").displacement),
            tip(&full.last().expect("path").displacement),
        );
        // The locking ladder RFC 0012 item 3 asked for: full integration
        // locks catastrophically (18x stiffer than the constraint here),
        // SRI recovers most but visibly not all of it (its documented
        // residual over-stiffening at nu -> 1/2 — the reason the mixed
        // path exists), and the mixed solution is the softest, with the
        // element-mean dilatation at 1 to the compliance scale.
        assert!(
            u_full < 0.2 * u_sri,
            "full integration must lock hard against SRI: {u_full} vs {u_sri}"
        );
        assert!(
            u_sri < u_mixed && u_mixed < 1.4 * u_sri,
            "SRI must over-stiffen modestly, not wildly: sri {u_sri} vs mixed {u_mixed}"
        );
        for jd in &mixed.mean_dilatation {
            assert!((jd - 1.0).abs() < 1.0e-6, "mixed J_bar {jd}");
        }
    }

    /// RFC 0012's verification item 4: no checkerboard pressure modes on
    /// the structured mesh - a uniformly loaded block carries a uniform
    /// pressure within the state's own variation.
    #[test]
    fn mixed_pressure_field_has_no_checkerboard() {
        let mesh: Hex8Mesh = hex_box(2, 2, 2, 1.0, 1.0, 1.0).expect("mesh");
        let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 1.0 });
        let lambda = 1.15;
        let mut dirichlet = Vec::new();
        for n in mesh.face_nodes(0, false) {
            dirichlet.push((3 * n, 0.0));
        }
        for n in mesh.face_nodes(1, false) {
            dirichlet.push((3 * n + 1, 0.0));
        }
        for n in mesh.face_nodes(2, false) {
            dirichlet.push((3 * n + 2, 0.0));
        }
        for n in mesh.face_nodes(0, true) {
            dirichlet.push((3 * n, lambda - 1.0));
        }
        let result = solve_mixed_static(
            &mesh,
            &model,
            &vec![0.0; mesh.dof_count()],
            &dirichlet,
            &MixedOptions::default(),
        )
        .expect("converges");
        let p_max = result
            .pressure
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let p_min = result
            .pressure
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        assert!(
            p_max - p_min < 1.0e-6 * p_max.abs().max(1.0),
            "checkerboard would separate the pressures: {p_min}..{p_max}"
        );
        for jd in &result.mean_dilatation {
            assert!((jd - 1.0).abs() < 1.0e-6);
        }
    }
}
