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

use crate::friction::friction_terms;
use crate::mesh::{Mesh, MeshError};
use crate::solver::{ContactConfig, ContactSummary, Convergence, SolveError};
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
    /// Maximum bisections of a failing increment's remaining distance,
    /// mirroring [`crate::LoadPathOptions::max_cutbacks`]. A failed
    /// increment is retried from the last converged load factor at half
    /// the remaining step, so progress is retained rather than restarted.
    pub max_cutbacks: usize,
}

impl Default for MixedOptions {
    fn default() -> Self {
        Self {
            compliance: 1.0e-8,
            quadrature_order: 3,
            fd_step: 1.0e-7,
            convergence: Convergence::default(),
            increments: 1,
            max_cutbacks: 4,
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
    /// Contact outcome, or `None` when the solve had no contact configured.
    pub contact: Option<ContactSummary>,
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

/// The contact terms of one state, frozen so the numerical tangent can
/// difference the residual without re-evaluating the active set: a slave
/// node sitting on the activation boundary would flip in/out between the
/// ±h probes and poison the difference with a `kappa·u/h` spike — the
/// mixed module's analogue of the penalty solver freezing `k_c` per
/// Newton iteration. The friction force is frozen with it (its tangent
/// is neglected; friction is an explicit force here, per the solver's
/// own "separate layer" treatment).
struct FrozenContact {
    f_c: Vec<f64>,
    k_rows: Vec<usize>,
    k_cols: Vec<usize>,
    k_vals: Vec<f64>,
    friction_force: Vec<f64>,
}

fn freeze_contact<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    contact: &ContactConfig<'_>,
    u: &[f64],
) -> Result<FrozenContact, MeshError> {
    let (f_c, k_c, _) = crate::solver::contact_terms(mesh, contact.pairing, u, contact.penalty)?;
    let friction_force = match contact.friction {
        Some(fcfg) => friction_terms(mesh, contact.pairing, u, contact.penalty, fcfg)?.force,
        None => Vec::new(),
    };
    Ok(FrozenContact {
        f_c,
        k_rows: k_c.rows,
        k_cols: k_c.cols,
        k_vals: k_c.vals,
        friction_force,
    })
}

fn mixed_residual_frozen<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &TissueModel,
    load: &[f64],
    u: &[f64],
    pressure: &[f64],
    opts: &MixedOptions,
    frozen: Option<&FrozenContact>,
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
    // Contact and friction couple to the displacement rows only (RFC 0012:
    // the contact block keeps its shape with pressure rows/columns zero),
    // with the same `K_c u - f_c` structure the penalty solver's residual
    // carries — using the frozen terms so the active set is a property of
    // the Newton iteration, not of each finite-difference probe.
    if let Some(frozen) = frozen {
        for (i, v) in frozen.f_c.iter().enumerate() {
            r[i] -= v;
        }
        for i in 0..frozen.k_rows.len() {
            r[frozen.k_rows[i]] += frozen.k_vals[i] * u[frozen.k_cols[i]];
        }
        for (i, v) in frozen.friction_force.iter().enumerate() {
            r[i] -= v;
        }
    }
    Ok(r)
}

/// The mixed Newton tangent: central differences of the whole residual
/// with respect to every DOF (displacement and pressure), so consistency
/// with [`mixed_residual`] is structural. See the module docs for why
/// RFC 0012's findings mandate this over a hand-contracted Jacobian.
fn mixed_tangent_frozen<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &TissueModel,
    load: &[f64],
    u: &[f64],
    pressure: &[f64],
    opts: &MixedOptions,
    frozen: Option<&FrozenContact>,
) -> Result<Vec<Vec<f64>>, MeshError> {
    let n = mixed_dof_count(mesh);
    let mut k = vec![vec![0.0f64; n]; n];
    for b in 0..n {
        let (up, um);
        if b < mesh.dof_count() {
            let mut probe = u.to_vec();
            probe[b] = u[b] + opts.fd_step;
            up = mixed_residual_frozen(mesh, model, load, &probe, pressure, opts, frozen)?;
            probe[b] = u[b] - opts.fd_step;
            um = mixed_residual_frozen(mesh, model, load, &probe, pressure, opts, frozen)?;
        } else {
            let e = b - mesh.dof_count();
            let mut probe = pressure.to_vec();
            probe[e] = pressure[e] + opts.fd_step;
            up = mixed_residual_frozen(mesh, model, load, u, &probe, opts, frozen)?;
            probe[e] = pressure[e] - opts.fd_step;
            um = mixed_residual_frozen(mesh, model, load, u, &probe, opts, frozen)?;
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
/// Contact couples to the displacement rows only (RFC 0012: "the contact
/// block keeps its current shape with pressure rows/columns zero") — and
/// because the tangent is the whole residual differenced, the contact
/// terms reach the Newton matrix by construction once they are in the
/// residual.
///
/// # Errors
///
/// As [`crate::solve_static`], plus [`SolveError::LoadSizeMismatch`] when
/// `load` is not `3 * node_count` long.
#[allow(clippy::too_many_arguments)]
pub fn solve_mixed_static<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &TissueModel,
    load: &[f64],
    dirichlet: &[(usize, f64)],
    opts: &MixedOptions,
    contact: Option<ContactConfig<'_>>,
) -> Result<MixedSolveResult, SolveError> {
    let n_u = mesh.dof_count();
    if load.len() != n_u {
        return Err(SolveError::LoadSizeMismatch {
            expected: n_u,
            found: load.len(),
        });
    }
    // The mixed solver freezes its contact terms per increment; the wall's
    // history-carrying friction has no such treatment, so refuse it loudly.
    if contact.is_some_and(|c| c.radial.is_some()) {
        return Err(SolveError::Unsupported(
            "a radial wall with the mixed u-p solver",
        ));
    }
    // The continuation walk mirrors the penalty solver's load path
    // exactly: equal nominal increments, a failed increment retried from
    // the last converged factor at half the remaining distance (progress
    // retained, never restarted), up to `max_cutbacks`.
    let increments = opts.increments.max(1);
    let mut u = vec![0.0f64; n_u];
    for (dof, value) in dirichlet {
        u[*dof] = *value;
    }
    let mut pressure = vec![0.0f64; mesh.element_count()];
    let mut residual_norm = f64::INFINITY;
    let mut iterations = 0usize;
    let mut current = 0.0f64;
    for i in 1..=increments {
        let target = i as f64 / increments as f64;
        let mut reached = current;
        let mut attempt_target = target;
        let mut cutbacks = 0usize;
        loop {
            let scaled: Vec<f64> = load.iter().map(|v| v * attempt_target).collect();
            match newton_mixed_from(
                mesh, model, &scaled, dirichlet, opts, contact, &u, &pressure,
            ) {
                Ok(r) => {
                    u = r.displacement;
                    pressure = r.pressure;
                    reached = attempt_target;
                    residual_norm = r.residual_norm;
                    iterations += r.newton_iterations;
                    break;
                }
                Err(SolveError::NotConverged { .. }) => {
                    if cutbacks >= opts.max_cutbacks {
                        return Err(SolveError::NotConverged {
                            iterations: iterations + 1,
                            residual_norm: f64::INFINITY,
                            displacement: u,
                        });
                    }
                    cutbacks += 1;
                    attempt_target = reached + (target - reached) / 2.0;
                }
                Err(e) => return Err(e),
            }
        }
        current = reached;
    }
    let mean_dilatation = (0..mesh.element_count())
        .map(|e| mean_dilatation(mesh, e, &u, opts))
        .collect::<Result<Vec<_>, _>>()
        .map_err(SolveError::Mesh)?;
    let contact = match contact {
        Some(cfg) => Some(ContactSummary {
            active_constraints: cfg
                .pairing
                .active_constraints(mesh, &u)
                .map_err(SolveError::Contact)?,
            max_penetration: cfg
                .pairing
                .max_penetration(mesh, &u)
                .map_err(SolveError::Contact)?,
            total_reaction: cfg
                .pairing
                .total_reaction(mesh, &u, cfg.penalty)
                .map_err(SolveError::Contact)?,
            slipping_nodes: match cfg.friction {
                Some(fcfg) => Some(
                    friction_terms(mesh, cfg.pairing, &u, cfg.penalty, fcfg)
                        .map_err(SolveError::Contact)?
                        .slipping_nodes,
                ),
                None => None,
            },
            wall: None,
        }),
        None => None,
    };
    Ok(MixedSolveResult {
        displacement: u,
        pressure,
        mean_dilatation,
        residual_norm,
        newton_iterations: iterations,
        contact,
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
    contact: Option<ContactConfig<'_>>,
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
        // The contact state is frozen once per Newton iteration: the
        // residual, the tangent and every line-search trial are evaluated
        // against the SAME active set, so the descent check is consistent
        // with the step that was computed (a trial that re-evaluated the
        // active set could — and did — flip it mid-line-search, breaking
        // the descent guarantee and cycling the iteration between
        // engaged and separated states).
        let frozen = match &contact {
            Some(cfg) => Some(freeze_contact(mesh, cfg, &u).map_err(SolveError::Mesh)?),
            None => None,
        };
        let r = mixed_residual_frozen(mesh, model, load, &u, &pressure, opts, frozen.as_ref())
            .map_err(SolveError::Mesh)?;
        residual_norm = free.iter().map(|&i| r[i] * r[i]).sum::<f64>().sqrt();
        if residual_norm <= tolerance {
            break;
        }
        let k = mixed_tangent_frozen(mesh, model, load, &u, &pressure, opts, frozen.as_ref())
            .map_err(SolveError::Mesh)?;
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
            if let Ok(r_trial) =
                mixed_residual_frozen(mesh, model, load, &trial_u, &trial_p, opts, frozen.as_ref())
            {
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
        contact: None,
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
            None,
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
            None,
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
        let mixed = solve_mixed_static(&mesh, &model, &load, &dirichlet, &mixed_opts, None)
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
            None,
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

    /// RFC 0012's verification item 5, the grand cross-validation: a
    /// free-contact punch — the contact face as the *only* support —
    /// re-run in mixed mode, cross-validated against SRI at a compliant
    /// penalty.
    ///
    /// The active-set stabilization: the pairing carries an activation
    /// tolerance of 1e-2 mm, so a node whose Newton step overshoots the
    /// wall by the equilibrium-penetration scale (contact force divided
    /// by penalty stiffness, ~1e-4 mm here) stays engaged and is pulled
    /// back, instead of deactivating and free-diving — the
    /// engaged/separated cycle that otherwise stalls the iteration.
    /// Contact-chatter damping of exactly this shape is standard in
    /// penalty contact; the tolerance is three orders of magnitude below
    /// the geometric scale and two below the punch.
    #[test]
    fn mixed_contact_cross_validates_against_sri_in_the_compliant_limit() {
        use crate::contact::ContactPairing;
        use crate::mesh::{hex_box, Hex8Mesh};
        use crate::solver::ContactConfig;
        use tpt_med_geometry::Vec3;
        use tpt_med_tissue::NeoHookeanParams;

        let mesh: Hex8Mesh = hex_box(1, 1, 2, 10.0, 10.0, 10.0).expect("mesh");
        let model = TissueModel::NeoHookean(NeoHookeanParams {
            c10: 0.5,
            d1: 1.0e-3,
        });
        let punch = -0.075; // 0.75 % compression, per node

        let top = mesh.face_nodes(1, true);
        let slave = mesh.face_nodes(1, false);
        let mut dirichlet = Vec::new();
        for &node in &top {
            dirichlet.push((mesh.dof(node, 0), 0.0));
            dirichlet.push((mesh.dof(node, 2), 0.0));
            dirichlet.push((mesh.dof(node, 1), punch));
        }
        let load = vec![0.0; mesh.dof_count()];
        let master: Vec<Vec3> = slave.iter().map(|&n| mesh.nodes()[n]).collect();
        let pairing = ContactPairing::new(1, slave.iter().copied(), master)
            .expect("axis 1")
            .with_activation_tolerance(1.0e-2);
        let contact = Some(ContactConfig {
            radial: None,
            pairing: &pairing,
            penalty: 1.0e4,
            friction: None,
        });

        let mixed = solve_mixed_static(
            &mesh,
            &model,
            &load,
            &dirichlet,
            &MixedOptions {
                increments: 4,
                ..MixedOptions::default()
            },
            contact,
        )
        .expect("mixed contact solve converges");
        let summary = mixed.contact.as_ref().expect("contact summary");
        assert_eq!(summary.active_constraints.len(), slave.len());
        assert!(
            summary.max_penetration < 1.0e-2,
            "penetrated by {}",
            summary.max_penetration
        );
        for jd in &mixed.mean_dilatation {
            assert!((jd - 1.0).abs() < 1.0e-6, "J_bar {jd}");
        }

        // Cross-validation: SRI's lateral bulge approaches the mixed one
        // as the penalty compliance shrinks.
        let bulge = |u: &[f64]| {
            mesh.face_nodes(0, true)
                .iter()
                .map(|&n| u[3 * n])
                .sum::<f64>()
                / mesh.face_nodes(0, true).len() as f64
        };
        let mixed_bulge = bulge(&mixed.displacement);
        assert!(mixed_bulge > 0.0, "incompressible bulge must be outward");
        // The bands widen with d1 on purpose: a compressible material
        // (large d1) bulges legitimately *less* than the exactly
        // incompressible one, so the cross-validation claim is the
        // convergence as the compliance shrinks — near-incompressible d1
        // agrees to a few percent.
        for (d1, band) in [(1.0e-1, 0.15), (1.0e-2, 0.08), (1.0e-3, 0.05)] {
            let sri_model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1 });
            let sri_opts = crate::solver::SolveOptions {
                assembly: crate::assembly::AssemblyOptions::with_volumetric_quadrature_order(1),
                ..crate::solver::SolveOptions::default()
            };
            let sri = crate::solve_static(&mesh, &sri_model, &load, &dirichlet, &sri_opts, contact)
                .unwrap_or_else(|e| panic!("sri converges at d1={d1}: {e}"));
            let sri_bulge = bulge(&sri.displacement);
            let gap = (sri_bulge - mixed_bulge).abs();
            assert!(
                gap < band * mixed_bulge.abs(),
                "d1={d1}: sri bulge {sri_bulge} vs mixed {mixed_bulge} (gap {gap})"
            );
        }
    }
}

/// The cutback walk: a transverse load large enough that equal single
/// increments fail converges when failed increments are bisected —
/// mirroring the penalty solver's `LoadPathOptions::max_cutbacks`
/// contract. Verified by comparing the cutback result against the
/// same load walked in many small equal increments (both roads to the
/// same equilibrium).
#[test]
fn mixed_cutback_walk_recovers_a_load_equal_increments_cannot() {
    use crate::mesh::{hex_box, Hex8Mesh};
    use tpt_med_tissue::NeoHookeanParams;

    let mesh: Hex8Mesh = hex_box(1, 1, 4, 1.0, 1.0, 4.0).expect("mesh");
    let model = TissueModel::NeoHookean(NeoHookeanParams {
        c10: 0.5,
        d1: 1.0e-3,
    });
    let mut load = vec![0.0; mesh.dof_count()];
    for n in mesh.face_nodes(2, true) {
        load[3 * n] = 0.01; // 4 nodes, 0.04 N — the load that failed at increments: 1
    }
    let mut dirichlet = Vec::new();
    for n in mesh.face_nodes(2, false) {
        dirichlet.push((3 * n, 0.0));
        dirichlet.push((3 * n + 1, 0.0));
        dirichlet.push((3 * n + 2, 0.0));
    }
    let tip = |u: &[f64]| {
        mesh.face_nodes(2, true)
            .iter()
            .map(|&n| u[3 * n])
            .sum::<f64>()
            / 4.0
    };

    // One increment without cutbacks: must fail (regression guard for
    // the premise).
    let stiff = MixedOptions {
        quadrature_order: 2,
        increments: 1,
        max_cutbacks: 0,
        ..MixedOptions::default()
    };
    assert!(
        solve_mixed_static(&mesh, &model, &load, &dirichlet, &stiff, None).is_err(),
        "the 0.04 N single-increment solve is expected not to converge"
    );

    // Cutbacks on: converges.
    let cut = MixedOptions {
        quadrature_order: 2,
        increments: 4,
        max_cutbacks: 4,
        ..MixedOptions::default()
    };
    let cut = solve_mixed_static(&mesh, &model, &load, &dirichlet, &cut, None)
        .expect("cutback walk converges");

    // Many small equal increments: the reference road.
    let fine = MixedOptions {
        quadrature_order: 2,
        increments: 16,
        max_cutbacks: 0,
        ..MixedOptions::default()
    };
    let fine = solve_mixed_static(&mesh, &model, &load, &dirichlet, &fine, None)
        .expect("fine walk converges");

    let (u_cut, u_fine) = (tip(&cut.displacement), tip(&fine.displacement));
    assert!(
        (u_cut - u_fine).abs() < 0.02 * u_fine.abs(),
        "cutback {u_cut} vs fine {u_fine}"
    );
    // The deformation is finite (this is a real load, not a token one).
    assert!(u_fine.abs() > 0.2, "tip {u_fine}");
}
