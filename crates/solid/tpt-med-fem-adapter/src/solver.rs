//! Nonlinear static equilibrium: Newton-Raphson on the assembled residual,
//! with optional unilateral contact re-evaluated at every iteration.
//!
//! The linear algebra is the substrate's: `tpt_fem_sparse::Coo` for the global
//! assembly and `tpt_fem_sparse::solve` for each Newton step. The iteration
//! itself is this crate's, and the reason is worth stating plainly.
//!
//! # Why not `tpt_fem_solve::newton`
//!
//! The substrate's Newton driver tests `||R(u)||_2` over the **whole** vector
//! against an **absolute** tolerance before condensing the Dirichlet DOFs. For
//! a displacement-controlled problem — which is how a hyperelastic uniaxial
//! test, and most hyperelastic analyses, are posed — the residual at a
//! prescribed DOF is the *reaction*, which is large and non-zero by
//! definition, so `||R||` has a floor that no iteration can reduce and the
//! driver always reports `MaxIterations` even when every free DOF is in
//! equilibrium to machine precision. RFC 0009 named `tpt-fem-solve::newton` as
//! the intended driver; the loop below keeps the same structure (condense the
//! essential DOFs, solve, update) and fixes the convergence measure, which is
//! the free-DOF residual norm against an absolute *plus relative* tolerance
//! scaled by the applied load.
//!
//! # Sign convention
//!
//! The residual is `R(u) = f_int(u) - f_ext + r_contact(u)`, where `r_contact`
//! is what `tpt_fem_contact::penalty_contact` contributes. That function
//! returns `(K, f)` such that `solve(K + penalty_terms, f + penalty_terms)`
//! holds the constrained DOF at `lower`; equivalently its contribution to a
//! residual written as `K u - f` is `kappa * (u_dof - lower)`. So the residual
//! gains `+kappa * (u_dof - lower)` and the Jacobian diagonal gains `+kappa`
//! for each active constraint — assembled here by calling `penalty_contact`
//! with an empty base matrix and a zero load vector, so the active-set
//! bookkeeping is the substrate's rather than a re-derivation of it.
//!
//! Because the active set is recomputed from the current geometry inside both
//! closures, a node that separates from the obstacle simply stops being
//! constrained on the next iteration; nothing has to be told to "release" it.
//!
//! # Friction
//!
//! Friction is a separate layer ([`crate::friction`]) and a separate field on
//! [`ContactConfig`], not a change to the loop above. The normal problem stays
//! exactly as it was, so a caller who leaves `friction: None` gets the
//! frictionless behaviour this module was verified against, unchanged.

use crate::assembly::{internal_force, tangent_stiffness, AssemblyOptions, Constitutive};
use crate::contact::ContactPairing;
use crate::friction::{friction_terms, FrictionConfig};
use crate::mesh::{Mesh, MeshError};
use crate::wall::{RadialWall, WallSummary};
use std::cell::Cell;
use std::collections::HashSet;
use tpt_fem_element::ReferenceElement;
use tpt_fem_sparse::{solve, Coo};

/// Maximum number of step halvings the line search will try.
const MAX_HALVINGS: usize = 12;

/// The mixed solver's line-search cap — same policy, separate module.
pub(crate) const MAX_HALVINGS_MIXED: usize = 12;

/// Errors returned by [`solve_static`].
#[derive(Debug)]
pub enum SolveError {
    /// A mesh or assembly input was malformed.
    Mesh(MeshError),
    /// The external load vector had the wrong length.
    LoadSizeMismatch {
        /// The length the solve requires (`3 * node_count`).
        expected: usize,
        /// The length supplied.
        found: usize,
    },
    /// A contact pairing could not be evaluated (bad DOF count, or a slave node
    /// outside the mesh).
    Contact(MeshError),
    /// The Newton iteration did not reach the requested tolerance.
    NotConverged {
        /// Iterations performed.
        iterations: usize,
        /// Final free-DOF residual norm.
        residual_norm: f64,
        /// The best displacement reached, so a caller can inspect how far the
        /// solve got instead of having to start over.
        displacement: Vec<f64>,
    },
    /// The condensed linear system was singular (a floating structure, or a
    /// Dirichlet set that leaves a rigid-body mode unconstrained).
    Singular,
    /// The requested combination is not supported by this solver (named in
    /// the message) — rejected rather than silently ignored.
    Unsupported(&'static str),
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolveError::Mesh(e) => write!(f, "{e}"),
            SolveError::LoadSizeMismatch { expected, found } => {
                write!(f, "load vector has {found} entries, expected {expected}")
            }
            SolveError::Contact(e) => write!(f, "contact evaluation failed: {e}"),
            SolveError::NotConverged {
                iterations,
                residual_norm,
                ..
            } => write!(
                f,
                "newton did not converge in {iterations} iterations \
                 (free-DOF residual {residual_norm:e})"
            ),
            SolveError::Singular => write!(f, "the condensed linear system is singular"),
            SolveError::Unsupported(what) => write!(f, "unsupported: {what}"),
        }
    }
}

impl std::error::Error for SolveError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            SolveError::Mesh(e) | SolveError::Contact(e) => Some(e),
            SolveError::LoadSizeMismatch { .. }
            | SolveError::NotConverged { .. }
            | SolveError::Singular
            | SolveError::Unsupported(_) => None,
        }
    }
}

impl From<MeshError> for SolveError {
    fn from(e: MeshError) -> Self {
        SolveError::Mesh(e)
    }
}

/// Convergence settings for the Newton iteration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Convergence {
    /// Absolute tolerance on the free-DOF residual norm. Guards the
    /// load-free case (`f_ext = 0`, e.g. a pure eigen- or residual check).
    pub abs_tol: f64,
    /// Relative tolerance on the free-DOF residual norm, scaled by the applied
    /// load norm. This is the criterion that makes a displacement-controlled
    /// problem converge: the reactions at the prescribed DOFs never enter it.
    pub rel_tol: f64,
    /// Maximum number of Newton iterations.
    pub max_iter: usize,
}

impl Default for Convergence {
    fn default() -> Self {
        Self {
            abs_tol: 1.0e-10,
            rel_tol: 1.0e-10,
            max_iter: 50,
        }
    }
}

/// Solver controls: convergence settings plus the assembly options.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SolveOptions {
    /// Newton convergence settings.
    pub convergence: Convergence,
    /// Quadrature order and finite-difference step for the assembly.
    pub assembly: AssemblyOptions,
}

/// Contact configuration for a solve: which pairing, and the penalty stiffness.
#[derive(Debug, Clone, Copy)]
pub struct ContactConfig<'a> {
    /// The pairing to enforce.
    pub pairing: &'a ContactPairing,
    /// Penalty stiffness added to each active constrained DOF. Should be
    /// several orders of magnitude above the structural stiffness, and large
    /// enough that the resulting penetration is below the caller's tolerance.
    pub penalty: f64,
    /// Optional Coulomb friction on the active contacts. `None` is exactly the
    /// previous frictionless behaviour, and costs nothing at solve time.
    pub friction: Option<FrictionConfig>,
    /// Optional rigid analytic cylindrical wall (a vessel lumen), enforced by
    /// the same penalty. `None` costs nothing. It composes with `pairing`:
    /// for a wall-only problem pass [`ContactPairing::inactive`] as the
    /// pairing. Wall friction lives on the [`RadialWall`] itself, because it
    /// carries committed history.
    pub radial: Option<&'a RadialWall>,
}

/// Contact outcome of a solve, for reporting and verification.
#[derive(Debug, Clone)]
pub struct ContactSummary {
    /// The active set at the converged configuration.
    pub active_constraints: Vec<crate::contact::Constraint>,
    /// Largest penetration depth at the converged configuration.
    pub max_penetration: f64,
    /// Total contact reaction magnitude at the converged configuration.
    pub total_reaction: f64,
    /// How many active nodes had saturated at the Coulomb bound, or `None` when
    /// the solve ran without friction.
    pub slipping_nodes: Option<usize>,
    /// The radial wall's outcome, or `None` when no wall was configured.
    pub wall: Option<WallSummary>,
}

/// The converged result of a static solve.
#[derive(Debug, Clone)]
pub struct SolveResult {
    /// Converged nodal displacement vector.
    pub displacement: Vec<f64>,
    /// Free-DOF residual norm at the returned displacement.
    pub residual_norm: f64,
    /// Number of Newton iterations taken.
    pub newton_iterations: usize,
    /// Contact outcome, or `None` when the solve had no contact configured.
    pub contact: Option<ContactSummary>,
}

/// The active-set contact terms `(load, stiffness, active)` at configuration `u`.
///
/// Passing `penalty_contact` an empty base matrix and a zero load vector makes
/// its return exactly the two pieces of the linearised penalty: a diagonal
/// stiffness of `kappa` per active DOF, and a load of `kappa * lower` per active
/// DOF — so the active-set bookkeeping is the substrate's, not a
/// re-derivation of it here.
pub(crate) fn contact_terms<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    pairing: &ContactPairing,
    u: &[f64],
    penalty: f64,
) -> Result<(Vec<f64>, Coo, Vec<crate::contact::Constraint>), MeshError> {
    let active = pairing.active_constraints(mesh, u)?;
    let n = mesh.dof_count();
    let (stiffness, load) =
        tpt_fem_contact::penalty_contact(&Coo::new(), &vec![0.0; n], &active, penalty);
    Ok((load, stiffness, active))
}

/// Assembles the residual `f_int(u) - f_ext + r_contact(u)` for a configuration.
fn residual_vector<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &dyn Constitutive,
    load: &[f64],
    contact: Option<ContactConfig<'_>>,
    assembly: &AssemblyOptions,
    u: &[f64],
) -> Result<Vec<f64>, MeshError> {
    let mut r = internal_force(mesh, model, u, assembly)?;
    for (i, v) in load.iter().enumerate() {
        r[i] -= v;
    }
    if let Some(cfg) = contact {
        // The penalty contact residual is `K_c u - f_c`, and both halves matter:
        // `penalty_contact` returns a diagonal stiffness `kappa` *and* a load
        // `kappa * lower` per active DOF, so the residual has to carry
        // `+kappa * u_dof` as well as `-kappa * lower`. Dropping the first
        // leaves a residual the Jacobian does not differentiate, which shows up
        // as a solve that creeps down by a fraction of a percent per iteration
        // instead of converging — the penalty acts on the Jacobian alone.
        let (f_c, k_c, _) = contact_terms(mesh, cfg.pairing, u, cfg.penalty)?;
        for (i, v) in f_c.iter().enumerate() {
            r[i] -= v;
        }
        for i in 0..k_c.len() {
            r[k_c.rows[i]] += k_c.vals[i] * u[k_c.cols[i]];
        }
        // Friction enters as a plain external force on the slave nodes, so it
        // is subtracted from the residual like any other applied load — it is
        // not a stiffness acting on `u`. Its own linearisation goes into the
        // Jacobian separately, below.
        if let Some(fcfg) = cfg.friction {
            let friction = friction_terms(mesh, cfg.pairing, u, cfg.penalty, fcfg)?;
            for (i, v) in friction.force.iter().enumerate() {
                r[i] -= v;
            }
        }
        if let Some(wall) = cfg.radial {
            let t = wall.terms(mesh, u, cfg.penalty)?;
            for (i, v) in t.force.iter().enumerate() {
                r[i] += v;
            }
        }
    }
    Ok(r)
}

/// Assembles `df_int/du` plus the active contact stiffness.
fn jacobian_matrix<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &dyn Constitutive,
    contact: Option<ContactConfig<'_>>,
    assembly: &AssemblyOptions,
    u: &[f64],
) -> Result<Coo, MeshError> {
    let mut k = tangent_stiffness(mesh, model, u, assembly)?;
    if let Some(cfg) = contact {
        let (_, k_c, _) = contact_terms(mesh, cfg.pairing, u, cfg.penalty)?;
        for i in 0..k_c.len() {
            k.push(k_c.rows[i], k_c.cols[i], k_c.vals[i]);
        }
        if let Some(fcfg) = cfg.friction {
            let friction = friction_terms(mesh, cfg.pairing, u, cfg.penalty, fcfg)?;
            for i in 0..friction.tangent.len() {
                k.push(
                    friction.tangent.rows[i],
                    friction.tangent.cols[i],
                    friction.tangent.vals[i],
                );
            }
        }
        if let Some(wall) = cfg.radial {
            let t = wall.terms(mesh, u, cfg.penalty)?;
            for i in 0..t.tangent.len() {
                k.push(t.tangent.rows[i], t.tangent.cols[i], t.tangent.vals[i]);
            }
        }
    }
    Ok(k)
}

/// Solves `f_int(u) = f_ext` with optional unilateral contact, by Newton.
///
/// `dirichlet` is a list of `(dof, value)` essential conditions, condensed out
/// of every linear solve and held fixed across iterations. The initial guess is
/// zero — for a large-deflection *load*-controlled path use
/// [`crate::solve_load_path`] instead, which seeds each increment from the last
/// and bisects a step that will not converge. For displacement control this
/// function is exact as it stands.
///
/// # Errors
///
/// [`SolveError::LoadSizeMismatch`] if `load` is not `3 * node_count` long,
/// [`SolveError::Mesh`] / [`SolveError::Contact`] if an assembly or contact
/// evaluation fails, [`SolveError::Singular`] if the condensed system is
/// singular, and [`SolveError::NotConverged`] if the free-DOF residual does not
/// reach `opts.convergence`.
pub fn solve_static<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &dyn Constitutive,
    load: &[f64],
    dirichlet: &[(usize, f64)],
    opts: &SolveOptions,
    contact: Option<ContactConfig<'_>>,
) -> Result<SolveResult, SolveError> {
    newton_from(mesh, model, load, dirichlet, opts, contact, None)
}

/// The Newton loop itself, starting from `initial` when given.
///
/// `solve_static` is this with no initial guess, which is why it is defined
/// separately rather than inlined: the load-path driver in [`solve_load_path`]
/// needs the identical iteration seeded from the previous converged increment,
/// and duplicating the loop to get that would be a way for the two to drift.
///
/// `initial` is validated against the Dirichlet set: any prescribed DOF is
/// overwritten by its Dirichlet value regardless of what `initial` holds, so a
/// caller cannot accidentally carry a stale essential condition forward.
///
/// `pub(crate)` rather than private because the continuation driver in
/// [`crate::loadpath`] must reuse this exact loop; a second copy would be free
/// to drift from the one verified against the closed forms.
pub(crate) fn newton_from<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &dyn Constitutive,
    load: &[f64],
    dirichlet: &[(usize, f64)],
    opts: &SolveOptions,
    contact: Option<ContactConfig<'_>>,
    initial: Option<&[f64]>,
) -> Result<SolveResult, SolveError> {
    let n = mesh.dof_count();
    if load.len() != n {
        return Err(SolveError::LoadSizeMismatch {
            expected: n,
            found: load.len(),
        });
    }
    if let Some(u0) = initial {
        if u0.len() != n {
            return Err(SolveError::LoadSizeMismatch {
                expected: n,
                found: u0.len(),
            });
        }
    }
    // Validate the pairing once, outside the iteration, so a bad pairing is
    // reported as an error rather than swallowed mid-Newton.
    if let Some(cfg) = contact {
        cfg.pairing
            .candidates(mesh, &vec![0.0; n])
            .map_err(SolveError::Contact)?;
    }

    let fixed: HashSet<usize> = dirichlet.iter().map(|(i, _)| *i).collect();
    let free: Vec<usize> = (0..n).filter(|i| !fixed.contains(i)).collect();
    let free_index: std::collections::HashMap<usize, usize> =
        free.iter().enumerate().map(|(k, &v)| (v, k)).collect();
    let load_norm = load.iter().map(|x| x * x).sum::<f64>().sqrt();
    let tolerance = opts.convergence.abs_tol + opts.convergence.rel_tol * load_norm;

    let mut u = match initial {
        Some(u0) => u0.to_vec(),
        None => vec![0.0; n],
    };
    for (dof, value) in dirichlet {
        u[*dof] = *value;
    }

    let mut iterations = 0usize;
    let mut residual_norm = f64::INFINITY;
    let jacobian_calls = Cell::new(0usize);
    for _ in 0..opts.convergence.max_iter {
        iterations += 1;
        let r = residual_vector(mesh, model, load, contact, &opts.assembly, &u)?;
        let r_free: Vec<f64> = free.iter().map(|&i| r[i]).collect();
        residual_norm = r_free.iter().map(|x| x * x).sum::<f64>().sqrt();
        if residual_norm <= tolerance {
            break;
        }
        jacobian_calls.set(jacobian_calls.get() + 1);
        let k = jacobian_matrix(mesh, model, contact, &opts.assembly, &u)?;
        // Condense the essential DOFs out of the global COO.
        let csr = k.to_csr();
        let mut condensed = Coo::new();
        for &row in &free {
            for c in csr.row_ptrs[row]..csr.row_ptrs[row + 1] {
                if let Some(&j) = free_index.get(&csr.col_ind[c]) {
                    condensed.push(free_index[&row], j, csr.values[c]);
                }
            }
        }
        // Diagonal equilibration before the linear solve. The substrate's
        // dense backend factors with partial pivoting and an *absolute*
        // singularity threshold (1e-12), which a matrix whose entries span four
        // orders of magnitude — a volumetric penalty diagonal next to a shear
        // off-diagonal — can trip on size alone, well before it is actually
        // ill-conditioned. Scaling to a unit diagonal puts every pivot on an
        // O(1) footing so that threshold means what it says. The scaling is
        // undone on the way out, so the step is exact.
        let mut scale = vec![1.0f64; free.len()];
        for k in 0..free.len() {
            let diag: f64 = condensed
                .rows
                .iter()
                .zip(&condensed.cols)
                .zip(&condensed.vals)
                .filter(|((&r, &c), _)| r == k && c == k)
                .map(|(_, &v)| v)
                .sum();
            if diag > 0.0 && diag.is_finite() {
                scale[k] = 1.0 / diag.sqrt();
            }
        }
        let mut scaled = Coo::with_capacity(condensed.len());
        for i in 0..condensed.len() {
            let (r, c) = (condensed.rows[i], condensed.cols[i]);
            scaled.push(r, c, condensed.vals[i] * scale[r] * scale[c]);
        }
        let rhs: Vec<f64> = r_free.iter().zip(&scale).map(|(v, s)| v * s).collect();
        let y = solve(&scaled, &rhs).map_err(|_| SolveError::Singular)?;
        let step: Vec<f64> = y.iter().zip(&scale).map(|(v, s)| v * s).collect();
        // Damped Newton. A hyperelastic tangent is only guaranteed to be
        // positive definite in the region where the material is stable, so a
        // full step from a poor initial guess (a zero start under a large
        // prescribed displacement) can leave the admissible region entirely and
        // invert an element. Halving until the free-DOF residual actually
        // decreases is the standard guard, and it costs one extra residual
        // evaluation per rejection.
        let mut alpha = 1.0f64;
        let mut accepted = false;
        for _ in 0..MAX_HALVINGS {
            let mut trial = u.clone();
            for (k, &dof) in free.iter().enumerate() {
                trial[dof] -= alpha * step[k];
            }
            match residual_vector(mesh, model, load, contact, &opts.assembly, &trial) {
                Ok(r) if r.iter().all(|v| v.is_finite()) => {
                    let norm: f64 = free.iter().map(|&i| r[i] * r[i]).sum::<f64>().sqrt();
                    if norm < residual_norm {
                        u = trial;
                        accepted = true;
                        break;
                    }
                }
                _ => {}
            }
            alpha *= 0.5;
        }
        if !accepted {
            // No damping rescued the step: take the full one and let the
            // convergence test at the top of the next iteration decide, so the
            // caller gets a NotConverged with a real residual rather than a
            // silent early exit.
            for (k, &dof) in free.iter().enumerate() {
                u[dof] -= step[k];
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
    let _ = jacobian_calls.get();

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
            wall: match cfg.radial {
                Some(w) => Some(
                    w.summary(mesh, &u, cfg.penalty)
                        .map_err(SolveError::Contact)?,
                ),
                None => None,
            },
        }),
        None => None,
    };
    Ok(SolveResult {
        displacement: u,
        residual_norm,
        newton_iterations: iterations,
        contact,
    })
}

/// The residual `f_int(u) - f_ext + r_contact(u)` at a configuration.
///
/// Exposed so a caller (and this crate's own tests) can check a converged
/// solution independently of the solve, and so the contact contribution can be
/// differenced against the Jacobian the solver uses.
///
/// # Errors
///
/// As [`solve_static`]'s assembly and contact evaluation.
pub fn residual<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &dyn Constitutive,
    load: &[f64],
    u: &[f64],
    opts: &SolveOptions,
    contact: Option<ContactConfig<'_>>,
) -> Result<Vec<f64>, SolveError> {
    residual_vector(mesh, model, load, contact, &opts.assembly, u).map_err(SolveError::from)
}
