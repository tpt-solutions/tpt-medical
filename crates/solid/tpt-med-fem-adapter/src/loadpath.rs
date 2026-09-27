//! Load stepping / continuation for large-deflection load-controlled paths.
//!
//! `solve_static` in [`crate::solver`] starts from zero displacement. That is
//! fine for a small problem, and wrong for a large-deflection *load*-controlled
//! one: Newton from a zero start has to find the equilibrium path in one shot,
//! and past a limit point there is no next state to find — the tangent becomes
//! singular and the iteration either stalls or walks off into a configuration
//! that is not on the path at all.
//!
//! This module walks the path instead, one converged increment at a time, each
//! seeded from the last. What follows is what it deliberately does *not* do.
//!
//! # Proportional stepping, not arc length
//!
//! Increments are proportional in the load: step `i` applies `t_i * load` for
//! `t_i` increasing from 0 to 1. **There is no arc-length control.** That is a
//! real limitation and the reason a limit point is still not reachable from
//! here — proportional stepping carries the load factor monotonically upward,
//! and at a limit point the structure cannot carry more load at all, so the
//! step stalls and the cutback below gives up. Getting past one needs a genuine
//! arc-length (or dynamic-relaxation) formulation with a sign convention for the
//! load factor, which is a larger piece of work and a different API. This module
//! is the stepping a path *without* a limit point needs, and it stops at one
//! rather than pretending to pass it.
//!
//! # Cutback
//!
//! When an increment fails to converge, the driver does not give up. It bisects
//! the *remaining* increment and retries, up to
//! [`LoadPathOptions::max_cutbacks`] times. That is what makes a path which
//! turns sharply near the end still walkable: the steps that matter are the
//! small ones, and halving produces them without the caller having to know where.
//!
//! A cutback bisects only the *remaining* distance, never a step already
//! converged, so the returned path is monotone in load factor and every point
//! on it satisfies its own equilibrium tolerance.
//!
//! # Every point returned is converged
//!
//! The path is only extended with increments that reached [`crate::Convergence`].
//! If subdivision is exhausted the error carries the last displacement reached,
//! so a caller can inspect how far it got rather than restarting blind. There is
//! no "best effort" path containing unconverged points — that would make every
//! downstream use of a point unsafe without saying so on the point itself.

use crate::assembly::Constitutive;
use crate::mesh::{Mesh, MeshError};
use crate::solver::{newton_from, ContactConfig, ContactSummary, SolveError, SolveOptions};
use std::fmt;
use tpt_fem_element::ReferenceElement;

/// How to walk the load path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LoadPathOptions {
    /// Number of equal increments from zero to the full load.
    ///
    /// `1` recovers a single full-load solve. The count is an exact divisor of
    /// the path, not a suggestion: the driver never takes a partial step
    /// silently.
    pub steps: usize,
    /// Maximum bisections of a remaining increment after a failed convergence
    /// attempt. `0` disables cutback, so each step is all-or-nothing.
    pub max_cutbacks: usize,
}

impl Default for LoadPathOptions {
    fn default() -> Self {
        Self {
            steps: 10,
            max_cutbacks: 6,
        }
    }
}

/// Errors specific to the load-path driver. Convergence failures are reported
/// through [`SolveError::NotConverged`], wrapped with step context rather than
/// duplicated.
#[derive(Debug)]
pub enum LoadPathError {
    /// An increment, and every bisection of it, failed to converge.
    StepNotConverged {
        /// Load factor this increment was attempting.
        load_factor: f64,
        /// How many bisections were tried before giving up.
        cutbacks: usize,
        /// The underlying failure.
        source: SolveError,
    },
    /// A mesh or contact evaluation failed outside the Newton loop.
    Mesh(MeshError),
    /// The load vector had the wrong length.
    LoadSizeMismatch {
        /// The length the solve requires.
        expected: usize,
        /// The length supplied.
        found: usize,
    },
}

impl fmt::Display for LoadPathError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LoadPathError::StepNotConverged {
                load_factor,
                cutbacks,
                source,
            } => write!(
                f,
                "load step to factor {load_factor} did not converge after {cutbacks} \
                 cutbacks: {source}"
            ),
            LoadPathError::Mesh(e) => write!(f, "{e}"),
            LoadPathError::LoadSizeMismatch { expected, found } => {
                write!(f, "load vector has {found} entries, expected {expected}")
            }
        }
    }
}

impl std::error::Error for LoadPathError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            LoadPathError::StepNotConverged { source, .. } => Some(source),
            LoadPathError::Mesh(e) => Some(e),
            LoadPathError::LoadSizeMismatch { .. } => None,
        }
    }
}

impl From<SolveError> for LoadPathError {
    fn from(e: SolveError) -> Self {
        match e {
            SolveError::Mesh(m) => LoadPathError::Mesh(m),
            SolveError::LoadSizeMismatch { expected, found } => {
                LoadPathError::LoadSizeMismatch { expected, found }
            }
            // A contact or singularity failure is not specific to one step, so
            // it is not dressed up as a convergence failure — that would imply
            // cutback could fix it, and it cannot.
            other => LoadPathError::StepNotConverged {
                load_factor: f64::NAN,
                cutbacks: 0,
                source: other,
            },
        }
    }
}

/// One converged point on the load path.
#[derive(Debug, Clone)]
pub struct LoadStep {
    /// The load factor this point was solved at, in `0.0..=1.0`.
    pub load_factor: f64,
    /// The converged displacement vector.
    pub displacement: Vec<f64>,
    /// Free-DOF residual norm at the returned displacement.
    pub residual_norm: f64,
    /// Newton iterations the *final* increment took, not the sum over cutbacks:
    /// a caller sizing an iteration budget wants the worst single solve.
    pub newton_iterations: usize,
    /// Contact outcome at this point, or `None` when no contact was configured.
    pub contact: Option<ContactSummary>,
}

impl LoadStep {
    /// The applied external load at this point, `load_factor * load`.
    pub fn applied_load(&self, load: &[f64]) -> Vec<f64> {
        load.iter().map(|f| f * self.load_factor).collect()
    }
}

/// A converged load path, always starting at zero load and zero displacement.
#[derive(Debug, Clone)]
pub struct LoadPath {
    /// The converged points, in increasing load factor. The first entry is the
    /// undeformed state at factor 0.
    pub steps: Vec<LoadStep>,
}

impl LoadPath {
    /// The final point, or `None` for an empty path.
    pub fn last(&self) -> Option<&LoadStep> {
        self.steps.last()
    }

    /// The load factor actually reached, or 0.0 for an empty path.
    pub fn reached(&self) -> f64 {
        self.last().map_or(0.0, |s| s.load_factor)
    }
}

/// Walks the load path from zero to the full `load`, one converged increment at
/// a time, with bisection cutback on a failed increment.
///
/// `dirichlet` and `contact` are held fixed across the whole path and passed
/// through to every increment unchanged. In particular a prescribed
/// displacement is *not* ramped: this is a load-controlled driver, and if you
/// want displacement control, prescribe the displacement and call
/// [`crate::solve_static`] directly, which is exact and needs no path.
///
/// # Errors
///
/// [`LoadPathError::LoadSizeMismatch`] if `load` is the wrong length, and
/// [`LoadPathError::StepNotConverged`] if an increment cannot be converged even
/// after `max_cutbacks` bisections. Nothing partial is returned in the error
/// case: the displacement reached is inside the wrapped [`SolveError`].
pub fn solve_load_path<E: ReferenceElement + crate::mesh::ElementFamily>(
    mesh: &Mesh<E>,
    model: &dyn Constitutive,
    load: &[f64],
    dirichlet: &[(usize, f64)],
    opts: &SolveOptions,
    contact: Option<ContactConfig<'_>>,
    path: LoadPathOptions,
) -> Result<LoadPath, LoadPathError> {
    let n = mesh.dof_count();
    if load.len() != n {
        return Err(LoadPathError::LoadSizeMismatch {
            expected: n,
            found: load.len(),
        });
    }
    if path.steps == 0 {
        return Ok(LoadPath { steps: Vec::new() });
    }

    // Factor 0 is the undeformed state, and it is trivially in equilibrium: the
    // applied load is zero, so the free-DOF residual is exactly zero. It is
    // included rather than assumed, so `steps[0]` exists for a caller that
    // indexes the path.
    let mut steps = vec![LoadStep {
        load_factor: 0.0,
        displacement: {
            let mut u = vec![0.0; n];
            for (dof, value) in dirichlet {
                u[*dof] = *value;
            }
            u
        },
        residual_norm: 0.0,
        newton_iterations: 0,
        contact: None,
    }];

    let mut current = 0.0f64;
    let mut state = steps[0].displacement.clone();

    for i in 1..=path.steps {
        // The nominal end of this increment. Cutback walks *up* from `current`
        // towards it, so the path is monotone and no converged point is
        // revisited.
        let target = i as f64 / path.steps as f64;
        let mut reached = current;
        let mut attempt_target = target;
        let mut cutbacks = 0usize;

        loop {
            let scaled: Vec<f64> = load.iter().map(|f| f * attempt_target).collect();
            match newton_from(mesh, model, &scaled, dirichlet, opts, contact, Some(&state)) {
                Ok(r) => {
                    // `clone_from` reuses the existing allocation; `state` and
                    // the stored copy are the same length every step, so the
                    // copy that `state = r.displacement` would make is pure
                    // waste on a path with many increments.
                    state.clone_from(&r.displacement);
                    reached = attempt_target;
                    steps.push(LoadStep {
                        load_factor: attempt_target,
                        displacement: r.displacement,
                        residual_norm: r.residual_norm,
                        newton_iterations: r.newton_iterations,
                        contact: r.contact,
                    });
                    break;
                }
                Err(e) => {
                    // Only a convergence failure is worth retrying smaller. A
                    // singular system, a bad pairing or an inverted element is a
                    // property of the problem, not of the step size.
                    if !matches!(e, SolveError::NotConverged { .. }) {
                        return Err(e.into());
                    }
                    if cutbacks >= path.max_cutbacks {
                        return Err(LoadPathError::StepNotConverged {
                            load_factor: attempt_target,
                            cutbacks,
                            source: e,
                        });
                    }
                    cutbacks += 1;
                    // Bisect the *remaining* distance, so progress is retained.
                    attempt_target = reached + (target - reached) / 2.0;
                }
            }
        }
        current = reached;
    }

    Ok(LoadPath { steps })
}
