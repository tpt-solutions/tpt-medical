//! # Tangential basis
//!
//! Friction acts in the plane orthogonal to the contact normal, which for this
//! crate's axis-aligned obstacle is exactly the two axes that are not
//! [`ContactPairing::axis`]. Slip is the tangential part of the separation from
//! the paired master point, so a node whose obstacle point starts at the same
//! tangential position has zero slip until it moves. No tangent basis needs
//! storing or rotating, because the normal is a coordinate axis by
//! construction.

use crate::contact::ContactPairing;
use crate::mesh::{HexMesh, MeshError};
use tpt_fem_element::ReferenceElement;
use tpt_fem_sparse::Coo;

/// Friction parameters for [`friction_terms`].
///
/// `mu = 0.0` is the frictionless case exactly, not approximately: every
/// tangential force and tangent is multiplied by `mu`, so a zero coefficient
/// reproduces normal-contact-only behaviour bit for bit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrictionConfig {
    /// Coulomb friction coefficient. Must be non-negative.
    pub mu: f64,
    /// Regularizing tangential stiffness `k_t` per active node, in the same
    /// units as the normal penalty.
    ///
    /// The stiffness of the "microslip" spring standing in for the stick
    /// branch. It must be large relative to the normal penalty for the
    /// regularization length `mu * f_n / k_t` to be small.
    pub tangential_stiffness: f64,
}

impl FrictionConfig {
    /// A friction configuration, rejecting a non-finite or negative parameter
    /// up front rather than producing a non-physical attraction later.
    ///
    /// # Errors
    ///
    /// [`FrictionError::InvalidMu`] or
    /// [`FrictionError::InvalidTangentialStiffness`].
    pub fn new(mu: f64, tangential_stiffness: f64) -> Result<Self, FrictionError> {
        if !mu.is_finite() || mu < 0.0 {
            return Err(FrictionError::InvalidMu(mu));
        }
        if !tangential_stiffness.is_finite() || tangential_stiffness < 0.0 {
            return Err(FrictionError::InvalidTangentialStiffness(
                tangential_stiffness,
            ));
        }
        Ok(Self {
            mu,
            tangential_stiffness,
        })
    }
}

/// Errors produced when configuring friction.
#[derive(Debug, Clone, PartialEq)]
pub enum FrictionError {
    /// The Coulomb coefficient was negative, infinite or NaN.
    InvalidMu(f64),
    /// The regularizing tangential stiffness was negative, infinite or NaN.
    InvalidTangentialStiffness(f64),
}

impl std::fmt::Display for FrictionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FrictionError::InvalidMu(mu) => {
                write!(
                    f,
                    "friction coefficient {mu} must be finite and non-negative"
                )
            }
            FrictionError::InvalidTangentialStiffness(k) => write!(
                f,
                "tangential friction stiffness {k} must be finite and non-negative"
            ),
        }
    }
}

impl std::error::Error for FrictionError {}

/// The friction state of one active node at one configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrictionState {
    /// The two tangential axes (everything except the contact normal).
    pub tangential_axes: [usize; 2],
    /// Tangential slip vector, in `tangential_axes` order.
    pub slip: [f64; 2],
    /// Its magnitude; zero when the node has not moved tangentially.
    pub slip_norm: f64,
    /// `true` when the force has saturated at `mu * f_n` (slipping), `false`
    /// when still on the `k_t * s` branch (sticking).
    pub slipping: bool,
    /// The force applied *to the slave node*, opposing `slip`, in
    /// `tangential_axes` order.
    pub force: [f64; 2],
    /// The `mu * f_n` bound this node is limited to.
    pub bound: f64,
}

/// Friction force and tangent for every currently active contact node.
pub struct FrictionTerms {
    /// One entry per active slave node.
    pub states: Vec<FrictionState>,
    /// Friction force on each node, in global DOF order. Only the tangential
    /// DOFs of active nodes are ever non-zero.
    pub force: Vec<f64>,
    /// The exact derivative `d(force)/d(u)` of [`FrictionTerms::force`].
    pub tangent: Coo,
    /// How many active nodes are saturated at the Coulomb bound.
    pub slipping_nodes: usize,
}

/// Evaluates regularized Coulomb friction for every currently active contact.
///
/// `penalty` is the *normal* penalty stiffness from [`crate::ContactConfig`],
/// needed because the Coulomb bound is `mu * f_n` and the normal penalty
/// reaction is `f_n = penalty * max(0, -gap)`.
///
/// A node is active — and so carries friction — on exactly the same test the
/// normal contact uses, [`ContactPairing::activation_tolerance`], so the two
/// layers can never disagree about who is touching.
///
/// # Errors
///
/// As [`ContactPairing::candidates`]: a DOF-count mismatch, or a slave node
/// outside the mesh.
pub fn friction_terms<E: ReferenceElement>(
    mesh: &HexMesh<E>,
    pairing: &ContactPairing,
    u: &[f64],
    penalty: f64,
    cfg: FrictionConfig,
) -> Result<FrictionTerms, MeshError> {
    let n = mesh.dof_count();
    let mut states = Vec::new();
    let mut force = vec![0.0; n];
    let mut tangent = Coo::new();
    let mut slipping_nodes = 0usize;

    // A zero coefficient or a zero regularizing stiffness is the frictionless
    // case *exactly* — no tangential force, zero tangent — rather than a `min`
    // that happens to evaluate to zero. Returning before any geometry work also
    // means a frictionless caller pays nothing.
    if cfg.mu == 0.0 || cfg.tangential_stiffness == 0.0 {
        return Ok(FrictionTerms {
            states,
            force,
            tangent,
            slipping_nodes,
        });
    }

    let normal = pairing.axis();
    let tangential_axes = match normal {
        0 => [1usize, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    let x = mesh.positions(u)?;
    let candidates = pairing.candidates(mesh, u)?;

    for c in &candidates {
        // A node that is separating is not in contact and must not be dragged.
        if c.gap > pairing.activation_tolerance() {
            continue;
        }
        let node = c.node;
        let p = x[node].to_array();
        let m = pairing.master_points()[c.master_index].to_array();

        // Tangential part of the separation from the paired obstacle point. The
        // normal component is excluded, so penetration itself never contributes
        // to slip and the friction force is genuinely tangential.
        let slip = [
            p[tangential_axes[0]] - m[tangential_axes[0]],
            p[tangential_axes[1]] - m[tangential_axes[1]],
        ];
        let slip_norm = slip[0].hypot(slip[1]);
        let dofs = [
            mesh.dof(node, tangential_axes[0]),
            mesh.dof(node, tangential_axes[1]),
        ];

        // The Coulomb bound, from the normal penalty reaction. `max(0, ...)` so
        // a node at exactly zero penetration — which the default activation
        // tolerance counts as active, deliberately, to keep the first linear
        // system non-singular — gets no friction bound rather than a spurious
        // one.
        let bound = cfg.mu * (penalty * (-c.gap).max(0.0));

        if slip_norm == 0.0 {
            // No slip, no force. The tangent is *not* zero: the node is held,
            // so a tangential perturbation of it must be resisted. Omitting
            // this leaves a stuck node free to creep tangentially in the linear
            // solve, which is exactly the artefact friction exists to prevent.
            for &d in &dofs {
                tangent.push(d, d, cfg.tangential_stiffness);
            }
            states.push(FrictionState {
                tangential_axes,
                slip,
                slip_norm,
                slipping: false,
                force: [0.0, 0.0],
                bound,
            });
            continue;
        }

        // The regularized law: `f = -min(mu * f_n, k_t * |s|) * s / |s|`, the
        // saturated magnitude along `-s`.
        let trial = cfg.tangential_stiffness * slip_norm;
        let slipping = trial >= bound;
        if slipping {
            slipping_nodes += 1;
        }
        let magnitude = if slipping { bound } else { trial };
        let f = [
            -magnitude * slip[0] / slip_norm,
            -magnitude * slip[1] / slip_norm,
        ];
        for (i, &d) in dofs.iter().enumerate() {
            force[d] += f[i];
        }

        // Exact derivative of the force w.r.t. the two tangential DOFs, holding
        // the normal gap — and hence the bound — fixed. The neglected coupling
        // through the gap is second order in the perturbation, and the
        // finite-difference Jacobian test covers the assembled result as a
        // whole.
        //
        // - Stick branch, `f = -k_t s`: derivative `-k_t I`, symmetric.
        // - Slip branch, `f = -bound * s/|s|`: derivative
        //   `-bound * (I/|s| - s s^T / |s|^3)`. Symmetric, and singular along
        //   the slip direction, which is correct — a slipping node carries no
        //   restoring force in the direction it is sliding.
        if slipping {
            let inv = 1.0 / slip_norm;
            for a in 0..2 {
                for b in 0..2 {
                    let outer = slip[a] * slip[b] * inv * inv * inv;
                    let delta = if a == b { 1.0 } else { 0.0 };
                    tangent.push(dofs[a], dofs[b], -bound * (delta - outer));
                }
            }
        } else {
            for &d in &dofs {
                tangent.push(d, d, -cfg.tangential_stiffness);
            }
        }

        states.push(FrictionState {
            tangential_axes,
            slip,
            slip_norm,
            slipping,
            force: f,
            bound,
        });
    }

    Ok(FrictionTerms {
        states,
        force,
        tangent,
        slipping_nodes,
    })
}
