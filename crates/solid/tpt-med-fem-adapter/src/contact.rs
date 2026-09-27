//! Unilateral contact between a deforming `Hex8` surface and a rigid obstacle,
//! layered onto the nonlinear Newton solve in [`crate::solver`].
//!
//! The pairing of slave nodes to obstacle points is
//! `tpt_fem_contact::contact_pairs` (octree-accelerated nearest-node search,
//! reused rather than reimplemented). The enforcement is
//! `tpt_fem_contact::penalty_contact`, whose `(stiffness, load)` pair is
//! exactly the linearisation of the penalty energy this adapter adds at each
//! Newton iteration — see [`crate::solver`] for the sign convention.
//!
//! # What the substrate does and does not give us
//!
//! `tpt-fem-contact` at 0.1.0 constrains a single global DOF to stay at or
//! above a scalar `lower` bound — exactly the form needed for a planar,
//! normal-aligned obstacle, which is what this adapter models. Two
//! consequences are named rather than hidden:
//!
//! - `contact_pairs` returns a **non-negative** nearest-node distance, which
//!   cannot itself express penetration. This adapter therefore uses the
//!   pairing only to choose the contact partner, and computes the signed
//!   normal gap itself as `x_slave[axis] - x_master[axis]`. A constraint is
//!   active when that signed gap is negative.
//! - There is **no friction** in the pinned 0.1.0 substrate source (RFC 0009's
//!   second unresolved question). Friction is added by this workspace in
//!   [`crate::friction`], as a regularized Coulomb layer on top of the
//!   primitives here, rather than upstream.
//!
//! # Re-evaluation, not a frozen active set
//!
//! [`ContactPairing::active_constraints`] is called from inside the Newton
//! residual and Jacobian closures, so the active set is recomputed from the
//! *current* geometry at every residual evaluation and every Jacobian
//! evaluation. That is the specific coupling design RFC 0009 listed as
//! unresolved: a fixed active set computed at the undeformed configuration is
//! simply wrong once the geometry moves, and freezing it makes a body that
//! lifts off a wall keep pushing against it.

use crate::mesh::{Mesh, MeshError};
use tpt_fem_contact::ContactConstraint;
use tpt_fem_element::ReferenceElement;
use tpt_med_geometry::Vec3;

pub use tpt_fem_contact::ContactConstraint as Constraint;

/// Errors produced when defining a contact pairing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContactError {
    /// The normal axis was not `0`, `1` or `2`.
    AxisOutOfRange(usize),
    /// The obstacle had no points, so no slave node could ever be paired.
    EmptyMaster,
}

impl std::fmt::Display for ContactError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ContactError::AxisOutOfRange(a) => {
                write!(
                    f,
                    "contact normal axis {a} is not one of 0 (x), 1 (y), 2 (z)"
                )
            }
            ContactError::EmptyMaster => write!(f, "the contact obstacle has no points"),
        }
    }
}

/// A paired slave node, the obstacle point it is closest to, and the signed
/// normal gap between them.
///
/// No `PartialEq` derive: the substrate's `ContactConstraint` does not
/// implement it, and the fields are public for callers that want to compare
/// them themselves.
#[derive(Debug, Clone, Copy)]
pub struct ContactCandidate {
    /// The constraint this pairing would impose (inactive unless the gap is
    /// negative).
    pub constraint: ContactConstraint,
    /// Signed gap along the normal axis: `x_slave[axis] - x_master[axis]`.
    /// Negative means penetration.
    pub gap: f64,
    /// Index into the obstacle point list of the paired master point.
    pub master_index: usize,
    /// The slave node this candidate belongs to.
    ///
    /// The substrate's `ContactConstraint` names only the DOF, which does not
    /// identify the node, so it is recorded here for callers that need the
    /// node's full displacement — the frictional layer needs the tangential
    /// components, not just the normal one.
    pub node: usize,
}

/// A frictionless unilateral contact pairing between deforming slave nodes and
/// a rigid obstacle, resolved along a single coordinate axis.
///
/// The pairing itself is purely geometric — it knows nothing about friction.
/// Friction is layered on top of it in [`crate::friction`], so that the normal
/// and tangential problems stay separable and the normal path is unchanged for
/// callers that do not ask for friction.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactPairing {
    axis: usize,
    slave: Vec<usize>,
    master: Vec<Vec3>,
    activation_tolerance: f64,
}

impl ContactPairing {
    /// Builds a pairing of `slave` node indices against `master` obstacle
    /// points, with contact resolved along `axis`.
    ///
    /// # Errors
    ///
    /// [`ContactError::AxisOutOfRange`] if `axis > 2`, or
    /// [`ContactError::EmptyMaster`] if the obstacle has no points.
    pub fn new(
        axis: usize,
        slave: impl IntoIterator<Item = usize>,
        master: impl IntoIterator<Item = Vec3>,
    ) -> Result<Self, ContactError> {
        if axis > 2 {
            return Err(ContactError::AxisOutOfRange(axis));
        }
        let master: Vec<Vec3> = master.into_iter().collect();
        if master.is_empty() {
            return Err(ContactError::EmptyMaster);
        }
        Ok(Self {
            axis,
            slave: slave.into_iter().collect(),
            master,
            activation_tolerance: 0.0,
        })
    }

    /// Sets the gap at or below which a contact is considered active.
    ///
    /// The default of `0.0` means *touching counts as active*, which is the
    /// choice that keeps the first Newton step well posed: a body resting
    /// exactly on the obstacle has no penetration to activate a constraint, yet
    /// without one its rigid-body mode along the normal is unrestrained and the
    /// first linear system is singular. Raising the value activates contact
    /// slightly before touching.
    pub fn with_activation_tolerance(mut self, activation_tolerance: f64) -> Self {
        self.activation_tolerance = activation_tolerance;
        self
    }

    /// The contact normal axis (`0` = x, `1` = y, `2` = z).
    pub fn axis(&self) -> usize {
        self.axis
    }

    /// Slave node indices.
    pub fn slave_nodes(&self) -> &[usize] {
        &self.slave
    }

    /// Rigid obstacle points.
    pub fn master_points(&self) -> &[Vec3] {
        &self.master
    }

    /// The gap at or below which a contact counts as active.
    ///
    /// Exposed so the frictional layer in [`crate::friction`] applies exactly
    /// the same activity test as the normal contact, rather than a second
    /// definition of "in contact" that could drift from it.
    pub fn activation_tolerance(&self) -> f64 {
        self.activation_tolerance
    }
}

impl std::error::Error for ContactError {}

impl ContactPairing {
    /// Every slave node paired to its nearest obstacle point, with its signed
    /// normal gap. One entry per slave node, in `slave_nodes()` order.
    ///
    /// # Errors
    ///
    /// [`MeshError::DofCountMismatch`] if `u` is not `3 * node_count` long, or
    /// [`MeshError::NodeIndexOutOfRange`] if a slave node is not in the mesh.
    pub fn candidates<E: ReferenceElement + crate::mesh::ElementFamily>(
        &self,
        mesh: &Mesh<E>,
        u: &[f64],
    ) -> Result<Vec<ContactCandidate>, MeshError> {
        let x = mesh.positions(u)?;
        for &node in &self.slave {
            if node >= mesh.node_count() {
                return Err(MeshError::NodeIndexOutOfRange {
                    element: usize::MAX,
                    local_node: node,
                    node,
                    node_count: mesh.node_count(),
                });
            }
        }
        let slave_points: Vec<(usize, Vec<f64>)> = self
            .slave
            .iter()
            .map(|&n| (n, x[n].to_array().to_vec()))
            .collect();
        let master_indexed: Vec<(usize, Vec<f64>)> = self
            .master
            .iter()
            .enumerate()
            .map(|(i, p)| (i, p.to_array().to_vec()))
            .collect();
        let pairs = tpt_fem_contact::contact_pairs(&slave_points, &master_indexed);
        let mut out = Vec::with_capacity(self.slave.len());
        for ((node, _), (_, nearest)) in slave_points.iter().zip(pairs) {
            // `contact_pairs` returns `None` only for an empty master list,
            // which `new` already rejected; the arm is kept rather than
            // indexing a sentinel that does not exist.
            let Some((master_index, _distance)) = nearest else {
                continue;
            };
            let lower = self.master[master_index].to_array()[self.axis];
            out.push(ContactCandidate {
                constraint: ContactConstraint {
                    dof: mesh.dof(*node, self.axis),
                    lower,
                },
                gap: x[*node].to_array()[self.axis] - lower,
                master_index,
                node: *node,
            });
        }
        Ok(out)
    }

    /// The active (penetrating) constraints at configuration `u`.
    ///
    /// # Errors
    ///
    /// As [`ContactPairing::candidates`].
    pub fn active_constraints<E: ReferenceElement + crate::mesh::ElementFamily>(
        &self,
        mesh: &Mesh<E>,
        u: &[f64],
    ) -> Result<Vec<ContactConstraint>, MeshError> {
        Ok(self
            .candidates(mesh, u)?
            .into_iter()
            .filter(|c| c.gap <= self.activation_tolerance)
            .map(|c| c.constraint)
            .collect())
    }

    /// Largest penetration depth over all slave nodes (`0.0` when clear).
    ///
    /// # Errors
    ///
    /// As [`ContactPairing::candidates`].
    pub fn max_penetration<E: ReferenceElement + crate::mesh::ElementFamily>(
        &self,
        mesh: &Mesh<E>,
        u: &[f64],
    ) -> Result<f64, MeshError> {
        Ok(self
            .candidates(mesh, u)?
            .iter()
            .map(|c| -c.gap)
            .fold(0.0f64, f64::max))
    }

    /// Total normal contact reaction magnitude
    /// `sum_i penalty * max(0, lower_i - x_i)` over the slave nodes: the force
    /// the obstacle exerts, positive *resisting* penetration (so it balances
    /// the magnitude of a load pressing the body into the obstacle).
    ///
    /// # Errors
    ///
    /// As [`ContactPairing::candidates`].
    pub fn total_reaction<E: ReferenceElement + crate::mesh::ElementFamily>(
        &self,
        mesh: &Mesh<E>,
        u: &[f64],
        penalty: f64,
    ) -> Result<f64, MeshError> {
        Ok(self
            .candidates(mesh, u)?
            .iter()
            .map(|c| penalty * (-c.gap).max(0.0))
            .sum())
    }
}
