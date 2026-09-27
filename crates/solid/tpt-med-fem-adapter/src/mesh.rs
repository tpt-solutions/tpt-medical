//! A 3-D trilinear hexahedral (`Hex8`) mesh: node coordinates, elements, and
//! the structured box builder the verification tests use.
//!
//! Element connectivity follows the reference-element node ordering of
//! `tpt_fem_element::Hex8` exactly, so `elements[e][k]` is the mesh node
//! sitting at reference coordinate `Hex8::nodes()[k]`. The structured builder
//! emits that ordering, and any mesh assembled by hand must too — the
//! isoparametric mapping in [`crate::assembly`] assumes it.
//!
//! # Example
//!
//! ```
//! use tpt_med_fem_adapter::{hex_box, Hex8Mesh, MeshError};
//! use tpt_med_geometry::Vec3;
//!
//! // Unit cube split into 2x2x2 hexahedra: 27 nodes, 8 elements, 81 DOFs.
//! let mesh = hex_box(2, 2, 2, 1.0, 1.0, 1.0).expect("non-degenerate box");
//! assert_eq!(mesh.node_count(), 27);
//! assert_eq!(mesh.element_count(), 8);
//! assert_eq!(mesh.dof_count(), 81);
//! // Every element has strictly positive volume (correct orientation).
//! for e in 0..mesh.element_count() {
//!     assert!(mesh.element_volume(e).expect("upright element") > 0.0);
//! }
//! // Out-of-range connectivity is rejected rather than read out of bounds.
//! let bad = Hex8Mesh::from_parts(
//!     vec![Vec3::ZERO; 8],
//!     vec![[0, 1, 2, 3, 4, 5, 6, 99]],
//! );
//! assert!(matches!(bad, Err(MeshError::NodeIndexOutOfRange { .. })));
//! ```

use tpt_fem_element::{hex_rule, Hex8, ReferenceElement};
use tpt_med_geometry::{Mat3, Vec3};

/// Errors produced when building or querying a [`Hex8Mesh`].
#[derive(Debug, Clone, PartialEq)]
pub enum MeshError {
    /// An element referenced a node index the mesh does not have.
    NodeIndexOutOfRange {
        /// The offending element index.
        element: usize,
        /// The offending local node position within the element.
        local_node: usize,
        /// The node index that was out of range.
        node: usize,
        /// Number of nodes the mesh actually has.
        node_count: usize,
    },
    /// A displacement vector's length did not equal `3 * node_count`.
    DofCountMismatch {
        /// The length the operation requires.
        expected: usize,
        /// The length that was supplied.
        found: usize,
    },
    /// An element's Jacobian determinant was non-positive (inverted or
    /// degenerate) when a quantity requiring a valid volume was requested.
    DegenerateElement {
        /// The offending element index.
        element: usize,
        /// The Jacobian determinant that was found.
        jacobian_determinant: f64,
    },
    /// An element's deformation gradient had a non-positive determinant at an
    /// integration point: the element has inverted (or degenerated) and its
    /// deformation gradient is not a valid state for the constitutive law.
    InvertedDeformation {
        /// The offending element index.
        element: usize,
        /// The determinant det F that was found.
        deformation_determinant: f64,
    },
    /// A structured-box dimension was zero.
    EmptyBox {
        /// Which extent was zero (`0` = x, `1` = y, `2` = z).
        axis: usize,
    },
}

impl std::fmt::Display for MeshError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MeshError::NodeIndexOutOfRange {
                element,
                local_node,
                node,
                node_count,
            } => write!(
                f,
                "element {element} references node {node} at local position \
                 {local_node}, but the mesh has only {node_count} nodes"
            ),
            MeshError::DofCountMismatch { expected, found } => {
                write!(f, "expected {expected} degrees of freedom, got {found}")
            }
            MeshError::DegenerateElement {
                element,
                jacobian_determinant,
            } => write!(
                f,
                "element {element} is degenerate or inverted \
                 (Jacobian determinant {jacobian_determinant})"
            ),
            MeshError::InvertedDeformation {
                element,
                deformation_determinant,
            } => write!(
                f,
                "element {element} has inverted: det F = {deformation_determinant} at an integration point, so its deformation gradient is not a valid state"
            ),
            MeshError::EmptyBox { axis } => {
                write!(f, "box extent along axis {axis} must be positive")
            }
        }
    }
}

impl std::error::Error for MeshError {}

/// A 3-D trilinear hexahedral mesh.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Hex8Mesh {
    nodes: Vec<Vec3>,
    elements: Vec<[usize; 8]>,
}

impl Hex8Mesh {
    /// An empty mesh.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            elements: Vec::new(),
        }
    }

    /// Builds a mesh from raw parts, validating every connectivity entry.
    pub fn from_parts(nodes: Vec<Vec3>, elements: Vec<[usize; 8]>) -> Result<Self, MeshError> {
        let mesh = Self { nodes, elements };
        mesh.validate()?;
        Ok(mesh)
    }

    /// Node coordinates in the reference configuration.
    pub fn nodes(&self) -> &[Vec3] {
        &self.nodes
    }

    /// Element connectivity, in `tpt_fem_element::Hex8` reference-node order.
    pub fn elements(&self) -> &[[usize; 8]] {
        &self.elements
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of elements.
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }

    /// Total degrees of freedom (three translational components per node).
    pub fn dof_count(&self) -> usize {
        3 * self.nodes.len()
    }

    /// Global DOF index of component `component` of `node`.
    pub fn dof(&self, node: usize, component: usize) -> usize {
        3 * node + component
    }

    /// Validates that every connectivity entry is in range.
    pub fn validate(&self) -> Result<(), MeshError> {
        for (e, element) in self.elements.iter().enumerate() {
            for (local_node, &node) in element.iter().enumerate() {
                if node >= self.nodes.len() {
                    return Err(MeshError::NodeIndexOutOfRange {
                        element: e,
                        local_node,
                        node,
                        node_count: self.nodes.len(),
                    });
                }
            }
        }
        Ok(())
    }

    /// Checked access to a displacement vector of the right length.
    fn check_dofs(&self, u: &[f64]) -> Result<(), MeshError> {
        if u.len() != self.dof_count() {
            return Err(MeshError::DofCountMismatch {
                expected: self.dof_count(),
                found: u.len(),
            });
        }
        Ok(())
    }

    /// Nodal positions in the current configuration, `X + u`.
    pub fn positions(&self, u: &[f64]) -> Result<Vec<Vec3>, MeshError> {
        self.check_dofs(u)?;
        Ok((0..self.node_count())
            .map(|n| {
                let d = Vec3::new(u[3 * n], u[3 * n + 1], u[3 * n + 2]);
                self.nodes[n] + d
            })
            .collect())
    }

    /// Element Jacobian `J = dX/dxi` at reference coordinate `xi`.
    ///
    /// `J[k][l] = sum_I X_I[k] dN_I/dxi_l`.
    pub fn jacobian(&self, element: usize, xi: &[f64; 3]) -> Mat3 {
        let grads = Hex8::grad(xi);
        let mut j = Mat3::ZERO;
        for local in 0..8 {
            let x = self.nodes[self.elements[element][local]].to_array();
            for k in 0..3 {
                for l in 0..3 {
                    let v = j.at(k, l) + x[k] * grads[local][l];
                    j.set(k, l, v);
                }
            }
        }
        j
    }

    /// Physical gradients `dN_I/dx_k` at reference coordinate `xi`, one row per
    /// element-local node.
    ///
    /// Returns `None` when the element's Jacobian determinant is non-positive
    /// (inverted or degenerate), which leaves `J^-T` undefined or singular.
    pub fn physical_gradients(&self, element: usize, xi: &[f64; 3]) -> Option<[[f64; 3]; 8]> {
        let j = self.jacobian(element, xi);
        if j.det() <= 0.0 {
            return None;
        }
        let j_inv = j.inverse()?;
        // dN/dx = J^-T dN/dxi  =>  dN/dx_k = sum_l (J^-1)_{l k} dN/dxi_l
        let grads = Hex8::grad(xi);
        let mut out = [[0.0f64; 3]; 8];
        for local in 0..8 {
            for k in 0..3 {
                let mut acc = 0.0;
                for l in 0..3 {
                    acc += j_inv.at(l, k) * grads[local][l];
                }
                out[local][k] = acc;
            }
        }
        Some(out)
    }

    /// Element volume, integrated with an order-2 (`2x2x2`) rule.
    ///
    /// # Errors
    ///
    /// [`MeshError::DegenerateElement`] if any integration point has a
    /// non-positive Jacobian determinant.
    pub fn element_volume(&self, element: usize) -> Result<f64, MeshError> {
        let rule = hex_rule(2);
        let mut volume = 0.0;
        for (xi, w) in rule.points.iter().zip(&rule.weights) {
            let det = self.jacobian(element, xi).det();
            if det <= 0.0 {
                return Err(MeshError::DegenerateElement {
                    element,
                    jacobian_determinant: det,
                });
            }
            volume += w * det;
        }
        Ok(volume)
    }

    /// Indices of elements whose current configuration is inverted or
    /// degenerate, evaluated on an order-2 rule.
    pub fn inverted_elements(&self, u: &[f64]) -> Result<Vec<usize>, MeshError> {
        self.check_dofs(u)?;
        let x = self.positions(u)?;
        let mut out = Vec::new();
        for e in 0..self.element_count() {
            if element_inverted(self.elements[e], &x) {
                out.push(e);
            }
        }
        Ok(out)
    }

    /// Node indices lying on the face perpendicular to `axis` at the minimum
    /// (`at_max == false`) or maximum (`at_max == true`) extent.
    ///
    /// Coordinates are compared with an absolute tolerance of `1e-9`.
    pub fn face_nodes(&self, axis: usize, at_max: bool) -> Vec<usize> {
        let mut target = if at_max {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
        for n in &self.nodes {
            let v = n.to_array()[axis];
            target = if at_max { target.max(v) } else { target.min(v) };
        }
        (0..self.node_count())
            .filter(|&n| (self.nodes[n].to_array()[axis] - target).abs() <= 1e-9)
            .collect()
    }
}

/// True when any point of the element has a non-positive Jacobian
/// determinant in the supplied nodal positions.
fn element_inverted(element: [usize; 8], x: &[Vec3]) -> bool {
    hex_rule(2).points.iter().any(|xi| {
        let grads = Hex8::grad(xi);
        let mut j = Mat3::ZERO;
        for local in 0..8 {
            let p = x[element[local]].to_array();
            for k in 0..3 {
                for l in 0..3 {
                    let v = j.at(k, l) + p[k] * grads[local][l];
                    j.set(k, l, v);
                }
            }
        }
        j.det() <= 0.0
    })
}

/// A structured `nx * ny * nz` hexahedral box spanning
/// `[0, lx] x [0, ly] x [0, lz]`.
///
/// Nodes are laid out with `x` slowest and `z` fastest. Each element's
/// connectivity follows `tpt_fem_element::Hex8` reference-node ordering, so all
/// elements have a strictly positive Jacobian determinant.
///
/// # Errors
///
/// [`MeshError::EmptyBox`] if any count or extent is zero or negative.
pub fn hex_box(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<Hex8Mesh, MeshError> {
    for (axis, n, l) in [(0usize, nx, lx), (1, ny, ly), (2, nz, lz)] {
        if n == 0 || l <= 0.0 {
            return Err(MeshError::EmptyBox { axis });
        }
    }
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;
    let dz = lz / nz as f64;
    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
    for i in 0..=nx {
        for j in 0..=ny {
            for k in 0..=nz {
                nodes.push(Vec3::new(i as f64 * dx, j as f64 * dy, k as f64 * dz));
            }
        }
    }
    let id = |i: usize, j: usize, k: usize| (i * (ny + 1) + j) * (nz + 1) + k;
    let mut elements = Vec::with_capacity(nx * ny * nz);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                elements.push([
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i + 1, j + 1, k),
                    id(i, j + 1, k),
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i + 1, j + 1, k + 1),
                    id(i, j + 1, k + 1),
                ]);
            }
        }
    }
    Hex8Mesh::from_parts(nodes, elements)
}
