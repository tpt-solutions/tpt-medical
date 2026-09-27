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
//!     vec![vec![0, 1, 2, 3, 4, 5, 6, 99]],
//! );
//! assert!(matches!(bad, Err(MeshError::NodeIndexOutOfRange { .. })));
//! // An element with the wrong node count is rejected too, rather than being
//! // read past the end of its own connectivity.
//! let short = Hex8Mesh::from_parts(vec![Vec3::ZERO; 8], vec![vec![0, 1, 2]]);
//! assert!(matches!(
//!     short,
//!     Err(MeshError::WrongElementNodeCount { .. })
//! ));
//! ```

use std::marker::PhantomData;
use tpt_fem_element::{hex_rule, Hex8, ReferenceElement};
use tpt_med_geometry::{Mat3, Vec3};

/// A 3-D hexahedral mesh whose element shape is chosen by the reference
/// element `E`.
///
/// # Why this is generic
///
/// The adapter's job is the *assembly* — the internal force, the tangent, the
/// contact coupling — and none of that is specific to a node count. Every
/// quantity here is already written in terms of `ReferenceElement::grad` and a
/// quadrature rule, so parameterising on the element costs nothing at the
/// formulation level and buys `Hex20`/`Hex27` for curved geometry without a
/// second assembly.
///
/// `Hex8Mesh` remains a type alias for `HexMesh<Hex8>`, so every existing
/// signature, caller and test is unchanged. New code should name
/// `HexMesh<Hex20>` explicitly rather than adding a parallel `Hex20Mesh` type:
/// a second type would be a second implementation to keep in step.
///
/// # Quadrature order
///
/// Quadratic elements get a higher default order than `Hex8`. A trilinear
/// element integrates its own Jacobian exactly at order 2, but a quadratic
/// one does not, and under-integrating a curved element flattens it. The
/// element type therefore carries its own floor via
/// [`HexMesh::default_quadrature_order`], which callers should prefer over a
/// hard-coded order.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct HexMesh<E: ReferenceElement> {
    nodes: Vec<Vec3>,
    elements: Vec<Vec<usize>>,
    _element: PhantomData<E>,
}

/// The trilinear hexahedral mesh this crate was originally scoped to.
///
/// An alias, not its own type: see [`HexMesh`].
pub type Hex8Mesh = HexMesh<Hex8>;

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
    /// An element's connectivity did not have the node count the reference
    /// element requires.
    WrongElementNodeCount {
        /// The offending element index.
        element: usize,
        /// Nodes the element actually listed.
        found: usize,
        /// Nodes the reference element requires.
        expected: usize,
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
            MeshError::WrongElementNodeCount {
                element,
                found,
                expected,
            } => write!(
                f,
                "element {element} has {found} nodes, but this element type needs {expected}"
            ),
        }
    }
}

impl std::error::Error for MeshError {}

impl<E: ReferenceElement> HexMesh<E> {
    /// Nodes per element, from the reference element.
    ///
    /// A function rather than an associated const: `E::NUM_NODES` cannot be
    /// used as an associated const on stable, and the places that need this
    /// number are all runtime (allocation sizes, loop bounds) anyway.
    pub fn nodes_per_element() -> usize {
        E::NUM_NODES
    }

    /// An empty mesh.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            elements: Vec::new(),
            _element: PhantomData,
        }
    }

    /// A quadrature order that integrates this element's own geometry exactly.
    ///
    /// `3` for a quadratic element, `2` for `Hex8`. Prefer this over a
    /// hard-coded order: under-integrating a curved element's Jacobian is a
    /// silent accuracy loss that shows up as a mesh that behaves more stiffly
    /// than it should, not as an error.
    pub fn default_quadrature_order() -> usize {
        if E::NUM_NODES > 8 {
            3
        } else {
            2
        }
    }

    /// Builds a mesh from raw parts, validating every connectivity entry.
    ///
    /// # Errors
    ///
    /// [`MeshError::WrongElementNodeCount`] if any element does not have
    /// `E::NUM_NODES` nodes, and [`MeshError::NodeIndexOutOfRange`] if any
    /// referenced node is outside the mesh.
    pub fn from_parts(nodes: Vec<Vec3>, elements: Vec<Vec<usize>>) -> Result<Self, MeshError> {
        let mesh = Self {
            nodes,
            elements,
            _element: PhantomData,
        };
        mesh.validate()?;
        Ok(mesh)
    }

    /// Node coordinates in the reference configuration.
    pub fn nodes(&self) -> &[Vec3] {
        &self.nodes
    }

    /// Element connectivity, in `E`'s reference-node order.
    pub fn elements(&self) -> &[Vec<usize>] {
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

    /// Validates every element's node count and every connectivity entry.
    pub fn validate(&self) -> Result<(), MeshError> {
        for (e, element) in self.elements.iter().enumerate() {
            // Checked before the indices: a short element would otherwise index
            // out of range inside the loop below, and report a confusing error.
            if element.len() != E::NUM_NODES {
                return Err(MeshError::WrongElementNodeCount {
                    element: e,
                    found: element.len(),
                    expected: E::NUM_NODES,
                });
            }
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
        let grads = E::grad(xi);
        let mut j = Mat3::ZERO;
        for (local, &node) in self.elements[element].iter().enumerate() {
            let x = self.nodes[node].to_array();
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
    /// element-local node (so `E::NUM_NODES` rows).
    ///
    /// Returns `None` when the element's Jacobian determinant is non-positive
    /// (inverted or degenerate), which leaves `J^-T` undefined or singular.
    pub fn physical_gradients(&self, element: usize, xi: &[f64; 3]) -> Option<Vec<[f64; 3]>> {
        let j = self.jacobian(element, xi);
        if j.det() <= 0.0 {
            return None;
        }
        let j_inv = j.inverse()?;
        // dN/dx = J^-T dN/dxi  =>  dN/dx_k = sum_l (J^-1)_{l k} dN/dxi_l
        let grads = E::grad(xi);
        // A `Vec` rather than a fixed-size array: `E::NUM_NODES` is an
        // associated const and cannot be used as an array length on stable, so
        // the gradient rows are allocated at the size the element declares.
        let mut out = vec![[0.0f64; 3]; E::NUM_NODES];
        for (local, grad) in grads.iter().enumerate().take(E::NUM_NODES) {
            for k in 0..3 {
                let mut acc = 0.0;
                for l in 0..3 {
                    acc += j_inv.at(l, k) * grad[l];
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
        let rule = hex_rule(Self::default_quadrature_order());
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
            if element_inverted::<E>(&self.elements[e], &x) {
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
fn element_inverted<E: ReferenceElement>(element: &[usize], x: &[Vec3]) -> bool {
    hex_rule(HexMesh::<E>::default_quadrature_order())
        .points
        .iter()
        .any(|xi| {
            let grads = E::grad(xi);
            let mut j = Mat3::ZERO;
            for (local, &node) in element.iter().enumerate().take(E::NUM_NODES) {
                let p = x[node].to_array();
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
/// `[0, lx] x [0, ly] x [0, lz]`, built from element type `E`.
///
/// Nodes are laid out with `x` slowest and `z` fastest. Each element's
/// connectivity follows `E`'s reference-node ordering, so all elements have a
/// strictly positive Jacobian determinant.
///
/// For a quadratic element the *global* node grid is refined twice per element
/// along each axis, which is what makes mid-edge nodes exist to connect. The
/// element count is still `nx * ny * nz` — `nx` is a count of elements, not of
/// cells in the node grid, so this does not change the meaning of the
/// arguments between `Hex8` and `Hex20`.
///
/// # Errors
///
/// [`MeshError::EmptyBox`] if any count or extent is zero or negative.
pub fn hex_box_of<E: ReferenceElement>(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<HexMesh<E>, MeshError> {
    for (axis, n, l) in [(0usize, nx, lx), (1, ny, ly), (2, nz, lz)] {
        if n == 0 || l <= 0.0 {
            return Err(MeshError::EmptyBox { axis });
        }
    }
    // Global node-grid subdivisions per element: 1 for a linear element, 2 for
    // a quadratic one (corner + midpoint).
    let s = if E::NUM_NODES > 8 { 2 } else { 1 };
    let (_gx, gy, gz) = (nx * s + 1, ny * s + 1, nz * s + 1);
    // The grid has `gx = nx * s + 1` nodes, so the spacing that makes the last
    // node land exactly on `lx` is `lx / (nx * s)` — dividing by the element
    // count would overshoot to `lx * s` (a Hex20 box of nominal side `l` coming
    // out with side `2l`), and dividing by the grid count would undershoot. The
    // arithmetic is the one place where a plausible-looking box can be
    // uniformly the wrong size, so it is pinned by the volume and extent tests.
    let (dx, dy, dz) = (
        lx / (nx * s) as f64,
        ly / (ny * s) as f64,
        lz / (nz * s) as f64,
    );
    // A serendipity element (`Hex20`) has no face-centre or body-centre node, so
    // on a twice-refined grid those grid positions are referenced by *no*
    // element. Leaving them in would give orphan nodes with an exactly zero
    // stiffness row — free DOFs that make the condensed system singular for a
    // reason that has nothing to do with the physics. So the connectivity is
    // built first and the node list is then compacted to the nodes actually
    // referenced, which is what a mesher would hand over anyway.
    let id = |i: usize, j: usize, k: usize| (i * gy + j) * gz + k;
    let map = |xi: f64, i: usize| i * s + (((xi + 1.0) * 0.5 * s as f64).round() as usize);
    let mut elements = Vec::with_capacity(nx * ny * nz);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let element: Vec<usize> = E::nodes()
                    .iter()
                    .map(|n| id(map(n[0], i), map(n[1], j), map(n[2], k)))
                    .collect();
                elements.push(element);
            }
        }
    }
    // Compact: grid position -> node index, assigned in first-referenced order
    // so the numbering is deterministic.
    let mut compact: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut nodes: Vec<Vec3> = Vec::new();
    for element in &mut elements {
        for slot in element.iter_mut() {
            let grid = *slot;
            let next = compact.len();
            let index = *compact.entry(grid).or_insert_with(|| {
                let (gi, rem) = (grid / (gy * gz), grid % (gy * gz));
                let (gj, gk) = (rem / gz, rem % gz);
                nodes.push(Vec3::new(gi as f64 * dx, gj as f64 * dy, gk as f64 * dz));
                next
            });
            *slot = index;
        }
    }
    HexMesh::<E>::from_parts(nodes, elements)
}

/// A structured `nx * ny * nz` trilinear hexahedral box.
///
/// The `Hex8` specialisation of [`hex_box_of`].
pub fn hex_box(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<Hex8Mesh, MeshError> {
    hex_box_of::<Hex8>(nx, ny, nz, lx, ly, lz)
}
