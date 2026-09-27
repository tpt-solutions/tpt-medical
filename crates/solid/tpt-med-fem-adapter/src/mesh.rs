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
use tpt_fem_element::{hex_rule, tet_rule, Hex8, ReferenceElement};
use tpt_fem_quadrature::{Quad3D, TetrahedronRule};
use tpt_med_geometry::{Mat3, Vec3};

/// The quadrature and subdivision behaviour that depends on an element's
/// *reference domain*, not on its polynomial order.
///
/// `ReferenceElement` alone cannot pick a rule. A hexahedron lives on
/// `[-1, 1]^3` and is integrated by a tensor-product Gauss rule; a tetrahedron
/// lives on the simplex `(0,0,0)..(1,0,0)` and needs a simplex rule, and the two
/// are not interchangeable. Without this trait the mesh would have to assume
/// "hexahedral", and a tet would silently integrate against the wrong reference
/// domain — a plausible-looking but wrong answer rather than an error.
///
/// Implemented for the element types themselves. There is deliberately no blanket
/// impl: choosing a rule is exactly the decision that has to be explicit.
pub trait ElementFamily {
    /// A quadrature rule exact for this element's own geometry at `order`.
    ///
    /// `order` is interpreted per family: Gauss points per axis for a
    /// hexahedron, and the rule *degree* for a tetrahedron, since a simplex rule
    /// has no tensor-product form to factor through.
    fn quadrature_rule(order: usize) -> Quad3D;

    /// The order [`Mesh::default_quadrature_order`] should use.
    fn natural_quadrature_order() -> usize;

    /// Global node-grid subdivisions per element: `1` for a linear element, `2`
    /// for a quadratic one (corner plus midpoint).
    fn grid_subdivisions() -> usize;

    /// Whether this element lives on the reference *simplex* (a tetrahedron)
    /// rather than the hexahedral cube.
    ///
    /// [`tet_box_of`] needs this to refuse a hexahedral element up front: it
    /// subdivides cells into tetrahedra, and the corner sets it produces are
    /// simply not a hexahedron's eight corners.
    const IS_TETRAHEDRON: bool;
}

/// Hexahedral elements: `[-1, 1]^3`, tensor-product Gauss.
pub trait HexFamily: ReferenceElement {}

/// Tetrahedral elements: the reference simplex, Keast rules.
pub trait TetFamily: ReferenceElement {}

impl HexFamily for Hex8 {}
impl HexFamily for tpt_fem_element::Hex20 {}
impl HexFamily for tpt_fem_element::Hex27 {}
impl TetFamily for tpt_fem_element::Tet4 {}
impl TetFamily for tpt_fem_element::Tet10 {}

/// The hexahedral rule, shared by every `HexFamily` element.
fn hex_element_rule(order: usize) -> Quad3D {
    hex_rule(order)
}

/// The tetrahedral rule, shared by every `TetFamily` element.
///
/// # Why Keast4 is not used
///
/// `tpt_fem_quadrature` 0.1.0's `TetrahedronRule::Keast4` is **wrong**, and using
/// it is not merely inaccurate — it is silently so. Its weights sum to the
/// reference volume 1/6 as they should, and every coordinate is positive, but
/// several of its eleven points have barycentric coordinates summing to more
/// than 1 (measured maximum 1.2607), which places them *outside* the reference
/// tetrahedron. Evaluating a shape function outside its own element is
/// undefined, and the resulting `J` and `det J` are garbage: a `Tet10` uniaxial
/// solve run with that rule came out **42x** the closed-form answer rather than
/// erroring.
///
/// Keast3 is the highest rule that is actually correct here, and it is
/// sufficient: for a `Tet10` the isoparametric map is quadratic, so `J` is
/// linear and `det J` is cubic, which a degree-3 rule integrates exactly.
/// Verified in `a_tet_element_uses_a_simplex_rule_not_a_tensor_product_one`.
///
/// This is worth reporting upstream, and worth re-checking whenever the
/// substrate version moves: if Keast4 is ever fixed, the natural order for a
/// quadratic tet can go back to 4.
fn tet_element_rule(order: usize) -> Quad3D {
    let rule = match order {
        0 | 1 => TetrahedronRule::Degree1,
        2 => TetrahedronRule::Degree2,
        // Clamped at Keast3: Keast4 is unusable, see above.
        _ => TetrahedronRule::Keast3,
    };
    tet_rule(rule)
}

// Explicit impls per element rather than a blanket one over `HexFamily` /
// `TetFamily`: a blanket impl for both would make the two overlap the moment a
// downstream type claimed both markers, and Rust's coherence rules would reject
// the crate outright. Spelling out five impls is not the cost.
macro_rules! hex_family {
    ($($t:ty),+ $(,)?) => {
        $(impl ElementFamily for $t {
            fn quadrature_rule(order: usize) -> Quad3D { hex_element_rule(order) }
            // Order 3 for a quadratic element, not 2: its Jacobian is not
            // integrated exactly at order 2, and under-integrating a curved
            // element's own geometry is a silent accuracy loss.
            fn natural_quadrature_order() -> usize { if Self::NUM_NODES > 8 { 3 } else { 2 } }
            fn grid_subdivisions() -> usize { if Self::NUM_NODES > 8 { 2 } else { 1 } }
            const IS_TETRAHEDRON: bool = false;
        })+
    };
}

macro_rules! tet_family {
    ($($t:ty),+ $(,)?) => {
        $(impl ElementFamily for $t {
            fn quadrature_rule(order: usize) -> Quad3D { tet_element_rule(order) }
            // Degree 3, i.e. Keast3: a `Tet10` map is quadratic, so `det J` is
            // cubic and degree 3 integrates it exactly. Degree 4 would be the
            // textbook choice but the substrate's Keast4 rule is broken — see
            // `tet_element_rule`.
            fn natural_quadrature_order() -> usize { if Self::NUM_NODES > 8 { 3 } else { 2 } }
            fn grid_subdivisions() -> usize { if Self::NUM_NODES > 8 { 2 } else { 1 } }
            const IS_TETRAHEDRON: bool = true;
        })+
    };
}

hex_family!(Hex8, tpt_fem_element::Hex20, tpt_fem_element::Hex27);
tet_family!(tpt_fem_element::Tet4, tpt_fem_element::Tet10);

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
/// `Hex8Mesh` remains a type alias for `Mesh<Hex8>`, so every existing
/// signature, caller and test is unchanged. New code should name
/// `Mesh<Hex20>` explicitly rather than adding a parallel `Hex20Mesh` type:
/// a second type would be a second implementation to keep in step.
///
/// # Quadrature order
///
/// Quadratic elements get a higher default order than `Hex8`. A trilinear
/// element integrates its own Jacobian exactly at order 2, but a quadratic
/// one does not, and under-integrating a curved element flattens it. The
/// element type therefore carries its own floor via
/// [`Mesh::default_quadrature_order`], which callers should prefer over a
/// hard-coded order.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Mesh<E: ReferenceElement + ElementFamily> {
    nodes: Vec<Vec3>,
    elements: Vec<Vec<usize>>,
    _element: PhantomData<E>,
}

/// The trilinear hexahedral mesh this crate was originally scoped to.
///
/// An alias, not its own type: see [`Mesh`].
pub type Hex8Mesh = Mesh<Hex8>;

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
    /// [`tet_box_of`] was asked for an element that is not a tetrahedron.
    UnsupportedTetBox {
        /// The node count of the element that was supplied, for diagnosis.
        nodes_per_element: usize,
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
            MeshError::UnsupportedTetBox { nodes_per_element } => write!(
                f,
                "tet_box_of needs a tetrahedral element, but this one has \
                 {nodes_per_element} nodes on a hexahedral domain"
            ),
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

impl<E: ReferenceElement + ElementFamily> Mesh<E> {
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
        E::natural_quadrature_order()
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
        let rule = E::quadrature_rule(Self::default_quadrature_order());
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
fn element_inverted<E: ReferenceElement + ElementFamily>(element: &[usize], x: &[Vec3]) -> bool {
    E::quadrature_rule(Mesh::<E>::default_quadrature_order())
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
pub fn hex_box_of<E: ReferenceElement + ElementFamily>(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<Mesh<E>, MeshError> {
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
    Mesh::<E>::from_parts(nodes, elements)
}

/// A structured `nx * ny * nz` trilinear hexahedral box.
///
/// The six tetrahedra of one hexahedral cell, in `Hex8` corner order.
///
/// The Freudenthal (Kuhn) subdivision: all six share the cube's main diagonal,
/// which in `Hex8` numbering is corner **0 to corner 6** — `(0,0,0)` to
/// `(1,1,1)`. The other six corners form the equator in the order `1, 2, 3, 7,
/// 4, 5`, walking once around the diagonal.
///
/// Getting that apex wrong is a silent failure, and worth stating because it is
/// easy to do: numbering the diagonal as `0..7` instead produces two
/// tetrahedra whose four vertices are coplanar — one lying entirely in an
/// `x = 0` face — so the mesh looks structurally fine while two elements have
/// zero volume. The orientation fix-up in [`tet_box_of`] cannot rescue a flat
/// tetrahedron.
const KUHN_TETS: [[usize; 4]; 6] = [
    [0, 1, 2, 6],
    [0, 2, 3, 6],
    [0, 3, 7, 6],
    [0, 7, 4, 6],
    [0, 4, 5, 6],
    [0, 5, 1, 6],
];

/// `Tet10`'s six mid-edge nodes as `(vertex_a, vertex_b)` pairs, in its own
/// reference ordering: mid 1-2, 2-3, 1-3, 1-4, 2-4, 3-4, taken from
/// `tpt_fem_element::Tet10::nodes`.
const TET_EDGES: [(usize, usize); 6] = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)];

/// A structured `nx * ny * nz` tetrahedral box spanning
/// `[0, lx] x [0, ly] x [0, lz]`, built from element type `E`.
///
/// Each cell becomes six tetrahedra sharing a main diagonal (see
/// [`KUHN_TETS`]), so `nx * ny * nz` keeps meaning *cells*, as in
/// [`hex_box_of`], and the element count is six times that.
///
/// # Quadratic elements
///
/// The mid-edge nodes are created per *tetrahedral edge* and shared by every
/// tetrahedron using it, keyed on the sorted corner pair. That includes edges
/// on internal faces and on the shared body diagonal, which is what makes the
/// quadratic mesh conforming too: an unshared mid-node on an internal face is a
/// crack in exactly the way a mismatched face split is.
///
/// The consequence worth stating plainly: a `Tet10` box has interior nodes on
/// the face and body diagonals, not only on the grid lines. That is inherent to
/// subdividing a hexahedron into tetrahedra rather than a defect here, but it
/// does mean a `Tet10` box is not a drop-in replacement for a `Hex20` box of the
/// same dimensions — the node counts differ substantially.
///
/// # Errors
///
/// [`MeshError::EmptyBox`] if any count or extent is zero or negative, and
/// [`MeshError::UnsupportedTetBox`] if `E` is not a tetrahedral element.
pub fn tet_box_of<E: ReferenceElement + ElementFamily>(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<Mesh<E>, MeshError> {
    if !E::IS_TETRAHEDRON {
        return Err(MeshError::UnsupportedTetBox {
            nodes_per_element: E::NUM_NODES,
        });
    }
    for (axis, n, l) in [(0usize, nx, lx), (1, ny, ly), (2, nz, lz)] {
        if n == 0 || l <= 0.0 {
            return Err(MeshError::EmptyBox { axis });
        }
    }
    let (dx, dy, dz) = (lx / nx as f64, ly / ny as f64, lz / nz as f64);
    let mut nodes: Vec<Vec3> = Vec::new();
    let mut elements: Vec<Vec<usize>> = Vec::with_capacity(6 * nx * ny * nz);
    // Corners keyed on grid position, so a node shared by several cells is made
    // once. Without this the mesh has duplicate coincident nodes and no
    // stiffness ever passes between cells.
    let mut corner_of: std::collections::HashMap<(usize, usize, usize), usize> =
        std::collections::HashMap::new();
    // Mid-edge nodes keyed on the sorted corner pair, so every tetrahedron on an
    // edge shares the one node.
    let mut mid_of: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();

    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let grid = [
                    (i, j, k),
                    (i + 1, j, k),
                    (i + 1, j + 1, k),
                    (i, j + 1, k),
                    (i, j, k + 1),
                    (i + 1, j, k + 1),
                    (i + 1, j + 1, k + 1),
                    (i, j + 1, k + 1),
                ];
                let mut cell = [0usize; 8];
                for (slot, key) in cell.iter_mut().zip(grid) {
                    *slot = match corner_of.get(&key) {
                        Some(&existing) => existing,
                        None => {
                            let index = nodes.len();
                            nodes.push(Vec3::new(
                                key.0 as f64 * dx,
                                key.1 as f64 * dy,
                                key.2 as f64 * dz,
                            ));
                            corner_of.insert(key, index);
                            index
                        }
                    };
                }
                for tet in KUHN_TETS {
                    let mut verts = [cell[tet[0]], cell[tet[1]], cell[tet[2]], cell[tet[3]]];
                    // Orient positively. A negative tetrahedron would later be
                    // reported as "degenerate", blaming the physics for a meshing
                    // choice; fixing it here keeps that error honest. The sign
                    // convention matches the reference tetrahedron
                    // `(0,0,0),(1,0,0),(0,1,0),(0,0,1)`, whose volume is +1/6.
                    if signed_volume(&nodes, &verts) < 0.0 {
                        verts.swap(1, 2);
                    }
                    let mut element: Vec<usize> = verts.to_vec();
                    if E::NUM_NODES > 4 {
                        for (a, b) in TET_EDGES {
                            let (p, q) = if verts[a] < verts[b] {
                                (verts[a], verts[b])
                            } else {
                                (verts[b], verts[a])
                            };
                            let index = match mid_of.get(&(p, q)) {
                                Some(&existing) => existing,
                                None => {
                                    let mid = (nodes[p] + nodes[q]) * 0.5;
                                    let existing = nodes.len();
                                    nodes.push(mid);
                                    mid_of.insert((p, q), existing);
                                    existing
                                }
                            };
                            element.push(index);
                        }
                    }
                    elements.push(element);
                }
            }
        }
    }
    Mesh::<E>::from_parts(nodes, elements)
}

/// Signed volume of a tetrahedron: one sixth of the scalar triple product.
fn signed_volume(nodes: &[Vec3], verts: &[usize; 4]) -> f64 {
    let p: Vec<[f64; 3]> = verts.iter().map(|&v| nodes[v].to_array()).collect();
    let (a, b, c) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
    let cross = [
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0],
    ];
    dot(a, cross) / 6.0
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// A structured `nx * ny * nz` tetrahedral box of `Tet4`.
pub fn tet_box(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<Mesh<tpt_fem_element::Tet4>, MeshError> {
    tet_box_of::<tpt_fem_element::Tet4>(nx, ny, nz, lx, ly, lz)
}

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
