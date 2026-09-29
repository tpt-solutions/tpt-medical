//! Global FEM model: assembly, boundary conditions, solve, post-processing.

use std::collections::BTreeSet;

use crate::hex::{centre_strain_displacement, check_element, trilinear_hex_stiffness};
use crate::results::{ElementStress, StressResult};
use crate::sparse::{conjugate_gradient, CsrMatrix};
use tpt_med_geometry::{Mat3, Vec3};

/// A linear-elastic biomechanics model over a hex mesh.
#[derive(Debug, Clone)]
pub struct BiomechanicsModel {
    /// Node positions (mm).
    pub nodes: Vec<Vec3>,
    /// Hex connectivity (node indices per element, 8 each).
    pub elements: Vec<[u32; 8]>,
    /// Young's modulus per element (MPa).
    pub element_modulus: Vec<f64>,
    /// Poisson's ratio per element.
    pub element_poisson: Vec<f64>,
}

/// Solver configuration error.
#[derive(Debug)]
pub enum SolverError {
    /// Model inconsistency (empty mesh, mismatched arrays, degenerate
    /// elements, singular system).
    Invalid(String),
}

impl core::fmt::Display for SolverError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SolverError::Invalid(why) => write!(f, "invalid model: {why}"),
        }
    }
}

impl std::error::Error for SolverError {}

/// Boundary conditions: fixed DOFs plus nodal forces (N).
#[derive(Debug, Clone, Default)]
pub struct BoundaryConditions {
    /// Fully constrained nodes (all 3 DOFs).
    pub fixed_nodes: BTreeSet<u32>,
    /// Nodal forces `(node, force)` in newtons. Summed on duplicates.
    pub nodal_forces: Vec<(u32, Vec3)>,
    /// Optionally prescribed nodal displacements (mm) — applied after
    /// `fixed_nodes`.
    pub prescribed_displacements: Vec<(u32, Vec3)>,
    /// Per-node partial constraints — one axis flag per DOF
    /// (`[x, y, z]`, `true` = constrained): symmetry planes and roller
    /// supports constrain individual axes, unlike all-3-DOFs
    /// [`Self::fixed_nodes`]. Folded into the constraint elimination after
    /// `fixed_nodes` and before prescribed displacements.
    pub constrained_dofs: Vec<(u32, [bool; 3])>,
}

/// Per-axis DOF constraint mask helpers.
impl BoundaryConditions {
    /// Constrains selected axes of a node: `[true, false, false]` pins x
    /// only (a symmetry plane normal to x); `[true, true, false]` is a
    /// roller in the x-y plane.
    pub fn constrain_dofs(
        &mut self,
        nodes: impl IntoIterator<Item = u32>,
        axes: [bool; 3],
    ) -> &mut Self {
        for n in nodes {
            self.constrained_dofs.push((n, axes));
        }
        self
    }
}

impl BoundaryConditions {
    /// Fixes all nodes of an element set (convenience for faces selected by
    /// coordinate tests).
    pub fn fix_nodes(&mut self, nodes: impl IntoIterator<Item = u32>) {
        self.fixed_nodes.extend(nodes);
    }

    /// Adds a nodal force (N).
    pub fn add_force(&mut self, node: u32, force: Vec3) {
        self.nodal_forces.push((node, force));
    }
}

impl BiomechanicsModel {
    /// Adapts a voxel mesh (node positions + per-element HU-derived
    /// moduli).
    pub fn from_voxel_mesh(mesh: &tpt_med_meshing::VoxelHexMesh) -> Self {
        Self {
            nodes: mesh.nodes.clone(),
            elements: mesh.elements.clone(),
            element_modulus: mesh.materials.iter().map(|m| m.youngs_modulus).collect(),
            element_poisson: mesh.materials.iter().map(|m| m.poissons_ratio).collect(),
        }
    }

    /// Builds a model from raw parts (e.g. synthetic beam meshes).
    pub fn from_parts(
        nodes: Vec<Vec3>,
        elements: Vec<[u32; 8]>,
        modulus: f64,
        poisson: f64,
    ) -> Self {
        let n_el = elements.len();
        Self {
            nodes,
            elements,
            element_modulus: vec![modulus; n_el],
            element_poisson: vec![poisson; n_el],
        }
    }

    /// Assembles and solves the static system, then post-processes element
    /// stresses.
    pub fn solve(
        &self,
        bc: &BoundaryConditions,
        tolerance: f64,
        max_iterations: usize,
    ) -> Result<StressResult, SolverError> {
        if self.nodes.is_empty() || self.elements.is_empty() {
            return Err(SolverError::Invalid("empty mesh".into()));
        }
        if self.element_modulus.len() != self.elements.len()
            || self.element_poisson.len() != self.elements.len()
        {
            return Err(SolverError::Invalid(
                "material array length mismatch".into(),
            ));
        }

        let n_nodes = self.nodes.len();
        let n_dofs = 3 * n_nodes;

        // Geometry sanity on every element.
        let mut geometry_errors = Vec::new();
        for (id, el) in self.elements.iter().enumerate() {
            let mut corner_nodes = [Vec3::ZERO; 8];
            for (a, &n) in el.iter().enumerate() {
                match self.nodes.get(n as usize) {
                    Some(&p) => corner_nodes[a] = p,
                    None => {
                        return Err(SolverError::Invalid(format!(
                            "element {id} references missing node {n}"
                        )))
                    }
                }
            }
            check_element(&corner_nodes, id, &mut geometry_errors);
        }
        if !geometry_errors.is_empty() {
            return Err(SolverError::Invalid(geometry_errors.join("; ")));
        }

        // ---- assembly -------------------------------------------------
        let mut triplets = Vec::new();
        let mut force = vec![0.0f64; n_dofs];
        for &(node, f) in &bc.nodal_forces {
            let i = node as usize;
            if i >= n_nodes {
                return Err(SolverError::Invalid(format!(
                    "force on missing node {node}"
                )));
            }
            force[3 * i] += f.x;
            force[3 * i + 1] += f.y;
            force[3 * i + 2] += f.z;
        }

        for (id, el) in self.elements.iter().enumerate() {
            let mut corner_nodes = [Vec3::ZERO; 8];
            for (a, &n) in el.iter().enumerate() {
                corner_nodes[a] = self.nodes[n as usize];
            }
            let ke = trilinear_hex_stiffness(
                &corner_nodes,
                self.element_modulus[id],
                self.element_poisson[id],
            );
            for (a, &na) in el.iter().enumerate() {
                for (b, &nb) in el.iter().enumerate() {
                    for i in 0..3 {
                        for j in 0..3 {
                            let v = ke[3 * a + i][3 * b + j];
                            if v != 0.0 {
                                triplets.push((3 * na as usize + i, 3 * nb as usize + j, v));
                            }
                        }
                    }
                }
            }
        }
        let n_triplets = triplets.len();
        let k = CsrMatrix::from_triplets(n_dofs, triplets);

        // ---- constrain: eliminate fixed DOFs by zeroing rows/cols and
        // putting 1 on the diagonal -----------------------------------
        let mut fixed = vec![false; n_dofs];
        for &n in &bc.fixed_nodes {
            if n as usize >= n_nodes {
                return Err(SolverError::Invalid(format!("fix on missing node {n}")));
            }
            for i in 0..3 {
                fixed[3 * n as usize + i] = true;
            }
        }
        for &(n, axes) in &bc.constrained_dofs {
            if n as usize >= n_nodes {
                return Err(SolverError::Invalid(format!(
                    "partial constraint on missing node {n}"
                )));
            }
            for (i, &on) in axes.iter().enumerate() {
                if on {
                    fixed[3 * n as usize + i] = true;
                }
            }
        }
        for &(n, _d) in &bc.prescribed_displacements {
            if n as usize >= n_nodes {
                return Err(SolverError::Invalid(format!(
                    "displacement on missing node {n}"
                )));
            }
            for i in 0..3 {
                fixed[3 * n as usize + i] = true;
            }
        }

        // K_c = K with fixed rows/columns zeroed and unit diagonal.
        let mut constrained_triplets = Vec::with_capacity(n_triplets);
        for i in 0..n_dofs {
            for kk in k.indptr[i]..k.indptr[i + 1] {
                let j = k.indices[kk];
                if !fixed[i] && !fixed[j] {
                    constrained_triplets.push((i, j, k.data[kk]));
                }
            }
        }
        for i in 0..n_dofs {
            if fixed[i] {
                constrained_triplets.push((i, i, 1.0));
            }
        }
        let kc = CsrMatrix::from_triplets(n_dofs, constrained_triplets);

        // Prescribed displacement load terms: f_mod = f − K_up · d_prescribed.
        for &(n, d) in &bc.prescribed_displacements {
            let i = n as usize;
            for c in 0..3 {
                if d.to_array()[c] == 0.0 {
                    continue;
                }
                let dof = 3 * i + c;
                // Coupling into free DOFs via row `dof` of the raw K.
                for kk in k.indptr[dof]..k.indptr[dof + 1] {
                    let jj = k.indices[kk];
                    if !fixed[jj] {
                        force[jj] -= k.data[kk] * d.to_array()[c];
                    }
                }
            }
        }
        for &n in &bc.fixed_nodes {
            for i in 0..3 {
                force[3 * n as usize + i] = 0.0;
            }
        }

        // ---- solve ----------------------------------------------------
        let (u, stats) = conjugate_gradient(&kc, &force, tolerance, max_iterations);
        if !stats.converged {
            return Err(SolverError::Invalid(format!(
                "CG did not converge: {} iterations, relative residual {:.3e}",
                stats.iterations, stats.relative_residual
            )));
        }

        // Prescribed displacements back into the solution vector.
        let mut displacements = Vec::with_capacity(n_nodes);
        for i in 0..n_nodes {
            displacements.push(Vec3::new(u[3 * i], u[3 * i + 1], u[3 * i + 2]));
        }
        for &(n, d) in &bc.prescribed_displacements {
            displacements[n as usize] = d;
        }

        // ---- post-process ----------------------------------------------
        let mut stresses = Vec::with_capacity(self.elements.len());
        for (id, el) in self.elements.iter().enumerate() {
            let mut corner_nodes = [Vec3::ZERO; 8];
            let mut ue = [0.0f64; 24];
            for (a, &n) in el.iter().enumerate() {
                corner_nodes[a] = self.nodes[n as usize];
                let d = displacements[n as usize];
                ue[3 * a] = d.x;
                ue[3 * a + 1] = d.y;
                ue[3 * a + 2] = d.z;
            }
            let (b_mat, _det) = centre_strain_displacement(&corner_nodes);
            let mut strain = [0.0f64; 6];
            for r in 0..6 {
                for c in 0..24 {
                    strain[r] += b_mat[r][c] * ue[c];
                }
            }
            let d_mat = crate::hex::isotropic_d(self.element_modulus[id], self.element_poisson[id]);
            let mut stress = [0.0f64; 6];
            for r in 0..6 {
                for c in 0..6 {
                    stress[r] += d_mat[r][c] * strain[c];
                }
            }
            stresses.push(ElementStress {
                element: id as u32,
                strain,
                stress,
            });
        }

        Ok(StressResult {
            displacements,
            stats,
            stresses,
        })
    }

    /// Elements whose centroid satisfies `keep(centroid)` — helper for
    /// selecting faces by coordinates.
    pub fn element_centroid(&self, id: usize) -> Vec3 {
        let el = &self.elements[id];
        let mut c = Vec3::ZERO;
        for &n in el.iter() {
            c += self.nodes[n as usize];
        }
        c / 8.0
    }

    /// Global stiffness matrix invariants used by tests/verification:
    /// rigid-body translations must produce zero force.
    pub fn rigid_mode_residual(&self) -> f64 {
        let mut worst = 0.0f64;
        for dof in 0..3 {
            let mut triplets = Vec::new();
            for (id, el) in self.elements.iter().enumerate() {
                let mut corner_nodes = [Vec3::ZERO; 8];
                for (a, &n) in el.iter().enumerate() {
                    corner_nodes[a] = self.nodes[n as usize];
                }
                let ke = trilinear_hex_stiffness(
                    &corner_nodes,
                    self.element_modulus[id],
                    self.element_poisson[id],
                );
                for (a, &na) in el.iter().enumerate() {
                    for (b, &nb) in el.iter().enumerate() {
                        triplets.push((
                            3 * na as usize + dof,
                            3 * nb as usize + dof,
                            ke[3 * a + dof][3 * b + dof],
                        ));
                    }
                }
            }
            let k = CsrMatrix::from_triplets(3 * self.nodes.len(), triplets);
            let ones = vec![1.0; 3 * self.nodes.len()];
            let mut y = vec![0.0; 3 * self.nodes.len()];
            k.mul_vec(&ones, &mut y);
            for v in y {
                worst = worst.max(v.abs());
            }
        }
        worst
    }

    /// Convenience access to the stiffness coupling used in verification.
    pub fn coupling(&self) -> Mat3 {
        Mat3::IDENTITY
    }
}
