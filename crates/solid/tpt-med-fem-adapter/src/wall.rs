//! Rigid analytic cylindrical wall: radial penalty contact with history-
//! carrying regularized Coulomb friction (RFC 0013's vessel).
//!
//! [`crate::contact::ContactPairing`] resolves gaps along a single coordinate
//! axis against sampled obstacle points, which cannot express a cylinder. A
//! stent deploys into a lumen, so this module adds the wall as an *analytic*
//! surface: for a slave node at transverse radius `ρ` about the cylinder
//! axis, the signed penetration is `ρ − R` (body inside the wall) or
//! `R − ρ` (body outside it), and the penalty energy is `½ κ pen²` for
//! `pen > 0`. Residual and tangent are the exact gradient and Hessian of that
//! energy, including the curvature term `κ·pen/ρ (I − r̂r̂ᵀ)`, so Newton stays
//! quadratic through the wall.
//!
//! # Friction
//!
//! The planar layer in [`crate::friction`] has no memory: its slip is the
//! offset to a fixed master point. A wall needs real slip, so this layer is
//! an elastic-predictor / Coulomb-corrector with explicit history:
//!
//! ```text
//! Δs   = (u − u_anchor) with its wall-normal part removed
//! trial = f_held,t − k_t Δs          (f_held projected on the tangent plane)
//! f     = trial                      if |trial| ≤ μ N     (stick)
//!       = μ N trial/|trial|          otherwise            (slip)
//! ```
//!
//! `u_anchor` and `f_held` are *committed* state, advanced only by
//! [`RadialWall::committed`] after a converged increment — exactly as the
//! superelastic field is — so a cut-back increment restarts from the same
//! friction state. The tangent includes the stick stiffness, the radial-
//! return term and the `∂(μN)/∂u` coupling; it neglects the derivative of
//! the tangent projector itself (an `O(|Δs|/ρ)` effect), which the Jacobian
//! finite-difference test bounds.
//!
//! **Stated limits:** a node that enters contact mid-increment carries its
//! tangential motion since the increment start as slip (an error of at most
//! one increment — refine the stage); a rigid wall has no compliance.

use crate::friction::FrictionConfig;
use crate::mesh::{Mesh, MeshError};
use tpt_fem_element::ReferenceElement;
use tpt_fem_sparse::Coo;

const TINY: f64 = 1e-300;

/// Which side of the wall the deforming body lives on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WallSide {
    /// The body is inside the cylinder (a stent in a lumen): penetration is
    /// `ρ − R`.
    Inside,
    /// The body is outside the cylinder (a sheath-crimped mandrel): penetration
    /// is `R − ρ`.
    Outside,
}

/// Errors defining a wall.
#[derive(Debug, Clone, PartialEq)]
pub enum WallError {
    /// The axis is not `0`, `1` or `2`.
    AxisOutOfRange(usize),
    /// The radius is not finite and positive.
    InvalidRadius(f64),
    /// The centre is not finite.
    InvalidCenter,
    /// A radius profile is not finite and strictly increasing in `z`.
    InvalidProfile,
}

impl std::fmt::Display for WallError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AxisOutOfRange(a) => write!(f, "wall axis {a} is not 0, 1 or 2"),
            Self::InvalidRadius(r) => write!(f, "wall radius {r} must be finite and positive"),
            Self::InvalidCenter => write!(f, "wall centre must be finite"),
            Self::InvalidProfile => {
                write!(f, "radius profile must be finite with increasing z")
            }
        }
    }
}

impl std::error::Error for WallError {}

/// A rigid analytic cylindrical wall with slave nodes.
#[derive(Debug, Clone, PartialEq)]
pub struct RadialWall {
    axis: usize,
    center: [f64; 2],
    radius: f64,
    side: WallSide,
    slave: Vec<usize>,
    friction: Option<FrictionConfig>,
    /// Axial radius offsets `(z, Δr)`, linearly interpolated and held flat
    /// beyond the ends; empty for a straight cylinder.
    profile: Vec<(f64, f64)>,
    /// Committed displacement of each slave node at the last commit.
    anchor: Vec<[f64; 3]>,
    /// Committed friction force on each slave node at the last commit.
    held: Vec<[f64; 3]>,
}

/// Wall contact terms at a configuration.
#[derive(Debug, Clone)]
pub struct WallTerms {
    /// The residual contribution (penalty gradient minus friction force),
    /// length `3 * node_count`.
    pub force: Vec<f64>,
    /// The Jacobian contribution.
    pub tangent: Coo,
    /// Slave nodes in contact.
    pub active: usize,
    /// Of those, how many are sliding (saturated at `μN`).
    pub slipping: usize,
    /// Largest penetration depth.
    pub max_penetration: f64,
    /// Sum of normal reaction magnitudes `κ·pen`.
    pub total_reaction: f64,
    /// Friction force per slave node (slave order), for committing history.
    pub friction_force: Vec<[f64; 3]>,
}

/// The contact outcome of a wall at a converged configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallSummary {
    /// Slave nodes in contact.
    pub active: usize,
    /// Sliding nodes (`0` when no friction is configured).
    pub slipping: usize,
    /// Largest penetration depth.
    pub max_penetration: f64,
    /// Sum of normal reaction magnitudes.
    pub total_reaction: f64,
}

fn transverse(axis: usize) -> (usize, usize) {
    match axis {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    }
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn project_tangent(v: &[f64; 3], n: &[f64; 3]) -> [f64; 3] {
    let c = dot(v, n);
    [v[0] - c * n[0], v[1] - c * n[1], v[2] - c * n[2]]
}

impl RadialWall {
    /// A frictionless wall about the line parallel to `axis` through the
    /// transverse point `center` (the two coordinates other than `axis`, in
    /// cyclic order: axis 0 → `(y, z)`, 1 → `(z, x)`, 2 → `(x, y)`).
    ///
    /// # Errors
    /// [`WallError`] for a bad axis, radius or centre.
    pub fn new(
        axis: usize,
        center: [f64; 2],
        radius: f64,
        side: WallSide,
        slave: impl IntoIterator<Item = usize>,
    ) -> Result<Self, WallError> {
        if axis > 2 {
            return Err(WallError::AxisOutOfRange(axis));
        }
        if !radius.is_finite() || radius <= 0.0 {
            return Err(WallError::InvalidRadius(radius));
        }
        if !center.iter().all(|c| c.is_finite()) {
            return Err(WallError::InvalidCenter);
        }
        let slave: Vec<usize> = slave.into_iter().collect();
        let n = slave.len();
        Ok(Self {
            axis,
            center,
            radius,
            side,
            slave,
            friction: None,
            profile: Vec::new(),
            anchor: vec![[0.0; 3]; n],
            held: vec![[0.0; 3]; n],
        })
    }

    /// Adds Coulomb friction (`μ`, tangential penalty stiffness).
    pub fn with_friction(mut self, friction: FrictionConfig) -> Self {
        self.friction = Some(friction);
        self
    }

    /// Gives the lumen an axial profile: `points` are `(z, Δr)` pairs, the
    /// radius at axial coordinate `z` being `radius + Δr(z)` by linear
    /// interpolation, flat beyond the ends. A stenosis, a taper or a step
    /// between two ring groups is a profile.
    ///
    /// The slope enters the contact gradient (the wall pushes axially where
    /// it is inclined) and the Jacobian's rank-one part; the profile is
    /// piecewise linear, so the radius' second derivative — an axial-axial
    /// Hessian term — is zero except at the kinks, where it is not modelled.
    ///
    /// # Errors
    /// [`WallError::InvalidProfile`] unless the `z` values are finite and
    /// strictly increasing and the offsets finite.
    pub fn with_profile(mut self, points: Vec<(f64, f64)>) -> Result<Self, WallError> {
        let finite = points.iter().all(|(z, d)| z.is_finite() && d.is_finite());
        let increasing = points.windows(2).all(|w| w[0].0 < w[1].0);
        if !finite || !increasing {
            return Err(WallError::InvalidProfile);
        }
        self.profile = points;
        Ok(self)
    }

    /// The wall radius at axial coordinate `z`.
    pub fn radius_at(&self, z: f64) -> f64 {
        self.radius + self.offset_and_slope(z).0
    }

    /// `dR/dz` at axial coordinate `z` (zero beyond the profile's ends).
    pub fn slope_at(&self, z: f64) -> f64 {
        self.offset_and_slope(z).1
    }

    fn offset_and_slope(&self, z: f64) -> (f64, f64) {
        let pts = &self.profile;
        match pts.len() {
            0 => (0.0, 0.0),
            1 => (pts[0].1, 0.0),
            _ => {
                if z <= pts[0].0 {
                    return (pts[0].1, 0.0);
                }
                if z >= pts[pts.len() - 1].0 {
                    return (pts[pts.len() - 1].1, 0.0);
                }
                let i = pts.partition_point(|p| p.0 <= z) - 1;
                let (z0, d0) = pts[i];
                let (z1, d1) = pts[i + 1];
                let slope = (d1 - d0) / (z1 - z0);
                (d0 + slope * (z - z0), slope)
            }
        }
    }

    /// The cylinder axis.
    pub fn axis(&self) -> usize {
        self.axis
    }

    /// The wall radius.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// The slave nodes.
    pub fn slave_nodes(&self) -> &[usize] {
        &self.slave
    }

    /// The same wall with a different radius and unchanged friction history —
    /// the vessel stage of a deployment moves the wall between stages.
    ///
    /// # Errors
    /// [`WallError::InvalidRadius`].
    pub fn with_radius(mut self, radius: f64) -> Result<Self, WallError> {
        if !radius.is_finite() || radius <= 0.0 {
            return Err(WallError::InvalidRadius(radius));
        }
        self.radius = radius;
        Ok(self)
    }

    /// Residual, Jacobian and statistics at configuration `u`, with normal
    /// penalty stiffness `penalty`.
    ///
    /// # Errors
    /// [`MeshError`] if `u` has the wrong length or a slave node is outside
    /// the mesh.
    pub fn terms<E: ReferenceElement + crate::mesh::ElementFamily>(
        &self,
        mesh: &Mesh<E>,
        u: &[f64],
        penalty: f64,
    ) -> Result<WallTerms, MeshError> {
        let x = mesh.positions(u)?;
        let (p, q) = transverse(self.axis);
        let mut out = WallTerms {
            force: vec![0.0; mesh.dof_count()],
            tangent: Coo::new(),
            active: 0,
            slipping: 0,
            max_penetration: 0.0,
            total_reaction: 0.0,
            friction_force: vec![[0.0; 3]; self.slave.len()],
        };
        let s = match self.side {
            WallSide::Inside => 1.0,
            WallSide::Outside => -1.0,
        };
        for (k, &node) in self.slave.iter().enumerate() {
            if node >= mesh.node_count() {
                return Err(MeshError::NodeIndexOutOfRange {
                    element: usize::MAX,
                    local_node: node,
                    node,
                    node_count: mesh.node_count(),
                });
            }
            let xa = x[node].to_array();
            let (rp, rq) = (xa[p] - self.center[0], xa[q] - self.center[1]);
            let rho = rp.hypot(rq);
            if rho < TINY {
                continue;
            }
            let z = xa[self.axis];
            let pen = s * (rho - self.radius_at(z));
            if pen <= 0.0 {
                continue;
            }
            // Gradient of the penetration: radial in the transverse plane, and
            // axial where the wall is inclined. Its norm is 1 for a straight
            // cylinder, which makes every expression below reduce to the
            // straight-wall one.
            let mut g = [0.0; 3];
            let mut rh = [0.0; 3];
            rh[p] = rp / rho;
            rh[q] = rq / rho;
            g[p] = s * rh[p];
            g[q] = s * rh[q];
            g[self.axis] = -s * self.slope_at(z);
            let gn = dot(&g, &g).sqrt();
            // Unit normal, for projecting tangents.
            let n = [g[0] / gn, g[1] / gn, g[2] / gn];
            let big_n = penalty * pen;
            let normal_force = big_n * gn;
            out.active += 1;
            out.max_penetration = out.max_penetration.max(pen);
            out.total_reaction += normal_force;

            let dofs = [mesh.dof(node, 0), mesh.dof(node, 1), mesh.dof(node, 2)];
            for i in 0..3 {
                out.force[dofs[i]] += big_n * g[i];
            }
            // H = κ [ g gᵀ + pen·s/ρ (I - r̂r̂ᵀ)_transverse ].
            for i in 0..3 {
                for j in 0..3 {
                    let curv = if (i == p || i == q) && (j == p || j == q) {
                        let id = if i == j { 1.0 } else { 0.0 };
                        s * pen / rho * (id - rh[i] * rh[j])
                    } else {
                        0.0
                    };
                    out.tangent
                        .push(dofs[i], dofs[j], penalty * (g[i] * g[j] + curv));
                }
            }

            let Some(fc) = self.friction else { continue };
            if fc.mu == 0.0 || fc.tangential_stiffness == 0.0 {
                continue;
            }
            let un = [u[dofs[0]], u[dofs[1]], u[dofs[2]]];
            let du = [
                un[0] - self.anchor[k][0],
                un[1] - self.anchor[k][1],
                un[2] - self.anchor[k][2],
            ];
            let ds = project_tangent(&du, &n);
            let held = project_tangent(&self.held[k], &n);
            let kt = fc.tangential_stiffness;
            let trial = [
                held[0] - kt * ds[0],
                held[1] - kt * ds[1],
                held[2] - kt * ds[2],
            ];
            let tn = dot(&trial, &trial).sqrt();
            let bound = fc.mu * normal_force;
            let slip = tn > bound;
            let f = if slip && tn > TINY {
                let sc = bound / tn;
                [trial[0] * sc, trial[1] * sc, trial[2] * sc]
            } else {
                trial
            };
            if slip {
                out.slipping += 1;
            }
            out.friction_force[k] = f;
            // Friction is an applied force on the node: residual -= f.
            for i in 0..3 {
                out.force[dofs[i]] -= f[i];
            }
            // r = ... - f(u), so the Jacobian carries -∂f/∂u.
            // P_t as a 3x3 on the tangent plane.
            let pt = |i: usize, j: usize| {
                let id = if i == j { 1.0 } else { 0.0 };
                id - n[i] * n[j]
            };
            if slip && tn > TINY {
                let t_hat = [trial[0] / tn, trial[1] / tn, trial[2] / tn];
                for i in 0..3 {
                    for j in 0..3 {
                        // -∂f/∂u = (bound/|trial|) k_t (P_t - t t^T) P_t ...
                        let rad = pt(i, j) - t_hat[i] * t_hat[j];
                        let mut v = bound / tn * kt * rad;
                        // ... minus the μN coupling: ∂f_i/∂u_j ∋ t_i μ κ n_j.
                        v -= t_hat[i] * fc.mu * penalty * gn * g[j];
                        out.tangent.push(dofs[i], dofs[j], v);
                    }
                }
            } else {
                for i in 0..3 {
                    for j in 0..3 {
                        out.tangent.push(dofs[i], dofs[j], kt * pt(i, j));
                    }
                }
            }
        }
        Ok(out)
    }

    /// The summary at a converged configuration.
    ///
    /// # Errors
    /// As [`terms`](Self::terms).
    pub fn summary<E: ReferenceElement + crate::mesh::ElementFamily>(
        &self,
        mesh: &Mesh<E>,
        u: &[f64],
        penalty: f64,
    ) -> Result<WallSummary, MeshError> {
        let t = self.terms(mesh, u, penalty)?;
        Ok(WallSummary {
            active: t.active,
            slipping: t.slipping,
            max_penetration: t.max_penetration,
            total_reaction: t.total_reaction,
        })
    }

    /// The wall with friction history advanced to the converged displacement
    /// `u`: anchors reset to the present displacements and the held friction
    /// force set to what this configuration produced (zero off contact).
    ///
    /// # Errors
    /// As [`terms`](Self::terms).
    pub fn committed<E: ReferenceElement + crate::mesh::ElementFamily>(
        &self,
        mesh: &Mesh<E>,
        u: &[f64],
        penalty: f64,
    ) -> Result<Self, MeshError> {
        let t = self.terms(mesh, u, penalty)?;
        let mut next = self.clone();
        for (k, &node) in self.slave.iter().enumerate() {
            next.anchor[k] = [
                u[mesh.dof(node, 0)],
                u[mesh.dof(node, 1)],
                u[mesh.dof(node, 2)],
            ];
            next.held[k] = t.friction_force[k];
        }
        Ok(next)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::hex_box;
    use crate::solver::{solve_static, ContactConfig, SolveOptions};
    use crate::ContactPairing;
    use tpt_med_tissue::{NeoHookeanParams, TissueModel};

    fn dense(coo: &Coo, n: usize) -> Vec<Vec<f64>> {
        let mut m = vec![vec![0.0; n]; n];
        for i in 0..coo.len() {
            m[coo.rows[i]][coo.cols[i]] += coo.vals[i];
        }
        m
    }

    #[test]
    fn constructor_validates() {
        assert_eq!(
            RadialWall::new(3, [0.0; 2], 1.0, WallSide::Inside, []).err(),
            Some(WallError::AxisOutOfRange(3))
        );
        assert!(matches!(
            RadialWall::new(2, [0.0; 2], 0.0, WallSide::Inside, []).err(),
            Some(WallError::InvalidRadius(_))
        ));
        assert_eq!(
            RadialWall::new(2, [f64::NAN, 0.0], 1.0, WallSide::Inside, []).err(),
            Some(WallError::InvalidCenter)
        );
    }

    #[test]
    fn active_set_follows_geometry_and_side() {
        let mesh = hex_box(1, 1, 1, 2.0, 2.0, 2.0).expect("box");
        let u = vec![0.0; mesh.dof_count()];
        let all: Vec<usize> = (0..mesh.node_count()).collect();
        // Axis z, centre at the origin: corner nodes sit at rho = 0, 2, 2, 2.83.
        let inside = RadialWall::new(2, [0.0; 2], 1.9, WallSide::Inside, all.clone()).unwrap();
        let t = inside.terms(&mesh, &u, 100.0).unwrap();
        // x=2,y=0 and x=0,y=2 (2 nodes each, two z levels) plus the diagonal pair.
        assert_eq!(t.active, 6);
        assert!((t.max_penetration - (8.0f64.sqrt() - 1.9)).abs() < 1e-12);
        // Outside: only the node on the axis (rho = 0) is skipped as degenerate;
        // all others are outside a wall of radius 5 and penetrate it.
        let outside = RadialWall::new(2, [0.0; 2], 5.0, WallSide::Outside, all).unwrap();
        let t = outside.terms(&mesh, &u, 100.0).unwrap();
        assert_eq!(t.active, 6);
    }

    /// Central-difference Jacobian of the wall force vs the assembled tangent.
    fn fd_check(wall: &RadialWall, mesh: &crate::Hex8Mesh, u: &[f64], penalty: f64, tol: f64) {
        let n = mesh.dof_count();
        let t = wall.terms(mesh, u, penalty).unwrap();
        let k = dense(&t.tangent, n);
        let h = 1e-7;
        let mut worst = 0.0f64;
        let mut scale = 0.0f64;
        for j in 0..n {
            let (mut up, mut um) = (u.to_vec(), u.to_vec());
            up[j] += h;
            um[j] -= h;
            let fp = wall.terms(mesh, &up, penalty).unwrap().force;
            let fm = wall.terms(mesh, &um, penalty).unwrap().force;
            for i in 0..n {
                let fd = (fp[i] - fm[i]) / (2.0 * h);
                worst = worst.max((fd - k[i][j]).abs());
                scale = scale.max(fd.abs());
            }
        }
        assert!(worst <= tol * scale.max(1.0), "worst {worst} scale {scale}");
    }

    fn perturbed(n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| 0.05 * ((i * 7 % 11) as f64 / 11.0 - 0.5))
            .collect()
    }

    #[test]
    fn normal_jacobian_matches_finite_differences() {
        let mesh = hex_box(1, 1, 1, 2.0, 2.0, 2.0).expect("box");
        let wall =
            RadialWall::new(2, [0.3, 0.2], 1.7, WallSide::Inside, 0..mesh.node_count()).unwrap();
        fd_check(&wall, &mesh, &perturbed(mesh.dof_count()), 250.0, 1e-6);
        let sheath =
            RadialWall::new(2, [0.3, 0.2], 3.0, WallSide::Outside, 0..mesh.node_count()).unwrap();
        fd_check(&sheath, &mesh, &perturbed(mesh.dof_count()), 250.0, 1e-6);
    }

    fn single_node(mesh: &crate::Hex8Mesh, at: [f64; 3]) -> usize {
        (0..mesh.node_count())
            .find(|&i| {
                let p = mesh.nodes()[i].to_array();
                (0..3).all(|k| (p[k] - at[k]).abs() < 1e-12)
            })
            .expect("node")
    }

    #[test]
    fn friction_sticks_then_saturates_at_the_coulomb_bound() {
        let mesh = hex_box(1, 1, 1, 2.0, 2.0, 2.0).expect("box");
        let node = single_node(&mesh, [2.0, 0.0, 0.0]);
        let fc = FrictionConfig::new(0.5, 1000.0).unwrap();
        let wall = RadialWall::new(2, [0.0; 2], 1.9, WallSide::Inside, [node])
            .unwrap()
            .with_friction(fc);
        let penalty = 100.0; // N = 100 * 0.1 = 10, bound = 5
        let mut u = vec![0.0; mesh.dof_count()];
        let dz = mesh.dof(node, 2);
        u[dz] = 0.002; // trial 2 < 5: stick
        let t = wall.terms(&mesh, &u, penalty).unwrap();
        assert_eq!(t.slipping, 0);
        assert!(
            (t.friction_force[0][2] + 2.0).abs() < 1e-9,
            "{:?}",
            t.friction_force
        );
        u[dz] = 0.01; // trial 10 > 5: slip, saturated
        let t = wall.terms(&mesh, &u, penalty).unwrap();
        assert_eq!(t.slipping, 1);
        assert!(
            (t.friction_force[0][2] + 5.0).abs() < 1e-9,
            "{:?}",
            t.friction_force
        );
    }

    #[test]
    fn friction_state_is_carried_across_commits() {
        let mesh = hex_box(1, 1, 1, 2.0, 2.0, 2.0).expect("box");
        let node = single_node(&mesh, [2.0, 0.0, 0.0]);
        let fc = FrictionConfig::new(0.5, 1000.0).unwrap();
        let wall = RadialWall::new(2, [0.0; 2], 1.9, WallSide::Inside, [node])
            .unwrap()
            .with_friction(fc);
        let mut u = vec![0.0; mesh.dof_count()];
        u[mesh.dof(node, 2)] = 0.002;
        let next = wall.committed(&mesh, &u, 100.0).unwrap();
        // At the same configuration the slip since the new anchor is zero, so
        // the force is exactly the held one: stick memory, not a relaxation.
        let t = next.terms(&mesh, &u, 100.0).unwrap();
        assert!((t.friction_force[0][2] + 2.0).abs() < 1e-9);
        // A further small slide adds to the held force.
        u[mesh.dof(node, 2)] = 0.003;
        let t = next.terms(&mesh, &u, 100.0).unwrap();
        assert!((t.friction_force[0][2] + 3.0).abs() < 1e-9);
    }

    #[test]
    fn friction_jacobian_matches_finite_differences() {
        let mesh = hex_box(1, 1, 1, 2.0, 2.0, 2.0).expect("box");
        let all = 0..mesh.node_count();
        let fc_stick = FrictionConfig::new(5.0, 20.0).unwrap();
        let stick = RadialWall::new(2, [0.3, 0.2], 1.7, WallSide::Inside, all.clone())
            .unwrap()
            .with_friction(fc_stick);
        let u0 = perturbed(mesh.dof_count());
        // Commit at u0, then move a little: stick branch.
        let stick = stick.committed(&mesh, &u0, 250.0).unwrap();
        let u1: Vec<f64> = u0
            .iter()
            .enumerate()
            .map(|(i, v)| v + 1e-3 * (i % 5) as f64)
            .collect();
        fd_check(&stick, &mesh, &u1, 250.0, 5e-3);
        // Slip branch: tiny mu so the bound is exceeded.
        let fc_slip = FrictionConfig::new(1e-3, 1.0e4).unwrap();
        let slip = RadialWall::new(2, [0.3, 0.2], 1.7, WallSide::Inside, all)
            .unwrap()
            .with_friction(fc_slip);
        let slip = slip.committed(&mesh, &u0, 250.0).unwrap();
        let t = slip.terms(&mesh, &u1, 250.0).unwrap();
        assert!(t.slipping > 0, "slip branch not exercised");
        fd_check(&slip, &mesh, &u1, 250.0, 5e-3);
    }

    #[test]
    fn block_pressed_into_the_wall_reaches_a_symmetric_equilibrium() {
        let l = 10.0;
        let mesh = hex_box(2, 2, 2, l, l, l).expect("box");
        let mut dirichlet = Vec::new();
        for n in mesh.face_nodes(0, false) {
            dirichlet.push((mesh.dof(n, 0), 0.0));
        }
        for n in mesh.face_nodes(1, false) {
            dirichlet.push((mesh.dof(n, 1), 0.0));
        }
        for n in mesh.face_nodes(2, false) {
            dirichlet.push((mesh.dof(n, 2), 0.0));
        }
        let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.49, d1: 0.5 });
        let wall_nodes: Vec<usize> = (0..mesh.node_count()).collect();
        let radius = 13.0; // the corner column sits at 14.14
        let wall = RadialWall::new(2, [0.0; 2], radius, WallSide::Inside, wall_nodes).unwrap();
        let pairing = ContactPairing::inactive();
        let solve = |penalty: f64| {
            solve_static(
                &mesh,
                &model,
                &vec![0.0; mesh.dof_count()],
                &dirichlet,
                &SolveOptions::default(),
                Some(ContactConfig {
                    pairing: &pairing,
                    penalty,
                    friction: None,
                    radial: Some(&wall),
                }),
            )
            .expect("converges")
        };
        let soft = solve(50.0);
        let stiff = solve(500.0);
        let s = soft.contact.as_ref().unwrap().wall.unwrap();
        let t = stiff.contact.as_ref().unwrap().wall.unwrap();
        assert!(s.active > 0 && t.active > 0);
        // Stiffer penalty -> smaller penetration, roughly inversely.
        assert!(t.max_penetration < 0.3 * s.max_penetration, "{s:?} {t:?}");
        // The wall acts on the corner column only (nodes with rho <= R inactive).
        let all =
            RadialWall::new(2, [0.0; 2], radius, WallSide::Inside, 0..mesh.node_count()).unwrap();
        let terms = all.terms(&mesh, &stiff.displacement, 500.0).unwrap();
        assert_eq!(terms.active, t.active);
        // x <-> y mirror symmetry of the solved displacement field.
        for n in 0..mesh.node_count() {
            let p = mesh.nodes()[n].to_array();
            let m = (0..mesh.node_count())
                .find(|&j| {
                    let q = mesh.nodes()[j].to_array();
                    (q[0] - p[1]).abs() < 1e-9
                        && (q[1] - p[0]).abs() < 1e-9
                        && (q[2] - p[2]).abs() < 1e-9
                })
                .unwrap();
            let (a, b) = (mesh.dof(n, 0), mesh.dof(m, 1));
            assert!(
                (stiff.displacement[a] - stiff.displacement[b]).abs() < 1e-6,
                "asymmetry at node {n}"
            );
        }
        // Normal reaction balances: the inward force the wall applies equals
        // the force the body transmits (residual of the free dofs is ~0 at
        // convergence, so the wall's total radial force on the body is
        // carried by the symmetry-plane reactions).
        assert!(stiff.residual_norm < 1e-6);
    }

    #[test]
    fn profile_interpolates_clamps_and_validates() {
        let w = RadialWall::new(2, [0.0; 2], 3.0, WallSide::Inside, [])
            .unwrap()
            .with_profile(vec![(0.0, 0.0), (2.0, -0.5), (3.0, -0.5)])
            .unwrap();
        assert!((w.radius_at(-1.0) - 3.0).abs() < 1e-15);
        assert!((w.radius_at(1.0) - 2.75).abs() < 1e-15);
        assert!((w.slope_at(1.0) + 0.25).abs() < 1e-15);
        assert!((w.radius_at(2.5) - 2.5).abs() < 1e-15);
        assert_eq!(w.slope_at(2.5), 0.0);
        assert!((w.radius_at(9.0) - 2.5).abs() < 1e-15);
        let bad = RadialWall::new(2, [0.0; 2], 3.0, WallSide::Inside, [])
            .unwrap()
            .with_profile(vec![(1.0, 0.0), (1.0, 0.1)]);
        assert_eq!(bad.err(), Some(WallError::InvalidProfile));
        let nan = RadialWall::new(2, [0.0; 2], 3.0, WallSide::Inside, [])
            .unwrap()
            .with_profile(vec![(0.0, f64::NAN)]);
        assert_eq!(nan.err(), Some(WallError::InvalidProfile));
    }

    #[test]
    fn inclined_wall_pushes_axially_and_its_jacobian_matches_finite_differences() {
        let mesh = hex_box(1, 1, 1, 2.0, 2.0, 2.0).expect("box");
        // A taper over the block's axial extent: the lumen narrows with z.
        let wall = RadialWall::new(2, [0.3, 0.2], 1.7, WallSide::Inside, 0..mesh.node_count())
            .unwrap()
            .with_profile(vec![(-1.0, 0.2), (3.0, -0.4)])
            .unwrap();
        let u = perturbed(mesh.dof_count());
        let t = wall.terms(&mesh, &u, 250.0).unwrap();
        // The slope (dR/dz < 0) makes the penetration gradient point toward
        // +z, so the penalty residual carries a positive axial component.
        let axial: f64 = (0..mesh.node_count())
            .map(|n| t.force[mesh.dof(n, 2)])
            .sum();
        assert!(axial > 0.0, "axial wall force {axial}");
        fd_check(&wall, &mesh, &u, 250.0, 1e-6);
        // A straight wall has no axial force from the normal contact.
        let straight =
            RadialWall::new(2, [0.3, 0.2], 1.7, WallSide::Inside, 0..mesh.node_count()).unwrap();
        let ts = straight.terms(&mesh, &u, 250.0).unwrap();
        let axial_s: f64 = (0..mesh.node_count())
            .map(|n| ts.force[mesh.dof(n, 2)])
            .sum();
        assert!(axial_s.abs() < 1e-12, "{axial_s}");
    }
}
