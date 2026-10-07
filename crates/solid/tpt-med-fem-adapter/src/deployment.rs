//! Staged, state-carrying deployment driver (RFC 0013, second slice).
//!
//! Off by default (`superelastic` feature). Where [`crate::solve_load_path`]
//! walks one stateless path, a superelastic deployment is a *sequence* of
//! stages over history-dependent material, so this driver owns what a stage
//! sequence needs and a single solve does not:
//!
//! - the per-quadrature-point committed state
//!   ([`SuperelasticField`](crate::superelastic::SuperelasticField)), advanced
//!   only after a converged increment so a cut-back increment restarts from
//!   untouched state;
//! - the current displacement and the set of prescribed DOFs;
//! - displacement-controlled stages with cutback ([`Deployment3D::prescribe`]
//!   — the **crimp** and any displacement-driven expansion);
//! - a force-controlled **release** ([`Deployment3D::release`]): constraints
//!   are removed and their reaction forces ramped to zero, so superelastic
//!   recovery drives the expansion rather than a prescribed path.
//!
//! Planar axis-aligned contact (the adapter's existing pairing, optionally
//! with friction) can ride along on any stage. **Not in this slice:** a
//! radial rigid-cylinder vessel — the pairing resolves gaps along one
//! coordinate axis, so a cylindrical wall needs a new radial pairing — and
//! the crown-ring / two-group fixtures. Both are the next RFC 0013 slices.

use crate::assembly::{internal_force, AssemblyOptions};
use crate::mesh::{ElementFamily, Mesh, MeshError};
use crate::solver::{newton_from, ContactConfig, SolveError, SolveOptions};
use crate::superelastic::{FieldError, SouzaAuricchio, SuperelasticField};
use std::collections::BTreeMap;
use tpt_fem_element::ReferenceElement;

/// Controls for a stage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StageOptions {
    /// Equal increments from the stage start to its end.
    pub steps: usize,
    /// Maximum bisections of one increment after a convergence failure.
    pub max_cutbacks: usize,
}

impl Default for StageOptions {
    fn default() -> Self {
        Self {
            steps: 10,
            max_cutbacks: 6,
        }
    }
}

/// Errors from the deployment driver.
#[derive(Debug)]
pub enum DeploymentError {
    /// An increment (and every bisection of it) failed to converge.
    StepNotConverged {
        /// Stage progress in `[0, 1]` that could not be reached.
        progress: f64,
        /// Bisections attempted.
        cutbacks: usize,
        /// The final solver failure.
        source: SolveError,
    },
    /// The solve failed for a reason retrying smaller cannot fix (singular
    /// system, bad contact pairing, inverted element).
    Solve(SolveError),
    /// Advancing the committed state failed.
    State(FieldError),
    /// Mesh-level failure evaluating reactions.
    Mesh(MeshError),
    /// A DOF index was outside the mesh.
    DofOutOfRange(usize),
    /// `release` named a DOF that is not currently prescribed.
    NotPrescribed(usize),
}

impl std::fmt::Display for DeploymentError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::StepNotConverged {
                progress,
                cutbacks,
                source,
            } => write!(
                f,
                "stage increment to {progress} did not converge after {cutbacks} cutbacks: {source}"
            ),
            Self::Solve(e) => write!(f, "{e}"),
            Self::State(e) => write!(f, "{e}"),
            Self::Mesh(e) => write!(f, "{e}"),
            Self::DofOutOfRange(d) => write!(f, "dof {d} is outside the mesh"),
            Self::NotPrescribed(d) => write!(f, "dof {d} is not prescribed"),
        }
    }
}

impl std::error::Error for DeploymentError {}

/// One converged increment of a stage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IncrementReport {
    /// Stage progress reached, in `(0, 1]`.
    pub progress: f64,
    /// Newton iterations the increment took.
    pub newton_iterations: usize,
    /// Smallest / largest / mean committed martensite fraction afterwards.
    pub martensite: (f64, f64, f64),
}

/// What a finished stage did.
#[derive(Debug, Clone, PartialEq)]
pub struct StageReport {
    /// Every converged increment, in order.
    pub increments: Vec<IncrementReport>,
    /// Total bisections taken across the stage.
    pub cutbacks: usize,
}

/// Staged 3-D superelastic deployment on a caller-supplied mesh.
#[derive(Debug)]
pub struct Deployment3D<'m, E: ReferenceElement + ElementFamily> {
    mesh: &'m Mesh<E>,
    field: SuperelasticField,
    opts: SolveOptions,
    u: Vec<f64>,
    prescribed: BTreeMap<usize, f64>,
}

impl<'m, E: ReferenceElement + ElementFamily> Deployment3D<'m, E> {
    /// A fresh deployment: undeformed, austenitic, with `fixed` DOFs held at
    /// the given values for the whole run unless later released.
    ///
    /// # Errors
    /// [`DeploymentError::DofOutOfRange`] for a DOF outside the mesh.
    pub fn new(
        mesh: &'m Mesh<E>,
        model: SouzaAuricchio,
        opts: SolveOptions,
        fixed: &[(usize, f64)],
    ) -> Result<Self, DeploymentError> {
        let n = mesh.dof_count();
        let mut prescribed = BTreeMap::new();
        let mut u = vec![0.0; n];
        for &(d, v) in fixed {
            if d >= n {
                return Err(DeploymentError::DofOutOfRange(d));
            }
            prescribed.insert(d, v);
            u[d] = v;
        }
        let field = SuperelasticField::new(model, mesh, opts.assembly.quadrature_order);
        Ok(Self {
            mesh,
            field,
            opts,
            u,
            prescribed,
        })
    }

    /// The current displacement.
    pub fn displacement(&self) -> &[f64] {
        &self.u
    }

    /// The committed per-point state.
    pub fn field(&self) -> &SuperelasticField {
        &self.field
    }

    /// Internal-force reaction summed over `dofs` at the current state.
    ///
    /// # Errors
    /// [`DeploymentError::Mesh`] if assembly fails; [`DeploymentError::DofOutOfRange`].
    pub fn reaction(&self, dofs: &[usize]) -> Result<f64, DeploymentError> {
        let f = self.internal_force()?;
        let mut sum = 0.0;
        for &d in dofs {
            sum += *f.get(d).ok_or(DeploymentError::DofOutOfRange(d))?;
        }
        Ok(sum)
    }

    fn internal_force(&self) -> Result<Vec<f64>, DeploymentError> {
        let asm: AssemblyOptions = self.opts.assembly;
        internal_force(self.mesh, &self.field, &self.u, &asm).map_err(DeploymentError::Mesh)
    }

    /// A displacement-controlled stage: each target DOF ramps linearly from
    /// its present value to the target; all other prescribed DOFs hold.
    /// Used for the crimp, and for any displacement-driven expansion.
    ///
    /// # Errors
    /// [`DeploymentError::DofOutOfRange`]; [`DeploymentError::StepNotConverged`]
    /// when cutback is exhausted; [`DeploymentError::Solve`] /
    /// [`DeploymentError::State`] for non-retryable failures.
    pub fn prescribe(
        &mut self,
        targets: &[(usize, f64)],
        stage: StageOptions,
        contact: Option<ContactConfig<'_>>,
    ) -> Result<StageReport, DeploymentError> {
        let n = self.mesh.dof_count();
        for &(d, _) in targets {
            if d >= n {
                return Err(DeploymentError::DofOutOfRange(d));
            }
        }
        let starts: Vec<(usize, f64, f64)> = targets
            .iter()
            .map(|&(d, v)| (d, self.prescribed.get(&d).copied().unwrap_or(self.u[d]), v))
            .collect();
        let zero = vec![0.0; n];
        self.march(stage, contact, |t, held| {
            for &(d, a, b) in &starts {
                held.insert(d, a + (b - a) * t);
            }
            zero.clone()
        })
    }

    /// A force-controlled **release**: the named prescribed DOFs are removed
    /// from the constraint set and replaced by their current reaction forces,
    /// which are ramped to zero — the body is let go and superelastic
    /// recovery moves it. This is the RFC's "remove the crimp constraints".
    ///
    /// # Errors
    /// [`DeploymentError::NotPrescribed`] for a DOF that is not held; the
    /// same stage errors as [`prescribe`](Self::prescribe).
    pub fn release(
        &mut self,
        dofs: &[usize],
        stage: StageOptions,
        contact: Option<ContactConfig<'_>>,
    ) -> Result<StageReport, DeploymentError> {
        for &d in dofs {
            if !self.prescribed.contains_key(&d) {
                return Err(DeploymentError::NotPrescribed(d));
            }
        }
        let f_int = self.internal_force()?;
        let mut reaction = vec![0.0; self.mesh.dof_count()];
        for &d in dofs {
            reaction[d] = f_int[d];
            self.prescribed.remove(&d);
        }
        self.march(stage, contact, |t, _| {
            reaction.iter().map(|r| r * (1.0 - t)).collect()
        })
    }

    /// Walk a stage: `setup(progress, held)` fills the Dirichlet values
    /// for progress `t` and returns the external load at that progress.
    fn march(
        &mut self,
        stage: StageOptions,
        contact: Option<ContactConfig<'_>>,
        mut setup: impl FnMut(f64, &mut BTreeMap<usize, f64>) -> Vec<f64>,
    ) -> Result<StageReport, DeploymentError> {
        let mut report = StageReport {
            increments: Vec::new(),
            cutbacks: 0,
        };
        if stage.steps == 0 {
            return Ok(report);
        }
        let dt = 1.0 / stage.steps as f64;
        let mut reached = 0.0f64;
        while reached < 1.0 - 1e-12 {
            let nominal = (reached + dt).min(1.0);
            let mut attempt = nominal;
            let mut cutbacks = 0usize;
            loop {
                let mut held = self.prescribed.clone();
                let load = setup(attempt, &mut held);
                let dirichlet: Vec<(usize, f64)> = held.iter().map(|(&d, &v)| (d, v)).collect();
                match newton_from(
                    self.mesh,
                    &self.field,
                    &load,
                    &dirichlet,
                    &self.opts,
                    contact,
                    Some(&self.u),
                ) {
                    Ok(r) => {
                        // Commit only a converged increment, onto a copy, so a
                        // failed commit leaves the committed state untouched.
                        let mut next = self.field.clone();
                        next.commit(self.mesh, &r.displacement)
                            .map_err(DeploymentError::State)?;
                        self.field = next;
                        self.u = r.displacement;
                        self.prescribed = held;
                        reached = attempt;
                        report.increments.push(IncrementReport {
                            progress: attempt,
                            newton_iterations: r.newton_iterations,
                            martensite: self.field.martensite_summary(),
                        });
                        break;
                    }
                    Err(e @ SolveError::NotConverged { .. }) => {
                        if cutbacks >= stage.max_cutbacks {
                            return Err(DeploymentError::StepNotConverged {
                                progress: attempt,
                                cutbacks,
                                source: e,
                            });
                        }
                        cutbacks += 1;
                        report.cutbacks += 1;
                        // Bisect the remaining distance of this increment so
                        // progress is retained; the next increment returns to
                        // the nominal size.
                        attempt = reached + (attempt - reached) / 2.0;
                    }
                    Err(e) => return Err(DeploymentError::Solve(e)),
                }
            }
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::{hex_box, Hex8Mesh};
    use crate::solver::SolveOptions;
    use crate::superelastic::{SuperelasticParams, SuperelasticState};
    use tpt_med_geometry::Mat3;

    const L: f64 = 10.0;

    fn model() -> SouzaAuricchio {
        SouzaAuricchio::new(SuperelasticParams::default()).expect("valid")
    }

    /// Symmetry-plane block under uniaxial control along y: the faces at
    /// x = 0, y = 0, z = 0 are slid shut, the top face is the prescribed one.
    /// Returns the deployment plus the top y-DOFs.
    fn block(
        mesh: &Hex8Mesh,
        opts: SolveOptions,
    ) -> (Deployment3D<'_, tpt_fem_element::Hex8>, Vec<usize>) {
        let mut fixed = Vec::new();
        for n in mesh.face_nodes(0, false) {
            fixed.push((mesh.dof(n, 0), 0.0));
        }
        for n in mesh.face_nodes(1, false) {
            fixed.push((mesh.dof(n, 1), 0.0));
        }
        for n in mesh.face_nodes(2, false) {
            fixed.push((mesh.dof(n, 2), 0.0));
        }
        let top: Vec<usize> = mesh
            .face_nodes(1, true)
            .iter()
            .map(|&n| mesh.dof(n, 1))
            .collect();
        for &d in &top {
            fixed.push((d, 0.0));
        }
        (
            Deployment3D::new(mesh, model(), opts, &fixed).expect("deployment"),
            top,
        )
    }

    fn opts() -> SolveOptions {
        SolveOptions::default()
    }

    /// Point-model nominal axial stress through a committed strain sequence,
    /// the reference the mesh result must reproduce.
    fn point_path(strains: &[f64]) -> Vec<f64> {
        let m = model();
        let mut st = SuperelasticState::default();
        strains
            .iter()
            .map(|&e| {
                let lam = 1.0 + e;
                let (mut a, mut b) = (0.5, 1.2);
                for _ in 0..100 {
                    let mid = 0.5 * (a + b);
                    let s = m
                        .update(&Mat3::diagonal([mid, lam, mid]), &st)
                        .expect("update")
                        .piola
                        .at(0, 0);
                    if s > 0.0 {
                        b = mid;
                    } else {
                        a = mid;
                    }
                }
                let mu = 0.5 * (a + b);
                let u = m
                    .update(&Mat3::diagonal([mu, lam, mu]), &st)
                    .expect("update");
                st = u.state;
                u.piola.at(1, 1)
            })
            .collect()
    }

    fn nominal(dep: &Deployment3D<'_, tpt_fem_element::Hex8>, top: &[usize]) -> f64 {
        dep.reaction(top).expect("reaction") / (L * L)
    }

    #[test]
    fn single_hex_loop_reproduces_the_point_model() {
        let mesh = hex_box(1, 1, 1, L, L, L).expect("box");
        let (mut dep, top) = block(&mesh, opts());
        let stage = StageOptions {
            steps: 6,
            max_cutbacks: 8,
        };
        let marks = [0.03, 0.07, 0.03, 0.0];
        let mut got = Vec::new();
        for &e in &marks {
            let targets: Vec<(usize, f64)> = top.iter().map(|&d| (d, e * L)).collect();
            dep.prescribe(&targets, stage, None).expect("stage");
            got.push(nominal(&dep, &top));
        }
        // The point reference must walk the same committed sub-steps; with a
        // rate-independent law only the kink placement can differ, so refine.
        let mut path = Vec::new();
        let mut at = 0.0;
        let mut ends = Vec::new();
        for &e in &marks {
            for k in 1..=6 {
                path.push(at + (e - at) * k as f64 / 6.0);
            }
            at = e;
            ends.push(path.len() - 1);
        }
        let reference = point_path(&path);
        for (g, &i) in got.iter().zip(&ends) {
            let r = reference[i];
            assert!(
                (g - r).abs() <= 2e-3 * r.abs().max(100.0),
                "mesh {g} vs point {r}"
            );
        }
        // Loop closed: no residual stress or martensite after unloading.
        assert!(got[3].abs() < 1e-3, "residual {}", got[3]);
        assert_eq!(dep.field().martensite_summary().1, 0.0);
        // Hysteresis through the real solver.
        assert!(got[2] < got[0] - 50.0, "load {} unload {}", got[0], got[2]);
    }

    #[test]
    fn refined_block_is_uniform_and_matches_the_single_element() {
        let one = hex_box(1, 1, 1, L, L, L).expect("box");
        let many = hex_box(2, 2, 2, L, L, L).expect("box");
        let (mut a, ta) = block(&one, opts());
        let (mut b, tb) = block(&many, opts());
        let stage = StageOptions {
            steps: 6,
            max_cutbacks: 8,
        };
        for e in [0.03, 0.07] {
            let t1: Vec<_> = ta.iter().map(|&d| (d, e * L)).collect();
            let t2: Vec<_> = tb.iter().map(|&d| (d, e * L)).collect();
            a.prescribe(&t1, stage, None).expect("one");
            b.prescribe(&t2, stage, None).expect("many");
        }
        let (sa, sb) = (nominal(&a, &ta), nominal(&b, &tb));
        assert!((sa - sb).abs() < 1e-3 * sa.abs(), "{sa} vs {sb}");
        let (lo, hi, _) = b.field().martensite_summary();
        assert!(hi - lo < 1e-6, "xi not uniform: {lo}..{hi}");
        assert!(hi > 0.99, "fully transformed at 7 %: {hi}");
    }

    #[test]
    fn cutback_reaches_the_same_state_as_fine_stepping() {
        let mesh = hex_box(1, 1, 1, L, L, L).expect("box");
        let (mut coarse, tc) = block(&mesh, opts());
        let (mut fine, tf) = block(&mesh, opts());
        let t1: Vec<_> = tc.iter().map(|&d| (d, 0.07 * L)).collect();
        let t2: Vec<_> = tf.iter().map(|&d| (d, 0.07 * L)).collect();
        let rep = coarse
            .prescribe(
                &t1,
                StageOptions {
                    steps: 1,
                    max_cutbacks: 10,
                },
                None,
            )
            .expect("coarse");
        fine.prescribe(
            &t2,
            StageOptions {
                steps: 28,
                max_cutbacks: 4,
            },
            None,
        )
        .expect("fine");
        let (sc, sf) = (nominal(&coarse, &tc), nominal(&fine, &tf));
        assert!((sc - sf).abs() < 5e-3 * sf.abs(), "{sc} vs {sf}");
        assert!(rep.increments.last().unwrap().progress == 1.0);
    }

    #[test]
    fn release_recovers_the_superelastic_strain() {
        let mesh = hex_box(2, 2, 2, L, L, L).expect("box");
        let (mut dep, top) = block(&mesh, opts());
        let t: Vec<_> = top.iter().map(|&d| (d, 0.07 * L)).collect();
        dep.prescribe(
            &t,
            StageOptions {
                steps: 6,
                max_cutbacks: 8,
            },
            None,
        )
        .expect("load");
        let loaded = nominal(&dep, &top);
        assert!(loaded > 500.0, "loaded stress {loaded}");

        // Let go: reaction forces ramp to zero and the body recovers.
        dep.release(
            &top,
            StageOptions {
                steps: 20,
                max_cutbacks: 8,
            },
            None,
        )
        .expect("release");
        let y_top = dep.displacement()[top[0]];
        assert!(
            y_top.abs() < 1e-4 * L,
            "top displacement after release {y_top}"
        );
        assert!(nominal(&dep, &top).abs() < 1e-3);
        assert_eq!(dep.field().martensite_summary().1, 0.0);
    }

    #[test]
    fn stage_input_is_validated() {
        let mesh = hex_box(1, 1, 1, L, L, L).expect("box");
        let (mut dep, top) = block(&mesh, opts());
        let n = mesh.dof_count();
        assert!(matches!(
            dep.prescribe(&[(n, 1.0)], StageOptions::default(), None),
            Err(DeploymentError::DofOutOfRange(_))
        ));
        // A free (never prescribed) DOF cannot be released.
        let free = (0..n)
            .find(|d| !top.contains(d) && !dep.prescribed.contains_key(d))
            .expect("a free dof");
        assert!(matches!(
            dep.release(&[free], StageOptions::default(), None),
            Err(DeploymentError::NotPrescribed(_))
        ));
        assert!(matches!(
            Deployment3D::new(&mesh, model(), opts(), &[(n + 1, 0.0)]),
            Err(DeploymentError::DofOutOfRange(_))
        ));
    }
}
