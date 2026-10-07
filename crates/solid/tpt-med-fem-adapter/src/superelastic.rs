//! Souza–Auricchio-style 3-D superelastic point model (RFC 0013, first slice).
//!
//! Off by default (`superelastic` feature). This is the constitutive core of
//! the stent Level 3 ladder: a per-quadrature-point internal-variable model
//! with an explicit return mapping. It is **not** yet a deployment driver —
//! per-quadrature state storage, the crown-ring fixtures and the contact
//! deployment stages are the later slices RFC 0013 sequences.
//!
//! # Formulation
//!
//! Moderate-strain, St. Venant–Kirchhoff elastic response in the
//! Green–Lagrange strain `E = ½(FᵀF − I)`, with the transformation strain
//! `E_tr = ξ·ε_L·N` carried as an internal variable:
//!
//! ```text
//! E_e = E − ξ ε_L N,    S = 2G(ξ) E_e + λ(ξ) tr(E_e) I,    P = F S
//! ```
//!
//! `N` is the deviatoric transformation direction, normalised so a uniaxial
//! state has `N = diag(1, −½, −½)` (its axial component is exactly `ε_L`
//! per unit `ξ`); it is fixed at activation from the deviatoric stress and
//! released when `ξ` returns to zero, as the RFC specifies. Moduli mix
//! linearly between austenite and martensite with `ξ`.
//!
//! **Departure from the RFC text, recorded here:** the RFC sketched
//! Neo-Hookean mixing. This slice uses the quadratic St. Venant–Kirchhoff
//! form instead, because it makes the uniaxial reduction track the 1-D stent
//! model's partition `ε = σ/E + ξ ε_L` (the verification identity the RFC
//! leads with) and nitinol's transformation strain is moderate (≲ 8 %). A
//! finite-strain (Neo-Hookean / Hencky) swap-in is a localised change to
//! `elastic_stress` should a validation case demand it.
//!
//! # Kinetics
//!
//! Rate-independent, stress-driven, sharing the 1-D plateau bounds. The drive
//! is the deviatoric stress projected on `N` (equal to the axial stress in a
//! uniaxial state, and signed, so unloading through zero is distinguished
//! from reverse loading); while no direction exists it is the von Mises
//! equivalent. Forward transformation runs on `σ_f(ξ)` (cosine interface
//! `σ_ms → σ_mf`), reverse on `σ_r(ξ)` (`σ_af → σ_as`); between the two the
//! response is elastic at the current `ξ`. The return mapping solves
//! `drive(ξ) = σ_bound(ξ)` by bisection on the bracket the monotone
//! structure guarantees, with saturation at `ξ = 0` / `ξ = 1`.
//!
//! One honest difference from the 1-D model: its reverse branch interpolates
//! *strain* linearly in `ξ`, whereas this model's reverse branch is
//! stress-consistent. The two loops agree on the plateau stresses, elastic
//! slopes and closure; the verification compares them with a stated
//! tolerance rather than claiming a bit-level identity.

use crate::assembly::Constitutive;
use tpt_med_geometry::Mat3;

/// Superelastic material parameters (MPa, dimensionless strains).
///
/// The kinetic bounds mirror `tpt-med-stents`' `NitinolParams`; `poisson` is
/// the 3-D-only elastic constant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SuperelasticParams {
    /// Austenite Young's modulus (MPa).
    pub e_austenite: f64,
    /// Martensite Young's modulus (MPa).
    pub e_martensite: f64,
    /// Uniaxial transformation strain ε_L.
    pub transformation_strain: f64,
    /// Forward start stress σ_ms (MPa).
    pub sigma_ms: f64,
    /// Forward finish stress σ_mf (MPa).
    pub sigma_mf: f64,
    /// Reverse start stress σ_as (MPa).
    pub sigma_as: f64,
    /// Reverse finish stress σ_af (MPa).
    pub sigma_af: f64,
    /// Poisson ratio (shared by both phases).
    pub poisson: f64,
}

impl Default for SuperelasticParams {
    /// The 1-D crate's defaults plus `ν = 0.3`.
    fn default() -> Self {
        Self {
            e_austenite: 55_000.0,
            e_martensite: 28_000.0,
            transformation_strain: 0.05,
            sigma_ms: 480.0,
            sigma_mf: 560.0,
            sigma_as: 380.0,
            sigma_af: 260.0,
            poisson: 0.3,
        }
    }
}

/// Errors from the superelastic model.
#[derive(Debug, Clone, PartialEq)]
pub enum SuperelasticError {
    /// The plateau bounds are not ordered `σ_af < σ_as < σ_ms < σ_mf`
    /// (a non-monotone or hysteresis-free plateau has no return map).
    NonMonotonePlateau,
    /// A parameter is non-finite or out of its physical range.
    InvalidParameter(&'static str),
    /// The deformation gradient contains a non-finite entry.
    NonFiniteInput,
    /// The return mapping found no admissible `ξ` update.
    ReturnMappingFailed {
        /// Drive residual at the bracket end that failed.
        drive: f64,
        /// The martensite fraction the update started from.
        xi: f64,
    },
}

impl core::fmt::Display for SuperelasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonMonotonePlateau => write!(
                f,
                "plateau bounds must satisfy sigma_af < sigma_as < sigma_ms < sigma_mf"
            ),
            Self::InvalidParameter(what) => write!(f, "invalid superelastic parameter: {what}"),
            Self::NonFiniteInput => write!(f, "deformation gradient is not finite"),
            Self::ReturnMappingFailed { drive, xi } => write!(
                f,
                "return mapping failed from xi = {xi} (drive residual {drive})"
            ),
        }
    }
}

impl std::error::Error for SuperelasticError {}

/// Per-quadrature-point internal variables, advanced by the caller (the
/// deployment driver) — the `Constitutive` trait itself is stateless.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SuperelasticState {
    /// Martensite volume fraction ξ ∈ [0, 1].
    pub xi: f64,
    /// Transformation direction `N` (zero while `ξ = 0`).
    pub direction: Mat3,
}

impl Default for SuperelasticState {
    fn default() -> Self {
        Self {
            xi: 0.0,
            direction: Mat3::ZERO,
        }
    }
}

/// Result of one state update.
#[derive(Debug, Clone, Copy)]
pub struct SuperelasticUpdate {
    /// First Piola–Kirchhoff stress at the updated state.
    pub piola: Mat3,
    /// The updated internal variables.
    pub state: SuperelasticState,
    /// Von Mises equivalent of the deviatoric second Piola stress (MPa).
    pub equivalent_stress: f64,
}

/// A validated Souza–Auricchio-style model.
#[derive(Debug, Clone, Copy)]
pub struct SouzaAuricchio {
    params: SuperelasticParams,
}

const BISECTION_ITERS: usize = 100;

fn dev(m: &Mat3) -> Mat3 {
    let t = m.trace() / 3.0;
    *m - Mat3::IDENTITY * t
}

fn double_dot(a: &Mat3, b: &Mat3) -> f64 {
    let mut s = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            s += a.at(i, j) * b.at(i, j);
        }
    }
    s
}

fn von_mises(s: &Mat3) -> f64 {
    let d = dev(s);
    (1.5 * double_dot(&d, &d)).sqrt()
}

/// Signed drive: the deviatoric stress projected on the transformation
/// direction (`σ` for a uniaxial state), so unloading through zero and into
/// reverse loading changes sign instead of being folded positive. Falls back
/// to the von Mises equivalent while no direction exists yet.
fn drive(s: &Mat3, n: &Mat3) -> f64 {
    if double_dot(n, n) > 0.0 {
        double_dot(&dev(s), n)
    } else {
        von_mises(s)
    }
}

/// Cosine-interface stress between plateau bounds, `ξ ∈ [0, 1]`.
fn interface(xi: f64, lo: f64, hi: f64) -> f64 {
    let x = (1.0 - 2.0 * xi.clamp(0.0, 1.0)).acos();
    lo + (hi - lo) * x / core::f64::consts::PI
}

impl SouzaAuricchio {
    /// Validate and build the model.
    ///
    /// # Errors
    /// [`SuperelasticError::NonMonotonePlateau`] unless
    /// `0 < σ_af < σ_as < σ_ms < σ_mf`; [`SuperelasticError::InvalidParameter`]
    /// for non-finite or non-positive moduli/strain or `ν ∉ (−1, ½)`.
    pub fn new(params: SuperelasticParams) -> Result<Self, SuperelasticError> {
        let p = &params;
        let all = [
            p.e_austenite,
            p.e_martensite,
            p.transformation_strain,
            p.sigma_ms,
            p.sigma_mf,
            p.sigma_as,
            p.sigma_af,
            p.poisson,
        ];
        if all.iter().any(|v| !v.is_finite()) {
            return Err(SuperelasticError::InvalidParameter("non-finite value"));
        }
        if p.e_austenite <= 0.0 || p.e_martensite <= 0.0 {
            return Err(SuperelasticError::InvalidParameter(
                "moduli must be positive",
            ));
        }
        if p.transformation_strain <= 0.0 {
            return Err(SuperelasticError::InvalidParameter(
                "transformation strain must be positive",
            ));
        }
        if !(p.poisson > -1.0 && p.poisson < 0.5) {
            return Err(SuperelasticError::InvalidParameter(
                "poisson must lie in (-1, 0.5)",
            ));
        }
        if !(0.0 < p.sigma_af
            && p.sigma_af < p.sigma_as
            && p.sigma_as < p.sigma_ms
            && p.sigma_ms < p.sigma_mf)
        {
            return Err(SuperelasticError::NonMonotonePlateau);
        }
        Ok(Self { params })
    }

    /// The validated parameters.
    pub fn params(&self) -> &SuperelasticParams {
        &self.params
    }

    fn young(&self, xi: f64) -> f64 {
        self.params.e_austenite + xi * (self.params.e_martensite - self.params.e_austenite)
    }

    /// Second Piola stress at strain `e`, fraction `xi`, direction `n`.
    fn elastic_stress(&self, e: &Mat3, xi: f64, n: &Mat3) -> Mat3 {
        let nu = self.params.poisson;
        let young = self.young(xi);
        let g = young / (2.0 * (1.0 + nu));
        let lam = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let ee = *e - *n * (xi * self.params.transformation_strain);
        ee * (2.0 * g) + Mat3::IDENTITY * (lam * ee.trace())
    }

    /// Advance the internal variables for a deformation gradient `f` from the
    /// committed `state`, returning the stress at the updated state.
    ///
    /// # Errors
    /// [`SuperelasticError::NonFiniteInput`] for a non-finite `f`;
    /// [`SuperelasticError::ReturnMappingFailed`] if the residual is
    /// non-finite at a bracket end.
    pub fn update(
        &self,
        f: &Mat3,
        state: &SuperelasticState,
    ) -> Result<SuperelasticUpdate, SuperelasticError> {
        if f.to_array().iter().any(|v| !v.is_finite()) {
            return Err(SuperelasticError::NonFiniteInput);
        }
        let p = &self.params;
        let e = (f.transpose().mul_mat(f) - Mat3::IDENTITY) * 0.5;
        let (xi_n, n_n) = (state.xi, state.direction);

        let trial = self.elastic_stress(&e, xi_n, &n_n);
        let drive_trial = drive(&trial, &n_n);

        let forward_trigger = xi_n < 1.0 && drive_trial > interface(xi_n, p.sigma_ms, p.sigma_mf);
        let reverse_trigger = xi_n > 0.0 && drive_trial < interface(xi_n, p.sigma_af, p.sigma_as);

        let (xi, n) = if forward_trigger {
            // A direction is needed to transform at all: take it from the
            // trial deviatoric stress when starting from the austenite.
            let n = if xi_n == 0.0 {
                let q = von_mises(&trial);
                if q > 0.0 {
                    dev(&trial) * (1.5 / q)
                } else {
                    Mat3::ZERO
                }
            } else {
                n_n
            };
            let g = |xi: f64| {
                drive(&self.elastic_stress(&e, xi, &n), &n) - interface(xi, p.sigma_ms, p.sigma_mf)
            };
            (self.solve(&g, xi_n, 1.0, xi_n)?, n)
        } else if reverse_trigger {
            let g = |xi: f64| {
                drive(&self.elastic_stress(&e, xi, &n_n), &n_n)
                    - interface(xi, p.sigma_af, p.sigma_as)
            };
            let xi = self.solve(&g, 0.0, xi_n, xi_n)?;
            (xi, if xi == 0.0 { Mat3::ZERO } else { n_n })
        } else {
            (xi_n, n_n)
        };

        let s = self.elastic_stress(&e, xi, &n);
        Ok(SuperelasticUpdate {
            piola: f.mul_mat(&s),
            state: SuperelasticState { xi, direction: n },
            equivalent_stress: von_mises(&s),
        })
    }

    /// Root of a residual that is non-increasing in `ξ` on `[lo, hi]`:
    /// the upper end if the residual is still non-negative there, the lower
    /// end if already non-positive, bisection otherwise.
    fn solve(
        &self,
        g: &dyn Fn(f64) -> f64,
        lo: f64,
        hi: f64,
        xi_from: f64,
    ) -> Result<f64, SuperelasticError> {
        let (g_lo, g_hi) = (g(lo), g(hi));
        if !g_lo.is_finite() || !g_hi.is_finite() {
            return Err(SuperelasticError::ReturnMappingFailed {
                drive: if g_lo.is_finite() { g_hi } else { g_lo },
                xi: xi_from,
            });
        }
        if g_hi >= 0.0 {
            return Ok(hi);
        }
        if g_lo <= 0.0 {
            return Ok(lo);
        }
        let (mut a, mut b) = (lo, hi);
        for _ in 0..BISECTION_ITERS {
            let mid = 0.5 * (a + b);
            if g(mid) > 0.0 {
                a = mid;
            } else {
                b = mid;
            }
        }
        Ok(0.5 * (a + b))
    }

    /// A [`Constitutive`] view of this model at a committed `state`: each
    /// call runs the full return mapping from `state` and discards the
    /// updated variables, so a central-difference
    /// [`material_tangent`](crate::material_tangent) is the consistent
    /// algorithmic tangent **away from the switching point**. At a committed
    /// state sitting exactly on a kinetic bound the response has a kink
    /// (plateau on one side, elastic on the other), and a central difference
    /// there averages the two slopes — a driver should difference one-sided
    /// in the loading direction, or commit states strictly inside a branch.
    /// The deployment driver commits states between
    /// load steps.
    pub fn at_state(&self, state: SuperelasticState) -> SuperelasticAt<'_> {
        SuperelasticAt { model: self, state }
    }
}

/// [`SouzaAuricchio`] frozen at a committed internal state.
///
/// `first_piola` cannot return a `Result`; a failed update yields a
/// NaN-filled stress so the Newton loop fails loudly instead of continuing
/// on a clamped value.
#[derive(Debug, Clone, Copy)]
pub struct SuperelasticAt<'a> {
    model: &'a SouzaAuricchio,
    state: SuperelasticState,
}

impl Constitutive for SuperelasticAt<'_> {
    fn first_piola(&self, f: &Mat3) -> Mat3 {
        match self.model.update(f, &self.state) {
            Ok(u) => u.piola,
            Err(_) => Mat3::from_array([f64::NAN; 9]),
        }
    }
}

/// Errors from advancing a [`SuperelasticField`].
#[derive(Debug)]
pub enum FieldError {
    /// The mesh could not supply a deformation gradient (wrong DOF count,
    /// degenerate element).
    Mesh(crate::mesh::MeshError),
    /// A point update failed.
    Point(SuperelasticError),
}

impl core::fmt::Display for FieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Mesh(e) => write!(f, "{e}"),
            Self::Point(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for FieldError {}

/// Committed internal variables at every quadrature point of a mesh, as a
/// [`Constitutive`] the assembly can use directly.
///
/// This is the per-point state storage RFC 0013 puts in the deployment
/// driver: the assembly asks for stress through
/// [`Constitutive::first_piola_at`], and each answer is the return mapping
/// from that point's **committed** state at the trial deformation. Nothing is
/// advanced during a Newton solve; [`commit`](Self::commit) advances every
/// point once, after a converged step, so a cut-back step restarts from
/// untouched state.
///
/// Used through the plain [`Constitutive::first_piola`] (no point identity)
/// it returns a NaN stress rather than guess a state.
#[derive(Debug, Clone)]
pub struct SuperelasticField {
    model: SouzaAuricchio,
    points: usize,
    quadrature_order: usize,
    states: Vec<SuperelasticState>,
}

impl SuperelasticField {
    /// A fresh (austenite, unloaded) field for `mesh`, indexed at
    /// `quadrature_order` — which must equal the assembly's.
    pub fn new<E: tpt_fem_element::ReferenceElement + crate::mesh::ElementFamily>(
        model: SouzaAuricchio,
        mesh: &crate::mesh::Mesh<E>,
        quadrature_order: usize,
    ) -> Self {
        let points = E::quadrature_rule(quadrature_order).points.len();
        Self {
            model,
            points,
            quadrature_order,
            states: vec![SuperelasticState::default(); mesh.element_count() * points],
        }
    }

    /// The model.
    pub fn model(&self) -> &SouzaAuricchio {
        &self.model
    }

    /// The committed state at `(element, point)`.
    pub fn state(&self, element: usize, point: usize) -> Option<&SuperelasticState> {
        if point >= self.points {
            return None;
        }
        self.states.get(element * self.points + point)
    }

    /// `(min, max, mean)` of the committed martensite fraction.
    pub fn martensite_summary(&self) -> (f64, f64, f64) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        let mut sum = 0.0;
        for s in &self.states {
            lo = lo.min(s.xi);
            hi = hi.max(s.xi);
            sum += s.xi;
        }
        if self.states.is_empty() {
            return (0.0, 0.0, 0.0);
        }
        (lo, hi, sum / self.states.len() as f64)
    }

    /// Advance every point's committed state to the displacement `u`.
    ///
    /// All-or-nothing: if any point fails, no state changes.
    ///
    /// # Errors
    /// [`FieldError`] if a deformation gradient cannot be formed or a point
    /// update fails.
    pub fn commit<E: tpt_fem_element::ReferenceElement + crate::mesh::ElementFamily>(
        &mut self,
        mesh: &crate::mesh::Mesh<E>,
        u: &[f64],
    ) -> Result<(), FieldError> {
        if u.len() != mesh.dof_count() {
            return Err(FieldError::Mesh(crate::mesh::MeshError::DofCountMismatch {
                expected: mesh.dof_count(),
                found: u.len(),
            }));
        }
        let rule = E::quadrature_rule(self.quadrature_order);
        let mut next = Vec::with_capacity(self.states.len());
        for e in 0..mesh.element_count() {
            for (q, xi) in rule.points.iter().enumerate() {
                let f = crate::assembly::element_deformation_gradient(mesh, e, u, xi).ok_or(
                    FieldError::Mesh(crate::mesh::MeshError::DegenerateElement {
                        element: e,
                        jacobian_determinant: mesh.jacobian(e, xi).det(),
                    }),
                )?;
                let committed = &self.states[e * self.points + q];
                let up = self
                    .model
                    .update(&f, committed)
                    .map_err(FieldError::Point)?;
                next.push(up.state);
            }
        }
        self.states = next;
        Ok(())
    }
}

impl Constitutive for SuperelasticField {
    fn first_piola(&self, _f: &Mat3) -> Mat3 {
        Mat3::from_array([f64::NAN; 9])
    }

    fn first_piola_at(&self, element: usize, point: usize, f: &Mat3) -> Mat3 {
        let Some(state) = self.state(element, point) else {
            return Mat3::from_array([f64::NAN; 9]);
        };
        match self.model.update(f, state) {
            Ok(u) => u.piola,
            Err(_) => Mat3::from_array([f64::NAN; 9]),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::material_tangent;
    use tpt_med_stents::{NitinolParams, SuperelasticState as OneD};

    fn model() -> SouzaAuricchio {
        SouzaAuricchio::new(SuperelasticParams::default()).expect("valid")
    }

    /// Lateral stretch that zeroes the lateral nominal stress at axial
    /// stretch `lam`, from a committed `state`.
    fn lateral_stretch(m: &SouzaAuricchio, lam: f64, state: &SuperelasticState) -> f64 {
        let (mut a, mut b) = (0.5, 1.2);
        for _ in 0..100 {
            let mid = 0.5 * (a + b);
            let s = m
                .update(&Mat3::diagonal([lam, mid, mid]), state)
                .expect("update")
                .piola
                .at(1, 1);
            // Lateral stress rises with the lateral stretch.
            if s > 0.0 {
                b = mid;
            } else {
                a = mid;
            }
        }
        0.5 * (a + b)
    }

    /// Uniaxial stretch `1 + e` with free sides, from a committed state.
    fn uniaxial(m: &SouzaAuricchio, e: f64, state: &SuperelasticState) -> SuperelasticUpdate {
        let lam = 1.0 + e;
        let mu = lateral_stretch(m, lam, state);
        m.update(&Mat3::diagonal([lam, mu, mu]), state)
            .expect("update")
    }

    fn sweep(m: &SouzaAuricchio, path: &[f64]) -> Vec<(f64, f64, f64)> {
        let mut st = SuperelasticState::default();
        path.iter()
            .map(|&e| {
                let u = uniaxial(m, e, &st);
                st = u.state;
                (e, u.piola.at(0, 0), u.state.xi)
            })
            .collect()
    }

    fn up_down(peak: f64, n: usize) -> Vec<f64> {
        let mut v: Vec<f64> = (0..=n).map(|i| peak * i as f64 / n as f64).collect();
        v.extend((0..n).rev().map(|i| peak * i as f64 / n as f64));
        v
    }

    fn loaded_to(m: &SouzaAuricchio, e: f64) -> SuperelasticState {
        let mut st = SuperelasticState::default();
        for k in 1..=100 {
            st = uniaxial(m, e * k as f64 / 100.0, &st).state;
        }
        st
    }

    #[test]
    fn validation_rejects_bad_parameters() {
        let ok = SuperelasticParams::default();
        let bad = |f: &dyn Fn(&mut SuperelasticParams)| {
            let mut p = ok;
            f(&mut p);
            SouzaAuricchio::new(p).err()
        };
        assert_eq!(
            bad(&|p| p.sigma_as = p.sigma_ms + 1.0),
            Some(SuperelasticError::NonMonotonePlateau)
        );
        assert_eq!(
            bad(&|p| p.sigma_mf = p.sigma_ms),
            Some(SuperelasticError::NonMonotonePlateau)
        );
        assert!(matches!(
            bad(&|p| p.poisson = 0.5),
            Some(SuperelasticError::InvalidParameter(_))
        ));
        assert!(matches!(
            bad(&|p| p.e_martensite = 0.0),
            Some(SuperelasticError::InvalidParameter(_))
        ));
        assert!(matches!(
            bad(&|p| p.sigma_ms = f64::NAN),
            Some(SuperelasticError::InvalidParameter(_))
        ));
    }

    #[test]
    fn non_finite_gradient_is_an_error() {
        let f = Mat3::diagonal([f64::NAN, 1.0, 1.0]);
        assert_eq!(
            model().update(&f, &SuperelasticState::default()).err(),
            Some(SuperelasticError::NonFiniteInput)
        );
    }

    #[test]
    fn elastic_austenite_slope_is_the_austenite_modulus() {
        let m = model();
        let r = sweep(&m, &[0.0, 0.001, 0.002]);
        let slope = (r[2].1 - r[1].1) / 0.001;
        assert!((slope - 55_000.0).abs() / 55_000.0 < 0.01, "slope {slope}");
        assert!(r.iter().all(|x| x.2 == 0.0));
    }

    #[test]
    fn return_mapping_hits_the_hand_computed_midplateau() {
        // xi = 1/2 sits at the cosine-interface midpoint, sigma_f = mean of
        // the forward bounds. Hand-build the strain that puts it there.
        let m = model();
        let p = m.params();
        let xi = 0.5;
        let sigma = 0.5 * (p.sigma_ms + p.sigma_mf);
        let e_axial_gl = sigma / m.young(xi) + xi * p.transformation_strain;
        let lam = (1.0 + 2.0 * e_axial_gl).sqrt();
        let st = loaded_to(&m, lam - 1.0);
        assert!((st.xi - xi).abs() < 5e-3, "xi {}", st.xi);
        let u = uniaxial(&m, lam - 1.0, &st);
        assert!(
            (u.equivalent_stress - sigma).abs() < 2.0,
            "{}",
            u.equivalent_stress
        );
    }

    #[test]
    fn loop_closes_and_unload_sits_below_load() {
        let m = model();
        let r = sweep(&m, &up_down(0.07, 140));
        let last = r.last().unwrap();
        assert!(last.1.abs() < 1e-6 && last.2 == 0.0, "residual {last:?}");
        // Hysteresis: at 3 % strain the unloading stress is below loading.
        let load = r
            .iter()
            .take(141)
            .find(|x| (x.0 - 0.03).abs() < 1e-9)
            .unwrap();
        let unload = r
            .iter()
            .skip(141)
            .find(|x| (x.0 - 0.03).abs() < 1e-9)
            .unwrap();
        assert!(
            unload.1 < load.1 - 50.0,
            "load {} unload {}",
            load.1,
            unload.1
        );
        assert!(unload.2 > 0.0, "still partly martensitic on unload at 3 %");
    }

    #[test]
    fn partial_unloading_is_elastic() {
        let m = model();
        let st = loaded_to(&m, 0.03);
        assert!(st.xi > 0.0 && st.xi < 1.0);
        // A small unload sits inside the elastic window between the bounds.
        let u = uniaxial(&m, 0.0295, &st);
        assert_eq!(u.state.xi, st.xi);
    }

    #[test]
    fn forward_plateau_tracks_the_one_d_model() {
        // The forward kinetics share the 1-D plateau bounds, so the 3-D
        // loading stress must follow the 1-D one on the plateau within the
        // documented mixture-modulus / strain-measure tolerance. The reverse
        // interface differs by design (strain-linear in 1-D) and is not
        // compared point-for-point.
        let m = model();
        let q = NitinolParams::default();
        let path: Vec<f64> = (0..=140).map(|i| 0.07 * i as f64 / 140.0).collect();
        let r = sweep(&m, &path);
        let mut one = OneD::new();
        let mut worst = 0.0f64;
        let mut compared = 0;
        for &(e, s3, _) in &r {
            let s1 = one.strain_to_stress(e, &q);
            if s1 > q.sigma_ms + 5.0 && s1 < q.sigma_mf - 5.0 {
                worst = worst.max((s3 - s1).abs() / s1);
                compared += 1;
            }
        }
        assert!(compared > 20, "plateau points compared: {compared}");
        assert!(worst < 0.05, "forward plateau deviates {worst}");
    }

    #[test]
    fn algorithmic_tangent_is_elastic_off_plateau_and_soft_on_it() {
        let m = model();
        let p = m.params();
        let nu = p.poisson;
        let c1111 = p.e_austenite * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let a = material_tangent(
            &m.at_state(SuperelasticState::default()),
            &Mat3::IDENTITY,
            1e-6,
        );
        assert!(
            (a[0][0][0][0] - c1111).abs() / c1111 < 1e-3,
            "{}",
            a[0][0][0][0]
        );

        // On the plateau the physically meaningful tangent is the free-side
        // uniaxial one (the fixed-lateral component A1111 also carries the
        // stiff volumetric response): differencing the converged uniaxial
        // response from the committed state is the algorithmic slope, and it
        // must sit near the plateau slope, far below the elastic modulus.
        let st = loaded_to(&m, 0.03);
        let h = 1e-5;
        let hi = uniaxial(&m, 0.03 + h, &st).piola.at(0, 0);
        // One-sided: the backward side of a committed plateau state is the
        // elastic unloading branch, so a central difference there averages
        // two different slopes (the kink recorded on `at_state`).
        let lo = uniaxial(&m, 0.03, &st).piola.at(0, 0);
        let slope = (hi - lo) / h;
        assert!(
            slope > 0.0 && slope < 0.1 * p.e_austenite,
            "plateau slope {slope}"
        );
    }
}
