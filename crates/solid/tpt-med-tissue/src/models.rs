//! Strain-energy models and stress computation.

use tpt_med_geometry::Mat3;

/// Neo-Hookean parameters (compressible penalty form).
///
/// `W = C10 (J^{-2/3} I1 − 3) + (1/D1) (J − 1)²`.
///
/// `D1 → ∞` recovers incompressibility in a penalty sense; for exact
/// incompressibility use `J = 1` deformation fields (as the verification
/// tests do).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NeoHookeanParams {
    /// Shear-term modulus.
    pub c10: f64,
    /// Volumetric penalty coefficient (`1/D1 (J−1)²`).
    pub d1: f64,
}

/// Mooney–Rivlin parameters.
///
/// `W = C10 (J^{-2/3} I1 − 3) + C01 (J^{-4/3} I2 − 3) + (1/D1) (J − 1)²`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MooneyRivlinParams {
    /// First-order `I1` modulus.
    pub c10: f64,
    /// `I2` modulus.
    pub c01: f64,
    /// Volumetric penalty coefficient.
    pub d1: f64,
}

/// Yeoh parameters (three-term in `Ī1 − 3`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct YeohParams {
    /// Linear coefficient.
    pub c1: f64,
    /// Quadratic coefficient.
    pub c2: f64,
    /// Cubic coefficient.
    pub c3: f64,
    /// Volumetric penalty coefficient.
    pub d1: f64,
}

/// Ogden parameters (deviatoric principal-stretch form).
///
/// `W = Σᵢ (2 μᵢ / αᵢ²) (λ̄1^{αᵢ} + λ̄2^{αᵢ} + λ̄3^{αᵢ} − 3) + (1/D1)(J−1)²`
/// with deviatoric stretches `λ̄ = λ J^{−1/3}`.
#[derive(Debug, Clone, PartialEq)]
pub struct OgdenParams {
    /// Moduli per term.
    pub mu: Vec<f64>,
    /// Exponents per term (same length as `mu`).
    pub alpha: Vec<f64>,
    /// Volumetric penalty coefficient.
    pub d1: f64,
}

/// A hyperelastic tissue model.
#[derive(Debug, Clone, PartialEq)]
pub enum TissueModel {
    /// Neo-Hookean (linear in `Ī1`).
    NeoHookean(NeoHookeanParams),
    /// Mooney–Rivlin (linear in `Ī1` and `Ī2`).
    MooneyRivlin(MooneyRivlinParams),
    /// Yeoh (cubic in `Ī1`).
    Yeoh(YeohParams),
    /// Ogden (principal-stretch power series).
    Ogden(OgdenParams),
    /// Holzapfel–Gasser–Ogden (fiber-reinforced arterial wall).
    HolzapfelGasserOgden(crate::HgoParams),
}

impl TissueModel {
    /// Strain energy density `W(F)`.
    pub fn strain_energy(&self, f: &Mat3) -> f64 {
        let f = *f;
        match self {
            TissueModel::NeoHookean(p) => {
                let j = f.det();
                let i1_bar = j.powf(-2.0 / 3.0) * invariant_i1(&f);
                p.c10 * (i1_bar - 3.0) + (j - 1.0).powi(2) / p.d1
            }
            TissueModel::MooneyRivlin(p) => {
                let j = f.det();
                let i1_bar = j.powf(-2.0 / 3.0) * invariant_i1(&f);
                let i2_bar = j.powf(-4.0 / 3.0) * invariant_i2(&f);
                p.c10 * (i1_bar - 3.0) + p.c01 * (i2_bar - 3.0) + (j - 1.0).powi(2) / p.d1
            }
            TissueModel::Yeoh(p) => {
                let j = f.det();
                let s = j.powf(-2.0 / 3.0) * invariant_i1(&f) - 3.0;
                p.c1 * s + p.c2 * s * s + p.c3 * s * s * s + (j - 1.0).powi(2) / p.d1
            }
            TissueModel::Ogden(p) => {
                let j = f.det();
                let stretches = principal_stretches(&f);
                let mut w = (j - 1.0).powi(2) / p.d1;
                let j13 = j.powf(1.0 / 3.0);
                for (mu, alpha) in p.mu.iter().zip(&p.alpha) {
                    let sum: f64 = stretches.iter().map(|l| (l / j13).powf(*alpha)).sum();
                    w += 2.0 * mu / (alpha * alpha) * (sum - 3.0);
                }
                w
            }
            TissueModel::HolzapfelGasserOgden(p) => p.strain_energy(&f),
        }
    }

    /// First Piola–Kirchhoff stress `P = ∂W/∂F` (analytic).
    pub fn first_piola(&self, f: &Mat3) -> Mat3 {
        let f = *f;
        let f_inv_t = f.inverse().map(|i| i.transpose()).unwrap_or(Mat3::ZERO);
        let j = f.det();
        let i1 = invariant_i1(&f);
        let vol_term = |d1: f64| {
            if j > EPS_F64 {
                2.0 * j * (j - 1.0) / d1 * f_inv_t
            } else {
                Mat3::ZERO
            }
        };

        match self {
            TissueModel::NeoHookean(p) => {
                2.0 * p.c10 * j.powf(-2.0 / 3.0) * (f - (i1 / 3.0) * f_inv_t) + vol_term(p.d1)
            }
            TissueModel::MooneyRivlin(p) => {
                let i2 = invariant_i2(&f);
                let c = f.transpose() * f;
                let iso1 = 2.0 * p.c10 * j.powf(-2.0 / 3.0) * (f - (i1 / 3.0) * f_inv_t);
                let iso2 = 2.0
                    * p.c01
                    * j.powf(-4.0 / 3.0)
                    * (i1 * f - f * c - (2.0 * i2 / 3.0) * f_inv_t);
                iso1 + iso2 + vol_term(p.d1)
            }
            TissueModel::Yeoh(p) => {
                let s = j.powf(-2.0 / 3.0) * i1 - 3.0;
                let dw = p.c1 + 2.0 * p.c2 * s + 3.0 * p.c3 * s * s;
                2.0 * j.powf(-2.0 / 3.0) * dw * (f - (i1 / 3.0) * f_inv_t) + vol_term(p.d1)
            }
            // Ogden/HGO: analytic forms require eigenvector derivatives; the
            // finite-difference reference below is the verified path (the
            // FD test locks the step so results are deterministic).
            TissueModel::Ogden(_) | TissueModel::HolzapfelGasserOgden(_) => {
                self.first_piola_numerical(&f)
            }
        }
    }

    /// Volumetric part of the first Piola–Kirchhoff stress:
    /// `d/dF [(J - 1)^2 / d1] = 2J(J - 1)/d1 * F^{-T}`.
    ///
    /// Split out from [`TissueModel::first_piola`] so a caller can integrate the
    /// volumetric response on its own — which is what selective reduced
    /// integration needs, and what makes the penalty visible as a separable term
    /// rather than fused into `P`.
    ///
    /// All five variants share the identical `(J - 1)^2 / d1` penalty, so the
    /// volumetric first Piola is the same closed form for every one of them. That
    /// is a convenience, not an accident of the derivation: it is what lets this
    /// be one method rather than five.
    ///
    /// Returns zero for `J <= 0`, matching the guard in
    /// [`TissueModel::first_piola`]. An inverted configuration has no valid
    /// volumetric response, and the adapter rejects inversion before asking.
    pub fn volumetric_first_piola(&self, f: &Mat3) -> Mat3 {
        let f = *f;
        let j = f.det();
        if j <= EPS_F64 {
            return Mat3::ZERO;
        }
        let d1 = match self {
            TissueModel::NeoHookean(p) => p.d1,
            TissueModel::MooneyRivlin(p) => p.d1,
            TissueModel::Yeoh(p) => p.d1,
            TissueModel::Ogden(p) => p.d1,
            TissueModel::HolzapfelGasserOgden(p) => p.d1,
        };
        let f_inv_t = f.inverse().map(|i| i.transpose()).unwrap_or(Mat3::ZERO);
        2.0 * j * (j - 1.0) / d1 * f_inv_t
    }

    /// Central-difference reference for `∂W/∂F` (used for Ogden/HGO and to
    /// verify all analytic derivatives).
    pub fn first_piola_numerical(&self, f: &Mat3) -> Mat3 {
        let f = *f;
        let f = &f;
        const H: f64 = 1.0e-6;
        let mut p = Mat3::ZERO;
        for r in 0..3 {
            for c in 0..3 {
                let mut fp = *f;
                fp.set(r, c, fp.at(r, c) + H);
                let mut fm = *f;
                fm.set(r, c, fm.at(r, c) - H);
                p.set(
                    r,
                    c,
                    (self.strain_energy(&fp) - self.strain_energy(&fm)) / (2.0 * H),
                );
            }
        }
        p
    }

    /// Second-order material tangent `A[i][j](k,l) = ∂P_ij/∂F_kl` — the
    /// tensor a nonlinear Newton solver assembles the element stiffness
    /// from.
    ///
    /// Delivery follows the same policy as [`Self::first_piola`]: analytic
    /// for the Neo-Hookean and Yeoh families (which share the
    /// `P_dev = 2β·q·G` structure) plus the model-independent volumetric
    /// penalty; central differences through [`Self::first_piola`] for
    /// Mooney–Rivlin, Ogden and HGO (whose analytic forms need eigenvector
    /// derivatives). The numerical path for those three carries FD-of-FD
    /// round-off (~1e-4 relative), which consumers should treat as the
    /// tangent's own tolerance.
    pub fn material_tangent(&self, f: &Mat3) -> MaterialTangent {
        let f = *f;
        match self {
            TissueModel::NeoHookean(p) => {
                let mut a = deviatoric_tangent_two_beta_q_g(f, p.c10, 0.0);
                add_volumetric_tangent(&mut a, f, p.d1);
                a
            }
            TissueModel::Yeoh(p) => {
                let j = f.det();
                let s = j.powf(-2.0 / 3.0) * invariant_i1(&f) - 3.0;
                let q = p.c1 + 2.0 * p.c2 * s + 3.0 * p.c3 * s * s;
                let q_prime = 2.0 * p.c2 + 6.0 * p.c3 * s;
                let mut a = deviatoric_tangent_two_beta_q_g(f, q, q_prime);
                add_volumetric_tangent(&mut a, f, p.d1);
                a
            }
            _ => self.material_tangent_numerical(&f),
        }
    }

    /// Central-difference reference for [`Self::material_tangent`]:
    /// `A[i][j](k,l) = ΔP_ij / ΔF_kl` with the same step as
    /// [`Self::first_piola_numerical`].
    pub fn material_tangent_numerical(&self, f: &Mat3) -> MaterialTangent {
        let f = *f;
        const H: f64 = 1.0e-6;
        let mut a = zero_tangent();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut fp = f;
                        fp.set(k, l, fp.at(k, l) + H);
                        let mut fm = f;
                        fm.set(k, l, fm.at(k, l) - H);
                        a[i][j].set(
                            k,
                            l,
                            (self.first_piola(&fp).at(i, j) - self.first_piola(&fm).at(i, j))
                                / (2.0 * H),
                        );
                    }
                }
            }
        }
        a
    }

    /// The volumetric part of the material tangent — `d/dF` of
    /// [`Self::volumetric_first_piola`]:
    ///
    /// ```text
    /// A_vol[ij,kl] = (2/d1)·J·[(2J−1)·F⁻¹_lk·F⁻¹_ji − (J−1)·F⁻¹_jk·F⁻¹_li]
    /// ```
    ///
    /// All five models share the `(J−1)²/d1` penalty, so — exactly as with
    /// `volumetric_first_piola` — this is one closed form for every model,
    /// which is what a mixed `u`-`p` formulation needs from the tangent
    /// side. Zero for `J <= 0`, matching the stress-side guard.
    pub fn volumetric_tangent(&self, f: &Mat3) -> MaterialTangent {
        let f = *f;
        let d1 = self.volumetric_d1();
        let mut a = zero_tangent();
        if f.det() <= EPS_F64 {
            return a;
        }
        add_volumetric_tangent(&mut a, f, d1);
        a
    }

    /// Deviatoric energy with the **mean dilatation** `j_bar` substituted
    /// for the pointwise `J` in the isochoric factors, and the volumetric
    /// penalty dropped — the constitutive half of a mixed `u`-`p` or
    /// mean-dilatation (`B-bar`) formulation.
    ///
    /// The point of the substitution: for the workspace's laws the
    /// isochoric factors are the only place the pointwise `J` enters
    /// (e.g. Neo-Hookean's `J^{-2/3} I1`), so replacing `J -> j_bar`
    /// exactly removes the pointwise volumetric coupling a penalty or a
    /// constraint must otherwise fight — while `j_bar = J` reproduces the
    /// plain deviatoric response, which is the identity
    /// `mean_dilatation_energy(f, f.det()) == strain_energy(f) -
    /// (f.det()-1)^2/d1` the test suite pins. The fiber terms of HGO
    /// carry no dilatation factor and pass through unchanged.
    ///
    /// Pair with [`Self::mean_dilatation_first_piola`], whose stress then
    /// needs the constraint term `p * cof(F)` supplied by the
    /// formulation.
    pub fn mean_dilatation_energy(&self, f: &Mat3, j_bar: f64) -> f64 {
        let f = *f;
        let c_mat = f.transpose() * f;
        let i1 = c_mat.trace();
        let j_bar_factor = j_bar.powf(-2.0 / 3.0);
        match self {
            TissueModel::NeoHookean(p) => p.c10 * (j_bar_factor * i1 - 3.0),
            TissueModel::MooneyRivlin(p) => {
                let i2 = invariant_i2(&f);
                p.c10 * (j_bar_factor * i1 - 3.0) + p.c01 * (j_bar_factor * j_bar_factor * i2 - 3.0)
            }
            TissueModel::Yeoh(p) => {
                let s = j_bar_factor * i1 - 3.0;
                p.c1 * s + p.c2 * s * s + p.c3 * s * s * s
            }
            TissueModel::Ogden(p) => {
                let stretches = principal_stretches(&f);
                let mut w = 0.0;
                let j_bar13 = j_bar.powf(1.0 / 3.0);
                for (mu, alpha) in p.mu.iter().zip(&p.alpha) {
                    let sum: f64 = stretches.iter().map(|l| (l / j_bar13).powf(*alpha)).sum();
                    w += 2.0 * mu / (alpha * alpha) * (sum - 3.0);
                }
                w
            }
            TissueModel::HolzapfelGasserOgden(p) => {
                // The ground-substance and fiber terms carry no (J^-2/3)
                // dilatation factor: the substitution is a no-op, the
                // volumetric penalty is simply dropped.
                crate::hgo::deviatoric_energy(&f, p)
            }
            #[allow(unreachable_patterns)]
            _ => unreachable!("all five variants matched"),
        }
    }

    /// First Piola of [`Self::mean_dilatation_energy`] — the derivative of
    /// the substituted energy with `j_bar` treated as the independent
    /// parameter it is in a mixed or B-bar formulation. Because the
    /// substituted energy drops the pointwise `J`-dependence, this stress
    /// carries **no** `F^{-T}` term — the hydrostatic contribution the
    /// classical deviatoric Piola carries through `dJ/dF` is exactly what
    /// the formulation's constraint stress `p * cof(F)` reinstates through
    /// the pressure field. Analytic for the invariant-based laws
    /// (Neo-Hookean, Mooney-Rivlin, Yeoh); central differences through the
    /// energy for Ogden and HGO, matching the crate's established
    /// analytic/numerical split in [`Self::first_piola`].
    ///
    /// Zero for `J <= 0`, matching every other stress guard here: an
    /// inverted configuration has no valid response.
    pub fn mean_dilatation_first_piola(&self, f: &Mat3, j_bar: f64) -> Mat3 {
        if f.det() <= EPS_F64 {
            return Mat3::ZERO;
        }
        match self {
            TissueModel::NeoHookean(p) => 2.0 * p.c10 * j_bar.powf(-2.0 / 3.0) * *f,
            TissueModel::MooneyRivlin(p) => {
                let f = *f;
                let i1 = invariant_i1(&f);
                let c_mat = f.transpose() * f;
                // d/dF [I2] = 2 (I1 F - F C).
                2.0 * p.c10 * j_bar.powf(-2.0 / 3.0) * f
                    + 2.0 * p.c01 * j_bar.powf(-4.0 / 3.0) * (i1 * f - f * c_mat)
            }
            TissueModel::Yeoh(p) => {
                let f = *f;
                let i1 = invariant_i1(&f);
                let s = j_bar.powf(-2.0 / 3.0) * i1 - 3.0;
                let dw = p.c1 + 2.0 * p.c2 * s + 3.0 * p.c3 * s * s;
                2.0 * j_bar.powf(-2.0 / 3.0) * dw * f
            }
            TissueModel::Ogden(_) | TissueModel::HolzapfelGasserOgden(_) => {
                // Eigenvector-derivative territory, exactly as in
                // `first_piola`: the finite-difference reference through the
                // energy is the crate's verified path for these two.
                const H: f64 = 1.0e-6;
                let f = *f;
                let mut p = Mat3::ZERO;
                for r in 0..3 {
                    for c in 0..3 {
                        let mut fp = f;
                        fp.set(r, c, fp.at(r, c) + H);
                        let mut fm = f;
                        fm.set(r, c, fm.at(r, c) - H);
                        p.set(
                            r,
                            c,
                            (self.mean_dilatation_energy(&fp, j_bar)
                                - self.mean_dilatation_energy(&fm, j_bar))
                                / (2.0 * H),
                        );
                    }
                }
                p
            }
        }
    }

    /// Linearization at `F = I`: the small-strain `(shear, bulk)` moduli
    /// (MPa) of the model — the values a linear-elastic solver needs to
    /// include soft tissue in a model that also carries linear bone
    /// (e.g. via `tpt-med-biomechanics`' mixed-material input).
    ///
    /// Deviatoric: `μ = 2C10` (Neo-Hookean), `2(C10+C01)` (Mooney–Rivlin),
    /// `2c1` (Yeoh), `Σμᵢ` (Ogden), `2c` (HGO ground substance).
    /// Volumetric: every variant shares the `(J−1)²/D1` penalty, so
    /// `K = 2/D1` throughout. The fiber terms of HGO are *not* active at
    /// `F = I` (the invariant `E_f = 0` there), so the HGO linearization
    /// is the ground substance alone — documented, not accidental.
    pub fn linearized_elastic_constants(&self) -> (f64, f64) {
        let shear = match self {
            TissueModel::NeoHookean(p) => 2.0 * p.c10,
            TissueModel::MooneyRivlin(p) => 2.0 * (p.c10 + p.c01),
            TissueModel::Yeoh(p) => 2.0 * p.c1,
            TissueModel::Ogden(p) => p.mu.iter().sum(),
            TissueModel::HolzapfelGasserOgden(p) => 2.0 * p.c,
        };
        (shear, 2.0 / self.volumetric_d1())
    }

    /// Engineering constants `(E, ν)` from the
    /// [`Self::linearized_elastic_constants`] moduli:
    /// `E = 9Kμ/(3K+μ)`, `ν = (3K−2μ)/(2(3K+μ))`. None for non-physical
    /// parameters (`μ ≤ 0`, `K ≤ 0`, or `3K ≤ 2μ`, which would give
    /// `ν ≥ ½`).
    pub fn linearized_engineering_constants(&self) -> Option<(f64, f64)> {
        let (mu, k) = self.linearized_elastic_constants();
        if mu <= 0.0 || k <= 0.0 || 3.0 * k <= 2.0 * mu {
            return None;
        }
        let e = 9.0 * k * mu / (3.0 * k + mu);
        let nu = (3.0 * k - 2.0 * mu) / (2.0 * (3.0 * k + mu));
        Some((e, nu))
    }

    /// The shared `d1` of the `(J−1)²/d1` penalty.
    fn volumetric_d1(&self) -> f64 {
        match self {
            TissueModel::NeoHookean(p) => p.d1,
            TissueModel::MooneyRivlin(p) => p.d1,
            TissueModel::Yeoh(p) => p.d1,
            TissueModel::Ogden(p) => p.d1,
            TissueModel::HolzapfelGasserOgden(p) => p.d1,
        }
    }
}

/// Second-order material tangent: `A[i][j](k,l) = ∂P_ij/∂F_kl`, i.e.
/// `A[i][j]` is the 3×3 matrix of derivatives of stress component `P_ij`
/// with respect to the nine deformation-gradient components.
pub type MaterialTangent = [[Mat3; 3]; 3];

/// An all-zero tangent.
pub(crate) fn zero_tangent() -> MaterialTangent {
    [
        [Mat3::ZERO, Mat3::ZERO, Mat3::ZERO],
        [Mat3::ZERO, Mat3::ZERO, Mat3::ZERO],
        [Mat3::ZERO, Mat3::ZERO, Mat3::ZERO],
    ]
}

/// Analytic deviatoric tangent for the family `P_dev = 2β·q·G` with
/// `β = J^{−2/3}`, `G = F − (Ī1/3)·F^{−T}`: `q = C10` (constant, `q′ = 0`)
/// is Neo-Hookean; `q = q(Ī̄1)` with `q′` its derivative is Yeoh.
///
/// ```text
/// A_dev[ij,kl] = 2β·[(−2/3·F⁻¹_lk·q + q′·(−2/3·Ī1·F⁻¹_lk + 2F_kl))·G_ij
///                    + q·(δ_ik δ_jl − 2/3·F_kl·F⁻¹_ji + Ī1/3·F⁻¹_jk·F⁻¹_li)]
/// ```
fn deviatoric_tangent_two_beta_q_g(f: Mat3, q: f64, q_prime: f64) -> MaterialTangent {
    let j = f.det();
    let mut a = zero_tangent();
    if j <= EPS_F64 {
        return a; // matches the stress-side guard for inverted configurations
    }
    let Some(finv) = f.inverse() else {
        return a;
    };
    let beta = j.powf(-2.0 / 3.0);
    let i1 = invariant_i1(&f);
    let g = f - (i1 / 3.0) * finv.transpose();
    for i in 0..3 {
        for j2 in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let delta_ik = if i == k { 1.0 } else { 0.0 };
                    let delta_jl = if j2 == l { 1.0 } else { 0.0 };
                    // ∂(βq)/∂F_kl = β·[−(2/3)·Fi_lk·q + β·q′·(−(2/3)·Ī1·Fi_lk
                    // + 2F_kl)]: the q′ term carries a second β because ∂s/∂F
                    // has its own (the outer β is factored out below).
                    let d_beta_q = -(2.0 / 3.0) * finv.at(l, k) * q
                        + beta * q_prime * (-(2.0 / 3.0) * i1 * finv.at(l, k) + 2.0 * f.at(k, l));
                    let d_g = delta_ik * delta_jl - (2.0 / 3.0) * f.at(k, l) * finv.at(j2, i)
                        + (i1 / 3.0) * finv.at(j2, k) * finv.at(l, i);
                    a[i][j2].set(k, l, 2.0 * beta * (d_beta_q * g.at(i, j2) + q * d_g));
                }
            }
        }
    }
    a
}

/// Adds the analytic volumetric tangent (shared by every model) in place:
/// `d/dF [2J(J−1)/d1·F^{−T}]`.
fn add_volumetric_tangent(a: &mut MaterialTangent, f: Mat3, d1: f64) {
    let j = f.det();
    if j <= EPS_F64 {
        return; // matches volumetric_first_piola's guard
    }
    let Some(finv) = f.inverse() else {
        return;
    };
    for i in 0..3 {
        for j2 in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let term = (2.0 / d1)
                        * j
                        * ((2.0 * j - 1.0) * finv.at(l, k) * finv.at(j2, i)
                            - (j - 1.0) * finv.at(j2, k) * finv.at(l, i));
                    a[i][j2].set(k, l, a[i][j2].at(k, l) + term);
                }
            }
        }
    }
}

/// `I1 = tr(C)`, `C = Fᵀ F`.
pub fn invariant_i1(f: &Mat3) -> f64 {
    let f = *f;
    (f.transpose() * f).trace()
}

/// `I2 = ½ (I1² − tr(C²))`.
pub fn invariant_i2(f: &Mat3) -> f64 {
    let f = *f;
    let c = f.transpose() * f;
    let i1 = c.trace();
    0.5 * (i1 * i1 - (c * c).trace())
}

/// Principal stretches (square roots of the eigenvalues of `C = Fᵀ F`).
pub fn principal_stretches(f: &Mat3) -> [f64; 3] {
    let f = *f;
    (f.transpose() * f)
        .symmetric_eigenvalues()
        .map(|l| l.max(0.0).sqrt())
}

/// Numeric tolerance used across the tissue crate tests.
pub(crate) const EPS_F64: f64 = 1.0e-12;

/// Reduced-dimensionality wrappers over the full 3×3 `F` interface: a
/// caller working in 2D supplies the four in-plane gradient components and
/// gets the in-plane first Piola back, with the out-of-plane stretch
/// solved per the plane condition — **plane strain** (`F₃₃ = 1`) or
/// **plane stress** (`P₃₃ = 0`, solved by robust bracketing on the scalar
/// unknown `F₃₃`).
#[derive(Debug, Clone)]
pub struct ReducedPlaneModel {
    model: TissueModel,
    condition: PlaneCondition,
}

/// The out-of-plane condition of a reduced model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlaneCondition {
    /// `F₃₃ = 1`: thick section, no out-of-plane deformation.
    PlaneStrain,
    /// `P₃₃ = 0`: thin section, traction-free through-thickness faces.
    PlaneStress,
}

/// The plane problem's solution: in-plane stress, the solved out-of-plane
/// stretch, and the full 3-D first Piola it came from.
#[derive(Debug, Clone, Copy)]
pub struct PlaneSolution {
    /// In-plane components of the first Piola, matching the input layout:
    /// `p[i][j]` = `P_ij` for i, j in {0, 1}.
    pub p_in_plane: [[f64; 2]; 2],
    /// The out-of-plane stretch (1.0 for plane strain, solved for plane
    /// stress).
    pub f33: f64,
    /// The full 3×3 first Piola at the completed gradient.
    pub p_full: Mat3,
}

impl ReducedPlaneModel {
    /// Wraps a full 3-D model under a plane condition.
    pub fn new(model: TissueModel, condition: PlaneCondition) -> Self {
        Self { model, condition }
    }

    /// The wrapped model.
    pub fn model(&self) -> &TissueModel {
        &self.model
    }

    /// The plane condition.
    pub fn condition(&self) -> PlaneCondition {
        self.condition
    }

    /// Solves the plane problem for the in-plane gradient
    /// `[[F11, F12], [F21, F22]]`.
    ///
    /// Plane stress brackets `F₃₃` geometrically (positive, bounded away
    /// from inversion, expanded until `P₃₃` changes sign) and bisects —
    /// `P₃₃` is monotone decreasing in `F₃₃` for every law in this crate,
    /// which is what makes the scalar bracket robust. Returns `None` when
    /// no bracket exists within the expanded bounds (an inadmissible
    /// in-plane gradient, e.g. volumetrically inverted).
    pub fn solve(&self, in_plane: [[f64; 2]; 2]) -> Option<PlaneSolution> {
        let f33 = match self.condition {
            PlaneCondition::PlaneStrain => 1.0,
            PlaneCondition::PlaneStress => self.solve_f33_plane_stress(&in_plane)?,
        };
        let mut full = Mat3::ZERO;
        for (i, row) in in_plane.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                full.set(i, j, v);
            }
        }
        full.set(2, 2, f33);
        let p = self.model.first_piola(&full);
        Some(PlaneSolution {
            p_in_plane: [[p.at(0, 0), p.at(0, 1)], [p.at(1, 0), p.at(1, 1)]],
            f33,
            p_full: p,
        })
    }

    /// Brackets and bisects the root of `P₃₃(F₃₃) = 0`.
    fn solve_f33_plane_stress(&self, in_plane: &[[f64; 2]; 2]) -> Option<f64> {
        let p33 = |f33: f64| {
            let mut f = Mat3::ZERO;
            for (i, row) in in_plane.iter().enumerate() {
                for (j, &v) in row.iter().enumerate() {
                    f.set(i, j, v);
                }
            }
            f.set(2, 2, f33);
            self.model.first_piola(&f).at(2, 2)
        };
        // Expand outward from a geometric bracket until P₃₃ changes sign.
        let mut lo = 0.05f64;
        let mut hi = 2.0f64;
        let mut flo = p33(lo);
        let mut fhi = p33(hi);
        let mut guard = 0;
        while flo * fhi > 0.0 {
            guard += 1;
            if guard > 60 {
                return None;
            }
            lo *= 0.5;
            hi *= 2.0;
            flo = p33(lo);
            fhi = p33(hi);
        }
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            let fm = p33(mid);
            if fm == 0.0 {
                return Some(mid);
            }
            if fm * flo < 0.0 {
                hi = mid;
            } else {
                lo = mid;
                flo = fm;
            }
        }
        Some(0.5 * (lo + hi))
    }
}

/// A named soft-tissue material: model + physical metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct SoftTissueMaterial {
    /// Constitutive model.
    pub model: TissueModel,
    /// Mass density (g/cm³), ~1.06 for most soft tissue.
    pub density: f64,
    /// True when the model is used under strictly incompressible kinematics.
    pub is_incompressible: bool,
}

impl SoftTissueMaterial {
    /// Strain energy (see [`TissueModel::strain_energy`]).
    pub fn strain_energy(&self, f: &Mat3) -> f64 {
        self.model.strain_energy(f)
    }

    /// First Piola–Kirchhoff stress (see [`TissueModel::first_piola`]).
    pub fn first_piola(&self, f: &Mat3) -> Mat3 {
        self.model.first_piola(f)
    }

    /// Material tangent (see [`TissueModel::material_tangent`]).
    pub fn material_tangent(&self, f: &Mat3) -> MaterialTangent {
        self.model.material_tangent(f)
    }
}
