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
}
