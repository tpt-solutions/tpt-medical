//! Verification suite: model stresses vs textbook closed-form solutions.
//!
//! Nearly-incompressible penalty formulations carry a model-internal
//! hydrostatic pressure at `J = 1` (the derivative of `J^{-2/3}`), so the
//! pressure-dependent components are not unique. The pressure-**independent**
//! quantity is the deviatoric Cauchy stress `s = σ − (tr σ/3) I`; closed-form
//! incompressible solutions define it exactly. All uniaxial verifications
//! below therefore compare `s`. Volumetric response is verified separately
//! under pure dilatation, and every analytic `P` is additionally checked
//! against a central-difference reference.

use crate::hgo::HgoParams;
use crate::models::{
    MooneyRivlinParams, NeoHookeanParams, OgdenParams, SoftTissueMaterial, TissueModel, YeohParams,
};
use tpt_med_geometry::{Mat3, Vec3};

/// Uniaxial stretch with incompressible lateral contraction: `J = 1`.
fn uniaxial_f(lam: f64) -> Mat3 {
    let s = 1.0 / lam.sqrt();
    Mat3::from_rows(
        Vec3::new(lam, 0.0, 0.0),
        Vec3::new(0.0, s, 0.0),
        Vec3::new(0.0, 0.0, s),
    )
}

/// Left Cauchy–Green tensor `B = F Fᵀ`.
fn left_cauchy_green(f: &Mat3) -> Mat3 {
    *f * f.transpose()
}

/// Deviatoric Cauchy stress `s = σ − (tr σ/3) I` with `σ = P Fᵀ / J`.
fn deviatoric_cauchy(model: &TissueModel, f: &Mat3) -> Mat3 {
    let j = f.det();
    let sigma = model.first_piola(f) * f.transpose() * (1.0 / j);
    let p = sigma.trace() / 3.0;
    sigma - Mat3::IDENTITY * p
}

/// Largest absolute entry difference between two matrices.
fn max_diff(a: &Mat3, b: &Mat3) -> f64 {
    let mut mx = 0.0f64;
    for r in 0..3 {
        for c in 0..3 {
            mx = mx.max((a.at(r, c) - b.at(r, c)).abs());
        }
    }
    mx
}

#[test]
fn neo_hookean_deviatoric_stress_matches_closed_form() {
    // Incompressible NH: σ = −pI + 2 C10 B  ⇒  s = dev(2 C10 B).
    let c10 = 0.49;
    let model = TissueModel::NeoHookean(NeoHookeanParams { c10, d1: 100.0 });
    for lam in [1.1, 1.3, 1.6] {
        let f = uniaxial_f(lam);
        let b = left_cauchy_green(&f);
        let tr_b = b.trace();
        let mut s_ref = Mat3::ZERO;
        for i in 0..3 {
            s_ref.set(i, i, 2.0 * c10 * (b.at(i, i) - tr_b / 3.0));
        }
        let s_code = deviatoric_cauchy(&model, &f);
        assert!(max_diff(&s_code, &s_ref) < 1e-12, "λ={lam}: {s_code:?}");
    }
}

#[test]
fn mooney_rivlin_deviatoric_stress_matches_closed_form() {
    // Incompressible MR: σ = −pI + 2 C10 B − 2 C01 B⁻¹
    // ⇒ s = dev(2 C10 B − 2 C01 B⁻¹).
    let (c10, c01) = (1.3, 0.7);
    let model = TissueModel::MooneyRivlin(MooneyRivlinParams {
        c10,
        c01,
        d1: 100.0,
    });
    for lam in [1.1, 1.35, 1.7] {
        let f = uniaxial_f(lam);
        let b = left_cauchy_green(&f);
        let binv = b.inverse().expect("B invertible on J=1");
        let mut s_ref = Mat3::ZERO;
        for i in 0..3 {
            s_ref.set(
                i,
                i,
                2.0 * c10 * (b.at(i, i) - b.trace() / 3.0)
                    - 2.0 * c01 * (binv.at(i, i) - binv.trace() / 3.0),
            );
        }
        let s_code = deviatoric_cauchy(&model, &f);
        assert!(max_diff(&s_code, &s_ref) < 1e-12, "λ={lam}: {s_code:?}");
    }
}

#[test]
fn yeoh_deviatoric_stress_matches_closed_form() {
    // Yeoh W = Σ cᵢ(Ī1−3)ⁱ: s = dev(2 dW/dI1 · B).
    let p = YeohParams {
        c1: 0.3,
        c2: -0.05,
        c3: 0.02,
        d1: 100.0,
    };
    let model = TissueModel::Yeoh(p);
    for lam in [1.05, 1.3, 1.8] {
        let f = uniaxial_f(lam);
        let b = left_cauchy_green(&f);
        let s_bar = b.trace() - 3.0; // J = 1 ⇒ Ī1 = I1
        let dw = p.c1 + 2.0 * p.c2 * s_bar + 3.0 * p.c3 * s_bar * s_bar;
        let mut s_ref = Mat3::ZERO;
        for i in 0..3 {
            s_ref.set(i, i, 2.0 * dw * (b.at(i, i) - b.trace() / 3.0));
        }
        let s_code = deviatoric_cauchy(&model, &f);
        assert!(max_diff(&s_code, &s_ref) < 1e-12, "λ={lam}: {s_code:?}");
    }
}

#[test]
fn ogden_alpha2_matches_neo_hookean_deviation() {
    // One-term Ogden with α = 2 equals NH with C10 = μ/2.
    let mu = 0.6;
    let ogden = TissueModel::Ogden(OgdenParams {
        mu: vec![mu],
        alpha: vec![2.0],
        d1: 100.0,
    });
    let nh = TissueModel::NeoHookean(NeoHookeanParams {
        c10: mu / 2.0,
        d1: 100.0,
    });
    for lam in [1.1, 1.4] {
        let f = uniaxial_f(lam);
        assert!(
            (ogden.strain_energy(&f) - nh.strain_energy(&f)).abs() < 1e-9,
            "λ={lam}"
        );
        assert!(max_diff(&deviatoric_cauchy(&ogden, &f), &deviatoric_cauchy(&nh, &f)) < 1e-9);
    }
}

#[test]
fn all_analytic_stresses_match_finite_difference() {
    // Every model with an analytic P is locked to its central-difference
    // reference (this catches sign/invariant errors independently of any
    // closed form).
    let models: Vec<TissueModel> = vec![
        TissueModel::NeoHookean(NeoHookeanParams {
            c10: 0.49,
            d1: 10.0,
        }),
        TissueModel::MooneyRivlin(MooneyRivlinParams {
            c10: 1.3,
            c01: 0.7,
            d1: 10.0,
        }),
        TissueModel::Yeoh(YeohParams {
            c1: 0.3,
            c2: -0.05,
            c3: 0.02,
            d1: 10.0,
        }),
        TissueModel::HolzapfelGasserOgden(HgoParams {
            c: 0.8,
            k1: 5.0,
            k2: 12.0,
            kappa: 0.2,
            fiber_directions: vec![Vec3::new(1.0, 1.0, 0.0), Vec3::new(-1.0, 1.0, 0.0)],
            d1: 100.0,
        }),
    ];
    // Simple shear + stretch states.
    let mut shear = Mat3::IDENTITY;
    shear.set(0, 1, 0.3);
    shear.set(1, 0, 0.1);
    for model in &models {
        for f in [uniaxial_f(1.2), shear] {
            let a = model.first_piola(&f);
            let n = model.first_piola_numerical(&f);
            let scale = [a.at(0, 0).abs(), a.at(1, 1).abs(), a.at(2, 2).abs()]
                .iter()
                .fold(1e-3f64, |m, &v| m.max(v));
            assert!(max_diff(&a, &n) < 1e-4 * scale, "{model:?}: {a:?} vs {n:?}");
        }
    }
}

#[test]
fn volumetric_term_produces_pressure_on_dilatation() {
    // Pure dilatation must produce purely hydrostatic P with magnitude
    // driven by (J−1)/D1 — the verified volumetric branch.
    let j = 1.2f64;
    let f = Mat3::IDENTITY * j.powf(1.0 / 3.0);
    let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 10.0 });
    let p = model.first_piola(&f);
    for i in 0..3 {
        for jj in 0..3 {
            if i != jj {
                assert!(p.at(i, jj).abs() < 1e-12);
            }
        }
    }
    assert!((p.at(0, 0) - p.at(1, 1)).abs() < 1e-12);
    assert!((p.at(0, 0) - p.at(2, 2)).abs() < 1e-12);
    // With I1 = 3 J^{2/3} the isochoric branch vanishes identically and
    // P11 collapses to 2 J (J−1) / D1 · J^{-1/3}.
    let expected = 2.0 * j * (j - 1.0) / 10.0 * j.powf(-1.0 / 3.0);
    assert!(
        (p.at(0, 0) - expected).abs() < 1e-9,
        "{} vs {expected}",
        p.at(0, 0)
    );
}

#[test]
fn soft_tissue_material_delegates() {
    let m = SoftTissueMaterial {
        model: TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 10.0 }),
        density: 1.06,
        is_incompressible: false,
    };
    let f = uniaxial_f(1.1);
    assert_eq!(m.strain_energy(&f), m.model.strain_energy(&f));
    assert_eq!(m.first_piola(&f), m.model.first_piola(&f));
}
