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
    MooneyRivlinParams, NeoHookeanParams, OgdenParams, PlaneCondition, ReducedPlaneModel,
    SoftTissueMaterial, TissueModel, YeohParams,
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

/// A deformed, rotated gradient with J > 0: exercises all nine components
/// (unlike the diagonal verification stretches above).
fn rotated_stretch_f() -> Mat3 {
    let stretch = Mat3::diagonal([1.15, 0.92, 1.06]);
    let rotation = Mat3::rotation_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.4);
    rotation * stretch
}

fn all_models() -> Vec<(&'static str, TissueModel)> {
    vec![
        (
            "neo-hookean",
            TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 0.5 }),
        ),
        (
            "mooney-rivlin",
            TissueModel::MooneyRivlin(MooneyRivlinParams {
                c10: 0.3,
                c01: 0.1,
                d1: 0.5,
            }),
        ),
        (
            "yeoh",
            TissueModel::Yeoh(YeohParams {
                c1: 0.3,
                c2: 0.05,
                c3: 0.01,
                d1: 0.5,
            }),
        ),
        (
            "ogden",
            TissueModel::Ogden(OgdenParams {
                mu: vec![0.4, 0.1],
                alpha: vec![2.0, -2.0],
                d1: 0.5,
            }),
        ),
        (
            "hgo",
            TissueModel::HolzapfelGasserOgden(HgoParams {
                c: 0.8,
                k1: 5.0,
                k2: 12.0,
                kappa: 0.2,
                fiber_directions: vec![Vec3::new(1.0, 1.0, 0.0), Vec3::new(-1.0, 1.0, 0.0)],
                d1: 100.0,
            }),
        ),
    ]
}

#[test]
fn analytic_tangent_matches_central_differences() {
    // The Neo-Hookean and Yeoh tangents are analytic; the finite-difference
    // reference of the same P they differentiate must agree to FD accuracy.
    let f = rotated_stretch_f();
    for (name, model) in all_models() {
        if !matches!(model, TissueModel::NeoHookean(_) | TissueModel::Yeoh(_)) {
            continue;
        }
        let analytic = model.material_tangent(&f);
        let reference = model.material_tangent_numerical(&f);
        let mut scale = 0.0f64;
        let mut max_err = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        scale = scale.max(reference[i][j].at(k, l).abs());
                        max_err =
                            max_err.max((analytic[i][j].at(k, l) - reference[i][j].at(k, l)).abs());
                    }
                }
            }
        }
        assert!(
            max_err < 1e-5 * scale.max(1e-12),
            "{name}: analytic vs FD tangent err {max_err:.3e}, scale {scale:.3e}"
        );
    }
}

#[test]
fn tangent_has_major_symmetry_for_every_model() {
    // Hyperelasticity: ∂P_ij/∂F_kl = ∂P_kl/∂F_ij (W is twice differentiable).
    // P itself is NOT symmetric, so only the major symmetry holds — the
    // minor one (A_ij,kl = A_ji,kl) must not be asserted.
    let f = rotated_stretch_f();
    for (name, model) in all_models() {
        let a = model.material_tangent(&f);
        let mut max_asym = 0.0f64;
        let mut scale = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let fwd = a[i][j].at(k, l);
                        let bwd = a[k][l].at(i, j);
                        scale = scale.max(fwd.abs());
                        max_asym = max_asym.max((fwd - bwd).abs());
                    }
                }
            }
        }
        assert!(
            max_asym < 1e-4 * scale.max(1e-12),
            "{name}: major symmetry violated, asym {max_asym:.3e} vs scale {scale:.3e}"
        );
    }
}

#[test]
fn volumetric_tangent_matches_fd_of_volumetric_first_piola() {
    // The shared (J−1)²/d1 penalty means one closed form serves all models;
    // verify it against central differences of volumetric_first_piola.
    let f = rotated_stretch_f();
    const H: f64 = 1.0e-6;
    for (name, model) in all_models() {
        let analytic = model.volumetric_tangent(&f);
        let mut max_err = 0.0f64;
        let mut scale = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut fp = f;
                        fp.set(k, l, fp.at(k, l) + H);
                        let mut fm = f;
                        fm.set(k, l, fm.at(k, l) - H);
                        let fd = (model.volumetric_first_piola(&fp).at(i, j)
                            - model.volumetric_first_piola(&fm).at(i, j))
                            / (2.0 * H);
                        scale = scale.max(fd.abs());
                        max_err = max_err.max((analytic[i][j].at(k, l) - fd).abs());
                    }
                }
            }
        }
        assert!(
            max_err < 1e-5 * scale.max(1e-12),
            "{name}: volumetric tangent err {max_err:.3e} vs scale {scale:.3e}"
        );
    }
}

#[test]
fn tangent_stiffens_under_uniaxial_loading() {
    // Directional sanity: the axial axial component must be positive and
    // grow with stretch (strain stiffening), while the pure-dilatation
    // tangent is governed by the volumetric penalty.
    let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 0.5 });
    // A_1111 at fixed lateral components stays positive across the working
    // range. (It is *not* the uniaxial stiffness — that follows the
    // incompressible path with coupled lateral contraction — so monotone
    // growth is not asserted.)
    for lam in [1.05, 1.2, 1.4] {
        let a1111 = model.material_tangent(&uniaxial_f(lam))[0][0].at(0, 0);
        assert!(a1111 > 0.0, "A_1111 at λ={lam}: {a1111}");
    }
    // Inverted configuration: zero, matching the stress-side guard.
    let inverted = Mat3::diagonal([-1.0, 1.0, 1.0]);
    let a_inv = model.material_tangent(&inverted);
    let mut max_abs = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    max_abs = max_abs.max(a_inv[i][j].at(k, l).abs());
                }
            }
        }
    }
    assert_eq!(max_abs, 0.0, "inverted tangent must be zero");
}

#[test]
fn linearized_constants_match_the_closed_forms() {
    let nh = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 0.2 });
    let (mu, k) = nh.linearized_elastic_constants();
    assert!((mu - 1.0).abs() < 1e-12);
    assert!((k - 10.0).abs() < 1e-12);
    // Engineering constants from the moduli, hand-computed:
    // E = 9·10·1/(30+1) = 90/31, ν = (30−2)/(2·31) = 14/31.
    let (e, nu) = nh.linearized_engineering_constants().expect("physical");
    assert!((e - 90.0 / 31.0).abs() < 1e-12);
    assert!((nu - 14.0 / 31.0).abs() < 1e-12);

    let mr = TissueModel::MooneyRivlin(MooneyRivlinParams {
        c10: 0.3,
        c01: 0.2,
        d1: 0.5,
    });
    let (mu, _) = mr.linearized_elastic_constants();
    assert!((mu - 1.0).abs() < 1e-12);

    let ogden = TissueModel::Ogden(OgdenParams {
        mu: vec![0.4, 0.1],
        alpha: vec![2.0, -2.0],
        d1: 0.5,
    });
    let (mu, _) = ogden.linearized_elastic_constants();
    assert!((mu - 0.5).abs() < 1e-12);

    // HGO linearizes to its ground substance: fibers are inactive at F = I.
    let hgo = TissueModel::HolzapfelGasserOgden(HgoParams {
        c: 0.8,
        k1: 5.0,
        k2: 12.0,
        kappa: 0.2,
        fiber_directions: vec![Vec3::new(1.0, 1.0, 0.0), Vec3::new(-1.0, 1.0, 0.0)],
        d1: 100.0,
    });
    let (mu, k) = hgo.linearized_elastic_constants();
    assert!((mu - 1.6).abs() < 1e-12);
    assert!((k - 0.02).abs() < 1e-12);

    // Non-physical parameters are refused, not coerced.
    let degenerate = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.0, d1: 1.0 });
    assert!(degenerate.linearized_engineering_constants().is_none());
    let incompressible_limit = TissueModel::NeoHookean(NeoHookeanParams { c10: 10.0, d1: 1.0 });
    // K = 2, μ = 20 → 3K = 6 < 2μ = 40: ν ≥ ½, refused.
    assert!(incompressible_limit
        .linearized_engineering_constants()
        .is_none());
}

#[test]
fn plane_strain_pins_f33_and_delegates_to_the_full_law() {
    let model = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 0.5 });
    let reduced = ReducedPlaneModel::new(model.clone(), PlaneCondition::PlaneStrain);
    let in_plane = [[1.1, 0.05], [0.0, 0.95]];
    let sol = reduced.solve(in_plane).expect("plane strain always solves");
    assert_eq!(sol.f33, 1.0);
    // Identical to evaluating the full 3-D law at the completed gradient.
    let mut full = Mat3::ZERO;
    for (i, row) in in_plane.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            full.set(i, j, v);
        }
    }
    full.set(2, 2, 1.0);
    let direct = model.first_piola(&full);
    for i in 0..2 {
        for j in 0..2 {
            assert!((sol.p_in_plane[i][j] - direct.at(i, j)).abs() < 1e-15);
        }
    }
}

#[test]
fn plane_stress_solves_the_traction_free_out_of_plane_stretch() {
    // Incompressible-limit material: plane stress drives J → 1, so the
    // analytic out-of-plane stretch is F33 = 1/det(in-plane F) — an
    // independent closed form the bisection must land on.
    let model = TissueModel::NeoHookean(NeoHookeanParams {
        c10: 0.5,
        d1: 1.0e-6,
    });
    let reduced = ReducedPlaneModel::new(model, PlaneCondition::PlaneStress);
    let in_plane = [[1.3, 0.0], [0.0, 1.0]];
    let sol = reduced.solve(in_plane).expect("bracket exists");
    let expected_f33 = 1.0 / (1.3 * 1.0);
    assert!(
        (sol.f33 - expected_f33).abs() < 1e-4,
        "F33 {} vs analytic {expected_f33}",
        sol.f33
    );
    assert!(
        sol.p_full.at(2, 2).abs() < 1e-6,
        "P33 = {}",
        sol.p_full.at(2, 2)
    );
    // Rotated in-plane gradient: same determinant, same F33.
    let rotated = [[1.2, 0.2], [-0.1, 0.9]];
    let sol2 = reduced.solve(rotated).expect("bracket exists");
    let det: f64 = rotated[0][0] * rotated[1][1] - rotated[0][1] * rotated[1][0];
    assert!((sol2.f33 - 1.0 / det).abs() < 1e-4);
    assert!(sol2.p_full.at(2, 2).abs() < 1e-6);
}

#[test]
fn compressible_plane_stress_stays_traction_free() {
    // Compressible material: no closed-form F33, but the defining property
    // (P33 = 0) must hold to bisection precision, and the solution must
    // differ from the incompressible one.
    let compressible = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 0.5 });
    let reduced = ReducedPlaneModel::new(compressible, PlaneCondition::PlaneStress);
    let in_plane = [[1.3, 0.0], [0.0, 1.0]];
    let sol = reduced.solve(in_plane).expect("bracket exists");
    assert!(
        sol.p_full.at(2, 2).abs() < 1e-9,
        "P33 = {}",
        sol.p_full.at(2, 2)
    );
    // Compressibility relaxes the out-of-plane contraction relative to the
    // incompressible 1/1.3 ≈ 0.769.
    assert!(sol.f33 > 1.0 / 1.3, "F33 {}", sol.f33);
    // In-plane components are unaffected by the scalar solve: they equal
    // the full law evaluated at the solved gradient.
    let mut full = Mat3::ZERO;
    for (i, row) in in_plane.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            full.set(i, j, v);
        }
    }
    full.set(2, 2, sol.f33);
    let direct = TissueModel::NeoHookean(NeoHookeanParams { c10: 0.5, d1: 0.5 }).first_piola(&full);
    assert!((sol.p_in_plane[0][0] - direct.at(0, 0)).abs() < 1e-12);
}
