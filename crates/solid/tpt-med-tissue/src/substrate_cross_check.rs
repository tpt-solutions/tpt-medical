//! Cross-checks the in-house closed-form uniaxial Neo-Hookean stress against
//! `tpt-fem-hyperelastic`'s independently-implemented 1-D bar Newton solve
//! (`rfcs/0009-nonlinear-fem-substrate-adapter.md`).
//!
//! This is the narrow, safe first slice of substrate integration that RFC
//! calls out explicitly: no new production API, no 3D assembly (which the
//! substrate does not yet provide at 0.1.0 — see the RFC's "actual gap"),
//! just two independently-coded implementations of the same incompressible
//! Neo-Hookean uniaxial-tension physics agreeing on the same answer.
//!
//! Deliberately duplicates (rather than shares) the `uniaxial_f`/deviatoric-
//! stress helpers already private to `tests.rs` — this module is gated
//! behind an off-by-default feature and should not create a dependency
//! edge, in either direction, on the always-on verification suite.

use crate::models::{NeoHookeanParams, TissueModel};
use tpt_fem_mesh::{CellType, MeshBuilder};
use tpt_med_geometry::{Mat3, Vec3};

/// Uniaxial stretch with incompressible lateral contraction: `J = 1`. Same
/// construction as `tests.rs`'s helper of the same name.
fn uniaxial_f(lam: f64) -> Mat3 {
    let s = 1.0 / lam.sqrt();
    Mat3::from_rows(
        Vec3::new(lam, 0.0, 0.0),
        Vec3::new(0.0, s, 0.0),
        Vec3::new(0.0, 0.0, s),
    )
}

/// `sigma11 - sigma22` of the Cauchy stress `sigma = P Fᵀ / J`. For an
/// isotropic incompressible material under uniaxial tension with
/// traction-free lateral sides (`sigma22 = sigma33 = 0`, which
/// `uniaxial_f`'s symmetric lateral contraction satisfies exactly), this
/// difference equals `sigma11` itself — the physically correct axial true
/// stress, independent of the model's own internal (and, for a penalty
/// formulation, otherwise arbitrary) pressure term. Same reasoning
/// `tests.rs`'s `deviatoric_cauchy` helper relies on, carried one step
/// further into a single scalar so it can be compared against the
/// substrate's nominal-stress formula below.
fn sigma11_minus_sigma22(model: &TissueModel, f: &Mat3) -> f64 {
    let j = f.det();
    let sigma = model.first_piola(f) * f.transpose() * (1.0 / j);
    sigma.at(0, 0) - sigma.at(1, 1)
}

/// A straight 1-D bar mesh of `n_elements` equal `Line2` segments spanning
/// `[0, length]`, matching `tpt_fem_hyperelastic::solve_hyperelastic_bar`'s
/// expected input shape.
fn bar_mesh(n_elements: usize, length: f64) -> tpt_fem_mesh::Mesh {
    let mut b = MeshBuilder::new();
    let nodes: Vec<_> = (0..=n_elements)
        .map(|i| b.add_node(vec![length * i as f64 / n_elements as f64]))
        .collect();
    for i in 0..n_elements {
        b.add_element(CellType::Line, vec![nodes[i], nodes[i + 1]]);
    }
    b.build()
}

#[test]
fn neo_hookean_uniaxial_matches_substrate_bar_solve() {
    let c10 = 0.49;
    // Standard incompressible-NH convention shared by both formulations:
    // nominal uniaxial stress `P11 = mu (lambda - lambda^-2)` with
    // `mu = 2 * c10` -- see rfcs/0009's "actual gap" item 2 for the
    // derivation of why the in-house penalty model's raw `first_piola`
    // does *not* equal this directly, and why `sigma11_minus_sigma22`
    // (not `first_piola` alone) is the quantity that does.
    let mu = 2.0 * c10;
    let length = 10.0;
    let mesh = bar_mesh(20, length);

    for &lam in &[1.1, 1.3, 1.6] {
        // Independent numerical solve: the substrate's own Newton
        // iteration on a 1-D bar mesh, prescribing the end displacement
        // that achieves stretch `lam`.
        let disp = tpt_fem_hyperelastic::solve_hyperelastic_bar(&mesh, 1.0, mu, lam)
            .expect("substrate Newton solve converges");
        let n = disp.len();
        let achieved_stretch = 1.0 + (disp[n - 1] - disp[0]) / length;
        assert!(
            (achieved_stretch - lam).abs() < 1e-9,
            "substrate solve did not achieve the prescribed stretch"
        );
        // The substrate's own nominal-stress formula, evaluated at its own
        // solved stretch (should reproduce the closed form it was built
        // from — this is the substrate's internal self-consistency, not
        // yet the cross-check).
        let substrate_p11 = tpt_fem_hyperelastic::neo_hookean_1d(achieved_stretch, mu);

        // In-house closed form, via the deviatoric-stress route.
        let model = TissueModel::NeoHookean(NeoHookeanParams { c10, d1: 1.0e6 });
        let f = uniaxial_f(lam);
        let sigma11 = sigma11_minus_sigma22(&model, &f);
        let in_house_p11 = sigma11 / lam;

        assert!(
            (in_house_p11 - substrate_p11).abs() / substrate_p11.abs() < 1e-9,
            "lam={lam}: in-house P11={in_house_p11}, substrate P11={substrate_p11}"
        );
    }
}
