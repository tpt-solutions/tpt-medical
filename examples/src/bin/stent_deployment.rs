//! Phase 5 milestone: stent deployment simulation with artery contact.
//!
//! Crimps a Nitinol ring stent onto a delivery system, expands it into a
//! compliant artery model, and reports radial force, contact pressure,
//! recoil, and dogboning across a range of balloon pressures (the
//! ASTM F2394-style bench metrics).
//!
//! ```console
//! cargo run -p tpt-med-examples --bin stent-deployment
//! ```

use tpt_med_stents::{simulate_deployment, NitinolParams, StentModel, SuperelasticState};
use tpt_med_units::Pressure;

fn main() {
    println!("=== tpt-medical stent deployment (Phase 5 milestone) ===");

    let nitinol = NitinolParams::default();
    let stent = StentModel {
        expanded_diameter: 6.0, // mm nominal
        crimped_diameter: 1.8,  // mm on the balloon
        n_crowns: 12,
        crown_stiffness: 0.5, // N/mm per crown
    };

    // Crimping check: the DIAMETER change is large (~70%), but the crown
    // bending strain scales with strut thickness over crown radius:
    // ε ≈ 2·(t/D)·(ΔD/D). With t = 0.16 mm on a 6 mm ring this lands
    // ~3.7% — inside the superelastic plateau. Sweep the material through
    // that cycle.
    let crimp_ratio = (stent.expanded_diameter - stent.crimped_diameter) / stent.expanded_diameter;
    let strut_strain = 2.0 * (0.16 / stent.expanded_diameter) * crimp_ratio;
    let mut state = SuperelasticState::new();
    let p = &nitinol;
    let mut peak_stress = 0.0f64;
    let n = 200;
    // Crown apexes locally exceed the plateau; sweep to 8% so the cycle
    // passes martensite finish before unloading.
    let cycle_strain = 0.08f64;
    for i in 0..=n {
        let e = cycle_strain * i as f64 / n as f64;
        peak_stress = peak_stress.max(state.strain_to_stress(e, p));
    }
    let mut stress_after_release = 0.0;
    for i in (0..=n).rev() {
        let e = cycle_strain * i as f64 / n as f64;
        stress_after_release = state.strain_to_stress(e, p);
    }
    println!(
        "crimp: ΔD = {:.1} mm ({:.0}% diameter change → ~{:.1}% strut strain),          peak strut stress {:.0} MPa, stress after release {:.0} MPa (superelastic recovery)",
        stent.expanded_diameter - stent.crimped_diameter,
        100.0 * crimp_ratio,
        100.0 * strut_strain,
        peak_stress,
        stress_after_release
    );

    // Deployment across clinically relevant balloon pressures; the vessel
    // pressure–diameter law: compliant artery, D = 4.6 + 6.0·p (mm at MPa)
    // with a stiff-plaque floor at 4.6 mm.
    let vessel = |p_mpa: f64| 4.6 + 6.0 * p_mpa;
    println!(
        "deployment vs balloon pressure (artery ≈ {:?} mm lumen range):",
        4.6..5.2
    );
    println!(
        "{:>8} {:>10} {:>10} {:>12} {:>10} {:>10}",
        "p[MPa]", "D_eq[mm]", "F_rad[N]", "P_contact[Pa]", "recoil[%]", "dogbone[%]"
    );
    for kpa in [60.0, 80.0, 100.0, 121.0, 141.0] {
        let balloon = Pressure::from_kpa(kpa);
        let r = simulate_deployment(&stent, &nitinol, vessel, balloon);
        println!(
            "{:>8.2} {:>10.3} {:>10.2} {:>12.0} {:>9.1}% {:>9.1}%",
            kpa / 10.0,
            r.diameter,
            r.radial_force,
            r.contact_pressure * 1.0e6, // MPa → Pa for readability
            100.0 * r.recoil,
            100.0 * r.dogboning
        );
    }

    // Chronic outward force (COF) and acute recoil summary at nominal
    // deployment pressure (~8 atm ≈ 0.81 MPa... contact law above uses lumen
    // at the given pressure).
    let nominal = simulate_deployment(&stent, &nitinol, vessel, Pressure::from_mpa(0.10));
    println!("--- summary at nominal deployment ---");
    println!(
        "  equilibrium diameter {:.2} mm (nominal {:.1} mm), recoil {:.1}%",
        nominal.diameter,
        stent.expanded_diameter,
        100.0 * nominal.recoil
    );
    println!("  chronic outward force: {:.2} N", nominal.radial_force);
    println!(
        "  vessel wall contact pressure: {:.1} kPa",
        nominal.contact_pressure * 1000.0
    );
}
