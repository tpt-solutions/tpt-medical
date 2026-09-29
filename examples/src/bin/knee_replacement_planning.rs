//! Phase 6 milestone: virtual total knee replacement planning.
//!
//! Builds a synthetic tibial-femoral voxel scene, performs a joint-line
//! osteotomy, moves the resected fragment to restore alignment, and sizes
//! femoral/tibial components from landmarks.
//!
//! ```console
//! cargo run -p tpt-med-examples --bin knee-replacement-planning
//! ```

use tpt_med_geometry::{Plane, Vec3};
use tpt_med_implant_sizing::{size_tka, KneeLandmarks, SizeChart, SizeEntry};
use tpt_med_surgical_planning::{
    DiscardedSide, FragmentTransform, OsteotomyCut, VirtualSurgery, VoxelModel,
};

fn main() {
    println!("=== tpt-medical virtual TKA planning (Phase 6 milestone) ===");

    // Synthetic distal femur: 40×40×30 grid, 1 mm voxels; distal condyles
    // modelled as a blunt block for the cut demo.
    let mut femur = VoxelModel {
        dims: (40, 40, 30),
        spacing: (1.0, 1.0, 1.0),
        origin: Vec3::new(-20.0, -20.0, -30.0),
        values: vec![0.0; 40 * 40 * 30],
    };
    for z in 0..30 {
        for y in 0..40 {
            for x in 0..40 {
                let xf = x as f64;
                let yf = y as f64;
                let condyle = z < 8 && (xf - 20.0).abs() < 14.0 && (yf - 20.0).abs() < 16.0;
                let shaft = z >= 8 && (xf - 20.0).abs() < 9.0 && (yf - 20.0).abs() < 9.0;
                if condyle || shaft {
                    let idx = femur.index(x, y, z).unwrap();
                    femur.values[idx] = 700.0; // cortical HU
                }
            }
        }
    }
    println!(
        "pre-op distal femur: {} voxels of bone",
        femur.count_above(200.0)
    );

    // Landmarks consistent with the synthetic scene (patient coords, mm).
    let landmarks = KneeLandmarks {
        medial_epicondyle: Vec3::new(-15.0, 20.0, -22.0),
        lateral_epicondyle: Vec3::new(15.0, 20.0, -22.0),
        trochlea_point: Vec3::new(0.0, 20.0, -30.0),
        posterior_condyle: Vec3::new(0.0, 20.0, -22.0),
        tibial_medial: Vec3::new(-14.0, 18.0, -52.0),
        tibial_lateral: Vec3::new(14.0, 18.0, -52.0),
        tibial_center: Vec3::new(0.0, 18.0, -52.0),
        tibial_tubercle: Vec3::new(3.0, 18.0, -50.0),
    };

    // Plan: 9 mm resection plane above the distal condyle level, then a
    // 2° correction rotation about the TEA and +2 mm distal release.
    let resection_z = -21.0;
    let mut plan = VirtualSurgery::new(femur.clone());
    plan.cut(OsteotomyCut {
        plane: Plane::from_point_normal(Vec3::new(0.0, 0.0, resection_z), Vec3::Z).expect("plane"),
        fragment_name: "distal_resection".into(),
        keep_positive: true,
        kerf_width: 0.0,
        discarded: DiscardedSide::Resect,
    })
    .expect("valid plan");
    plan.move_fragment(FragmentTransform {
        rotation_axis: Vec3::new(-1.0, 0.0, 0.0), // TEA (medial→lateral = +x; rotate about it)
        rotation_angle: 2.0f64.to_radians(),
        pivot: Vec3::new(0.0, 20.0, -22.0),
        translation: Vec3::new(0.0, 0.0, -2.0),
    });
    let (operated, log) = plan.execute();
    println!("plan steps executed: {log:?}");
    println!(
        "post-op bone voxels: {} (resected {})",
        operated.count_above(200.0),
        femur.count_above(200.0) - operated.count_above(200.0)
    );

    // Component sizing from the (pre-operative) landmarks.
    let femoral_chart = SizeChart {
        family: "synthetic-tka-femoral".into(),
        entries: (1..=8)
            .map(|i| SizeEntry {
                label: i,
                nominal: 50.0 + 4.0 * i as f64,
            })
            .collect(),
    };
    let tibial_chart = SizeChart {
        family: "synthetic-tka-tibial".into(),
        entries: (1..=6)
            .map(|i| SizeEntry {
                label: i,
                nominal: 45.0 + 5.0 * i as f64,
            })
            .collect(),
    };
    match size_tka(&landmarks, &femoral_chart, &tibial_chart) {
        Some(s) => {
            println!("--- component sizing ---");
            println!(
                "  measurements: TEA {:.1} mm, AP {:.1} mm, plateau {:.1} mm",
                s.tea_width_mm, s.ap_depth_mm, s.plateau_width_mm
            );
            println!(
                "  femoral component: size {} | tibial component: size {}",
                s.femoral_size, s.tibial_size
            );
            println!(
                "  alignment: femorotibial {:.1}°, tibial slope {:.1}°",
                s.femorotibial_angle_deg, s.tibial_slope_deg
            );
        }
        None => println!("  sizing failed: chart coverage"),
    }
}
