//! Automated implant sizing from anatomical landmarks.
//!
//! Manufacturer sizing is modelled as measurement → size interpolation:
//! landmark-driven measurements (e.g. transepicondylar width) are mapped
//! through a size chart with linear interpolation between chart sizes, and
//! the alignment targets (distal femoral valgus angle, tibial slope) are
//! computed from landmark geometry following ISB conventions.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use tpt_med_geometry::Vec3;

/// One size in a manufacturer chart.
#[derive(Debug, Clone, Copy)]
pub struct SizeEntry {
    /// Size label (manufacturer-specific, e.g. 1–8 or S/M/L).
    pub label: u32,
    /// Chart measurement value (mm) — typically transepicondylar width or
    /// AP depth for femoral charts, plateau width for tibial.
    pub nominal: f64,
}

/// A sizing chart over an implant family.
#[derive(Debug, Clone)]
pub struct SizeChart {
    /// Implant family identifier (e.g. "tka-femoral-vectorzen").
    pub family: String,
    /// Ordered entries (ascending `nominal`).
    pub entries: Vec<SizeEntry>,
}

impl SizeChart {
    /// Selects the smallest size whose nominal ≥ measurement; when the
    /// measurement is between sizes, `SizeChoice` reports the interpolated
    /// position so callers can decide to size up/down.
    pub fn select(&self, measurement: f64) -> Option<SizeChoice> {
        if self.entries.is_empty() {
            return None;
        }
        let first = self.entries.first()?;
        let last = self.entries.last()?;
        if measurement <= first.nominal {
            return Some(SizeChoice {
                label: first.label,
                below_chart: true,
                above_chart: false,
                interpolation: 0.0,
            });
        }
        if measurement >= last.nominal {
            return Some(SizeChoice {
                label: last.label,
                below_chart: false,
                above_chart: true,
                interpolation: 0.0,
            });
        }
        for w in self.entries.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            if measurement >= a.nominal && measurement <= b.nominal {
                let t = (measurement - a.nominal) / (b.nominal - a.nominal);
                // Round-half-up to the nearer size, ties go up.
                let label = if t < 0.5 { a.label } else { b.label };
                return Some(SizeChoice {
                    label,
                    below_chart: false,
                    above_chart: false,
                    interpolation: t,
                });
            }
        }
        None
    }
}

/// Result of a chart lookup.
#[derive(Debug, Clone, Copy)]
pub struct SizeChoice {
    /// Selected size label.
    pub label: u32,
    /// Measurement below the smallest chart entry.
    pub below_chart: bool,
    /// Measurement above the largest chart entry.
    pub above_chart: bool,
    /// Fractional position between the bracketing sizes [0, 1).
    pub interpolation: f64,
}

/// Landmark-driven measurements for total knee arthroplasty sizing and
/// alignment.
#[derive(Debug, Clone)]
pub struct KneeLandmarks {
    /// Medial femoral epicondyle.
    pub medial_epicondyle: Vec3,
    /// Lateral femoral epicondyle.
    pub lateral_epicondyle: Vec3,
    /// Femoral AP-most distal point (trochlea bottom).
    pub trochlea_point: Vec3,
    /// Posterior condyle reference (most posterior femoral point).
    pub posterior_condyle: Vec3,
    /// Medial tibial plateau edge.
    pub tibial_medial: Vec3,
    /// Lateral tibial plateau edge.
    pub tibial_lateral: Vec3,
    /// Tibial midpoint at the resection level.
    pub tibial_center: Vec3,
    /// Most proximal point of tibial tubercle.
    pub tibial_tubercle: Vec3,
}

impl KneeLandmarks {
    /// Surgical transepicondylar axis (TEA) length (mm).
    pub fn tea_width(&self) -> f64 {
        (self.lateral_epicondyle - self.medial_epicondyle).norm()
    }

    /// Anterior–posterior depth: trochlea to posterior condyle projected on
    /// the TEA-normal direction (mm).
    pub fn ap_depth(&self) -> f64 {
        let tea = (self.lateral_epicondyle - self.medial_epicondyle).normalize();
        let d = self.posterior_condyle - self.trochlea_point;
        (d - tea * d.dot(tea)).norm()
    }

    /// Tibial plateau width (mm).
    pub fn plateau_width(&self) -> f64 {
        (self.tibial_lateral - self.tibial_medial).norm()
    }

    /// Femorotibial angle (screening alignment proxy): angle between the
    /// distal femoral axis (posterior condyle → trochlea) and the tibial
    /// axis (trochlea → tibial centre), in degrees. Neutral alignment ≈
    /// 0–10°; larger values indicate varus/valgus deformity.
    pub fn femorotibial_angle_deg(&self) -> f64 {
        let fem = (self.trochlea_point - self.posterior_condyle).normalize();
        let tib = (self.tibial_center - self.trochlea_point).normalize();
        (fem.dot(tib).clamp(-1.0, 1.0)).acos().to_degrees()
    }

    /// Posterior tibial slope: angle of the plateau line (medial→lateral is
    /// the ML axis; slope from tubercle-center offset), degrees.
    pub fn tibial_slope_deg(&self) -> f64 {
        let ml = (self.tibial_lateral - self.tibial_medial).normalize();
        let to_tubercle = self.tibial_tubercle - self.tibial_center;
        let along_ml = to_tubercle.dot(ml);
        let perp = (to_tubercle - ml * along_ml).norm();
        let base = along_ml.abs().max(1e-9);
        (perp / base).atan().to_degrees().clamp(0.0, 30.0)
    }
}

/// Sizing recommendation bundle for a TKA case.
#[derive(Debug, Clone, Copy)]
pub struct TkaSizing {
    /// Selected femoral component size label.
    pub femoral_size: u32,
    /// Selected tibial component size label.
    pub tibial_size: u32,
    /// TEA width measurement (mm).
    pub tea_width_mm: f64,
    /// AP depth measurement (mm).
    pub ap_depth_mm: f64,
    /// Plateau width measurement (mm).
    pub plateau_width_mm: f64,
    /// Femorotibial alignment angle (degrees).
    pub femorotibial_angle_deg: f64,
    /// Posterior tibial slope (degrees).
    pub tibial_slope_deg: f64,
}

/// Computes TKA sizing from landmarks and two charts.
pub fn size_tka(
    landmarks: &KneeLandmarks,
    femoral_chart: &SizeChart,
    tibial_chart: &SizeChart,
) -> Option<TkaSizing> {
    let tea = landmarks.tea_width();
    let ap = landmarks.ap_depth();
    let plateau = landmarks.plateau_width();
    let femoral = femoral_chart.select(0.6 * tea + 0.4 * ap)?.label;
    let tibial = tibial_chart.select(plateau)?.label;
    Some(TkaSizing {
        femoral_size: femoral,
        tibial_size: tibial,
        tea_width_mm: tea,
        ap_depth_mm: ap,
        plateau_width_mm: plateau,
        femorotibial_angle_deg: landmarks.femorotibial_angle_deg(),
        tibial_slope_deg: landmarks.tibial_slope_deg(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chart() -> SizeChart {
        SizeChart {
            family: "tka-femoral-test".into(),
            entries: (1..=8)
                .map(|i| SizeEntry {
                    label: i,
                    nominal: 50.0 + 4.0 * i as f64, // 54..82 mm
                })
                .collect(),
        }
    }

    fn landmarks() -> KneeLandmarks {
        KneeLandmarks {
            medial_epicondyle: Vec3::new(-35.0, 0.0, 0.0),
            lateral_epicondyle: Vec3::new(35.0, 0.0, 0.0), // TEA 70 mm
            trochlea_point: Vec3::new(0.0, 10.0, -60.0),
            posterior_condyle: Vec3::new(0.0, 10.0, 0.0),
            tibial_medial: Vec3::new(-30.0, 0.0, -80.0),
            tibial_lateral: Vec3::new(30.0, 0.0, -80.0), // 60 mm plateau
            tibial_center: Vec3::new(0.0, 0.0, -80.0),
            tibial_tubercle: Vec3::new(5.0, 0.0, -78.0),
        }
    }

    #[test]
    fn chart_selects_and_interpolates() {
        let c = chart();
        let mid = c.select(57.0).expect("in range"); // t = 0.75 → size 2
        assert_eq!(mid.label, 2);
        let between = c.select(55.0).expect("in range"); // t = 0.25 → size 1
        assert_eq!(between.label, 1);
        let tie = c.select(56.0).expect("in range"); // t = 0.5 → up-size 2
        assert_eq!(tie.label, 2);
        assert!(c.select(40.0).unwrap().below_chart);
        assert!(c.select(90.0).unwrap().above_chart);
        assert!(SizeChart {
            family: "e".into(),
            entries: vec![]
        }
        .select(50.0)
        .is_none());
    }

    #[test]
    fn landmark_measurements() {
        let l = landmarks();
        assert!((l.tea_width() - 70.0).abs() < 1e-9);
        assert!((l.plateau_width() - 60.0).abs() < 1e-9);
        assert!((l.ap_depth() - 60.0).abs() < 1e-9); // straight AP offset
    }

    #[test]
    fn tka_sizing_end_to_end() {
        let fem = chart();
        let tib = SizeChart {
            family: "tka-tibial-test".into(),
            entries: (1..=6)
                .map(|i| SizeEntry {
                    label: i,
                    nominal: 45.0 + 5.0 * i as f64, // 50..75
                })
                .collect(),
        };
        let s = size_tka(&landmarks(), &fem, &tib).expect("sizes");
        assert_eq!(s.tea_width_mm, 70.0);
        // Femoral index: 0.6·70 + 0.4·60 = 66 → chart nominal 66 is
        // exactly label 4 (50 + 4·4).
        assert_eq!(s.femoral_size, 4);
        // Plateau 60 → between 55(2) and 60(3) → label 3.
        assert_eq!(s.tibial_size, 3);
    }

    #[test]
    fn alignment_angles_in_clinical_bands() {
        let l = landmarks();
        let s = size_tka(&l, &chart(), &chart()).unwrap();
        // Screening geometry: angles finite and within plausible bands.
        assert!((0.0..=45.0).contains(&s.femorotibial_angle_deg));
        assert!((0.0..=30.0).contains(&s.tibial_slope_deg));
    }
}
