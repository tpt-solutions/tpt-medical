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

/// Chart validation failure, with the offending entry so a
/// mis-transcribed chart is caught rather than silently producing
/// recommendations.
#[derive(Debug, Clone, PartialEq)]
pub enum ChartError {
    /// No entries.
    Empty,
    /// A nominal value is not finite or not positive.
    NonPositiveNominal {
        /// Index of the offending entry.
        index: usize,
    },
    /// Nominal values are not strictly ascending.
    NotAscending {
        /// Index of the entry that broke the ordering.
        index: usize,
    },
    /// Duplicate size labels.
    DuplicateLabel {
        /// The duplicated label.
        label: u32,
    },
}

impl core::fmt::Display for ChartError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ChartError::Empty => write!(f, "chart has no entries"),
            ChartError::NonPositiveNominal { index } => {
                write!(f, "entry {index} has a non-finite or non-positive nominal")
            }
            ChartError::NotAscending { index } => {
                write!(f, "entry {index} breaks ascending nominal order")
            }
            ChartError::DuplicateLabel { label } => write!(f, "duplicate label {label}"),
        }
    }
}

impl std::error::Error for ChartError {}

impl SizeChart {
    /// Validates the chart: non-empty, finite positive nominals in strict
    /// ascending order, unique labels. [`Self::select`] implicitly assumes
    /// all of this; call this once when a chart is transcribed or loaded.
    pub fn validate(&self) -> Result<(), ChartError> {
        if self.entries.is_empty() {
            return Err(ChartError::Empty);
        }
        let mut seen = std::collections::BTreeSet::new();
        for (i, e) in self.entries.iter().enumerate() {
            if !e.nominal.is_finite() || e.nominal <= 0.0 {
                return Err(ChartError::NonPositiveNominal { index: i });
            }
            if i > 0 && e.nominal <= self.entries[i - 1].nominal {
                return Err(ChartError::NotAscending { index: i });
            }
            if !seen.insert(e.label) {
                return Err(ChartError::DuplicateLabel { label: e.label });
            }
        }
        Ok(())
    }

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

    /// Posterior condylar offset (screening proxy): the posterior condyle's
    /// perpendicular distance from the TEA line (mm). The surgical
    /// definition references the femoral anatomical axis, which these
    /// landmarks do not include — this proxy is the same *kind* of
    /// posterior-reach measurement and is documented as a proxy, not the
    /// clinical quantity.
    pub fn posterior_condylar_offset(&self) -> f64 {
        let tea = (self.lateral_epicondyle - self.medial_epicondyle).normalize();
        let from_medial = self.posterior_condyle - self.medial_epicondyle;
        let along = from_medial.dot(tea);
        (from_medial - tea * along).norm()
    }

    /// The three femoral measurements as sizing inputs with the canonical
    /// precedence order (TEA, then AP depth, then posterior condylar
    /// offset) — a common vendor convention, overridable by the caller.
    pub fn femoral_measurements(&self) -> [MeasurementInput; 3] {
        [
            MeasurementInput {
                name: "tea_width",
                value_mm: self.tea_width(),
                precedence: 1,
            },
            MeasurementInput {
                name: "ap_depth",
                value_mm: self.ap_depth(),
                precedence: 2,
            },
            MeasurementInput {
                name: "posterior_condylar_offset",
                value_mm: self.posterior_condylar_offset(),
                precedence: 3,
            },
        ]
    }
}

/// One measurement feeding a multi-measurement sizing decision.
#[derive(Debug, Clone, Copy)]
pub struct MeasurementInput {
    /// Stable measurement name for the audit record.
    pub name: &'static str,
    /// The measurement (mm).
    pub value_mm: f64,
    /// Precedence rank: **lower wins** when measurements disagree on a
    /// size; equal ranks resolve to the larger size (conservative
    /// up-size).
    pub precedence: u32,
}

/// How a multi-measurement decision resolved disagreement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Resolution {
    /// Every measurement selected the same size.
    Consensus,
    /// The lowest-precedence-rank measurement decided (named).
    Precedence(&'static str),
    /// Equal-precedence disagreement, resolved to the larger size.
    UpsizeTieBreak,
}

/// One measurement's vote, recorded in the decision.
#[derive(Debug, Clone, Copy)]
pub struct MeasurementVote {
    /// The measurement's name.
    pub name: &'static str,
    /// Its value (mm).
    pub value_mm: f64,
    /// The size it alone would select.
    pub label: u32,
    /// Its value fell below the chart.
    pub below_chart: bool,
    /// Its value rose above the chart.
    pub above_chart: bool,
}

/// A multi-measurement sizing decision: every measurement maps through the
/// chart, and an **explicit vendor precedence rule** resolves disagreement
/// — the point where a single blended index (like [`size_tka`]'s fixed
/// TEA/AP mix) hides which measurement actually drove the choice.
#[derive(Debug, Clone)]
pub struct MultiMeasurementDecision {
    /// The resolved size label.
    pub label: u32,
    /// How disagreement (if any) was resolved.
    pub resolved_by: Resolution,
    /// Every measurement's individual vote, in input order.
    pub votes: Vec<MeasurementVote>,
}

impl MultiMeasurementDecision {
    /// True when every measurement agreed on the selected size.
    pub fn is_unanimous(&self) -> bool {
        matches!(self.resolved_by, Resolution::Consensus)
    }
}

/// Sizes one component from several measurements jointly, resolving
/// disagreement by the explicit precedence rule: lowest precedence rank
/// wins; equal ranks up-size to the larger label.
pub fn size_from_measurements(
    chart: &SizeChart,
    inputs: &[MeasurementInput],
) -> Option<MultiMeasurementDecision> {
    let mut votes = Vec::with_capacity(inputs.len());
    for input in inputs {
        let Some(choice) = chart.select(input.value_mm) else {
            return None; // empty chart: nothing to vote on
        };
        votes.push(MeasurementVote {
            name: input.name,
            value_mm: input.value_mm,
            label: choice.label,
            below_chart: choice.below_chart,
            above_chart: choice.above_chart,
        });
    }
    if votes.is_empty() {
        return None;
    }

    let resolved_by = if votes.iter().all(|v| v.label == votes[0].label) {
        Resolution::Consensus
    } else {
        let best_rank = inputs
            .iter()
            .map(|i| i.precedence)
            .min()
            .expect("non-empty");
        let top: Vec<&MeasurementVote> = votes
            .iter()
            .filter(|v| {
                inputs
                    .iter()
                    .find(|i| i.name == v.name)
                    .is_some_and(|i| i.precedence == best_rank)
            })
            .collect();
        if top.len() == 1 {
            Resolution::Precedence(top[0].name)
        } else {
            Resolution::UpsizeTieBreak
        }
    };

    let label = match resolved_by {
        Resolution::Consensus => votes[0].label,
        Resolution::Precedence(name) => votes.iter().find(|v| v.name == name).expect("voted").label,
        Resolution::UpsizeTieBreak => votes.iter().map(|v| v.label).max().expect("non-empty"),
    };

    Some(MultiMeasurementDecision {
        label,
        resolved_by,
        votes,
    })
}

/// The resection and component-thickness inputs of a planned total knee:
/// how much bone each cut removes and how much implant each surface
/// replaces. All values in millimetres.
#[derive(Debug, Clone, Copy)]
pub struct ResectionPlan {
    /// Distal femoral resection (mm) — bone removed from the distal femur.
    pub distal_femoral_resection: f64,
    /// Posterior femoral resection (mm, per condyle) — bone removed from
    /// the posterior condyles.
    pub posterior_femoral_resection: f64,
    /// Proximal tibial resection (mm).
    pub tibial_resection: f64,
    /// Distal femoral component thickness (mm) — what the implant adds
    /// back at the distal femur.
    pub femoral_distal_thickness: f64,
    /// Posterior femoral condyle thickness (mm).
    pub femoral_posterior_thickness: f64,
    /// Total tibial component thickness (mm, insert + tray).
    pub tibial_component_thickness: f64,
}

/// The flexion/extension gap assessment for a [`ResectionPlan`].
#[derive(Debug, Clone, Copy)]
pub struct GapReport {
    /// Extension gap (mm): bone resected minus implant thickness at the
    /// distal femur and tibia — the space left after seating the
    /// components in extension. Negative means the components are
    /// **overstuffed** (thicker than the bone removed).
    pub extension_gap_mm: f64,
    /// Flexion gap (mm): the same balance at the posterior condyles and
    /// tibia.
    pub flexion_gap_mm: f64,
    /// `|flexion − extension|` (mm) — the imbalance a surgeon levels.
    pub imbalance_mm: f64,
    /// True when both gaps are non-negative (nothing overstuffed) and the
    /// imbalance is within the tolerance.
    pub is_balanced: bool,
}

/// Checks that a selected size leaves acceptable **gap balancing**: the
/// screening arithmetic of TKA mechanics — the extension gap is the distal
/// femoral and tibial resections minus the corresponding component
/// thicknesses, the flexion gap the posterior resections minus theirs.
/// A plan balances when neither gap went negative (an overstuffed
/// component, which lifts the joint line and tightens the collateral) and
/// the two gaps agree within `tolerance_mm` (a flexion/extension mismatch
/// that a soft-tissue release should not have to paper over).
///
/// This is the *checkable geometric half* of soft-tissue assessment:
/// actual ligament tension and stability need the soft-tissue structures
/// themselves, which this crate does not model.
pub fn check_gap_balance(plan: &ResectionPlan, tolerance_mm: f64) -> GapReport {
    let extension_gap = (plan.distal_femoral_resection - plan.femoral_distal_thickness)
        + (plan.tibial_resection - plan.tibial_component_thickness);
    let flexion_gap = (plan.posterior_femoral_resection - plan.femoral_posterior_thickness)
        + (plan.tibial_resection - plan.tibial_component_thickness);
    let imbalance = (flexion_gap - extension_gap).abs();
    GapReport {
        extension_gap_mm: extension_gap,
        flexion_gap_mm: flexion_gap,
        imbalance_mm: imbalance,
        is_balanced: extension_gap >= 0.0 && flexion_gap >= 0.0 && imbalance <= tolerance_mm,
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

/// Femoral-side landmarks for hip stem sizing: the canal geometry, the
/// head centre, and the lesser trochanter. All the stem-sizing
/// measurements derive from these geometrically; the size chart itself
/// stays caller-supplied (vendor data with citations), exactly as for
/// the knee.
#[derive(Debug, Clone, Copy)]
pub struct HipLandmarks {
    /// Centre of femoral head rotation.
    pub head_center: Vec3,
    /// Tip of the lesser trochanter (the femoral-side leg-length
    /// reference).
    pub lesser_trochanter: Vec3,
    /// Medial endosteal edge of the canal at the isthmus (the narrowest
    /// point).
    pub isthmus_medial: Vec3,
    /// Lateral endosteal edge of the canal at the isthmus.
    pub isthmus_lateral: Vec3,
    /// A point on the proximal canal axis (piriformis fossa / entry
    /// region) — anchors the canal axis the offset is measured against.
    pub canal_entry: Vec3,
}

impl HipLandmarks {
    /// The femoral canal axis: unit line through `canal_entry` toward the
    /// isthmus midpoint. The offset and the leg-length projection are
    /// measured against it.
    pub fn canal_axis(&self) -> Vec3 {
        let mid = (self.isthmus_medial + self.isthmus_lateral) * 0.5;
        (mid - self.canal_entry).normalize()
    }

    /// Endosteal canal width at the isthmus (mm): the medial–lateral
    /// isthmus edge separation. This drives the distal stem size.
    pub fn canal_width(&self) -> f64 {
        (self.isthmus_lateral - self.isthmus_medial).norm()
    }

    /// Femoral offset (mm): the perpendicular distance from the head
    /// centre to the canal axis — the standard geometric definition, and
    /// the measurement that decides the proximal body / head offset.
    pub fn femoral_offset(&self) -> f64 {
        let axis = self.canal_axis();
        let d = self.head_center - self.canal_entry;
        (d - axis * d.dot(axis)).norm()
    }

    /// Femoral leg length (mm): the head centre's axial projection past
    /// the lesser trochanter along the canal axis. A femoral-side proxy
    /// for the leg-length measurement (the clinical definition references
    /// a pelvis landmark these landmarks do not include) — documented as
    /// a proxy.
    pub fn head_to_lesser_trochanter(&self) -> f64 {
        let axis = self.canal_axis();
        let d = self.head_center - self.lesser_trochanter;
        d.dot(axis).abs()
    }

    /// The three stem-sizing measurements with the canonical precedence
    /// order (canal width, then offset, then leg length).
    pub fn stem_measurements(&self) -> [MeasurementInput; 3] {
        [
            MeasurementInput {
                name: "canal_width",
                value_mm: self.canal_width(),
                precedence: 1,
            },
            MeasurementInput {
                name: "femoral_offset",
                value_mm: self.femoral_offset(),
                precedence: 2,
            },
            MeasurementInput {
                name: "head_to_lesser_trochanter",
                value_mm: self.head_to_lesser_trochanter(),
                precedence: 3,
            },
        ]
    }
}

/// Sizes a femoral stem from hip landmarks against a caller-supplied
/// stem chart: the same precedence-resolved multi-measurement decision
/// as the knee path, over the canal-width / offset / leg-length triple.
///
/// # Errors
///
/// Returns `None` when the chart is empty or a measurement falls outside
/// every chart entry's reach (mirroring [`size_from_measurements`]).
pub fn size_hip(chart: &SizeChart, landmarks: &HipLandmarks) -> Option<MultiMeasurementDecision> {
    size_from_measurements(chart, &landmarks.stem_measurements())
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
    fn chart_validation_catches_transcription_errors() {
        let ok = SizeChart {
            family: "ok".into(),
            entries: vec![
                SizeEntry {
                    label: 1,
                    nominal: 50.0,
                },
                SizeEntry {
                    label: 2,
                    nominal: 55.0,
                },
            ],
        };
        assert!(ok.validate().is_ok());

        let empty = SizeChart {
            family: "e".into(),
            entries: vec![],
        };
        assert_eq!(empty.validate(), Err(ChartError::Empty));

        let zero = SizeChart {
            family: "z".into(),
            entries: vec![SizeEntry {
                label: 1,
                nominal: 0.0,
            }],
        };
        assert_eq!(
            zero.validate(),
            Err(ChartError::NonPositiveNominal { index: 0 })
        );

        let unsorted = SizeChart {
            family: "u".into(),
            entries: vec![
                SizeEntry {
                    label: 1,
                    nominal: 55.0,
                },
                SizeEntry {
                    label: 2,
                    nominal: 50.0,
                },
            ],
        };
        assert_eq!(
            unsorted.validate(),
            Err(ChartError::NotAscending { index: 1 })
        );

        let dup = SizeChart {
            family: "d".into(),
            entries: vec![
                SizeEntry {
                    label: 3,
                    nominal: 50.0,
                },
                SizeEntry {
                    label: 3,
                    nominal: 55.0,
                },
            ],
        };
        assert_eq!(dup.validate(), Err(ChartError::DuplicateLabel { label: 3 }));
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
    fn gap_balance_flags_overstuffing_and_flexion_mismatch() {
        // Resections equal to component thicknesses: both gaps zero,
        // balanced.
        let plan = ResectionPlan {
            distal_femoral_resection: 9.0,
            posterior_femoral_resection: 10.0,
            tibial_resection: 10.0,
            femoral_distal_thickness: 9.0,
            femoral_posterior_thickness: 10.0,
            tibial_component_thickness: 10.0,
        };
        let ok = check_gap_balance(&plan, 2.0);
        assert!((ok.extension_gap_mm - 0.0).abs() < 1e-12);
        assert!((ok.flexion_gap_mm - 0.0).abs() < 1e-12);
        assert!(ok.is_balanced);

        // A thicker tibial insert stuffs both gaps negative.
        let stuffed = ResectionPlan {
            tibial_component_thickness: 13.0,
            ..plan
        };
        let r = check_gap_balance(&stuffed, 2.0);
        assert!(r.extension_gap_mm < 0.0 && r.flexion_gap_mm < 0.0);
        assert!(!r.is_balanced, "overstuffed must fail");

        // A posterior-heavy femoral resection opens the flexion gap
        // relative to extension: flagged by the imbalance tolerance.
        let flexion_open = ResectionPlan {
            posterior_femoral_resection: 14.0,
            ..plan
        };
        let r = check_gap_balance(&flexion_open, 2.0);
        assert!((r.imbalance_mm - 4.0).abs() < 1e-12);
        assert!(!r.is_balanced);
        // …and accepted at a looser tolerance.
        assert!(check_gap_balance(&flexion_open, 4.0).is_balanced);
    }

    #[test]
    fn posterior_condylar_offset_is_the_tea_perpendicular() {
        let l = landmarks();
        // posterior_condyle is (0,10,0); the TEA runs along x through
        // (±35,0,0), so the perpendicular offset is sqrt(0² + 10²) = 10.
        assert!((l.posterior_condylar_offset() - 10.0).abs() < 1e-9);
        // The convenience inputs carry the canonical precedence order.
        let inputs = l.femoral_measurements();
        assert_eq!(inputs[0].name, "tea_width");
        assert!(inputs[0].precedence < inputs[1].precedence);
        assert!(inputs[1].precedence < inputs[2].precedence);
    }

    #[test]
    fn unanimous_measurements_resolve_by_consensus() {
        let c = chart();
        // 58, 58.5, 57.5 all sit inside size 2's bracket.
        let decision = size_from_measurements(
            &c,
            &[
                MeasurementInput {
                    name: "tea_width",
                    value_mm: 58.0,
                    precedence: 1,
                },
                MeasurementInput {
                    name: "ap_depth",
                    value_mm: 58.5,
                    precedence: 2,
                },
                MeasurementInput {
                    name: "posterior_condylar_offset",
                    value_mm: 57.5,
                    precedence: 3,
                },
            ],
        )
        .expect("decision");
        assert_eq!(decision.label, 2);
        assert_eq!(decision.resolved_by, Resolution::Consensus);
        assert!(decision.is_unanimous());
        assert_eq!(decision.votes.len(), 3);
    }

    #[test]
    fn disagreement_resolves_by_precedence_then_upsize() {
        let c = chart();
        // TEA (rank 1) votes label 5, AP (rank 2) votes label 7, PCO
        // (rank 3) votes label 1: the TEA decides.
        let by_precedence = size_from_measurements(
            &c,
            &[
                MeasurementInput {
                    name: "tea_width",
                    value_mm: 70.0,
                    precedence: 1,
                },
                MeasurementInput {
                    name: "ap_depth",
                    value_mm: 78.0,
                    precedence: 2,
                },
                MeasurementInput {
                    name: "posterior_condylar_offset",
                    value_mm: 54.0,
                    precedence: 3,
                },
            ],
        )
        .expect("decision");
        assert_eq!(by_precedence.label, 5);
        assert_eq!(
            by_precedence.resolved_by,
            Resolution::Precedence("tea_width")
        );
        assert!(!by_precedence.is_unanimous());
        // Same rank on the two disagreeing inputs: conservative up-size.
        let tie = size_from_measurements(
            &c,
            &[
                MeasurementInput {
                    name: "tea_width",
                    value_mm: 70.0,
                    precedence: 1,
                },
                MeasurementInput {
                    name: "ap_depth",
                    value_mm: 90.0,
                    precedence: 1,
                },
            ],
        )
        .expect("decision");
        assert_eq!(tie.label, 8);
        assert_eq!(tie.resolved_by, Resolution::UpsizeTieBreak);
        // Out-of-chart votes are recorded, not silently clamped.
        assert!(tie.votes.iter().any(|v| v.above_chart));
    }

    #[test]
    fn alignment_angles_in_clinical_bands() {
        let l = landmarks();
        let s = size_tka(&l, &chart(), &chart()).unwrap();
        // Screening geometry: angles finite and within plausible bands.
        assert!((0.0..=45.0).contains(&s.femorotibial_angle_deg));
        assert!((0.0..=30.0).contains(&s.tibial_slope_deg));
    }

    #[test]
    fn hip_offset_and_canal_width_are_exact_geometric_measurements() {
        // Canal axis along +z through (0, 0), head offset 4 mm in +x,
        // isthmus width 12 mm in x, lesser trochanter 3 mm below the head
        // along the axis.
        let hip = HipLandmarks {
            head_center: Vec3::new(4.0, 0.0, 5.0),
            lesser_trochanter: Vec3::new(0.0, 0.0, 2.0),
            isthmus_medial: Vec3::new(-6.0, 0.0, 12.0),
            isthmus_lateral: Vec3::new(6.0, 0.0, 12.0),
            canal_entry: Vec3::new(0.0, 0.0, 0.0),
        };
        assert!((hip.canal_width() - 12.0).abs() < 1e-12);
        // The offset is the |x| distance: the head sits 4 mm off the axis.
        assert!((hip.femoral_offset() - 4.0).abs() < 1e-12);
        // The axis direction is +z: axial head-to-trochanter = 3 mm.
        assert!((hip.head_to_lesser_trochanter() - 3.0).abs() < 1e-12);
        // The canal axis is unit and along +z.
        let axis = hip.canal_axis();
        assert!((axis.norm() - 1.0).abs() < 1e-12);
        assert!(axis.z > 0.99);
    }

    #[test]
    fn hip_sizing_resolves_by_precedence_like_the_knee() {
        let chart = SizeChart {
            family: "hip-stem-test".into(),
            entries: vec![
                SizeEntry {
                    label: 1,
                    nominal: 10.0,
                },
                SizeEntry {
                    label: 2,
                    nominal: 12.0,
                },
                SizeEntry {
                    label: 3,
                    nominal: 14.0,
                },
                SizeEntry {
                    label: 4,
                    nominal: 16.0,
                },
            ],
        };
        chart.validate().expect("chart valid");

        // Consensus: all three measurements land inside the chart near
        // label 3 (nominal 14). The head sits 14 mm off the canal axis,
        // the isthmus is 14 mm wide, and the head trails the lesser
        // trochanter by 13 mm along the axis — offsets the chart can see.
        let hip = HipLandmarks {
            head_center: Vec3::new(14.0, 0.0, 5.0),
            lesser_trochanter: Vec3::new(0.0, 0.0, -8.0),
            isthmus_medial: Vec3::new(-7.0, 0.0, 12.0),
            isthmus_lateral: Vec3::new(7.0, 0.0, 12.0),
            canal_entry: Vec3::new(0.0, 0.0, 0.0),
        };
        let decision = size_hip(&chart, &hip).expect("sizes");
        assert!(decision.is_unanimous());
        assert_eq!(decision.label, 3);

        // Disagreement: canal width votes low (label 1 at 10.5 mm),
        // offset votes high (label 3) — the canal width wins by
        // precedence.
        let narrow = HipLandmarks {
            isthmus_medial: Vec3::new(-5.25, 0.0, 12.0),
            isthmus_lateral: Vec3::new(5.25, 0.0, 12.0),
            ..hip
        };
        let decision = size_hip(&chart, &narrow).expect("sizes");
        assert_eq!(decision.label, 1);
        assert!(matches!(
            decision.resolved_by,
            Resolution::Precedence("canal_width")
        ));
        assert_eq!(decision.votes.len(), 3);
    }

    #[test]
    fn hip_sizing_is_monotone_in_canal_width() {
        let chart = SizeChart {
            family: "hip-stem-test".into(),
            entries: vec![
                SizeEntry {
                    label: 1,
                    nominal: 10.0,
                },
                SizeEntry {
                    label: 2,
                    nominal: 12.0,
                },
                SizeEntry {
                    label: 3,
                    nominal: 14.0,
                },
                SizeEntry {
                    label: 4,
                    nominal: 16.0,
                },
            ],
        };
        let base = HipLandmarks {
            head_center: Vec3::new(4.0, 0.0, 5.0),
            lesser_trochanter: Vec3::new(0.0, 0.0, 2.0),
            isthmus_medial: Vec3::new(-6.0, 0.0, 12.0),
            isthmus_lateral: Vec3::new(6.0, 0.0, 12.0),
            canal_entry: Vec3::new(0.0, 0.0, 0.0),
        };
        // Widening the canal never decreases the selected stem size.
        let mut previous = 0;
        for width in [10.5, 11.5, 12.5, 13.5, 14.5, 15.5] {
            let half = width / 2.0;
            let hip = HipLandmarks {
                isthmus_medial: Vec3::new(-half, 0.0, 12.0),
                isthmus_lateral: Vec3::new(half, 0.0, 12.0),
                ..base
            };
            let label = size_hip(&chart, &hip).expect("sizes").label;
            assert!(label >= previous, "width {width}: {label} after {previous}");
            previous = label;
        }
    }
}
