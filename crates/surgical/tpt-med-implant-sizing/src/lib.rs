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
}
