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

/// Humeral landmarks for shoulder stem sizing: the canal geometry, the
/// head centre and anatomical neck, and the elbow epicondyles that anchor
/// the distal (transepicondylar) retroversion reference. All stem-sizing
/// measurements derive from these geometrically; the chart itself stays
/// caller-supplied (vendor data with citations), exactly as for the knee
/// and hip.
#[derive(Debug, Clone, Copy)]
pub struct ShoulderLandmarks {
    /// Centre of the humeral head.
    pub head_center: Vec3,
    /// Centre of the anatomical neck (the base of the head — the
    /// resection-plane reference).
    pub anatomical_neck_center: Vec3,
    /// Medial endosteal edge of the diaphyseal canal at the isthmus (the
    /// narrowest point).
    pub isthmus_medial: Vec3,
    /// Lateral endosteal edge of the canal at the isthmus.
    pub isthmus_lateral: Vec3,
    /// A point on the proximal canal axis (entry region) — anchors the
    /// canal axis the offset and angles are measured against.
    pub canal_entry: Vec3,
    /// Medial humeral epicondyle — one end of the distal retroversion
    /// reference.
    pub medial_epicondyle: Vec3,
    /// Lateral humeral epicondyle.
    pub lateral_epicondyle: Vec3,
}

impl ShoulderLandmarks {
    /// The humeral canal axis: unit line through `canal_entry` toward the
    /// isthmus midpoint (pointing distal). The offset and both alignment
    /// angles are measured against it.
    pub fn canal_axis(&self) -> Vec3 {
        let mid = (self.isthmus_medial + self.isthmus_lateral) * 0.5;
        (mid - self.canal_entry).normalize()
    }

    /// Endosteal canal width at the isthmus (mm). This drives the stem
    /// size.
    pub fn canal_width(&self) -> f64 {
        (self.isthmus_lateral - self.isthmus_medial).norm()
    }

    /// The neck axis: unit direction from the anatomical neck centre up
    /// to the head centre (pointing proximal).
    pub fn neck_axis(&self) -> Vec3 {
        (self.head_center - self.anatomical_neck_center).normalize()
    }

    /// Head height (mm): the head centre's axial height above the
    /// anatomical neck centre along the canal axis. A proxy for the
    /// calcar-to-apex head height (the clinical definition references
    /// the calcar and the head apex, which these landmarks do not
    /// include) — documented as a proxy.
    pub fn head_height(&self) -> f64 {
        (self.head_center - self.anatomical_neck_center)
            .dot(self.canal_axis())
            .abs()
    }

    /// Head offset (mm): the perpendicular distance from the head centre
    /// to the canal axis — the medial–lateral head position the implant
    /// restores, the standard geometric definition.
    pub fn head_offset(&self) -> f64 {
        let axis = self.canal_axis();
        let d = self.head_center - self.canal_entry;
        (d - axis * d.dot(axis)).norm()
    }

    /// Neck–shaft angle (degrees): the angle between the neck axis
    /// (pointing proximal, toward the head) and the distal canal axis.
    /// The normal band is ≈ 130–140°; a screening alignment output, not
    /// a sizing input.
    pub fn neck_shaft_angle_deg(&self) -> f64 {
        self.neck_axis()
            .dot(self.canal_axis())
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    }

    /// Humeral retroversion (degrees): the angle, in the plane
    /// perpendicular to the canal axis, between the projected neck axis
    /// and the projected transepicondylar axis — the standard distal
    /// reference for humeral version. The normal band is ≈ 20–30°;
    /// reported as a magnitude (screening output, not a sizing input).
    pub fn retroversion_deg(&self) -> f64 {
        let axis = self.canal_axis();
        let project = |v: Vec3| v - axis * v.dot(axis);
        let neck = project(self.neck_axis());
        let tea = project(self.lateral_epicondyle - self.medial_epicondyle);
        neck.normalize()
            .dot(tea.normalize())
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    }

    /// The three stem-sizing measurements with the canonical precedence
    /// order (canal width, then head offset, then head height).
    pub fn stem_measurements(&self) -> [MeasurementInput; 3] {
        [
            MeasurementInput {
                name: "canal_width",
                value_mm: self.canal_width(),
                precedence: 1,
            },
            MeasurementInput {
                name: "head_offset",
                value_mm: self.head_offset(),
                precedence: 2,
            },
            MeasurementInput {
                name: "head_height",
                value_mm: self.head_height(),
                precedence: 3,
            },
        ]
    }
}

/// Shoulder sizing recommendation bundle: the precedence-resolved stem
/// decision plus the measurements and both alignment outputs.
#[derive(Debug, Clone)]
pub struct ShoulderSizing {
    /// The stem decision (votes, resolution, label).
    pub stem: MultiMeasurementDecision,
    /// Canal width at the isthmus (mm).
    pub canal_width_mm: f64,
    /// Head offset from the canal axis (mm).
    pub head_offset_mm: f64,
    /// Head height above the neck (mm).
    pub head_height_mm: f64,
    /// Neck–shaft angle (degrees).
    pub neck_shaft_angle_deg: f64,
    /// Humeral retroversion relative to the transepicondylar axis
    /// (degrees).
    pub retroversion_deg: f64,
}

/// Sizes a humeral stem from shoulder landmarks against a
/// caller-supplied stem chart: the same precedence-resolved
/// multi-measurement decision as the knee and hip paths, over the
/// canal-width / head-offset / head-height triple, with the neck–shaft
/// angle and retroversion reported alongside.
///
/// # Errors
///
/// Returns `None` when the chart is empty or a measurement falls outside
/// every chart entry's reach (mirroring [`size_from_measurements`]).
pub fn size_shoulder(chart: &SizeChart, landmarks: &ShoulderLandmarks) -> Option<ShoulderSizing> {
    Some(ShoulderSizing {
        stem: size_from_measurements(chart, &landmarks.stem_measurements())?,
        canal_width_mm: landmarks.canal_width(),
        head_offset_mm: landmarks.head_offset(),
        head_height_mm: landmarks.head_height(),
        neck_shaft_angle_deg: landmarks.neck_shaft_angle_deg(),
        retroversion_deg: landmarks.retroversion_deg(),
    })
}

/// Landmarks for total ankle arthroplasty sizing: the tibial plafond
/// edges (medial/lateral for width, anterior/posterior for depth), the
/// talar dome edges and centre, and a proximal point on the tibial
/// anatomical axis for the alignment measurement. The charts stay
/// caller-supplied (vendor data with citations), exactly as for the
/// knee, hip and shoulder.
#[derive(Debug, Clone, Copy)]
pub struct AnkleLandmarks {
    /// Medial edge of the tibial plafond (distal tibia articular
    /// surface).
    pub plafond_medial: Vec3,
    /// Lateral edge of the tibial plafond.
    pub plafond_lateral: Vec3,
    /// Anterior edge of the plafond at its mid-width.
    pub plafond_anterior: Vec3,
    /// Posterior edge of the plafond at its mid-width.
    pub plafond_posterior: Vec3,
    /// Medial edge of the talar dome.
    pub talar_medial: Vec3,
    /// Lateral edge of the talar dome.
    pub talar_lateral: Vec3,
    /// Centre of the talar dome.
    pub talar_center: Vec3,
    /// A point proximal on the tibial anatomical axis — anchors the
    /// tibial axis the alignment angle is measured from.
    pub tibial_axis_point: Vec3,
}

impl AnkleLandmarks {
    /// Tibial plafond medial–lateral width (mm) — drives the tibial
    /// component size.
    pub fn plafond_width(&self) -> f64 {
        (self.plafond_lateral - self.plafond_medial).norm()
    }

    /// Tibial plafond anterior–posterior depth (mm) — the second tibial
    /// sizing measurement.
    pub fn plafond_depth(&self) -> f64 {
        (self.plafond_posterior - self.plafond_anterior).norm()
    }

    /// Talar dome medial–lateral width (mm) — drives the talar
    /// component size.
    pub fn talar_width(&self) -> f64 {
        (self.talar_lateral - self.talar_medial).norm()
    }

    /// Tibiotalar alignment deviation (degrees): the angle between the
    /// tibial anatomical axis and the talar dome line's perpendicular —
    /// zero when the axis is perpendicular to the dome (neutral), larger
    /// with varus/valgus tilt. A frontal-plane screening proxy computed
    /// from the dome line, not the full clinical tibiotalar angle.
    pub fn tibiotalar_angle_deg(&self) -> f64 {
        let plafond_center = (self.plafond_medial + self.plafond_lateral) * 0.5;
        let axis = (plafond_center - self.tibial_axis_point).normalize();
        let dome = (self.talar_lateral - self.talar_medial).normalize();
        axis.dot(dome).clamp(-1.0, 1.0).abs().asin().to_degrees()
    }

    /// The tibial component's two sizing measurements with the canonical
    /// precedence order (plafond width, then plafond depth).
    pub fn tibial_measurements(&self) -> [MeasurementInput; 2] {
        [
            MeasurementInput {
                name: "plafond_width",
                value_mm: self.plafond_width(),
                precedence: 1,
            },
            MeasurementInput {
                name: "plafond_depth",
                value_mm: self.plafond_depth(),
                precedence: 2,
            },
        ]
    }

    /// The talar component's sizing measurement (dome width).
    pub fn talar_measurements(&self) -> [MeasurementInput; 1] {
        [MeasurementInput {
            name: "talar_width",
            value_mm: self.talar_width(),
            precedence: 1,
        }]
    }
}

/// Ankle sizing recommendation bundle: precedence-resolved decisions for
/// both components plus the measurements and the alignment deviation.
#[derive(Debug, Clone)]
pub struct AnkleSizing {
    /// The tibial component decision (votes, resolution, label).
    pub tibial: MultiMeasurementDecision,
    /// The talar component decision (votes, resolution, label).
    pub talar: MultiMeasurementDecision,
    /// Plafond medial–lateral width (mm).
    pub plafond_width_mm: f64,
    /// Plafond anterior–posterior depth (mm).
    pub plafond_depth_mm: f64,
    /// Talar dome width (mm).
    pub talar_width_mm: f64,
    /// Tibiotalar alignment deviation (degrees).
    pub tibiotalar_angle_deg: f64,
}

/// Sizes a tibial and talar component from ankle landmarks against two
/// caller-supplied charts: the tibial decision weighs plafond width and
/// depth (width takes precedence), the talar decision maps the dome
/// width, and the tibiotalar deviation is reported alongside.
///
/// # Errors
///
/// Returns `None` when either chart is empty or a measurement falls
/// outside every chart entry's reach (mirroring
/// [`size_from_measurements`]).
pub fn size_ankle(
    tibial_chart: &SizeChart,
    talar_chart: &SizeChart,
    landmarks: &AnkleLandmarks,
) -> Option<AnkleSizing> {
    Some(AnkleSizing {
        tibial: size_from_measurements(tibial_chart, &landmarks.tibial_measurements())?,
        talar: size_from_measurements(talar_chart, &landmarks.talar_measurements())?,
        plafond_width_mm: landmarks.plafond_width(),
        plafond_depth_mm: landmarks.plafond_depth(),
        talar_width_mm: landmarks.talar_width(),
        tibiotalar_angle_deg: landmarks.tibiotalar_angle_deg(),
    })
}

/// Which side of the joint a collateral ligament spans — the axis of the
/// mediolateral imbalance the balance screen reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CollateralSide {
    /// Medial (tibial) collateral.
    Medial,
    /// Lateral (fibular) collateral.
    Lateral,
}

/// One collateral ligament as a **1-D spring with a slack range** — the
/// minimal soft-tissue structure the gap-balance check alone was missing:
/// origin and insertion (patient space, mm), the length below which the
/// fibre bundle is slack, and the stiffness beyond it. The tension law is
/// deliberately the screening piecewise-linear one — zero below
/// `slack_length_mm`, `stiffness_n_per_mm` times elongation above — not a
/// viscoelastic or fibre-bundle model.
#[derive(Debug, Clone)]
pub struct LigamentModel {
    /// Attachment on the proximal bone (femoral origin), patient space.
    pub origin: Vec3,
    /// Attachment on the distal bone (tibial insertion), patient space.
    pub insertion: Vec3,
    /// Unloaded (anatomical) length (mm): with no resection or component
    /// placement the ligament sits at this length — by definition the
    /// pre-operative state.
    pub slack_length_mm: f64,
    /// Stiffness beyond slack (N/mm).
    pub stiffness_n_per_mm: f64,
    /// Which side of the joint this collateral spans.
    pub side: CollateralSide,
}

impl LigamentModel {
    /// The anatomical (attachment-to-attachment) length (mm).
    pub fn anatomical_length(&self) -> f64 {
        (self.insertion - self.origin).norm()
    }

    /// The ligament's length (mm) when the plan has moved the attachments
    /// — a resection lowers the distal bone, a component thickens the
    /// joint line, a fragment move relocates an attachment. Pass the
    /// planned positions of both attachments; the unplanned call is
    /// [`Self::anatomical_length`].
    pub fn planned_length(&self, origin: Vec3, insertion: Vec3) -> f64 {
        (insertion - origin).norm()
    }

    /// Tension (N) at `length_mm`: zero while slack, linear beyond.
    pub fn tension(&self, length_mm: f64) -> f64 {
        let elongation = length_mm - self.slack_length_mm;
        if elongation <= 0.0 {
            0.0
        } else {
            self.stiffness_n_per_mm * elongation
        }
    }

    /// True when the ligament carries tension at `length_mm`.
    pub fn is_taut(&self, length_mm: f64) -> bool {
        length_mm > self.slack_length_mm
    }
}

/// The balance verdict for one collateral: anatomical vs planned length,
/// the elongation, the tension, and whether it is taut at all.
#[derive(Debug, Clone)]
pub struct LigamentTension {
    /// Anatomical length (mm).
    pub anatomical_length_mm: f64,
    /// Planned length (mm).
    pub planned_length_mm: f64,
    /// `planned − anatomical` (mm): positive lengthens (and tensions) the
    /// collateral, negative slackens it.
    pub elongation_mm: f64,
    /// Tension (N) at the planned length.
    pub tension_n: f64,
    /// True when the planned length exceeds slack.
    pub is_taut: bool,
}

/// The whole-joint ligament balance: every collateral's verdict plus the
/// **mediolateral imbalance** — `|Σ medial tensions − Σ lateral
/// tensions|`, the net unbalanced force a surgeon levels with releases or
/// insert thickness.
#[derive(Debug, Clone)]
pub struct LigamentBalanceReport {
    /// Per-ligament verdicts, in input order.
    pub ligaments: Vec<LigamentTension>,
    /// `|Σ medial − Σ lateral|` tensions (N).
    pub mediolateral_imbalance_n: f64,
    /// True when the imbalance is within `tolerance_n`. A slack collateral
    /// carries no tension, so an over-released side balances numerically
    /// while leaving the joint loose — pair this with the per-ligament
    /// `is_taut` flags rather than reading the flag alone.
    pub is_balanced: bool,
}

/// The ligament-balance assessment proper — the soft-tissue half the
/// geometric [`check_gap_balance`] deliberately left out. Each entry pairs
/// a [`LigamentModel`] with the **planned positions of its two
/// attachments** (what the resection depths, component thicknesses and
/// fragment moves of the plan do to them); the screen reports the length
/// change, the spring tension, and the mediolateral imbalance against
/// `tolerance_n`.
pub fn assess_ligament_balance(
    ligaments: &[(LigamentModel, Vec3, Vec3)],
    tolerance_n: f64,
) -> Option<LigamentBalanceReport> {
    if ligaments.is_empty() || !tolerance_n.is_finite() || tolerance_n < 0.0 {
        return None;
    }
    let mut verdicts = Vec::with_capacity(ligaments.len());
    let mut medial = 0.0f64;
    let mut lateral = 0.0f64;
    for (lig, origin, insertion) in ligaments {
        let planned = lig.planned_length(*origin, *insertion);
        let anatomical = lig.anatomical_length();
        let tension = lig.tension(planned);
        verdicts.push(LigamentTension {
            anatomical_length_mm: anatomical,
            planned_length_mm: planned,
            elongation_mm: planned - anatomical,
            tension_n: tension,
            is_taut: lig.is_taut(planned),
        });
        match lig.side {
            CollateralSide::Medial => medial += tension,
            CollateralSide::Lateral => lateral += tension,
        }
    }
    let imbalance = (medial - lateral).abs();
    Some(LigamentBalanceReport {
        ligaments: verdicts,
        mediolateral_imbalance_n: imbalance,
        is_balanced: imbalance <= tolerance_n,
    })
}

/// The **stability screen**: the tension balance of
/// [`assess_ligament_balance`] converted into a varus/valgus moment
/// statement. Each collateral pulls along its own line of action, so a
/// tension imbalance only *stabilises* the joint if the lines' moment
/// arms about the joint centre make it — two collaterals with equal
/// tension but unequal arms leave a net moment, and a well-armed pair
/// can carry an asymmetric load. The screen reports each ligament's
/// moment arm (the perpendicular distance from the joint centre to the
/// ligament's line of action) and the net moment
/// `Σ (medial arm × medial tension − lateral arm × lateral tension)` at
/// the plan's attachment positions.
///
/// Still deliberately out of scope: the *response* to a prescribed
/// stress — the kinematics of how the joint opens under a moment, which
/// needs the articulating surfaces this crate does not model. This screen
/// says whether the *plan's* ligament configuration leaves a net moment;
/// it does not simulate the opening that would result.
pub struct StabilityScreen;

impl StabilityScreen {
    /// A collateral's moment arm (mm): the perpendicular distance from
    /// `joint_center` to the ligament's line of action through its
    /// planned attachments.
    pub fn moment_arm_mm(
        ligament: &LigamentModel,
        joint_center: Vec3,
        planned_origin: Vec3,
        planned_insertion: Vec3,
    ) -> f64 {
        let axis = (planned_insertion - planned_origin)
            * (1.0 / ligament.planned_length(planned_origin, planned_insertion));
        let w = joint_center - planned_origin;
        (w.cross(axis)).norm()
    }

    /// The net varus/valgus moment (N·mm) at the plan's attachment
    /// positions: positive when the medial side's moment contribution
    /// (its arm times its planned tension) dominates — the convention is
    /// documented as a magnitude with sign by side pairing, since the
    /// clinical varus/valgus sign depends on the joint's left/right
    /// orientation. Arms and tensions both come from the plan geometry.
    pub fn net_moment_nmm(
        ligaments: &[(LigamentModel, Vec3, Vec3)],
        joint_center: Vec3,
    ) -> Option<f64> {
        if ligaments.is_empty() {
            return None;
        }
        let mut moment = 0.0f64;
        for (lig, origin, insertion) in ligaments {
            let arm = Self::moment_arm_mm(lig, joint_center, *origin, *insertion);
            let tension = lig.tension(lig.planned_length(*origin, *insertion));
            let signed = match lig.side {
                CollateralSide::Medial => arm * tension,
                CollateralSide::Lateral => -(arm * tension),
            };
            moment += signed;
        }
        Some(moment)
    }

    /// Full screen: per-ligament arms and the net moment.
    pub fn assess(
        ligaments: &[(LigamentModel, Vec3, Vec3)],
        joint_center: Vec3,
    ) -> Option<StabilityReport> {
        if ligaments.is_empty() {
            return None;
        }
        let mut arms = Vec::with_capacity(ligaments.len());
        for (lig, origin, insertion) in ligaments {
            arms.push(Self::moment_arm_mm(lig, joint_center, *origin, *insertion));
        }
        let net = Self::net_moment_nmm(ligaments, joint_center)?;
        Some(StabilityReport {
            moment_arms_mm: arms,
            net_moment_nmm: net,
        })
    }
}

/// The stability screen's report: per-ligament moment arms (input order)
/// and the net varus/valgus moment (N·mm, medial-positive by side
/// pairing).
#[derive(Debug, Clone)]
pub struct StabilityReport {
    /// Moment arm (mm) per ligament, input order.
    pub moment_arms_mm: Vec<f64>,
    /// Net moment (N·mm): `Σ ±arm·tension`, medial-positive.
    pub net_moment_nmm: f64,
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

    /// Canal along +z through the origin, isthmus 12 mm wide at z = 12;
    /// neck axis 10 mm long at 135° to the distal axis, its in-plane
    /// component at `retroversion` degrees from the +x TEA.
    fn shoulder_landmarks(retroversion_deg: f64) -> ShoulderLandmarks {
        let s = 45.0f64.to_radians().sin(); // in-plane magnitude of the neck axis
        let rv = retroversion_deg.to_radians();
        // Unit neck axis: in-plane component at angle rv from the TEA
        // (+x), axial component −s so the neck–shaft angle is 135°.
        let neck_axis = Vec3::new(s * rv.cos(), s * rv.sin(), -s);
        ShoulderLandmarks {
            anatomical_neck_center: Vec3::new(0.0, 0.0, 2.0),
            head_center: Vec3::new(0.0, 0.0, 2.0) + neck_axis * 10.0,
            isthmus_medial: Vec3::new(-6.0, 0.0, 12.0),
            isthmus_lateral: Vec3::new(6.0, 0.0, 12.0),
            canal_entry: Vec3::new(0.0, 0.0, 0.0),
            medial_epicondyle: Vec3::new(-20.0, 0.0, -10.0),
            lateral_epicondyle: Vec3::new(20.0, 0.0, -10.0),
        }
    }

    #[test]
    fn shoulder_measurements_are_exact_geometry() {
        // Neck axis at 135° to the distal canal axis, in-plane component
        // along +x (TEA along x): retroversion is exactly zero.
        let sh = shoulder_landmarks(0.0);
        assert!((sh.canal_width() - 12.0).abs() < 1e-12);
        // The neck vector is 10 mm at 45° to the axis: axial height and
        // in-plane offset are both 10·sin45°.
        let in_plane = 10.0 * std::f64::consts::FRAC_1_SQRT_2;
        assert!((sh.head_height() - in_plane).abs() < 1e-12);
        assert!((sh.head_offset() - in_plane).abs() < 1e-12);
        assert!((sh.neck_shaft_angle_deg() - 135.0).abs() < 1e-9);
        assert!((sh.retroversion_deg() - 0.0).abs() < 1e-9);
        // The canal axis is unit and distal.
        assert!((sh.canal_axis().norm() - 1.0).abs() < 1e-12);
        assert!(sh.canal_axis().z > 0.99);

        // Rotating the in-plane neck component 25° off the TEA reads as
        // 25° of retroversion, without touching any other measurement.
        let rv = shoulder_landmarks(25.0);
        assert!((rv.retroversion_deg() - 25.0).abs() < 1e-9);
        assert!((rv.neck_shaft_angle_deg() - 135.0).abs() < 1e-9);
        assert!((rv.head_height() - sh.head_height()).abs() < 1e-9);
        assert!((rv.head_offset() - sh.head_offset()).abs() < 1e-9);
    }

    #[test]
    fn shoulder_sizing_resolves_by_precedence_like_the_hip() {
        let chart = SizeChart {
            family: "shoulder-stem-test".into(),
            entries: (1..=4)
                .map(|i| SizeEntry {
                    label: i,
                    nominal: 8.0 + 2.0 * i as f64, // 10, 12, 14, 16 mm
                })
                .collect(),
        };
        chart.validate().expect("chart valid");

        // Consensus: canal width exactly 14 (entry 3), and a neck long
        // enough that head offset AND head height both read 14.5 mm
        // (t = 0.25 toward entry 4 → label 3).
        let s = 45.0f64.to_radians().sin();
        let rv = 20.0f64.to_radians();
        let neck_axis = Vec3::new(s * rv.cos(), s * rv.sin(), -s);
        let wide = ShoulderLandmarks {
            head_center: Vec3::new(0.0, 0.0, 2.0) + neck_axis * (14.5 / s),
            anatomical_neck_center: Vec3::new(0.0, 0.0, 2.0),
            isthmus_medial: Vec3::new(-7.0, 0.0, 12.0),
            isthmus_lateral: Vec3::new(7.0, 0.0, 12.0),
            canal_entry: Vec3::new(0.0, 0.0, 0.0),
            medial_epicondyle: Vec3::new(-20.0, 0.0, -10.0),
            lateral_epicondyle: Vec3::new(20.0, 0.0, -10.0),
        };
        assert!((wide.head_offset() - 14.5).abs() < 1e-9);
        assert!((wide.head_height() - 14.5).abs() < 1e-9);
        let sizing = size_shoulder(&chart, &wide).expect("sizes");
        assert!(sizing.stem.is_unanimous());
        assert_eq!(sizing.stem.label, 3);
        assert!((sizing.neck_shaft_angle_deg - 135.0).abs() < 1e-9);
        assert!((sizing.retroversion_deg - 20.0).abs() < 1e-9);

        // Disagreement: narrowing the canal to 10.5 mm votes label 1
        // while offset and height still vote 3 — the canal width wins by
        // precedence.
        let narrow = ShoulderLandmarks {
            isthmus_medial: Vec3::new(-5.25, 0.0, 12.0),
            isthmus_lateral: Vec3::new(5.25, 0.0, 12.0),
            ..wide
        };
        let sizing = size_shoulder(&chart, &narrow).expect("sizes");
        assert_eq!(sizing.stem.label, 1);
        assert!(matches!(
            sizing.stem.resolved_by,
            Resolution::Precedence("canal_width")
        ));
        assert_eq!(sizing.stem.votes.len(), 3);
    }

    #[test]
    fn shoulder_sizing_is_monotone_in_canal_width() {
        let chart = SizeChart {
            family: "shoulder-stem-test".into(),
            entries: (1..=4)
                .map(|i| SizeEntry {
                    label: i,
                    nominal: 8.0 + 2.0 * i as f64,
                })
                .collect(),
        };
        let base = shoulder_landmarks(20.0);
        let mut previous = 0;
        for width in [10.5, 11.5, 12.5, 13.5, 14.5, 15.5] {
            let half = width / 2.0;
            let sh = ShoulderLandmarks {
                isthmus_medial: Vec3::new(-half, 0.0, 12.0),
                isthmus_lateral: Vec3::new(half, 0.0, 12.0),
                ..base
            };
            let label = size_shoulder(&chart, &sh).expect("sizes").stem.label;
            assert!(label >= previous, "width {width}: {label} after {previous}");
            previous = label;
        }
    }

    /// Plafond 30 mm wide (±15 in x) and 24 mm deep (±12 in y) at z = −30;
    /// dome 28 mm wide at z = −32; tibial axis point placed so the axis
    /// tilts exactly `tilt_deg` from the dome normal (30° tilt with a
    /// 50 mm axis arm: `sin 30° = 0.5`).
    fn ankle_landmarks(tilt_deg: f64) -> AnkleLandmarks {
        let arm = 50.0f64;
        let tilt = tilt_deg.to_radians();
        // Axis direction (distal) tilted `tilt` into −x from −z.
        let axis = Vec3::new(-tilt.sin(), 0.0, -tilt.cos());
        AnkleLandmarks {
            plafond_medial: Vec3::new(-15.0, 0.0, -30.0),
            plafond_lateral: Vec3::new(15.0, 0.0, -30.0),
            plafond_anterior: Vec3::new(0.0, 12.0, -30.0),
            plafond_posterior: Vec3::new(0.0, -12.0, -30.0),
            talar_medial: Vec3::new(-14.0, 0.0, -32.0),
            talar_lateral: Vec3::new(14.0, 0.0, -32.0),
            talar_center: Vec3::new(0.0, 0.0, -32.0),
            tibial_axis_point: Vec3::new(0.0, 0.0, -30.0) - axis * arm,
        }
    }

    #[test]
    fn ankle_measurements_are_exact_geometry() {
        let neutral = ankle_landmarks(0.0);
        assert!((neutral.plafond_width() - 30.0).abs() < 1e-12);
        assert!((neutral.plafond_depth() - 24.0).abs() < 1e-12);
        assert!((neutral.talar_width() - 28.0).abs() < 1e-12);
        // Axis perpendicular to the dome: zero deviation.
        assert!(neutral.tibiotalar_angle_deg().abs() < 1e-12);

        // A 30° axis tilt reads exactly 30° of deviation (sin 30° = 0.5).
        let tilted = ankle_landmarks(30.0);
        assert!((tilted.tibiotalar_angle_deg() - 30.0).abs() < 1e-9);
        // The width measurements are tilt-independent.
        assert!((tilted.plafond_width() - 30.0).abs() < 1e-12);
        assert!((tilted.talar_width() - 28.0).abs() < 1e-12);
    }

    #[test]
    fn ankle_sizing_bundles_both_components() {
        let tibial_chart = SizeChart {
            family: "taa-tibial-test".into(),
            entries: (1..=5)
                .map(|i| SizeEntry {
                    label: i,
                    nominal: 21.0 + 5.0 * i as f64, // 26, 31, 36, 41, 46 mm
                })
                .collect(),
        };
        let talar_chart = SizeChart {
            family: "taa-talar-test".into(),
            entries: (1..=3)
                .map(|i| SizeEntry {
                    label: i,
                    nominal: 20.0 + 4.0 * i as f64, // 24, 28, 32 mm
                })
                .collect(),
        };
        tibial_chart.validate().expect("tibial chart valid");
        talar_chart.validate().expect("talar chart valid");

        let sizing = size_ankle(&tibial_chart, &talar_chart, &ankle_landmarks(0.0)).expect("sizes");
        // Plafond width 30 → between 26 and 31, t = 0.8 → label 2;
        // depth 24 falls below the chart and votes label 1 — width wins
        // by precedence.
        assert_eq!(sizing.tibial.label, 2);
        assert!(matches!(
            sizing.tibial.resolved_by,
            Resolution::Precedence("plafond_width")
        ));
        assert!(sizing.tibial.votes.iter().any(|v| v.below_chart));
        // Dome width 28 sits exactly on talar entry 2: consensus (a
        // single vote is unanimous by construction).
        assert!(sizing.talar.is_unanimous());
        assert_eq!(sizing.talar.label, 2);
        assert!((sizing.plafond_width_mm - 30.0).abs() < 1e-12);
        assert!((sizing.tibiotalar_angle_deg - 0.0).abs() < 1e-12);
    }

    #[test]
    fn ankle_sizing_is_monotone_in_plafond_width() {
        let tibial_chart = SizeChart {
            family: "taa-tibial-test".into(),
            entries: (1..=5)
                .map(|i| SizeEntry {
                    label: i,
                    nominal: 21.0 + 5.0 * i as f64,
                })
                .collect(),
        };
        let talar_chart = tibial_chart.clone();
        let base = ankle_landmarks(0.0);
        let mut previous = 0;
        for width in [26.0, 30.0, 34.0, 38.0, 42.0, 46.0] {
            let half = width / 2.0;
            let ankle = AnkleLandmarks {
                plafond_medial: Vec3::new(-half, 0.0, -30.0),
                plafond_lateral: Vec3::new(half, 0.0, -30.0),
                ..base
            };
            let label = size_ankle(&tibial_chart, &talar_chart, &ankle)
                .expect("sizes")
                .tibial
                .label;
            assert!(label >= previous, "width {width}: {label} after {previous}");
            previous = label;
        }
    }

    fn collateral(side: CollateralSide, slack: f64, stiffness: f64) -> LigamentModel {
        LigamentModel {
            origin: Vec3::new(0.0, 0.0, 0.0),
            insertion: Vec3::new(0.0, 0.0, -40.0),
            slack_length_mm: slack,
            stiffness_n_per_mm: stiffness,
            side,
        }
    }

    #[test]
    fn ligament_tension_is_zero_below_slack_and_linear_above() {
        let lig = collateral(CollateralSide::Medial, 40.0, 25.0);
        assert!((lig.anatomical_length() - 40.0).abs() < 1e-12);
        // Slack range: no tension however unloaded.
        assert_eq!(lig.tension(40.0), 0.0);
        assert_eq!(lig.tension(35.0), 0.0);
        assert!(!lig.is_taut(40.0));
        // Linear beyond: 2 mm of elongation under 25 N/mm → 50 N.
        assert!((lig.tension(42.0) - 50.0).abs() < 1e-12);
        assert!(lig.is_taut(42.0));
    }

    #[test]
    fn planned_length_follows_moved_attachments() {
        let lig = collateral(CollateralSide::Lateral, 40.0, 25.0);
        // Lowering the tibial side by the resection + component offset
        // (4 mm distal along −z) lengthens the collateral to 44.
        assert!((lig.planned_length(lig.origin, Vec3::new(0.0, 0.0, -44.0)) - 44.0).abs() < 1e-12);
        // A distalized origin shortens it instead.
        assert!((lig.planned_length(Vec3::new(0.0, 0.0, 2.0), lig.insertion) - 42.0).abs() < 1e-12);
    }

    #[test]
    fn mediolateral_imbalance_sums_both_sides() {
        let medial = collateral(CollateralSide::Medial, 40.0, 25.0);
        let lateral = collateral(CollateralSide::Lateral, 40.0, 25.0);
        // A tight medial side (4 mm lengthening → 100 N) against a
        // relaxed lateral (1 mm → 25 N): imbalance 75 N.
        let report = assess_ligament_balance(
            &[
                (medial.clone(), medial.origin, Vec3::new(0.0, 0.0, -44.0)),
                (lateral.clone(), lateral.origin, Vec3::new(0.0, 0.0, -41.0)),
            ],
            20.0,
        )
        .expect("non-empty");
        assert!((report.ligaments[0].tension_n - 100.0).abs() < 1e-12);
        assert!((report.ligaments[1].tension_n - 25.0).abs() < 1e-12);
        assert!((report.mediolateral_imbalance_n - 75.0).abs() < 1e-12);
        assert!(!report.is_balanced);
        // Within a looser tolerance it balances numerically — and the
        // per-ligament flags still show both sides taut.
        let loose = assess_ligament_balance(
            &[
                (
                    collateral(CollateralSide::Medial, 40.0, 25.0),
                    Vec3::ZERO,
                    Vec3::new(0.0, 0.0, -44.0),
                ),
                (
                    collateral(CollateralSide::Lateral, 40.0, 25.0),
                    Vec3::ZERO,
                    Vec3::new(0.0, 0.0, -41.0),
                ),
            ],
            80.0,
        )
        .expect("non-empty");
        assert!(loose.is_balanced);
        assert!(loose.ligaments.iter().all(|l| l.is_taut));
    }

    #[test]
    fn an_over_released_side_balances_numerically_but_reads_slack() {
        // Over-release the lateral collateral past its slack range: it
        // carries no tension, so the imbalance vanishes — the trap the
        // docs warn about, asserted here.
        let report = assess_ligament_balance(
            &[
                (
                    collateral(CollateralSide::Medial, 40.0, 25.0),
                    Vec3::ZERO,
                    Vec3::new(0.0, 0.0, -41.0),
                ),
                (
                    collateral(CollateralSide::Lateral, 40.0, 25.0),
                    Vec3::ZERO,
                    Vec3::new(0.0, 0.0, -35.0),
                ),
            ],
            5.0,
        )
        .expect("non-empty");
        assert!((report.mediolateral_imbalance_n - 25.0).abs() < 1e-12);
        assert!(!report.is_balanced);
        let released = assess_ligament_balance(
            &[
                (
                    collateral(CollateralSide::Medial, 40.0, 25.0),
                    Vec3::ZERO,
                    Vec3::new(0.0, 0.0, -41.0),
                ),
                (
                    collateral(CollateralSide::Lateral, 40.0, 25.0),
                    Vec3::ZERO,
                    Vec3::new(0.0, 0.0, -30.0),
                ),
            ],
            5.0,
        )
        .expect("non-empty");
        assert_eq!(released.ligaments[1].tension_n, 0.0);
        assert!(
            !released.ligaments[1].is_taut,
            "over-released collateral is slack"
        );
    }

    #[test]
    fn empty_input_is_rejected() {
        assert!(assess_ligament_balance(&[], 1.0).is_none());
        assert!(assess_ligament_balance(&[], -1.0).is_none());
    }

    #[test]
    fn moment_arm_is_the_point_line_distance() {
        let lig = collateral(CollateralSide::Medial, 40.0, 25.0);
        // Ligament line: the z-axis through (0, 0). A joint centre 5 mm
        // off-axis has arm exactly 5; one ON the line has arm 0.
        assert!(
            (StabilityScreen::moment_arm_mm(
                &lig,
                Vec3::new(5.0, 0.0, 0.0),
                lig.origin,
                lig.insertion
            ) - 5.0)
                .abs()
                < 1e-12
        );
        assert!(
            (StabilityScreen::moment_arm_mm(
                &lig,
                Vec3::new(0.0, 0.0, -20.0),
                lig.origin,
                lig.insertion
            ))
            .abs()
                < 1e-12
        );
        // A slanted ligament: the arm is the point-line distance to the
        // line of action. Line direction (0, 3, −4)/5 through the origin
        // has constant x = 0, so a centre at (6, 0, 0) sits exactly 6 mm
        // off the line — past the insertion's projection, which is fine:
        // the arm follows the line of action, not the segment.
        let slant = LigamentModel {
            origin: Vec3::new(0.0, 0.0, 0.0),
            insertion: Vec3::new(0.0, 30.0, -40.0),
            ..lig
        };
        let arm = StabilityScreen::moment_arm_mm(
            &slant,
            Vec3::new(6.0, 0.0, 0.0),
            slant.origin,
            slant.insertion,
        );
        assert!((arm - 6.0).abs() < 1e-12, "{arm}");
    }

    #[test]
    fn net_moment_sums_signed_arm_tension_products() {
        // Symmetric collaterals 5 mm either side of the joint centre,
        // both taut at equal tension: the moments cancel exactly.
        let medial = LigamentModel {
            origin: Vec3::new(0.0, 5.0, 0.0),
            insertion: Vec3::new(0.0, 5.0, -40.0),
            slack_length_mm: 40.0,
            stiffness_n_per_mm: 25.0,
            side: CollateralSide::Medial,
        };
        let lateral = LigamentModel {
            origin: Vec3::new(0.0, -5.0, 0.0),
            insertion: Vec3::new(0.0, -5.0, -40.0),
            slack_length_mm: 40.0,
            stiffness_n_per_mm: 25.0,
            side: CollateralSide::Lateral,
        };
        let center = Vec3::new(0.0, 0.0, -20.0);
        let balanced = &[
            (medial.clone(), medial.origin, medial.insertion),
            (lateral.clone(), lateral.origin, lateral.insertion),
        ];
        let report = StabilityScreen::assess(balanced, center).expect("non-empty");
        assert!(
            report
                .moment_arms_mm
                .iter()
                .all(|a| (a - 5.0).abs() < 1e-12),
            "arms {:?}",
            report.moment_arms_mm
        );
        assert!(report.net_moment_nmm.abs() < 1e-12);

        // A taut medial against a slack lateral: the net moment is the
        // medial contribution alone (1 mm elongation → 25 N × 5 mm).
        let stretched = Vec3::new(0.0, 5.0, -41.0);
        let report = StabilityScreen::assess(
            &[
                (medial.clone(), medial.origin, stretched),
                (lateral.clone(), lateral.origin, lateral.insertion),
            ],
            center,
        )
        .expect("non-empty");
        assert!(
            (report.net_moment_nmm - 125.0).abs() < 1e-9,
            "{}",
            report.net_moment_nmm
        );
    }

    #[test]
    fn stability_screen_rejects_empty_input() {
        assert!(StabilityScreen::assess(&[], Vec3::ZERO).is_none());
        assert!(StabilityScreen::net_moment_nmm(&[], Vec3::ZERO).is_none());
    }
}
