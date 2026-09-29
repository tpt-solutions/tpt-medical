//! Osteotomy cuts and virtual surgery on voxel anatomical models.
//!
//! A virtual surgery is an ordered list of [`PlanStep`]s applied to a
//! voxel model: plane-based [`OsteotomyCut`]s split and discard/keep
//! fragments, [`FragmentTransform`]s reposition them (e.g. a tibial
//! tubercle distalization or a Le Fort advancement), and the final
//! `VirtualSurgery` produces the operated model plus an audit trail of
//! every step (consumed by `tpt-med-fda` in regulated pipelines).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::collections::BTreeMap;

use tpt_med_geometry::{Plane, Vec3};

/// A voxel scalar model with geometry (the surgical substrate).
#[derive(Debug, Clone)]
pub struct VoxelModel {
    /// Volume dimensions (nx, ny, nz).
    pub dims: (usize, usize, usize),
    /// Voxel pitch (mm) along (x, y, z).
    pub spacing: (f64, f64, f64),
    /// Patient-space position of voxel (0,0,0).
    pub origin: Vec3,
    /// Scalar values per voxel (HU or mask label; index `(z·ny+y)·nx+x`).
    pub values: Vec<f64>,
}

impl VoxelModel {
    /// Flat index of a voxel (None out of bounds).
    pub fn index(&self, x: usize, y: usize, z: usize) -> Option<usize> {
        let (nx, ny, nz) = self.dims;
        if x < nx && y < ny && z < nz {
            Some((z * ny + y) * nx + x)
        } else {
            None
        }
    }

    /// Patient-space center of a voxel.
    pub fn center(&self, x: usize, y: usize, z: usize) -> Vec3 {
        self.origin
            + Vec3::new(
                x as f64 * self.spacing.0,
                y as f64 * self.spacing.1,
                z as f64 * self.spacing.2,
            )
    }

    /// Count of voxels above `threshold`.
    pub fn count_above(&self, threshold: f64) -> usize {
        self.values.iter().filter(|&&v| v >= threshold).count()
    }
}

/// A plane cut: keeps the voxels on the `keep_side` side of `plane`
/// (`signed_distance ≥ 0`), with an optional **saw-kerf width** — material
/// within `kerf_width / 2` of the plane is removed on both sides, so a
/// kerf of 4 mm shifts the kept boundary 2 mm outward and the resection
/// measurement grows by the slab.
#[derive(Debug, Clone)]
pub struct OsteotomyCut {
    /// Cutting plane.
    pub plane: Plane,
    /// Fragment name produced from the kept side (audit label).
    pub fragment_name: String,
    /// Keep the positive (normal-side) region; if false keep the negative.
    pub keep_positive: bool,
    /// Saw-kerf width (mm): material within `kerf_width / 2` of the plane
    /// is discarded from both sides. 0.0 (the default) is a pure plane cut.
    pub kerf_width: f64,
}

impl OsteotomyCut {
    /// Applies the cut to a model; returns the kept fragment. Voxels on the
    /// discarded side are set to `f64::NAN` and compacted by `bounding_box`
    /// of the kept region.
    pub fn apply(&self, model: &VoxelModel) -> VoxelModel {
        self.apply_measured(model).0
    }

    /// [`Self::apply`] plus the step's measurements: the **resection
    /// volume** (discarded non-empty voxels × voxel volume) and the **cut
    /// depth** (deepest discarded voxel centre below the plane, 0 when
    /// nothing was discarded).
    pub fn apply_measured(&self, model: &VoxelModel) -> (VoxelModel, CutMeasurement) {
        let mut out = model.clone();
        let (nx, ny, nz) = model.dims;
        let voxel_volume = model.spacing.0 * model.spacing.1 * model.spacing.2;
        let mut min = [usize::MAX; 3];
        let mut max = [0usize; 3];
        let mut resected = 0usize;
        let mut max_depth = 0.0f64;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let Some(idx) = model.index(x, y, z) else {
                        continue;
                    };
                    if model.values[idx].is_nan() {
                        continue;
                    }
                    let c = model.center(x, y, z);
                    let d = self.plane.signed_distance(c);
                    let half_kerf = self.kerf_width * 0.5;
                    let keep = if self.keep_positive {
                        d >= half_kerf
                    } else {
                        d <= -half_kerf
                    };
                    if keep {
                        min[0] = min[0].min(x);
                        min[1] = min[1].min(y);
                        min[2] = min[2].min(z);
                        max[0] = max[0].max(x);
                        max[1] = max[1].max(y);
                        max[2] = max[2].max(z);
                    } else {
                        out.values[idx] = f64::NAN;
                        resected += 1;
                        // Depth past the kerf face on the kept side, not
                        // past the plane (with kerf 0 this is the distance
                        // past the plane on the discard side).
                        let depth = if self.keep_positive {
                            half_kerf - d
                        } else {
                            d + half_kerf
                        };
                        max_depth = max_depth.max(depth);
                    }
                }
            }
        }
        let measurement = CutMeasurement {
            fragment_name: self.fragment_name.clone(),
            resection_volume_mm3: resected as f64 * voxel_volume,
            max_depth_mm: max_depth,
        };
        if max[0] == 0 {
            return (out, measurement); // empty cut result; leave as-is
        }
        compact(&mut out, min, max);
        (out, measurement)
    }
}

/// Measurements of one [`OsteotomyCut`], for the surgical report.
#[derive(Debug, Clone)]
pub struct CutMeasurement {
    /// The cut's fragment label.
    pub fragment_name: String,
    /// Discarded tissue volume (mm³).
    pub resection_volume_mm3: f64,
    /// Deepest discarded voxel centre below the cutting plane (mm).
    pub max_depth_mm: f64,
}

/// Crops the model to the voxel box `[min, max]` (inclusive).
fn compact(model: &mut VoxelModel, min: [usize; 3], max: [usize; 3]) {
    let (nx, ny, _) = model.dims;
    let (sx, sy, sz) = model.spacing;
    let new_dims = (
        max[0] - min[0] + 1,
        max[1] - min[1] + 1,
        max[2] - min[2] + 1,
    );
    let mut values = Vec::with_capacity(new_dims.0 * new_dims.1 * new_dims.2);
    for z in min[2]..=max[2] {
        for y in min[1]..=max[1] {
            for x in min[0]..=max[0] {
                values.push(model.values[(z * ny + y) * nx + x]);
            }
        }
    }
    model.dims = new_dims;
    model.values = values;
    model.origin += Vec3::new(min[0] as f64 * sx, min[1] as f64 * sy, min[2] as f64 * sz);
}

/// A **multi-plane wedge**: the closed wedge removed by two intersecting
/// osteotomy planes — material on the discarded side of *both* planes
/// (`d_a ≤ 0 ∧ d_b ≤ 0`, each offset by a half-kerf) is removed, and
/// everything else is kept. Two sequential [`OsteotomyCut`]s cannot express
/// this (each keeps only one side of one plane), which is why the wedge is
/// its own step.
#[derive(Debug, Clone)]
pub struct WedgeCut {
    /// First wedge plane (its negative side is removed where it overlaps
    /// the second).
    pub plane_a: Plane,
    /// Second wedge plane.
    pub plane_b: Plane,
    /// Fragment name for the kept material (audit label).
    pub fragment_name: String,
    /// Saw-kerf width (mm), offsetting both planes outward into the kept
    /// region symmetrically.
    pub kerf_width: f64,
}

impl WedgeCut {
    /// Applies the wedge; returns the kept model plus the step's
    /// measurements (resection volume; cut depth measured to the nearer of
    /// the two planes).
    pub fn apply_measured(&self, model: &VoxelModel) -> (VoxelModel, CutMeasurement) {
        let mut out = model.clone();
        let (nx, ny, nz) = model.dims;
        let voxel_volume = model.spacing.0 * model.spacing.1 * model.spacing.2;
        let half_kerf = self.kerf_width * 0.5;
        let mut min = [usize::MAX; 3];
        let mut max = [0usize; 3];
        let mut resected = 0usize;
        let mut max_depth = 0.0f64;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let Some(idx) = model.index(x, y, z) else {
                        continue;
                    };
                    if model.values[idx].is_nan() {
                        continue;
                    }
                    let c = model.center(x, y, z);
                    let da = self.plane_a.signed_distance(c);
                    let db = self.plane_b.signed_distance(c);
                    // Keep whenever either plane's kept side is reached.
                    let keep = da >= half_kerf || db >= half_kerf;
                    if keep {
                        min[0] = min[0].min(x);
                        min[1] = min[1].min(y);
                        min[2] = min[2].min(z);
                        max[0] = max[0].max(x);
                        max[1] = max[1].max(y);
                        max[2] = max[2].max(z);
                    } else {
                        out.values[idx] = f64::NAN;
                        resected += 1;
                        // Depth into the wedge = distance past the nearer
                        // kept face.
                        max_depth = max_depth.max(half_kerf - da.max(db));
                    }
                }
            }
        }
        let measurement = CutMeasurement {
            fragment_name: self.fragment_name.clone(),
            resection_volume_mm3: resected as f64 * voxel_volume,
            max_depth_mm: max_depth,
        };
        if max[0] == 0 {
            return (out, measurement);
        }
        compact(&mut out, min, max);
        (out, measurement)
    }

    /// Applies the wedge without measurements.
    pub fn apply(&self, model: &VoxelModel) -> VoxelModel {
        self.apply_measured(model).0
    }
}

/// A rigid fragment transform: rotation (axis-angle) about a pivot plus/// A rigid fragment transform: rotation (axis-angle) about a pivot plus
/// translation (both in patient coordinates, mm / radians).
#[derive(Debug, Clone)]
pub struct FragmentTransform {
    /// Rotation axis (normalized internally).
    pub rotation_axis: Vec3,
    /// Rotation angle (radians, right-hand rule).
    pub rotation_angle: f64,
    /// Rotation pivot point (patient coordinates).
    pub pivot: Vec3,
    /// Translation after rotation (mm).
    pub translation: Vec3,
}

impl FragmentTransform {
    /// Patient-space position of a point after the transform.
    pub fn apply_to_point(&self, p: Vec3) -> Vec3 {
        let r =
            tpt_med_geometry::Mat3::rotation_axis_angle(self.rotation_axis, self.rotation_angle);
        r * (p - self.pivot) + self.pivot + self.translation
    }

    /// Applies the transform to all non-NaN voxels of a fragment model.
    /// The output grid is grown to cover the transformed extent so no
    /// voxels are lost (nearest-neighbour scatter by voxel centre).
    pub fn apply_to_model(&self, fragment: &VoxelModel) -> VoxelModel {
        let (nx, ny, nz) = fragment.dims;
        // Pass 1: destination index-space bounds.
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let Some(idx) = fragment.index(x, y, z) else {
                        continue;
                    };
                    if fragment.values[idx].is_nan() {
                        continue;
                    }
                    let dst = self.apply_to_point(fragment.center(x, y, z));
                    let rel = dst - fragment.origin;
                    let g = [
                        rel.x / fragment.spacing.0,
                        rel.y / fragment.spacing.1,
                        rel.z / fragment.spacing.2,
                    ];
                    for a in 0..3 {
                        lo[a] = lo[a].min(g[a]);
                        hi[a] = hi[a].max(g[a]);
                    }
                }
            }
        }
        if hi[0] < lo[0] {
            return fragment.clone(); // nothing to move
        }
        let min = [
            lo[0].floor() as i64,
            lo[1].floor() as i64,
            lo[2].floor() as i64,
        ];
        let max = [
            hi[0].ceil() as i64,
            hi[1].ceil() as i64,
            hi[2].ceil() as i64,
        ];
        let new_dims = (
            (max[0] - min[0] + 1) as usize,
            (max[1] - min[1] + 1) as usize,
            (max[2] - min[2] + 1) as usize,
        );
        let mut new_values = vec![f64::NAN; new_dims.0 * new_dims.1 * new_dims.2];
        let to_local = |g: [i64; 3]| {
            (
                (g[0] - min[0]) as usize,
                (g[1] - min[1]) as usize,
                (g[2] - min[2]) as usize,
            )
        };
        // Pass 2: scatter.
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let Some(idx) = fragment.index(x, y, z) else {
                        continue;
                    };
                    if fragment.values[idx].is_nan() {
                        continue;
                    }
                    let dst = self.apply_to_point(fragment.center(x, y, z));
                    let rel = dst - fragment.origin;
                    let g = [
                        (rel.x / fragment.spacing.0).round() as i64,
                        (rel.y / fragment.spacing.1).round() as i64,
                        (rel.z / fragment.spacing.2).round() as i64,
                    ];
                    let (gx, gy, gz) = to_local(g);
                    let didx = (gz * new_dims.1 + gy) * new_dims.0 + gx;
                    new_values[didx] = fragment.values[idx];
                }
            }
        }
        VoxelModel {
            dims: new_dims,
            spacing: fragment.spacing,
            origin: fragment.origin
                + Vec3::new(
                    min[0] as f64 * fragment.spacing.0,
                    min[1] as f64 * fragment.spacing.1,
                    min[2] as f64 * fragment.spacing.2,
                ),
            values: new_values,
        }
    }
}

/// One recorded step in a virtual surgery plan.
#[derive(Debug, Clone)]
pub enum PlanStep {
    /// Plane osteotomy.
    Cut(OsteotomyCut),
    /// Two-plane closed-wedge osteotomy.
    Wedge(WedgeCut),
    /// Rigid fragment reposition.
    Move(FragmentTransform),
}

/// Measurement of one [`PlanStep::Move`]: what the transform was prescribed
/// to do to the fragment centroid vs what the voxel scatter achieved —
/// their difference is the **achieved alignment error** (the quantisation
/// of a sub-voxel plan onto the grid).
#[derive(Debug, Clone)]
pub struct MoveMeasurement {
    /// Centroid displacement the prescribed transform implies (mm).
    pub prescribed_centroid_translation: Vec3,
    /// Centroid displacement actually realised by the scatter (mm).
    pub achieved_centroid_translation: Vec3,
    /// `|prescribed − achieved|` (mm).
    pub alignment_error_mm: f64,
}

/// One measured step outcome, index-aligned with the step log.
#[derive(Debug, Clone)]
pub enum StepMeasurement {
    /// Plane osteotomy measurements.
    Cut(CutMeasurement),
    /// Fragment reposition measurements.
    Move(MoveMeasurement),
}

/// Measurements recorded alongside the audit log by
/// [`VirtualSurgery::execute_with_report`].
#[derive(Debug, Clone)]
pub struct SurgeryReport {
    /// Step descriptions in execution order (the audit log).
    pub step_log: Vec<String>,
    /// Per-step measurements, same order and length as `step_log`.
    pub measurements: Vec<StepMeasurement>,
    /// Summed resection volume over all cut steps (mm³).
    pub total_resection_volume_mm3: f64,
}

/// A virtual surgery plan over a base model.
#[derive(Debug, Clone)]
pub struct VirtualSurgery {
    base: VoxelModel,
    steps: Vec<PlanStep>,
    /// Labels of fragments produced/modified, in order.
    fragment_log: Vec<String>,
}

impl VirtualSurgery {
    /// Starts a plan on the base (pre-operative) model.
    pub fn new(base: VoxelModel) -> Self {
        Self {
            base,
            steps: Vec::new(),
            fragment_log: Vec::new(),
        }
    }

    /// Appends an osteotomy step.
    pub fn cut(&mut self, cut: OsteotomyCut) -> &mut Self {
        self.fragment_log.push(cut.fragment_name.clone());
        self.steps.push(PlanStep::Cut(cut));
        self
    }

    /// Appends a two-plane closed-wedge step.
    pub fn wedge(&mut self, cut: WedgeCut) -> &mut Self {
        self.fragment_log
            .push(format!("wedge:{}", cut.fragment_name));
        self.steps.push(PlanStep::Wedge(cut));
        self
    }

    /// Appends a fragment transform step.
    pub fn move_fragment(&mut self, transform: FragmentTransform) -> &mut Self {
        self.fragment_log.push(format!(
            "rotate({:.3} rad) + translate({:.1}, {:.1}, {:.1}) mm",
            transform.rotation_angle,
            transform.translation.x,
            transform.translation.y,
            transform.translation.z
        ));
        self.steps.push(PlanStep::Move(transform));
        self
    }

    /// Executes the plan, returning the operated model and the audit log
    /// (step descriptions in execution order).
    pub fn execute(&self) -> (VoxelModel, Vec<String>) {
        let (model, report) = self.execute_with_report();
        (model, report.step_log)
    }

    /// [`Self::execute`] with the measurement report: resection volumes and
    /// cut depths per cut, achieved-vs-prescribed centroid displacement per
    /// move (the alignment error), and the total resection volume. The
    /// report is index-aligned with the audit log, so a submission bundle
    /// can attach the numbers to the steps they belong to.
    pub fn execute_with_report(&self) -> (VoxelModel, SurgeryReport) {
        let mut model = self.base.clone();
        let mut log = Vec::with_capacity(self.steps.len());
        let mut measurements = Vec::with_capacity(self.steps.len());
        let mut total_resection = 0.0f64;
        for step in &self.steps {
            match step {
                PlanStep::Cut(cut) => {
                    let (next, m) = cut.apply_measured(&model);
                    total_resection += m.resection_volume_mm3;
                    log.push(format!("cut:{}", cut.fragment_name));
                    measurements.push(StepMeasurement::Cut(m));
                    model = next;
                }
                PlanStep::Wedge(cut) => {
                    let (next, m) = cut.apply_measured(&model);
                    total_resection += m.resection_volume_mm3;
                    log.push(format!("wedge:{}", cut.fragment_name));
                    measurements.push(StepMeasurement::Cut(m));
                    model = next;
                }
                PlanStep::Move(m) => {
                    let before = fragment_centroid(&model);
                    let next = m.apply_to_model(&model);
                    let after = fragment_centroid(&next);
                    let achieved = after - before;
                    let prescribed = m.apply_to_point(before) - before;
                    log.push(format!("move:{}", format_args!("{:?}", m.translation)));
                    measurements.push(StepMeasurement::Move(MoveMeasurement {
                        prescribed_centroid_translation: prescribed,
                        achieved_centroid_translation: achieved,
                        alignment_error_mm: (prescribed - achieved).norm(),
                    }));
                    model = next;
                }
            }
        }
        (
            model,
            SurgeryReport {
                step_log: log,
                measurements,
                total_resection_volume_mm3: total_resection,
            },
        )
    }

    /// The pre-operative base model (for side-by-side planning views).
    pub fn base_model(&self) -> &VoxelModel {
        &self.base
    }

    /// Fragment labels recorded so far, deduplicated and ordered.
    pub fn fragments(&self) -> BTreeMap<usize, String> {
        self.fragment_log
            .iter()
            .enumerate()
            .map(|(i, s)| (i, s.clone()))
            .collect()
    }
}

/// Centroid of the non-empty (non-NaN) voxels of a model, in patient
/// coordinates. `Vec3::ZERO` for an empty model.
fn fragment_centroid(model: &VoxelModel) -> Vec3 {
    let (nx, ny, nz) = model.dims;
    let mut sum = Vec3::ZERO;
    let mut count = 0usize;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let Some(idx) = model.index(x, y, z) else {
                    continue;
                };
                if model.values[idx].is_nan() {
                    continue;
                }
                sum += model.center(x, y, z);
                count += 1;
            }
        }
    }
    if count == 0 {
        Vec3::ZERO
    } else {
        sum / count as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 10×10×10 model, values 100 inside a centred 6×6×6 cube.
    fn cube_model() -> VoxelModel {
        let mut m = VoxelModel {
            dims: (10, 10, 10),
            spacing: (1.0, 1.0, 1.0),
            origin: Vec3::new(-5.0, -5.0, -5.0),
            values: vec![0.0; 1000],
        };
        for z in 2..8 {
            for y in 2..8 {
                for x in 2..8 {
                    let i = m.index(x, y, z).unwrap();
                    m.values[i] = 100.0;
                }
            }
        }
        m
    }

    #[test]
    fn plane_cut_keeps_positive_side() {
        let model = cube_model();
        let plane = Plane::from_point_normal(Vec3::new(0.0, 0.0, 0.0), Vec3::Z).unwrap();
        let cut = OsteotomyCut {
            plane,
            fragment_name: "proximal".into(),
            keep_positive: true,
            kerf_width: 0.0,
        };
        let out = cut.apply(&model);
        // Positive-z side of the cube: z ∈ 0..=7 slices kept from the cube
        // region (centre z = -5+z ≥ 0 → z ≥ 5): z ∈ 5..8 → 3 slices of 36.
        assert_eq!(out.count_above(50.0), 3 * 36);
        // Grid compacts to the kept bounding box: z ∈ 5..=9 (cube slices
        // 5..7 plus zero-valued tissue slices 8..9 above the plane).
        assert_eq!(out.dims.2, 5);
    }

    #[test]
    fn cut_then_move_shifts_fragment() {
        let model = cube_model();
        let plane = Plane::from_point_normal(Vec3::new(0.0, 0.0, 0.0), Vec3::Z).unwrap();
        let mut plan = VirtualSurgery::new(model);
        plan.cut(OsteotomyCut {
            plane,
            fragment_name: "distal".into(),
            keep_positive: false,
            kerf_width: 0.0,
        });
        // Negative-z fragment: z ∈ 2..5 (3 slices). Distract +2 mm in z.
        plan.move_fragment(FragmentTransform {
            rotation_axis: Vec3::Z,
            rotation_angle: 0.0,
            pivot: Vec3::ZERO,
            translation: Vec3::new(0.0, 0.0, 2.0),
        });
        let (operated, log) = plan.execute();
        assert_eq!(log.len(), 2);
        assert!(log[0].starts_with("cut:distal"));
        // Negative-z fragment keeps cube slices z ∈ 2..=5 (4 × 36 voxels);
        // the transform's output grid grows, so none are lost.
        assert_eq!(
            operated.count_above(50.0),
            4 * 36,
            "voxel count preserved by translation"
        );
    }

    #[test]
    fn rotation_about_z_maps_points() {
        let t = FragmentTransform {
            rotation_axis: Vec3::Z,
            rotation_angle: core::f64::consts::FRAC_PI_2,
            pivot: Vec3::ZERO,
            translation: Vec3::ZERO,
        };
        let q = t.apply_to_point(Vec3::new(1.0, 0.0, 0.0));
        assert!((q.x - 0.0).abs() < 1e-9);
        assert!((q.y - 1.0).abs() < 1e-9);
    }

    #[test]
    fn fragment_bookkeeping() {
        let mut plan = VirtualSurgery::new(cube_model());
        assert!(plan.fragments().is_empty());
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::X).unwrap(),
            fragment_name: "lateral".into(),
            keep_positive: true,
            kerf_width: 0.0,
        });
        assert_eq!(plan.fragments().len(), 1);
        assert_eq!(plan.base_model().dims, (10, 10, 10));
    }

    #[test]
    fn cut_measurements_report_volume_and_depth() {
        let model = cube_model();
        // Keep z ≥ 0: every voxel with centre below the plane is resected —
        // including zero-valued tissue (only NaN is empty) — 5 slices of
        // 100 voxels, centres down to z = −5 mm below the plane.
        let cut = OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
            fragment_name: "distal".into(),
            keep_positive: true,
            kerf_width: 0.0,
        };
        let (out, m) = cut.apply_measured(&model);
        assert_eq!(out.count_above(50.0), 3 * 36);
        assert_eq!(m.fragment_name, "distal");
        assert!(
            (m.resection_volume_mm3 - 5.0 * 100.0).abs() < 1e-9,
            "{}",
            m.resection_volume_mm3
        );
        assert!((m.max_depth_mm - 5.0).abs() < 1e-9, "{}", m.max_depth_mm);
        // `apply` is the un-measured form of the same operation.
        assert_eq!(cut.apply(&model).count_above(50.0), out.count_above(50.0));
    }

    #[test]
    fn saw_kerf_removes_a_symmetric_slab_and_grows_the_resection() {
        let model = cube_model();
        // Zero-kerf baseline: keep z ≥ 0.
        let mut base = OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
            fragment_name: "d".into(),
            keep_positive: true,
            kerf_width: 0.0,
        };
        let (kept0, m0) = base.apply_measured(&model);
        assert_eq!(kept0.count_above(50.0), 3 * 36);
        assert!((m0.max_depth_mm - 5.0).abs() < 1e-9);
        // A 2 mm kerf removes the slab (−1, +1) around the plane: the
        // kept boundary moves from centre 0 to centre +1, costing the
        // slice at centre 0 (36 voxels of cube). Depth is measured past
        // the kept face at +1 mm, so the deepest voxel (centre −5) sits
        // 6 mm from it.
        base.kerf_width = 2.0;
        let (kept2, m2) = base.apply_measured(&model);
        assert_eq!(kept2.count_above(50.0), 2 * 36);
        assert!((m2.max_depth_mm - 6.0).abs() < 1e-9, "{}", m2.max_depth_mm);
        // The slab is resected tissue: volume grows by the slice.
        assert!((m2.resection_volume_mm3 - m0.resection_volume_mm3 - 100.0).abs() < 1e-9);
        // The negative side with kerf: the slab is removed from there too.
        base.keep_positive = false;
        let (kept_neg, _) = base.apply_measured(&model);
        assert_eq!(kept_neg.count_above(50.0), 3 * 36);
    }

    #[test]
    fn wedge_removes_exactly_the_two_plane_intersection() {
        let model = cube_model();
        // Two planes through the origin, normals +x and +z: the wedge
        // region is {x < 0 ∧ z < 0} — voxel centres with x ≤ −1 and z ≤ −1.
        let wedge = WedgeCut {
            plane_a: Plane::from_point_normal(Vec3::ZERO, Vec3::X).unwrap(),
            plane_b: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
            fragment_name: "wedge".into(),
            kerf_width: 0.0,
        };
        let (out, m) = wedge.apply_measured(&model);
        // Discarded: 5 x-slices × 5 z-slices of 10 = 250 voxels.
        assert!(
            (m.resection_volume_mm3 - 250.0).abs() < 1e-9,
            "{}",
            m.resection_volume_mm3
        );
        // A voxel in the kept region on plane A's positive side but plane
        // B's negative side survives — the property two sequential cuts
        // cannot express.
        let survived = out.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(survived, 1000 - 250);
        // Depth into the wedge is measured to the nearer face.
        assert!(m.max_depth_mm > 0.0 && m.max_depth_mm <= 5.0 + 1e-9);
        // The wedge integrates through the plan with its own log label.
        let mut plan = VirtualSurgery::new(model);
        plan.wedge(wedge);
        let (_, report) = plan.execute_with_report();
        assert!(report.step_log[0].starts_with("wedge:wedge"));
        assert!((report.total_resection_volume_mm3 - 250.0).abs() < 1e-9);
    }

    #[test]
    fn wedge_kerf_offsets_both_planes() {
        let model = cube_model();
        let wedge = WedgeCut {
            plane_a: Plane::from_point_normal(Vec3::ZERO, Vec3::X).unwrap(),
            plane_b: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
            fragment_name: "w".into(),
            kerf_width: 2.0,
        };
        let (_, m) = wedge.apply_measured(&model);
        // Half-kerf offset: discard now includes centres with x ≤ 0 or
        // z ≤ 0 in the overlap — strictly more than the zero-kerf wedge.
        assert!(m.resection_volume_mm3 > 250.0);
    }

    #[test]
    fn move_measurements_report_alignment_error() {
        let model = cube_model();
        let mut plan = VirtualSurgery::new(model);
        // Grid-aligned translation: the scatter achieves it exactly.
        plan.move_fragment(FragmentTransform {
            rotation_axis: Vec3::Z,
            rotation_angle: 0.0,
            pivot: Vec3::ZERO,
            translation: Vec3::new(0.0, 0.0, 2.0),
        });
        // Sub-voxel translation: rounds to one voxel (1 mm) — the 0.3 mm
        // shortfall is the alignment error the report exists to surface.
        plan.move_fragment(FragmentTransform {
            rotation_axis: Vec3::Z,
            rotation_angle: 0.0,
            pivot: Vec3::ZERO,
            translation: Vec3::new(0.3, 0.0, 0.0),
        });
        let (operated, report) = plan.execute_with_report();
        assert_eq!(report.measurements.len(), 2);
        assert_eq!(report.step_log.len(), 2);
        match &report.measurements[0] {
            StepMeasurement::Move(m) => {
                assert!((m.prescribed_centroid_translation.z - 2.0).abs() < 1e-9);
                assert!((m.achieved_centroid_translation.z - 2.0).abs() < 1e-9);
                assert!(m.alignment_error_mm < 1e-9, "grid-aligned move is exact");
            }
            other => panic!("expected a move measurement, got {other:?}"),
        }
        match &report.measurements[1] {
            StepMeasurement::Move(m) => {
                assert!((m.prescribed_centroid_translation.x - 0.3).abs() < 1e-9);
                let achieved = m.achieved_centroid_translation.x;
                assert!(
                    achieved.abs() < 1e-9,
                    "a 0.3 mm shift rounds back onto the same voxel grid: {achieved}"
                );
                assert!(
                    (m.alignment_error_mm - 0.3).abs() < 1e-9,
                    "{}",
                    m.alignment_error_mm
                );
            }
            other => panic!("expected a move measurement, got {other:?}"),
        }
        // No cut in this plan: the full 6³ cube survives both moves.
        assert_eq!(operated.count_above(50.0), 216);
    }

    #[test]
    fn report_totals_resection_across_cuts() {
        let model = cube_model();
        let mut plan = VirtualSurgery::new(model);
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
            fragment_name: "a".into(),
            keep_positive: true,
            kerf_width: 0.0,
        });
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::X).unwrap(),
            fragment_name: "b".into(),
            keep_positive: false,
            kerf_width: 0.0,
        });
        let (_, report) = plan.execute_with_report();
        // First cut removes 5 slices of 100 (zero-valued tissue included);
        // the second keeps x ≤ 0 and discards the 4 positive-x columns of
        // the 10×10×5 compacted remainder (200 voxels).
        let cut_a = match &report.measurements[0] {
            StepMeasurement::Cut(m) => m.resection_volume_mm3,
            other => panic!("{other:?}"),
        };
        let cut_b = match &report.measurements[1] {
            StepMeasurement::Cut(m) => m.resection_volume_mm3,
            other => panic!("{other:?}"),
        };
        assert!((cut_a - 500.0).abs() < 1e-9, "{cut_a}");
        assert!((cut_b - 200.0).abs() < 1e-9, "{cut_b}");
        assert!((report.total_resection_volume_mm3 - (cut_a + cut_b)).abs() < 1e-9);
        // Measurements are index-aligned with the step log.
        assert!(report.step_log[0].starts_with("cut:a"));
        assert!(report.step_log[1].starts_with("cut:b"));
        assert!(matches!(report.measurements[1], StepMeasurement::Cut(_)));
    }
}
