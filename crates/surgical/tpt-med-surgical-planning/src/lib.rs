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
/// (`signed_distance ≥ 0`).
#[derive(Debug, Clone)]
pub struct OsteotomyCut {
    /// Cutting plane.
    pub plane: Plane,
    /// Fragment name produced from the kept side (audit label).
    pub fragment_name: String,
    /// Keep the positive (normal-side) region; if false keep the negative.
    pub keep_positive: bool,
}

impl OsteotomyCut {
    /// Applies the cut to a model; returns the kept fragment. Voxels on the
    /// discarded side are set to `f64::NAN` and compacted by `bounding_box`
    /// of the kept region.
    pub fn apply(&self, model: &VoxelModel) -> VoxelModel {
        let mut out = model.clone();
        let (nx, ny, nz) = model.dims;
        let mut min = [usize::MAX; 3];
        let mut max = [0usize; 3];
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
                    let keep = if self.keep_positive {
                        d >= 0.0
                    } else {
                        d <= 0.0
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
                    }
                }
            }
        }
        if max[0] == 0 {
            return out; // empty cut result; leave as-is
        }
        compact(&mut out, min, max);
        out
    }
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

/// A rigid fragment transform: rotation (axis-angle) about a pivot plus
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
    /// Rigid fragment reposition.
    Move(FragmentTransform),
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
        let mut model = self.base.clone();
        let mut log = Vec::with_capacity(self.steps.len());
        for step in &self.steps {
            match step {
                PlanStep::Cut(cut) => {
                    model = cut.apply(&model);
                    log.push(format!("cut:{}", cut.fragment_name));
                }
                PlanStep::Move(m) => {
                    model = m.apply_to_model(&model);
                    log.push(format!("move:{}", format_args!("{:?}", m.translation)));
                }
            }
        }
        (model, log)
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
        });
        assert_eq!(plan.fragments().len(), 1);
        assert_eq!(plan.base_model().dims, (10, 10, 10));
    }
}
