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
/// What happens to a cut's discarded side (`rfcs/0011`'s resolution):
/// resected — the executor's original semantics, and the default — or
/// **retained as a named fragment**, so the two pieces of an osteotomy
/// can be addressed and repositioned independently.
#[derive(Debug, Clone, PartialEq)]
pub enum DiscardedSide {
    /// The discarded side is removed from the model (original behaviour).
    Resect,
    /// The discarded side is kept as a named, independently addressable
    /// fragment.
    RetainAs {
        /// The retained fragment's name (audit label; must be non-empty).
        name: String,
    },
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
    /// The discarded side's fate. `Resect` (the default) is the original
    /// behaviour; `RetainAs` keeps it as a second named fragment.
    pub discarded: DiscardedSide,
}

impl OsteotomyCut {
    /// Applies the cut to a model; returns the kept fragment. Voxels on the
    /// discarded side are set to `f64::NAN` and compacted by `bounding_box`
    /// of the kept region.
    pub fn apply(&self, model: &VoxelModel) -> VoxelModel {
        self.apply_measured(model).0
    }

    /// [`Self::apply`] plus the step's measurements: the **resection
    /// volume** (resected non-empty voxels × voxel volume — a retained
    /// side is not resection) and the **cut depth** (deepest discarded
    /// voxel centre below the plane, 0 when nothing was discarded).
    /// Returns `(kept model, discarded model when retained, measurement)`.
    pub fn apply_measured(
        &self,
        model: &VoxelModel,
    ) -> (VoxelModel, Option<VoxelModel>, CutMeasurement) {
        let mut out = model.clone();
        let (nx, ny, nz) = model.dims;
        let voxel_volume = model.spacing.0 * model.spacing.1 * model.spacing.2;
        let retaining = matches!(self.discarded, DiscardedSide::RetainAs { .. });
        // The retained side keeps its own copy of the grid (a fragment is
        // an independent model; moves later diverge the grids).
        let mut retained = if retaining { Some(model.clone()) } else { None };
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
                        if let Some(r) = &mut retained {
                            r.values[idx] = f64::NAN;
                        }
                    } else {
                        out.values[idx] = f64::NAN;
                        if !retaining {
                            resected += 1;
                        }
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
            return (out, None, measurement); // empty cut result; leave as-is
        }
        compact(&mut out, min, max);
        let retained_model = retained.and_then(compact_retained);
        (out, retained_model, measurement)
    }
}

/// Measurements of one cut step (plane [`OsteotomyCut`], [`WedgeCut`] or
/// [`CylindricalCut`]), for the surgical report.
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

/// A **curved (cylindrical) resection**: the cylinder of radius `radius`
/// about an axis (`axis_origin`, `axis_direction`), keeping one side of
/// the wall — the core within `radius` (a reaming or core decompression)
/// or the annulus outside it. This is the first curved resection surface
/// in the crate: a cylinder is the shape a reamer or a burr actually
/// leaves, and the surface a rotational (derotation) osteotomy swings
/// about. An optional **saw-kerf width** removes the radial slab within
/// `kerf_width / 2` of the wall on both sides, exactly as
/// [`OsteotomyCut`]'s kerf does for a plane. The discarded side follows
/// the `rfcs/0011` resolution ([`DiscardedSide`]): resected by default,
/// or retained as a named fragment (e.g. a cylindrical core kept for
/// grafting).
#[derive(Debug, Clone)]
pub struct CylindricalCut {
    /// A point on the cylinder's axis (patient coordinates).
    pub axis_origin: Vec3,
    /// The cylinder's axis direction (normalized internally; a zero
    /// vector degrades the distance to "distance from `axis_origin`", so
    /// always pass a real axis).
    pub axis_direction: Vec3,
    /// Cylinder radius (mm).
    pub radius: f64,
    /// Fragment name produced from the kept side (audit label).
    pub fragment_name: String,
    /// Keep the material within `radius` of the axis (the core); `false`
    /// keeps the annulus outside and resects the core.
    pub keep_inside: bool,
    /// Saw-kerf width (mm): material within `kerf_width / 2` of the
    /// cylinder wall is discarded from both sides. 0.0 (the default) is a
    /// pure cylindrical surface.
    pub kerf_width: f64,
    /// The discarded side's fate, as for [`OsteotomyCut`].
    pub discarded: DiscardedSide,
}

impl CylindricalCut {
    /// Perpendicular distance from the cylinder's axis (mm).
    pub fn radial_distance(&self, p: Vec3) -> f64 {
        let dir = self.axis_direction.normalize();
        let w = p - self.axis_origin;
        let along = w.dot(dir);
        (w - dir * along).norm()
    }

    /// Applies the cut; returns the kept fragment (see
    /// [`OsteotomyCut::apply`]).
    pub fn apply(&self, model: &VoxelModel) -> VoxelModel {
        self.apply_measured(model).0
    }

    /// [`Self::apply`] plus the step's measurements, with the same
    /// conventions as [`OsteotomyCut::apply_measured`]: the resection
    /// volume counts discarded non-empty voxels, and the cut depth is the
    /// deepest discarded voxel centre past the kept face — the cylinder
    /// wall shifted by half the kerf.
    pub fn apply_measured(
        &self,
        model: &VoxelModel,
    ) -> (VoxelModel, Option<VoxelModel>, CutMeasurement) {
        let mut out = model.clone();
        let (nx, ny, nz) = model.dims;
        let voxel_volume = model.spacing.0 * model.spacing.1 * model.spacing.2;
        let retaining = matches!(self.discarded, DiscardedSide::RetainAs { .. });
        let mut retained = if retaining { Some(model.clone()) } else { None };
        let mut min = [usize::MAX; 3];
        let mut max = [0usize; 3];
        let mut resected = 0usize;
        let mut max_depth = 0.0f64;
        let half_kerf = self.kerf_width * 0.5;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let Some(idx) = model.index(x, y, z) else {
                        continue;
                    };
                    if model.values[idx].is_nan() {
                        continue;
                    }
                    let d = self.radial_distance(model.center(x, y, z));
                    let keep = if self.keep_inside {
                        d <= self.radius - half_kerf
                    } else {
                        d >= self.radius + half_kerf
                    };
                    if keep {
                        min[0] = min[0].min(x);
                        min[1] = min[1].min(y);
                        min[2] = min[2].min(z);
                        max[0] = max[0].max(x);
                        max[1] = max[1].max(y);
                        max[2] = max[2].max(z);
                        if let Some(r) = &mut retained {
                            r.values[idx] = f64::NAN;
                        }
                    } else {
                        out.values[idx] = f64::NAN;
                        if !retaining {
                            resected += 1;
                        }
                        let depth = if self.keep_inside {
                            d - (self.radius - half_kerf)
                        } else {
                            (self.radius + half_kerf) - d
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
            return (out, None, measurement); // empty cut result; leave as-is
        }
        compact(&mut out, min, max);
        let retained_model = retained.and_then(compact_retained);
        (out, retained_model, measurement)
    }
}

/// Compacts a retained fragment to its own bounding box; `None` when the
/// side retained nothing non-empty. Shared by every cut type that can
/// retain ([`OsteotomyCut`], [`CylindricalCut`]) so the two cannot drift.
fn compact_retained(retained: VoxelModel) -> Option<VoxelModel> {
    let mut rmin = [usize::MAX; 3];
    let mut rmax = [0usize; 3];
    for z in 0..retained.dims.2 {
        for y in 0..retained.dims.1 {
            for x in 0..retained.dims.0 {
                if !retained.values[retained.index(x, y, z).expect("in-bounds")].is_nan() {
                    rmin[0] = rmin[0].min(x);
                    rmin[1] = rmin[1].min(y);
                    rmin[2] = rmin[2].min(z);
                    rmax[0] = rmax[0].max(x);
                    rmax[1] = rmax[1].max(y);
                    rmax[2] = rmax[2].max(z);
                }
            }
        }
    }
    let mut r = retained;
    if rmax[0] != 0 {
        compact(&mut r, rmin, rmax);
    }
    (r.values.iter().any(|v| !v.is_nan())).then_some(r)
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
    /// Identity transform — no rotation, no translation.
    pub fn no_op() -> Self {
        Self {
            rotation_axis: Vec3::Z,
            rotation_angle: 0.0,
            pivot: Vec3::ZERO,
            translation: Vec3::ZERO,
        }
    }

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
    /// Curved (cylindrical) resection.
    Cylinder(CylindricalCut),
    /// Rigid reposition of **every** fragment (the original semantics).
    Move(FragmentTransform),
    /// Rigid reposition of one **named** fragment
    /// (`rfcs/0011-per-fragment-addressing.md`, first slice). The name
    /// must be a fragment that exists at this point in the plan.
    MoveNamed {
        /// The target fragment's name.
        fragment: String,
        /// The transform to apply.
        transform: FragmentTransform,
    },
}

/// Plan-structure rejection, raised at **build time** (when a step is
/// appended) so an invalid plan cannot be recorded, and `execute` stays
/// infallible.
#[derive(Debug, Clone, PartialEq)]
pub enum PlanError {
    /// A [`VirtualSurgery::move_fragment_named`] named a fragment that
    /// does not exist at that point in the plan.
    UnknownFragment {
        /// The unknown name.
        name: String,
    },
    /// A cut or wedge was appended after a named-fragment move. Once a
    /// fragment has moved onto its own grid, a later cut cannot be applied
    /// without a lossy re-composition (`rfcs/0011` v0 rejects rather than
    /// risks a silently wrong osteotomy). Whole-model moves do not block
    /// cuts when no fragment has been split.
    CutAfterMove {
        /// The fragment that moved before the cut.
        fragment: String,
    },
    /// A retaining cut ([`DiscardedSide::RetainAs`]) was appended when the
    /// model already holds more than one fragment: one cut cannot name a
    /// discarded side per intersected fragment.
    ModelAlreadySplit {
        /// The fragment count at the point of the cut.
        fragments: usize,
    },
    /// A [`DiscardedSide::RetainAs`] carried an empty name.
    EmptyRetainedName,
}

impl core::fmt::Display for PlanError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PlanError::UnknownFragment { name } => {
                write!(f, "no fragment named {name:?} exists at this point in the plan")
            }
            PlanError::CutAfterMove { fragment } => write!(
                f,
                "cut after a named move of {fragment:?}: re-composition across moved grids is not supported in v0"
            ),
            PlanError::ModelAlreadySplit { fragments } => write!(
                f,
                "retaining cut on a model of {fragments} fragments: one cut cannot name a discarded side per fragment"
            ),
            PlanError::EmptyRetainedName => {
                write!(f, "retained side must carry a non-empty fragment name")
            }
        }
    }
}

impl std::error::Error for PlanError {}

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

/// A named, independently addressable piece of the operated model
/// (`rfcs/0011` first slice).
#[derive(Debug, Clone)]
struct NamedFragment {
    name: Option<String>,
    model: VoxelModel,
}

/// A virtual surgery plan over a base model.
///
/// Steps are validated **as they are appended** (fragment existence, cut
/// sequencing), so `execute` is infallible: every rejection this crate can
/// detect happens at plan-building time with a [`PlanError`].
#[derive(Debug, Clone)]
pub struct VirtualSurgery {
    base: VoxelModel,
    steps: Vec<PlanStep>,
    /// Labels of fragments produced/modified, in order.
    fragment_log: Vec<String>,
    /// Build-time fragment-name state ("" = the unnamed base).
    fragment_names: Vec<String>,
    /// Build-time fragment count (mirrors `fragment_names.len()`).
    fragment_count: usize,
    /// The last named-moved fragment, if any.
    last_named_move: Option<String>,
}

impl VirtualSurgery {
    /// Starts a plan on the base (pre-operative) model.
    pub fn new(base: VoxelModel) -> Self {
        Self {
            base,
            steps: Vec::new(),
            fragment_log: Vec::new(),
            fragment_names: vec![String::new()],
            fragment_count: 1,
            last_named_move: None,
        }
    }

    /// Appends an osteotomy step.
    ///
    /// # Errors
    ///
    /// [`PlanError::CutAfterMove`] when a named-fragment move precedes the
    /// cut, [`PlanError::ModelAlreadySplit`] when the model already holds
    /// several fragments, [`PlanError::EmptyRetainedName`] for an unnamed
    /// retained side.
    pub fn cut(&mut self, cut: OsteotomyCut) -> Result<&mut Self, PlanError> {
        self.validate_cut()?;
        if let DiscardedSide::RetainAs { name } = &cut.discarded {
            if name.trim().is_empty() {
                return Err(PlanError::EmptyRetainedName);
            }
        }
        let retained = matches!(cut.discarded, DiscardedSide::RetainAs { .. });
        self.fragment_log.push(format!("cut:{}", cut.fragment_name));
        self.steps.push(PlanStep::Cut(cut));
        self.fragment_names = vec![self.latest_kept_name().to_string()];
        if retained {
            if let Some(PlanStep::Cut(c)) = self.steps.last() {
                if let DiscardedSide::RetainAs { name } = &c.discarded {
                    self.fragment_names.push(name.clone());
                }
            }
        }
        self.fragment_count = self.fragment_names.len();
        Ok(self)
    }

    /// Appends a two-plane closed-wedge step.
    ///
    /// # Errors
    ///
    /// As [`Self::cut`] (wedges always resect their discarded region).
    pub fn wedge(&mut self, cut: WedgeCut) -> Result<&mut Self, PlanError> {
        self.validate_cut()?;
        self.fragment_log
            .push(format!("wedge:{}", cut.fragment_name));
        self.steps.push(PlanStep::Wedge(cut));
        self.fragment_names = vec![self.latest_kept_name().to_string()];
        self.fragment_count = 1;
        Ok(self)
    }

    /// Appends a curved (cylindrical) resection step.
    ///
    /// # Errors
    ///
    /// As [`Self::cut`] (including retention via
    /// [`DiscardedSide::RetainAs`]).
    pub fn cylinder(&mut self, cut: CylindricalCut) -> Result<&mut Self, PlanError> {
        self.validate_cut()?;
        if let DiscardedSide::RetainAs { name } = &cut.discarded {
            if name.trim().is_empty() {
                return Err(PlanError::EmptyRetainedName);
            }
        }
        let retained = matches!(cut.discarded, DiscardedSide::RetainAs { .. });
        self.fragment_log
            .push(format!("cylinder:{}", cut.fragment_name));
        self.steps.push(PlanStep::Cylinder(cut));
        self.fragment_names = vec![self.latest_kept_name().to_string()];
        if retained {
            if let Some(PlanStep::Cylinder(c)) = self.steps.last() {
                if let DiscardedSide::RetainAs { name } = &c.discarded {
                    self.fragment_names.push(name.clone());
                }
            }
        }
        self.fragment_count = self.fragment_names.len();
        Ok(self)
    }

    /// Appends a whole-model fragment transform — every fragment moves
    /// together, preserving their relative positions (the original
    /// semantics; existing plans are unaffected).
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

    /// Appends a reposition of **one named fragment** — the per-fragment
    /// addressing of `rfcs/0011`. The fragment must exist at this point in
    /// the plan (a cut's kept side, or a retained discarded side).
    ///
    /// # Errors
    ///
    /// [`PlanError::UnknownFragment`] when no fragment carries `name`.
    pub fn move_fragment_named(
        &mut self,
        name: impl Into<String>,
        transform: FragmentTransform,
    ) -> Result<&mut Self, PlanError> {
        let name = name.into();
        if !self.fragment_names.iter().any(|n| n == &name) {
            return Err(PlanError::UnknownFragment { name });
        }
        self.fragment_log.push(format!(
            "move_named:{name}:rotate({:.3} rad) + translate({:.1}, {:.1}, {:.1}) mm",
            transform.rotation_angle,
            transform.translation.x,
            transform.translation.y,
            transform.translation.z
        ));
        self.steps.push(PlanStep::MoveNamed {
            fragment: name.clone(),
            transform,
        });
        self.last_named_move = Some(name);
        Ok(self)
    }

    /// Shared cut/wedge sequencing checks.
    fn validate_cut(&self) -> Result<(), PlanError> {
        if let Some(fragment) = &self.last_named_move {
            return Err(PlanError::CutAfterMove {
                fragment: fragment.clone(),
            });
        }
        if self.fragment_count > 1 {
            return Err(PlanError::ModelAlreadySplit {
                fragments: self.fragment_count,
            });
        }
        Ok(())
    }

    /// The kept-side name of the most recent cut/wedge step.
    fn latest_kept_name(&self) -> &str {
        for step in self.steps.iter().rev() {
            match step {
                PlanStep::Cut(c) => return &c.fragment_name,
                PlanStep::Wedge(w) => return &w.fragment_name,
                PlanStep::Cylinder(c) => return &c.fragment_name,
                _ => continue,
            }
        }
        ""
    }

    /// Executes the plan, returning the operated model and the audit log
    /// (step descriptions in execution order). Infallible: the plan was
    /// validated as it was built.
    pub fn execute(&self) -> (VoxelModel, Vec<String>) {
        let (model, report) = self.execute_with_report();
        (model, report.step_log)
    }

    /// [`Self::execute`] with the measurement report: resection volumes and
    /// cut depths per cut, achieved-vs-prescribed centroid displacement per
    /// move (the alignment error), and the total resection volume. The
    /// report is index-aligned with the audit log, so a submission bundle
    /// can attach the numbers to the steps they belong to.
    ///
    /// Execution tracks one [`NamedFragment`] per piece. The operated model
    /// is the composition of all fragments — for plans that never split
    /// the model (every existing plan) this is exactly the single-model
    /// executor this crate has always had, byte for byte.
    pub fn execute_with_report(&self) -> (VoxelModel, SurgeryReport) {
        let mut fragments: Vec<NamedFragment> = vec![NamedFragment {
            name: None,
            model: self.base.clone(),
        }];
        let mut log = Vec::with_capacity(self.steps.len());
        let mut measurements = Vec::with_capacity(self.steps.len());
        let mut total_resection = 0.0f64;
        for step in &self.steps {
            match step {
                PlanStep::Cut(cut) => {
                    // Validated at build: exactly one fragment.
                    let current = &fragments[0].model;
                    let (kept, retained, m) = cut.apply_measured(current);
                    total_resection += m.resection_volume_mm3;
                    log.push(format!("cut:{}", cut.fragment_name));
                    measurements.push(StepMeasurement::Cut(m));
                    fragments[0] = NamedFragment {
                        name: Some(cut.fragment_name.clone()),
                        model: kept,
                    };
                    if let (Some(r), DiscardedSide::RetainAs { name }) = (retained, &cut.discarded)
                    {
                        fragments.push(NamedFragment {
                            name: Some(name.clone()),
                            model: r,
                        });
                    }
                }
                PlanStep::Wedge(cut) => {
                    let current = &fragments[0].model;
                    let (next, m) = cut.apply_measured(current);
                    total_resection += m.resection_volume_mm3;
                    log.push(format!("wedge:{}", cut.fragment_name));
                    measurements.push(StepMeasurement::Cut(m));
                    fragments[0] = NamedFragment {
                        name: Some(cut.fragment_name.clone()),
                        model: next,
                    };
                }
                PlanStep::Cylinder(cut) => {
                    // Validated at build: exactly one fragment.
                    let current = &fragments[0].model;
                    let (kept, retained, m) = cut.apply_measured(current);
                    total_resection += m.resection_volume_mm3;
                    log.push(format!("cylinder:{}", cut.fragment_name));
                    measurements.push(StepMeasurement::Cut(m));
                    fragments[0] = NamedFragment {
                        name: Some(cut.fragment_name.clone()),
                        model: kept,
                    };
                    if let (Some(r), DiscardedSide::RetainAs { name }) = (retained, &cut.discarded)
                    {
                        fragments.push(NamedFragment {
                            name: Some(name.clone()),
                            model: r,
                        });
                    }
                }
                PlanStep::Move(m) => {
                    let before = union_centroid(&fragments);
                    for f in &mut fragments {
                        f.model = m.apply_to_model(&f.model);
                    }
                    let after = union_centroid(&fragments);
                    let achieved = after - before;
                    let prescribed = m.apply_to_point(before) - before;
                    log.push(format!("move:{}", format_args!("{:?}", m.translation)));
                    measurements.push(StepMeasurement::Move(MoveMeasurement {
                        prescribed_centroid_translation: prescribed,
                        achieved_centroid_translation: achieved,
                        alignment_error_mm: (prescribed - achieved).norm(),
                    }));
                }
                PlanStep::MoveNamed {
                    fragment,
                    transform,
                } => {
                    let Some(f) = fragments
                        .iter_mut()
                        .find(|f| f.name.as_deref() == Some(fragment.as_str()))
                    else {
                        unreachable!("validated at build time");
                    };
                    let before = fragment_centroid(&f.model);
                    f.model = transform.apply_to_model(&f.model);
                    let after = fragment_centroid(&f.model);
                    let achieved = after - before;
                    let prescribed = transform.apply_to_point(before) - before;
                    log.push(format!(
                        "move_named:{fragment}:{}",
                        format_args!("{:?}", transform.translation)
                    ));
                    measurements.push(StepMeasurement::Move(MoveMeasurement {
                        prescribed_centroid_translation: prescribed,
                        achieved_centroid_translation: achieved,
                        alignment_error_mm: (prescribed - achieved).norm(),
                    }));
                }
            }
        }
        let model = compose_fragments(&fragments, &self.base);
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

/// Scatters every fragment into one grid covering the union extent
/// (nearest-neighbour by voxel centre; all fragments share the base
/// spacing). A single-fragment composition is the identity, which is what
/// keeps every pre-existing plan byte-identical.
fn compose_fragments(fragments: &[NamedFragment], base: &VoxelModel) -> VoxelModel {
    if fragments.len() == 1 {
        return fragments[0].model.clone();
    }
    let spacing = base.spacing;
    let mut min_corner = fragments[0].model.origin;
    let mut max_corner = fragments[0].model.origin;
    for f in fragments {
        let m = &f.model;
        let far = m.origin
            + Vec3::new(
                m.dims.0 as f64 * spacing.0,
                m.dims.1 as f64 * spacing.1,
                m.dims.2 as f64 * spacing.2,
            );
        min_corner = Vec3::new(
            min_corner.x.min(m.origin.x),
            min_corner.y.min(m.origin.y),
            min_corner.z.min(m.origin.z),
        );
        max_corner = Vec3::new(
            max_corner.x.max(far.x),
            max_corner.y.max(far.y),
            max_corner.z.max(far.z),
        );
    }
    let dims = (
        ((max_corner.x - min_corner.x) / spacing.0).round() as usize + 1,
        ((max_corner.y - min_corner.y) / spacing.1).round() as usize + 1,
        ((max_corner.z - min_corner.z) / spacing.2).round() as usize + 1,
    );
    let mut values = vec![f64::NAN; dims.0 * dims.1 * dims.2];
    for f in fragments {
        let m = &f.model;
        for z in 0..m.dims.2 {
            for y in 0..m.dims.1 {
                for x in 0..m.dims.0 {
                    let Some(idx) = m.index(x, y, z) else {
                        continue;
                    };
                    let v = m.values[idx];
                    if v.is_nan() {
                        continue;
                    }
                    let c = m.center(x, y, z);
                    let gx = ((c.x - min_corner.x) / spacing.0).round();
                    let gy = ((c.y - min_corner.y) / spacing.1).round();
                    let gz = ((c.z - min_corner.z) / spacing.2).round();
                    if gx < 0.0 || gy < 0.0 || gz < 0.0 {
                        continue;
                    }
                    let (gx, gy, gz) = (gx as usize, gy as usize, gz as usize);
                    if gx < dims.0 && gy < dims.1 && gz < dims.2 {
                        values[(gz * dims.1 + gy) * dims.0 + gx] = v;
                    }
                }
            }
        }
    }
    VoxelModel {
        dims,
        spacing,
        origin: min_corner,
        values,
    }
}

/// Centroid over every fragment of a tracked state.
fn union_centroid(fragments: &[NamedFragment]) -> Vec3 {
    let mut sum = Vec3::ZERO;
    let mut count = 0usize;
    for f in fragments {
        let (s, c) = centroid_sum(&f.model);
        sum += s;
        count += c;
    }
    if count == 0 {
        Vec3::ZERO
    } else {
        sum / count as f64
    }
}

/// Centroid of the non-empty (non-NaN) voxels of a model, in patient
/// coordinates. `Vec3::ZERO` for an empty model.
fn fragment_centroid(model: &VoxelModel) -> Vec3 {
    let (sum, count) = centroid_sum(model);
    if count == 0 {
        Vec3::ZERO
    } else {
        sum / count as f64
    }
}

/// `(coordinate sum, voxel count)` of a model's non-empty voxels.
fn centroid_sum(model: &VoxelModel) -> (Vec3, usize) {
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
    (sum, count)
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
            discarded: DiscardedSide::Resect,
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
            discarded: DiscardedSide::Resect,
        })
        .expect("valid plan");
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
            discarded: DiscardedSide::Resect,
        })
        .expect("valid plan");
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
            discarded: DiscardedSide::Resect,
        };
        let (out, _, m) = cut.apply_measured(&model);
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
            discarded: DiscardedSide::Resect,
        };
        let (kept0, _, m0) = base.apply_measured(&model);
        assert_eq!(kept0.count_above(50.0), 3 * 36);
        assert!((m0.max_depth_mm - 5.0).abs() < 1e-9);
        // A 2 mm kerf removes the slab (−1, +1) around the plane: the
        // kept boundary moves from centre 0 to centre +1, costing the
        // slice at centre 0 (36 voxels of cube). Depth is measured past
        // the kept face at +1 mm, so the deepest voxel (centre −5) sits
        // 6 mm from it.
        base.kerf_width = 2.0;
        let (kept2, _, m2) = base.apply_measured(&model);
        assert_eq!(kept2.count_above(50.0), 2 * 36);
        assert!((m2.max_depth_mm - 6.0).abs() < 1e-9, "{}", m2.max_depth_mm);
        // The slab is resected tissue: volume grows by the slice.
        assert!((m2.resection_volume_mm3 - m0.resection_volume_mm3 - 100.0).abs() < 1e-9);
        // The negative side with kerf: the slab is removed from there too.
        base.keep_positive = false;
        let (kept_neg, _, _) = base.apply_measured(&model);
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
        plan.wedge(wedge).expect("valid plan");
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

    /// The cube grid's column count with centres x² + y² within a
    /// radius² bound (x, y centres −5..4, a 10-wide grid): {r² ≤ 9} →
    /// 29 columns, {r² ≤ 4} → 13, the complement → 71.
    fn columns_within(radius_squared: f64) -> usize {
        let mut n = 0;
        for x in -5i32..5 {
            for y in -5i32..5 {
                if (x * x + y * y) as f64 <= radius_squared {
                    n += 1;
                }
            }
        }
        n
    }

    fn core_cut(
        radius: f64,
        kerf: f64,
        keep_inside: bool,
        discarded: DiscardedSide,
    ) -> CylindricalCut {
        CylindricalCut {
            axis_origin: Vec3::ZERO,
            axis_direction: Vec3::Z,
            radius,
            fragment_name: "core".into(),
            keep_inside,
            kerf_width: kerf,
            discarded,
        }
    }

    #[test]
    fn cylindrical_cut_keeps_exactly_the_core() {
        let model = cube_model();
        let cut = core_cut(3.0, 0.0, true, DiscardedSide::Resect);
        let (out, _, m) = cut.apply_measured(&model);
        // Columns within r = 3 exist in every z slice (the axis runs the
        // full grid), zero-valued tissue included — 29 × 10.
        let cols = columns_within(9.0);
        assert_eq!(cols, 29);
        // Cube voxels in the core: the cube's x/y centres run −3..2, so
        // the (3,0) column and the four (x,3) columns are outside it —
        // 27 of the 29 core columns are cube.
        assert_eq!(out.count_above(50.0), 27 * 6, "cube voxels in the core");
        let kept_total = out.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(kept_total, cols * 10);
        // Resection volume: everything not kept (tissue + cube corners).
        assert!(
            (m.resection_volume_mm3 - (1000.0 - (cols * 10) as f64)).abs() < 1e-9,
            "{}",
            m.resection_volume_mm3
        );
        // Depth past the wall: the farthest discarded centre is a grid
        // corner at r = √50.
        assert!((m.max_depth_mm - (50.0f64.sqrt() - 3.0)).abs() < 1e-9);
        // The kept bounding box compacts to the core's extent.
        assert_eq!(out.dims, (7, 7, 10));
        assert!((out.origin.x - (-3.0)).abs() < 1e-12);
    }

    #[test]
    fn cylindrical_cut_keep_outside_resects_the_core() {
        let model = cube_model();
        let cut = core_cut(3.0, 0.0, false, DiscardedSide::Resect);
        let (out, _, m) = cut.apply_measured(&model);
        // Kept: the annulus d ≥ 3. The four r = 3 boundary columns are
        // kept by *both* settings (the core keeps d ≤ 3) — the same
        // measure-zero boundary overlap the plane cut has with its
        // d ≥ 0 / d ≤ 0 convention — so the annulus is 71 + 4 = 75 of
        // the grid's 100 columns × 10 slices.
        let annulus_cols = 75;
        let kept_total = out.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(kept_total, annulus_cols * 10);
        // Cube voxels in the annulus: 36 cube columns minus the 25 the
        // resection took (the 27 core∩cube columns less the two r = 3
        // boundary columns it keeps).
        assert_eq!(out.count_above(50.0), (36 - 25) * 6);
        // Resected: the strict core interior.
        assert!(
            (m.resection_volume_mm3 - (25.0 * 10.0)).abs() < 1e-9,
            "{}",
            m.resection_volume_mm3
        );
        // Depth into the core: the axis column (r = 0) is the deepest.
        assert!((m.max_depth_mm - 3.0).abs() < 1e-9);
    }

    #[test]
    fn cylindrical_kerf_removes_a_radial_slab() {
        let model = cube_model();
        // Zero-kerf baseline.
        let (_, _, m0) = core_cut(3.0, 0.0, true, DiscardedSide::Resect).apply_measured(&model);
        // A 2 mm kerf keeps r ≤ 2 only: 13 columns × 10 slices, so the
        // resection grows by exactly the (29 − 13)-column shell, and the
        // depth is measured past the kept face at r = 2.
        let cut = core_cut(3.0, 2.0, true, DiscardedSide::Resect);
        let (out, _, m2) = cut.apply_measured(&model);
        let cols = columns_within(4.0);
        assert_eq!(cols, 13);
        let kept_total = out.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(kept_total, cols * 10);
        assert!(
            (m2.resection_volume_mm3 - (m0.resection_volume_mm3 + 160.0)).abs() < 1e-9,
            "{} vs {}",
            m2.resection_volume_mm3,
            m0.resection_volume_mm3
        );
        assert!(m2.resection_volume_mm3 > m0.resection_volume_mm3);
        assert!((m2.max_depth_mm - (50.0f64.sqrt() - 2.0)).abs() < 1e-9);
    }

    #[test]
    fn cylindrical_axis_direction_and_offset_are_honoured() {
        let model = cube_model();
        // A non-unit axis direction and an offset axis: the cut is the
        // radius-1.5 cylinder about the line (x=3, y=0). Columns within
        // 1.5 of (3, 0): dy=0 → dx ∈ {−1,0,1} → 3; dy=±1 → dx² + 1 ≤ 2.25
        // → dx ∈ {−1,0,1} → 3 each. 9 columns, all inside the grid.
        let cut = CylindricalCut {
            axis_origin: Vec3::new(3.0, 0.0, 0.0),
            axis_direction: Vec3::new(0.0, 0.0, 7.0),
            radius: 1.5,
            fragment_name: "core".into(),
            keep_inside: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        };
        let (out, _, m) = cut.apply_measured(&model);
        let kept_total = out.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(kept_total, 9 * 10, "9 columns of the offset cylinder");
        // Cube columns in the core: only centre x = 2 is inside the
        // cube's x centres −3..2 (the axis sits at x = 3), y ∈ {−1,0,1}
        // → 3 per slice.
        assert_eq!(out.count_above(50.0), 3 * 6);
        assert_eq!(m.fragment_name, "core");
    }

    #[test]
    fn cylindrical_cut_integrates_with_retention_and_named_moves() {
        let model = cube_model();
        let mut plan = VirtualSurgery::new(model);
        plan.cylinder(core_cut(
            3.0,
            0.0,
            true,
            DiscardedSide::RetainAs {
                name: "graft".into(),
            },
        ))
        .expect("single-fragment model accepts a retaining cylinder");
        // Distract the kept core +2 mm axially; the retained annulus
        // stays. A z-translation keeps every core voxel in its own
        // column, and the core/annulus columns are disjoint, so the
        // composition conserves every voxel exactly.
        plan.move_fragment_named(
            "core",
            FragmentTransform {
                rotation_axis: Vec3::Z,
                rotation_angle: 0.0,
                pivot: Vec3::ZERO,
                translation: Vec3::new(0.0, 0.0, 2.0),
            },
        )
        .expect("core exists");
        let (operated, report) = plan.execute_with_report();
        assert!(report.step_log[0].starts_with("cylinder:core"));
        assert!(report.step_log[1].starts_with("move_named:core:"));
        // Retention is not resection.
        assert!((report.total_resection_volume_mm3 - 0.0).abs() < 1e-9);
        // Both sides survive the composition: core 29×10 + annulus
        // 71×10 = the full 1000 voxels.
        assert_eq!(operated.values.iter().filter(|v| !v.is_nan()).count(), 1000);
        assert_eq!(operated.count_above(50.0), 36 * 6);
        // The composed grid grew axially to hold the distracted core.
        let far_z = operated.origin.z + operated.dims.2 as f64 * operated.spacing.2;
        assert!(far_z > 5.0, "composed z extent {far_z}");
    }

    #[test]
    fn cylindrical_cut_respects_plan_validation() {
        let model = cube_model();
        // A cylinder after a plain cut is a legal second step.
        let mut plan = VirtualSurgery::new(model.clone());
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
            fragment_name: "proximal".into(),
            keep_positive: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        })
        .expect("valid");
        assert!(plan
            .cylinder(core_cut(2.0, 0.0, true, DiscardedSide::Resect))
            .is_ok());
        // A cylinder cannot follow a named move.
        let mut plan2 = VirtualSurgery::new(model.clone());
        plan2
            .cut(OsteotomyCut {
                plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
                fragment_name: "a".into(),
                keep_positive: true,
                kerf_width: 0.0,
                discarded: DiscardedSide::RetainAs { name: "b".into() },
            })
            .expect("valid");
        plan2
            .move_fragment_named("b", FragmentTransform::no_op())
            .expect("b exists");
        assert_eq!(
            plan2
                .cylinder(core_cut(2.0, 0.0, true, DiscardedSide::Resect))
                .err(),
            Some(PlanError::CutAfterMove {
                fragment: "b".into()
            })
        );
        // Nor a retaining cylinder on an already-split model.
        let mut plan3 = VirtualSurgery::new(model);
        plan3
            .cut(OsteotomyCut {
                plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
                fragment_name: "a".into(),
                keep_positive: true,
                kerf_width: 0.0,
                discarded: DiscardedSide::RetainAs { name: "b".into() },
            })
            .expect("valid");
        assert_eq!(
            plan3
                .cylinder(core_cut(
                    2.0,
                    0.0,
                    true,
                    DiscardedSide::RetainAs { name: "d".into() }
                ))
                .err(),
            Some(PlanError::ModelAlreadySplit { fragments: 2 })
        );
    }

    #[test]
    fn retaining_cut_keeps_both_sides_as_addressable_fragments() {
        let model = cube_model();
        let mut plan = VirtualSurgery::new(model);
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
            fragment_name: "distal".into(),
            keep_positive: false,
            kerf_width: 0.0,
            discarded: DiscardedSide::RetainAs {
                name: "proximal".into(),
            },
        })
        .expect("single-fragment model accepts a retaining cut");
        // Move only the distal piece 3 mm laterally — a collision-free
        // distraction (moving it through the retained piece would overlap
        // voxels, which the composition records last-write-wins).
        plan.move_fragment_named(
            "distal",
            FragmentTransform {
                rotation_axis: Vec3::Z,
                rotation_angle: 0.0,
                pivot: Vec3::ZERO,
                translation: Vec3::new(3.0, 0.0, 0.0),
            },
        )
        .expect("distal exists");
        let (operated, log) = plan.execute();
        assert!(log[0].starts_with("cut:distal"));
        assert!(log[1].starts_with("move_named:distal:"));
        // Both pieces survive: the full 6³ cube (216 voxels ≥ 50 HU) is
        // present in the composed model.
        assert_eq!(operated.count_above(50.0), 216);
        // And the composed grid grew laterally to hold the moved piece:
        // its x extent now exceeds the original cube's.
        let max_x = operated.origin.x + operated.dims.0 as f64 * operated.spacing.0;
        assert!(max_x > 5.0, "composed x extent {max_x}");
    }

    #[test]
    fn retained_side_is_not_counted_as_resection() {
        let model = cube_model();
        let mut plan = VirtualSurgery::new(model);
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
            fragment_name: "kept".into(),
            keep_positive: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::RetainAs {
                name: "held".into(),
            },
        })
        .expect("valid");
        let (_, report) = plan.execute_with_report();
        // Nothing is resected: the discarded side was retained.
        assert!((report.total_resection_volume_mm3 - 0.0).abs() < 1e-9);
    }

    #[test]
    fn plan_errors_are_raised_at_build_time() {
        let model = cube_model();
        // Unknown fragment name.
        let mut plan = VirtualSurgery::new(model.clone());
        let err = plan
            .move_fragment_named("ghost", FragmentTransform::no_op())
            .unwrap_err();
        assert_eq!(
            err,
            PlanError::UnknownFragment {
                name: "ghost".into()
            }
        );
        // Cut after a named move.
        let mut plan = VirtualSurgery::new(model.clone());
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
            fragment_name: "a".into(),
            keep_positive: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::RetainAs { name: "b".into() },
        })
        .expect("valid");
        plan.move_fragment_named("b", FragmentTransform::no_op())
            .expect("b exists");
        let err = plan
            .cut(OsteotomyCut {
                plane: Plane::from_point_normal(Vec3::ZERO, Vec3::X).unwrap(),
                fragment_name: "c".into(),
                keep_positive: true,
                kerf_width: 0.0,
                discarded: DiscardedSide::Resect,
            })
            .unwrap_err();
        assert_eq!(
            err,
            PlanError::CutAfterMove {
                fragment: "b".into()
            }
        );
        // Retaining cut on an already-split model.
        let mut plan2 = VirtualSurgery::new(model);
        plan2
            .cut(OsteotomyCut {
                plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
                fragment_name: "a".into(),
                keep_positive: true,
                kerf_width: 0.0,
                discarded: DiscardedSide::RetainAs { name: "b".into() },
            })
            .expect("valid");
        let err = plan2
            .cut(OsteotomyCut {
                plane: Plane::from_point_normal(Vec3::ZERO, Vec3::X).unwrap(),
                fragment_name: "c".into(),
                keep_positive: true,
                kerf_width: 0.0,
                discarded: DiscardedSide::RetainAs { name: "d".into() },
            })
            .unwrap_err();
        assert_eq!(err, PlanError::ModelAlreadySplit { fragments: 2 });
        // Empty retained name.
        let mut plan3 = VirtualSurgery::new(cube_model());
        assert_eq!(
            plan3
                .cut(OsteotomyCut {
                    plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap(),
                    fragment_name: "a".into(),
                    keep_positive: true,
                    kerf_width: 0.0,
                    discarded: DiscardedSide::RetainAs { name: "  ".into() },
                })
                .err(),
            Some(PlanError::EmptyRetainedName)
        );
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
            discarded: DiscardedSide::Resect,
        })
        .expect("valid plan");
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::X).unwrap(),
            fragment_name: "b".into(),
            keep_positive: false,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        })
        .expect("valid plan");
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
