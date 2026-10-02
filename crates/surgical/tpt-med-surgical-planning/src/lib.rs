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

/// A **freeform (anatomically contoured) resection**: the cut surface is a
/// caller-supplied **closed triangle mesh** — the patient-matched
/// implant/back-of-the-condyle surface a planer or a patient-specific
/// jig actually follows, which no plane, wedge or cylinder can express.
/// The kept region is the mesh's interior (`keep_inside`) or its
/// exterior, with the same optional saw-kerf and `rfcs/0011`
/// discarded-side conventions as the plane and cylindrical cuts.
///
/// Side determination is the standard signed-distance construction: the
/// unsigned distance to the closest triangle, signed by a +x ray-parity
/// inside test — which is only meaningful for a **watertight** mesh, so
/// [`Self::validate`] checks closure (every undirected edge shared by
/// exactly two triangles) and consistent winding (each directed edge
/// exactly once) and is called at plan-append time. Rays are cast along
/// the grid's x axis; a mesh whose edges deliberately thread every voxel
/// centre row is outside this model's scope.
#[derive(Debug, Clone)]
pub struct MeshCut {
    /// The contour surface as a triangle soup (three vertices per
    /// triangle, patient coordinates). Must be closed with consistent
    /// winding — [`Self::validate`] enforces both.
    pub triangles: Vec<[Vec3; 3]>,
    /// Fragment name produced from the kept side (audit label).
    pub fragment_name: String,
    /// Keep the mesh's interior (the contoured region itself — e.g. a
    /// graft block shaped to a defect); `false` keeps the exterior and
    /// resects the interior (a contoured resection cap).
    pub keep_inside: bool,
    /// Saw-kerf width (mm): material within `kerf_width / 2` of the
    /// surface on either side is discarded. 0.0 is a pure contoured cut.
    pub kerf_width: f64,
    /// The discarded side's fate, as for [`OsteotomyCut`].
    pub discarded: DiscardedSide,
}

/// Mesh-contour rejection, with the offending edge when one is at fault.
#[derive(Debug, Clone, PartialEq)]
pub enum MeshCutError {
    /// No triangles.
    Empty,
    /// A triangle with zero area (repeated or collinear vertices).
    DegenerateTriangle {
        /// The triangle's index.
        index: usize,
    },
    /// The mesh is not watertight: some undirected edge is shared by a
    /// number of triangles other than two.
    OpenSurface {
        /// The offending edge's two endpoints.
        edge: (Vec3, Vec3),
    },
    /// The mesh is closed but its winding is inconsistent: some directed
    /// edge appears more than once (a flipped triangle).
    InconsistentWinding {
        /// The offending edge's two endpoints.
        edge: (Vec3, Vec3),
    },
}

impl core::fmt::Display for MeshCutError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let fmt_edge = |e: &(Vec3, Vec3)| {
            format!(
                "({:.3},{:.3},{:.3})-({:.3},{:.3},{:.3})",
                e.0.x, e.0.y, e.0.z, e.1.x, e.1.y, e.1.z
            )
        };
        match self {
            MeshCutError::Empty => write!(f, "contour mesh has no triangles"),
            MeshCutError::DegenerateTriangle { index } => {
                write!(f, "contour mesh triangle {index} is degenerate")
            }
            MeshCutError::OpenSurface { edge } => {
                write!(
                    f,
                    "contour mesh is not watertight at edge {}",
                    fmt_edge(edge)
                )
            }
            MeshCutError::InconsistentWinding { edge } => write!(
                f,
                "contour mesh winding is inconsistent at edge {}",
                fmt_edge(edge)
            ),
        }
    }
}

impl std::error::Error for MeshCutError {}

/// Quantizes a vertex to the edge-identity key (1 µm grid): vertices
/// closer than this are the same corner of the surface.
fn vertex_key(p: Vec3) -> [i64; 3] {
    const SCALE: f64 = 1.0e6;
    [
        (p.x * SCALE).round() as i64,
        (p.y * SCALE).round() as i64,
        (p.z * SCALE).round() as i64,
    ]
}

/// Squared distance from a point to a triangle (vertex/edge/face regions,
/// the standard closest-point-on-triangle decomposition).
fn point_triangle_distance2(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;
    let d1 = ab.dot(ap);
    let d2 = ac.dot(ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return ap.norm_squared();
    }
    let bp = p - b;
    let d3 = ab.dot(bp);
    let d4 = ac.dot(bp);
    if d3 >= 0.0 && d4 <= d3 {
        return bp.norm_squared();
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let t = d1 / (d1 - d3);
        let q = a + ab * t;
        return (p - q).norm_squared();
    }
    let cp = p - c;
    let d5 = ab.dot(cp);
    let d6 = ac.dot(cp);
    if d6 >= 0.0 && d5 <= d6 {
        return cp.norm_squared();
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let t = d2 / (d2 - d6);
        let q = a + ac * t;
        return (p - q).norm_squared();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let q = b + (c - b) * t;
        return (p - q).norm_squared();
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let q = a + ab * v + ac * w;
    (p - q).norm_squared()
}

/// Signed distance from `p` to the closed triangle mesh: the unsigned
/// distance to the closest triangle, negative inside by +x ray parity.
fn mesh_signed_distance(p: Vec3, triangles: &[[Vec3; 3]]) -> f64 {
    let mut best2 = f64::INFINITY;
    for t in triangles {
        best2 = best2.min(point_triangle_distance2(p, t[0], t[1], t[2]));
    }
    // Ray parity along +x (Möller–Trumbore).
    let dir = Vec3::X;
    let mut crossings = 0usize;
    for t in triangles {
        let e1 = t[1] - t[0];
        let e2 = t[2] - t[0];
        let h = dir.cross(e2);
        let a = e1.dot(h);
        if a.abs() < 1e-12 {
            continue;
        }
        let f = 1.0 / a;
        let s = p - t[0];
        let u = f * s.dot(h);
        if !(-1e-9..=1.0 + 1e-9).contains(&u) {
            continue;
        }
        let q = s.cross(e1);
        let v = f * dir.dot(q);
        if v < -1e-9 || u + v > 1.0 + 1e-9 {
            continue;
        }
        let t_hit = f * e2.dot(q);
        if t_hit > 1e-9 {
            crossings += 1;
        }
    }
    let inside = crossings % 2 == 1;
    let unsigned = best2.sqrt();
    if inside {
        -unsigned
    } else {
        unsigned
    }
}

/// Structural checks for a closed triangle mesh (the contract
/// [`MeshCut`], [`ImplantPlacement`] and [`GraftReconstruction`] all
/// need): non-empty, non-degenerate triangles, watertight (every
/// undirected edge shared by exactly two triangles) and consistently
/// wound (every directed edge exactly once — a flipped triangle is
/// caught here, before a parity test could silently invert a region).
fn validate_closed_mesh(triangles: &[[Vec3; 3]]) -> Result<(), MeshCutError> {
    {
        if triangles.is_empty() {
            return Err(MeshCutError::Empty);
        }
        // (undirected key) -> (count, directed: a->b count, b->a count)
        let mut edges: std::collections::BTreeMap<[i64; 6], (usize, usize, usize)> =
            std::collections::BTreeMap::new();
        for (i, t) in triangles.iter().enumerate() {
            let n = (t[1] - t[0]).cross(t[2] - t[0]);
            if !n.norm_squared().is_finite() || n.norm_squared() <= 1e-20 {
                return Err(MeshCutError::DegenerateTriangle { index: i });
            }
            let vs = [vertex_key(t[0]), vertex_key(t[1]), vertex_key(t[2])];
            for k in 0..3 {
                let (a, b) = (vs[k], vs[(k + 1) % 3]);
                let (mut lo, mut hi) = (a, b);
                if hi < lo {
                    std::mem::swap(&mut lo, &mut hi);
                }
                let key = [lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]];
                let e = edges.entry(key).or_insert((0, 0, 0));
                e.0 += 1;
                if b > a {
                    e.1 += 1;
                } else {
                    e.2 += 1;
                }
            }
        }
        for (key, (count, forward, backward)) in edges {
            let a = Vec3::new(
                key[0] as f64 * 1e-6,
                key[1] as f64 * 1e-6,
                key[2] as f64 * 1e-6,
            );
            let b = Vec3::new(
                key[3] as f64 * 1e-6,
                key[4] as f64 * 1e-6,
                key[5] as f64 * 1e-6,
            );
            if count != 2 {
                return Err(MeshCutError::OpenSurface { edge: (a, b) });
            }
            if forward != 1 || backward != 1 {
                return Err(MeshCutError::InconsistentWinding { edge: (a, b) });
            }
        }
        Ok(())
    }
}

impl MeshCut {
    /// Checks the contour mesh's structural contract — see
    /// [`validate_closed_mesh`].
    pub fn validate(&self) -> Result<(), MeshCutError> {
        validate_closed_mesh(&self.triangles)
    }

    /// Applies the cut; returns the kept fragment. Panics on a mesh that
    /// fails [`Self::validate`] — `VirtualSurgery::mesh` validates at
    /// plan-append time, and a direct caller should too.
    pub fn apply(&self, model: &VoxelModel) -> VoxelModel {
        self.apply_measured(model).0
    }
    /// [`Self::apply`] plus the step's measurements, with the same
    /// conventions as [`OsteotomyCut::apply_measured`]: the resection
    /// volume counts discarded non-empty voxels, and the cut depth is the
    /// deepest discarded voxel centre past the kept face — the contour
    /// surface shifted by half the kerf.
    pub fn apply_measured(
        &self,
        model: &VoxelModel,
    ) -> (VoxelModel, Option<VoxelModel>, CutMeasurement) {
        if let Err(e) = self.validate() {
            panic!("mesh cut: {e}");
        }
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
                    let d = mesh_signed_distance(model.center(x, y, z), &self.triangles);
                    let keep = if self.keep_inside {
                        d <= -half_kerf
                    } else {
                        d >= half_kerf
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
                            -half_kerf - d
                        } else {
                            d - half_kerf
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

/// Measurements of an [`ImplantPlacement`] step: the bone the implant
/// replaced, the implant's own occupied volume, and the **bone–implant
/// interface area** — the shared-face area between implant voxels and
/// surviving bone, the quantity osseointegration screening is about.
#[derive(Debug, Clone)]
pub struct ImplantMeasurement {
    /// The implant's audit label.
    pub fragment_name: String,
    /// Bone volume the implant replaced (mm³).
    pub resection_volume_mm3: f64,
    /// Implant occupied volume (mm³).
    pub implant_volume_mm3: f64,
    /// Shared-face area between implant and surviving bone (mm²).
    pub interface_area_mm2: f64,
}

/// Measurements of a [`GraftReconstruction`] step: the graft's occupied
/// volume and its interface area with the surrounding bone.
#[derive(Debug, Clone)]
pub struct GraftMeasurement {
    /// The graft's audit label.
    pub fragment_name: String,
    /// Defect volume the graft filled (mm³).
    pub graft_volume_mm3: f64,
    /// Shared-face area between graft and bone (mm²).
    pub interface_area_mm2: f64,
}

/// **Implant component placement**: a closed triangle mesh defines the
/// component's shape and pose (bake the pose into the vertices — this
/// crate does not move implants independently), and on execution every
/// voxel whose centre lies inside the mesh becomes the implant (value
/// `marker_value`), the bone it replaces counting as resected. The
/// measurement reports the resection, the implant volume, and the
/// bone–implant **interface area** (shared faces with surviving bone).
///
/// Sequencing follows the cut rules ([`PlanError::CutAfterMove`],
/// [`PlanError::ModelAlreadySplit`]): placement edits the operated
/// fragment's grid. There is deliberately **no stem/cement mechanics
/// here** — this is the geometric placement and interface bookkeeping, a
/// placement plan, not a fixation-strength model.
#[derive(Debug, Clone)]
pub struct ImplantPlacement {
    /// The component's closed contour surface, pose baked in. Must pass
    /// [`MeshCut::validate`]'s contract (watertight, consistently
    /// wound); checked at plan-append time.
    pub triangles: Vec<[Vec3; 3]>,
    /// Audit label.
    pub fragment_name: String,
    /// Scalar value written into implant voxels (e.g. a metal HU marker
    /// for downstream meshing/auditing).
    pub marker_value: f64,
}

impl ImplantPlacement {
    /// Mesh structural check.
    pub fn validate(&self) -> Result<(), MeshCutError> {
        validate_closed_mesh(&self.triangles)
    }
}

/// **Bone graft / defect reconstruction**: a closed triangle mesh
/// defines the graft's shape and pose; on execution every *empty*
/// (`NaN`) voxel whose centre lies inside the mesh is filled with
/// `value`, and the measurement reports the filled volume and the
/// graft–bone interface area. Existing (bone) voxels inside the mesh are
/// left untouched — a graft fills a defect, it does not resect.
#[derive(Debug, Clone)]
pub struct GraftReconstruction {
    /// The graft's closed contour surface, pose baked in. Same mesh
    /// contract as [`ImplantPlacement`].
    pub triangles: Vec<[Vec3; 3]>,
    /// Audit label.
    pub fragment_name: String,
    /// Scalar value written into grafted voxels.
    pub value: f64,
}

impl GraftReconstruction {
    /// Mesh structural check.
    pub fn validate(&self) -> Result<(), MeshCutError> {
        validate_closed_mesh(&self.triangles)
    }
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
    /// Freeform (anatomically contoured) resection against a closed
    /// triangle mesh.
    Mesh(MeshCut),
    /// Implant component placement (bone replaced by the component).
    Implant(ImplantPlacement),
    /// Bone graft / defect reconstruction (empty voxels filled).
    Graft(GraftReconstruction),
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
    /// A [`MeshCut`]'s contour mesh failed its structural check (open
    /// surface, inconsistent winding, degenerate or missing triangles).
    InvalidMesh {
        /// What the mesh check found.
        reason: String,
    },
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
            PlanError::InvalidMesh { reason } => {
                write!(f, "contour-mesh cut rejected: {reason}")
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
    /// Implant placement measurements.
    Implant(ImplantMeasurement),
    /// Graft reconstruction measurements.
    Graft(GraftMeasurement),
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
    /// Structures-at-risk screening: one entry per watched structure
    /// × per cut step, in (structure, step) order. Empty when no
    /// structures were registered (see
    /// [`VirtualSurgery::watch_structures`]).
    pub structures_at_risk: Vec<StructureRisk>,
}

/// A **soft-tissue structure at risk**: a capsule (a segment with a
/// radius) standing in for a neurovascular bundle, ligament or tendon,
/// registered on the plan and screened against every cut surface. The
/// capsule is a screening stand-in, not deformable soft tissue — the
/// question answered is *does the plan's cut field come within the
/// structure's envelope*, not how the tissue moves.
#[derive(Debug, Clone)]
pub struct SoftTissueStructure {
    /// Audit label (e.g. `"common-peroneal-nerve"`).
    pub name: String,
    /// Capsule segment start (patient space).
    pub start: Vec3,
    /// Capsule segment end (patient space).
    pub end: Vec3,
    /// Capsule radius (mm) — the structure's envelope.
    pub radius_mm: f64,
}

/// One watched structure's clearance against one cut: the smallest
/// signed distance from the capsule envelope to the cut's removed
/// region, negative when the envelope is breached.
#[derive(Debug, Clone)]
pub struct StructureRisk {
    /// The structure's audit label.
    pub structure: String,
    /// The cut step's audit label.
    pub cut: String,
    /// Clearance (mm): min over capsule samples of (signed distance to
    /// the removed region) − envelope radius. ≤ 0 means the cut
    /// intersects the structure's envelope.
    pub clearance_mm: f64,
    /// True when the cut breaches the envelope.
    pub breached: bool,
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
    /// Soft-tissue structures watched for collateral damage.
    watched: Vec<SoftTissueStructure>,
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
            watched: Vec::new(),
        }
    }

    /// Registers soft-tissue structures (capsule stand-ins) to screen
    /// against every cut surface at execution time; the verdicts land in
    /// [`SurgeryReport::structures_at_risk`].
    pub fn watch_structures(
        &mut self,
        structures: impl IntoIterator<Item = SoftTissueStructure>,
    ) -> &mut Self {
        self.watched.extend(structures);
        self
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

    /// Appends a freeform (anatomically contoured) resection step against
    /// a closed triangle mesh.
    ///
    /// # Errors
    ///
    /// As [`Self::cut`] (including retention via
    /// [`DiscardedSide::RetainAs`]), plus
    /// [`PlanError::InvalidMesh`] when the contour mesh fails
    /// [`MeshCut::validate`].
    pub fn mesh(&mut self, cut: MeshCut) -> Result<&mut Self, PlanError> {
        self.validate_cut()?;
        if let Err(e) = cut.validate() {
            return Err(PlanError::InvalidMesh {
                reason: e.to_string(),
            });
        }
        if let DiscardedSide::RetainAs { name } = &cut.discarded {
            if name.trim().is_empty() {
                return Err(PlanError::EmptyRetainedName);
            }
        }
        let retained = matches!(cut.discarded, DiscardedSide::RetainAs { .. });
        self.fragment_log
            .push(format!("mesh:{}", cut.fragment_name));
        self.steps.push(PlanStep::Mesh(cut));
        self.fragment_names = vec![self.latest_kept_name().to_string()];
        if retained {
            if let Some(PlanStep::Mesh(c)) = self.steps.last() {
                if let DiscardedSide::RetainAs { name } = &c.discarded {
                    self.fragment_names.push(name.clone());
                }
            }
        }
        self.fragment_count = self.fragment_names.len();
        Ok(self)
    }

    /// Appends an **implant component placement** step.
    ///
    /// # Errors
    ///
    /// As [`Self::cut`], plus [`PlanError::InvalidMesh`] when the
    /// component's contour mesh fails [`MeshCut::validate`].
    pub fn place_implant(&mut self, implant: ImplantPlacement) -> Result<&mut Self, PlanError> {
        self.validate_cut()?;
        if let Err(e) = implant.validate() {
            return Err(PlanError::InvalidMesh {
                reason: e.to_string(),
            });
        }
        self.fragment_log
            .push(format!("implant:{}", implant.fragment_name));
        self.steps.push(PlanStep::Implant(implant));
        Ok(self)
    }

    /// Appends a **bone graft / defect reconstruction** step.
    ///
    /// # Errors
    ///
    /// As [`Self::cut`], plus [`PlanError::InvalidMesh`] when the
    /// graft's contour mesh fails [`MeshCut::validate`].
    pub fn add_graft(&mut self, graft: GraftReconstruction) -> Result<&mut Self, PlanError> {
        self.validate_cut()?;
        if let Err(e) = graft.validate() {
            return Err(PlanError::InvalidMesh {
                reason: e.to_string(),
            });
        }
        self.fragment_log
            .push(format!("graft:{}", graft.fragment_name));
        self.steps.push(PlanStep::Graft(graft));
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
                PlanStep::Mesh(m) => return &m.fragment_name,
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
                PlanStep::Mesh(cut) => {
                    // Validated at build: exactly one fragment, valid mesh.
                    let current = &fragments[0].model;
                    let (kept, retained, m) = cut.apply_measured(current);
                    total_resection += m.resection_volume_mm3;
                    log.push(format!("mesh:{}", cut.fragment_name));
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
                PlanStep::Implant(implant) => {
                    // Validated at build: exactly one fragment, valid mesh.
                    let current = &fragments[0].model;
                    let mut next = current.clone();
                    let (nx, ny, nz) = current.dims;
                    let (sx, sy, sz) = current.spacing;
                    // Pass 1: classify the implant region once, so the
                    // interface count cannot depend on the scan order.
                    let mut implant_mask = vec![false; current.values.len()];
                    for z in 0..nz {
                        for y in 0..ny {
                            for x in 0..nx {
                                let Some(idx) = current.index(x, y, z) else {
                                    continue;
                                };
                                let centre = current.center(x, y, z);
                                implant_mask[idx] =
                                    mesh_signed_distance(centre, &implant.triangles) < 0.0;
                            }
                        }
                    }
                    let mut resected = 0usize;
                    let mut implant_voxels = 0usize;
                    let mut interface_faces = 0usize;
                    for z in 0..nz {
                        for y in 0..ny {
                            for x in 0..nx {
                                let Some(idx) = current.index(x, y, z) else {
                                    continue;
                                };
                                if !implant_mask[idx] {
                                    continue;
                                }
                                if !current.values[idx].is_nan() {
                                    resected += 1;
                                }
                                implant_voxels += 1;
                                next.values[idx] = implant.marker_value;
                                // Interface: shared faces with surviving
                                // (non-implant, non-empty) neighbours,
                                // classified from the finished mask.
                                let neigh = [
                                    current.index(x + 1, y, z),
                                    x.checked_sub(1).and_then(|v| current.index(v, y, z)),
                                    current.index(x, y + 1, z),
                                    y.checked_sub(1).and_then(|v| current.index(x, v, z)),
                                    current.index(x, y, z + 1),
                                    z.checked_sub(1).and_then(|v| current.index(x, y, v)),
                                ];
                                for n in neigh.iter().flatten() {
                                    if !implant_mask[*n] && !current.values[*n].is_nan() {
                                        interface_faces += 1;
                                    }
                                }
                            }
                        }
                    }
                    let voxel_volume = sx * sy * sz;
                    // Mean voxel face area: each of the 6 faces carries a
                    // different pairing of spacings; a uniform grid has
                    // h². Counting per-axis would be exact but the mean
                    // keeps the bookkeeping honest for anisotropic grids
                    // too (documented; the tests use uniform pitch).
                    let face_area = (sx * sy + sy * sz + sx * sz) / 3.0;
                    let m = ImplantMeasurement {
                        fragment_name: implant.fragment_name.clone(),
                        resection_volume_mm3: resected as f64 * voxel_volume,
                        implant_volume_mm3: implant_voxels as f64 * voxel_volume,
                        interface_area_mm2: interface_faces as f64 * face_area,
                    };
                    total_resection += m.resection_volume_mm3;
                    log.push(format!("implant:{}", implant.fragment_name));
                    measurements.push(StepMeasurement::Implant(m));
                    fragments[0] = NamedFragment {
                        name: fragments[0].name.clone(),
                        model: next,
                    };
                }
                PlanStep::Graft(graft) => {
                    let current = &fragments[0].model;
                    let mut next = current.clone();
                    let (nx, ny, nz) = current.dims;
                    let (sx, sy, sz) = current.spacing;
                    // Pass 1: classify the graft region once (a graft
                    // fills empty voxels only).
                    let mut graft_mask = vec![false; current.values.len()];
                    for z in 0..nz {
                        for y in 0..ny {
                            for x in 0..nx {
                                let Some(idx) = current.index(x, y, z) else {
                                    continue;
                                };
                                if current.values[idx].is_nan() {
                                    let centre = current.center(x, y, z);
                                    graft_mask[idx] =
                                        mesh_signed_distance(centre, &graft.triangles) < 0.0;
                                }
                            }
                        }
                    }
                    let mut filled = 0usize;
                    let mut interface_faces = 0usize;
                    for z in 0..nz {
                        for y in 0..ny {
                            for x in 0..nx {
                                let Some(idx) = current.index(x, y, z) else {
                                    continue;
                                };
                                if !graft_mask[idx] {
                                    continue;
                                }
                                filled += 1;
                                next.values[idx] = graft.value;
                                let neigh = [
                                    current.index(x + 1, y, z),
                                    x.checked_sub(1).and_then(|v| current.index(v, y, z)),
                                    current.index(x, y + 1, z),
                                    y.checked_sub(1).and_then(|v| current.index(x, v, z)),
                                    current.index(x, y, z + 1),
                                    z.checked_sub(1).and_then(|v| current.index(x, y, v)),
                                ];
                                for n in neigh.iter().flatten() {
                                    if !graft_mask[*n] && !current.values[*n].is_nan() {
                                        interface_faces += 1;
                                    }
                                }
                            }
                        }
                    }
                    let voxel_volume = sx * sy * sz;
                    let face_area = (sx * sy + sy * sz + sx * sz) / 3.0;
                    let m = GraftMeasurement {
                        fragment_name: graft.fragment_name.clone(),
                        graft_volume_mm3: filled as f64 * voxel_volume,
                        interface_area_mm2: interface_faces as f64 * face_area,
                    };
                    log.push(format!("graft:{}", graft.fragment_name));
                    measurements.push(StepMeasurement::Graft(m));
                    fragments[0] = NamedFragment {
                        name: fragments[0].name.clone(),
                        model: next,
                    };
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
        let structures_at_risk = self.screen_structures(&self.steps, &self.watched);
        (
            model,
            SurgeryReport {
                step_log: log,
                measurements,
                total_resection_volume_mm3: total_resection,
                structures_at_risk,
            },
        )
    }

    /// The pre-operative base model (for side-by-side planning views).
    pub fn base_model(&self) -> &VoxelModel {
        &self.base
    }

    /// Screens the watched structures against every cut surface in the
    /// plan. The signed distance to a cut's removed region is closed
    /// form for planes, wedges and cylinders; mesh cuts sample the
    /// structure segment at a stride of at most half the capsule radius
    /// (documented sampling error, no closed-form point-mesh distance
    /// along a segment).
    fn screen_structures(
        &self,
        steps: &[PlanStep],
        structures: &[SoftTissueStructure],
    ) -> Vec<StructureRisk> {
        // Collect (label, removed-region signed distance) per cut step.
        let mut cuts: Vec<(String, Box<dyn Fn(Vec3) -> f64>)> = Vec::new();
        for step in steps {
            match step {
                PlanStep::Cut(c) => {
                    let plane = c.plane;
                    let keep_positive = c.keep_positive;
                    let half = c.kerf_width * 0.5;
                    cuts.push((
                        format!("cut:{}", c.fragment_name),
                        Box::new(move |p: Vec3| {
                            let d = plane.signed_distance(p);
                            if keep_positive {
                                d + half
                            } else {
                                -d + half
                            }
                        }),
                    ));
                }
                PlanStep::Wedge(w) => {
                    let (a, b) = (w.plane_a, w.plane_b);
                    let half = w.kerf_width * 0.5;
                    cuts.push((
                        format!("wedge:{}", w.fragment_name),
                        Box::new(move |p: Vec3| {
                            (a.signed_distance(p) + half).max(b.signed_distance(p) + half)
                        }),
                    ));
                }
                PlanStep::Cylinder(c) => {
                    let axis_origin = c.axis_origin;
                    let dir = c.axis_direction.normalize();
                    let radius = c.radius;
                    let half = c.kerf_width * 0.5;
                    let keep_inside = c.keep_inside;
                    cuts.push((
                        format!("cylinder:{}", c.fragment_name),
                        Box::new(move |p: Vec3| {
                            let w = p - axis_origin;
                            let along = w.dot(dir);
                            let radial = (w - dir * along).norm();
                            if keep_inside {
                                // Removed region: outside the wall.
                                (radius + half) - radial
                            } else {
                                // Removed region: the core.
                                radial - (radius - half)
                            }
                        }),
                    ));
                }
                PlanStep::Mesh(c) => {
                    let triangles = c.triangles.clone();
                    let half = c.kerf_width * 0.5;
                    let keep_inside = c.keep_inside;
                    cuts.push((
                        format!("mesh:{}", c.fragment_name),
                        Box::new(move |p: Vec3| {
                            let d = mesh_signed_distance(p, &triangles);
                            if keep_inside {
                                // Removed: the exterior beyond the kerf
                                // face (d > −half).
                                -(d + half)
                            } else {
                                // Removed: the interior inside the kerf
                                // face (d < +half).
                                d - half
                            }
                        }),
                    ));
                }
                _ => continue,
            }
        }
        let mut risks = Vec::new();
        for structure in structures {
            // Sample the capsule axis at a stride of at most half the
            // envelope radius (bounds the sampling error of the mesh
            // path; the closed-form paths are exact at every sample).
            let len = (structure.end - structure.start).norm();
            let stride = (structure.radius_mm.max(0.5)) * 0.5;
            let samples = (1usize + (len / stride).ceil() as usize).max(2);
            for (label, sdf) in &cuts {
                let mut min_signed = f64::INFINITY;
                for s in 0..samples {
                    let t = s as f64 / (samples - 1) as f64;
                    let p = structure.start * (1.0 - t) + structure.end * t;
                    min_signed = min_signed.min(sdf(p));
                }
                let clearance = min_signed - structure.radius_mm;
                risks.push(StructureRisk {
                    structure: structure.name.clone(),
                    cut: label.clone(),
                    clearance_mm: clearance,
                    breached: clearance <= 0.0,
                });
            }
        }
        risks
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

    /// A closed, consistently-wound box surface: 12 outward-facing
    /// triangles between `lo` and `hi`.
    fn box_mesh(lo: Vec3, hi: Vec3) -> Vec<[Vec3; 3]> {
        let p = [
            Vec3::new(lo.x, lo.y, lo.z),
            Vec3::new(hi.x, lo.y, lo.z),
            Vec3::new(hi.x, hi.y, lo.z),
            Vec3::new(lo.x, hi.y, lo.z),
            Vec3::new(lo.x, lo.y, hi.z),
            Vec3::new(hi.x, lo.y, hi.z),
            Vec3::new(hi.x, hi.y, hi.z),
            Vec3::new(lo.x, hi.y, hi.z),
        ];
        vec![
            [p[0], p[4], p[7]],
            [p[0], p[7], p[3]], // −x
            [p[5], p[1], p[2]],
            [p[5], p[2], p[6]], // +x
            [p[0], p[1], p[5]],
            [p[0], p[5], p[4]], // −y
            [p[3], p[7], p[6]],
            [p[3], p[6], p[2]], // +y
            [p[0], p[3], p[2]],
            [p[0], p[2], p[1]], // −z
            [p[4], p[5], p[6]],
            [p[4], p[6], p[7]], // +z
        ]
    }

    /// A closed icosphere (subdivided icosahedron), outward-facing.
    fn icosphere(center: Vec3, radius: f64, subdivisions: usize) -> Vec<[Vec3; 3]> {
        let t = (1.0 + 5.0f64.sqrt()) / 2.0;
        let mut v: Vec<Vec3> = [
            [-1.0, t, 0.0],
            [1.0, t, 0.0],
            [-1.0, -t, 0.0],
            [1.0, -t, 0.0],
            [0.0, -1.0, t],
            [0.0, 1.0, t],
            [0.0, -1.0, -t],
            [0.0, 1.0, -t],
            [t, 0.0, -1.0],
            [t, 0.0, 1.0],
            [-t, 0.0, -1.0],
            [-t, 0.0, 1.0],
        ]
        .iter()
        .map(|a| Vec3::new(a[0], a[1], a[2]))
        .collect();
        let mut faces: Vec<[usize; 3]> = vec![
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ];
        for _ in 0..subdivisions {
            let mut next = Vec::with_capacity(faces.len() * 4);
            let mut cache: std::collections::BTreeMap<(usize, usize), usize> =
                std::collections::BTreeMap::new();
            let mut midpoint = |a: usize, b: usize| -> usize {
                let key = (a.min(b), a.max(b));
                if let Some(&m) = cache.get(&key) {
                    return m;
                }
                v.push((v[a] + v[b]) * 0.5);
                let m = v.len() - 1;
                cache.insert(key, m);
                m
            };
            for f in &faces {
                let a = midpoint(f[0], f[1]);
                let b = midpoint(f[1], f[2]);
                let c = midpoint(f[2], f[0]);
                next.push([f[0], a, c]);
                next.push([f[1], b, a]);
                next.push([f[2], c, b]);
                next.push([a, b, c]);
            }
            faces = next;
        }
        faces
            .iter()
            .map(|f| {
                f.map(|i| {
                    let d = v[i] * (1.0 / v[i].norm());
                    center + d * radius
                })
            })
            .collect()
    }

    #[test]
    fn mesh_validate_catches_open_flipped_and_degenerate_meshes() {
        let ok = MeshCut {
            triangles: box_mesh(Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)),
            fragment_name: "m".into(),
            keep_inside: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        };
        assert!(ok.validate().is_ok());

        let empty = MeshCut {
            triangles: vec![],
            ..ok.clone()
        };
        assert!(matches!(empty.validate(), Err(MeshCutError::Empty)));

        // Degenerate: two vertices coincide.
        let degenerate = MeshCut {
            triangles: vec![[Vec3::ZERO, Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)]],
            ..ok.clone()
        };
        assert!(matches!(
            degenerate.validate(),
            Err(MeshCutError::DegenerateTriangle { index: 0 })
        ));

        // Open: a single quad (its boundary edges appear once).
        let open = MeshCut {
            triangles: vec![[
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
            ]],
            ..ok.clone()
        };
        assert!(matches!(
            open.validate(),
            Err(MeshCutError::OpenSurface { .. })
        ));

        // Flipped: reverse one triangle of the closed box — one directed
        // edge doubles up.
        let mut flipped_triangles = ok.triangles.clone();
        flipped_triangles[0].reverse();
        let flipped = MeshCut {
            triangles: flipped_triangles,
            ..ok.clone()
        };
        assert!(matches!(
            flipped.validate(),
            Err(MeshCutError::InconsistentWinding { .. })
        ));
    }

    #[test]
    fn point_triangle_distance_has_hand_checked_regions() {
        let a = Vec3::ZERO;
        let b = Vec3::new(4.0, 0.0, 0.0);
        let c = Vec3::new(0.0, 3.0, 0.0);
        // Face interior (1, 1, 5): closest point (1, 1, 0) → 25.
        assert!((point_triangle_distance2(Vec3::new(1.0, 1.0, 5.0), a, b, c) - 25.0).abs() < 1e-12);
        // Vertex region at a: (−1, −1, 0) → 2.
        assert!(
            (point_triangle_distance2(Vec3::new(-1.0, -1.0, 0.0), a, b, c) - 2.0).abs() < 1e-12
        );
        // Edge ab interior: (2, −2, 0) → closest (2, 0, 0) → 4.
        assert!((point_triangle_distance2(Vec3::new(2.0, -2.0, 0.0), a, b, c) - 4.0).abs() < 1e-12);
        // Hypotenuse edge bc: (3, 4, 0) → closest on x+y=4 → (1.8, 2.2)? The
        // projection of (3,4) onto the segment from (4,0) to (0,3): the
        // closest point is (1.96, 1.52)... verified numerically below.
        let d2 = point_triangle_distance2(Vec3::new(3.0, 4.0, 0.0), a, b, c);
        let mut best = f64::INFINITY;
        for i in 0..=1000 {
            let t = i as f64 / 1000.0;
            let q = b * (1.0 - t) + c * t;
            best = best.min((Vec3::new(3.0, 4.0, 0.0) - q).norm_squared());
        }
        assert!((d2 - best).abs() < 1e-3, "{d2} vs {best}");
    }

    #[test]
    fn mesh_cut_keeps_a_box_interior_exactly() {
        let model = cube_model();
        let cut = MeshCut {
            // Off the exact diagonal planes of the face triangulation, so
            // no +x ray from a voxel centre grazes a shared edge (the
            // documented alignment limit of the parity test).
            triangles: box_mesh(
                Vec3::new(-3.169, -3.15, -3.111),
                Vec3::new(2.231, 2.25, 2.208),
            ),
            fragment_name: "contoured".into(),
            keep_inside: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        };
        let kept = cut.apply(&model);
        let count = kept.values.iter().filter(|v| !v.is_nan()).count();
        // Voxel centres at −3..2 in each axis: exactly 6³ lie inside.
        assert_eq!(count, 216);
    }

    #[test]
    fn mesh_cut_exterior_with_retention_conserves_voxels() {
        let model = cube_model();
        let cut = MeshCut {
            triangles: box_mesh(
                Vec3::new(-3.169, -3.15, -3.111),
                Vec3::new(2.231, 2.25, 2.208),
            ),
            fragment_name: "cap".into(),
            keep_inside: false,
            kerf_width: 0.0,
            discarded: DiscardedSide::RetainAs {
                name: "contoured-core".into(),
            },
        };
        let (kept, retained, m) = cut.apply_measured(&model);
        let kept_n = kept.values.iter().filter(|v| !v.is_nan()).count();
        let core = retained.expect("retained side kept");
        let core_n = core.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(kept_n + core_n, 1000);
        // The resection measurement counts nothing: everything was kept
        // somewhere.
        assert_eq!(m.resection_volume_mm3, 0.0);
        assert_eq!(core_n, 216);
    }

    #[test]
    fn mesh_cut_volume_and_kerf_match_the_continuum_on_a_fine_grid() {
        let mut model = VoxelModel {
            dims: (24, 24, 24),
            spacing: (0.25, 0.25, 0.25),
            origin: Vec3::new(-3.0, -3.0, -3.0),
            values: vec![100.0; 24 * 24 * 24],
        };
        let _ = &mut model;
        let make_cut = |kerf: f64| MeshCut {
            triangles: box_mesh(
                Vec3::new(-2.031, -2.012, -2.05),
                Vec3::new(1.969, 1.981, 1.95),
            ),
            fragment_name: "b".into(),
            keep_inside: true,
            kerf_width: kerf,
            discarded: DiscardedSide::Resect,
        };
        // Centres at ±(0.125 + 0.25k): no centre lies on either side of a
        // kerf-shifted face, so both counts are exact.
        let plain = make_cut(0.0).apply(&model);
        let plain_n = plain.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(plain_n, 16 * 16 * 16);
        let plain_vol = plain_n as f64 * 0.25_f64.powi(3);
        assert!((plain_vol - 64.0).abs() < 1e-9);

        let kerfed = make_cut(0.5).apply(&model);
        let kerfed_n = kerfed.values.iter().filter(|v| !v.is_nan()).count();
        // Kept face: 2.0 − 0.25 = 1.75 per side → 14 centres per axis.
        assert_eq!(kerfed_n, 14 * 14 * 14);
        let _ = &mut model;
    }

    #[test]
    fn mesh_cut_sphere_volume_matches_analytic() {
        // 80-face icosphere, radius 3, on a 0.3 grid. The voxelisation
        // error is bounded by the surface shell (A·h ≈ 34 mm³ worst case);
        // assert the midpoint statistics stay well inside that.
        let h = 0.3;
        let n = 34;
        let mut model = VoxelModel {
            dims: (n, n, n),
            spacing: (h, h, h),
            origin: Vec3::new(0.1, 0.1, 0.1),
            values: vec![100.0; n * n * n],
        };
        let _ = &mut model;
        let cut = MeshCut {
            triangles: icosphere(Vec3::new(5.0, 5.0, 5.0), 3.0, 1),
            fragment_name: "sphere".into(),
            keep_inside: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        };
        cut.validate().expect("icosphere is closed and wound");
        let kept = cut.apply(&model);
        let vol = kept.values.iter().filter(|v| !v.is_nan()).count() as f64 * h.powi(3);
        // The reference is the MESH's own enclosed volume — the signed
        // tetrahedron sum over the (outward-wound) faces — not the ideal
        // sphere: an 80-face icosphere's face planes sit apothem-deep
        // inside the sphere, a few percent below (4/3)πR³.
        let mesh_volume = (1.0 / 6.0
            * cut
                .triangles
                .iter()
                .map(|t| t[0].dot(t[1].cross(t[2])))
                .sum::<f64>())
        .abs();
        let sphere = 4.0 / 3.0 * core::f64::consts::PI * 27.0;
        assert!(
            (mesh_volume - sphere).abs() < 0.15 * sphere,
            "subdiv-1 icosphere sanity: {mesh_volume} vs {sphere}"
        );
        // Voxelisation error is bounded by the surface shell (A·h ≈ 34
        // mm³ worst case); 5 % is comfortably inside the midpoint
        // statistics.
        assert!(
            (vol - mesh_volume).abs() < 0.05 * mesh_volume,
            "voxelised {vol} vs mesh {mesh_volume}"
        );
        let _ = &mut model;
    }

    #[test]
    fn virtual_surgery_runs_a_mesh_cut_and_rejects_bad_meshes() {
        let mut plan = VirtualSurgery::new(cube_model());
        // An open surface is rejected at build time.
        let open = MeshCut {
            triangles: vec![[
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
            ]],
            fragment_name: "open".into(),
            keep_inside: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        };
        assert!(matches!(
            plan.mesh(open),
            Err(PlanError::InvalidMesh { .. })
        ));

        let valid = MeshCut {
            triangles: box_mesh(
                Vec3::new(-3.169, -3.15, -3.111),
                Vec3::new(2.231, 2.25, 2.208),
            ),
            fragment_name: "contour".into(),
            keep_inside: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        };
        plan.mesh(valid).expect("valid mesh cut");
        let (model, report) = plan.execute_with_report();
        assert!(
            report
                .step_log
                .iter()
                .any(|s| s.starts_with("mesh:contour")),
            "audit log records the mesh step"
        );
        let kept = model.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(kept, 216);

        // A mesh cut after a named move is rejected like any cut.
        let mut plan = VirtualSurgery::new(cube_model());
        plan.move_fragment_named("", FragmentTransform::no_op())
            .expect("named move");
        let later = MeshCut {
            triangles: box_mesh(
                Vec3::new(-3.169, -3.15, -3.111),
                Vec3::new(2.231, 2.25, 2.208),
            ),
            fragment_name: "late".into(),
            keep_inside: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        };
        assert!(matches!(
            plan.mesh(later),
            Err(PlanError::CutAfterMove { .. })
        ));
    }

    #[test]
    fn implant_placement_resects_and_measures_the_interface() {
        // A 2×2×2 component embedded in the bone block: every outward
        // face of the corner block touches bone → 24 interface faces.
        let mut plan = VirtualSurgery::new(cube_model());
        let implant = ImplantPlacement {
            triangles: box_mesh(Vec3::new(-0.6, -0.7, -0.55), Vec3::new(1.4, 1.3, 1.45)),
            fragment_name: "tka-femoral".into(),
            marker_value: 2000.0,
        };
        implant.validate().expect("valid component");
        plan.place_implant(implant).expect("valid placement");
        let (model, report) = plan.execute_with_report();
        let marker_count = model.values.iter().filter(|&&v| v == 2000.0).count();
        assert_eq!(marker_count, 8, "implant occupies its 2×2×2 region");
        let StepMeasurement::Implant(m) = &report.measurements[0] else {
            panic!("implant measurement expected");
        };
        assert_eq!(m.fragment_name, "tka-femoral");
        assert!((m.resection_volume_mm3 - 8.0).abs() < 1e-9);
        assert!((m.implant_volume_mm3 - 8.0).abs() < 1e-9);
        assert!((m.interface_area_mm2 - 24.0).abs() < 1e-9);
        assert!((report.total_resection_volume_mm3 - 8.0).abs() < 1e-9);
        assert!(report.step_log[0].starts_with("implant:tka-femoral"));

        // Sequencing: an open component mesh is rejected at build time,
        // and placement after a named move follows the cut rules.
        let mut plan = VirtualSurgery::new(cube_model());
        let open = ImplantPlacement {
            triangles: vec![[
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
            ]],
            fragment_name: "open".into(),
            marker_value: 2000.0,
        };
        assert!(matches!(
            plan.place_implant(open),
            Err(PlanError::InvalidMesh { .. })
        ));
        plan.move_fragment_named("", FragmentTransform::no_op())
            .expect("named move");
        let later = ImplantPlacement {
            triangles: box_mesh(Vec3::new(-0.6, -0.7, -0.55), Vec3::new(1.4, 1.3, 1.45)),
            fragment_name: "late".into(),
            marker_value: 2000.0,
        };
        assert!(matches!(
            plan.place_implant(later),
            Err(PlanError::CutAfterMove { .. })
        ));
    }

    #[test]
    fn graft_fills_a_contoured_defect_and_conserves_the_model() {
        // A contoured cavity is cut out of the block (exterior kept, so
        // the compacted grid still surrounds the NaN pocket), then the
        // same surface grafts it back: the voxel count is restored.
        let surface = box_mesh(Vec3::new(-0.6, -0.7, -0.55), Vec3::new(1.4, 1.3, 1.45));
        let mut plan = VirtualSurgery::new(cube_model());
        plan.mesh(MeshCut {
            triangles: surface.clone(),
            fragment_name: "defect".into(),
            keep_inside: false,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        })
        .expect("cavity cut");
        let (defected, _) = plan.execute();
        let defect_count = defected.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(defect_count, 1000 - 8);

        let mut plan = VirtualSurgery::new(cube_model());
        plan.mesh(MeshCut {
            triangles: surface.clone(),
            fragment_name: "defect".into(),
            keep_inside: false,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        })
        .expect("cavity cut");
        plan.add_graft(GraftReconstruction {
            triangles: surface,
            fragment_name: "ibg".into(),
            value: 150.0,
        })
        .expect("graft");
        let (reconstructed, report) = plan.execute_with_report();
        let restored = reconstructed.values.iter().filter(|v| !v.is_nan()).count();
        assert_eq!(restored, 1000, "graft restores the resected voxels");
        let StepMeasurement::Graft(m) = &report.measurements[1] else {
            panic!("graft measurement expected");
        };
        assert_eq!(m.fragment_name, "ibg");
        assert!((m.graft_volume_mm3 - 8.0).abs() < 1e-9);
        // Each of the 8 cavity voxels faces bone on its 3 outward sides.
        assert!((m.interface_area_mm2 - 24.0).abs() < 1e-9);
        assert!(report.step_log[1].starts_with("graft:ibg"));
    }

    fn structure(name: &str, start: Vec3, end: Vec3, radius: f64) -> SoftTissueStructure {
        SoftTissueStructure {
            name: name.into(),
            start,
            end,
            radius_mm: radius,
        }
    }

    #[test]
    fn plane_cut_screen_breaches_and_clears_exactly() {
        // A plane at z = 0 keeping the positive side. A nerve crossing
        // the plane is breached; one parallel at a known distance is
        // cleared by exactly that distance minus the envelope radius
        // (the plane path is closed form along the whole segment, so no
        // sampling slop).
        let mut plan = VirtualSurgery::new(cube_model());
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).expect("plane"),
            fragment_name: "distal".into(),
            keep_positive: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        })
        .expect("cut");
        plan.watch_structures([
            structure(
                "crossing",
                Vec3::new(0.0, 0.0, -2.0),
                Vec3::new(0.0, 0.0, 2.0),
                0.5,
            ),
            structure(
                "parallel",
                Vec3::new(2.0, 0.0, 1.0),
                Vec3::new(2.0, 0.0, 3.0),
                0.5,
            ),
        ]);
        let (_, report) = plan.execute_with_report();
        assert_eq!(report.structures_at_risk.len(), 2);
        let crossing = &report.structures_at_risk[0];
        assert_eq!(crossing.structure, "crossing");
        assert!(crossing.breached, "a crossing structure must breach");
        // The segment's far end sits 2 mm inside the removed side; the
        // minimum over samples is exactly at the end: clearance =
        // −(2 + 0.5).
        assert!(
            (crossing.clearance_mm + 2.5).abs() < 1e-9,
            "{}",
            crossing.clearance_mm
        );
        let parallel = &report.structures_at_risk[1];
        assert!(!parallel.breached);
        // Both endpoints on the kept side, the nearer 1 mm from the
        // plane: clearance = 1 − 0.5.
        assert!(
            (parallel.clearance_mm - 0.5).abs() < 1e-9,
            "{}",
            parallel.clearance_mm
        );
    }

    #[test]
    fn mesh_and_cylinder_cuts_screen_their_removed_regions() {
        // A box-contoured resection (interior removed) with a nerve
        // through the middle: breached. A cylinder cut keeping its core
        // with a nerve outside the wall: cleared by the radial margin.
        let mut plan = VirtualSurgery::new(cube_model());
        plan.mesh(MeshCut {
            triangles: box_mesh(Vec3::new(-0.6, -0.7, -0.55), Vec3::new(1.4, 1.3, 1.45)),
            fragment_name: "contour".into(),
            keep_inside: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        })
        .expect("mesh cut");
        plan.cylinder(CylindricalCut {
            axis_origin: Vec3::new(0.0, 0.0, 0.0),
            axis_direction: Vec3::Z,
            radius: 2.0,
            fragment_name: "core".into(),
            keep_inside: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        })
        .expect("cylinder cut");
        plan.watch_structures([
            structure(
                "in-box",
                Vec3::new(0.0, 0.0, -0.4),
                Vec3::new(1.0, 0.5, 1.2),
                0.3,
            ),
            structure(
                "in-core",
                Vec3::new(1.2, 0.0, -1.0),
                Vec3::new(1.2, 0.0, 1.0),
                0.5,
            ),
        ]);
        let (_, report) = plan.execute_with_report();
        // (structure, cut) order: in-box×mesh, in-box×cylinder,
        // in-core×mesh, in-core×cylinder.
        assert_eq!(report.structures_at_risk.len(), 4);
        let in_box_mesh = &report.structures_at_risk[0];
        assert!(in_box_mesh.breached);
        assert!(in_box_mesh.cut.starts_with("mesh:"));
        // The nerve is inside the KEPT contour, but its envelope
        // (0.3 mm) pokes through the nearest face (the shallowest
        // sample sits 0.15 mm inside): clearance = 0.15 − 0.3.
        assert!(
            (in_box_mesh.clearance_mm + 0.15).abs() < 0.05,
            "{}",
            in_box_mesh.clearance_mm
        );
        // The in-core nerve sits inside the kept cylinder core (radial
        // 1.2 < 2.0): the cylinder screen clears it by the radial
        // margin minus the envelope (closed form along the segment).
        let in_core_cyl = &report.structures_at_risk[3];
        assert!(in_core_cyl.cut.starts_with("cylinder:"));
        assert!(!in_core_cyl.breached);
        assert!((in_core_cyl.clearance_mm - 0.3).abs() < 1e-9);
        // The same nerve crosses the mesh cut's removed exterior (the
        // box only spans z ∈ (−0.55, 1.45)), so the mesh screen
        // breaches it — two cuts, two verdicts, one structure.
        let in_core_mesh = &report.structures_at_risk[2];
        assert!(in_core_mesh.cut.starts_with("mesh:"));
        assert!(in_core_mesh.breached);
    }

    #[test]
    fn no_structures_means_no_screen_entries() {
        let mut plan = VirtualSurgery::new(cube_model());
        plan.cut(OsteotomyCut {
            plane: Plane::from_point_normal(Vec3::ZERO, Vec3::Z).expect("plane"),
            fragment_name: "distal".into(),
            keep_positive: true,
            kerf_width: 0.0,
            discarded: DiscardedSide::Resect,
        })
        .expect("cut");
        let (_, report) = plan.execute_with_report();
        assert!(report.structures_at_risk.is_empty());
    }
}
