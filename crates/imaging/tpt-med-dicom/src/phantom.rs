//! Calibration-phantom rod sampling (`rfcs/0008-phantom-rod-sampling.md`).
//!
//! Turns a CT scan of a calibration phantom into the `(HU, known_value)`
//! pairs [`crate::QctCalibration::fit`] (and, via
//! `rfcs/0007-bmd-apparent-density-conversion.md`,
//! [`crate::BmdToAshDensity`]) already consume. Splits the problem into a
//! manufacturer-agnostic half ([`locate_phantom_centroid`]: any solid
//! phantom is denser than the air/table around it, so its cross-section
//! centroid can be found by thresholding alone) and a manufacturer-specific
//! half ([`PhantomModel`]: rod layout and known values, which this crate
//! does not guess at and requires a citation for). Rotation and slice
//! selection are caller-supplied in v0 — see the RFC for why blind detection
//! of either is not attempted.

use crate::error::{DicomError, Result};
use crate::series::DicomSlice;

/// One phantom rod: its offset from the phantom's cross-section centroid (in
/// the phantom's own frame, at a stated reference angle), physical radius,
/// and known calibration value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhantomRod {
    /// Radial distance from centroid, mm.
    pub radial_offset_mm: f64,
    /// Angle from the phantom's reference axis, radians.
    pub angle_rad: f64,
    /// Rod radius, mm.
    pub radius_mm: f64,
    /// The value to pair with this rod's measured HU (density, BMD,
    /// whatever the caller's downstream calibration expects).
    pub known_value: f64,
}

/// A named, cited rod layout for one phantom model.
///
/// Cannot be constructed without a `source` citation (datasheet, revision,
/// lot if known) — the manufacturer-specific data this crate refuses to
/// bake in silently, mirroring `BmdToAshDensity`/`AshFraction`.
#[derive(Debug, Clone, PartialEq)]
pub struct PhantomModel {
    rods: Vec<PhantomRod>,
    source: String,
}

impl PhantomModel {
    /// `rods` must be non-empty; `source` must be non-empty; every rod's
    /// fields must be finite, and `radius_mm` must be positive.
    pub fn new(rods: Vec<PhantomRod>, source: impl Into<String>) -> Result<Self> {
        if rods.is_empty() {
            return Err(DicomError::Phantom(
                "a phantom model needs at least one rod".into(),
            ));
        }
        let source = source.into();
        if source.trim().is_empty() {
            return Err(DicomError::Phantom(
                "a citation (source) is required and cannot be empty — see \
                 rfcs/0008-phantom-rod-sampling.md"
                    .into(),
            ));
        }
        for rod in &rods {
            if !rod.radial_offset_mm.is_finite()
                || !rod.angle_rad.is_finite()
                || !rod.radius_mm.is_finite()
                || !rod.known_value.is_finite()
            {
                return Err(DicomError::Phantom(
                    "phantom rod fields must be finite".into(),
                ));
            }
            if rod.radius_mm <= 0.0 {
                return Err(DicomError::Phantom(format!(
                    "phantom rod radius must be positive, got {}",
                    rod.radius_mm
                )));
            }
        }
        Ok(Self { rods, source })
    }

    /// The rod layout, in the order [`sample_phantom_rods`] returns results.
    pub fn rods(&self) -> &[PhantomRod] {
        &self.rods
    }

    /// The model's citation.
    pub fn source(&self) -> &str {
        &self.source
    }
}

/// Locates a solid phantom's cross-section centroid in one slice by
/// intensity thresholding — the manufacturer-agnostic half of
/// `rfcs/0008-phantom-rod-sampling.md`.
///
/// `background_max_hu` separates phantom material from surrounding air/
/// table. Pixels strictly above it are candidate phantom material; the
/// largest 4-connected component of such pixels is taken as the phantom
/// (rejecting small, disconnected noise regions), and `None` is returned if
/// no component reaches `min_area_px` — i.e. no phantom-sized object was
/// found in this slice.
pub fn locate_phantom_centroid(
    slice: &DicomSlice,
    background_max_hu: f64,
    min_area_px: usize,
) -> Option<(f64, f64)> {
    let rows = slice.rows as usize;
    let cols = slice.columns as usize;
    if rows == 0 || cols == 0 {
        return None;
    }

    let hu = slice.hu_plane();
    let above = |idx: usize| hu[idx] > background_max_hu;

    let mut visited = vec![false; rows * cols];
    let mut best: Vec<usize> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();

    for start in 0..(rows * cols) {
        if visited[start] || !above(start) {
            continue;
        }
        let mut component = Vec::new();
        stack.clear();
        stack.push(start);
        visited[start] = true;
        while let Some(idx) = stack.pop() {
            component.push(idx);
            let r = idx / cols;
            let c = idx % cols;
            let neighbours = [
                (r > 0).then(|| idx - cols),
                (r + 1 < rows).then(|| idx + cols),
                (c > 0).then(|| idx - 1),
                (c + 1 < cols).then(|| idx + 1),
            ];
            for n in neighbours.into_iter().flatten() {
                if !visited[n] && above(n) {
                    visited[n] = true;
                    stack.push(n);
                }
            }
        }
        if component.len() > best.len() {
            best = component;
        }
    }

    if best.len() < min_area_px {
        return None;
    }

    let n = best.len() as f64;
    let (sum_r, sum_c) = best.iter().fold((0.0, 0.0), |(sr, sc), &idx| {
        (sr + (idx / cols) as f64, sc + (idx % cols) as f64)
    });
    Some((sum_r / n, sum_c / n))
}

/// Samples `model`'s rod layout across `slices`, given a located centroid
/// and the phantom's rotation, returning one `(mean_hu, known_value)` pair
/// per rod in `model.rods()` order.
///
/// `rotation_rad` is added to each rod's `angle_rad` before converting to
/// image-space offsets — `0.0` matches a phantom scanned at its documented
/// reference orientation. Each rod's ROI is a disk of radius
/// `rod.radius_mm * roi_fraction` centred on its computed image-space
/// position, deliberately smaller than the rod's true radius so
/// partial-volume pixels at the rod's own edge do not bias the mean.
/// Per-rod means are averaged uniformly across `slices`.
pub fn sample_phantom_rods(
    slices: &[DicomSlice],
    model: &PhantomModel,
    centroid: (f64, f64),
    rotation_rad: f64,
    roi_fraction: f64,
) -> Result<Vec<(f64, f64)>> {
    if slices.is_empty() {
        return Err(DicomError::Phantom("at least one slice is required".into()));
    }
    if !(roi_fraction > 0.0 && roi_fraction <= 1.0) {
        return Err(DicomError::Phantom(format!(
            "roi_fraction must be in (0.0, 1.0], got {roi_fraction}"
        )));
    }

    let first = &slices[0];
    for slice in slices {
        if slice.rows != first.rows
            || slice.columns != first.columns
            || slice.pixel_spacing != first.pixel_spacing
        {
            return Err(DicomError::Phantom(
                "all slices must share rows/columns/pixel_spacing".into(),
            ));
        }
    }

    let rows = first.rows as usize;
    let cols = first.columns as usize;
    let (row_spacing, col_spacing) = first.pixel_spacing;
    if row_spacing <= 0.0 || col_spacing <= 0.0 {
        return Err(DicomError::Phantom("pixel spacing must be positive".into()));
    }

    let mut results = Vec::with_capacity(model.rods().len());
    for rod in model.rods() {
        let theta = rod.angle_rad + rotation_rad;
        let offset_col_mm = rod.radial_offset_mm * theta.cos();
        let offset_row_mm = rod.radial_offset_mm * theta.sin();
        let centre_row = centroid.0 + offset_row_mm / row_spacing;
        let centre_col = centroid.1 + offset_col_mm / col_spacing;

        let radius_mm = rod.radius_mm * roi_fraction;
        let radius_row_px = radius_mm / row_spacing;
        let radius_col_px = radius_mm / col_spacing;

        let r_min = (centre_row - radius_row_px).floor();
        let r_max = (centre_row + radius_row_px).ceil();
        let c_min = (centre_col - radius_col_px).floor();
        let c_max = (centre_col + radius_col_px).ceil();

        if r_min < 0.0 || c_min < 0.0 || r_max >= rows as f64 || c_max >= cols as f64 {
            return Err(DicomError::Phantom(format!(
                "rod ROI at (row={centre_row:.2}, col={centre_col:.2}, r={radius_mm:.2}mm) \
                 falls outside the {rows}x{cols} slice bounds"
            )));
        }

        let r0 = r_min.max(0.0) as usize;
        let r1 = (r_max.min(rows as f64 - 1.0)) as usize;
        let c0 = c_min.max(0.0) as usize;
        let c1 = (c_max.min(cols as f64 - 1.0)) as usize;

        let mut slice_means = Vec::with_capacity(slices.len());
        for slice in slices {
            let hu = slice.hu_plane();
            let mut sum = 0.0;
            let mut count = 0usize;
            for r in r0..=r1 {
                for c in c0..=c1 {
                    let dr = (r as f64 - centre_row) / radius_row_px;
                    let dc = (c as f64 - centre_col) / radius_col_px;
                    if dr * dr + dc * dc <= 1.0 {
                        sum += hu[r * cols + c];
                        count += 1;
                    }
                }
            }
            if count == 0 {
                return Err(DicomError::Phantom(
                    "rod ROI contained no pixels (radius too small relative to pixel spacing)"
                        .into(),
                ));
            }
            slice_means.push(sum / count as f64);
        }

        let mean_hu = slice_means.iter().sum::<f64>() / slice_means.len() as f64;
        results.push((mean_hu, rod.known_value));
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn slice_with(
        rows: u16,
        cols: u16,
        spacing: (f64, f64),
        hu: impl Fn(usize, usize) -> f64,
    ) -> DicomSlice {
        let mut pixel_data = Vec::with_capacity(rows as usize * cols as usize);
        for r in 0..rows as usize {
            for c in 0..cols as usize {
                pixel_data.push(hu(r, c) as i32);
            }
        }
        DicomSlice {
            rows,
            columns: cols,
            pixel_data,
            rescale_slope: 1.0,
            rescale_intercept: 0.0,
            pixel_spacing: spacing,
            ..DicomSlice::default()
        }
    }

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn locates_a_centred_disk() {
        let slice = slice_with(40, 40, (1.0, 1.0), |r, c| {
            let dr = r as f64 - 20.0;
            let dc = c as f64 - 20.0;
            if dr * dr + dc * dc <= 100.0 {
                200.0
            } else {
                -1000.0
            }
        });
        let centroid = locate_phantom_centroid(&slice, 0.0, 50).expect("found");
        assert!(close(centroid.0, 20.0, 0.6));
        assert!(close(centroid.1, 20.0, 0.6));
    }

    #[test]
    fn locates_an_off_centre_disk() {
        let slice = slice_with(50, 50, (1.0, 1.0), |r, c| {
            let dr = r as f64 - 15.0;
            let dc = c as f64 - 35.0;
            if dr * dr + dc * dc <= 64.0 {
                300.0
            } else {
                -1000.0
            }
        });
        let centroid = locate_phantom_centroid(&slice, 0.0, 50).expect("found");
        assert!(close(centroid.0, 15.0, 0.6));
        assert!(close(centroid.1, 35.0, 0.6));
    }

    #[test]
    fn returns_none_when_nothing_reaches_min_area() {
        let slice = slice_with(20, 20, (1.0, 1.0), |_, _| -1000.0);
        assert!(locate_phantom_centroid(&slice, 0.0, 10).is_none());
    }

    #[test]
    fn rejects_small_disconnected_noise() {
        let slice = slice_with(20, 20, (1.0, 1.0), |r, c| {
            if r == 5 && c == 5 {
                500.0 // single isolated hot pixel, not a phantom
            } else {
                -1000.0
            }
        });
        assert!(locate_phantom_centroid(&slice, 0.0, 5).is_none());
    }

    fn two_rod_model() -> PhantomModel {
        PhantomModel::new(
            vec![
                PhantomRod {
                    radial_offset_mm: 10.0,
                    angle_rad: 0.0,
                    radius_mm: 3.0,
                    known_value: 100.0,
                },
                PhantomRod {
                    radial_offset_mm: 10.0,
                    angle_rad: std::f64::consts::PI,
                    radius_mm: 3.0,
                    known_value: 200.0,
                },
            ],
            "test fixture",
        )
        .expect("constructs")
    }

    #[test]
    fn samples_rods_at_known_offsets() {
        // Rod 0 at (row=30, col=40), rod 1 at (row=30, col=20) given
        // centroid (30, 30), spacing 1mm/px, rotation 0.
        let slice = slice_with(60, 60, (1.0, 1.0), |r, c| {
            let d0 = ((r as f64 - 30.0).powi(2) + (c as f64 - 40.0).powi(2)).sqrt();
            let d1 = ((r as f64 - 30.0).powi(2) + (c as f64 - 20.0).powi(2)).sqrt();
            if d0 <= 3.0 {
                111.0
            } else if d1 <= 3.0 {
                222.0
            } else {
                0.0
            }
        });
        let model = two_rod_model();
        let results =
            sample_phantom_rods(&[slice], &model, (30.0, 30.0), 0.0, 0.8).expect("samples");
        assert_eq!(results.len(), 2);
        assert!(close(results[0].0, 111.0, 1e-9));
        assert_eq!(results[0].1, 100.0);
        assert!(close(results[1].0, 222.0, 1e-9));
        assert_eq!(results[1].1, 200.0);
    }

    #[test]
    fn averages_across_multiple_slices() {
        let make = |hu0: f64, hu1: f64| {
            slice_with(60, 60, (1.0, 1.0), move |r, c| {
                let d0 = ((r as f64 - 30.0).powi(2) + (c as f64 - 40.0).powi(2)).sqrt();
                let d1 = ((r as f64 - 30.0).powi(2) + (c as f64 - 20.0).powi(2)).sqrt();
                if d0 <= 3.0 {
                    hu0
                } else if d1 <= 3.0 {
                    hu1
                } else {
                    0.0
                }
            })
        };
        let slices = vec![make(100.0, 200.0), make(120.0, 220.0)];
        let model = two_rod_model();
        let results =
            sample_phantom_rods(&slices, &model, (30.0, 30.0), 0.0, 0.8).expect("samples");
        assert!(close(results[0].0, 110.0, 1e-9));
        assert!(close(results[1].0, 210.0, 1e-9));
    }

    #[test]
    fn rejects_roi_outside_bounds() {
        let slice = slice_with(20, 20, (1.0, 1.0), |_, _| 0.0);
        let model = two_rod_model();
        // Centroid near the edge pushes rod ROIs out of bounds.
        let err = sample_phantom_rods(&[slice], &model, (1.0, 1.0), 0.0, 0.8).unwrap_err();
        assert!(matches!(err, DicomError::Phantom(_)));
    }

    #[test]
    fn rejects_empty_slices() {
        let model = two_rod_model();
        let err = sample_phantom_rods(&[], &model, (30.0, 30.0), 0.0, 0.8).unwrap_err();
        assert!(matches!(err, DicomError::Phantom(_)));
    }

    #[test]
    fn rejects_invalid_roi_fraction() {
        let slice = slice_with(60, 60, (1.0, 1.0), |_, _| 0.0);
        let model = two_rod_model();
        assert!(
            sample_phantom_rods(std::slice::from_ref(&slice), &model, (30.0, 30.0), 0.0, 0.0)
                .is_err()
        );
        assert!(sample_phantom_rods(&[slice], &model, (30.0, 30.0), 0.0, 1.5).is_err());
    }

    #[test]
    fn rejects_mismatched_slice_dims() {
        let a = slice_with(60, 60, (1.0, 1.0), |_, _| 0.0);
        let b = slice_with(40, 40, (1.0, 1.0), |_, _| 0.0);
        let model = two_rod_model();
        let err = sample_phantom_rods(&[a, b], &model, (30.0, 30.0), 0.0, 0.8).unwrap_err();
        assert!(matches!(err, DicomError::Phantom(_)));
    }

    #[test]
    fn phantom_model_rejects_empty_rods() {
        assert!(PhantomModel::new(vec![], "test").is_err());
    }

    #[test]
    fn phantom_model_rejects_empty_source() {
        let rods = vec![PhantomRod {
            radial_offset_mm: 1.0,
            angle_rad: 0.0,
            radius_mm: 1.0,
            known_value: 1.0,
        }];
        assert!(PhantomModel::new(rods, "").is_err());
    }

    #[test]
    fn phantom_model_rejects_non_positive_radius() {
        let rods = vec![PhantomRod {
            radial_offset_mm: 1.0,
            angle_rad: 0.0,
            radius_mm: 0.0,
            known_value: 1.0,
        }];
        assert!(PhantomModel::new(rods, "test").is_err());
    }
}
