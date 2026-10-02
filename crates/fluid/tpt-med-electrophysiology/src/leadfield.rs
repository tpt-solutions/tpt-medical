//! Pseudo-bidomain lead-field projection — the first slice of
//! `rfcs/0005-cardiac-electrophysiology.md` **Stage 2** (the ECG/EGM
//! forward problem).
//!
//! The extracellular potential an electrode sees is the tissue's
//! transmembrane-potential gradient projected onto the lead field of an
//! unbounded homogeneous conductor (Geselowitz 1967; the standard
//! pseudo-ECG reduction: the monodomain tissue embedded in a conductive
//! bath, no torso-boundary heterogeneity):
//!
//! ```text
//! V_e = −(σ_i/σ_e)·(V_m,peak/4π) · ∫_H ∇u_m · ∇(1/|x − x_e|) dV
//! ```
//!
//! with `u_m` the dimensionless monodomain field this crate's Stage 1
//! stepper produces, `H` the tissue, `σ_i/σ_e` the intra/extracellular
//! conductivity ratio, and `V_m,peak` the rest-to-peak transmembrane
//! amplitude that turns the dimensionless field into volts. The bath
//! carries no transmembrane source, so the integral runs over the tissue
//! only.
//!
//! Discretely the integral is evaluated **face-wise** over the faces
//! between adjacent solid voxel pairs:
//!
//! ```text
//! ∫_H ∇u_m · ∇K dV ≈ Σ_faces (u⁺ − u⁻)·(K⁺ − K⁻)·h_b·h_c/h_a
//! ```
//!
//! Each face contributes its jump weighted by the kernel difference
//! across it — the finite-volume form of the same integral, exact for
//! anisotropic spacing. There are deliberately **no bath faces**: the
//! transmembrane field does not continue into the bath, and adding a
//! fictitious `u = 0` jump at the tissue boundary would (by the
//! integration-by-parts identity `−∫∇u·∇K = −∮u∂K/∂n`) double-count the
//! surface term that the interior faces already carry. That identity is
//! also the physics this projection inherits from the reduction, and the
//! tests pin both sides of it: a **closed** polarization front (u stepping
//! and stepping back, wholly inside the tissue) integrates to ~zero —
//! its facing surfaces cancel — while an **open** front that terminates
//! at the tissue boundary produces the classical solid-angle signal with
//! a dipolar far field. Real depolarization wavefronts are open surfaces
//! (they end on the heart's own boundary), which is why the reduction
//! has a signal at all.
//!
//! # Scope, stated plainly
//!
//! This is the **source integral in an unbounded bath**: no torso
//! geometry, no conductivity heterogeneity, no electrode transfer
//! impedances — the pieces a real 12-lead comparison would need and
//! which the RFC deliberately leaves to a follow-up. What it *does*
//! give, at screening fidelity, is morphology and timing of the
//! potential at any caller-placed electrode (epicardial electrograms
//! included), which is the quantitative half of the Stage 2 question
//! that does not require a torso model.
use crate::error::{EpError, Result};
use tpt_med_geometry::Vec3;
use tpt_med_meshing::SegmentationMask;

/// A lead-field projection over a segmentation mask's voxel grid: the
/// geometry half of the pseudo-ECG integral, reusable across any number
/// of electrodes and time samples of the Stage 1 field.
#[derive(Debug, Clone)]
pub struct LeadFieldProjection {
    dims: (usize, usize, usize),
    spacing: (f64, f64, f64),
    origin: Vec3,
    row_dir: Vec3,
    col_dir: Vec3,
    slice_dir: Vec3,
    solid: Vec<bool>,
}

impl LeadFieldProjection {
    /// Builds the projection geometry from the same mask the Stage 1
    /// tissue was built from. Errs for a mask with no solid voxels (a
    /// projection onto empty tissue is meaningless, and the tissue's own
    /// `from_mask` rejects it too).
    pub fn from_mask(mask: &SegmentationMask) -> Result<Self> {
        if !mask.voxels.iter().any(|&s| s) {
            return Err(EpError::EmptyMask);
        }
        Ok(Self {
            dims: mask.dims,
            spacing: mask.spacing,
            origin: mask.origin,
            row_dir: mask.row_dir,
            col_dir: mask.col_dir,
            slice_dir: mask.slice_dir,
            solid: mask.voxels.clone(),
        })
    }

    /// Patient-space position of a voxel centre (mm): the mask origin
    /// plus one spacing step along each image basis per index — the same
    /// convention `tpt-med-meshing` uses everywhere.
    fn voxel_position(&self, x: f64, y: f64, z: f64) -> Vec3 {
        self.origin
            + self.row_dir * (x * self.spacing.0)
            + self.col_dir * (y * self.spacing.1)
            + self.slice_dir * (z * self.spacing.2)
    }

    fn index(&self, x: usize, y: usize, z: usize) -> Option<usize> {
        let (nx, ny, nz) = self.dims;
        (x < nx && y < ny && z < nz).then(|| (z * ny + y) * nx + x)
    }

    /// The discretized `∫_H ∇u_m · ∇K dV`: the face-wise sum over every
    /// face between two adjacent **solid** voxels. A face along axis `a`
    /// contributes `(u⁺ − u⁻)·(K⁺ − K⁻)·h_b·h_c/h_a` — the two
    /// face-centred derivatives times the face's volume weight, exact
    /// for anisotropic spacing too.
    fn source_integral(&self, electrode: Vec3, v: &[f64]) -> f64 {
        let (nx, ny, nz) = self.dims;
        let n = [nx, ny, nz];
        let h = [self.spacing.0, self.spacing.1, self.spacing.2];
        let k = |p: Vec3| 1.0 / (p - electrode).norm();
        let mut sum = 0.0f64;

        for axis in 0..3 {
            let (span_b, span_c) = match axis {
                0 => (1usize, 2usize),
                1 => (0usize, 2usize),
                _ => (0usize, 1usize),
            };
            let weight = h[span_b] * h[span_c] / h[axis];
            for a in 0..n[axis] {
                // Interior faces only: layer a+1 must exist on the grid.
                let Some(hi_layer) = a.checked_add(1) else {
                    break;
                };
                if hi_layer >= n[axis] {
                    break;
                }
                let mut face = 0.0f64;
                for b in 0..n[span_b] {
                    for c in 0..n[span_c] {
                        let mut lo = [0isize; 3];
                        let mut hi = [0isize; 3];
                        for (arr, layer) in [(&mut lo, a as isize), (&mut hi, hi_layer as isize)] {
                            arr[axis] = layer;
                            arr[span_b] = b as isize;
                            arr[span_c] = c as isize;
                        }
                        let lo_idx = self
                            .index(lo[0] as usize, lo[1] as usize, lo[2] as usize)
                            .expect("layer a < n, spans in-bounds");
                        let hi_idx = self
                            .index(hi[0] as usize, hi[1] as usize, hi[2] as usize)
                            .expect("layer a+1 < n, spans in-bounds");
                        if !self.solid[lo_idx] || !self.solid[hi_idx] {
                            continue;
                        }
                        let p_lo = self.voxel_position(lo[0] as f64, lo[1] as f64, lo[2] as f64);
                        let p_hi = self.voxel_position(hi[0] as f64, hi[1] as f64, hi[2] as f64);
                        face += (v[hi_idx] - v[lo_idx]) * (k(p_hi) - k(p_lo));
                    }
                }
                sum += face * weight;
            }
        }
        sum
    }

    /// The pseudo-ECG potential (mV) at `electrode` (patient-space, mm)
    /// for a Stage 1 field `v` (length `dims.0·dims.1·dims.2`).
    ///
    /// `transmembrane_amplitude_mv` is the rest-to-peak `V_m` the
    /// dimensionless field represents (a cited action-potential
    /// amplitude, caller-supplied); `conductivity_ratio` = `σ_i/σ_e`
    /// (dimensionless, caller-supplied — the screening knob of the
    /// unbounded-bath reduction). Both scales enter linearly, so timing
    /// and morphology are independent of them; only the amplitude is
    /// not.
    ///
    /// # Errors
    ///
    /// `Err` when `v`'s length is wrong, a scale is not finite and
    /// positive, or the electrode lies within one minimum spacing of a
    /// solid voxel centre — the unbounded-medium kernel is singular at
    /// the source, so electrodes at or inside the tissue are rejected
    /// rather than softened silently.
    pub fn potential_at(
        &self,
        electrode: Vec3,
        v: &[f64],
        transmembrane_amplitude_mv: f64,
        conductivity_ratio: f64,
    ) -> Result<f64> {
        if v.len() != self.solid.len() {
            return Err(EpError::InvalidLeadField {
                reason: format!(
                    "field length {} does not match the grid's {} voxels",
                    v.len(),
                    self.solid.len()
                ),
            });
        }
        if !transmembrane_amplitude_mv.is_finite() || transmembrane_amplitude_mv <= 0.0 {
            return Err(EpError::InvalidLeadField {
                reason: format!(
                    "transmembrane amplitude must be finite and positive, got \
                     {transmembrane_amplitude_mv}"
                ),
            });
        }
        if !conductivity_ratio.is_finite() || conductivity_ratio <= 0.0 {
            return Err(EpError::InvalidLeadField {
                reason: format!(
                    "conductivity ratio must be finite and positive, got {conductivity_ratio}"
                ),
            });
        }
        if !electrode.x.is_finite() || !electrode.y.is_finite() || !electrode.z.is_finite() {
            return Err(EpError::InvalidLeadField {
                reason: "electrode position must be finite".into(),
            });
        }

        let (nx, ny, nz) = self.dims;
        let min_spacing = self.spacing.0.min(self.spacing.1).min(self.spacing.2);
        let min_spacing2 = min_spacing * min_spacing;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let idx = self.index(x, y, z).expect("in-bounds by construction");
                    if self.solid[idx]
                        && (electrode - self.voxel_position(x as f64, y as f64, z as f64))
                            .norm_squared()
                            < min_spacing2
                    {
                        return Err(EpError::InvalidLeadField {
                            reason: "electrode lies within one spacing of a solid voxel; \
                                     the unbounded-medium kernel is singular there"
                                .into(),
                        });
                    }
                }
            }
        }

        let integral = self.source_integral(electrode, v);
        Ok(
            -conductivity_ratio * transmembrane_amplitude_mv / (4.0 * core::f64::consts::PI)
                * integral,
        )
    }

    /// The potentials at several electrodes (a lead set) for one field
    /// snapshot, in electrode order.
    ///
    /// # Errors
    ///
    /// As [`Self::potential_at`]; the first offending electrode rejects
    /// the whole set.
    pub fn potentials_at(
        &self,
        electrodes: &[Vec3],
        v: &[f64],
        transmembrane_amplitude_mv: f64,
        conductivity_ratio: f64,
    ) -> Result<Vec<f64>> {
        electrodes
            .iter()
            .map(|&e| self.potential_at(e, v, transmembrane_amplitude_mv, conductivity_ratio))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_solid(nx: usize, ny: usize, nz: usize, h: f64) -> SegmentationMask {
        SegmentationMask {
            dims: (nx, ny, nz),
            origin: Vec3::ZERO,
            row_dir: Vec3::X,
            col_dir: Vec3::Y,
            slice_dir: Vec3::Z,
            spacing: (h, h, h),
            voxels: vec![true; nx * ny * nz],
            hu: vec![0.0; nx * ny * nz],
        }
    }

    #[test]
    fn a_two_voxel_step_matches_the_hand_computed_face_sum() {
        // u = [1, -1, 0] on a 3x1x1 grid of unit spacing. The interior
        // solid-solid faces are exactly two (along x): 0|1 with jump -2
        // and 1|2 with jump +1; the 1-voxel width leaves no y/z faces.
        // K is 1/r from the electrode at (1, 0, 4).
        let proj = LeadFieldProjection::from_mask(&all_solid(3, 1, 1, 1.0)).expect("builds");
        let v = vec![1.0, -1.0, 0.0];
        let e = Vec3::new(1.0, 0.0, 4.0);
        let k = |x: f64| 1.0 / ((1.0 - x).powi(2) + 16.0).sqrt();
        let expected = -1.0 / (4.0 * core::f64::consts::PI)
            * ((-2.0) * (k(1.0) - k(0.0)) + 1.0 * (k(2.0) - k(1.0)));
        let got = proj.potential_at(e, &v, 1.0, 1.0).expect("well outside");
        assert!((got - expected).abs() < 1e-12, "{got} vs {expected}");
    }

    #[test]
    fn sign_and_scales_are_linear() {
        let proj = LeadFieldProjection::from_mask(&all_solid(3, 1, 1, 1.0)).expect("builds");
        let v = vec![1.0, -1.0, 0.0];
        let e = Vec3::new(1.0, 0.0, 4.0);
        let base = proj.potential_at(e, &v, 1.0, 1.0).expect("ok");
        let flipped = proj
            .potential_at(e, &[-1.0, 1.0, 0.0], 1.0, 1.0)
            .expect("ok");
        assert!(
            (flipped + base).abs() < 1e-12,
            "reversing the field flips the sign"
        );
        let doubled = proj.potential_at(e, &v, 2.0, 1.0).expect("ok");
        assert!(
            (doubled - 2.0 * base).abs() < 1e-12,
            "amplitude scales linearly"
        );
        let ratio = proj.potential_at(e, &v, 1.0, 3.0).expect("ok");
        assert!(
            (ratio - 3.0 * base).abs() < 1e-12,
            "conductivity ratio scales linearly"
        );
    }

    #[test]
    fn a_closed_front_cancels_and_an_open_front_carries_the_signal() {
        // The defining identity of the reduction: −∫∇u·∇K = −∮u∂K/∂n.
        // A closed front (steps up, steps back, wholly inside) has no
        // boundary term and integrates to ~zero; the same front opened
        // at the tissue boundary carries the classical solid-angle
        // signal, with the dipolar 1/R² far field of an open patch.
        let h = 0.5;
        let n = 31;
        let proj = LeadFieldProjection::from_mask(&all_solid(n, n, n, h)).expect("builds");
        let center = (n as f64 - 1.0) * h / 2.0;
        let w = 1.0;

        let mut closed = vec![0.0f64; n * n * n];
        let mut open = vec![0.0f64; n * n * n];
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let x0 = x as f64 * h - center;
                    let y0 = y as f64 * h - center;
                    let z0 = z as f64 * h - center;
                    let idx = (z * n + y) * n + x;
                    // A front crossing the grid's lower half: u steps up
                    // through a tanh band at z0 = −3 and (for the closed
                    // case) back down at z0 = +3 — 4.5 w from the z
                    // boundary — and is windowed in x/y so it is compact
                    // in every direction (u ~ 5e-3 at every grid face;
                    // an unwindowed band would stay open at the sides).
                    let up = ((z0 + 3.0) / w).tanh();
                    let down = ((3.0 - z0) / w).tanh();
                    let window = |t: f64| 0.5 * (1.0 - ((t.abs() - 5.0) / w).tanh());
                    closed[idx] = 0.25 * (1.0 + up) * (1.0 + down) * window(x0) * window(y0);
                    // The same front opened: it never steps back and is
                    // not windowed — the polarization runs to the grid's
                    // top and side boundaries.
                    open[idx] = 0.5 * (1.0 + up);
                }
            }
        }

        let axial = |dist: f64| Vec3::new(center, center, center + dist);
        let closed_near = proj
            .potential_at(axial(10.0), &closed, 1.0, 1.0)
            .expect("ok");
        let open_near = proj.potential_at(axial(10.0), &open, 1.0, 1.0).expect("ok");
        let open_far = proj.potential_at(axial(80.0), &open, 1.0, 1.0).expect("ok");
        assert!(open_near.abs() > 1e-3, "an open front must carry signal");
        assert!(
            closed_near.abs() < 0.05 * open_near.abs(),
            "a closed front must cancel: {closed_near} vs open {open_near}"
        );
        let ratio = open_far.abs() / open_near.abs();
        // Dipole far field: 8× the distance quarters the signal's
        // square-law fall (2³). At R = 80·h from a ~15·h front the
        // expansion is good to a few tens of percent; only a band is
        // honest.
        assert!(
            (ratio - 1.0 / 64.0).abs() < 0.05,
            "open-front far-field ratio {ratio}, dipole expects ≈ {}",
            1.0 / 64.0
        );
    }

    #[test]
    fn uniform_open_front_converges_to_the_solid_angle_integral() {
        // u = 1 below a plane, running to the grid boundary: the open
        // front. The integral collapses to (1/4π)·Ω, the solid angle the
        // grid's cross-section rectangle subtends at the electrode.
        // Independent reference: brute-force midpoint quadrature of
        // ∫∫ d/(x²+y²+d²)^{3/2} dA over the effective rectangle (voxel
        // centres ± h/2) — no memorised closed form.
        let omega = |nx: usize, ny: usize, h: f64, d: f64| {
            let half_x = nx as f64 * h / 2.0;
            let half_y = ny as f64 * h / 2.0;
            let m = 400;
            let mx = 2.0 * half_x / m as f64;
            let my = 2.0 * half_y / m as f64;
            let mut acc = 0.0;
            for i in 0..m {
                for j in 0..m {
                    let px = -half_x + (i as f64 + 0.5) * mx;
                    let py = -half_y + (j as f64 + 0.5) * my;
                    acc += d / (px * px + py * py + d * d).powf(1.5) * mx * my;
                }
            }
            acc
        };

        // u = 1 in layers 0..=6, 0 above; the interface plane sits at
        // z = 6.5·h, `d` below the electrode, on the axis.
        let run = |nx: usize, ny: usize, h: f64| {
            let nz = 12;
            let proj = LeadFieldProjection::from_mask(&all_solid(nx, ny, nz, h)).expect("builds");
            let mut v = vec![0.0f64; nx * ny * nz];
            for z in 0..7 {
                for y in 0..ny {
                    for x in 0..nx {
                        v[(z * ny + y) * nx + x] = 1.0;
                    }
                }
            }
            let d = 7.0;
            let electrode = Vec3::new(
                (nx as f64 - 1.0) * h / 2.0,
                (ny as f64 - 1.0) * h / 2.0,
                6.5 * h + d,
            );
            let got = proj.potential_at(electrode, &v, 1.0, 1.0).expect("ok");
            let reference = omega(nx, ny, h, d) / (4.0 * core::f64::consts::PI);
            (got, reference)
        };

        let (coarse, ref_coarse) = run(9, 9, 1.0);
        // nx·h = 9 in both runs: the same physical cross-section.
        let (fine, ref_fine) = run(18, 18, 0.5);
        // Same physical geometry at both resolutions: the references
        // agree, then the discrete error must shrink with refinement.
        assert!(
            (ref_coarse - ref_fine).abs() < 1e-4,
            "{ref_coarse} vs {ref_fine}"
        );
        let err_coarse = (coarse - ref_coarse).abs();
        let err_fine = (fine - ref_fine).abs();
        assert!(
            err_fine < err_coarse * 0.7,
            "refinement must converge: {err_coarse} -> {err_fine}"
        );
    }

    #[test]
    fn rejects_invalid_inputs() {
        let proj = LeadFieldProjection::from_mask(&all_solid(4, 4, 4, 1.0)).expect("builds");
        let v = vec![0.0; 64];
        // Wrong field length.
        assert!(proj
            .potential_at(Vec3::new(10.0, 0.0, 0.0), &v[..63], 1.0, 1.0)
            .is_err());
        // Non-finite or non-positive scales.
        assert!(proj
            .potential_at(Vec3::new(10.0, 0.0, 0.0), &v, 0.0, 1.0)
            .is_err());
        assert!(proj
            .potential_at(Vec3::new(10.0, 0.0, 0.0), &v, 1.0, -1.0)
            .is_err());
        assert!(proj
            .potential_at(Vec3::new(10.0, 0.0, 0.0), &v, f64::NAN, 1.0)
            .is_err());
        // Electrode at a voxel centre and within one spacing of one.
        assert!(proj
            .potential_at(Vec3::new(2.0, 1.0, 1.0), &v, 1.0, 1.0)
            .is_err());
        assert!(proj
            .potential_at(Vec3::new(2.0, 1.0, 1.9), &v, 1.0, 1.0)
            .is_err());
        // Non-finite electrode.
        assert!(proj
            .potential_at(Vec3::new(f64::NAN, 1.0, 1.0), &v, 1.0, 1.0)
            .is_err());
        // A well-placed electrode is accepted.
        assert!(proj
            .potential_at(Vec3::new(2.0, 1.0, 4.0), &v, 1.0, 1.0)
            .is_ok());
        // potentials_at shares the same contract.
        assert!(proj
            .potentials_at(
                &[Vec3::new(2.0, 1.0, 4.0), Vec3::new(2.0, 1.0, 1.0)],
                &v,
                1.0,
                1.0
            )
            .is_err());
    }

    #[test]
    fn an_empty_mask_is_rejected() {
        let mut mask = all_solid(2, 2, 2, 1.0);
        mask.voxels.fill(false);
        assert!(matches!(
            LeadFieldProjection::from_mask(&mask),
            Err(EpError::EmptyMask)
        ));
    }

    #[test]
    fn stage1_to_stage2_end_to_end() {
        // A real Stage 1 march: stimulate one end of a slab, let the
        // depolarization wave cross, then project the field at two
        // electrodes beyond the far face. The signal is nonzero, grows
        // toward the approaching wavefront, and decays with distance.
        use crate::{MitchellSchaefferParams, MonodomainTissue};
        let h = 0.5;
        let mask = all_solid(12, 12, 12, h);
        let params = MitchellSchaefferParams::human_ventricular_default();
        let mut tissue = MonodomainTissue::from_mask(&mask, 0.05, params).expect("builds");
        let dt = tissue.max_stable_dt() * 0.9;
        let stim_steps = (1.0 / dt).ceil() as u64;
        for _ in 0..stim_steps {
            tissue.stimulate(|x, _, _| x == 0, 3.0);
            tissue.step(dt).expect("stable step");
        }
        for _ in 0..(30.0 / dt).ceil() as u64 {
            tissue.step(dt).expect("stable step");
        }

        let proj = LeadFieldProjection::from_mask(&mask).expect("builds");
        let mut v = Vec::with_capacity(12 * 12 * 12);
        for z in 0..12 {
            for y in 0..12 {
                for x in 0..12 {
                    v.push(tissue.v_at(x, y, z).unwrap());
                }
            }
        }
        // The wave started at x = 0; after ~31 ms it has crossed much of
        // the 5.5 mm slab (CV ~ 0.4 mm/ms here), so a front is en route.
        let e_near = Vec3::new(12.0 * h + 3.0, 2.75, 2.75);
        let e_far = Vec3::new(12.0 * h + 6.0, 2.75, 2.75);
        let potentials = proj
            .potentials_at(&[e_near, e_far], &v, 100.0, 1.0)
            .expect("electrodes well outside");
        let (near, far) = (potentials[0], potentials[1]);
        assert!(near.is_finite() && far.is_finite());
        assert!(
            near.abs() > 1e-3,
            "a crossing wavefront must be visible: {near}"
        );
        assert!(
            far.abs() < near.abs(),
            "the signal must decay with distance: {far} vs {near}"
        );
    }
}
