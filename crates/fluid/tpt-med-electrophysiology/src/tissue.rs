//! Monodomain reaction-diffusion on a voxel grid
//! (`rfcs/0005-cardiac-electrophysiology.md` Stage 1).
//!
//! ```text
//! dV/dt = ∇·(D ∇V) + J_in(V, h) + J_out(V) + J_stim
//! ```
//!
//! Isotropic, homogeneous diffusivity `D` by default; caller-supplied
//! per-voxel anisotropy (fiber direction with longitudinal/transverse
//! conductivities) via [`MonodomainTissue::set_anisotropy`] — the crate
//! *accepts* a fiber field but still has no *source* for one (atlas or DTI
//! derivation stays external; see the RFC's Unresolved Questions). The
//! grid itself may be anisotropic in spacing (a real CT/MRI voxel is
//! rarely cubic); the discretisation and stability bound below account for
//! that directly rather than assuming `dx = dy = dz`.

use crate::error::{EpError, Result};
use crate::kinetics::MitchellSchaefferParams;
use tpt_med_geometry::Vec3;
use tpt_med_meshing::SegmentationMask;

/// Per-voxel anisotropic conductivity: a fiber direction with
/// longitudinal (along-fiber) and transverse conductivities about it —
/// the axisymmetric (transversely isotropic) reduction used when the
/// second fiber family is not resolved.
///
/// The face conductivity the 7-point stencil sees along a grid axis `a` is
/// the tensor projection `D_t + (D_l − D_t)·(d̂·â)²`, so a fiber aligned
/// with the grid conducts at `D_l` along itself and `D_t` across.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiberConductivity {
    /// Fiber direction (normalized internally; must be finite and
    /// non-zero).
    pub direction: Vec3,
    /// Longitudinal (along-fiber) diffusivity, `mm²/ms`.
    pub longitudinal: f64,
    /// Transverse diffusivity, `mm²/ms` (≤ longitudinal).
    pub transverse: f64,
}

impl FiberConductivity {
    fn projected(&self, axis: usize) -> f64 {
        let d = self.direction.normalize();
        let e = [d.x, d.y, d.z][axis];
        self.transverse + (self.longitudinal - self.transverse) * e * e
    }

    fn max_diffusivity(&self) -> f64 {
        self.longitudinal.max(self.transverse)
    }
}

/// A monodomain tissue on a voxel grid, built directly from a
/// `tpt-med-meshing` segmentation mask so the same CT/MRI-derived geometry
/// that feeds structural meshing feeds electrophysiology.
#[derive(Debug, Clone)]
pub struct MonodomainTissue {
    dims: (usize, usize, usize),
    spacing: (f64, f64, f64),
    solid: Vec<bool>,
    v: Vec<f64>,
    h: Vec<f64>,
    stim: Vec<f64>,
    diffusivity: f64,
    fibers: Option<Vec<FiberConductivity>>,
    params: MitchellSchaefferParams,
    time: f64,
    activation: Vec<Option<f64>>,
}

impl MonodomainTissue {
    /// Builds a resting tissue (`V = 0`, `h = 1` at every solid voxel) from
    /// `mask`. `diffusivity` is in `mm^2/ms`, matching the mask's spacing
    /// (mm) and the kinetics' time unit (ms). Errs if the mask has no solid
    /// voxels, or if `diffusivity` is not finite and positive.
    pub fn from_mask(
        mask: &SegmentationMask,
        diffusivity: f64,
        params: MitchellSchaefferParams,
    ) -> Result<Self> {
        if !diffusivity.is_finite() || diffusivity <= 0.0 {
            return Err(EpError::InvalidDiffusivity(diffusivity));
        }
        if !mask.voxels.iter().any(|&s| s) {
            return Err(EpError::EmptyMask);
        }
        let n = mask.voxels.len();
        Ok(Self {
            dims: mask.dims,
            spacing: mask.spacing,
            solid: mask.voxels.clone(),
            v: vec![0.0; n],
            h: vec![1.0; n],
            stim: vec![0.0; n],
            diffusivity,
            fibers: None,
            params,
            time: 0.0,
            activation: vec![None; n],
        })
    }

    /// Replaces the isotropic diffusivity with a caller-supplied
    /// **anisotropic conductivity field**: one [`FiberConductivity`] per
    /// solid voxel (the closure is evaluated in `(x, y, z)` order over the
    /// grid; non-solid voxels ignore their value). This is the acceptance
    /// point for an external fiber-field source (atlas or DTI derivation)
    /// — the crate does not source fibers itself.
    ///
    /// Validation per solid voxel: finite positive conductivities with
    /// `transverse ≤ longitudinal`, and a finite non-zero direction
    /// (normalized internally). Any violation rejects the whole field with
    /// [`EpError::InvalidAnisotropy`], leaving the tissue isotropic.
    pub fn set_anisotropy(
        &mut self,
        fiber_at: impl Fn(usize, usize, usize) -> Option<FiberConductivity>,
    ) -> Result<()> {
        let (nx, ny, nz) = self.dims;
        let mut fibers = vec![
            FiberConductivity {
                direction: Vec3::Z,
                longitudinal: self.diffusivity,
                transverse: self.diffusivity,
            };
            self.v.len()
        ];
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let idx = self.index(x, y, z).expect("in-bounds by construction");
                    if !self.solid[idx] {
                        continue;
                    }
                    let Some(f) = fiber_at(x, y, z) else {
                        return Err(EpError::InvalidAnisotropy {
                            reason: "missing fiber entry at a solid voxel",
                        });
                    };
                    if !f.longitudinal.is_finite()
                        || !f.transverse.is_finite()
                        || f.longitudinal <= 0.0
                        || f.transverse <= 0.0
                    {
                        return Err(EpError::InvalidAnisotropy {
                            reason: "conductivities must be finite and positive",
                        });
                    }
                    if f.transverse > f.longitudinal {
                        return Err(EpError::InvalidAnisotropy {
                            reason: "transverse diffusivity above longitudinal",
                        });
                    }
                    let d2 = f.direction.norm_squared();
                    if !d2.is_finite() || d2 <= 0.0 {
                        return Err(EpError::InvalidAnisotropy {
                            reason: "fiber direction must be finite and non-zero",
                        });
                    }
                    fibers[idx] = f;
                }
            }
        }
        self.fibers = Some(fibers);
        Ok(())
    }

    /// True when an anisotropic conductivity field is in effect.
    pub fn is_anisotropic(&self) -> bool {
        self.fibers.is_some()
    }

    /// The conductivity governing the explicit diffusion stability bound:
    /// the largest longitudinal/transverse value in the field when
    /// anisotropic, the isotropic diffusivity otherwise.
    fn max_diffusivity(&self) -> f64 {
        match &self.fibers {
            None => self.diffusivity,
            Some(fibers) => fibers
                .iter()
                .zip(&self.solid)
                .filter(|(_, &s)| s)
                .map(|(f, _)| f.max_diffusivity())
                .fold(0.0f64, f64::max),
        }
    }

    /// Voxel grid dimensions `(nx, ny, nz)`.
    pub fn dims(&self) -> (usize, usize, usize) {
        self.dims
    }

    /// Elapsed simulated time, ms.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// `V` at `(x, y, z)`, or `None` outside the grid.
    pub fn v_at(&self, x: usize, y: usize, z: usize) -> Option<f64> {
        self.index(x, y, z).map(|i| self.v[i])
    }

    /// `h` at `(x, y, z)`, or `None` outside the grid.
    pub fn h_at(&self, x: usize, y: usize, z: usize) -> Option<f64> {
        self.index(x, y, z).map(|i| self.h[i])
    }

    /// Sets the active stimulus for the next [`Self::step`] call: `region`
    /// is evaluated at every solid voxel, and matching voxels get
    /// `amplitude` (dimensionless current, added to `dV/dt`) for that one
    /// step. Replaces whatever stimulus was previously set — call this
    /// again every step a pulse should remain active (e.g. for a
    /// multi-step-long S1 stimulus), and stop calling it once the pulse
    /// ends.
    pub fn stimulate(&mut self, region: impl Fn(usize, usize, usize) -> bool, amplitude: f64) {
        self.stim.fill(0.0);
        let (nx, ny, nz) = self.dims;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let idx = self.index(x, y, z).expect("in-bounds by construction");
                    if self.solid[idx] && region(x, y, z) {
                        self.stim[idx] = amplitude;
                    }
                }
            }
        }
    }

    /// The largest `dt` (ms) this tissue's explicit RK2 update can take
    /// stably: the smaller of two independent bounds.
    ///
    /// - **Diffusion CFL**: `1 / (2 D_max (1/dx² + 1/dy² + 1/dz²))`, where
    ///   `D_max` is the largest conductivity in effect (isotropic
    ///   diffusivity, or the largest longitudinal/transverse value of an
    ///   anisotropic field) — the general (possibly anisotropic-spacing)
    ///   form of the isotropic `dx²/(6D)` bound for an explicit
    ///   7-point-stencil update.
    /// - **Reaction stiffness**: the upstroke term `J_in` acts on the
    ///   timescale `tau_in`, the fastest process in the kinetics (typically
    ///   far faster than any diffusion-limited step on a realistic voxel
    ///   grid) — an explicit step much larger than `tau_in` overshoots the
    ///   upstroke's stable manifold instead of saturating at it, diverging
    ///   rather than producing an action potential. Bounded here at
    ///   `0.2 * tau_in`, a standard conservative margin for explicit
    ///   integration of a stiff excitable ODE.
    pub fn max_stable_dt(&self) -> f64 {
        let (dx, dy, dz) = self.spacing;
        let diffusion_bound = 1.0
            / (2.0
                * self.max_diffusivity()
                * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
        let reaction_bound = 0.2 * self.params.tau_in;
        diffusion_bound.min(reaction_bound)
    }

    /// One explicit RK2 (Heun) step of size `dt` (ms). Errs
    /// (`EpError::UnstableTimeStep`) rather than stepping if `dt` is not
    /// finite/positive or exceeds [`Self::max_stable_dt`].
    pub fn step(&mut self, dt: f64) -> Result<()> {
        let max_stable_dt = self.max_stable_dt();
        if !dt.is_finite() || dt <= 0.0 || dt > max_stable_dt {
            return Err(EpError::UnstableTimeStep { dt, max_stable_dt });
        }

        let n = self.v.len();
        let (k1v, k1h) = self.rhs(&self.v, &self.h);

        let mut v1 = self.v.clone();
        let mut h1 = self.h.clone();
        for i in 0..n {
            if self.solid[i] {
                v1[i] += dt * k1v[i];
                h1[i] += dt * k1h[i];
            }
        }

        let (k2v, k2h) = self.rhs(&v1, &h1);
        for i in 0..n {
            if self.solid[i] {
                self.v[i] += dt * 0.5 * (k1v[i] + k2v[i]);
                self.h[i] += dt * 0.5 * (k1h[i] + k2h[i]);
            }
        }

        self.time += dt;
        let v_gate = self.params.v_gate;
        for i in 0..n {
            if self.solid[i] && self.activation[i].is_none() && self.v[i] >= v_gate {
                self.activation[i] = Some(self.time);
            }
        }

        self.stim.fill(0.0);
        Ok(())
    }

    /// First `V >= v_gate` crossing time (ms) per voxel, `None` if the
    /// voxel is not solid or has never activated.
    pub fn activation_map(&self) -> Vec<Option<f64>> {
        self.activation.clone()
    }

    fn index(&self, x: usize, y: usize, z: usize) -> Option<usize> {
        let (nx, ny, nz) = self.dims;
        if x < nx && y < ny && z < nz {
            Some((z * ny + y) * nx + x)
        } else {
            None
        }
    }

    fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        self.index(x, y, z).is_some_and(|i| self.solid[i])
    }

    /// No-flux (Neumann) discrete Laplacian at `(x, y, z)`: sums
    /// `D·(v_neighbour - v_self) / spacing^2` over every neighbour that both
    /// exists on the grid and is solid tissue, contributing zero flux
    /// across any grid edge or tissue/non-tissue boundary. `D` is the
    /// conductivity projected on the face normal — the isotropic
    /// diffusivity by default, or the cell's own fiber projection
    /// `D_t + (D_l − D_t)(d̂·â)²` when anisotropy is set.
    fn laplacian_at(&self, v: &[f64], x: usize, y: usize, z: usize) -> f64 {
        let (dx, dy, dz) = self.spacing;
        let idx = self.index(x, y, z).expect("in-bounds by construction");
        let vc = v[idx];
        let face_d = |axis: usize| match &self.fibers {
            None => self.diffusivity,
            Some(fibers) => fibers[idx].projected(axis),
        };
        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let dz2 = dz * dz;
        let mut lap = 0.0;

        if self.is_solid(x + 1, y, z) {
            lap += face_d(0) * (v[self.index(x + 1, y, z).unwrap()] - vc) / dx2;
        }
        if x > 0 && self.is_solid(x - 1, y, z) {
            lap += face_d(0) * (v[self.index(x - 1, y, z).unwrap()] - vc) / dx2;
        }
        if self.is_solid(x, y + 1, z) {
            lap += face_d(1) * (v[self.index(x, y + 1, z).unwrap()] - vc) / dy2;
        }
        if y > 0 && self.is_solid(x, y - 1, z) {
            lap += face_d(1) * (v[self.index(x, y - 1, z).unwrap()] - vc) / dy2;
        }
        if self.is_solid(x, y, z + 1) {
            lap += face_d(2) * (v[self.index(x, y, z + 1).unwrap()] - vc) / dz2;
        }
        if z > 0 && self.is_solid(x, y, z - 1) {
            lap += face_d(2) * (v[self.index(x, y, z - 1).unwrap()] - vc) / dz2;
        }
        lap
    }

    fn rhs(&self, v: &[f64], h: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = v.len();
        let mut dv = vec![0.0; n];
        let mut dh = vec![0.0; n];
        let (nx, ny, nz) = self.dims;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let idx = self.index(x, y, z).expect("in-bounds by construction");
                    if !self.solid[idx] {
                        continue;
                    }
                    let lap = self.laplacian_at(v, x, y, z);
                    dv[idx] =
                        self.diffusivity * lap + self.params.dv_dt(v[idx], h[idx], self.stim[idx]);
                    dh[idx] = self.params.dh_dt(v[idx], h[idx]);
                }
            }
        }
        (dv, dh)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_med_geometry::Vec3;

    fn uniform_mask(nx: usize, ny: usize, nz: usize, spacing: (f64, f64, f64)) -> SegmentationMask {
        SegmentationMask {
            dims: (nx, ny, nz),
            origin: Vec3::ZERO,
            row_dir: Vec3::new(1.0, 0.0, 0.0),
            col_dir: Vec3::new(0.0, 1.0, 0.0),
            slice_dir: Vec3::new(0.0, 0.0, 1.0),
            spacing,
            voxels: vec![true; nx * ny * nz],
            hu: vec![0.0; nx * ny * nz],
        }
    }

    #[test]
    fn rejects_empty_mask() {
        let mut mask = uniform_mask(4, 4, 1, (1.0, 1.0, 1.0));
        mask.voxels.fill(false);
        let params = MitchellSchaefferParams::human_ventricular_default();
        assert!(matches!(
            MonodomainTissue::from_mask(&mask, 0.1, params),
            Err(EpError::EmptyMask)
        ));
    }

    #[test]
    fn rejects_invalid_diffusivity() {
        let mask = uniform_mask(4, 4, 1, (1.0, 1.0, 1.0));
        let params = MitchellSchaefferParams::human_ventricular_default();
        assert!(MonodomainTissue::from_mask(&mask, 0.0, params).is_err());
        assert!(MonodomainTissue::from_mask(&mask, -1.0, params).is_err());
        assert!(MonodomainTissue::from_mask(&mask, f64::NAN, params).is_err());
    }

    #[test]
    fn rejects_a_dt_over_the_stability_bound() {
        let mask = uniform_mask(10, 10, 10, (0.5, 0.5, 0.5));
        let params = MitchellSchaefferParams::human_ventricular_default();
        let mut tissue = MonodomainTissue::from_mask(&mask, 0.1, params).expect("builds");
        let bound = tissue.max_stable_dt();
        assert!(tissue.step(bound * 1.5).is_err());
    }

    #[test]
    fn accepts_a_dt_under_the_stability_bound_and_remains_bounded() {
        let mask = uniform_mask(10, 10, 10, (0.5, 0.5, 0.5));
        let params = MitchellSchaefferParams::human_ventricular_default();
        let mut tissue = MonodomainTissue::from_mask(&mask, 0.1, params).expect("builds");
        let dt = tissue.max_stable_dt() * 0.5;
        tissue.stimulate(|x, _, _| x == 5, 2.0);
        for _ in 0..20 {
            tissue.step(dt).expect("stable step");
        }
        for x in 0..10 {
            for y in 0..10 {
                for z in 0..10 {
                    let v = tissue.v_at(x, y, z).unwrap();
                    assert!(v.is_finite() && (0.0..=1.5).contains(&v));
                }
            }
        }
    }

    #[test]
    fn stimulated_voxel_activates_and_diffuses_to_neighbours() {
        let mask = uniform_mask(20, 1, 1, (0.25, 1.0, 1.0));
        let params = MitchellSchaefferParams::human_ventricular_default();
        let mut tissue = MonodomainTissue::from_mask(&mask, 0.05, params).expect("builds");
        let dt = tissue.max_stable_dt() * 0.9;

        // Stimulate voxel 0 for a 1ms burst, then let it propagate for
        // 100ms — comfortably longer than the ~6ms conduction time this
        // grid/diffusivity implies (CV ~ sqrt(D/tau_in) ~ 0.4 mm/ms, 2.5mm
        // to voxel 10).
        let stim_steps = (1.0 / dt).ceil() as u64;
        let run_steps = (100.0 / dt).ceil() as u64;
        for _ in 0..stim_steps {
            tissue.stimulate(|x, _, _| x == 0, 3.0);
            tissue.step(dt).expect("stable step");
        }
        for _ in 0..run_steps {
            tissue.step(dt).expect("stable step");
        }

        let map = tissue.activation_map();
        assert!(map[0].is_some(), "stimulated voxel must activate");
        assert!(
            map[10].is_some(),
            "diffusion must eventually activate a downstream voxel"
        );
        // Activation time is monotonically non-decreasing with distance
        // from the stimulus site — a conduction wave, not simultaneous
        // activation.
        assert!(map[0].unwrap() < map[5].unwrap());
        assert!(map[5].unwrap() < map[10].unwrap());
    }

    #[test]
    fn fiber_projection_gives_longitudinal_along_and_transverse_across() {
        let along_x = FiberConductivity {
            direction: Vec3::new(2.0, 0.0, 0.0), // normalized internally
            longitudinal: 0.20,
            transverse: 0.02,
        };
        assert!(
            (along_x.projected(0) - 0.20).abs() < 1e-12,
            "along the fiber"
        );
        assert!(
            (along_x.projected(1) - 0.02).abs() < 1e-12,
            "across the fiber"
        );
        assert!((along_x.projected(2) - 0.02).abs() < 1e-12);
        // 45° in-plane fiber: the x and y faces each see the mean.
        let diagonal = FiberConductivity {
            direction: Vec3::new(1.0, 1.0, 0.0),
            longitudinal: 0.20,
            transverse: 0.02,
        };
        let expected = 0.5 * (0.20 + 0.02);
        assert!((diagonal.projected(0) - expected).abs() < 1e-12);
        assert!((diagonal.projected(1) - expected).abs() < 1e-12);
        assert!((diagonal.projected(2) - 0.02).abs() < 1e-12);
    }

    #[test]
    fn an_isotropic_fiber_field_reproduces_isotropic_conduction() {
        let mask = uniform_mask(16, 1, 1, (0.25, 1.0, 1.0));
        let params = MitchellSchaefferParams::human_ventricular_default();
        let d = 0.05;
        let mut plain = MonodomainTissue::from_mask(&mask, d, params).expect("builds");
        let mut fibered = MonodomainTissue::from_mask(&mask, d, params).expect("builds");
        fibered
            .set_anisotropy(|_, _, _| {
                Some(FiberConductivity {
                    direction: Vec3::X,
                    longitudinal: d,
                    transverse: d,
                })
            })
            .expect("valid field");
        assert!(fibered.is_anisotropic());

        for tissue in [&mut plain, &mut fibered] {
            let dt = tissue.max_stable_dt() * 0.9;
            for _ in 0..(1.0 / dt).ceil() as u64 {
                tissue.stimulate(|x, _, _| x == 0, 3.0);
                tissue.step(dt).expect("stable step");
            }
            for _ in 0..(80.0 / dt).ceil() as u64 {
                tissue.step(dt).expect("stable step");
            }
        }
        for x in 0..16 {
            let a = plain.activation_map()[x];
            let b = fibered.activation_map()[x];
            assert_eq!(a.is_some(), b.is_some(), "voxel {x}");
            if let (Some(ta), Some(tb)) = (a, b) {
                assert!((ta - tb).abs() < 1e-9, "voxel {x}: {ta} vs {tb}");
            }
        }
    }

    #[test]
    fn the_stability_bound_honours_the_fiber_maximum() {
        // Fine spacing so the diffusion CFL (not the 0.2·tau_in reaction
        // bound) is the governing limit for both fields.
        let mask = uniform_mask(6, 6, 1, (0.2, 0.2, 0.2));
        let params = MitchellSchaefferParams::human_ventricular_default();
        let mut tissue = MonodomainTissue::from_mask(&mask, 0.05, params).expect("builds");
        let isotropic_bound = tissue.max_stable_dt();
        tissue
            .set_anisotropy(|_, _, _| {
                Some(FiberConductivity {
                    direction: Vec3::X,
                    longitudinal: 0.20,
                    transverse: 0.02,
                })
            })
            .expect("valid field");
        let anisotropic_bound = tissue.max_stable_dt();
        assert!(
            anisotropic_bound < isotropic_bound,
            "a 0.20 longitudinal fiber must tighten the bound: {anisotropic_bound} vs {isotropic_bound}"
        );
        assert!(tissue.step(anisotropic_bound * 1.5).is_err());
        assert!(tissue.step(anisotropic_bound * 0.5).is_ok());
    }

    #[test]
    fn invalid_fiber_fields_are_rejected_wholesale() {
        let mask = uniform_mask(4, 4, 4, (1.0, 1.0, 1.0));
        let params = MitchellSchaefferParams::human_ventricular_default();
        type FiberField = Box<dyn Fn(usize, usize, usize) -> Option<FiberConductivity>>;
        let cases: Vec<(&str, FiberField)> = vec![
            (
                "transverse above longitudinal",
                Box::new(|_, _, _| {
                    Some(FiberConductivity {
                        direction: Vec3::X,
                        longitudinal: 0.02,
                        transverse: 0.20,
                    })
                }),
            ),
            (
                "zero direction",
                Box::new(|_, _, _| {
                    Some(FiberConductivity {
                        direction: Vec3::ZERO,
                        longitudinal: 0.20,
                        transverse: 0.02,
                    })
                }),
            ),
            ("missing entry", Box::new(|_, _, _| None)),
        ];
        for (why, fiber_at) in cases {
            let mut tissue = MonodomainTissue::from_mask(&mask, 0.1, params).expect("builds");
            let err = tissue.set_anisotropy(fiber_at).expect_err(why);
            assert!(
                matches!(err, EpError::InvalidAnisotropy { .. }),
                "{why}: {err:?}"
            );
            assert!(!tissue.is_anisotropic(), "{why}: field must stay unset");
        }
    }

    #[test]
    fn unstimulated_tissue_stays_at_rest() {
        let mask = uniform_mask(5, 5, 1, (1.0, 1.0, 1.0));
        let params = MitchellSchaefferParams::human_ventricular_default();
        let mut tissue = MonodomainTissue::from_mask(&mask, 0.1, params).expect("builds");
        let dt = tissue.max_stable_dt() * 0.5;
        for _ in 0..50 {
            tissue.step(dt).expect("stable step");
        }
        for x in 0..5 {
            for y in 0..5 {
                assert!(tissue.v_at(x, y, 0).unwrap().abs() < 1e-9);
            }
        }
        assert!(tissue.activation_map().iter().all(|a| a.is_none()));
    }
}
