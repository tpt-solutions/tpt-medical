//! Monodomain reaction-diffusion on a voxel grid
//! (`rfcs/0005-cardiac-electrophysiology.md` Stage 1).
//!
//! ```text
//! dV/dt = D * laplacian(V) + J_in(V, h) + J_out(V) + J_stim
//! ```
//!
//! Isotropic, homogeneous diffusivity `D` for v0 — anisotropic (fiber-
//! direction) conductivity needs a fiber-direction field this crate has no
//! source for (see the RFC's Unresolved Questions). The grid itself may be
//! anisotropic in spacing (a real CT/MRI voxel is rarely cubic); the
//! discretisation and stability bound below account for that directly
//! rather than assuming `dx = dy = dz`.

use crate::error::{EpError, Result};
use crate::kinetics::MitchellSchaefferParams;
use tpt_med_meshing::SegmentationMask;

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
            params,
            time: 0.0,
            activation: vec![None; n],
        })
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
    /// - **Diffusion CFL**: `1 / (2 D (1/dx² + 1/dy² + 1/dz²))`, the general
    ///   (possibly anisotropic-spacing) form of the isotropic `dx²/(6D)`
    ///   bound for an explicit 7-point-stencil update.
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
        let diffusion_bound =
            1.0 / (2.0 * self.diffusivity * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
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
    /// `(v_neighbour - v_self) / spacing^2` over every neighbour that both
    /// exists on the grid and is solid tissue, contributing zero flux
    /// across any grid edge or tissue/non-tissue boundary.
    fn laplacian_at(&self, v: &[f64], x: usize, y: usize, z: usize) -> f64 {
        let (dx, dy, dz) = self.spacing;
        let idx = self.index(x, y, z).expect("in-bounds by construction");
        let vc = v[idx];
        let mut lap = 0.0;

        if self.is_solid(x + 1, y, z) {
            lap += (v[self.index(x + 1, y, z).unwrap()] - vc) / (dx * dx);
        }
        if x > 0 && self.is_solid(x - 1, y, z) {
            lap += (v[self.index(x - 1, y, z).unwrap()] - vc) / (dx * dx);
        }
        if self.is_solid(x, y + 1, z) {
            lap += (v[self.index(x, y + 1, z).unwrap()] - vc) / (dy * dy);
        }
        if y > 0 && self.is_solid(x, y - 1, z) {
            lap += (v[self.index(x, y - 1, z).unwrap()] - vc) / (dy * dy);
        }
        if self.is_solid(x, y, z + 1) {
            lap += (v[self.index(x, y, z + 1).unwrap()] - vc) / (dz * dz);
        }
        if z > 0 && self.is_solid(x, y, z - 1) {
            lap += (v[self.index(x, y, z - 1).unwrap()] - vc) / (dz * dz);
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
