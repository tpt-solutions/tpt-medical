//! Voxel fluid domain with boundary classification.

/// Axis-aligned voxel fluid domain.
///
/// `mask[i][j][k]` (flattened) marks fluid cells. The domain drives the
/// MAC grid extents; boundary conditions (inlet velocity, outlet pressure
/// reference, walls elsewhere) are supplied to the solver.
#[derive(Debug, Clone)]
pub struct FluidDomain {
    /// Cell counts `(nx, ny, nz)`.
    pub dims: (usize, usize, usize),
    /// Uniform cell pitch (mm) — `(dx, dy, dz)`.
    pub spacing: (f64, f64, f64),
    /// Fluid flag per cell, index `(i*ny + j)*nz + k`.
    pub mask: Vec<bool>,
    /// Which axis the through-flow runs along (0=x, 1=y, 2=z).
    pub flow_axis: usize,
    /// True if flow enters at the low face of `flow_axis` (vs the high
    /// face).
    pub inlet_low: bool,
}

impl FluidDomain {
    /// Builds a domain from a mask.
    pub fn from_mask(
        dims: (usize, usize, usize),
        spacing: (f64, f64, f64),
        mask: Vec<bool>,
        flow_axis: usize,
        inlet_low: bool,
    ) -> Self {
        assert_eq!(mask.len(), dims.0 * dims.1 * dims.2, "mask length mismatch");
        Self {
            dims,
            spacing,
            mask,
            flow_axis,
            inlet_low,
        }
    }

    /// Cylinder of radius `r` (cells) along the flow axis — the canonical
    /// verification domain.
    pub fn cylinder(
        n_axial: usize,
        n_radius: usize,
        radius_cells: f64,
        spacing: f64,
        flow_axis: usize,
    ) -> Self {
        let dims = match flow_axis {
            0 => (n_axial, n_radius, n_radius),
            2 => (n_radius, n_radius, n_axial),
            _ => (n_radius, n_axial, n_radius),
        };
        let (nx, ny, nz) = dims;
        // Radial-plane centers must come from the RADIAL dimensions.
        let (n1, n2) = match flow_axis {
            0 => (ny, nz),
            1 => (nx, nz),
            _ => (nx, ny),
        };
        let c1 = (n1 as f64 - 1.0) / 2.0;
        let c2 = (n2 as f64 - 1.0) / 2.0;
        let mut mask = vec![false; nx * ny * nz];
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let (a, b) = match flow_axis {
                        0 => (j as f64 - c1, k as f64 - c2),
                        1 => (i as f64 - c1, k as f64 - c2),
                        _ => (i as f64 - c1, j as f64 - c2),
                    };
                    let radial = (a * a + b * b).sqrt();
                    mask[(i * ny + j) * nz + k] = radial <= radius_cells;
                }
            }
        }
        Self::from_mask(dims, (spacing, spacing, spacing), mask, flow_axis, true)
    }

    /// Flat cell index.
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        let (_nx, ny, nz) = self.dims;
        (i * ny + j) * nz + k
    }

    /// Fluid flag (out-of-bounds → wall).
    pub fn is_fluid(&self, i: i64, j: i64, k: i64) -> bool {
        let (nx, ny, nz) = self.dims;
        if i < 0 || j < 0 || k < 0 || i >= nx as i64 || j >= ny as i64 || k >= nz as i64 {
            return false;
        }
        self.mask[self.index(i as usize, j as usize, k as usize)]
    }

    /// Cell count and volume of the fluid region (mm³).
    pub fn fluid_stats(&self) -> (usize, f64) {
        let n = self.mask.iter().filter(|&&f| f).count();
        let cell_vol = self.spacing.0 * self.spacing.1 * self.spacing.2;
        (n, n as f64 * cell_vol)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cylinder_mask_counts() {
        let d = FluidDomain::cylinder(16, 12, 4.5, 0.5, 0);
        let (n, vol) = d.fluid_stats();
        assert!(n > 0);
        // Radius 4.5 cells → area ≈ π·4.5² ≈ 63.6 cells; ×16 slices.
        assert!(n > 16 * 40 && n < 16 * 72, "cells {n}");
        let cell = 0.5f64.powi(3);
        assert!((vol - n as f64 * cell).abs() < 1e-9);
    }

    #[test]
    fn bounds_checked_fluid_flag() {
        let d = FluidDomain::cylinder(8, 8, 3.0, 1.0, 2);
        assert!(!d.is_fluid(-1, 0, 0));
        assert!(!d.is_fluid(0, 0, 0) || d.is_fluid(0, 0, 0)); // no panic
        assert!(d.is_fluid(4, 4, 4), "center of axial cylinder");
    }
}
