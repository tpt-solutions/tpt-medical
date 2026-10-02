//! Reduced fluid–structure coupling — the wall-compliance slice of the
//! hemodynamics roadmap.
//!
//! The vessel wall is an **axisymmetric membrane**: each axial section's
//! radius relaxes first-order toward a pressure-set equilibrium
//!
//! ```text
//! r_eq = r₀·(1 + C·(p − p_ext)),   dr/dt = (r_eq − r)/τ
//! ```
//!
//! with `C` the caller-cited fractional compliance (radius change per
//! MPa) and `τ` the wall's relaxation time — and the flow domain's mask
//! is rebuilt from the radii every step, so the fluid feels the wall
//! moving and the wall feels the fluid's pressure. This is the reduced
//! "compliant tube" coupling (pressure → area via the tube law, area →
//! mask), not an immersed-boundary or ALE solve: the grid is fixed and
//! the stair-stepped mask moves through it, newly-wetted cells start at
//! rest, and the wall never moves *inflow* of the order of a cell per
//! step — `τ` should keep the per-step radius change well inside a
//! cell. Full FSI stays with the substrate upgrade path.
//!
//! Pressures are read from the projection potential in the same units
//! the Windkessel anchor writes: gauge-relative differences convert with
//! `p[MPa] = ΔΠ·ρ/1e6`, so only transmural-relevant *differences*
//! (section vs outlet) drive the wall — the absolute gauge never does.
use crate::domain::FluidDomain;
use crate::solver::HemodynamicsSolver;

/// The compliant membrane wall law (see the module docs).
#[derive(Debug, Clone)]
pub struct MembraneWall {
    /// Base (unloaded) wall radius (mm).
    pub base_radius_mm: f64,
    /// Fractional radius change per MPa of transmural pressure —
    /// caller-cited (a compliance of 0.05/MPa means 5 % radius growth
    /// per MPa above `external_pressure_mpa`).
    pub compliance_per_mpa: f64,
    /// Wall relaxation time τ (s); each coupled step moves the radius a
    /// fraction `dt/τ` toward its equilibrium.
    pub relaxation_time_s: f64,
    /// External (transmural reference) pressure (MPa), typically ~0.01
    /// for tissue-scale screening — caller-supplied.
    pub external_pressure_mpa: f64,
}

impl MembraneWall {
    /// Validates; `Err` for a non-positive radius, relaxation time, or a
    /// negative compliance.
    pub fn new(
        base_radius_mm: f64,
        compliance_per_mpa: f64,
        relaxation_time_s: f64,
        external_pressure_mpa: f64,
    ) -> Result<Self, String> {
        if !base_radius_mm.is_finite() || base_radius_mm <= 0.0 {
            return Err(format!(
                "base radius must be positive and finite, got {base_radius_mm}"
            ));
        }
        if !compliance_per_mpa.is_finite() || compliance_per_mpa < 0.0 {
            return Err(format!(
                "compliance must be non-negative and finite, got {compliance_per_mpa}"
            ));
        }
        if !relaxation_time_s.is_finite() || relaxation_time_s <= 0.0 {
            return Err(format!(
                "relaxation time must be positive and finite, got {relaxation_time_s}"
            ));
        }
        if !external_pressure_mpa.is_finite() {
            return Err("external pressure must be finite".into());
        }
        Ok(Self {
            base_radius_mm,
            compliance_per_mpa,
            relaxation_time_s,
            external_pressure_mpa,
        })
    }

    /// The pressure-set equilibrium radius (mm) at gauge-transmural
    /// pressure `p_rel_mpa`.
    pub fn equilibrium_radius(&self, p_rel_mpa: f64) -> f64 {
        self.base_radius_mm
            * (1.0 + self.compliance_per_mpa * (p_rel_mpa - self.external_pressure_mpa))
    }
}

/// The membrane wall's per-section state: one radius per axial section
/// (mm), persisted across the coupled march.
#[derive(Debug, Clone)]
pub struct MembraneWallState {
    /// Wall radius (mm) per axial section.
    pub radii_mm: Vec<f64>,
}

/// The outcome of one coupled step: the gauge-transmural pressure (MPa)
/// each section's wall felt, index-aligned with the state's radii.
#[derive(Debug, Clone)]
pub struct CoupledWallStep {
    /// Section-mean gauge-transmural pressure (MPa) driving the wall.
    pub section_pressure_mpa: Vec<f64>,
}

/// One coupled fluid–wall step: advance the flow, read each section's
/// gauge-transmural pressure, relax the membrane radii toward their
/// pressure-set equilibrium, and rebuild the domain mask from the
/// radii. The grid dimensions never change (the mask moves through the
/// fixed grid); the solver's arrays are untouched in size, so newly
/// wetted cells simply start at rest.
///
/// A zero-compliance wall is an exact no-op on the geometry: the coupled
/// march reproduces the fixed-domain run step for step (asserted in the
/// tests).
pub fn step_coupled(
    solver: &mut HemodynamicsSolver,
    wall: &MembraneWall,
    state: &mut MembraneWallState,
) -> CoupledWallStep {
    let (nx, ny, nz) = solver.domain.dims;
    assert_eq!(state.radii_mm.len(), nx, "one radius per axial section");
    let (_dx, dy, dz) = solver.domain.spacing;

    // 1. Flow.
    solver.step();

    // 2. Section-mean gauge pressure: the accumulated projection
    //    potential is gauge-relative to the outlet anchor, so a
    //    section's mean relative to the outlet layer is the transmural
    //    driver. Internal units convert as p[MPa] = ΔΠ·ρ/1e6 — the exact
    //    inverse of the Windkessel anchor's write convention.
    let rho = solver.config.density;
    let mut section_pressure = vec![0.0f64; nx];
    for i in 0..nx {
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for j in 0..ny {
            for k in 0..nz {
                let idx = solver.domain.index(i, j, k);
                if solver.domain.mask[idx] {
                    sum += solver.p[idx];
                    count += 1;
                }
            }
        }
        section_pressure[i] = if count > 0 {
            sum / count as f64
        } else {
            f64::NAN
        };
    }
    // Gauge reference: the last (outlet) section; fall back to the first
    // fluid section behind it if the outlet is occluded.
    let outlet_mean = section_pressure[nx - 1];
    let gauge = if outlet_mean.is_finite() {
        outlet_mean
    } else {
        section_pressure
            .iter()
            .rev()
            .find(|p| p.is_finite())
            .copied()
            .unwrap_or(0.0)
    };
    for p in section_pressure.iter_mut() {
        *p = if p.is_finite() {
            (*p - gauge) * rho / 1.0e6
        } else {
            // No fluid in this section: the wall sees the last known
            // driving pressure (keep the previous equilibrium).
            0.0
        };
    }

    // 3. Relax the radii toward the pressure-set equilibrium, clamped to
    //    the grid's radial half-extent.
    let relax = (solver.config.dt / wall.relaxation_time_s).min(1.0);
    let grid_half = 0.5 * dy.min(dz) * (ny.min(nz) as f64 - 1.0);
    for (r, &p) in state.radii_mm.iter_mut().zip(&section_pressure) {
        let r_eq = wall
            .equilibrium_radius(p)
            .clamp(0.5 * dy.min(dz), grid_half);
        *r += (r_eq - *r) * relax;
        *r = r.clamp(0.5 * dy.min(dz), grid_half);
    }

    // 4. Rebuild the mask from the radii (the cylinder convention:
    //    radial centre at ((n−1)/2, (n−1)/2) of the cross-section),
    //    measuring radial distance in mm so anisotropic radial spacing
    //    stays correct.
    let (c1, c2) = ((ny as f64 - 1.0) / 2.0, (nz as f64 - 1.0) / 2.0);
    let mut mask = vec![false; solver.domain.mask.len()];
    for i in 0..nx {
        let r = state.radii_mm[i];
        for j in 0..ny {
            for k in 0..nz {
                let a = (j as f64 - c1) * dy;
                let b = (k as f64 - c2) * dz;
                mask[solver.domain.index(i, j, k)] = (a * a + b * b).sqrt() <= r;
            }
        }
    }
    solver.domain.mask = mask;

    CoupledWallStep {
        section_pressure_mpa: section_pressure,
    }
}

/// The membrane state for a [`FluidDomain::cylinder`]-shaped start: one
/// base radius per axial section.
pub fn initial_state(domain: &FluidDomain, wall: &MembraneWall) -> MembraneWallState {
    MembraneWallState {
        radii_mm: vec![wall.base_radius_mm; domain.dims.0],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blood::BloodModel;
    use crate::solver::SolverConfig;

    fn fsi_solver() -> (HemodynamicsSolver, MembraneWall, MembraneWallState) {
        let domain = FluidDomain::cylinder(16, 9, 2.0, 0.5, 0);
        let wall = MembraneWall::new(1.0, 0.05, 1.0e-3, 0.01).expect("valid wall");
        let state = initial_state(&domain, &wall);
        let solver = HemodynamicsSolver::new(
            domain,
            BloodModel::Newtonian { viscosity: 0.0035 },
            50.0,
            SolverConfig::default(),
        );
        (solver, wall, state)
    }

    #[test]
    fn a_rigid_wall_reproduces_the_fixed_domain_run_exactly() {
        let (mut coupled, wall, mut state) = fsi_solver();
        // Zero compliance: the equilibrium radius never leaves the base.
        let rigid = MembraneWall {
            compliance_per_mpa: 0.0,
            ..wall
        };
        let (mut plain, ..) = {
            let domain = FluidDomain::cylinder(16, 9, 2.0, 0.5, 0);
            let solver = HemodynamicsSolver::new(
                domain,
                BloodModel::Newtonian { viscosity: 0.0035 },
                50.0,
                SolverConfig::default(),
            );
            (solver, 0, 0)
        };
        for _ in 0..12 {
            step_coupled(&mut coupled, &rigid, &mut state);
            plain.step();
        }
        // Bitwise-identical flow fields: the driver's geometry path is a
        // no-op at zero compliance.
        assert_eq!(coupled.u, plain.u);
        assert_eq!(coupled.v, plain.v);
        assert_eq!(coupled.w, plain.w);
        assert_eq!(coupled.p, plain.p);
        assert_eq!(coupled.domain.mask, plain.domain.mask);
    }

    #[test]
    fn compliant_wall_settles_on_the_pressure_set_equilibrium() {
        let (mut solver, wall, mut state) = fsi_solver();
        for _ in 0..400 {
            step_coupled(&mut solver, &wall, &mut state);
        }
        // Recompute the sections' gauge pressures from the final state
        // and check every radius against the law it claims to follow.
        let (nx, ny, nz) = solver.domain.dims;
        let rho = solver.config.density;
        let mut means = vec![0.0f64; nx];
        for i in 0..nx {
            let mut sum = 0.0;
            let mut count = 0;
            for j in 0..ny {
                for k in 0..nz {
                    let idx = solver.domain.index(i, j, k);
                    if solver.domain.mask[idx] {
                        sum += solver.p[idx];
                        count += 1;
                    }
                }
            }
            means[i] = if count > 0 {
                sum / count as f64
            } else {
                means[i - 1]
            };
        }
        let gauge = means[nx - 1];
        for i in 0..nx {
            let p_rel = (means[i] - gauge) * rho / 1.0e6;
            let expected = wall.equilibrium_radius(p_rel);
            assert!(
                (state.radii_mm[i] - expected).abs() < 1e-6,
                "section {i}: radius {} vs equilibrium {expected}",
                state.radii_mm[i]
            );
        }
        // The wall actually responded: at least one section left its
        // base radius (the flow's pressure field is not identically the
        // external reference).
        assert!(state
            .radii_mm
            .iter()
            .any(|&r| (r - wall.base_radius_mm).abs() > 1e-4));
    }

    #[test]
    fn the_mask_tracks_the_radii() {
        let (mut solver, wall, mut state) = fsi_solver();
        for _ in 0..60 {
            step_coupled(&mut solver, &wall, &mut state);
        }
        let (nx, ny, nz) = solver.domain.dims;
        let (_dx, dy, dz) = solver.domain.spacing;
        let (c1, c2) = ((ny as f64 - 1.0) / 2.0, (nz as f64 - 1.0) / 2.0);
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let a = (j as f64 - c1) * dy;
                    let b = (k as f64 - c2) * dz;
                    let expected = (a * a + b * b).sqrt() <= state.radii_mm[i];
                    assert_eq!(solver.domain.mask[solver.domain.index(i, j, k)], expected);
                }
            }
        }
    }

    #[test]
    fn wall_parameters_are_validated() {
        assert!(MembraneWall::new(0.0, 0.05, 1.0e-3, 0.01).is_err());
        assert!(MembraneWall::new(1.0, -0.05, 1.0e-3, 0.01).is_err());
        assert!(MembraneWall::new(1.0, 0.05, 0.0, 0.01).is_err());
        assert!(MembraneWall::new(1.0, 0.05, 1.0e-3, f64::NAN).is_err());
        assert!(MembraneWall::new(1.0, 0.05, 1.0e-3, 0.01).is_ok());
    }
}
