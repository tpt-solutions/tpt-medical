//! Incompressible Navier–Stokes CFD for blood flow on voxel domains.
//!
//! A staggered-grid (MAC) projection method: explicit advection/viscous
//! sub-step, then a pressure-Poisson projection enforcing incompressibility
//! inside a voxel fluid mask (stair-step walls). Blood rheology supports
//! Newtonian, Carreau–Yasuda, and Casson models. Post-processing computes
//! wall shear stress (WSS) and the oscillatory shear index (OSI).
//!
//! This is a screening-grade solver: no body-fitted meshes, no turbulence
//! modelling, laminar range only (the relevant regime for most arterial
//! screening at Re < 2000). Unstructured/high-fidelity CFD is the
//! documented upgrade path via `tpt-sci-cfd-core` / `tpt-sci-hemodynamics`
//! (pinned in the workspace manifest).
//!
//! Units: lengths mm, velocity mm/s, time s, pressure/viscosity in SI-derived
//! consistent set (μ in Pa·s, stress in Pa).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod blood;
pub mod domain;
pub mod fsi;
pub mod heat;
pub mod solver;
pub mod wss;

pub use blood::BloodModel;
pub use domain::FluidDomain;
pub use fsi::{initial_state, step_coupled, CoupledWallStep, MembraneWall, MembraneWallState};
pub use heat::{bulk, conjugate_step, stable_time_step, step, ScalarWall};
pub use solver::{HemodynamicsSolver, PressureSolver, SolverConfig, SteadyStats};
pub use wss::{
    extract_wss, wss_resolution_study, OsiAccumulator, WssField, WssResolutionLevel,
    WssResolutionStudy,
};

#[cfg(test)]
mod tests;
