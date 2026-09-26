//! Monodomain cardiac electrophysiology on voxel geometry
//! (`rfcs/0005-cardiac-electrophysiology.md` Stage 1).
//!
//! Ionic kinetics: **Mitchell–Schaeffer** (Mitchell & Schaeffer, 2003), a
//! two-variable model chosen over ten-plus-variable ionic detail
//! (Luo–Rudy, ten Tusscher) for screening fidelity with tractable parameter
//! hygiene. Tissue: isotropic, homogeneous-diffusivity monodomain
//! reaction-diffusion on the same voxel grid `tpt-med-meshing` already
//! builds from CT/MRI (`SegmentationMask`), stepped by explicit RK2.
//!
//! **Question of interest:** given a tissue mask, a conductivity, and a
//! pacing protocol, what is the local activation time at each point in the
//! tissue? **Model risk:** low — no output here feeds a device-sizing or
//! surgical-planning decision. **Model influence:** "supporting" (ASME
//! V&V 40): illustrative conduction behaviour with population-level,
//! literature-cited parameters, no patient-specific calibration. See the
//! RFC for the full V&V 40 discussion and what Stages 2/3 would need before
//! that influence level could honestly rise.
//!
//! # Examples
//!
//! ```
//! use tpt_med_electrophysiology::{MitchellSchaefferParams, MonodomainTissue};
//! use tpt_med_geometry::Vec3;
//! use tpt_med_meshing::SegmentationMask;
//!
//! let mask = SegmentationMask {
//!     dims: (10, 10, 1),
//!     origin: Vec3::ZERO,
//!     row_dir: Vec3::new(1.0, 0.0, 0.0),
//!     col_dir: Vec3::new(0.0, 1.0, 0.0),
//!     slice_dir: Vec3::new(0.0, 0.0, 1.0),
//!     spacing: (0.5, 0.5, 0.5),
//!     voxels: vec![true; 100],
//!     hu: vec![0.0; 100],
//! };
//!
//! let params = MitchellSchaefferParams::human_ventricular_default();
//! let mut tissue = MonodomainTissue::from_mask(&mask, 0.1, params)?;
//! let dt = tissue.max_stable_dt() * 0.9;
//!
//! tissue.stimulate(|x, _, _| x == 0, 2.0);
//! tissue.step(dt)?;
//! # Ok::<(), tpt_med_electrophysiology::EpError>(())
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod kinetics;
pub mod restitution;
pub mod tissue;

pub use error::EpError;
pub use kinetics::MitchellSchaefferParams;
pub use restitution::S1S2Protocol;
pub use tissue::MonodomainTissue;

/// Crate result alias.
pub type Result<T> = core::result::Result<T, EpError>;
