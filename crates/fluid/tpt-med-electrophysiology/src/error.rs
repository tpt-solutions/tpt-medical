//! Error types for monodomain electrophysiology.

/// Errors produced while building or stepping a [`crate::MonodomainTissue`],
/// or running an [`crate::S1S2Protocol`].
#[derive(Debug, Clone, PartialEq)]
pub enum EpError {
    /// A [`crate::tissue::MonodomainTissue::from_mask`] was given a mask
    /// with no solid (tissue) voxels — nothing to build a tissue on.
    EmptyMask,
    /// A [`crate::MitchellSchaefferParams`] field was NaN or infinite.
    NonFiniteParameter {
        /// The field's name.
        name: &'static str,
    },
    /// `step(dt)` was called with a `dt` that violates the explicit
    /// diffusion stability bound for this tissue's grid spacing and
    /// diffusivity.
    UnstableTimeStep {
        /// The requested step.
        dt: f64,
        /// The largest `dt` this tissue's grid/diffusivity can take
        /// explicitly.
        max_stable_dt: f64,
    },
    /// A diffusivity was not finite or not positive.
    InvalidDiffusivity(f64),
    /// A fiber-conductivity field was rejected (bad direction, transverse
    /// diffusivity above the longitudinal one, or a non-finite value).
    InvalidAnisotropy {
        /// What was wrong.
        reason: &'static str,
    },
    /// An [`crate::S1S2Protocol`] failed to establish even the baseline S1
    /// pacing train (no measurable action potential was detected at all —
    /// usually a cycle length shorter than the model's own refractory
    /// period, or a stimulus too weak to reach `v_gate`).
    RestitutionFailed(String),
    /// A [`crate::leadfield::LeadFieldProjection`] input was invalid: a
    /// potential field of the wrong length, a non-finite scale, or an
    /// electrode at or inside the tissue, where the unbounded-medium
    /// kernel this projection uses is singular.
    InvalidLeadField {
        /// What was wrong.
        reason: String,
    },
}

impl core::fmt::Display for EpError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            EpError::EmptyMask => write!(f, "segmentation mask has no solid (tissue) voxels"),
            EpError::NonFiniteParameter { name } => {
                write!(f, "parameter {name} must be finite")
            }
            EpError::UnstableTimeStep { dt, max_stable_dt } => write!(
                f,
                "dt={dt} exceeds the explicit diffusion stability bound \
                 (max stable dt={max_stable_dt}) for this grid/diffusivity"
            ),
            EpError::InvalidDiffusivity(d) => {
                write!(f, "diffusivity must be finite and positive, got {d}")
            }
            EpError::InvalidAnisotropy { reason } => {
                write!(f, "fiber conductivity field rejected: {reason}")
            }
            EpError::RestitutionFailed(why) => write!(f, "S1-S2 restitution protocol: {why}"),
            EpError::InvalidLeadField { reason } => {
                write!(f, "lead-field projection rejected: {reason}")
            }
        }
    }
}

impl std::error::Error for EpError {}

/// Crate result alias.
pub type Result<T> = core::result::Result<T, EpError>;
