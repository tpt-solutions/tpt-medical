//! Pure-Rust DICOM parsing and Hounsfield Unit (HU) mapping for CT-based
//! patient-specific modelling.
//!
//! Scope of the v0 parser (see `rfcs/0001-dicom-ingestion.md` for the full
//! roadmap):
//!
//! - Transfer syntaxes: **implicit VR little endian**
//!   (`1.2.840.10008.1.2`) and **explicit VR little endian**
//!   (`1.2.840.10008.1.2.1`) — the two uncompressed syntaxes carried by the
//!   overwhelming majority of archive exports.
//! - Single-frame CT/MR slices; sequences are parsed and skipped.
//! - Encapsulated/compressed pixel data (JPEG, JPEG-LS, JPEG 2000, RLE) is
//!   **not** decoded; [`DicomError::CompressedPixelData`] is returned.
//!
//! Pixel values are stored as `i32` raw stored values; HU are obtained via
//! `stored × RescaleSlope + RescaleIntercept` per slice.
//!
//! # Examples
//!
//! ```no_run
//! use std::path::Path;
//! use tpt_med_dicom::DicomSeries;
//!
//! let series = DicomSeries::load_from_dir(Path::new("test-data/dicom/synthetic_ct"))?;
//! println!("{} slices, spacing {:?}", series.slices.len(), series.pixel_spacing);
//! # Ok::<(), tpt_med_dicom::DicomError>(())
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bmd;
pub mod error;
pub mod hounsfield;
#[cfg(feature = "jpeg")]
pub mod jpeg;
#[cfg(feature = "jpeg2000")]
pub mod jpeg2000;
#[cfg(feature = "jpeg-ls")]
pub mod jpeg_ls;
pub mod parser;
pub mod phantom;
#[cfg(feature = "rle")]
pub mod rle;
pub mod series;
pub mod synthetic;
pub mod tags;

/// Crate result alias.
pub type Result<T> = core::result::Result<T, error::DicomError>;

pub use bmd::{AshFraction, BmdConvention, BmdToApparentDensity, BmdToAshDensity};
pub use error::DicomError;
pub use hounsfield::{BoneRegion, HounsfieldMapper, QctCalibration};
pub use parser::{DicomElement, DicomParser};
pub use phantom::{locate_phantom_centroid, sample_phantom_rods, PhantomModel, PhantomRod};
pub use series::{DicomSeries, DicomSlice};
pub use tags::{TransferSyntax, Vr};

pub(crate) const DICM_MAGIC: &[u8; 4] = b"DICM";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn synthetic_series_roundtrip() {
        let series = synthetic::SyntheticCtBuilder::new()
            .cols(8)
            .rows(8)
            .slices(3)
            .spacing(1.0, 1.0)
            .slice_thickness(1.0)
            .hu_fn(|_x, _y, _z| 0.0)
            .build("1.2.826.0.1.3680043.8.498.123.1");
        assert_eq!(series.slices.len(), 3);
        let re = DicomParser::parse_bytes(&series.slices[0].bytes).expect("reparse");
        assert_eq!(re.rows, 8);
        assert_eq!(re.columns, 8);
    }
}
