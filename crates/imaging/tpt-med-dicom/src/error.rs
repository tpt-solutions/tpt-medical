//! Error types for DICOM ingestion.

use std::path::PathBuf;

/// Errors produced while parsing DICOM files or assembling series.
#[derive(Debug)]
pub enum DicomError {
    /// File I/O failure.
    Io(std::io::Error),
    /// Not a DICOM Part-10 file (missing preamble magic).
    NotDicom(PathBuf),
    /// The byte stream ended before a complete element was read.
    UnexpectedEof {
        /// Offset at which the stream ended.
        offset: usize,
        /// What was being read.
        while_reading: &'static str,
    },
    /// A transfer syntax that is recognized but unsupported (compressed
    /// pixel data).
    CompressedPixelData(String),
    /// An unrecognised transfer syntax UID.
    UnknownTransferSyntax(String),
    /// A VR the parser cannot interpret in this position.
    UnsupportedVr {
        /// Tag the VR was found on.
        tag: (u16, u16),
        /// The two-character VR.
        vr: String,
    },
    /// A scalar element (US/IS/DS/…) could not be decoded.
    BadValue {
        /// Tag being decoded.
        tag: (u16, u16),
        /// Why decoding failed.
        reason: String,
    },
    /// A series is internally inconsistent (mixed geometry, empty dir, …).
    InconsistentSeries(String),
    /// A QCT phantom calibration could not be fitted from the supplied
    /// points (too few points, or a degenerate/collinear HU spread).
    Calibration(String),
}

impl core::fmt::Display for DicomError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DicomError::Io(e) => write!(f, "io error: {e}"),
            DicomError::NotDicom(p) => write!(f, "not a DICOM Part 10 file: {}", p.display()),
            DicomError::UnexpectedEof {
                offset,
                while_reading,
            } => {
                write!(
                    f,
                    "unexpected EOF at byte {offset} while reading {while_reading}"
                )
            }
            DicomError::CompressedPixelData(uid) => write!(
                f,
                "compressed transfer syntax {uid} not supported; decompress before ingestion \
                 (see rfcs/0001-dicom-ingestion.md)"
            ),
            DicomError::UnknownTransferSyntax(uid) => write!(f, "unknown transfer syntax {uid}"),
            DicomError::UnsupportedVr { tag, vr } => {
                write!(
                    f,
                    "unsupported VR {vr:?} on tag ({:04x},{:04x})",
                    tag.0, tag.1
                )
            }
            DicomError::BadValue { tag, reason } => {
                write!(
                    f,
                    "bad value on tag ({:04x},{:04x}): {reason}",
                    tag.0, tag.1
                )
            }
            DicomError::InconsistentSeries(why) => write!(f, "inconsistent series: {why}"),
            DicomError::Calibration(why) => write!(f, "QCT calibration: {why}"),
        }
    }
}

impl std::error::Error for DicomError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            DicomError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for DicomError {
    fn from(e: std::io::Error) -> Self {
        DicomError::Io(e)
    }
}

/// Crate result alias.
pub type Result<T> = core::result::Result<T, DicomError>;
