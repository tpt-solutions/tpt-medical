//! Error types for NIfTI-1 ingestion.

use std::path::PathBuf;

/// Errors produced while parsing a NIfTI-1 volume.
#[derive(Debug)]
pub enum NiftiError {
    /// File I/O failure.
    Io(std::io::Error),
    /// Not a NIfTI-1 file: neither little- nor big-endian interpretation of
    /// `sizeof_hdr` gives 348, or the magic bytes are not `"n+1"`.
    NotNifti(PathBuf),
    /// The input is gzip-compressed (`.nii.gz`, or a gzipped `.hdr`/`.img`
    /// part, detected by its `1F 8B` magic bytes) and this build has no
    /// decompressor: the `gzip` feature is off. Enable `gzip`, or
    /// decompress before ingestion; see `rfcs/0006-nifti-ingestion.md`.
    Gzipped(PathBuf),
    /// The stream announced itself as gzip (`1F 8B` magic bytes) but could
    /// not be decompressed: its gzip header, deflate stream, or checksum is
    /// corrupt. Only reachable with the `gzip` feature enabled — without it
    /// `Gzipped` is produced instead.
    CorruptGzip(std::io::Error),
    /// A `datatype` code this crate does not decode.
    UnsupportedDatatype(i16),
    /// The header's `bitpix` disagrees with the size implied by `datatype`.
    BitpixMismatch {
        /// The `datatype` field.
        datatype: i16,
        /// The `bitpix` field actually present.
        bitpix: i16,
    },
    /// The byte stream ended before a complete header or the expected voxel
    /// data was read.
    UnexpectedEof {
        /// Offset at which the stream ended.
        offset: usize,
        /// What was being read.
        while_reading: &'static str,
    },
    /// A header field was present but out of range (zero/negative
    /// dimension, `vox_offset` before the end of the header, …).
    BadValue {
        /// Why the value is rejected.
        reason: String,
    },
}

impl core::fmt::Display for NiftiError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            NiftiError::Io(e) => write!(f, "io error: {e}"),
            NiftiError::NotNifti(p) => write!(f, "not a NIfTI-1 file: {}", p.display()),
            NiftiError::Gzipped(p) => write!(
                f,
                "{} is gzip-compressed (.nii.gz); enable this crate's `gzip` \
                 feature or decompress before ingestion (see \
                 rfcs/0006-nifti-ingestion.md)",
                p.display()
            ),
            NiftiError::CorruptGzip(e) => write!(f, "corrupt gzip stream: {e}"),
            NiftiError::UnsupportedDatatype(dt) => {
                write!(f, "unsupported NIfTI datatype code {dt}")
            }
            NiftiError::BitpixMismatch { datatype, bitpix } => write!(
                f,
                "bitpix {bitpix} does not match datatype {datatype}'s expected size"
            ),
            NiftiError::UnexpectedEof {
                offset,
                while_reading,
            } => write!(
                f,
                "unexpected EOF at byte {offset} while reading {while_reading}"
            ),
            NiftiError::BadValue { reason } => write!(f, "bad header value: {reason}"),
        }
    }
}

impl std::error::Error for NiftiError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            NiftiError::Io(e) | NiftiError::CorruptGzip(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for NiftiError {
    fn from(e: std::io::Error) -> Self {
        NiftiError::Io(e)
    }
}

/// Crate result alias.
pub type Result<T> = core::result::Result<T, NiftiError>;
