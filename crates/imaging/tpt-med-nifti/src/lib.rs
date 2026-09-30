//! Pure-Rust NIfTI volume parsing (versions 1 and 2) for research-space
//! imaging.
//!
//! Reads a NIfTI volume into [`NiftiVolume`]: voxel dimensions, physical
//! spacing, RAS orientation, and per-voxel scaled values, with the parsed
//! [`NiftiVersion`] recorded on the volume. Both format versions are
//! accepted transparently at the same entry points — the 348-byte
//! NIfTI-1 header (`n+1`/`ni1`) and the 540-byte NIfTI-2 header
//! (`n+2`/`ni2`, with its `i64` dims and `f64` geometry) normalize into
//! one shared decode path. Three layouts are accepted — the single-file
//! `.nii`, the dual-file `.hdr`/`.img` pair (via
//! [`NiftiVolume::parse_dual_file`] / [`NiftiVolume::parse_dual_bytes`]),
//! and, behind the off-by-default `gzip` feature, a gzip-compressed
//! `.nii.gz` (or gzipped `.hdr`/`.img` parts). It does not decode
//! acquisition metadata (no `descrip`/`aux_file`/`intent_*`/slice-timing
//! fields) — see `rfcs/0006-nifti-ingestion.md` for the full scope.
//!
//! # Why
//!
//! `tpt-med-dicom` covers DICOM; this crate exists because research tooling
//! (`dcm2niix`, FSL, FreeSurfer, ANTs, most public imaging datasets) works in
//! NIfTI, not DICOM, once a scan leaves the archive. Converting back to
//! DICOM just to enter this pipeline is the gap this crate closes.
//!
//! # Verification
//!
//! The qform quaternion → rotation formula is checked against a hand-derived
//! closed-form case and an orthonormality property test (see `geometry.rs`).
//! The header/voxel-data decode itself is checked by round-tripping a
//! `SyntheticNiftiBuilder` through the parser, plus one hand-built minimal
//! header verified field-by-field against values computed by hand from the
//! bytes — not against this crate's own writer. The dual-file path is
//! round-tripped against the single-file one on identical geometry and data,
//! and (with `gzip` on) a compressed volume is asserted to parse to the same
//! result as the plain bytes it wraps. There is no golden dataset entry (no
//! simulation output originates here) and no test against a real-world
//! writer's output (`dcm2niix`, FSL, …) — see the RFC's "What remains
//! unverified".

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
mod geometry;
mod header;
pub mod synthetic;

use std::borrow::Cow;
use std::path::{Path, PathBuf};

pub use error::{NiftiError, Result};
pub use header::Datatype;
use header::{
    parse_header, parse_header2, parse_header2_dual, parse_header_dual, HEADER2_LEN, HEADER_LEN,
};
use tpt_med_geometry::{Mat3, Vec3};

/// The NIfTI format version of a parsed volume. The two versions carry
/// identical semantics for everything this crate reads — geometry, units,
/// scaling, datatypes — and differ in header layout only (`nifti2.h`:
/// 540 bytes, `i64` dims, `f64` geometry/scaling, `i64` `vox_offset`),
/// so every accessor is version-agnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NiftiVersion {
    /// NIfTI-1 (348-byte header, magic `n+1`/`ni1`).
    V1,
    /// NIfTI-2 (540-byte header, magic `n+2`/`ni2`).
    V2,
}

/// A parsed NIfTI volume (either format version): dimensions, physical
/// geometry (RAS), and scaled per-voxel values.
#[derive(Debug, Clone)]
pub struct NiftiVolume {
    /// Voxel counts `(nx, ny, nz)`.
    pub dims: (usize, usize, usize),
    /// Voxel spacing in mm `(dx, dy, dz)`.
    pub voxel_spacing: (f64, f64, f64),
    /// Position of voxel `(0, 0, 0)`'s center, in RAS mm.
    pub origin: Vec3,
    /// Unit-column rotation: voxel i/j/k axis directions in RAS. Voxel
    /// spacing is *not* baked in here — see [`Self::voxel_position`].
    pub rotation: Mat3,
    /// The on-disk voxel datatype, for informational purposes (e.g.
    /// deciding whether the source precision can represent a downstream
    /// computation exactly). `values` is always `f64` regardless.
    pub datatype: Datatype,
    /// Per-voxel values, already scaled by `scl_slope`/`scl_inter` (or
    /// unscaled if `scl_slope == 0`, per the NIfTI convention). Row-major,
    /// x-fastest: `index = i + nx * (j + ny * k)`.
    pub values: Vec<f64>,
    /// Which NIfTI format version the volume was parsed from.
    pub version: NiftiVersion,
}

impl NiftiVolume {
    /// Parses a NIfTI volume (version 1 or 2) from a file on disk.
    pub fn parse_file(path: &Path) -> Result<Self> {
        let bytes = std::fs::read(path)?;
        Self::parse_bytes(&bytes).map_err(|e| match e {
            NiftiError::NotNifti(_) => NiftiError::NotNifti(path.to_path_buf()),
            NiftiError::Gzipped(_) => NiftiError::Gzipped(path.to_path_buf()),
            other => other,
        })
    }

    /// Parses a NIfTI volume (version 1 or 2) from an in-memory byte
    /// buffer.
    ///
    /// Accepts the single-file `.nii` layout (magic `n+1` or `n+2`); the
    /// parsed volume's [`NiftiVolume::version`] says which. A
    /// gzip-compressed buffer is decompressed when the `gzip` feature is on
    /// and named with [`NiftiError::Gzipped`] when it is off. A dual-file
    /// `.hdr` buffer (magic `ni1`/`ni2`) is rejected with a pointer to
    /// [`Self::parse_dual_bytes`], since its voxel data lives in a separate
    /// `.img`.
    pub fn parse_bytes(bytes: &[u8]) -> Result<Self> {
        if is_gzip(bytes) {
            return gzip_input(bytes);
        }
        // NIfTI-2 announces itself through sizeof_hdr = 540 (and carries
        // its magic at offset 4, not 344), so the version dispatch reads
        // the size field first.
        let is_v2 = bytes.len() >= 4
            && (bytes[0..4] == 540u32.to_le_bytes() || bytes[0..4] == 540u32.to_be_bytes());
        if is_v2 {
            if bytes.len() >= 7 && &bytes[4..7] == b"ni2" {
                return Err(NiftiError::BadValue {
                    reason: "magic is \"ni2\", a dual-file .hdr/.img header — \
                             read it with NiftiVolume::parse_dual_file or \
                             NiftiVolume::parse_dual_bytes"
                        .into(),
                });
            }
            let header = match parse_header2(bytes) {
                Ok(h) => h,
                Err(NiftiError::BadValue { .. }) => {
                    return Err(NiftiError::NotNifti(PathBuf::from("<memory>")))
                }
                Err(e) => return Err(e),
            };
            let fields = decode_header(&header.common(), Layout::Single, NiftiVersion::V2)?;
            if bytes.len() < fields.data_end {
                return Err(NiftiError::UnexpectedEof {
                    offset: bytes.len(),
                    while_reading: "voxel data",
                });
            }
            let values = decode_values(&bytes[fields.data_start..fields.data_end], &fields)?;
            return Ok(Self::from_decoded(fields, values));
        }
        if bytes.len() < HEADER_LEN {
            return Err(NiftiError::NotNifti(PathBuf::from("<memory>")));
        }
        if &bytes[344..347] == b"ni1" {
            return Err(NiftiError::BadValue {
                reason: "magic is \"ni1\", a dual-file .hdr/.img header — \
                         read it with NiftiVolume::parse_dual_file or \
                         NiftiVolume::parse_dual_bytes"
                    .into(),
            });
        }
        let header = match parse_header(bytes) {
            Ok(h) => h,
            Err(NiftiError::BadValue { .. }) => {
                return Err(NiftiError::NotNifti(PathBuf::from("<memory>")))
            }
            Err(e) => return Err(e),
        };
        let fields = decode_header(&header.common(), Layout::Single, NiftiVersion::V1)?;
        if bytes.len() < fields.data_end {
            return Err(NiftiError::UnexpectedEof {
                offset: bytes.len(),
                while_reading: "voxel data",
            });
        }
        let values = decode_values(&bytes[fields.data_start..fields.data_end], &fields)?;
        Ok(Self::from_decoded(fields, values))
    }

    /// Parses a dual-file NIfTI pair (`.hdr` + `.img`, magic `ni1` or
    /// `ni2`) from in-memory buffers.
    ///
    /// `vox_offset` is honoured as an offset into `img_bytes` (0 is the
    /// norm for the two-file layout). With the `gzip` feature on, either
    /// part may be gzip-compressed (`.hdr.gz`/`.img.gz`, which the NIfTI
    /// spec allows); without it a gzipped part is named with
    /// [`NiftiError::Gzipped`].
    pub fn parse_dual_bytes(hdr_bytes: &[u8], img_bytes: &[u8]) -> Result<Self> {
        let hdr = expand_part(hdr_bytes)?;
        let img = expand_part(img_bytes)?;
        let is_v2 = hdr.len() >= 4
            && (hdr[0..4] == 540u32.to_le_bytes() || hdr[0..4] == 540u32.to_be_bytes());
        if !is_v2 && hdr.len() < HEADER_LEN {
            return Err(NiftiError::NotNifti(PathBuf::from("<memory>")));
        }
        if is_v2 {
            if hdr.len() >= 7 && &hdr[4..7] == b"n+2" {
                return Err(NiftiError::BadValue {
                    reason: "magic is \"n+2\", a single-file .nii header — read \
                             it with NiftiVolume::parse_file or \
                             NiftiVolume::parse_bytes"
                        .into(),
                });
            }
            let header = match parse_header2_dual(&hdr) {
                Ok(h) => h,
                Err(NiftiError::BadValue { .. }) => {
                    return Err(NiftiError::NotNifti(PathBuf::from("<memory>")))
                }
                Err(e) => return Err(e),
            };
            let fields = decode_header(&header.common(), Layout::Dual, NiftiVersion::V2)?;
            if img.len() < fields.data_end {
                return Err(NiftiError::UnexpectedEof {
                    offset: img.len(),
                    while_reading: "dual-file voxel data (.img)",
                });
            }
            let values = decode_values(&img[fields.data_start..fields.data_end], &fields)?;
            return Ok(Self::from_decoded(fields, values));
        }
        if &hdr[344..347] == b"n+1" {
            return Err(NiftiError::BadValue {
                reason: "magic is \"n+1\", a single-file .nii header — read \
                         it with NiftiVolume::parse_file or \
                         NiftiVolume::parse_bytes"
                    .into(),
            });
        }
        let header = match parse_header_dual(&hdr) {
            Ok(h) => h,
            Err(NiftiError::BadValue { .. }) => {
                return Err(NiftiError::NotNifti(PathBuf::from("<memory>")))
            }
            Err(e) => return Err(e),
        };
        let fields = decode_header(&header.common(), Layout::Dual, NiftiVersion::V1)?;
        if img.len() < fields.data_end {
            return Err(NiftiError::UnexpectedEof {
                offset: img.len(),
                while_reading: "dual-file voxel data (.img)",
            });
        }
        let values = decode_values(&img[fields.data_start..fields.data_end], &fields)?;
        Ok(Self::from_decoded(fields, values))
    }

    /// Parses a dual-file NIfTI-1 pair from disk: `hdr_path` is the `.hdr`
    /// header file, `img_path` the `.img` voxel-data file. Read errors name
    /// the file that caused them; without the `gzip` feature a compressed
    /// part is reported as [`NiftiError::Gzipped`] against its own path.
    pub fn parse_dual_file(hdr_path: &Path, img_path: &Path) -> Result<Self> {
        let hdr = std::fs::read(hdr_path)?;
        let img = std::fs::read(img_path)?;
        // Without the `gzip` feature, name which file is compressed — the
        // in-buffer fallback inside `parse_dual_bytes` can only say
        // "<memory>".
        #[cfg(not(feature = "gzip"))]
        {
            if is_gzip(&hdr) {
                return Err(NiftiError::Gzipped(hdr_path.to_path_buf()));
            }
            if is_gzip(&img) {
                return Err(NiftiError::Gzipped(img_path.to_path_buf()));
            }
        }
        Self::parse_dual_bytes(&hdr, &img).map_err(|e| match e {
            NiftiError::NotNifti(_) => NiftiError::NotNifti(hdr_path.to_path_buf()),
            other => other,
        })
    }

    /// Decompresses a gzip-wrapped single-file volume and parses it.
    ///
    /// Inflation happens in two phases so the header's own geometry caps
    /// the allocation: phase 1 inflates exactly the 348-byte header, phase
    /// 2 inflates only up to the declared `data_end` (`vox_offset` +
    /// dimensions x datatype) and stops there. A corrupt or hostile stream
    /// can therefore never make this reader allocate more than the file's
    /// own header says the volume needs — the same geometry-capped posture
    /// `tpt-med-dicom`'s RLE decoder takes against an expanding run — and a
    /// stream shorter than its header claims fails as
    /// [`NiftiError::UnexpectedEof`] rather than being padded.
    #[cfg(feature = "gzip")]
    fn parse_gzip(bytes: &[u8]) -> Result<Self> {
        use std::io::Read;

        let mut decoder = flate2::read::GzDecoder::new(bytes);

        // Phase 1: the header, sized by its own sizeof_hdr (348 or 540),
        // read incrementally so the version is known before its header
        // length is demanded.
        let mut buf = vec![0u8; 4];
        let mut filled = 0usize;
        while filled < 4 {
            match decoder.read(&mut buf[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(NiftiError::CorruptGzip(e)),
            }
        }
        let header_len = {
            let le = i32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
            let be = i32::from_be_bytes([buf[0], buf[1], buf[2], buf[3]]);
            // An unrecognised size is not an error here: defaulting to
            // the NIfTI-1 length lets the version-1 parse run and classify
            // the payload as NotNifti, exactly as the uncompressed path
            // does for garbage input.
            match (le, be) {
                (540, _) | (_, 540) => HEADER2_LEN,
                _ => HEADER_LEN,
            }
        };
        buf.resize(header_len, 0);
        while filled < header_len {
            match decoder.read(&mut buf[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(NiftiError::CorruptGzip(e)),
            }
        }
        if filled < header_len {
            return Err(NiftiError::UnexpectedEof {
                offset: filled,
                while_reading: "gzip-compressed NIfTI header",
            });
        }
        let is_v2 = header_len == HEADER2_LEN;
        if &buf[344..347] == b"ni1" {
            return Err(NiftiError::BadValue {
                reason: "gzip stream holds a dual-file .hdr — read both \
                         parts with NiftiVolume::parse_dual_file or \
                         NiftiVolume::parse_dual_bytes (with the `gzip` \
                         feature for .hdr.gz/.img.gz)"
                    .into(),
            });
        }
        if is_v2 && buf.len() >= 7 && &buf[4..7] == b"ni2" {
            return Err(NiftiError::BadValue {
                reason: "gzip stream holds a dual-file .hdr — read both                          parts with NiftiVolume::parse_dual_file or                          NiftiVolume::parse_dual_bytes (with the `gzip`                          feature for .hdr.gz/.img.gz)"
                    .into(),
            });
        }
        let fields = if is_v2 {
            let header = match parse_header2(&buf) {
                Ok(h) => h,
                Err(NiftiError::BadValue { .. }) => {
                    return Err(NiftiError::NotNifti(PathBuf::from("<memory>")))
                }
                Err(e) => return Err(e),
            };
            decode_header(&header.common(), Layout::Single, NiftiVersion::V2)?
        } else {
            let header = match parse_header(&buf) {
                Ok(h) => h,
                Err(NiftiError::BadValue { .. }) => {
                    return Err(NiftiError::NotNifti(PathBuf::from("<memory>")))
                }
                Err(e) => return Err(e),
            };
            decode_header(&header.common(), Layout::Single, NiftiVersion::V1)?
        };

        // Phase 2: the rest, grown chunk by chunk up to `data_end` — the
        // buffer follows the bytes actually inflated, never the (header-
        // chosen) declared extent.
        let mut chunk = [0u8; 64 * 1024];
        while buf.len() < fields.data_end {
            let want = (fields.data_end - buf.len()).min(chunk.len());
            match decoder.read(&mut chunk[..want]) {
                Ok(0) => {
                    return Err(NiftiError::UnexpectedEof {
                        offset: buf.len(),
                        while_reading: "gzip-compressed voxel data",
                    })
                }
                Ok(n) => buf.extend_from_slice(&chunk[..n]),
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                    return Err(NiftiError::UnexpectedEof {
                        offset: buf.len(),
                        while_reading: "gzip-compressed voxel data",
                    })
                }
                Err(e) => return Err(NiftiError::CorruptGzip(e)),
            }
        }
        let values = decode_values(&buf[fields.data_start..fields.data_end], &fields)?;
        Ok(Self::from_decoded(fields, values))
    }

    /// Assembles a volume from decoded header fields and the value vector
    /// every parse path ends with.
    fn from_decoded(fields: DecodedHeader, values: Vec<f64>) -> Self {
        Self {
            dims: fields.dims,
            voxel_spacing: fields.spacing,
            origin: fields.origin,
            rotation: fields.rotation,
            datatype: fields.datatype,
            values,
            version: fields.version,
        }
    }

    /// The value at voxel `(i, j, k)`, or `None` if out of bounds.
    pub fn value_at(&self, i: usize, j: usize, k: usize) -> Option<f64> {
        let (nx, ny, nz) = self.dims;
        if i >= nx || j >= ny || k >= nz {
            return None;
        }
        self.values.get(i + nx * (j + ny * k)).copied()
    }

    /// RAS-mm position of voxel `(i, j, k)`'s center.
    pub fn voxel_position(&self, i: usize, j: usize, k: usize) -> Vec3 {
        self.origin
            + self.rotation.col(0) * (i as f64 * self.voxel_spacing.0)
            + self.rotation.col(1) * (j as f64 * self.voxel_spacing.1)
            + self.rotation.col(2) * (k as f64 * self.voxel_spacing.2)
    }
}

/// Whether `bytes` announce themselves as gzip (`1F 8B` magic).
fn is_gzip(bytes: &[u8]) -> bool {
    bytes.len() >= 2 && bytes[0] == 0x1F && bytes[1] == 0x8B
}

/// Handles a gzip-magic input: decompress and parse when the `gzip`
/// feature is on, produce the typed [`NiftiError::Gzipped`] when it is off.
fn gzip_input(bytes: &[u8]) -> Result<NiftiVolume> {
    #[cfg(feature = "gzip")]
    {
        NiftiVolume::parse_gzip(bytes)
    }
    #[cfg(not(feature = "gzip"))]
    {
        let _ = bytes; // only the typed "enable `gzip`" error exists here
        Err(NiftiError::Gzipped(PathBuf::from("<memory>")))
    }
}

/// Returns `bytes` decompressed when they are gzip, or as-is when not.
/// With the `gzip` feature off, a gzip input produces
/// [`NiftiError::Gzipped`]. This is the dual-file path's helper: both parts
/// are already in memory, so unlike the single-file path there is no
/// streaming here — a `.hdr` is a few hundred bytes, and a `.img.gz` is
/// bounded by DEFLATE's maximum 1032:1 expansion ratio against the buffer
/// the caller passed in (whose size the caller already chose to pay for).
fn expand_part(bytes: &[u8]) -> Result<Cow<'_, [u8]>> {
    if !is_gzip(bytes) {
        return Ok(Cow::Borrowed(bytes));
    }
    #[cfg(feature = "gzip")]
    {
        use std::io::Read;
        let mut out = Vec::new();
        flate2::read::GzDecoder::new(bytes)
            .read_to_end(&mut out)
            .map_err(NiftiError::CorruptGzip)?;
        Ok(Cow::Owned(out))
    }
    #[cfg(not(feature = "gzip"))]
    {
        Err(NiftiError::Gzipped(PathBuf::from("<memory>")))
    }
}

/// Which file layout the header came from — it changes only what
/// `vox_offset` means: an offset into the same buffer (single-file), or an
/// offset into a separate `.img` buffer, where 0 is the norm (dual-file).
#[derive(Clone, Copy)]
enum Layout {
    Single,
    Dual,
}

/// Everything the header says about the volume except the values
/// themselves — shared by the single-file, dual-file, and gzip paths so the
/// three cannot drift apart on dimension, datatype, geometry, or offset
/// checks.
struct DecodedHeader {
    dims: (usize, usize, usize),
    datatype: Datatype,
    spacing: (f64, f64, f64),
    origin: Vec3,
    rotation: Mat3,
    /// First voxel-data byte in its buffer.
    data_start: usize,
    /// One past the last voxel-data byte in its buffer.
    data_end: usize,
    big_endian: bool,
    scl_slope: f64,
    scl_inter: f64,
    version: NiftiVersion,
}

/// Decodes the geometry- and layout-relevant fields of a parsed header.
/// The header values both NIfTI-1 and NIfTI-2 decode from. The two
/// versions' fields carry identical semantics at different offsets and
/// widths (NIfTI-2 widens `dim` to `i64` and the geometry/scaling to
/// `f64`); normalising once here means the decode logic below — checks,
/// geometry recovery, extent arithmetic — is written once and cannot
/// drift between the versions.
pub(crate) struct CommonHeader {
    pub dim: [i64; 8],
    pub datatype_code: i16,
    pub bitpix: i16,
    pub pixdim: [f64; 8],
    pub vox_offset: f64,
    pub scl_slope: f64,
    pub scl_inter: f64,
    pub qform_code: i32,
    pub sform_code: i32,
    pub quatern_bcd: (f64, f64, f64),
    pub qoffset: (f64, f64, f64),
    pub srow: [[f64; 4]; 3],
    pub big_endian: bool,
}

fn decode_header(
    header: &CommonHeader,
    layout: Layout,
    version: NiftiVersion,
) -> Result<DecodedHeader> {
    let dim = header.dim;
    let ndim = dim[0].max(0) as usize;
    if ndim < 3 {
        return Err(NiftiError::BadValue {
            reason: format!("dim[0] = {ndim}, need at least 3 spatial dimensions"),
        });
    }
    let (nx, ny, nz) = (dim[1], dim[2], dim[3]);
    if nx <= 0 || ny <= 0 || nz <= 0 {
        return Err(NiftiError::BadValue {
            reason: format!("non-positive dimension in dim[1..4] = {nx},{ny},{nz}"),
        });
    }
    let (nx, ny, nz) = (nx as usize, ny as usize, nz as usize);
    // Checked product: NIfTI-1's i16 dimensions fit usize on 64-bit but can
    // overflow wasm32, where an unchecked product would wrap and then have
    // its (wrapped) data_end validated against a real buffer length.
    let n_voxels = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(nz))
        .ok_or_else(|| NiftiError::BadValue {
            reason: format!("dimensions {nx}x{ny}x{nz} overflow usize"),
        })?;

    let datatype_code = header.datatype_code;
    let datatype =
        Datatype::from_code(datatype_code).ok_or(NiftiError::UnsupportedDatatype(datatype_code))?;
    let bitpix = header.bitpix;
    if bitpix != datatype.bitpix() {
        return Err(NiftiError::BitpixMismatch {
            datatype: datatype_code,
            bitpix,
        });
    }

    let pixdim = header.pixdim;
    let qfac = if pixdim[0] == 0.0 { 1.0 } else { pixdim[0] };
    let (spacing, origin, rotation) = if header.sform_code > 0 {
        geometry_from_srow(&header.srow)?
    } else if header.qform_code > 0 {
        geometry_from_qform(header, qfac)
    } else {
        (
            (pixdim[1], pixdim[2], pixdim[3]),
            Vec3::ZERO,
            Mat3::IDENTITY,
        )
    };

    let vox_offset = header.vox_offset;
    let header_len = match version {
        NiftiVersion::V1 => HEADER_LEN,
        NiftiVersion::V2 => HEADER2_LEN,
    };
    let data_start = match layout {
        Layout::Single => {
            if !vox_offset.is_finite() || vox_offset < header_len as f64 {
                return Err(NiftiError::BadValue {
                    reason: format!(
                        "vox_offset {vox_offset} is before the end of the \
                         {header_len}-byte header"
                    ),
                });
            }
            vox_offset as usize
        }
        Layout::Dual => {
            if !vox_offset.is_finite() || vox_offset < 0.0 {
                return Err(NiftiError::BadValue {
                    reason: format!(
                        "vox_offset {vox_offset} is invalid; in the dual-file \
                         layout it is an offset into the .img file (0 is the norm)"
                    ),
                });
            }
            vox_offset as usize
        }
    };
    let data_end = n_voxels
        .checked_mul(datatype.size_bytes())
        .and_then(|n| data_start.checked_add(n))
        .ok_or_else(|| NiftiError::BadValue {
            reason: "voxel data extent overflows usize".into(),
        })?;

    Ok(DecodedHeader {
        dims: (nx, ny, nz),
        datatype,
        spacing,
        origin,
        rotation,
        data_start,
        data_end,
        big_endian: header.big_endian,
        scl_slope: header.scl_slope,
        scl_inter: header.scl_inter,
        version,
    })
}

/// Decodes `voxel_bytes` — sliced by the caller to start at
/// `fields.data_start` — into scaled `f64` values.
fn decode_values(voxel_bytes: &[u8], fields: &DecodedHeader) -> Result<Vec<f64>> {
    let (nx, ny, nz) = fields.dims;
    let n_voxels = nx * ny * nz;
    let sample = fields.datatype.size_bytes();
    // Callers validate against `data_end` first; this keeps the decode loop
    // itself panic-free regardless.
    if voxel_bytes.len() < n_voxels * sample {
        return Err(NiftiError::UnexpectedEof {
            offset: voxel_bytes.len(),
            while_reading: "voxel data",
        });
    }
    let mut values = Vec::with_capacity(n_voxels);
    for i in 0..n_voxels {
        let start = i * sample;
        let raw = fields
            .datatype
            .decode(&voxel_bytes[start..start + sample], fields.big_endian);
        let value = if fields.scl_slope != 0.0 {
            raw * fields.scl_slope + fields.scl_inter
        } else {
            raw
        };
        values.push(value);
    }
    Ok(values)
}

/// Builds `(spacing, origin, rotation)` from the sform affine (`srow_x/y/z`):
/// each row directly gives `ras_axis = row[0..3]·(i,j,k) + row[3]`. The
/// affine mixes rotation and spacing together (unlike DICOM's separate
/// direction-cosine and spacing tags), so spacing is recovered as each
/// column's norm and the rotation as that column normalised — an explicit
/// no-shear assumption, see RFC 0006 Drawbacks.
fn geometry_from_srow(srow: &[[f64; 4]; 3]) -> Result<((f64, f64, f64), Vec3, Mat3)> {
    let (sx, sy, sz) = (srow[0], srow[1], srow[2]);
    let origin = Vec3::new(sx[3], sy[3], sz[3]);

    let cols = [
        Vec3::new(sx[0], sy[0], sz[0]),
        Vec3::new(sx[1], sy[1], sz[1]),
        Vec3::new(sx[2], sy[2], sz[2]),
    ];
    let mut spacing = [0.0f64; 3];
    let mut units = [Vec3::ZERO; 3];
    for (i, c) in cols.iter().enumerate() {
        let n = c.norm();
        if n < 1e-9 {
            return Err(NiftiError::BadValue {
                reason: format!("sform column {i} has near-zero norm; degenerate affine"),
            });
        }
        spacing[i] = n;
        units[i] = *c / n;
    }
    Ok((
        (spacing[0], spacing[1], spacing[2]),
        origin,
        Mat3::from_cols(units[0], units[1], units[2]),
    ))
}

/// Builds `(spacing, origin, rotation)` from the qform quaternion.
fn geometry_from_qform(header: &CommonHeader, qfac: f64) -> ((f64, f64, f64), Vec3, Mat3) {
    let (b, c, d) = header.quatern_bcd;
    let rotation = geometry::quaternion_to_rotation(b, c, d);
    let (ox, oy, oz) = header.qoffset;
    let pixdim = header.pixdim;

    // qfac flips the k (third) column's handedness; magnitude spacing stays
    // positive, matching pixdim's own always-positive convention.
    let k_col = if qfac < 0.0 {
        rotation.col(2) * -1.0
    } else {
        rotation.col(2)
    };
    let rotation = Mat3::from_cols(rotation.col(0), rotation.col(1), k_col);

    (
        (pixdim[1], pixdim[2], pixdim[3]),
        Vec3::new(ox, oy, oz),
        rotation,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synthetic::SyntheticNiftiBuilder;

    #[test]
    fn round_trips_sform_uint8_volume() {
        let (nx, ny, nz) = (4usize, 3, 2);
        let values: Vec<f64> = (0..nx * ny * nz).map(|i| (i % 256) as f64).collect();
        let bytes = SyntheticNiftiBuilder::new(nx, ny, nz, (1.0, 1.0, 2.0))
            .with_values_u8(&values)
            .build();
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        assert_eq!(vol.dims, (nx, ny, nz));
        assert_eq!(vol.values, values);
    }

    #[test]
    fn round_trips_float32_with_scaling() {
        let (nx, ny, nz) = (2usize, 2, 2);
        let raw: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let bytes = SyntheticNiftiBuilder::new(nx, ny, nz, (0.5, 0.5, 0.5))
            .with_scaling(2.0, -1000.0)
            .with_values_f32(&raw)
            .build();
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        let expected: Vec<f64> = raw.iter().map(|&v| v * 2.0 - 1000.0).collect();
        for (got, want) in vol.values.iter().zip(&expected) {
            assert!((got - want).abs() < 1e-4);
        }
    }

    #[test]
    fn sform_geometry_recovers_spacing_and_origin() {
        let bytes = SyntheticNiftiBuilder::new(2, 2, 2, (1.5, 2.5, 3.5))
            .with_origin(Vec3::new(10.0, -20.0, 30.0))
            .with_values_u8(&[0.0; 8])
            .build();
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        assert!((vol.voxel_spacing.0 - 1.5).abs() < 1e-4);
        assert!((vol.voxel_spacing.1 - 2.5).abs() < 1e-4);
        assert!((vol.voxel_spacing.2 - 3.5).abs() < 1e-4);
        assert!((vol.origin - Vec3::new(10.0, -20.0, 30.0)).norm() < 1e-4);
    }

    #[cfg(not(feature = "gzip"))]
    #[test]
    fn gzip_magic_is_rejected_with_named_error() {
        let err = NiftiVolume::parse_bytes(&[0x1F, 0x8B, 0x08, 0x00]).unwrap_err();
        assert!(matches!(err, NiftiError::Gzipped(_)));
    }

    #[cfg(not(feature = "gzip"))]
    #[test]
    fn gzipped_dual_file_part_is_named_without_the_feature() {
        let err = NiftiVolume::parse_dual_bytes(&[0x1F, 0x8B, 0x08, 0x00], &[]).unwrap_err();
        assert!(matches!(err, NiftiError::Gzipped(_)));
    }

    #[test]
    fn garbage_is_rejected_not_panicking() {
        let err = NiftiVolume::parse_bytes(&[0u8; 400]).unwrap_err();
        assert!(matches!(err, NiftiError::NotNifti(_)));
    }

    #[test]
    fn truncated_header_is_rejected() {
        let err = NiftiVolume::parse_bytes(&[0u8; 100]).unwrap_err();
        assert!(matches!(err, NiftiError::NotNifti(_)));
    }

    #[test]
    fn unsupported_datatype_is_rejected() {
        let mut bytes = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_u8(&[0.0; 8])
            .build();
        // Overwrite datatype (offset 70) with DT_COMPLEX64 (32), keep bitpix
        // matching so this exercises the datatype check specifically.
        bytes[70..72].copy_from_slice(&32i16.to_le_bytes());
        bytes[72..74].copy_from_slice(&64i16.to_le_bytes());
        let err = NiftiVolume::parse_bytes(&bytes).unwrap_err();
        assert!(matches!(err, NiftiError::UnsupportedDatatype(32)));
    }

    #[test]
    fn bitpix_mismatch_is_rejected() {
        let mut bytes = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_u8(&[0.0; 8])
            .build();
        // datatype stays DT_UINT8 (expects bitpix 8); corrupt bitpix to 16.
        bytes[72..74].copy_from_slice(&16i16.to_le_bytes());
        let err = NiftiVolume::parse_bytes(&bytes).unwrap_err();
        assert!(matches!(
            err,
            NiftiError::BitpixMismatch {
                datatype: 2,
                bitpix: 16
            }
        ));
    }

    #[test]
    fn truncated_voxel_data_is_rejected_not_padded() {
        let bytes = SyntheticNiftiBuilder::new(4, 4, 4, (1.0, 1.0, 1.0))
            .with_values_u8(&vec![0.0; 64])
            .build();
        let truncated = &bytes[..bytes.len() - 10];
        let err = NiftiVolume::parse_bytes(truncated).unwrap_err();
        assert!(matches!(err, NiftiError::UnexpectedEof { .. }));
    }

    #[test]
    fn qform_path_is_used_when_sform_absent() {
        // 180-degree rotation about k (b=0,c=0,d=1), qfac=1: rotation should
        // come out as diag(-1,-1,1) per geometry.rs's hand-derived case.
        let bytes = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_origin(Vec3::new(5.0, 6.0, 7.0))
            .with_qform((0.0, 0.0, 1.0), 1.0)
            .with_values_u8(&[0.0; 8])
            .build();
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        assert!((vol.origin - Vec3::new(5.0, 6.0, 7.0)).norm() < 1e-4);
        assert_eq!(vol.rotation, Mat3::diagonal([-1.0, -1.0, 1.0]));
    }

    #[test]
    fn qform_qfac_negative_flips_k_column() {
        let bytes = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_qform((0.0, 0.0, 0.0), -1.0) // identity quaternion, qfac=-1
            .with_values_u8(&[0.0; 8])
            .build();
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        assert_eq!(vol.rotation, Mat3::diagonal([1.0, 1.0, -1.0]));
    }

    #[test]
    fn falls_back_to_axis_aligned_when_no_form_present() {
        // Neither with_qform nor a real sform_code -- overwrite sform_code
        // (offset 254) to 0 after building, so both codes are absent.
        let mut bytes = SyntheticNiftiBuilder::new(2, 2, 2, (2.0, 3.0, 4.0))
            .with_values_u8(&[0.0; 8])
            .build();
        bytes[254..256].copy_from_slice(&0i16.to_le_bytes());
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        assert_eq!(vol.rotation, Mat3::IDENTITY);
        assert!((vol.origin - Vec3::ZERO).norm() < 1e-9);
        assert!((vol.voxel_spacing.0 - 2.0).abs() < 1e-4);
        assert!((vol.voxel_spacing.2 - 4.0).abs() < 1e-4);
    }

    #[test]
    fn datatype_field_reflects_on_disk_type() {
        let bytes = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_f32(&[0.0; 8])
            .build();
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        assert_eq!(vol.datatype, Datatype::Float32);
    }

    #[test]
    fn voxel_position_uses_spacing_and_origin() {
        let bytes = SyntheticNiftiBuilder::new(2, 2, 2, (2.0, 3.0, 4.0))
            .with_origin(Vec3::new(0.0, 0.0, 0.0))
            .with_values_u8(&[0.0; 8])
            .build();
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        let p = vol.voxel_position(1, 1, 1);
        assert!((p - Vec3::new(2.0, 3.0, 4.0)).norm() < 1e-4);
    }

    /// Field-by-field comparison with float tolerance (the fixtures round-
    /// trip through `f32` storage, hence 1e-4 on values).
    fn assert_volumes_eq(a: &NiftiVolume, b: &NiftiVolume) {
        assert_eq!(a.dims, b.dims);
        assert_eq!(a.datatype, b.datatype);
        assert!((a.voxel_spacing.0 - b.voxel_spacing.0).abs() < 1e-6);
        assert!((a.voxel_spacing.1 - b.voxel_spacing.1).abs() < 1e-6);
        assert!((a.voxel_spacing.2 - b.voxel_spacing.2).abs() < 1e-6);
        assert!((a.origin - b.origin).norm() < 1e-6);
        for c in 0..3 {
            assert!((a.rotation.col(c) - b.rotation.col(c)).norm() < 1e-6);
        }
        assert_eq!(a.values.len(), b.values.len());
        for (x, y) in a.values.iter().zip(&b.values) {
            assert!((x - y).abs() < 1e-4, "{x} != {y}");
        }
    }

    #[test]
    fn round_trips_dual_file_layout_like_the_single_file_one() {
        let (nx, ny, nz) = (4usize, 3, 2);
        let values: Vec<f64> = (0..nx * ny * nz).map(|i| (i % 256) as f64).collect();
        let build = || {
            SyntheticNiftiBuilder::new(nx, ny, nz, (1.5, 2.5, 3.5))
                .with_origin(Vec3::new(10.0, -20.0, 30.0))
                .with_values_u8(&values)
        };
        let single = NiftiVolume::parse_bytes(&build().build()).expect("single-file parses");
        let (hdr, img) = build().build_dual();
        let dual = NiftiVolume::parse_dual_bytes(&hdr, &img).expect("dual-file parses");
        assert_eq!(dual.dims, (nx, ny, nz));
        assert_eq!(dual.values, values);
        assert_volumes_eq(&single, &dual);
    }

    #[test]
    fn dual_file_honours_a_nonzero_vox_offset_into_the_img() {
        let (mut hdr, img) = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_u8(&[7.0; 8])
            .build_dual();
        // vox_offset lives at byte 108: 4 bytes of padding before the data.
        hdr[108..112].copy_from_slice(&4.0f32.to_le_bytes());
        let padded = [vec![0u8; 4], img].concat();
        let vol = NiftiVolume::parse_dual_bytes(&hdr, &padded).expect("parses");
        assert_eq!(vol.values, vec![7.0; 8]);
    }

    #[test]
    fn dual_file_header_shorter_than_the_header_is_rejected() {
        let err = NiftiVolume::parse_dual_bytes(&[0u8; 100], &[]).unwrap_err();
        assert!(matches!(err, NiftiError::NotNifti(_)));
    }

    #[test]
    fn dual_file_img_shorter_than_its_geometry_claims_is_rejected() {
        let (hdr, img) = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_u8(&[7.0; 8])
            .build_dual();
        let err = NiftiVolume::parse_dual_bytes(&hdr, &img[..img.len() - 3]).unwrap_err();
        assert!(matches!(err, NiftiError::UnexpectedEof { .. }));
    }

    #[test]
    fn single_file_entry_point_points_at_the_dual_one() {
        let (hdr, _) = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_u8(&[0.0; 8])
            .build_dual();
        match NiftiVolume::parse_bytes(&hdr).unwrap_err() {
            NiftiError::BadValue { reason } => {
                assert!(reason.contains("parse_dual"), "reason: {reason}");
            }
            other => panic!("expected BadValue, got {other:?}"),
        }
    }

    #[test]
    fn dual_file_entry_point_points_at_the_single_file_one() {
        let bytes = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_u8(&[0.0; 8])
            .build();
        match NiftiVolume::parse_dual_bytes(&bytes, &[]).unwrap_err() {
            NiftiError::BadValue { reason } => {
                assert!(reason.contains("parse_bytes"), "reason: {reason}");
            }
            other => panic!("expected BadValue, got {other:?}"),
        }
    }

    // ---- gzip feature (off by default; these need `--features gzip`) ----

    #[cfg(feature = "gzip")]
    fn gzip(data: &[u8]) -> Vec<u8> {
        use std::io::Write;
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(data).expect("gzip write");
        encoder.finish().expect("gzip finish")
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn gzip_volume_parses_to_the_same_result_as_the_plain_one() {
        let values: Vec<f64> = (1..=24).map(|i| i as f64).collect();
        let bytes = SyntheticNiftiBuilder::new(4, 3, 2, (1.5, 2.5, 3.5))
            .with_origin(Vec3::new(10.0, -20.0, 30.0))
            .with_scaling(2.0, -1000.0)
            .with_values_f32(&values)
            .build();
        let plain = NiftiVolume::parse_bytes(&bytes).expect("plain parses");
        let from_gz = NiftiVolume::parse_bytes(&gzip(&bytes)).expect("gzip parses");
        assert_volumes_eq(&plain, &from_gz);
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn gzip_of_a_truncated_volume_is_rejected_not_padded() {
        let bytes = SyntheticNiftiBuilder::new(4, 4, 4, (1.0, 1.0, 1.0))
            .with_values_u8(&vec![0.0; 64])
            .build();
        let err = NiftiVolume::parse_bytes(&gzip(&bytes[..bytes.len() - 10])).unwrap_err();
        assert!(
            matches!(err, NiftiError::UnexpectedEof { .. }),
            "got {err:?}"
        );
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn corrupt_gzip_stream_is_named_not_misparsed() {
        // Valid magic, reserved FLG bits set (0xF8): not a decodable gzip.
        let err =
            NiftiVolume::parse_bytes(&[0x1F, 0x8B, 0x08, 0xF8, 0, 0, 0, 0, 0, 0]).unwrap_err();
        assert!(matches!(err, NiftiError::CorruptGzip(_)), "got {err:?}");
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn gzip_of_a_non_nifti_payload_is_rejected_not_misparsed() {
        let err = NiftiVolume::parse_bytes(&gzip(&[0u8; 400])).unwrap_err();
        assert!(matches!(err, NiftiError::NotNifti(_)), "got {err:?}");
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn dual_file_parts_may_be_gzipped() {
        let (hdr, img) = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_u8(&[5.0; 8])
            .build_dual();
        let vol = NiftiVolume::parse_dual_bytes(&gzip(&hdr), &gzip(&img)).expect("parses");
        assert_eq!(vol.values, vec![5.0; 8]);
    }
}

#[cfg(test)]
mod nifti2_tests {
    use super::*;
    use crate::synthetic::SyntheticNiftiBuilder;

    /// A non-f32-exact origin: NIfTI-2 stores geometry as f64, so 1/3
    /// round-trips exactly where the NIfTI-1 path (f32 srows) could not.
    const THIRD: f64 = 1.0 / 3.0;

    #[test]
    fn round_trips_nifti2_sform_volume() {
        let (nx, ny, nz) = (4usize, 3, 2);
        let values: Vec<f64> = (0..nx * ny * nz).map(|i| (i % 256) as f64).collect();
        let bytes = SyntheticNiftiBuilder::new(nx, ny, nz, (1.0, 1.0, 2.0))
            .with_origin(Vec3::new(THIRD, -2.0, 5.0))
            .with_values_u8(&values)
            .build2();
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        assert_eq!(vol.version, NiftiVersion::V2);
        assert_eq!(vol.dims, (nx, ny, nz));
        assert_eq!(vol.datatype, Datatype::Uint8);
        assert_eq!(vol.voxel_spacing, (1.0, 1.0, 2.0));
        // The f64 srow recovers 1/3 exactly; a f32 header could not.
        assert_eq!(vol.origin.x, THIRD);
        assert_eq!(vol.origin.y, -2.0);
        assert_eq!(vol.origin.z, 5.0);
        assert_eq!(vol.values, values);
        assert_eq!(
            vol.value_at(nx - 1, ny - 1, nz - 1),
            Some(&values[values.len() - 1]).copied()
        );
    }

    #[test]
    fn round_trips_nifti2_dual_and_qform() {
        let values: Vec<f64> = (0..24).map(|i| i as f64).collect();
        let (hdr, img) = SyntheticNiftiBuilder::new(4, 3, 2, (1.5, 1.5, 2.5))
            .with_qform((0.0, 0.0, 0.0), 1.0)
            .with_scaling(2.0, -10.0)
            .with_values_f32(&values)
            .build2_dual();
        let vol = NiftiVolume::parse_dual_bytes(&hdr, &img).expect("parses");
        assert_eq!(vol.version, NiftiVersion::V2);
        // Scaling applies through the f64 scl fields: parsed = raw·2 − 10.
        let expected: Vec<f64> = values
            .iter()
            .map(|&v| (v as f32 * 2.0 - 10.0) as f64)
            .collect();
        assert_eq!(vol.values, expected);
        assert_eq!(vol.voxel_spacing, (1.5, 1.5, 2.5));
    }

    #[test]
    fn nifti2_i64_dimensions_exceed_the_nifti1_range() {
        // dim as i16 (NIfTI-1) tops out at 32767; NIfTI-2's i64 dims take
        // 40000 voxels along x. A u8 volume keeps the fixture at 160 KB.
        let (nx, ny, nz) = (40000usize, 2, 2);
        let values: Vec<f64> = (0..nx * ny * nz).map(|i| (i % 251) as f64).collect();
        let bytes = SyntheticNiftiBuilder::new(nx, ny, nz, (1.0, 1.0, 1.0))
            .with_values_u8(&values)
            .build2();
        let vol = NiftiVolume::parse_bytes(&bytes).expect("parses");
        assert_eq!(vol.version, NiftiVersion::V2);
        assert_eq!(vol.dims, (nx, ny, nz));
        assert_eq!(vol.value_at(nx - 1, 0, 0), Some((39999 % 251) as f64));
    }

    #[test]
    fn nifti2_wrong_entry_points_point_at_the_right_one() {
        let bytes = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_u8(&[0.0; 8])
            .build2();
        let (hdr, img) = SyntheticNiftiBuilder::new(2, 2, 2, (1.0, 1.0, 1.0))
            .with_values_u8(&[0.0; 8])
            .build2_dual();
        // A single-file NIfTI-2 buffer through the dual entry point...
        let err =
            NiftiVolume::parse_dual_bytes(&bytes, &img).expect_err("n+2 through parse_dual_bytes");
        assert!(
            matches!(&err, NiftiError::BadValue { reason } if reason.contains("parse_bytes")),
            "unexpected: {err:?}"
        );
        // ...and a dual-file NIfTI-2 header through the single entry point.
        let err = NiftiVolume::parse_bytes(&hdr).expect_err("ni2 through parse_bytes");
        assert!(
            matches!(&err, NiftiError::BadValue { reason } if reason.contains("parse_dual")),
            "unexpected: {err:?}"
        );
    }
}
