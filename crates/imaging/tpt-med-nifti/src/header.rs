//! NIfTI-1 fixed 348-byte header layout (`nifti1.h`), byte-offset access.
//!
//! Only the fields this crate actually uses are decoded — see
//! `rfcs/0006-nifti-ingestion.md` for the full field table and why the rest
//! (`descrip`, `aux_file`, `cal_min/max`, `intent_*`, slice timing, …) is
//! deliberately skipped.

use crate::error::{NiftiError, Result};

/// Header size in bytes; also `sizeof_hdr`'s required value.
pub const HEADER_LEN: usize = 348;

/// A little/big-endian `i16`/`i32`/`f32` reader over a fixed-size header
/// buffer, with the endianness resolved once by [`parse_header`] /
/// [`parse_header_dual`] from `sizeof_hdr`.
pub struct Header<'a> {
    bytes: &'a [u8],
    big_endian: bool,
}

/// NIfTI datatype codes this crate decodes, and their expected `bitpix`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Datatype {
    /// `DT_UINT8` (2), 8 bits.
    Uint8,
    /// `DT_INT16` (4), 16 bits.
    Int16,
    /// `DT_INT32` (8), 32 bits.
    Int32,
    /// `DT_FLOAT32` (16), 32 bits.
    Float32,
    /// `DT_FLOAT64` (64), 64 bits.
    Float64,
    /// `DT_INT8` (256), 8 bits.
    Int8,
    /// `DT_UINT16` (512), 16 bits.
    Uint16,
    /// `DT_UINT32` (768), 32 bits.
    Uint32,
}

impl Datatype {
    /// Decodes the raw `datatype` code, or `None` for one this crate does
    /// not support.
    pub fn from_code(code: i16) -> Option<Self> {
        match code {
            2 => Some(Self::Uint8),
            4 => Some(Self::Int16),
            8 => Some(Self::Int32),
            16 => Some(Self::Float32),
            64 => Some(Self::Float64),
            256 => Some(Self::Int8),
            512 => Some(Self::Uint16),
            768 => Some(Self::Uint32),
            _ => None,
        }
    }

    /// The raw NIfTI datatype code.
    pub fn code(self) -> i16 {
        match self {
            Self::Uint8 => 2,
            Self::Int16 => 4,
            Self::Int32 => 8,
            Self::Float32 => 16,
            Self::Float64 => 64,
            Self::Int8 => 256,
            Self::Uint16 => 512,
            Self::Uint32 => 768,
        }
    }

    /// Expected `bitpix` for this datatype.
    pub fn bitpix(self) -> i16 {
        match self {
            Self::Uint8 | Self::Int8 => 8,
            Self::Int16 | Self::Uint16 => 16,
            Self::Int32 | Self::Float32 | Self::Uint32 => 32,
            Self::Float64 => 64,
        }
    }

    /// Bytes per voxel.
    pub fn size_bytes(self) -> usize {
        self.bitpix() as usize / 8
    }

    /// Decodes one voxel's raw bytes (already sliced to `size_bytes()`) into
    /// `f64`, honouring endianness. Widening to `f64` here (rather than a
    /// separate integer path per type) keeps `NiftiVolume::values` a single
    /// uniform `Vec<f64>`, matching `tpt-med-dicom`'s own `Vec<i32>` stored-
    /// value convention of "one numeric type downstream code can rely on".
    pub fn decode(self, bytes: &[u8], big_endian: bool) -> f64 {
        match self {
            Self::Uint8 => bytes[0] as f64,
            Self::Int8 => bytes[0] as i8 as f64,
            Self::Int16 => {
                let b = [bytes[0], bytes[1]];
                (if big_endian {
                    i16::from_be_bytes(b)
                } else {
                    i16::from_le_bytes(b)
                }) as f64
            }
            Self::Uint16 => {
                let b = [bytes[0], bytes[1]];
                (if big_endian {
                    u16::from_be_bytes(b)
                } else {
                    u16::from_le_bytes(b)
                }) as f64
            }
            Self::Int32 => {
                let b = [bytes[0], bytes[1], bytes[2], bytes[3]];
                (if big_endian {
                    i32::from_be_bytes(b)
                } else {
                    i32::from_le_bytes(b)
                }) as f64
            }
            Self::Uint32 => {
                let b = [bytes[0], bytes[1], bytes[2], bytes[3]];
                (if big_endian {
                    u32::from_be_bytes(b)
                } else {
                    u32::from_le_bytes(b)
                }) as f64
            }
            Self::Float32 => {
                let b = [bytes[0], bytes[1], bytes[2], bytes[3]];
                (if big_endian {
                    f32::from_be_bytes(b)
                } else {
                    f32::from_le_bytes(b)
                }) as f64
            }
            Self::Float64 => {
                let b = [
                    bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
                ];
                if big_endian {
                    f64::from_be_bytes(b)
                } else {
                    f64::from_le_bytes(b)
                }
            }
        }
    }
}

impl<'a> Header<'a> {
    fn i16_at(&self, offset: usize) -> i16 {
        let b = [self.bytes[offset], self.bytes[offset + 1]];
        if self.big_endian {
            i16::from_be_bytes(b)
        } else {
            i16::from_le_bytes(b)
        }
    }

    fn f32_at(&self, offset: usize) -> f32 {
        let b = [
            self.bytes[offset],
            self.bytes[offset + 1],
            self.bytes[offset + 2],
            self.bytes[offset + 3],
        ];
        if self.big_endian {
            f32::from_be_bytes(b)
        } else {
            f32::from_le_bytes(b)
        }
    }

    /// `dim[0..8]` (byte offset 40).
    pub fn dim(&self) -> [i16; 8] {
        let mut d = [0i16; 8];
        for (i, v) in d.iter_mut().enumerate() {
            *v = self.i16_at(40 + 2 * i);
        }
        d
    }

    /// `datatype` (byte offset 70).
    pub fn datatype_code(&self) -> i16 {
        self.i16_at(70)
    }

    /// `bitpix` (byte offset 72).
    pub fn bitpix(&self) -> i16 {
        self.i16_at(72)
    }

    /// `pixdim[0..8]` (byte offset 76).
    pub fn pixdim(&self) -> [f32; 8] {
        let mut p = [0f32; 8];
        for (i, v) in p.iter_mut().enumerate() {
            *v = self.f32_at(76 + 4 * i);
        }
        p
    }

    /// `vox_offset` (byte offset 108).
    pub fn vox_offset(&self) -> f32 {
        self.f32_at(108)
    }

    /// `scl_slope` (byte offset 112).
    pub fn scl_slope(&self) -> f32 {
        self.f32_at(112)
    }

    /// `scl_inter` (byte offset 116).
    pub fn scl_inter(&self) -> f32 {
        self.f32_at(116)
    }

    /// `qform_code` (byte offset 252).
    pub fn qform_code(&self) -> i16 {
        self.i16_at(252)
    }

    /// `sform_code` (byte offset 254).
    pub fn sform_code(&self) -> i16 {
        self.i16_at(254)
    }

    /// `(quatern_b, quatern_c, quatern_d)` (byte offsets 256, 260, 264).
    pub fn quatern_bcd(&self) -> (f32, f32, f32) {
        (self.f32_at(256), self.f32_at(260), self.f32_at(264))
    }

    /// `(qoffset_x, qoffset_y, qoffset_z)` (byte offsets 268, 272, 276).
    pub fn qoffset(&self) -> (f32, f32, f32) {
        (self.f32_at(268), self.f32_at(272), self.f32_at(276))
    }

    /// `srow_x[0..4]` (byte offset 280).
    pub fn srow_x(&self) -> [f32; 4] {
        [
            self.f32_at(280),
            self.f32_at(284),
            self.f32_at(288),
            self.f32_at(292),
        ]
    }

    /// `srow_y[0..4]` (byte offset 296).
    pub fn srow_y(&self) -> [f32; 4] {
        [
            self.f32_at(296),
            self.f32_at(300),
            self.f32_at(304),
            self.f32_at(308),
        ]
    }

    /// `srow_z[0..4]` (byte offset 312).
    pub fn srow_z(&self) -> [f32; 4] {
        [
            self.f32_at(312),
            self.f32_at(316),
            self.f32_at(320),
            self.f32_at(324),
        ]
    }

    /// Whether the header is big-endian.
    pub fn is_big_endian(&self) -> bool {
        self.big_endian
    }
}

/// NIfTI-2 header size in bytes; also `sizeof_hdr`'s required value.
/// The 64-bit update of NIfTI-1 (540 = 348 struct + padding), with `dim`
/// widened to `i64`, `pixdim`/scaling/quaternions/srows widened to `f64`,
/// `vox_offset` to `i64`, and the magic moved to offset 4 as an 8-byte
/// value (`"n+2\0\r\n\x1a\n"` single-file, `"ni2..."` dual-file).
pub const HEADER2_LEN: usize = 540;

/// A little/big-endian `i16`/`i32`/`i64`/`f32`/`f64` reader over a fixed
/// NIfTI-2 header buffer, with the endianness resolved once by
/// [`parse_header2`] / [`parse_header2_dual`] from `sizeof_hdr`.
pub struct Header2<'a> {
    bytes: &'a [u8],
    big_endian: bool,
}

impl<'a> Header2<'a> {
    fn i16_at(&self, offset: usize) -> i16 {
        let b = [self.bytes[offset], self.bytes[offset + 1]];
        if self.big_endian {
            i16::from_be_bytes(b)
        } else {
            i16::from_le_bytes(b)
        }
    }

    fn i32_at(&self, offset: usize) -> i32 {
        let b = [
            self.bytes[offset],
            self.bytes[offset + 1],
            self.bytes[offset + 2],
            self.bytes[offset + 3],
        ];
        if self.big_endian {
            i32::from_be_bytes(b)
        } else {
            i32::from_le_bytes(b)
        }
    }

    fn i64_at(&self, offset: usize) -> i64 {
        let mut b = [0u8; 8];
        b.copy_from_slice(&self.bytes[offset..offset + 8]);
        if self.big_endian {
            i64::from_be_bytes(b)
        } else {
            i64::from_le_bytes(b)
        }
    }

    fn f64_at(&self, offset: usize) -> f64 {
        let mut b = [0u8; 8];
        b.copy_from_slice(&self.bytes[offset..offset + 8]);
        if self.big_endian {
            f64::from_be_bytes(b)
        } else {
            f64::from_le_bytes(b)
        }
    }

    /// `dim[0..8]` (byte offset 16; `i64` in NIfTI-2).
    pub fn dim(&self) -> [i64; 8] {
        let mut d = [0i64; 8];
        for (i, v) in d.iter_mut().enumerate() {
            *v = self.i64_at(16 + 8 * i);
        }
        d
    }

    /// `datatype` (byte offset 12).
    pub fn datatype_code(&self) -> i16 {
        self.i16_at(12)
    }

    /// `bitpix` (byte offset 14).
    pub fn bitpix(&self) -> i16 {
        self.i16_at(14)
    }

    /// `pixdim[0..8]` (byte offset 104; `f64` in NIfTI-2).
    pub fn pixdim(&self) -> [f64; 8] {
        let mut p = [0f64; 8];
        for (i, v) in p.iter_mut().enumerate() {
            *v = self.f64_at(104 + 8 * i);
        }
        p
    }

    /// `vox_offset` (byte offset 168; `i64` in NIfTI-2).
    pub fn vox_offset(&self) -> i64 {
        self.i64_at(168)
    }

    /// `scl_slope` (byte offset 176).
    pub fn scl_slope(&self) -> f64 {
        self.f64_at(176)
    }

    /// `scl_inter` (byte offset 184).
    pub fn scl_inter(&self) -> f64 {
        self.f64_at(184)
    }

    /// `qform_code` (byte offset 344; `i32` in NIfTI-2).
    pub fn qform_code(&self) -> i32 {
        self.i32_at(344)
    }

    /// `sform_code` (byte offset 348).
    pub fn sform_code(&self) -> i32 {
        self.i32_at(348)
    }

    /// `(quatern_b, quatern_c, quatern_d)` (byte offsets 352, 360, 368).
    pub fn quatern_bcd(&self) -> (f64, f64, f64) {
        (self.f64_at(352), self.f64_at(360), self.f64_at(368))
    }

    /// `(qoffset_x, qoffset_y, qoffset_z)` (byte offsets 376, 384, 392).
    pub fn qoffset(&self) -> (f64, f64, f64) {
        (self.f64_at(376), self.f64_at(384), self.f64_at(392))
    }

    /// `srow_x[0..4]` (byte offset 400; `f64` elements, 8-byte stride).
    pub fn srow_x(&self) -> [f64; 4] {
        [
            self.f64_at(400),
            self.f64_at(408),
            self.f64_at(416),
            self.f64_at(424),
        ]
    }

    /// `srow_y[0..4]` (byte offset 432).
    pub fn srow_y(&self) -> [f64; 4] {
        [
            self.f64_at(432),
            self.f64_at(440),
            self.f64_at(448),
            self.f64_at(456),
        ]
    }

    /// `srow_z[0..4]` (byte offset 464).
    pub fn srow_z(&self) -> [f64; 4] {
        [
            self.f64_at(464),
            self.f64_at(472),
            self.f64_at(480),
            self.f64_at(488),
        ]
    }

    /// Whether the header is big-endian.
    pub fn is_big_endian(&self) -> bool {
        self.big_endian
    }
}

/// Detects endianness from `sizeof_hdr` (byte offset 0) and returns a
/// [`Header2`] view over `bytes`, or an error if neither byte order gives
/// 540, or the magic (byte offset 4, `"n+2"`) is not a single-file NIfTI-2
/// signature.
///
/// `bytes` must be at least [`HEADER2_LEN`] long.
pub fn parse_header2(bytes: &[u8]) -> Result<Header2<'_>> {
    parse_header2_with_magic(bytes, b"n+2")
}

/// Like [`parse_header2`], but requires the dual-file (`.hdr`/`.img`)
/// magic `"ni2"` at byte offset 4 instead of the single-file `"n+2"`.
pub fn parse_header2_dual(bytes: &[u8]) -> Result<Header2<'_>> {
    parse_header2_with_magic(bytes, b"ni2")
}

fn parse_header2_with_magic<'a>(bytes: &'a [u8], magic: &[u8; 3]) -> Result<Header2<'a>> {
    if bytes.len() < HEADER2_LEN {
        return Err(NiftiError::UnexpectedEof {
            offset: bytes.len(),
            while_reading: "NIfTI-2 header",
        });
    }
    let sizeof_hdr_le = i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
    let sizeof_hdr_be = i32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
    let big_endian = if sizeof_hdr_le == HEADER2_LEN as i32 {
        false
    } else if sizeof_hdr_be == HEADER2_LEN as i32 {
        true
    } else {
        return Err(NiftiError::BadValue {
            reason: format!(
                "sizeof_hdr is neither {HEADER2_LEN} (LE) nor {HEADER2_LEN} (BE): \
                 got {sizeof_hdr_le} (LE) / {sizeof_hdr_be} (BE)"
            ),
        });
    };

    // Magic (byte 4): the first three bytes identify the version and
    // layout; bytes 4..12 are the 8-byte `magic[8]` field whose tail is
    // the NIfTI control sequence CR LF SUB LF, which real writers set
    // consistently but which is not worth refusing a file over.
    if &bytes[4..7] != magic {
        return Err(NiftiError::BadValue {
            reason: format!(
                "magic is {:?}, expected {:?}",
                core::str::from_utf8(&bytes[4..7]).unwrap_or("<not utf-8>"),
                core::str::from_utf8(magic).unwrap_or("<not utf-8>")
            ),
        });
    }

    Ok(Header2 { bytes, big_endian })
}

impl<'a> Header<'a> {
    /// Widens this NIfTI-1 header into the shared normalization both
    /// versions' decode path consumes (`f32`/`i16` values are
    /// precision-preservingly widened).
    pub(crate) fn common(&self) -> super::CommonHeader {
        super::CommonHeader {
            dim: self.dim().map(i64::from),
            datatype_code: self.datatype_code(),
            bitpix: self.bitpix(),
            pixdim: self.pixdim().map(f64::from),
            vox_offset: self.vox_offset() as f64,
            scl_slope: self.scl_slope() as f64,
            scl_inter: self.scl_inter() as f64,
            qform_code: i32::from(self.qform_code()),
            sform_code: i32::from(self.sform_code()),
            quatern_bcd: {
                let (b, c, d) = self.quatern_bcd();
                (b as f64, c as f64, d as f64)
            },
            qoffset: {
                let (x, y, z) = self.qoffset();
                (x as f64, y as f64, z as f64)
            },
            srow: [
                self.srow_x().map(f64::from),
                self.srow_y().map(f64::from),
                self.srow_z().map(f64::from),
            ],
            big_endian: self.is_big_endian(),
        }
    }
}

impl<'a> Header2<'a> {
    /// Normalizes this NIfTI-2 header into the shared shape.
    pub(crate) fn common(&self) -> super::CommonHeader {
        super::CommonHeader {
            dim: self.dim(),
            datatype_code: self.datatype_code(),
            bitpix: self.bitpix(),
            pixdim: self.pixdim(),
            vox_offset: self.vox_offset() as f64,
            scl_slope: self.scl_slope(),
            scl_inter: self.scl_inter(),
            qform_code: self.qform_code(),
            sform_code: self.sform_code(),
            quatern_bcd: self.quatern_bcd(),
            qoffset: self.qoffset(),
            srow: [self.srow_x(), self.srow_y(), self.srow_z()],
            big_endian: self.is_big_endian(),
        }
    }
}

/// Detects endianness from `sizeof_hdr` (byte offset 0) and returns a
/// [`Header`] view over `bytes`, or an error if neither byte order gives
/// 348, or the magic (byte offset 344, `"n+1"`) is not a single-file NIfTI-1
/// signature.
///
/// `bytes` must be at least [`HEADER_LEN`] long.
pub fn parse_header(bytes: &[u8]) -> Result<Header<'_>> {
    parse_header_with_magic(bytes, "n+1")
}

/// Like [`parse_header`], but requires the dual-file (`.hdr`/`.img`) magic
/// `"ni1"` at byte offset 344 instead of the single-file `"n+1"`. Which
/// magic is *correct* depends on which entry point the caller used —
/// `NiftiVolume::parse_bytes` takes `n+1`, `NiftiVolume::parse_dual_bytes`
/// takes `ni1` — so the check lives here with the magic as a parameter and
/// each entry point asks for its own.
///
/// `bytes` must be at least [`HEADER_LEN`] long.
pub fn parse_header_dual(bytes: &[u8]) -> Result<Header<'_>> {
    parse_header_with_magic(bytes, "ni1")
}

fn parse_header_with_magic<'a>(bytes: &'a [u8], magic: &str) -> Result<Header<'a>> {
    if bytes.len() < HEADER_LEN {
        return Err(NiftiError::UnexpectedEof {
            offset: bytes.len(),
            while_reading: "NIfTI-1 header",
        });
    }
    let sizeof_hdr_le = i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
    let sizeof_hdr_be = i32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
    let big_endian = if sizeof_hdr_le == HEADER_LEN as i32 {
        false
    } else if sizeof_hdr_be == HEADER_LEN as i32 {
        true
    } else {
        return Err(NiftiError::BadValue {
            reason: format!(
                "sizeof_hdr is neither {HEADER_LEN} (LE) nor {HEADER_LEN} (BE): \
                 got {sizeof_hdr_le} (LE) / {sizeof_hdr_be} (BE)"
            ),
        });
    };

    // Magic (byte 344): only the first 3 bytes are checked -- some
    // real-world writers pad the 4th byte inconsistently (a known quirk,
    // not something worth refusing a file over).
    if &bytes[344..347] != magic.as_bytes() {
        return Err(NiftiError::BadValue {
            reason: format!(
                "magic is {:?}, expected {magic:?}",
                core::str::from_utf8(&bytes[344..347]).unwrap_or("<not utf-8>")
            ),
        });
    }

    Ok(Header { bytes, big_endian })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn datatype_bitpix_table_is_internally_consistent() {
        for &dt in &[
            Datatype::Uint8,
            Datatype::Int16,
            Datatype::Int32,
            Datatype::Float32,
            Datatype::Float64,
            Datatype::Int8,
            Datatype::Uint16,
            Datatype::Uint32,
        ] {
            assert_eq!(Datatype::from_code(dt.code()), Some(dt));
            assert_eq!(dt.bitpix() as usize, dt.size_bytes() * 8);
        }
    }

    #[test]
    fn unknown_datatype_code_is_none() {
        assert_eq!(Datatype::from_code(32), None); // DT_COMPLEX64, unsupported
        assert_eq!(Datatype::from_code(2304), None); // DT_RGBA32, unsupported
    }
}
