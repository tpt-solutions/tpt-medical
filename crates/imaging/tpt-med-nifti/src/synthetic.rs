//! A minimal, spec-correct NIfTI-1 writer for tests: round-tripping a known
//! volume through this writer and [`crate::NiftiVolume::parse_bytes`]
//! exercises the real byte layout, the same role
//! `tpt_med_dicom::SyntheticCtBuilder` plays for the DICOM parser. Not a
//! general-purpose NIfTI encoder — no dtype other than the two variants
//! tests need, no PHI (there is none to have).

use tpt_med_geometry::Vec3;

use crate::header::{HEADER2_LEN, HEADER_LEN};

/// Byte offset voxel data starts at when no header extensions are written:
/// the 348-byte header plus the mandatory 4-byte extension-flag field.
const VOX_OFFSET: usize = HEADER_LEN + 4;

enum Values {
    U8(Vec<f64>),
    F32(Vec<f64>),
}

impl Values {
    /// `(datatype, bitpix)` for this payload.
    fn bits(&self) -> (i16, i16) {
        match self {
            Values::U8(_) => (2, 8),
            Values::F32(_) => (16, 32),
        }
    }

    /// Appends the raw little-endian voxel samples to `buf`.
    fn append_to(&self, buf: &mut Vec<u8>) {
        match self {
            Values::U8(v) => {
                for x in v {
                    buf.push(*x as u8);
                }
            }
            Values::F32(v) => {
                for x in v {
                    buf.extend_from_slice(&(*x as f32).to_le_bytes());
                }
            }
        }
    }
}

/// Qform parameters, when the builder is asked to write a qform instead of
/// (the default) sform.
struct Qform {
    quatern_bcd: (f32, f32, f32),
    qfac: f32,
}

/// Builds a minimal single-file NIfTI-1 (`.nii`) byte buffer.
pub struct SyntheticNiftiBuilder {
    dims: (usize, usize, usize),
    spacing: (f64, f64, f64),
    origin: Vec3,
    scl_slope: f32,
    scl_inter: f32,
    qform: Option<Qform>,
    values: Option<Values>,
}

impl SyntheticNiftiBuilder {
    /// Starts a builder for a `dims` volume with axis-aligned spacing
    /// `spacing` (mm) and no scaling (`scl_slope = 1`).
    pub fn new(nx: usize, ny: usize, nz: usize, spacing: (f64, f64, f64)) -> Self {
        Self {
            dims: (nx, ny, nz),
            spacing,
            origin: Vec3::ZERO,
            scl_slope: 1.0,
            scl_inter: 0.0,
            qform: None,
            values: None,
        }
    }

    /// Sets the origin (RAS mm) — the sform translation, or the qform
    /// `qoffset` if [`Self::with_qform`] is used instead.
    pub fn with_origin(mut self, origin: Vec3) -> Self {
        self.origin = origin;
        self
    }

    /// Writes a qform instead of the default sform: `quatern_bcd` are the
    /// stored quaternion components, `qfac` is `pixdim[0]` (must be `1.0` or
    /// `-1.0`). The origin set via [`Self::with_origin`] becomes `qoffset`.
    pub fn with_qform(mut self, quatern_bcd: (f32, f32, f32), qfac: f32) -> Self {
        self.qform = Some(Qform { quatern_bcd, qfac });
        self
    }

    /// Sets `scl_slope`/`scl_inter`.
    pub fn with_scaling(mut self, slope: f32, inter: f32) -> Self {
        self.scl_slope = slope;
        self.scl_inter = inter;
        self
    }

    /// Writes `values` as `DT_UINT8`. Each value is truncated to `u8` — the
    /// caller is responsible for supplying values already in range.
    pub fn with_values_u8(mut self, values: &[f64]) -> Self {
        self.values = Some(Values::U8(values.to_vec()));
        self
    }

    /// Writes `values` as `DT_FLOAT32`.
    pub fn with_values_f32(mut self, values: &[f64]) -> Self {
        self.values = Some(Values::F32(values.to_vec()));
        self
    }

    /// Builds a single-file `.nii` byte buffer: the 348-byte header (magic
    /// `n+1`), the mandatory 4-byte extension flag, then the voxel data at
    /// `vox_offset`. Panics if no values were set (a synthetic fixture with
    /// no data is a test-writer bug, not a runtime condition).
    pub fn build(mut self) -> Vec<u8> {
        let values = self
            .values
            .take()
            .expect("SyntheticNiftiBuilder needs values");
        let (datatype, bitpix) = values.bits();
        let mut buf = self.write_header(b"n+1", VOX_OFFSET as f32, datatype, bitpix);
        buf.extend_from_slice(&[0u8; 4]); // mandatory extension-flag field
        debug_assert_eq!(buf.len(), VOX_OFFSET);

        let (nx, ny, nz) = self.dims;
        let expected_len = VOX_OFFSET + nx * ny * nz * (bitpix as usize / 8);
        values.append_to(&mut buf);
        debug_assert_eq!(buf.len(), expected_len);
        buf
    }

    /// Builds a dual-file (`.hdr`/`.img`) pair: a 348-byte `.hdr` carrying
    /// the magic `ni1` and `vox_offset = 0` (the spec's two-file convention
    /// — the voxel data starts at byte 0 of the `.img`), and an `.img`
    /// holding nothing but the voxel data. The halves parse through
    /// `NiftiVolume::parse_dual_bytes` exactly as `build`'s output parses
    /// through `NiftiVolume::parse_bytes`.
    pub fn build_dual(mut self) -> (Vec<u8>, Vec<u8>) {
        let values = self
            .values
            .take()
            .expect("SyntheticNiftiBuilder needs values");
        let (datatype, bitpix) = values.bits();
        let hdr = self.write_header(b"ni1", 0.0, datatype, bitpix);
        let mut img = Vec::new();
        values.append_to(&mut img);
        (hdr, img)
    }

    /// Builds a NIfTI-2 single-file (`.nii`, magic `n+2`) buffer: the
    /// 540-byte `nifti_2_header` (i64 dims, f64 geometry and scaling,
    /// i64 `vox_offset` — see `header::HEADER2_LEN`), then the extension
    /// flag and the voxel data. Parses through
    /// `NiftiVolume::parse_bytes` with `version == NiftiVersion::V2`.
    pub fn build2(mut self) -> Vec<u8> {
        let values = self
            .values
            .take()
            .expect("SyntheticNiftiBuilder needs values");
        let (datatype, bitpix) = values.bits();
        // "n+2" NUL CR LF SUB LF — the 8-byte NIfTI-2 magic.
        let magic: [u8; 8] = [b'n', b'+', b'2', 0, 0x0D, 0x0A, 0x1A, 0x0A];
        let mut buf = self.write_header2(&magic, VOX2_OFFSET, datatype, bitpix);
        buf.extend_from_slice(&[0u8; 4]); // mandatory extension-flag field
        debug_assert_eq!(buf.len(), VOX2_OFFSET);

        let (nx, ny, nz) = self.dims;
        let expected_len = VOX2_OFFSET + nx * ny * nz * (bitpix as usize / 8);
        values.append_to(&mut buf);
        debug_assert_eq!(buf.len(), expected_len);
        buf
    }

    /// Builds a NIfTI-2 dual-file (`.hdr`/`.img`, magic `ni2`) pair, the
    /// NIfTI-2 counterpart of [`Self::build_dual`]: `vox_offset = 0` into
    /// the `.img`.
    pub fn build2_dual(mut self) -> (Vec<u8>, Vec<u8>) {
        let values = self
            .values
            .take()
            .expect("SyntheticNiftiBuilder needs values");
        let (datatype, bitpix) = values.bits();
        let magic: [u8; 8] = [b'n', b'i', b'2', 0, 0x0D, 0x0A, 0x1A, 0x0A];
        let hdr = self.write_header2(&magic, 0, datatype, bitpix);
        let mut img = Vec::new();
        values.append_to(&mut img);
        (hdr, img)
    }

    /// The fixed 540-byte NIfTI-2 header: same semantics as
    /// [`Self::write_header`]'s NIfTI-1 layout at the `nifti_2_header`
    /// offsets and widths.
    fn write_header2(
        &self,
        magic: &[u8; 8],
        vox_offset: usize,
        datatype: i16,
        bitpix: i16,
    ) -> Vec<u8> {
        let mut buf = vec![0u8; HEADER2_LEN];
        write_i32(&mut buf, 0, HEADER2_LEN as i32); // sizeof_hdr
        buf[4..12].copy_from_slice(magic);

        write_i16(&mut buf, 12, datatype);
        write_i16(&mut buf, 14, bitpix);
        // dim[0..8] as i64.
        write_i64(&mut buf, 16, 3);
        write_i64(&mut buf, 24, self.dims.0 as i64);
        write_i64(&mut buf, 32, self.dims.1 as i64);
        write_i64(&mut buf, 40, self.dims.2 as i64);

        // pixdim[0..8] as f64.
        write_f64(
            &mut buf,
            104,
            self.qform.as_ref().map_or(1.0, |q| q.qfac as f64),
        );
        write_f64(&mut buf, 112, self.spacing.0);
        write_f64(&mut buf, 120, self.spacing.1);
        write_f64(&mut buf, 128, self.spacing.2);

        write_i64(&mut buf, 168, vox_offset as i64);
        write_f64(&mut buf, 176, self.scl_slope as f64);
        write_f64(&mut buf, 184, self.scl_inter as f64);

        match &self.qform {
            Some(q) => {
                write_i32(&mut buf, 344, 1); // qform_code = SCANNER_ANAT
                write_f64(&mut buf, 352, q.quatern_bcd.0 as f64);
                write_f64(&mut buf, 360, q.quatern_bcd.1 as f64);
                write_f64(&mut buf, 368, q.quatern_bcd.2 as f64);
                write_f64(&mut buf, 376, self.origin.x);
                write_f64(&mut buf, 384, self.origin.y);
                write_f64(&mut buf, 392, self.origin.z);
            }
            None => {
                // sform: identity rotation, spacing on the diagonal.
                write_i32(&mut buf, 348, 1); // sform_code = SCANNER_ANAT
                write_f64(&mut buf, 400, self.spacing.0); // srow_x[0]
                write_f64(&mut buf, 424, self.origin.x); // srow_x[3]
                write_f64(&mut buf, 440, self.spacing.1); // srow_y[1]
                write_f64(&mut buf, 456, self.origin.y); // srow_y[3]
                write_f64(&mut buf, 480, self.spacing.2); // srow_z[2]
                write_f64(&mut buf, 488, self.origin.z); // srow_z[3]
            }
        }
        buf
    }

    /// The fixed 348-byte header with `magic` at offset 344 and
    /// `vox_offset` at offset 108; dimensions, geometry, datatype and
    /// scaling come from the builder's own fields.
    fn write_header(
        &self,
        magic: &[u8; 3],
        vox_offset: f32,
        datatype: i16,
        bitpix: i16,
    ) -> Vec<u8> {
        let mut buf = vec![0u8; HEADER_LEN];
        write_i32(&mut buf, 0, HEADER_LEN as i32); // sizeof_hdr

        // dim[0..8]: dim[0]=3, dim[1..4]=nx,ny,nz, rest unused (0).
        write_i16(&mut buf, 40, 3);
        write_i16(&mut buf, 42, self.dims.0 as i16);
        write_i16(&mut buf, 44, self.dims.1 as i16);
        write_i16(&mut buf, 46, self.dims.2 as i16);

        write_i16(&mut buf, 70, datatype);
        write_i16(&mut buf, 72, bitpix);

        // pixdim[1..4]=spacing; pixdim[0]=qfac (only meaningful for qform).
        write_f32(&mut buf, 80, self.spacing.0 as f32);
        write_f32(&mut buf, 84, self.spacing.1 as f32);
        write_f32(&mut buf, 88, self.spacing.2 as f32);

        write_f32(&mut buf, 108, vox_offset);
        write_f32(&mut buf, 112, self.scl_slope);
        write_f32(&mut buf, 116, self.scl_inter);

        match &self.qform {
            Some(q) => {
                write_f32(&mut buf, 76, q.qfac);
                write_i16(&mut buf, 252, 1); // qform_code = NIFTI_XFORM_SCANNER_ANAT
                write_f32(&mut buf, 256, q.quatern_bcd.0);
                write_f32(&mut buf, 260, q.quatern_bcd.1);
                write_f32(&mut buf, 264, q.quatern_bcd.2);
                write_f32(&mut buf, 268, self.origin.x as f32);
                write_f32(&mut buf, 272, self.origin.y as f32);
                write_f32(&mut buf, 276, self.origin.z as f32);
            }
            None => {
                write_f32(&mut buf, 76, 1.0);
                // sform: identity rotation, so srow_i has spacing.i on the
                // diagonal.
                write_i16(&mut buf, 254, 1); // sform_code = NIFTI_XFORM_SCANNER_ANAT
                write_f32(&mut buf, 280, self.spacing.0 as f32); // srow_x
                write_f32(&mut buf, 284, 0.0);
                write_f32(&mut buf, 288, 0.0);
                write_f32(&mut buf, 292, self.origin.x as f32);
                write_f32(&mut buf, 296, 0.0); // srow_y
                write_f32(&mut buf, 300, self.spacing.1 as f32);
                write_f32(&mut buf, 304, 0.0);
                write_f32(&mut buf, 308, self.origin.y as f32);
                write_f32(&mut buf, 312, 0.0); // srow_z
                write_f32(&mut buf, 316, 0.0);
                write_f32(&mut buf, 320, self.spacing.2 as f32);
                write_f32(&mut buf, 324, self.origin.z as f32);
            }
        }

        buf[344..347].copy_from_slice(magic);
        buf[347] = 0; // 4th magic byte padded; the parser checks only 3
        buf
    }
}

/// Byte offset of the first voxel in a NIfTI-2 single file: the
/// 540-byte header plus the 4-byte extension flag.
const VOX2_OFFSET: usize = crate::header::HEADER2_LEN + 4;

fn write_i64(buf: &mut [u8], offset: usize, v: i64) {
    buf[offset..offset + 8].copy_from_slice(&v.to_le_bytes());
}

fn write_f64(buf: &mut [u8], offset: usize, v: f64) {
    buf[offset..offset + 8].copy_from_slice(&v.to_le_bytes());
}

fn write_i16(buf: &mut [u8], offset: usize, v: i16) {
    buf[offset..offset + 2].copy_from_slice(&v.to_le_bytes());
}

fn write_i32(buf: &mut [u8], offset: usize, v: i32) {
    buf[offset..offset + 4].copy_from_slice(&v.to_le_bytes());
}

fn write_f32(buf: &mut [u8], offset: usize, v: f32) {
    buf[offset..offset + 4].copy_from_slice(&v.to_le_bytes());
}
