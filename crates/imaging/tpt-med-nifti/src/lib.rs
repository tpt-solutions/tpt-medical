//! Pure-Rust NIfTI-1 volume parsing for research-space imaging.
//!
//! Reads an uncompressed, single-file NIfTI-1 (`.nii`) volume into
//! [`NiftiVolume`]: voxel dimensions, physical spacing, RAS orientation, and
//! per-voxel scaled values. It does not decode acquisition metadata (no
//! `descrip`/`aux_file`/`intent_*`/slice-timing fields), does not decompress
//! `.nii.gz`, and does not read the dual-file `.hdr`/`.img` form — see
//! `rfcs/0006-nifti-ingestion.md` for the full scope and why.
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
//! synthetic writer (`write_to_bytes`) through the parser, plus one
//! hand-built minimal header verified field-by-field against values computed
//! by hand from the bytes — not against this crate's own writer. There is no
//! golden dataset entry (no simulation output originates here) and no test
//! against a real-world writer's output (`dcm2niix`, FSL, …) — see the RFC's
//! "What remains unverified".

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
mod geometry;
mod header;
pub mod synthetic;

use std::path::Path;

pub use error::{NiftiError, Result};
pub use header::Datatype;
use header::{parse_header, HEADER_LEN};
use tpt_med_geometry::{Mat3, Vec3};

/// A parsed NIfTI-1 volume: dimensions, physical geometry (RAS), and scaled
/// per-voxel values.
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
}

impl NiftiVolume {
    /// Parses a NIfTI-1 volume from a file on disk.
    pub fn parse_file(path: &Path) -> Result<Self> {
        let bytes = std::fs::read(path)?;
        Self::parse_bytes(&bytes).map_err(|e| match e {
            NiftiError::NotNifti(_) => NiftiError::NotNifti(path.to_path_buf()),
            NiftiError::Gzipped(_) => NiftiError::Gzipped(path.to_path_buf()),
            other => other,
        })
    }

    /// Parses a NIfTI-1 volume from an in-memory byte buffer.
    pub fn parse_bytes(bytes: &[u8]) -> Result<Self> {
        if bytes.len() >= 2 && bytes[0] == 0x1F && bytes[1] == 0x8B {
            return Err(NiftiError::Gzipped(std::path::PathBuf::from("<memory>")));
        }
        if bytes.len() < HEADER_LEN {
            return Err(NiftiError::NotNifti(std::path::PathBuf::from("<memory>")));
        }
        let header = match parse_header(bytes) {
            Ok(h) => h,
            Err(NiftiError::BadValue { .. }) => {
                return Err(NiftiError::NotNifti(std::path::PathBuf::from("<memory>")))
            }
            Err(e) => return Err(e),
        };

        let dim = header.dim();
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

        let datatype_code = header.datatype_code();
        let datatype = Datatype::from_code(datatype_code)
            .ok_or(NiftiError::UnsupportedDatatype(datatype_code))?;
        let bitpix = header.bitpix();
        if bitpix != datatype.bitpix() {
            return Err(NiftiError::BitpixMismatch {
                datatype: datatype_code,
                bitpix,
            });
        }

        let pixdim = header.pixdim();
        let qfac = if pixdim[0] == 0.0 {
            1.0
        } else {
            pixdim[0] as f64
        };
        let (spacing, origin, rotation) = if header.sform_code() > 0 {
            geometry_from_sform(&header)?
        } else if header.qform_code() > 0 {
            geometry_from_qform(&header, qfac)
        } else {
            (
                (pixdim[1] as f64, pixdim[2] as f64, pixdim[3] as f64),
                Vec3::ZERO,
                Mat3::IDENTITY,
            )
        };

        let vox_offset = header.vox_offset();
        if !vox_offset.is_finite() || vox_offset < HEADER_LEN as f32 {
            return Err(NiftiError::BadValue {
                reason: format!(
                    "vox_offset {vox_offset} is before the end of the {HEADER_LEN}-byte header"
                ),
            });
        }
        let data_start = vox_offset as usize;

        let n_voxels = nx * ny * nz;
        let sample_bytes = datatype.size_bytes();
        let data_end =
            data_start
                .checked_add(n_voxels * sample_bytes)
                .ok_or(NiftiError::BadValue {
                    reason: "voxel data extent overflows usize".into(),
                })?;
        if bytes.len() < data_end {
            return Err(NiftiError::UnexpectedEof {
                offset: bytes.len(),
                while_reading: "voxel data",
            });
        }

        let big_endian = header.is_big_endian();
        let scl_slope = header.scl_slope() as f64;
        let scl_inter = header.scl_inter() as f64;
        let mut values = Vec::with_capacity(n_voxels);
        for i in 0..n_voxels {
            let start = data_start + i * sample_bytes;
            let raw = datatype.decode(&bytes[start..start + sample_bytes], big_endian);
            let value = if scl_slope != 0.0 {
                raw * scl_slope + scl_inter
            } else {
                raw
            };
            values.push(value);
        }

        Ok(Self {
            dims: (nx, ny, nz),
            voxel_spacing: spacing,
            origin,
            rotation,
            datatype,
            values,
        })
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

/// Builds `(spacing, origin, rotation)` from the sform affine (`srow_x/y/z`):
/// each row directly gives `ras_axis = row[0..3]·(i,j,k) + row[3]`. The
/// affine mixes rotation and spacing together (unlike DICOM's separate
/// direction-cosine and spacing tags), so spacing is recovered as each
/// column's norm and the rotation as that column normalised — an explicit
/// no-shear assumption, see RFC 0006 Drawbacks.
fn geometry_from_sform(header: &header::Header<'_>) -> Result<((f64, f64, f64), Vec3, Mat3)> {
    let sx = header.srow_x();
    let sy = header.srow_y();
    let sz = header.srow_z();
    let origin = Vec3::new(sx[3] as f64, sy[3] as f64, sz[3] as f64);

    let cols = [
        Vec3::new(sx[0] as f64, sy[0] as f64, sz[0] as f64),
        Vec3::new(sx[1] as f64, sy[1] as f64, sz[1] as f64),
        Vec3::new(sx[2] as f64, sy[2] as f64, sz[2] as f64),
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
fn geometry_from_qform(header: &header::Header<'_>, qfac: f64) -> ((f64, f64, f64), Vec3, Mat3) {
    let (b, c, d) = header.quatern_bcd();
    let rotation = geometry::quaternion_to_rotation(b as f64, c as f64, d as f64);
    let (ox, oy, oz) = header.qoffset();
    let pixdim = header.pixdim();

    // qfac flips the k (third) column's handedness; magnitude spacing stays
    // positive, matching pixdim's own always-positive convention.
    let k_col = if qfac < 0.0 {
        rotation.col(2) * -1.0
    } else {
        rotation.col(2)
    };
    let rotation = Mat3::from_cols(rotation.col(0), rotation.col(1), k_col);

    (
        (pixdim[1] as f64, pixdim[2] as f64, pixdim[3] as f64),
        Vec3::new(ox as f64, oy as f64, oz as f64),
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

    #[test]
    fn gzip_magic_is_rejected_with_named_error() {
        let err = NiftiVolume::parse_bytes(&[0x1F, 0x8B, 0x08, 0x00]).unwrap_err();
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
}
