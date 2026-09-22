//! CT series assembly: individual slices, sorted volumes, and HU access.

use std::path::Path;

use crate::error::{DicomError, Result};
use crate::parser::{DicomElement, DicomParser};
use crate::tags::{
    BITS_ALLOCATED, COLUMNS, IMAGE_ORIENTATION_PATIENT, IMAGE_POSITION_PATIENT, INSTANCE_NUMBER,
    MODALITY, PATIENT_ID, PIXEL_DATA, PIXEL_REPRESENTATION, PIXEL_SPACING, RESCALE_INTERCEPT,
    RESCALE_SLOPE, ROWS, SERIES_INSTANCE_UID, SLICE_THICKNESS, STUDY_DATE,
};
use tpt_med_geometry::{ImageFrame, Vec3};

/// One parsed CT/MR slice with geometry and raw stored pixel values.
#[derive(Debug, Clone)]
pub struct DicomSlice {
    /// Instance Number (0020,0013).
    pub instance_number: u32,
    /// Image Position (Patient): position of the first transmitted voxel.
    pub position: Vec3,
    /// Image Orientation (Patient): `[row_dir, col_dir]` direction cosines.
    pub orientation: [[f64; 3]; 2],
    /// Rows.
    pub rows: u16,
    /// Columns.
    pub columns: u16,
    /// Raw stored pixel values (`stored × slope + intercept = HU`).
    pub pixel_data: Vec<i32>,
    /// Rescale slope.
    pub rescale_slope: f64,
    /// Rescale intercept.
    pub rescale_intercept: f64,
    /// In-plane pixel spacing `(row_spacing, col_spacing)` in mm.
    pub pixel_spacing: (f64, f64),
    /// Slice thickness in mm.
    pub slice_thickness: f64,
    /// Modality, if present.
    pub modality: Option<String>,
    /// Patient ID tag, if present. Treat as PHI: never log.
    pub patient_id: Option<String>,
    /// Study date, if present. Semi-identifying: never log.
    pub study_date: Option<String>,
    /// Raw series UID, if present.
    pub series_uid: Option<String>,
}

impl Default for DicomSlice {
    fn default() -> Self {
        Self {
            instance_number: 0,
            position: Vec3::ZERO,
            orientation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            rows: 0,
            columns: 0,
            pixel_data: Vec::new(),
            rescale_slope: 1.0,
            rescale_intercept: 0.0,
            pixel_spacing: (1.0, 1.0),
            slice_thickness: 1.0,
            modality: None,
            patient_id: None,
            study_date: None,
            series_uid: None,
        }
    }
}

impl DicomSlice {
    /// Hounsfield Unit value at (row, col).
    pub fn hu_at(&self, row: usize, col: usize) -> Option<f64> {
        let idx = row * self.columns as usize + col;
        self.pixel_data
            .get(idx)
            .map(|&v| v as f64 * self.rescale_slope + self.rescale_intercept)
    }

    /// Full HU plane, row-major.
    pub fn hu_plane(&self) -> Vec<f64> {
        self.pixel_data
            .iter()
            .map(|&v| v as f64 * self.rescale_slope + self.rescale_intercept)
            .collect()
    }

    /// Slice normal (row × column direction cosines).
    pub fn normal(&self) -> Vec3 {
        let row = Vec3::new(
            self.orientation[0][0],
            self.orientation[0][1],
            self.orientation[0][2],
        );
        let col = Vec3::new(
            self.orientation[1][0],
            self.orientation[1][1],
            self.orientation[1][2],
        );
        row.cross(col).normalize()
    }

    /// Image frame geometry for this slice.
    pub fn frame(&self) -> ImageFrame {
        ImageFrame {
            origin: self.position,
            row_dir: Vec3::new(
                self.orientation[0][0],
                self.orientation[0][1],
                self.orientation[0][2],
            ),
            col_dir: Vec3::new(
                self.orientation[1][0],
                self.orientation[1][1],
                self.orientation[1][2],
            ),
        }
    }
}

/// Incrementally assembles a [`DicomSlice`] from parsed elements.
#[derive(Debug, Default)]
pub(crate) struct SliceBuilder {
    slice: DicomSlice,
    bits_allocated: u16,
    pixel_representation: u16,
}

impl SliceBuilder {
    pub(crate) fn absorb(&mut self, el: DicomElement) -> Result<()> {
        match el.tag {
            INSTANCE_NUMBER => self.slice.instance_number = el.as_is()? as u32,
            IMAGE_POSITION_PATIENT => {
                let v = el.as_ds_vec()?;
                if v.len() >= 3 {
                    self.slice.position = Vec3::new(v[0], v[1], v[2]);
                }
            }
            IMAGE_ORIENTATION_PATIENT => {
                let v = el.as_ds_vec()?;
                if v.len() >= 6 {
                    self.slice.orientation = [[v[0], v[1], v[2]], [v[3], v[4], v[5]]];
                }
            }
            ROWS => self.slice.rows = el.as_us()?,
            COLUMNS => self.slice.columns = el.as_us()?,
            BITS_ALLOCATED => self.bits_allocated = el.as_us()?,
            PIXEL_REPRESENTATION => self.pixel_representation = el.as_us()?,
            RESCALE_SLOPE => self.slice.rescale_slope = el.as_ds_first()?,
            RESCALE_INTERCEPT => self.slice.rescale_intercept = el.as_ds_first()?,
            PIXEL_SPACING => {
                let v = el.as_ds_vec()?;
                if v.len() >= 2 {
                    self.slice.pixel_spacing = (v[0], v[1]);
                }
            }
            SLICE_THICKNESS => self.slice.slice_thickness = el.as_ds_first()?,
            MODALITY => self.slice.modality = Some(el.as_text()),
            PATIENT_ID => self.slice.patient_id = Some(el.as_text()),
            STUDY_DATE => self.slice.study_date = Some(el.as_text()),
            SERIES_INSTANCE_UID => self.slice.series_uid = Some(el.as_text()),
            PIXEL_DATA => {
                // Defined-length native pixel data only; encapsulated data
                // arrives as an empty value after undefined-length skipping.
                self.slice.pixel_data = decode_pixels(
                    &el.value,
                    self.bits_allocated,
                    self.pixel_representation,
                    self.slice.rows as usize * self.slice.columns as usize,
                    el.tag,
                )?;
            }
            _ => {}
        }
        Ok(())
    }

    pub(crate) fn build(self) -> Result<DicomSlice> {
        Ok(self.slice)
    }
}

fn decode_pixels(
    bytes: &[u8],
    bits: u16,
    signed: u16,
    expected_len: usize,
    tag: crate::tags::Tag,
) -> Result<Vec<i32>> {
    match (bits, signed) {
        (16, 1) => {
            if bytes.len() < expected_len * 2 {
                return Ok(Vec::new()); // truncated/absent pixel payload
            }
            Ok(bytes[..expected_len * 2]
                .chunks_exact(2)
                .map(|c| i16::from_le_bytes([c[0], c[1]]) as i32)
                .collect())
        }
        (16, 0) => {
            if bytes.len() < expected_len * 2 {
                return Ok(Vec::new());
            }
            Ok(bytes[..expected_len * 2]
                .chunks_exact(2)
                .map(|c| u16::from_le_bytes([c[0], c[1]]) as i32)
                .collect())
        }
        (8, s) => {
            let expected = if expected_len == 0 {
                bytes.len()
            } else {
                expected_len
            };
            if bytes.len() < expected {
                return Ok(Vec::new());
            }
            Ok(bytes[..expected]
                .iter()
                .map(|&b| if s == 1 { b as i8 as i32 } else { b as i32 })
                .collect())
        }
        _ => Err(DicomError::BadValue {
            tag,
            reason: format!("unsupported pixel layout bits={bits} signed={signed}"),
        }),
    }
}

/// An assembled CT series: sorted slices plus shared geometry.
#[derive(Debug, Clone)]
pub struct DicomSeries {
    /// Series Instance UID.
    pub series_uid: String,
    /// Slices ordered along the series normal (ascending projection).
    pub slices: Vec<DicomSlice>,
    /// In-plane pixel spacing `(row_spacing, col_spacing)` in mm.
    pub pixel_spacing: (f64, f64),
    /// Slice thickness in mm.
    pub slice_thickness: f64,
    /// Modality of the series.
    pub modality: String,
}

impl DicomSeries {
    /// Volume dimensions: `(cols, rows, slices)`.
    pub fn dims(&self) -> (usize, usize, usize) {
        (
            self.slices.first().map_or(0, |s| s.columns as usize),
            self.slices.first().map_or(0, |s| s.rows as usize),
            self.slices.len(),
        )
    }

    /// HU volume, index order `(slice, row, col)` flattened row-major:
    /// `idx = (z * rows + y) * cols + x`.
    pub fn hu_volume(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.slices.iter().map(|s| s.pixel_data.len()).sum());
        for slice in &self.slices {
            out.extend(slice.hu_plane());
        }
        out
    }

    /// HU value at voxel `(x=col, y=row, z=slice_index_in_sorted_order)`.
    pub fn hu_at(&self, x: usize, y: usize, z: usize) -> Option<f64> {
        self.slices.get(z)?.hu_at(y, x)
    }

    /// Loads every DICOM Part-10 file in a directory and assembles the
    /// (single) series they belong to.
    pub fn load_from_dir(dir: &Path) -> Result<Self> {
        let mut entries: Vec<_> = std::fs::read_dir(dir)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.is_file())
            .collect();
        entries.sort();
        if entries.is_empty() {
            return Err(DicomError::InconsistentSeries(format!(
                "no files in {}",
                dir.display()
            )));
        }

        let mut slices = Vec::new();
        for path in entries {
            // Files without the Part-10 magic (READMEs etc.) are skipped so
            // documentation can live alongside data.
            match std::fs::read(&path) {
                Ok(bytes) if bytes.len() >= 132 && &bytes[128..132] == crate::DICM_MAGIC => {
                    slices.push(DicomParser::parse_bytes(&bytes).map_err(|e| {
                        DicomError::InconsistentSeries(format!("{}: {e}", path.display()))
                    })?);
                }
                _ => continue,
            }
        }
        if slices.is_empty() {
            return Err(DicomError::InconsistentSeries(format!(
                "no DICOM files in {}",
                dir.display()
            )));
        }
        Self::from_slices(slices)
    }

    /// Assembles a series from parsed slices, sorting them along their
    /// common normal (cross product of row/column cosines).
    pub fn from_slices(mut slices: Vec<DicomSlice>) -> Result<Self> {
        if slices.is_empty() {
            return Err(DicomError::InconsistentSeries("no slices".into()));
        }
        let normal = slices[0].normal();
        slices.sort_by(|a, b| {
            a.position
                .dot(normal)
                .partial_cmp(&b.position.dot(normal))
                .unwrap_or(core::cmp::Ordering::Equal)
        });

        let dims_consistent = slices
            .iter()
            .all(|s| s.rows == slices[0].rows && s.columns == slices[0].columns);
        if !dims_consistent {
            return Err(DicomError::InconsistentSeries(
                "slices have differing rows/columns".into(),
            ));
        }

        let series_uid = slices
            .iter()
            .find_map(|s| s.series_uid.clone())
            .unwrap_or_default();
        let modality = slices
            .iter()
            .find_map(|s| s.modality.clone())
            .unwrap_or_else(|| "CT".to_string());
        let pixel_spacing = slices[0].pixel_spacing;
        let slice_thickness = slices[0].slice_thickness;

        Ok(Self {
            series_uid,
            slices,
            pixel_spacing,
            slice_thickness,
            modality,
        })
    }
}
