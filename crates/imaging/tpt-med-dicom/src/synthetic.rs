//! Synthetic CT series generation for tests and demos.
//!
//! Produces valid DICOM Part-10 files (explicit VR little endian) from a
//! per-voxel HU function. All repository test data is generated here — real
//! patient data must never be committed.

use std::io::Write;
use std::path::Path;

use crate::parser::encode_element_explicit;
use crate::tags::{
    Tag, Vr, BITS_ALLOCATED, COLUMNS, IMAGE_ORIENTATION_PATIENT, IMAGE_POSITION_PATIENT,
    INSTANCE_NUMBER, MODALITY, PATIENT_ID, PIXEL_DATA, PIXEL_REPRESENTATION, PIXEL_SPACING,
    RESCALE_INTERCEPT, RESCALE_SLOPE, ROWS, SERIES_INSTANCE_UID, SLICE_THICKNESS, STUDY_DATE,
    TRANSFER_SYNTAX_UID,
};

const EXPLICIT_VR_LE: &str = "1.2.840.10008.1.2.1";

/// Builder for a synthetic axial CT series.
pub struct SyntheticCtBuilder {
    cols: usize,
    rows: usize,
    nslices: usize,
    spacing: (f64, f64),
    slice_thickness: f64,
    origin: [f64; 3],
    patient_id: String,
    study_date: String,
    hu_fn: Box<dyn Fn(usize, usize, usize) -> f64>,
}

impl Default for SyntheticCtBuilder {
    fn default() -> Self {
        Self {
            cols: 16,
            rows: 16,
            nslices: 8,
            spacing: (1.0, 1.0),
            slice_thickness: 1.0,
            origin: [-8.0, -8.0, 0.0],
            patient_id: "SYNTHETIC".to_string(),
            study_date: "20000101".to_string(),
            hu_fn: Box::new(|_, _, _| -1000.0),
        }
    }
}

impl SyntheticCtBuilder {
    /// Builder with 16×16×8 defaults at 1 mm isotropic spacing.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets column count.
    pub fn cols(mut self, n: usize) -> Self {
        self.cols = n;
        self
    }

    /// Sets row count.
    pub fn rows(mut self, n: usize) -> Self {
        self.rows = n;
        self
    }

    /// Sets slice count.
    pub fn slices(mut self, n: usize) -> Self {
        self.nslices = n;
        self
    }

    /// Sets in-plane spacing (row, col) in mm.
    pub fn spacing(mut self, row: f64, col: f64) -> Self {
        self.spacing = (row, col);
        self
    }

    /// Sets slice thickness in mm (also used as slice pitch).
    pub fn slice_thickness(mut self, t: f64) -> Self {
        self.slice_thickness = t;
        self
    }

    /// Sets the patient-space position of voxel (0,0,0).
    pub fn origin(mut self, origin: [f64; 3]) -> Self {
        self.origin = origin;
        self
    }

    /// Sets the PatientID tag value (synthetic only!).
    pub fn patient_id(mut self, id: &str) -> Self {
        self.patient_id = id.to_string();
        self
    }

    /// Sets the HU function of voxel indices `(x=col, y=row, z=slice)`.
    pub fn hu_fn(mut self, f: impl Fn(usize, usize, usize) -> f64 + 'static) -> Self {
        self.hu_fn = Box::new(f);
        self
    }

    /// Builds the in-memory series.
    pub fn build(self, series_uid: &str) -> SyntheticCtSeries {
        let mut slices = Vec::with_capacity(self.nslices);
        for z in 0..self.nslices {
            let mut hu = Vec::with_capacity(self.cols * self.rows);
            for y in 0..self.rows {
                for x in 0..self.cols {
                    hu.push((self.hu_fn)(x, y, z));
                }
            }
            let stored: Vec<i16> = hu
                .iter()
                .map(|&h| h.clamp(-32768.0, 32767.0).round() as i16)
                .collect();
            let bytes = Self::encode_slice(
                &stored,
                self.cols,
                self.rows,
                z as u32 + 1,
                self.spacing,
                self.slice_thickness,
                self.origin,
                series_uid,
                &self.patient_id,
                &self.study_date,
            );
            slices.push(SyntheticSlice {
                file_name: format!("slice_{:04}.dcm", z + 1),
                bytes,
            });
        }
        SyntheticCtSeries {
            series_uid: series_uid.to_string(),
            slices,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn encode_slice(
        pixels: &[i16],
        cols: usize,
        rows: usize,
        instance: u32,
        spacing: (f64, f64),
        thickness: f64,
        origin: [f64; 3],
        series_uid: &str,
        patient_id: &str,
        study_date: &str,
    ) -> Vec<u8> {
        let mut buf: Vec<u8> = vec![0u8; 128];
        buf.extend_from_slice(crate::DICM_MAGIC);

        // Meta group (always explicit VR LE).
        buf.extend_from_slice(&encode_element_explicit(
            TRANSFER_SYNTAX_UID,
            Vr::Ui,
            format!("{EXPLICIT_VR_LE}\0").as_bytes(),
        ));

        let ds = |t: Tag, vr: Vr, val: &[u8]| encode_element_explicit(t, vr, val);
        // DICOM values must have even length; pad per VR padding rules.
        let ds_str = |t: Tag, vr: Vr, s: &str| {
            let pad = if vr == Vr::Ui { '\0' } else { ' ' };
            let mut owned = s.to_string();
            if owned.len() % 2 == 1 {
                owned.push(pad);
            }
            encode_element_explicit(t, vr, owned.as_bytes())
        };
        let ds_us = |t: Tag, v: u16| encode_element_explicit(t, Vr::Us, &v.to_le_bytes());
        let ds_f64s = |t: Tag, vals: &[f64]| {
            let mut s = vals
                .iter()
                .map(|v| format!("{v:.8}"))
                .collect::<Vec<_>>()
                .join("\\");
            if s.len() % 2 == 1 {
                s.push(' ');
            }
            encode_element_explicit(t, Vr::Ds, s.as_bytes())
        };

        buf.extend_from_slice(&ds_str(STUDY_DATE, Vr::Da, study_date));
        buf.extend_from_slice(&ds_str(MODALITY, Vr::Cs, "CT"));
        buf.extend_from_slice(&ds_str(PATIENT_ID, Vr::Lo, patient_id));
        buf.extend_from_slice(&ds_f64s(SLICE_THICKNESS, &[thickness]));
        buf.extend_from_slice(&ds_str(
            SERIES_INSTANCE_UID,
            Vr::Ui,
            &format!("{series_uid}\0"),
        ));
        buf.extend_from_slice(&ds(
            INSTANCE_NUMBER,
            Vr::Is,
            format!("{instance} ").as_bytes(),
        ));
        let pos = format!(
            "{:.8}\\{:.8}\\{:.8} ",
            origin[0],
            origin[1],
            origin[2] + instance as f64 * thickness
        );
        buf.extend_from_slice(&ds_str(IMAGE_POSITION_PATIENT, Vr::Ds, &pos));
        buf.extend_from_slice(&ds_str(
            IMAGE_ORIENTATION_PATIENT,
            Vr::Ds,
            "1\\0\\0\\0\\1\\0 ", // axial: rows +x, cols +y
        ));
        buf.extend_from_slice(&ds_us(ROWS, rows as u16));
        buf.extend_from_slice(&ds_us(COLUMNS, cols as u16));
        buf.extend_from_slice(&ds_f64s(PIXEL_SPACING, &[spacing.0, spacing.1]));
        buf.extend_from_slice(&ds_us(BITS_ALLOCATED, 16));
        buf.extend_from_slice(&ds_us(PIXEL_REPRESENTATION, 1)); // signed
        buf.extend_from_slice(&ds_f64s(RESCALE_INTERCEPT, &[0.0]));
        buf.extend_from_slice(&ds_f64s(RESCALE_SLOPE, &[1.0]));
        let mut pixel_bytes = Vec::with_capacity(pixels.len() * 2);
        for p in pixels {
            pixel_bytes.extend_from_slice(&p.to_le_bytes());
        }
        buf.extend_from_slice(&ds(PIXEL_DATA, Vr::Ow, &pixel_bytes));
        buf
    }
}

/// One encoded synthetic slice.
#[derive(Debug, Clone)]
pub struct SyntheticSlice {
    /// File name to use when writing to disk.
    pub file_name: String,
    /// Complete DICOM Part-10 file bytes.
    pub bytes: Vec<u8>,
}

/// An in-memory synthetic series.
#[derive(Debug, Clone)]
pub struct SyntheticCtSeries {
    /// Series UID used for all slices.
    pub series_uid: String,
    /// Encoded slices in ascending instance order.
    pub slices: Vec<SyntheticSlice>,
}

impl SyntheticCtSeries {
    /// Parses the encoded slices back into a [`DicomSeries`].
    pub fn parse(&self) -> crate::Result<crate::DicomSeries> {
        let slices = self
            .slices
            .iter()
            .map(|s| crate::DicomParser::parse_bytes(&s.bytes))
            .collect::<crate::Result<Vec<_>>>()?;
        crate::DicomSeries::from_slices(slices)
    }

    /// Writes the series into a directory (created if missing), together
    /// with a `manifest.json` listing the slice files (used by the browser
    /// demos to fetch an exact file list without probing for 404s).
    pub fn write_to_dir(&self, dir: &Path) -> std::io::Result<()> {
        std::fs::create_dir_all(dir)?;
        for slice in &self.slices {
            let mut f = std::fs::File::create(dir.join(&slice.file_name))?;
            f.write_all(&slice.bytes)?;
        }
        let names = self
            .slices
            .iter()
            .map(|s| format!("    \"{}\"", s.file_name))
            .collect::<Vec<_>>()
            .join(",\n");
        let manifest = format!(
            "{{\n  \"series_uid\": \"{}\",\n  \"slices\": [\n{}\n  ]\n}}\n",
            self.series_uid, names
        );
        std::fs::write(dir.join("manifest.json"), manifest)?;
        Ok(())
    }
}

/// Generates the canonical synthetic proximal-femur-like phantom used by the
/// CLI milestone example: a cortical shell (700 HU) around a trabecular
/// core (150 HU) on a soft-tissue background (50 HU), all in air (−1000).
pub fn femur_phantom(cols: usize, rows: usize, nslices: usize) -> SyntheticCtSeries {
    let cx = (cols as f64 - 1.0) / 2.0;
    let cy = (rows as f64 - 1.0) / 2.0;
    let cz = (nslices as f64 - 1.0) / 2.0;
    let r_outer = (cols.min(rows) as f64) * 0.40;
    let r_inner = r_outer * 0.62;

    SyntheticCtBuilder::new()
        .cols(cols)
        .rows(rows)
        .slices(nslices)
        .spacing(1.0, 1.0)
        .slice_thickness(1.0)
        .origin([
            -(cols as f64 - 1.0) / 2.0,
            -(rows as f64 - 1.0) / 2.0,
            -(nslices as f64 - 1.0) / 2.0,
        ])
        .patient_id("SYNTHETIC-FEMUR")
        .hu_fn(move |x, y, z| {
            let dx = x as f64 - cx;
            let dy = y as f64 - cy;
            let dz = z as f64 - cz;
            // Slightly tapered cylinder along z
            let taper = 1.0 - 0.25 * (dz / cz.abs().max(1.0)).abs();
            let r = (dx * dx + dy * dy).sqrt();
            if r <= r_inner * taper {
                150.0 // trabecular core
            } else if r <= r_outer * taper {
                700.0 // cortical shell
            } else {
                50.0 // soft tissue background
            }
        })
        .build("1.2.826.0.1.3680043.8.498.9001.1")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DicomParser;

    #[test]
    fn phantom_parses_with_expected_hu() {
        let series = femur_phantom(16, 16, 8);
        assert_eq!(series.slices.len(), 8);
        let parsed = DicomParser::parse_bytes(&series.slices[0].bytes).expect("parse");
        assert_eq!((parsed.columns, parsed.rows), (16, 16));
        assert_eq!(parsed.pixel_data.len(), 256);
        // Centre voxel is trabecular core (150 HU); corner is soft tissue.
        let center = parsed.hu_at(8, 8).expect("center");
        assert!((center - 150.0).abs() < 1e-6, "got {center}");
        let corner = parsed.hu_at(0, 0).expect("corner");
        assert!((corner - 50.0).abs() < 1e-6, "got {corner}");
        assert_eq!(parsed.modality.as_deref(), Some("CT"));
        assert_eq!(parsed.instance_number, 1);
    }

    #[test]
    fn positions_increase_along_z() {
        let series = femur_phantom(8, 8, 4);
        let mut zs = Vec::new();
        for s in &series.slices {
            let parsed = DicomParser::parse_bytes(&s.bytes).expect("parse");
            zs.push(parsed.position.z);
        }
        assert!(zs.windows(2).all(|w| w[1] > w[0]), "{zs:?}");
    }

    #[test]
    fn write_and_reload_roundtrip() {
        let series = femur_phantom(8, 8, 3);
        let dir = std::env::temp_dir().join("tpt-med-synthetic-test");
        series.write_to_dir(&dir).expect("write");
        let reloaded = crate::DicomSeries::load_from_dir(&dir).expect("reload");
        assert_eq!(reloaded.slices.len(), 3);
        assert_eq!(reloaded.modality, "CT");
        std::fs::remove_dir_all(&dir).ok();
    }
}
