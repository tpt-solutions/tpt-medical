//! CT series assembly: individual slices, sorted volumes, and HU access.

use std::path::Path;

use crate::error::{DicomError, Result};
use crate::parser::{DicomElement, DicomParser};
use crate::tags::{
    TransferSyntax, BITS_ALLOCATED, COLUMNS, IMAGE_ORIENTATION_PATIENT, IMAGE_POSITION_PATIENT,
    INSTANCE_NUMBER, MODALITY, PATIENT_ID, PIXEL_DATA, PIXEL_REPRESENTATION, PIXEL_SPACING,
    RESCALE_INTERCEPT, RESCALE_SLOPE, ROWS, SERIES_INSTANCE_UID, SLICE_THICKNESS, STUDY_DATE,
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
    /// The dataset's transfer syntax, needed to decode encapsulated pixel data.
    /// `None` only in `Default`; the parser always sets it.
    transfer_syntax: Option<TransferSyntax>,
}

impl SliceBuilder {
    pub(crate) fn new(transfer_syntax: TransferSyntax) -> Self {
        Self {
            transfer_syntax: Some(transfer_syntax),
            ..Self::default()
        }
    }

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
                self.slice.pixel_data = self.decode_pixel_data(&el)?;
            }
            _ => {}
        }
        Ok(())
    }

    pub(crate) fn build(self) -> Result<DicomSlice> {
        Ok(self.slice)
    }

    /// Decodes a PixelData element according to the dataset's transfer syntax.
    ///
    /// For a compressed syntax the value holds the concatenated encapsulated
    /// fragments collected by the parser, so it is handed to the matching
    /// decoder. When that decoder is behind a cargo feature that is not
    /// enabled, the error names the feature rather than reporting a
    /// decode failure: the data is fine, the build just cannot read it.
    fn decode_pixel_data(&self, el: &DicomElement) -> Result<Vec<i32>> {
        let expected = self.slice.rows as usize * self.slice.columns as usize;

        macro_rules! feature_gated_decode {
            ($feature:literal, $module:ident, $uid:literal, $name:literal) => {{
                #[cfg(feature = $feature)]
                {
                    crate::$module::decode_frame(
                        &el.value,
                        self.slice.rows,
                        self.slice.columns,
                        self.bits_allocated,
                        self.pixel_representation,
                        el.tag,
                    )
                }
                #[cfg(not(feature = $feature))]
                {
                    let _ = (el, expected);
                    Err(DicomError::CompressedPixelData(
                        concat!($uid, " (", $name, ")").to_string(),
                    ))
                }
            }};
        }

        match self.transfer_syntax {
            Some(TransferSyntax::RleLossless) => {
                feature_gated_decode!("rle", rle, "1.2.840.10008.1.2.5", "RLE Lossless")
            }
            Some(TransferSyntax::JpegBaseline) => {
                feature_gated_decode!("jpeg", jpeg, "1.2.840.10008.1.2.4.50", "JPEG Baseline")
            }
            Some(TransferSyntax::JpegExtended) => {
                feature_gated_decode!("jpeg", jpeg, "1.2.840.10008.1.2.4.51", "JPEG Extended")
            }
            Some(TransferSyntax::JpegLossless) => {
                feature_gated_decode!("jpeg", jpeg, "1.2.840.10008.1.2.4.57", "JPEG Lossless")
            }
            Some(TransferSyntax::JpegLosslessSv1) => {
                feature_gated_decode!("jpeg", jpeg, "1.2.840.10008.1.2.4.70", "JPEG Lossless, SV1")
            }
            Some(TransferSyntax::JpegLsLossless) => {
                feature_gated_decode!(
                    "jpeg-ls",
                    jpeg_ls,
                    "1.2.840.10008.1.2.4.80",
                    "JPEG-LS Lossless"
                )
            }
            Some(TransferSyntax::JpegLsNearLossless) => {
                feature_gated_decode!(
                    "jpeg-ls",
                    jpeg_ls,
                    "1.2.840.10008.1.2.4.81",
                    "JPEG-LS Near-Lossless"
                )
            }
            Some(TransferSyntax::Jpeg2000Lossless) => {
                feature_gated_decode!(
                    "jpeg2000",
                    jpeg2000,
                    "1.2.840.10008.1.2.4.90",
                    "JPEG 2000 Lossless"
                )
            }
            Some(TransferSyntax::Jpeg2000) => {
                feature_gated_decode!("jpeg2000", jpeg2000, "1.2.840.10008.1.2.4.91", "JPEG 2000")
            }
            Some(TransferSyntax::Jpeg2000Part2MultiComponentLossless) => {
                feature_gated_decode!(
                    "jpeg2000",
                    jpeg2000,
                    "1.2.840.10008.1.2.4.92",
                    "JPEG 2000 Part 2 Multi-component Lossless"
                )
            }
            Some(TransferSyntax::Jpeg2000Part2MultiComponent) => {
                feature_gated_decode!(
                    "jpeg2000",
                    jpeg2000,
                    "1.2.840.10008.1.2.4.93",
                    "JPEG 2000 Part 2 Multi-component"
                )
            }
            _ => decode_pixels(
                &el.value,
                self.bits_allocated,
                self.pixel_representation,
                expected,
                el.tag,
            ),
        }
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

#[cfg(all(test, feature = "jpeg-ls"))]
mod encapsulated_pixel_data_tests {
    use crate::parser::{encode_element_explicit, DicomParser};
    use crate::tags::{
        Vr, BITS_ALLOCATED, COLUMNS, PIXEL_DATA, PIXEL_REPRESENTATION, ROWS, TRANSFER_SYNTAX_UID,
    };

    fn explicit_us(tag: crate::tags::Tag, v: u16) -> Vec<u8> {
        encode_element_explicit(tag, Vr::Us, &v.to_le_bytes())
    }

    /// Builds an encapsulated PixelData element (PS3.5 Annex A.4): tag, `OB`
    /// long-form header with undefined length, an empty Basic Offset Table
    /// item (single frame), one fragment item holding the whole compressed
    /// frame, then the sequence delimiter.
    fn encapsulated_pixel_data(frame: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&PIXEL_DATA.0.to_le_bytes());
        out.extend_from_slice(&PIXEL_DATA.1.to_le_bytes());
        out.extend_from_slice(b"OB");
        out.extend_from_slice(&[0, 0]); // reserved
        out.extend_from_slice(&0xFFFF_FFFFu32.to_le_bytes()); // undefined length
                                                              // Basic Offset Table: empty (single frame).
        out.extend_from_slice(&0xFFFEu16.to_le_bytes());
        out.extend_from_slice(&0xE000u16.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        // Fragment.
        out.extend_from_slice(&0xFFFEu16.to_le_bytes());
        out.extend_from_slice(&0xE000u16.to_le_bytes());
        out.extend_from_slice(&(frame.len() as u32).to_le_bytes());
        out.extend_from_slice(frame);
        // Sequence delimiter.
        out.extend_from_slice(&0xFFFEu16.to_le_bytes());
        out.extend_from_slice(&0xE0DDu16.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        out
    }

    fn part10(transfer_syntax_uid: &[u8], dataset: &[u8]) -> Vec<u8> {
        let mut buf = vec![0u8; 128];
        buf.extend_from_slice(b"DICM");
        buf.extend_from_slice(&encode_element_explicit(
            TRANSFER_SYNTAX_UID,
            Vr::Ui,
            transfer_syntax_uid,
        ));
        buf.extend_from_slice(dataset);
        buf
    }

    /// The parser and `SliceBuilder` are exercised together here, unlike
    /// `jpeg_ls::decode_frame`'s own unit tests: this proves the transfer
    /// syntax is actually recognised from its UID, the encapsulated
    /// fragment stream is actually collected by the parser (not just handed
    /// to the decoder pre-assembled), and the decoded values reach
    /// `DicomSlice::pixel_data` end to end.
    #[test]
    fn jpeg_ls_lossless_end_to_end() {
        let pixels: Vec<u16> = vec![100, 4095, 0, 2048];
        let (w, h) = (2u32, 2u32);
        let mut frame = Vec::new();
        jpegls::encode(&pixels, w, h, &mut frame).expect("encodes");

        let mut dataset = Vec::new();
        dataset.extend_from_slice(&explicit_us(ROWS, h as u16));
        dataset.extend_from_slice(&explicit_us(COLUMNS, w as u16));
        dataset.extend_from_slice(&explicit_us(BITS_ALLOCATED, 16));
        dataset.extend_from_slice(&explicit_us(PIXEL_REPRESENTATION, 0));
        dataset.extend_from_slice(&encapsulated_pixel_data(&frame));

        let buf = part10(b"1.2.840.10008.1.2.4.80 ", &dataset);
        let slice = DicomParser::parse_bytes(&buf).expect("parses");
        assert_eq!(slice.rows, 2);
        assert_eq!(slice.columns, 2);
        assert_eq!(
            slice.pixel_data,
            pixels.iter().map(|&v| v as i32).collect::<Vec<_>>()
        );
    }
}

#[cfg(all(test, feature = "jpeg2000"))]
mod jpeg2000_encapsulated_pixel_data_tests {
    use crate::error::DicomError;
    use crate::jpeg2000::{MINIMAL_J2C_2X2, SSIZ0_OFFSET};
    use crate::parser::{encode_element_explicit, DicomParser};
    use crate::tags::{
        Vr, BITS_ALLOCATED, COLUMNS, PIXEL_DATA, PIXEL_REPRESENTATION, ROWS, TRANSFER_SYNTAX_UID,
    };
    use crate::DicomSlice;

    fn explicit_us(tag: crate::tags::Tag, v: u16) -> Vec<u8> {
        encode_element_explicit(tag, Vr::Us, &v.to_le_bytes())
    }

    /// Encapsulated PixelData element, PS3.5 Annex A.4: `OB` undefined-length
    /// header, empty Basic Offset Table (single frame), one fragment item, then
    /// the sequence delimiter. Mirrors the same helper in the sibling
    /// `jpeg-ls` test module, which documents the layout in full.
    fn encapsulated_pixel_data(frame: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&PIXEL_DATA.0.to_le_bytes());
        out.extend_from_slice(&PIXEL_DATA.1.to_le_bytes());
        out.extend_from_slice(b"OB");
        out.extend_from_slice(&[0, 0]); // reserved
        out.extend_from_slice(&0xFFFF_FFFFu32.to_le_bytes()); // undefined length
                                                              // Basic Offset Table: empty (single frame).
        out.extend_from_slice(&0xFFFEu16.to_le_bytes());
        out.extend_from_slice(&0xE000u16.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        // Fragment.
        out.extend_from_slice(&0xFFFEu16.to_le_bytes());
        out.extend_from_slice(&0xE000u16.to_le_bytes());
        out.extend_from_slice(&(frame.len() as u32).to_le_bytes());
        out.extend_from_slice(frame);
        // Sequence delimiter.
        out.extend_from_slice(&0xFFFEu16.to_le_bytes());
        out.extend_from_slice(&0xE0DDu16.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        out
    }

    fn part10(transfer_syntax_uid: &[u8], dataset: &[u8]) -> Vec<u8> {
        let mut buf = vec![0u8; 128];
        buf.extend_from_slice(b"DICM");
        buf.extend_from_slice(&encode_element_explicit(
            TRANSFER_SYNTAX_UID,
            Vr::Ui,
            transfer_syntax_uid,
        ));
        buf.extend_from_slice(dataset);
        buf
    }

    /// A Part 10 byte stream whose PixelData holds `frame`, declared under
    /// `transfer_syntax_uid`. `bits_allocated`/`pixel_representation` are
    /// written to match the codestream, so a mismatch is never what makes a
    /// case below fail.
    fn stream_with(
        transfer_syntax_uid: &str,
        frame: &[u8],
        bits_allocated: u16,
        pixel_representation: u16,
    ) -> DicomSlice {
        let mut dataset = Vec::new();
        dataset.extend_from_slice(&explicit_us(ROWS, 2));
        dataset.extend_from_slice(&explicit_us(COLUMNS, 2));
        dataset.extend_from_slice(&explicit_us(BITS_ALLOCATED, bits_allocated));
        dataset.extend_from_slice(&explicit_us(PIXEL_REPRESENTATION, pixel_representation));
        dataset.extend_from_slice(&encapsulated_pixel_data(frame));

        // UI values are space-padded to an even length by PS3.5.
        let uid = format!("{transfer_syntax_uid} ");
        let buf = part10(uid.as_bytes(), &dataset);
        DicomParser::parse_bytes(&buf).expect("parses")
    }

    /// The JPEG 2000 counterpart of the `jpeg-ls` module's
    /// `jpeg_ls_lossless_end_to_end`, run once per transfer syntax this crate
    /// routes to `jpeg2000::decode_frame`. As with that test, the point is the
    /// whole path: the transfer syntax is recognised from its UID, the parser
    /// collects the encapsulated fragment stream itself (not handed a
    /// pre-assembled frame), and the decoded samples reach
    /// `DicomSlice::pixel_data`.
    ///
    /// The shared fixture is a 2x2 8-bit unsigned tile with all coefficients
    /// zero, so every sample decodes to the unsigned DC level-shift midpoint,
    /// `2^(8-1) = 128` — see `jpeg2000::MINIMAL_J2C_2X2`.
    #[test]
    fn every_jpeg2000_transfer_syntax_decodes_end_to_end() {
        for uid in [
            "1.2.840.10008.1.2.4.90", // JPEG 2000 Lossless Only
            "1.2.840.10008.1.2.4.91", // JPEG 2000
            "1.2.840.10008.1.2.4.92", // JPEG 2000 Part 2 Multi-component Lossless Only
            "1.2.840.10008.1.2.4.93", // JPEG 2000 Part 2 Multi-component
        ] {
            let slice = stream_with(uid, MINIMAL_J2C_2X2, 8, 0);
            assert_eq!(slice.rows, 2, "{uid}");
            assert_eq!(slice.columns, 2, "{uid}");
            assert_eq!(slice.pixel_data, vec![128, 128, 128, 128], "{uid}");
        }
    }

    /// `.92`/`.93` carry "multi-component" in their *transfer syntax* name,
    /// not a promise about this crate's decoding: `decode_frame` is still
    /// single-component, and a `SamplesPerPixel = 1` Part 2 file (a Part 2
    /// file using Part 2 extensions for some other reason, e.g. an
    /// alternative wavelet kernel) decodes normally. This pins that, so the
    /// routing cannot later start rejecting `.92`/`.93` outright.
    #[test]
    fn part2_single_component_files_decode() {
        for uid in ["1.2.840.10008.1.2.4.92", "1.2.840.10008.1.2.4.93"] {
            let slice = stream_with(uid, MINIMAL_J2C_2X2, 8, 0);
            assert_eq!(slice.pixel_data, vec![128, 128, 128, 128], "{uid}");
        }
    }

    /// The signed-bit workaround (`jpeg2000::component0_signed` re-reading
    /// `Ssiz` and undoing the codec's unconditional unsigned level shift) has
    /// to survive the full parse path too, not just `decode_frame`'s unit
    /// tests: the codestream declares signed *and* the dataset says
    /// `PixelRepresentation = 1`, so the level shift is undone and the
    /// zero-coefficient tile reads as 0 rather than 128.
    #[test]
    fn signed_pixel_representation_survives_the_full_parse_path() {
        let mut stream = MINIMAL_J2C_2X2.to_vec();
        assert_eq!(stream[SSIZ0_OFFSET], 0x07, "fixture layout changed");
        stream[SSIZ0_OFFSET] = 0x87; // Ssiz bit 7: component 0 is signed

        let slice = stream_with("1.2.840.10008.1.2.4.90", &stream, 8, 1);
        assert_eq!(slice.pixel_data, vec![0, 0, 0, 0]);
    }

    /// A non-conformant file whose codestream signedness disagrees with the
    /// dataset's `PixelRepresentation` is rejected at parse time, not
    /// silently decoded to a shifted image. Confirms the rejection is a
    /// `DicomError` rather than a panic on the parse path.
    #[test]
    fn disagreeing_signedness_is_rejected_at_parse_time() {
        let err = DicomParser::parse_bytes(&{
            let mut dataset = Vec::new();
            dataset.extend_from_slice(&explicit_us(ROWS, 2));
            dataset.extend_from_slice(&explicit_us(COLUMNS, 2));
            dataset.extend_from_slice(&explicit_us(BITS_ALLOCATED, 8));
            dataset.extend_from_slice(&explicit_us(PIXEL_REPRESENTATION, 1));
            dataset.extend_from_slice(&encapsulated_pixel_data(MINIMAL_J2C_2X2));
            part10(b"1.2.840.10008.1.2.4.90 ", &dataset)
        })
        .expect_err("unsigned codestream with signed PixelRepresentation must be rejected");
        assert!(
            matches!(err, DicomError::BadValue { .. }),
            "unexpected error: {err:?}"
        );
    }
}
