//! CT series assembly: individual slices, sorted volumes, and HU access.

use std::path::Path;

use crate::error::{DicomError, Result};
use crate::parser::{DicomElement, DicomParser};
use crate::tags::{
    TransferSyntax, BITS_ALLOCATED, COLUMNS, IMAGE_ORIENTATION_PATIENT, IMAGE_POSITION_PATIENT,
    INSTANCE_NUMBER, MODALITY, NUMBER_OF_FRAMES, PATIENT_ID, PER_FRAME_FUNCTIONAL_GROUPS,
    PIXEL_DATA, PIXEL_REPRESENTATION, PIXEL_SPACING, PIXEL_VALUE_TRANSFORMATION_SEQUENCE,
    PLANE_ORIENTATION_SEQUENCE, PLANE_POSITION_SEQUENCE, RESCALE_INTERCEPT, RESCALE_SLOPE, ROWS,
    SERIES_INSTANCE_UID, SHARED_FUNCTIONAL_GROUPS, SLICE_THICKNESS, STUDY_DATE,
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

/// Incrementally assembles [`DicomSlice`]s from parsed elements — one
/// slice for a single-frame object, one per frame for a multi-frame
/// object (PS3.3 C.7.6.6).
#[derive(Debug, Default)]
pub(crate) struct SliceBuilder {
    slice: DicomSlice,
    bits_allocated: u16,
    pixel_representation: u16,
    /// The dataset's transfer syntax, needed to decode encapsulated pixel data.
    /// `None` only in `Default`; the parser always sets it.
    transfer_syntax: Option<TransferSyntax>,
    /// Number of Frames (0028,0008); 0 = tag absent = single-frame.
    n_frames: u32,
    /// The per-frame functional-group items, in frame order.
    per_frame_items: Vec<Vec<DicomElement>>,
    /// The raw PixelData element value — native little-endian samples for
    /// the uncompressed syntaxes, concatenated fragments for the
    /// encapsulated ones. Decoded once, at `build` time, because the
    /// frame count (which sizes the decode) may be tabulated after the
    /// PixelData position would otherwise have decoded it.
    pixel_value: Vec<u8>,
    /// Whether a PixelData element was present at all; without it the
    /// slice keeps empty pixels (a geometry-only parse).
    pixel_seen: bool,
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
            NUMBER_OF_FRAMES => self.n_frames = el.as_is()? as u32,
            SHARED_FUNCTIONAL_GROUPS => {
                // Re-absorb the shared macros' leaves through the very same
                // match arms as top-level tags, so a shared Image Orientation
                // and a top-level one cannot be decoded differently.
                for leaf in functional_group_leaves(&el.items) {
                    self.absorb(leaf)?;
                }
            }
            PER_FRAME_FUNCTIONAL_GROUPS => self.per_frame_items = el.items,
            PIXEL_DATA => {
                self.pixel_value = el.value;
                self.pixel_seen = true;
            }
            _ => {}
        }
        Ok(())
    }

    pub(crate) fn build(mut self) -> Result<Vec<DicomSlice>> {
        let frames = self.n_frames.max(1) as usize;
        let samples_per_frame = self.slice.rows as usize * self.slice.columns as usize;
        let pixels = if self.pixel_seen {
            self.decode_pixel_data(samples_per_frame * frames)?
        } else {
            Vec::new()
        };

        if frames == 1 {
            self.slice.pixel_data = pixels;
            return Ok(vec![self.slice]);
        }

        // Multi-frame: the frames are cut out of the one native PixelData
        // payload and given their per-frame functional-group geometry. An
        // encapsulated multi-frame payload decodes only as one image here,
        // so it is rejected rather than silently presenting frame 0 N times.
        if self
            .transfer_syntax
            .is_some_and(TransferSyntax::is_encapsulated)
        {
            return Err(DicomError::BadValue {
                tag: PIXEL_DATA,
                reason: "multi-frame encapsulated pixel data is not supported; \
                         use an uncompressed transfer syntax"
                    .into(),
            });
        }
        if self.per_frame_items.len() != frames {
            return Err(DicomError::BadValue {
                tag: PER_FRAME_FUNCTIONAL_GROUPS,
                reason: format!(
                    "NumberOfFrames = {} but {} per-frame functional-group items present",
                    frames,
                    self.per_frame_items.len()
                ),
            });
        }

        let mut out = Vec::with_capacity(frames);
        for (i, item) in self.per_frame_items.iter().enumerate() {
            let mut s = self.slice.clone();
            // The frame index stands in for the (per-frame meaningless)
            // instance number; ordering is positional either way.
            s.instance_number = i as u32;
            // Plane Position is the one geometry a frame cannot inherit:
            // it is what distinguishes the frames (PS3.3 C.7.6.16 2).
            let pos = seq_leaf(item, PLANE_POSITION_SEQUENCE, IMAGE_POSITION_PATIENT).ok_or_else(
                || DicomError::BadValue {
                    tag: PLANE_POSITION_SEQUENCE,
                    reason: format!("frame {i} carries no Plane Position Sequence"),
                },
            )?;
            let v = pos.as_ds_vec()?;
            if v.len() >= 3 {
                s.position = Vec3::new(v[0], v[1], v[2]);
            }
            if let Some(o) = seq_leaf(item, PLANE_ORIENTATION_SEQUENCE, IMAGE_ORIENTATION_PATIENT) {
                let v = o.as_ds_vec()?;
                if v.len() >= 6 {
                    s.orientation = [[v[0], v[1], v[2]], [v[3], v[4], v[5]]];
                }
            }
            if let Some(slope) = seq_leaf(item, PIXEL_VALUE_TRANSFORMATION_SEQUENCE, RESCALE_SLOPE)
            {
                s.rescale_slope = slope.as_ds_first()?;
            }
            if let Some(intercept) =
                seq_leaf(item, PIXEL_VALUE_TRANSFORMATION_SEQUENCE, RESCALE_INTERCEPT)
            {
                s.rescale_intercept = intercept.as_ds_first()?;
            }
            s.pixel_data = pixels[i * samples_per_frame..(i + 1) * samples_per_frame].to_vec();
            out.push(s);
        }
        Ok(out)
    }

    /// Decodes the retained PixelData payload according to the dataset's
    /// transfer syntax, expecting `expected` stored samples in total.
    ///
    /// For a compressed syntax the payload holds the concatenated
    /// encapsulated fragments collected by the parser, so it is handed to
    /// the matching decoder (single-frame only — a multi-frame compressed
    /// object is rejected in [`Self::build`] before this is reached). When
    /// that decoder is behind a cargo feature that is not enabled, the
    /// error names the feature rather than reporting a decode failure: the
    /// data is fine, the build just cannot read it.
    fn decode_pixel_data(&self, expected: usize) -> Result<Vec<i32>> {
        macro_rules! feature_gated_decode {
            ($feature:literal, $module:ident, $uid:literal, $name:literal) => {{
                #[cfg(feature = $feature)]
                {
                    crate::$module::decode_frame(
                        &self.pixel_value,
                        self.slice.rows,
                        self.slice.columns,
                        self.bits_allocated,
                        self.pixel_representation,
                        PIXEL_DATA,
                    )
                }
                #[cfg(not(feature = $feature))]
                {
                    let _ = expected;
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
                &self.pixel_value,
                self.bits_allocated,
                self.pixel_representation,
                expected,
                PIXEL_DATA,
            ),
        }
    }
}

/// The leaf elements carried by a functional-groups sequence's items —
/// each macro sequence (`Plane Orientation Sequence`, ...) holds one item
/// of top-level-shaped attributes. Shared items come first, so a later
/// top-level tag still wins, matching how the single-frame reader treats
/// repeated tags.
fn functional_group_leaves(items: &[Vec<DicomElement>]) -> Vec<DicomElement> {
    let mut out = Vec::new();
    for item in items {
        for el in item {
            for sub_item in &el.items {
                out.extend(sub_item.iter().cloned());
            }
        }
    }
    out
}

/// The `leaf`-tagged element inside `item`'s `seq`-tagged sequence's
/// first item — the shape every functional-group macro this crate reads
/// shares.
fn seq_leaf(
    item: &[DicomElement],
    seq: crate::tags::Tag,
    leaf: crate::tags::Tag,
) -> Option<&DicomElement> {
    item.iter()
        .find(|el| el.tag == seq)?
        .items
        .first()?
        .iter()
        .find(|el| el.tag == leaf)
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
                    // A multi-frame file contributes all of its frames to
                    // the same series slice list (RFC 0001 v1 item 2).
                    slices.extend(DicomParser::parse_bytes_all(&bytes).map_err(|e| {
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

#[cfg(test)]
mod multiframe_tests {
    use crate::parser::{encode_element_explicit, DicomParser};
    use crate::tags::{
        Tag, Vr, BITS_ALLOCATED, COLUMNS, IMAGE_ORIENTATION_PATIENT, IMAGE_POSITION_PATIENT,
        NUMBER_OF_FRAMES, PER_FRAME_FUNCTIONAL_GROUPS, PIXEL_DATA, PIXEL_MEASURES_SEQUENCE,
        PIXEL_REPRESENTATION, PIXEL_SPACING, PIXEL_VALUE_TRANSFORMATION_SEQUENCE,
        PLANE_ORIENTATION_SEQUENCE, PLANE_POSITION_SEQUENCE, RESCALE_INTERCEPT, RESCALE_SLOPE,
        ROWS, SHARED_FUNCTIONAL_GROUPS, SLICE_THICKNESS, TRANSFER_SYNTAX_UID,
    };

    /// How the functional-group sequences in a fixture are encoded: both
    /// length forms are legal PS3.5, and real writers differ (dcm4che tends
    /// to undefined length, GDCM to defined), so both must work.
    #[derive(Clone, Copy, PartialEq)]
    enum Style {
        Defined,
        Undefined,
    }

    fn item(content: &[u8]) -> Vec<u8> {
        let mut out = vec![0xFE, 0xFF, 0x00, 0xE0]; // (FFFE,E000)
        out.extend_from_slice(&(content.len() as u32).to_le_bytes());
        out.extend_from_slice(content);
        out
    }

    fn item_undefined(content: &[u8]) -> Vec<u8> {
        let mut out = vec![0xFE, 0xFF, 0x00, 0xE0]; // (FFFE,E000)
        out.extend_from_slice(&0xFFFF_FFFFu32.to_le_bytes());
        out.extend_from_slice(content);
        out.extend_from_slice(&[0xFE, 0xFF, 0x0D, 0xE0]); // (FFFE,E00D)
        out.extend_from_slice(&0u32.to_le_bytes());
        out
    }

    fn seq(tag: Tag, implicit: bool, style: Style, items: &[Vec<u8>]) -> Vec<u8> {
        let body: Vec<u8> = items.concat();
        let mut out = Vec::new();
        out.extend_from_slice(&tag.0.to_le_bytes());
        out.extend_from_slice(&tag.1.to_le_bytes());
        if !implicit {
            // Implicit VR has no VR field: the header is tag + 4-byte
            // length (PS3.5 §6.2.2), so the VR code is written only in
            // the explicit encoding.
            out.extend_from_slice(b"SQ");
            out.extend_from_slice(&[0, 0]);
        }
        match style {
            Style::Defined => out.extend_from_slice(&(body.len() as u32).to_le_bytes()),
            Style::Undefined => out.extend_from_slice(&0xFFFF_FFFFu32.to_le_bytes()),
        }
        out.extend_from_slice(&body);
        if style == Style::Undefined {
            out.extend_from_slice(&[0xFE, 0xFF, 0xDD, 0xE0]); // (FFFE,E0DD)
            out.extend_from_slice(&0u32.to_le_bytes());
        }
        out
    }

    /// A `rows=2, cols=2`, 16-bit, `n`-frame Enhanced-style dataset in the
    /// given encoding. Frame `i` sits at z = 1.5·i with stored values
    /// `1000·(i+1) + k`; frame 1 carries its own rescale (slope 2,
    /// intercept −100).
    fn dataset(style: Style, implicit: bool, declared: u32, n: usize) -> Vec<u8> {
        let padded = |s: &str| {
            let mut v = s.as_bytes().to_vec();
            if v.len() % 2 == 1 {
                v.push(b' ');
            }
            v
        };
        // One element in the dataset's encoding, with the fixture tag's VR.
        let el = |tag: Tag, vr: Vr, value: &[u8]| -> Vec<u8> {
            if implicit {
                implicit_el(tag, value)
            } else {
                encode_element_explicit(tag, vr, value)
            }
        };
        let ds = |tag: Tag, s: &str| el(tag, Vr::Ds, &padded(s));
        let us = |tag: Tag, v: u16| el(tag, Vr::Us, &v.to_le_bytes());
        let pixel = |value: &[u8]| -> Vec<u8> {
            if implicit {
                implicit_el(PIXEL_DATA, value)
            } else {
                encode_element_explicit(PIXEL_DATA, Vr::Ob, value)
            }
        };

        let mut d = Vec::new();
        d.extend(us(ROWS, 2));
        d.extend(us(COLUMNS, 2));
        d.extend(us(BITS_ALLOCATED, 16));
        d.extend(us(PIXEL_REPRESENTATION, 0));
        d.extend(ds(PIXEL_SPACING, "0.8\u{5c}0.8"));
        d.extend(ds(SLICE_THICKNESS, "1.5"));
        d.extend(ds(
            IMAGE_ORIENTATION_PATIENT,
            "1\u{5c}0\u{5c}0\u{5c}0\u{5c}1\u{5c}0",
        ));
        d.extend(el(NUMBER_OF_FRAMES, Vr::Is, &padded(&declared.to_string())));

        let shared = seq(
            SHARED_FUNCTIONAL_GROUPS,
            implicit,
            style,
            &[
                item(&seq(
                    PIXEL_MEASURES_SEQUENCE,
                    implicit,
                    style,
                    &[item(&{
                        let mut v = ds(PIXEL_SPACING, "0.8\u{5c}0.8");
                        v.extend(ds(SLICE_THICKNESS, "1.5"));
                        v
                    })],
                )),
                item(&seq(
                    PLANE_ORIENTATION_SEQUENCE,
                    implicit,
                    style,
                    &[item(&ds(
                        IMAGE_ORIENTATION_PATIENT,
                        "1\u{5c}0\u{5c}0\u{5c}0\u{5c}1\u{5c}0",
                    ))],
                )),
            ],
        );

        let mut per_frame = Vec::new();
        for i in 0..n {
            let mut one = seq(
                PLANE_POSITION_SEQUENCE,
                implicit,
                style,
                &[item(&ds(
                    IMAGE_POSITION_PATIENT,
                    &format!("0\u{5c}0\u{5c}{}", i as f64 * 1.5),
                ))],
            );
            if i == 1 {
                one.extend(seq(
                    PIXEL_VALUE_TRANSFORMATION_SEQUENCE,
                    implicit,
                    style,
                    &[item(&{
                        let mut v = ds(RESCALE_SLOPE, "2");
                        v.extend(ds(RESCALE_INTERCEPT, "-100"));
                        v
                    })],
                ));
            }
            per_frame.push(if style == Style::Defined {
                item(&one)
            } else {
                item_undefined(&one)
            });
        }
        d.extend(shared);
        d.extend(seq(
            PER_FRAME_FUNCTIONAL_GROUPS,
            implicit,
            style,
            &per_frame,
        ));

        let mut px = Vec::new();
        for i in 0..n {
            for k in 0..4u16 {
                px.extend_from_slice(&(1000 * (i as u16 + 1) + k).to_le_bytes());
            }
        }
        d.extend(pixel(&px));
        d
    }

    fn implicit_el(tag: Tag, value: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(8 + value.len());
        out.extend_from_slice(&tag.0.to_le_bytes());
        out.extend_from_slice(&tag.1.to_le_bytes());
        out.extend_from_slice(&(value.len() as u32).to_le_bytes());
        out.extend_from_slice(value);
        out
    }

    fn part10(dataset: &[u8], transfer_syntax: &str) -> Vec<u8> {
        let mut buf = vec![0u8; 128];
        buf.extend_from_slice(b"DICM");
        let mut uid = transfer_syntax.as_bytes().to_vec();
        if uid.len() % 2 == 1 {
            uid.push(0);
        }
        buf.extend(&encode_element_explicit(TRANSFER_SYNTAX_UID, Vr::Ui, &uid));
        buf.extend_from_slice(dataset);
        buf
    }

    #[test]
    fn multiframe_maps_functional_groups_onto_the_slice_list() {
        let slices = DicomParser::parse_bytes_all(&part10(
            &dataset(Style::Defined, false, 3, 3),
            "1.2.840.10008.1.2.1",
        ))
        .expect("parses");
        assert_eq!(slices.len(), 3);
        for (i, s) in slices.iter().enumerate() {
            assert_eq!(s.rows, 2);
            assert_eq!(s.columns, 2);
            // Per-frame plane position.
            assert!((s.position.z - 1.5 * i as f64).abs() < 1e-9, "frame {i}");
            // Shared orientation and measures.
            assert_eq!(s.orientation, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
            assert_eq!(s.pixel_spacing, (0.8, 0.8));
            assert!((s.slice_thickness - 1.5).abs() < 1e-12);
            // Per-frame pixel payload cut from the one native value.
            let expected: Vec<i32> = (0..4u16)
                .map(|k| (1000 * (i as u16 + 1) + k) as i32)
                .collect();
            assert_eq!(s.pixel_data, expected, "frame {i}");
            assert_eq!(s.instance_number, i as u32);
        }
        // Frame 1's own Pixel Value Transformation (slope 2, intercept
        // -100) applies to that frame alone; the others keep 1/0.
        assert!((slices[0].hu_at(1, 1).unwrap() - 1003.0).abs() < 1e-9);
        assert!((slices[1].hu_at(1, 1).unwrap() - (2003.0 * 2.0 - 100.0)).abs() < 1e-9);
        assert!((slices[2].hu_at(1, 1).unwrap() - 3003.0).abs() < 1e-9);
    }

    #[test]
    fn undefined_length_functional_groups_parse_identically() {
        let defined = DicomParser::parse_bytes_all(&part10(
            &dataset(Style::Defined, false, 3, 3),
            "1.2.840.10008.1.2.1",
        ))
        .expect("defined");
        let undefined = DicomParser::parse_bytes_all(&part10(
            &dataset(Style::Undefined, false, 3, 3),
            "1.2.840.10008.1.2.1",
        ))
        .expect("undefined");
        assert_eq!(defined.len(), undefined.len());
        for (a, b) in defined.iter().zip(&undefined) {
            assert_eq!(a.position, b.position);
            assert_eq!(a.pixel_data, b.pixel_data);
            assert_eq!(a.rescale_slope, b.rescale_slope);
        }
    }

    #[test]
    fn implicit_vr_multiframe_parses() {
        let slices = DicomParser::parse_bytes_all(&part10(
            &dataset(Style::Defined, true, 3, 3),
            "1.2.840.10008.1.2",
        ))
        .expect("parses");
        assert_eq!(slices.len(), 3);
        assert!((slices[2].position.z - 3.0).abs() < 1e-9);
        assert_eq!(slices[1].pixel_data, vec![2000, 2001, 2002, 2003]);
    }

    #[test]
    fn implicit_vr_single_frame_parses() {
        // The implicit-VR dataset encoding at the Part-10 level: the meta
        // group must hand over to the implicit reader at the boundary, not
        // read the dataset's first element header as an explicit VR.
        let slice = DicomParser::parse_bytes(&part10(
            &dataset(Style::Defined, true, 1, 1),
            "1.2.840.10008.1.2",
        ))
        .expect("parses");
        assert_eq!(slice.rows, 2);
        assert_eq!(slice.pixel_data, vec![1000, 1001, 1002, 1003]);
    }

    #[test]
    fn parse_bytes_names_parse_bytes_all_for_multiframe() {
        let err = DicomParser::parse_bytes(&part10(
            &dataset(Style::Defined, false, 3, 3),
            "1.2.840.10008.1.2.1",
        ))
        .expect_err("multi-frame must not silently truncate");
        let msg = match &err {
            crate::DicomError::InconsistentSeries(m) => m.clone(),
            other => panic!("unexpected error: {other:?}"),
        };
        assert!(msg.contains("parse_bytes_all"), "{msg}");
    }

    #[test]
    fn declared_frame_count_must_match_per_frame_items() {
        let err = DicomParser::parse_bytes_all(&part10(
            &dataset(Style::Defined, false, 3, 2),
            "1.2.840.10008.1.2.1",
        ))
        .expect_err("declared 3 frames with 2 items");
        assert!(
            matches!(&err, crate::DicomError::BadValue { reason, .. } if reason.contains("per-frame")),
            "unexpected error: {err:?}"
        );
    }

    #[test]
    fn single_frame_declared_count_still_parses_bytes() {
        let slice = DicomParser::parse_bytes(&part10(
            &dataset(Style::Defined, false, 1, 1),
            "1.2.840.10008.1.2.1",
        ))
        .expect("single-frame entry point");
        assert_eq!(slice.pixel_data, vec![1000, 1001, 1002, 1003]);
    }

    #[test]
    fn multiframe_assembles_into_a_series() {
        let slices = DicomParser::parse_bytes_all(&part10(
            &dataset(Style::Defined, false, 3, 3),
            "1.2.840.10008.1.2.1",
        ))
        .expect("parses");
        let series = crate::DicomSeries::from_slices(slices).expect("series");
        assert_eq!(series.dims(), (2, 2, 3));
        assert_eq!(series.hu_volume().len(), 12);
        // Sorted ascending along the normal (z).
        assert!(series.slices[0].position.z < series.slices[2].position.z);
    }
}
