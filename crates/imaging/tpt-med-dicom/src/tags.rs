//! Well-known DICOM tag constants and implicit-VR lookups.
//!
//! Only the subset the pipeline consumes is tabulated; unknown tags are
//! still parsed (skipped) correctly.

/// (group, element) for a DICOM tag.
pub type Tag = (u16, u16);

/// Meta group: Transfer Syntax UID (0002,0010).
pub const TRANSFER_SYNTAX_UID: Tag = (0x0002, 0x0010);
/// Study Date (0008,0020).
pub const STUDY_DATE: Tag = (0x0008, 0x0020);
/// Modality (0008,0060).
pub const MODALITY: Tag = (0x0008, 0x0060);
/// Patient ID (0010,0020).
pub const PATIENT_ID: Tag = (0x0010, 0x0020);
/// Slice Thickness (0018,0050).
pub const SLICE_THICKNESS: Tag = (0x0018, 0x0050);
/// Series Instance UID (0020,000E).
pub const SERIES_INSTANCE_UID: Tag = (0x0020, 0x000E);
/// Instance Number (0020,0013).
pub const INSTANCE_NUMBER: Tag = (0x0020, 0x0013);
/// Image Position (Patient) (0020,0032).
pub const IMAGE_POSITION_PATIENT: Tag = (0x0020, 0x0032);
/// Image Orientation (Patient) (0020,0037).
pub const IMAGE_ORIENTATION_PATIENT: Tag = (0x0020, 0x0037);
/// Rows (0028,0010).
pub const ROWS: Tag = (0x0028, 0x0010);
/// Columns (0028,0011).
pub const COLUMNS: Tag = (0x0028, 0x0011);
/// Pixel Spacing (0028,0030).
pub const PIXEL_SPACING: Tag = (0x0028, 0x0030);
/// Bits Allocated (0028,0100).
pub const BITS_ALLOCATED: Tag = (0x0028, 0x0100);
/// Pixel Representation (0028,0103) — 0 unsigned, 1 two's-complement.
pub const PIXEL_REPRESENTATION: Tag = (0x0028, 0x0103);
/// Rescale Intercept (0028,1052).
pub const RESCALE_INTERCEPT: Tag = (0x0028, 0x1052);
/// Rescale Slope (0028,1053).
pub const RESCALE_SLOPE: Tag = (0x0028, 0x1053);
/// Pixel Data (7FE0,0010).
pub const PIXEL_DATA: Tag = (0x7FE0, 0x0010);

/// Value Representation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Vr {
    /// Application Entity.
    Ae,
    /// Age String.
    As,
    /// Attribute Tag.
    At,
    /// Code String.
    Cs,
    /// Date.
    Da,
    /// Decimal String.
    Ds,
    /// Date Time.
    Dt,
    /// Floating Point Single.
    Fl,
    /// Floating Point Double.
    Fd,
    /// Integer String.
    Is,
    /// Long String.
    Lo,
    /// Other Byte.
    Ob,
    /// Other Double.
    Od,
    /// Other Float.
    Of,
    /// Other Long.
    Ol,
    /// Other Very Long.
    Ov,
    /// Other Word.
    Ow,
    /// Person Name.
    Pn,
    /// Short String.
    Sh,
    /// Signed Long.
    Sl,
    /// Sequence of Items.
    Sq,
    /// Short Text.
    St,
    /// Signed Short.
    Ss,
    /// Short Vr, undefined-length capable.
    Sv,
    /// Time.
    Tm,
    /// Unlimited Characters.
    Uc,
    /// Unique Identifier (UID).
    Ui,
    /// Unsigned Long.
    Ul,
    /// Unknown.
    Un,
    /// Universal Resource Identifier.
    Ur,
    /// Unsigned Short.
    Us,
    /// Unlimited Text.
    Ut,
    /// Very Long.
    Uv,
}

impl Vr {
    /// Parses a two-byte ASCII VR code.
    pub fn from_ascii(b: [u8; 2]) -> Option<Vr> {
        let s = core::str::from_utf8(&b).ok()?;
        Some(match s {
            "AE" => Vr::Ae,
            "AS" => Vr::As,
            "AT" => Vr::At,
            "CS" => Vr::Cs,
            "DA" => Vr::Da,
            "DS" => Vr::Ds,
            "DT" => Vr::Dt,
            "FD" => Vr::Fd,
            "FL" => Vr::Fl,
            "IS" => Vr::Is,
            "LO" => Vr::Lo,
            "OB" => Vr::Ob,
            "OD" => Vr::Od,
            "OF" => Vr::Of,
            "OL" => Vr::Ol,
            "OV" => Vr::Ov,
            "OW" => Vr::Ow,
            "PN" => Vr::Pn,
            "SH" => Vr::Sh,
            "SL" => Vr::Sl,
            "SQ" => Vr::Sq,
            "ST" => Vr::St,
            "SS" => Vr::Ss,
            "SV" => Vr::Sv,
            "TM" => Vr::Tm,
            "UC" => Vr::Uc,
            "UI" => Vr::Ui,
            "UL" => Vr::Ul,
            "UN" => Vr::Un,
            "UR" => Vr::Ur,
            "US" => Vr::Us,
            "UT" => Vr::Ut,
            "UV" => Vr::Uv,
            _ => return None,
        })
    }

    /// The two-byte ASCII code for this VR.
    pub fn ascii(self) -> &'static [u8; 2] {
        match self {
            Vr::Ae => b"AE",
            Vr::As => b"AS",
            Vr::At => b"AT",
            Vr::Cs => b"CS",
            Vr::Da => b"DA",
            Vr::Ds => b"DS",
            Vr::Dt => b"DT",
            Vr::Fd => b"FD",
            Vr::Fl => b"FL",
            Vr::Is => b"IS",
            Vr::Lo => b"LO",
            Vr::Ob => b"OB",
            Vr::Od => b"OD",
            Vr::Of => b"OF",
            Vr::Ol => b"OL",
            Vr::Ov => b"OV",
            Vr::Ow => b"OW",
            Vr::Pn => b"PN",
            Vr::Sh => b"SH",
            Vr::Sl => b"SL",
            Vr::Sq => b"SQ",
            Vr::St => b"ST",
            Vr::Ss => b"SS",
            Vr::Sv => b"SV",
            Vr::Tm => b"TM",
            Vr::Uc => b"UC",
            Vr::Ui => b"UI",
            Vr::Ul => b"UL",
            Vr::Un => b"UN",
            Vr::Ur => b"UR",
            Vr::Us => b"US",
            Vr::Ut => b"UT",
            Vr::Uv => b"UV",
        }
    }

    /// True if this VR uses a 32-bit length field in explicit encoding
    /// (PS3.5 §7.1.2).
    pub fn uses_long_length(self) -> bool {
        matches!(
            self,
            Vr::Ob
                | Vr::Od
                | Vr::Of
                | Vr::Ol
                | Vr::Ov
                | Vr::Ow
                | Vr::Sq
                | Vr::Un
                | Vr::Ut
                | Vr::Uc
                | Vr::Sv
                | Vr::Uv
        )
    }
}

/// The transfer syntax a dataset is encoded with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransferSyntax {
    /// Implicit VR little endian (`1.2.840.10008.1.2`).
    ImplicitVrLittleEndian,
    /// Explicit VR little endian (`1.2.840.10008.1.2.1`).
    ExplicitVrLittleEndian,
}

impl TransferSyntax {
    /// UID → transfer syntax, rejecting compressed/headerless syntaxes with
    /// actionable errors.
    pub fn from_uid(uid: &str) -> Result<Self, String> {
        match uid.trim_end_matches('\0') {
            "1.2.840.10008.1.2" => Ok(Self::ImplicitVrLittleEndian),
            "1.2.840.10008.1.2.1" => Ok(Self::ExplicitVrLittleEndian),
            other
                if other.starts_with("1.2.840.10008.1.2.4.")
                    || other.starts_with("1.2.840.10008.1.2.5.") =>
            {
                Err(format!(
                    "compressed transfer syntax {other} not supported; \
                     decompress before ingestion (see rfcs/0001-dicom-ingestion.md)"
                ))
            }
            other => Err(format!("unknown transfer syntax {other}")),
        }
    }
}

/// VR for a tag under implicit VR little endian, for the tabulated subset.
/// `PixelData` reports `Ow`; `UN`-equivalent unknown tags default to raw
/// bytes which the caller never decodes.
pub fn implicit_vr(tag: Tag) -> Vr {
    match tag {
        TRANSFER_SYNTAX_UID | SERIES_INSTANCE_UID => Vr::Ui,
        STUDY_DATE => Vr::Da,
        MODALITY => Vr::Cs,
        PATIENT_ID => Vr::Lo,
        SLICE_THICKNESS
        | IMAGE_POSITION_PATIENT
        | IMAGE_ORIENTATION_PATIENT
        | PIXEL_SPACING
        | RESCALE_INTERCEPT
        | RESCALE_SLOPE => Vr::Ds,
        INSTANCE_NUMBER => Vr::Is,
        ROWS | COLUMNS | BITS_ALLOCATED | PIXEL_REPRESENTATION => Vr::Us,
        PIXEL_DATA => Vr::Ow,
        _ => Vr::Un,
    }
}
