//! JPEG 2000 (ISO/IEC 15444-1) encapsulated pixel data decoding, via the
//! `pdfluent-jpeg2000` crate (`hayro-jpeg2000`). Enabled by the `jpeg2000`
//! cargo feature.
//!
//! Covers two DICOM transfer syntaxes with one codec:
//!
//! - **JPEG 2000 Lossless Only** (`1.2.840.10008.1.2.4.90`): exact
//!   reconstruction (reversible 5/3 wavelet).
//! - **JPEG 2000** (`1.2.840.10008.1.2.4.91`): lossless *or* lossy depending
//!   on how the codestream was encoded (irreversible 9/7 wavelet plus
//!   quantization for the lossy case). The transfer syntax UID alone does
//!   not say which, and this module does not attempt to tell — callers must
//!   not assume `.91` is exact.
//!
//! # Working around a real limitation: the crate always applies an unsigned
//! # DC level shift
//!
//! `pdfluent-jpeg2000` reads a component's SIZ "signed" bit
//! (ISO/IEC 15444-1 A.5.1) but — by its own admission in its source
//! (`_is_signed` is read and discarded with the comment "No idea how to
//! process signed images, but ... let's do the same [as OpenJPEG]") —
//! applies the unsigned DC level shift (`+ 2^(precision-1)`) to every
//! component unconditionally, regardless of that bit. For a codestream
//! whose component actually is signed (already zero-centered, needing no
//! shift), that silently adds an unwanted offset to every sample.
//!
//! Since the offset is a fixed, known constant (`2^(precision-1)`, not
//! data-dependent), it is exactly compensable *if* we know independently
//! whether the component was really declared signed. The crate discards
//! that bit, so [`component0_signed`] re-reads it directly from the raw
//! SIZ marker bytes — the same kind of low-level parsing `rle.rs` and
//! `parser.rs` already do elsewhere in this crate — and `decode_frame`
//! subtracts the shift back out when the stream says signed.
//!
//! DICOM ties this to `PixelRepresentation`: per PS3.5, a conformant encoder
//! sets the JPEG 2000 component signed exactly when `PixelRepresentation ==
//! 1`. A file where the codestream's own signed bit and the dataset's
//! `PixelRepresentation` disagree is non-conformant, and there is no safe
//! way to resolve that disagreement — such a file is rejected rather than
//! guessed at.

use crate::error::{DicomError, Result};
use crate::tags::Tag;
use pdfluent_jpeg2000::{DecodeSettings, DecoderContext, Image};

/// Raw codestream signature (SOC + SIZ marker codes), matching the constant
/// `pdfluent-jpeg2000` itself checks `Image::new` against. DICOM
/// encapsulated pixel data is always a raw codestream, never a JP2-boxed
/// file (PS3.5 Annex A.4.4), so this is required, not merely preferred.
const CODESTREAM_MAGIC: [u8; 4] = [0xFF, 0x4F, 0xFF, 0x51];

/// Byte offset of the first component's `Ssiz` field within a raw J2C
/// codestream (ISO/IEC 15444-1 A.5.1), counting from the start of `SOC`:
///
/// `SOC`(2) + SIZ marker(2) + `Lsiz`(2) + `Rsiz`(2) + `Xsiz`(4) + `Ysiz`(4) +
/// `XOsiz`(4) + `YOsiz`(4) + `XTsiz`(4) + `YTsiz`(4) + `XTOsiz`(4) +
/// `YTOsiz`(4) + `Csiz`(2) = 42.
const SSIZ0_OFFSET: usize = 42;

/// Reads whether the first SIZ component is declared signed
/// (ISO/IEC 15444-1 Table A.11, `Ssiz` bit 7), directly from the codestream
/// bytes. Returns `None` if `fragment` is not a raw J2C codestream, or is
/// too short to contain a `Csiz >= 1` SIZ marker.
fn component0_signed(fragment: &[u8]) -> Option<bool> {
    if !fragment.starts_with(&CODESTREAM_MAGIC) {
        return None;
    }
    fragment.get(SSIZ0_OFFSET).map(|&b| (b & 0x80) != 0)
}

/// Decodes one JPEG 2000 frame into raw stored values.
///
/// `rows`, `columns` and `bits_allocated` come from the dataset and are
/// cross-checked against the codestream's own dimensions and precision
/// rather than trusted blindly. `signed` must agree with the codestream's
/// own declared signedness — see the module-level docs above for why a
/// mismatch is rejected rather than resolved by guessing.
pub fn decode_frame(
    fragment: &[u8],
    rows: u16,
    columns: u16,
    bits_allocated: u16,
    signed: u16,
    tag: Tag,
) -> Result<Vec<i32>> {
    let stream_signed = component0_signed(fragment).ok_or_else(|| DicomError::BadValue {
        tag,
        reason: "JPEG 2000 fragment is not a raw codestream, or is too short to contain \
                 a SIZ marker"
            .into(),
    })?;
    if stream_signed != (signed == 1) {
        return Err(DicomError::BadValue {
            tag,
            reason: format!(
                "JPEG 2000 codestream declares component 0 {}, but the dataset's \
                 PixelRepresentation says {} — non-conformant file, refusing to guess",
                if stream_signed { "signed" } else { "unsigned" },
                if signed == 1 { "signed" } else { "unsigned" }
            ),
        });
    }

    let settings = DecodeSettings::default();
    let image = Image::new(fragment, &settings).map_err(|e| DicomError::BadValue {
        tag,
        reason: format!("JPEG 2000 header parse failed: {e}"),
    })?;

    if image.width() != columns as u32 || image.height() != rows as u32 {
        return Err(DicomError::BadValue {
            tag,
            reason: format!(
                "JPEG 2000 frame is {}x{}, dataset declares {columns}x{rows}",
                image.width(),
                image.height()
            ),
        });
    }

    let mut ctx = DecoderContext::default();
    let decoded = image.decode(&mut ctx).map_err(|e| DicomError::BadValue {
        tag,
        reason: format!("JPEG 2000 decode failed: {e}"),
    })?;

    let components = decoded.components();
    if components.len() != 1 {
        return Err(DicomError::BadValue {
            tag,
            reason: format!(
                "JPEG 2000 frame has {} components; only single-component \
                 grayscale is supported",
                components.len()
            ),
        });
    }
    let component = &components[0];
    // BitsAllocated is the DICOM pixel container width (8 or 16); the
    // codestream's own precision (BitsStored-equivalent) may be narrower
    // than the container, but must not be wider, or values would not fit.
    let bit_depth = u16::from(component.bit_depth());
    if bit_depth > bits_allocated {
        return Err(DicomError::BadValue {
            tag,
            reason: format!(
                "JPEG 2000 component precision is {bit_depth} bits, wider than the \
                 dataset's BitsAllocated {bits_allocated}"
            ),
        });
    }

    let samples = component.samples();
    let expected = rows as usize * columns as usize;
    if samples.len() != expected {
        return Err(DicomError::BadValue {
            tag,
            reason: format!(
                "JPEG 2000 frame decoded {} samples, expected {expected} ({rows}x{columns})",
                samples.len()
            ),
        });
    }

    if stream_signed {
        // Undo the crate's unconditional unsigned level shift: the true
        // signed sample is `decoded - 2^(bit_depth-1)`, per the module docs.
        let shift = (1i64 << (bit_depth - 1)) as f32;
        let min = -(1i64 << (bit_depth - 1));
        let max = (1i64 << (bit_depth - 1)) - 1;
        Ok(samples
            .iter()
            .map(|&s| ((s - shift).round() as i64).clamp(min, max) as i32)
            .collect())
    } else {
        let max = (1i64 << bit_depth) - 1;
        Ok(samples
            .iter()
            .map(|&s| (s.round() as i64).clamp(0, max) as i32)
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TAG: Tag = (0x7FE0, 0x0010);

    // Minimal valid JPEG 2000 raw codestream (J2C) for a 2x2 greyscale image,
    // all coefficients zero (a single empty packet), so every decoded sample
    // is exactly the unsigned DC level-shift midpoint: 2^(8-1) = 128.
    //
    // Layout: SOC, SIZ (8-bit unsigned single component), COD (no MCT,
    // reversible 5/3, 1 layer), QCD (no quantization), then one tile with a
    // zero-length packet and EOC. This mirrors the fixture
    // `pdfluent-jpeg2000` ships in its own `lib.rs` tests, since a
    // hand-built J2C header is otherwise easy to get subtly wrong.
    #[rustfmt::skip]
    const MINIMAL_J2C_2X2: &[u8] = &[
        // SOC
        0xFF, 0x4F,
        // SIZ (FF 51), Lsiz=41
        0xFF, 0x51,
        0x00, 0x29,
        0x00, 0x00,             // Rsiz
        0x00, 0x00, 0x00, 0x02, // Xsiz
        0x00, 0x00, 0x00, 0x02, // Ysiz
        0x00, 0x00, 0x00, 0x00, // XOsiz
        0x00, 0x00, 0x00, 0x00, // YOsiz
        0x00, 0x00, 0x00, 0x02, // XTsiz
        0x00, 0x00, 0x00, 0x02, // YTsiz
        0x00, 0x00, 0x00, 0x00, // XTOsiz
        0x00, 0x00, 0x00, 0x00, // YTOsiz
        0x00, 0x01,             // Csiz = 1
        0x07, 0x01, 0x01,       // Ssiz=7 (8-bit unsigned), XRsiz=1, YRsiz=1
        // COD (FF 52), Lcod=12
        0xFF, 0x52,
        0x00, 0x0C,
        0x00,                   // Scod
        0x00,                   // progression order (LRCP)
        0x00, 0x01,             // num layers
        0x00,                   // MCT = 0
        0x00,                   // decomposition levels = 0
        0x00,                   // code block width
        0x00,                   // code block height
        0x00,                   // code block style
        0x01,                   // transform = 1 (reversible 5/3)
        // QCD (FF 5C), Lqcd=4
        0xFF, 0x5C,
        0x00, 0x04,
        0x00,                   // Sqcd = NoQuantization
        0x80,                   // step-size[0]
        // SOT (FF 90), Lsot=10
        0xFF, 0x90,
        0x00, 0x0A,
        0x00, 0x00,             // Isot
        0x00, 0x00, 0x00, 0x0F, // Psot
        0x00,                   // TPsot
        0x01,                   // TNsot
        // SOD
        0xFF, 0x93,
        0x00,                   // one empty packet
        // EOC
        0xFF, 0xD9,
    ];

    /// Same fixture as [`MINIMAL_J2C_2X2`], with `Ssiz` bit 7 set so
    /// component 0 declares itself signed (byte 42: `0x07` -> `0x87`).
    fn signed_variant() -> Vec<u8> {
        let mut v = MINIMAL_J2C_2X2.to_vec();
        assert_eq!(v[SSIZ0_OFFSET], 0x07, "fixture layout changed");
        v[SSIZ0_OFFSET] = 0x87;
        v
    }

    #[test]
    fn decodes_zero_coefficient_tile_to_midpoint() {
        let decoded = decode_frame(MINIMAL_J2C_2X2, 2, 2, 8, 0, TAG).expect("decodes");
        assert_eq!(decoded, vec![128, 128, 128, 128]);
    }

    /// A signed component's zero-coefficient tile decodes (via the crate) to
    /// the same midpoint 128 the unsigned case does, since the crate applies
    /// the same unconditional shift either way. After `decode_frame`
    /// subtracts that shift back out for a signed component, the correct
    /// result is 0 -- the actual reconstructed signed sample.
    #[test]
    fn decodes_signed_component_with_shift_correction() {
        let stream = signed_variant();
        let decoded = decode_frame(&stream, 2, 2, 8, 1, TAG).expect("decodes");
        assert_eq!(decoded, vec![0, 0, 0, 0]);
    }

    #[test]
    fn unsigned_stream_with_signed_pixel_representation_is_rejected() {
        let err = decode_frame(MINIMAL_J2C_2X2, 2, 2, 8, 1, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }

    #[test]
    fn signed_stream_with_unsigned_pixel_representation_is_rejected() {
        let stream = signed_variant();
        let err = decode_frame(&stream, 2, 2, 8, 0, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }

    #[test]
    fn dimension_mismatch_is_rejected() {
        let err = decode_frame(MINIMAL_J2C_2X2, 3, 3, 8, 0, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }

    #[test]
    fn garbage_is_rejected_not_panicking() {
        let err = decode_frame(&[0u8; 8], 2, 2, 8, 0, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }
}
