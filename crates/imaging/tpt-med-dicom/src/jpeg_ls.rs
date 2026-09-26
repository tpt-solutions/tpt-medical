//! JPEG-LS (ISO/IEC 14495-1 / ITU-T T.87) encapsulated pixel data decoding,
//! via the `pure_jpegls` crate. Enabled by the `jpeg-ls` cargo feature.
//!
//! Covers two DICOM transfer syntaxes with one codec:
//!
//! - **JPEG-LS Lossless** (`1.2.840.10008.1.2.4.80`): exact reconstruction.
//! - **JPEG-LS Near-Lossless** (`1.2.840.10008.1.2.4.81`): the encoder may
//!   bound each sample's error by a `NEAR` value greater than zero, so a
//!   decoded value is only guaranteed to be within `NEAR` of the original —
//!   never treated as an exact stored value by this module or its callers.
//!
//! Single-component (grayscale) only, matching `pure_jpegls`'s conformance
//! scope (`Nf = 1`, `ILV = 0`); DICOM CT/MR pixel data is always
//! single-component.

use crate::error::{DicomError, Result};
use crate::tags::Tag;

/// Decodes one JPEG-LS frame into raw stored values.
///
/// `rows`, `columns`, `bits_allocated` and `signed` come from the dataset;
/// the decoded frame's own dimensions are cross-checked against `rows` and
/// `columns` rather than trusted blindly.
pub fn decode_frame(
    fragment: &[u8],
    rows: u16,
    columns: u16,
    bits_allocated: u16,
    signed: u16,
    tag: Tag,
) -> Result<Vec<i32>> {
    if bits_allocated != 8 && bits_allocated != 16 {
        return Err(DicomError::BadValue {
            tag,
            reason: format!(
                "JPEG-LS decoding supports BitsAllocated 8 or 16, got {bits_allocated}"
            ),
        });
    }

    let (samples, width, height) =
        jpegls::decode(fragment, columns as u32, rows as u32).map_err(|e| {
            DicomError::BadValue {
                tag,
                reason: format!("JPEG-LS decode failed: {e}"),
            }
        })?;

    if width != columns as u32 || height != rows as u32 {
        return Err(DicomError::BadValue {
            tag,
            reason: format!("JPEG-LS frame is {width}x{height}, dataset declares {columns}x{rows}"),
        });
    }

    Ok(samples
        .into_iter()
        .map(|v| {
            if bits_allocated == 8 {
                if signed == 1 {
                    v as u8 as i8 as i32
                } else {
                    v as u8 as i32
                }
            } else if signed == 1 {
                v as i16 as i32
            } else {
                v as i32
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    const TAG: Tag = (0x7FE0, 0x0010);

    #[test]
    fn roundtrips_lossless_16bit() {
        let pixels: Vec<u16> = vec![0, 1, 65535, 32768, 100, 4095, 60000, 12];
        let (w, h) = (4u32, 2u32);
        let mut encoded = Vec::new();
        jpegls::encode(&pixels, w, h, &mut encoded).expect("encodes");

        let decoded = decode_frame(&encoded, h as u16, w as u16, 16, 0, TAG).expect("decodes");
        assert_eq!(
            decoded,
            pixels.iter().map(|&v| v as i32).collect::<Vec<_>>()
        );
    }

    #[test]
    fn roundtrips_lossless_signed_16bit() {
        // Encode the bit pattern for [-100, 1, -1, 32767]; JPEG-LS itself is
        // sample-value agnostic (it just codes 16-bit magnitudes), so the
        // sign convention is applied on the way back out, exactly as for the
        // native and RLE pixel-data paths.
        let signed_values: [i16; 4] = [-100, 1, -1, 32767];
        let pixels: Vec<u16> = signed_values.iter().map(|&v| v as u16).collect();
        let (w, h) = (2u32, 2u32);
        let mut encoded = Vec::new();
        jpegls::encode(&pixels, w, h, &mut encoded).expect("encodes");

        let decoded = decode_frame(&encoded, h as u16, w as u16, 16, 1, TAG).expect("decodes");
        assert_eq!(
            decoded,
            signed_values.iter().map(|&v| v as i32).collect::<Vec<_>>()
        );
    }

    #[test]
    fn roundtrips_lossless_8bit() {
        let pixels: Vec<u16> = vec![0, 127, 128, 255, 1, 254];
        let (w, h) = (3u32, 2u32);
        let mut opts = jpegls::EncodeOptions::default();
        opts.precision = Some(8);
        let mut encoded = Vec::new();
        jpegls::encode_with_options(&pixels, w, h, &opts, &mut encoded).expect("encodes");

        let decoded = decode_frame(&encoded, h as u16, w as u16, 8, 0, TAG).expect("decodes");
        assert_eq!(
            decoded,
            pixels.iter().map(|&v| v as i32).collect::<Vec<_>>()
        );
    }

    #[test]
    fn dimension_mismatch_is_rejected() {
        let pixels: Vec<u16> = vec![1, 2, 3, 4];
        let mut encoded = Vec::new();
        jpegls::encode(&pixels, 2, 2, &mut encoded).expect("encodes");

        let err = decode_frame(&encoded, 3, 3, 16, 0, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }

    #[test]
    fn unsupported_bits_allocated_is_rejected() {
        let pixels: Vec<u16> = vec![1, 2, 3, 4];
        let mut encoded = Vec::new();
        jpegls::encode(&pixels, 2, 2, &mut encoded).expect("encodes");

        let err = decode_frame(&encoded, 2, 2, 12, 0, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }

    #[test]
    fn garbage_is_rejected_not_panicking() {
        let err = decode_frame(&[0u8; 8], 2, 2, 16, 0, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }
}
