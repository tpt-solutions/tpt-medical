//! Classic ITU-T T.81 JPEG family encapsulated pixel data decoding, via the
//! `jpeg-decoder` crate (image-rs). Enabled by the `jpeg` cargo feature.
//!
//! One dependency covers two different DICOM transfer syntaxes because
//! `jpeg-decoder` implements both coding processes T.81 defines:
//!
//! - **Baseline / Extended DCT** (Process 1; Process 2 & 4 — DICOM UIDs
//!   `.50` / `.51`). Lossy: the decoded pixel values are only an
//!   approximation of the originals, because DCT quantization discards
//!   information by design. Fine for a thumbnail; not for a value a stress
//!   calculation should trust as an exact stored value.
//! - **Lossless, Process 14 and Process 14 Selection Value 1** (DICOM UIDs
//!   `.57` / `.70`, the latter being the "default lossless JPEG" transfer
//!   syntax). Exact: differential coding (DPCM) with Huffman entropy coding,
//!   not the DCT path, so the decoded stored values are bit-exact.
//!
//! Both single-component (grayscale) frames only — DICOM CT/MR pixel data is
//! always single-component; a JPEG frame with more components (e.g. an RGB
//! secondary-capture image) is rejected rather than reinterpreted.

use std::io::Cursor;

use crate::error::{DicomError, Result};
use crate::tags::Tag;

/// Decodes one JPEG (Baseline/Extended or Lossless Process 14/SV1) frame
/// into raw stored values.
///
/// `rows`, `columns`, `bits_allocated` and `signed` come from the dataset,
/// exactly as for the `rle` and `jpeg-ls` decoders, and the decoded frame's
/// own dimensions/precision are cross-checked against them rather than
/// trusted blindly.
pub fn decode_frame(
    fragment: &[u8],
    rows: u16,
    columns: u16,
    bits_allocated: u16,
    signed: u16,
    tag: Tag,
) -> Result<Vec<i32>> {
    let mut decoder = jpeg_decoder::Decoder::new(Cursor::new(fragment));
    let pixels = decoder.decode().map_err(|e| DicomError::BadValue {
        tag,
        reason: format!("JPEG decode failed: {e}"),
    })?;
    let info = decoder.info().ok_or_else(|| DicomError::BadValue {
        tag,
        reason: "JPEG stream produced no frame header".into(),
    })?;

    if info.width != columns || info.height != rows {
        return Err(DicomError::BadValue {
            tag,
            reason: format!(
                "JPEG frame is {}x{}, dataset declares {columns}x{rows}",
                info.width, info.height
            ),
        });
    }

    match info.pixel_format {
        jpeg_decoder::PixelFormat::L8 => {
            if bits_allocated != 8 {
                return Err(DicomError::BadValue {
                    tag,
                    reason: format!(
                        "JPEG frame decoded 8-bit samples, dataset declares \
                         BitsAllocated {bits_allocated}"
                    ),
                });
            }
            Ok(pixels
                .into_iter()
                .map(|b| {
                    if signed == 1 {
                        b as i8 as i32
                    } else {
                        b as i32
                    }
                })
                .collect())
        }
        jpeg_decoder::PixelFormat::L16 => {
            if bits_allocated != 16 {
                return Err(DicomError::BadValue {
                    tag,
                    reason: format!(
                        "JPEG frame decoded >8-bit samples, dataset declares \
                         BitsAllocated {bits_allocated}"
                    ),
                });
            }
            // `jpeg-decoder` returns >8-bit samples as native-endian u16
            // pairs (documented on `compute_image_lossless`/DCT output
            // conversion), matching this crate's own native pixel-data path.
            Ok(pixels
                .chunks_exact(2)
                .map(|c| {
                    let v = u16::from_ne_bytes([c[0], c[1]]);
                    if signed == 1 {
                        v as i16 as i32
                    } else {
                        v as i32
                    }
                })
                .collect())
        }
        other => Err(DicomError::BadValue {
            tag,
            reason: format!(
                "JPEG frame has pixel format {other:?}; only single-component \
                 grayscale (L8/L16) is supported"
            ),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TAG: Tag = (0x7FE0, 0x0010);

    /// Builds a minimal JPEG Lossless (SOF3, predictor Ra) bitstream for a
    /// 2x2 8-bit image whose samples are `[128, 128, 128, 128 + last_diff]`.
    ///
    /// Restricting to two DC categories (0, for a zero difference, and 1,
    /// for a +/-1 difference) keeps the Huffman table to two length-1 codes
    /// and the "extend" arithmetic (T.81 F.2.2.1) unambiguous, so the
    /// bitstream can be hand-verified rather than round-tripped through the
    /// same encoder the decoder is being tested against.
    ///
    /// Per T.81 H.1.2.1 / `lossless.rs`'s `predict`: pixel (0,0) predicts
    /// from `1 << (precision-1)` = 128; pixel (1,0) (rest of first row)
    /// predicts from its left neighbour; pixel (0,1) (start of a later row)
    /// predicts from the pixel above; pixel (1,1) uses the selected
    /// predictor (Ss=1, Ra = left neighbour). With all base samples at 128,
    /// every prediction is 128, so only the last sample's difference is
    /// non-zero.
    fn lossless_8bit_2x2(last_diff: i8) -> Vec<u8> {
        assert!(last_diff == 0 || last_diff == 1 || last_diff == -1);

        // SOI
        let mut out = vec![0xFF, 0xD8];

        // SOF3: precision=8, height=2, width=2, 1 component (id=1, sampling
        // 1x1, quant table 0).
        out.extend_from_slice(&[
            0xFF, 0xC3, 0x00, 0x0B, 0x08, 0x00, 0x02, 0x00, 0x02, 0x01, 0x01, 0x11, 0x00,
        ]);

        // DHT: one DC table (id 0), two length-1 codes.
        // BITS[1] = 2; HUFFVAL = [category 0, category 1] -> canonical codes
        // "0" -> category 0 (zero difference), "1" -> category 1 (+/-1).
        let mut dht = vec![0xFF, 0xC4];
        let mut payload = vec![0x00u8]; // Tc|Th = 0 (DC, table 0)
        let mut bits = [0u8; 16];
        bits[0] = 2; // BITS[1] (two codes of length 1)
        payload.extend_from_slice(&bits);
        payload.extend_from_slice(&[0x00, 0x01]); // HUFFVAL
        let len = (payload.len() + 2) as u16;
        dht.extend_from_slice(&len.to_be_bytes());
        dht.extend_from_slice(&payload);
        out.extend_from_slice(&dht);

        // SOS: Ns=1, (Cs=1, Td|Ta=0x00), Ss=1 (predictor Ra), Se=0, Ah|Al=0.
        out.extend_from_slice(&[0xFF, 0xDA, 0x00, 0x08, 0x01, 0x01, 0x00, 0x01, 0x00, 0x00]);

        // Entropy scan: three zero-category codes ("0"), then either another
        // "0" (last_diff == 0) or a category-1 code ("1") plus its 1-bit
        // extra ("1" -> +1, "0" -> -1 per the extend() mapping).
        let mut bits_out: Vec<u8> = vec![0, 0, 0];
        if last_diff == 0 {
            bits_out.push(0);
        } else {
            bits_out.push(1);
            bits_out.push(if last_diff == 1 { 1 } else { 0 });
        }
        let mut byte = 0u8;
        let mut nbits = 0u8;
        let mut scan = Vec::new();
        for b in bits_out {
            byte = (byte << 1) | b;
            nbits += 1;
            if nbits == 8 {
                scan.push(byte);
                if byte == 0xFF {
                    scan.push(0x00); // byte-stuffing
                }
                byte = 0;
                nbits = 0;
            }
        }
        if nbits > 0 {
            byte <<= 8 - nbits;
            byte |= (1u8 << (8 - nbits)) - 1; // pad with 1s per T.81 F.1.2.3
            scan.push(byte);
            if byte == 0xFF {
                scan.push(0x00);
            }
        }
        out.extend_from_slice(&scan);
        out.extend_from_slice(&[0xFF, 0xD9]); // EOI
        out
    }

    #[test]
    fn decodes_lossless_8bit_exact_zero_diff() {
        let stream = lossless_8bit_2x2(0);
        let decoded = decode_frame(&stream, 2, 2, 8, 0, TAG).expect("decodes");
        assert_eq!(decoded, vec![128, 128, 128, 128]);
    }

    #[test]
    fn decodes_lossless_8bit_exact_positive_diff() {
        let stream = lossless_8bit_2x2(1);
        let decoded = decode_frame(&stream, 2, 2, 8, 0, TAG).expect("decodes");
        assert_eq!(decoded, vec![128, 128, 128, 129]);
    }

    #[test]
    fn decodes_lossless_8bit_exact_negative_diff() {
        let stream = lossless_8bit_2x2(-1);
        let decoded = decode_frame(&stream, 2, 2, 8, 0, TAG).expect("decodes");
        assert_eq!(decoded, vec![128, 128, 128, 127]);
    }

    #[test]
    fn dimension_mismatch_is_rejected() {
        let stream = lossless_8bit_2x2(0);
        let err = decode_frame(&stream, 3, 3, 8, 0, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }

    #[test]
    fn garbage_is_rejected_not_panicking() {
        let err = decode_frame(&[0u8; 8], 2, 2, 8, 0, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }
}
