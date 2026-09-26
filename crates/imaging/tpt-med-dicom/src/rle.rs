//! RLE Lossless (1.2.840.10008.1.2.5) encapsulated pixel data decoding.
//!
//! Enabled by the `rle` cargo feature. RLE is the only compressed syntax
//! this crate decodes, and it is the only one that needs no external
//! dependency: the pixel data is PackBits-compressed, per PS3.5 Annex G.
//! JPEG, JPEG-LS and JPEG 2000 all require a codec crate and are still
//! rejected (see `rfcs/0001-dicom-ingestion.md`).
//!
//! # Why RLE gets decoded and JPEG 2000 does not
//!
//! RLE is a byte-level format this crate can implement and verify against
//! the standard directly. A lossless JPEG decoder is a large piece of
//! numerical code whose correctness would need its own verification
//! evidence, and a subtly wrong decoder here would be worse than an honest
//! error: it would feed plausible-looking wrong numbers into a stress field
//! that a surgeon might act on. So JPEG 2000 stays a rejection until it can
//! be brought in behind a vetted dependency.
//!
//! # Format
//!
//! Each encapsulated fragment is one frame:
//!
//! 1. A 64-byte header: a `u32` segment count, then 15 `u32` byte offsets
//!    from the start of the frame to each segment. 16-bit pixel data uses
//!    segment 0 for the high byte and segment 1 for the low byte; the
//!    remaining segments must be zero-length.
//! 2. Each segment is PackBits: a control byte `n` in `0..=127` copies the
//!    next `n + 1` bytes literally, `n` in `129..=255` repeats the next byte
//!    `257 - n` times, and `n == 128` is a no-op.

use crate::error::{DicomError, Result};
use crate::tags::Tag;

/// The offset table addresses 15 segments (PS3.5 Annex G.1).
const MAX_SEGMENTS: usize = 15;
/// Header length: one segment count plus 15 offsets.
const HEADER_LEN: usize = 4 + 4 * MAX_SEGMENTS;

/// Decodes one RLE-compressed frame into raw stored values.
///
/// `rows`, `columns` and `bits_allocated` come from the dataset and bound the
/// output: the frame decodes into exactly `rows * columns` samples, and
/// anything shorter is an error. That bound is what keeps a malformed or
/// hostile fragment from expanding without limit.
pub fn decode_frame(
    frame: &[u8],
    rows: u16,
    columns: u16,
    bits_allocated: u16,
    signed: u16,
    tag: Tag,
) -> Result<Vec<i32>> {
    if bits_allocated != 16 && bits_allocated != 8 {
        return Err(DicomError::BadValue {
            tag,
            reason: format!("RLE decoding supports BitsAllocated 8 or 16, got {bits_allocated}"),
        });
    }
    if frame.len() < HEADER_LEN {
        return Err(DicomError::BadValue {
            tag,
            reason: format!(
                "RLE frame is {} bytes, shorter than the {HEADER_LEN}-byte header",
                frame.len()
            ),
        });
    }

    let expected = rows as usize * columns as usize;
    if expected == 0 {
        return Ok(Vec::new());
    }

    let segment_count = u32::from_le_bytes([frame[0], frame[1], frame[2], frame[3]]) as usize;
    if segment_count == 0 || segment_count > MAX_SEGMENTS {
        return Err(DicomError::BadValue {
            tag,
            reason: format!(
                "RLE header declares {segment_count} segments, expected 1..={MAX_SEGMENTS}"
            ),
        });
    }

    let frame_len = frame.len();
    let segment_offset = |i: usize| -> Result<usize> {
        let p = 4 + 4 * i;
        let o = u32::from_le_bytes([frame[p], frame[p + 1], frame[p + 2], frame[p + 3]]) as usize;
        if o < HEADER_LEN || o > frame_len {
            return Err(DicomError::BadValue {
                tag,
                reason: format!("RLE segment {i} offset {o} is outside the {frame_len}-byte frame"),
            });
        }
        Ok(o)
    };

    // Each segment ends where the next one begins, so a run may not read into
    // a neighbouring segment (nor into the next frame's bytes in the
    // concatenated fragment buffer). Only the final segment runs to the end.
    let segment_end = |i: usize| -> usize {
        if i + 1 < segment_count {
            segment_offset(i + 1).unwrap_or(frame_len)
        } else {
            frame_len
        }
    };

    if bits_allocated == 8 {
        let start = segment_offset(0)?;
        let bytes = unpack_bits(frame, start..segment_end(0), expected, 0, tag)?;
        return Ok(bytes
            .into_iter()
            .map(|b| {
                if signed == 1 {
                    b as i8 as i32
                } else {
                    b as i32
                }
            })
            .collect());
    }

    // 16-bit: most significant byte plane first, then least significant.
    let h_start = segment_offset(0)?;
    let h_end = segment_end(0);
    let l_start = segment_offset(1)?;
    let l_end = segment_end(1);
    let high = unpack_bits(frame, h_start..h_end, expected, 0, tag)?;
    let low = unpack_bits(frame, l_start..l_end, expected, 1, tag)?;
    Ok(high
        .iter()
        .zip(&low)
        .map(|(&h, &l)| {
            let v = ((h as u16) << 8 | l as u16) as i16;
            if signed == 1 {
                v as i32
            } else {
                v as u16 as i32
            }
        })
        .collect())
}

/// PackBits-decodes one segment over `bounds` (`start..end` of the frame).
///
/// `expected` caps the output: a run that would overshoot is truncated to
/// reach it exactly, and trailing bytes in the segment are ignored. Without
/// that cap a corrupt run length could expand to gigabytes before the
/// mismatch with `rows * columns` was ever noticed.
fn unpack_bits(
    frame: &[u8],
    bounds: core::ops::Range<usize>,
    expected: usize,
    segment: usize,
    tag: Tag,
) -> Result<Vec<u8>> {
    // Destructured as a struct, not via `into_iter`: a `Range` used as an
    // iterator yields every element it spans, not its start and end.
    // `Range::into_inner` is unstable, so this is the MSRV-safe equivalent.
    let core::ops::Range { start, end } = bounds;
    let mut out = Vec::with_capacity(expected);
    let mut i = start;
    while out.len() < expected {
        let short = || DicomError::BadValue {
            tag,
            reason: format!(
                "RLE segment {segment} ended after {} of {expected} bytes",
                out.len()
            ),
        };
        let ctrl = *frame.get(i).filter(|_| i < end).ok_or_else(short)?;
        i += 1;
        match ctrl {
            0..=127 => {
                let n = ctrl as usize + 1;
                if i + n > end {
                    return Err(DicomError::BadValue {
                        tag,
                        reason: format!("RLE segment {segment} literal run runs past the segment"),
                    });
                }
                let take = n.min(expected - out.len());
                out.extend_from_slice(&frame[i..i + take]);
                i += n;
            }
            // A control byte of 128 is a documented no-op, not a 129-byte run.
            128 => {}
            _ => {
                let n = 1 - (ctrl as i8) as isize;
                let b = *frame
                    .get(i)
                    .filter(|_| i < end)
                    .ok_or_else(|| DicomError::BadValue {
                        tag,
                        reason: format!("RLE segment {segment} repeat is missing its value byte"),
                    })?;
                i += 1;
                let take = n.min(expected as isize - out.len() as isize);
                if take > 0 {
                    out.resize(out.len() + take as usize, b);
                }
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TAG: Tag = (0x7FE0, 0x0010);

    /// Builds an RLE frame from PackBits-encoded segments.
    fn frame(segments: &[Vec<u8>]) -> Vec<u8> {
        let mut offsets = Vec::new();
        let mut body = Vec::new();
        for seg in segments {
            offsets.push((HEADER_LEN + body.len()) as u32);
            body.extend_from_slice(seg);
        }
        let mut out = vec![0u8; HEADER_LEN];
        out[0..4].copy_from_slice(&(segments.len() as u32).to_le_bytes());
        for (i, off) in offsets.iter().enumerate() {
            let p = 4 + 4 * i;
            out[p..p + 4].copy_from_slice(&off.to_le_bytes());
        }
        out.extend_from_slice(&body);
        out
    }

    /// PackBits literal run of `data`.
    fn literal(data: &[u8]) -> Vec<u8> {
        let mut out = vec![(data.len() as u8) - 1];
        out.extend_from_slice(data);
        out
    }

    /// PackBits repeat of `b` exactly `n` times (2 <= n <= 128). A single-copy
    /// repeat has no encoding: control 0x00 means a 1-byte *literal*.
    fn repeat(b: u8, n: usize) -> Vec<u8> {
        assert!((2..=128).contains(&n));
        vec![(257 - n) as u8, b]
    }

    fn decode16(segments: &[Vec<u8>], rows: u16, cols: u16) -> Result<Vec<i32>> {
        decode_frame(&frame(segments), rows, cols, 16, 1, TAG)
    }

    #[test]
    fn roundtrips_16_bit_literals() {
        // 2x2: high bytes 0x00,0x12,0x34,0xFF low bytes 0x10,0x34,0x56,0xFF
        let high = literal(&[0x00, 0x12, 0x34, 0xFF]);
        let low = literal(&[0x10, 0x34, 0x56, 0xFF]);
        let got = decode16(&[high, low], 2, 2).expect("decodes");
        assert_eq!(got, vec![0x0010, 0x1234, 0x3456, -1]);
    }

    #[test]
    fn unpacks_mixed_runs_and_literals() {
        // 2x2. Control 0xFE is -2, so it repeats the next byte 3 times;
        // the trailing 0x00 is a 1-byte literal. Both byte planes mix the two.
        let high = vec![0xFE, 0x7F, 0x00, 0xAB];
        let low = vec![0xFE, 0xFF, 0x00, 0xCD];
        let got = decode16(&[high, low], 2, 2).expect("decodes");
        // decode16 decodes as signed, so 0xABCD sign-extends to a negative i32.
        assert_eq!(got, vec![0x7FFF, 0x7FFF, 0x7FFF, 0xABCDu16 as i16 as i32]);
    }

    #[test]
    fn repeat_control_off_by_one() {
        // Control 0xFF is -1, so the next byte repeats twice, not once. A
        // decoder that read it as 1 would leave the output short.
        let got = decode16(&[repeat(0xAA, 2), repeat(0xBB, 2)], 1, 2).expect("decodes");
        assert_eq!(got, vec![0xAABBu16 as i16 as i32, 0xAABBu16 as i16 as i32]);
    }

    #[test]
    fn noop_control_does_not_consume_output() {
        // A leading 0x80 no-op must not eat a byte of the literal that follows.
        let seg = vec![0x80, 0x01, 0xAA, 0xBB];
        let got = decode16(&[seg, vec![0x01, 0xCC, 0xDD]], 1, 2).expect("decodes");
        assert_eq!(got, vec![0xAACCu16 as i16 as i32, 0xBBDDu16 as i16 as i32]);
    }

    #[test]
    fn unsigned_16_bit_is_not_sign_extended() {
        let got = decode_frame(
            &frame(&[literal(&[0xFF, 0xFF]), literal(&[0xFF, 0xFE])]),
            1,
            2,
            16,
            0,
            TAG,
        )
        .expect("decodes");
        assert_eq!(got, vec![65535, 65534]);
    }

    #[test]
    fn eight_bit_signed_roundtrip() {
        let got = decode_frame(
            &frame(&[literal(&[0x00, 0x7F, 0x80, 0xFF])]),
            2,
            2,
            8,
            1,
            TAG,
        )
        .expect("decodes");
        assert_eq!(got, vec![0, 127, -128, -1]);
    }

    #[test]
    fn short_frame_is_rejected() {
        let err = decode_frame(&[0u8; 16], 2, 2, 16, 1, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }

    #[test]
    fn truncated_segment_is_rejected_not_padded() {
        // Geometry claims 4 samples; the run only supplies 2. Padding with
        // zeroes here would silently invent attenuation values.
        let seg = vec![0x01, 0xAA];
        let err = decode16(&[seg.clone(), seg], 2, 2).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }

    #[test]
    fn oversized_run_is_capped_to_geometry() {
        // A 128-byte run against a 1x1 geometry. The cap must stop at one
        // sample; without it the run would expand unbounded on a corrupt frame.
        let got = decode16(&[repeat(0x7F, 128), repeat(0x00, 128)], 1, 1).expect("decodes");
        assert_eq!(got, vec![0x7F00]);
    }

    #[test]
    fn zero_segment_count_is_rejected() {
        let mut f = frame(&[literal(&[1, 2, 3, 4])]);
        f[0..4].copy_from_slice(&0u32.to_le_bytes());
        assert!(decode_frame(&f, 2, 2, 16, 1, TAG).is_err());
    }

    #[test]
    fn offset_outside_frame_is_rejected() {
        let mut f = frame(&[literal(&[1, 2, 3, 4]), literal(&[5, 6, 7, 8])]);
        f[4..8].copy_from_slice(&9_999u32.to_le_bytes());
        assert!(decode_frame(&f, 2, 2, 16, 1, TAG).is_err());
    }

    #[test]
    fn unsupported_bits_allocated_is_rejected() {
        let f = frame(&[literal(&[0; 16])]);
        let err = decode_frame(&f, 2, 2, 12, 0, TAG).unwrap_err();
        assert!(matches!(err, DicomError::BadValue { .. }));
    }
}
