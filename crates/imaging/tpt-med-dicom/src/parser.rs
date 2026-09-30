//! Byte-level DICOM dataset parser (uncompressed little-endian syntaxes).

use std::path::Path;

use crate::error::{DicomError, Result};
use crate::tags::{implicit_vr, Tag, TransferSyntax, Vr};
use crate::{tags, DICM_MAGIC};

/// The resolved structure of encapsulated pixel data (PS3.5 Annex A.4):
/// the fragment list and the Basic Offset Table's frame offsets.
#[derive(Debug, Clone, PartialEq)]
pub struct PixelFragments {
    /// The fragment items after the Basic Offset Table, in stream order.
    pub fragments: Vec<Vec<u8>>,
    /// The Basic Offset Table's byte offsets - one per frame into the
    /// concatenated fragment stream. Empty when the table itself was
    /// empty (single-frame objects, and multi-frame objects whose writer
    /// left the table empty).
    pub basic_offset_table: Vec<u32>,
}

/// One decoded dataset element: tag, VR, raw value bytes, and - for a
/// sequence the parser resolved - its parsed items.
#[derive(Debug, Clone)]
pub struct DicomElement {
    /// Tag (group, element).
    pub tag: Tag,
    /// Value representation.
    pub vr: Vr,
    /// Raw value bytes (empty for a sequence whose items were resolved;
    /// for encapsulated pixel data, the concatenated fragments).
    pub value: Vec<u8>,
    /// For encapsulated pixel data, the resolved fragment structure.
    pub pixel_fragments: Option<PixelFragments>,
    /// For a resolved sequence, the parsed items in order (each item is a
    /// flat element list; nested sequences appear as elements with their
    /// own `items`). Empty for non-sequence elements, and for sequences
    /// the dataset encodes with an unknown VR - whose content stays raw
    /// in `value` (implicit VR cannot identify an untabulated sequence).
    pub items: Vec<Vec<DicomElement>>,
}

impl DicomElement {
    /// Decodes a `US` (unsigned short) value.
    pub fn as_us(&self) -> Result<u16> {
        if self.value.len() < 2 {
            return Err(DicomError::BadValue {
                tag: self.tag,
                reason: "US value shorter than 2 bytes".into(),
            });
        }
        Ok(u16::from_le_bytes([self.value[0], self.value[1]]))
    }

    /// Decodes the first number of an `IS` (integer string) value.
    pub fn as_is(&self) -> Result<i64> {
        let s = self.trimmed_ascii();
        s.parse::<i64>().map_err(|e| DicomError::BadValue {
            tag: self.tag,
            reason: format!("IS parse: {e}"),
        })
    }

    /// Decodes the first number of a `DS` (decimal string) value; `DS` may
    /// hold multiple backslash-separated values.
    pub fn as_ds_first(&self) -> Result<f64> {
        let s = self.trimmed_ascii();
        let first = s.split('\\').next().unwrap_or_default();
        first
            .trim()
            .parse::<f64>()
            .map_err(|e| DicomError::BadValue {
                tag: self.tag,
                reason: format!("DS parse: {e}"),
            })
    }

    /// Decodes all values of a `DS` decimal string.
    pub fn as_ds_vec(&self) -> Result<Vec<f64>> {
        let s = self.trimmed_ascii();
        s.split('\\')
            .map(|p| {
                p.trim().parse::<f64>().map_err(|e| DicomError::BadValue {
                    tag: self.tag,
                    reason: format!("DS parse: {e}"),
                })
            })
            .collect()
    }

    /// Decodes a text value (`UI`, `LO`, `CS`, `DA`, …), trimming padding.
    pub fn as_text(&self) -> String {
        self.trimmed_ascii().to_string()
    }

    fn trimmed_ascii(&self) -> &str {
        let end = self
            .value
            .iter()
            .rposition(|&b| b != 0 && !b.is_ascii_whitespace())
            .map_or(0, |i| i + 1);
        core::str::from_utf8(&self.value[..end]).unwrap_or("")
    }
}

/// Parser entry point for DICOM Part-10 files.
#[derive(Debug, Default)]
pub struct DicomParser;

impl DicomParser {
    /// Parses one DICOM Part-10 file from disk.
    ///
    /// Single-frame entry point: a multi-frame object is an
    /// [`DicomError::InconsistentSeries`] naming
    /// [`DicomParser::parse_file_all`], never a silently truncated read.
    pub fn parse_file(path: &Path) -> Result<crate::series::DicomSlice> {
        let slices = Self::parse_file_all(path)?;
        Self::single(&slices)
    }

    /// Parses one DICOM Part-10 file from disk into every frame it holds:
    /// one [`DicomSlice`](crate::series::DicomSlice) per frame, with a
    /// multi-frame object's per-frame functional groups mapped onto the
    /// slice list (PS3.3 C.7.6.6, RFC 0001 v1 item 2).
    pub fn parse_file_all(path: &Path) -> Result<Vec<crate::series::DicomSlice>> {
        let bytes = std::fs::read(path)?;
        Self::parse_bytes_all(&bytes)
    }

    /// Parses one DICOM Part-10 file from memory (single-frame; see
    /// [`Self::parse_file`]).
    pub fn parse_bytes(bytes: &[u8]) -> Result<crate::series::DicomSlice> {
        let slices = Self::parse_bytes_all(bytes)?;
        Self::single(&slices)
    }

    /// Parses one DICOM Part-10 file from memory into every frame it
    /// holds (see [`Self::parse_file_all`]).
    pub fn parse_bytes_all(bytes: &[u8]) -> Result<Vec<crate::series::DicomSlice>> {
        if bytes.len() < 132 || &bytes[128..132] != DICM_MAGIC {
            return Err(DicomError::NotDicom(Path::new("<memory>").to_path_buf()));
        }
        let mut cursor = Cursor {
            data: bytes,
            pos: 132,
            encapsulated: false,
            dataset_ts: TransferSyntax::ExplicitVrLittleEndian,
        };

        // File Meta Information is always explicit VR little endian
        // (PS3.10 §7.1); its last element hands us the transfer syntax.
        // The meta/dataset boundary is found by peeking the next tag's
        // group rather than by reading the element explicitly: an
        // implicit-VR dataset's first element header is a 32-bit length,
        // which an explicit read would misinterpret as a VR and reject.
        let mut transfer_syntax = None;
        while cursor.has_more() && cursor.peek_tag()?.0 == 0x0002 {
            let el = cursor.next_explicit()?;
            if el.tag == tags::TRANSFER_SYNTAX_UID {
                transfer_syntax = Some(
                    TransferSyntax::from_uid(&el.as_text())
                        .map_err(DicomError::UnknownTransferSyntax)?,
                );
            }
        }
        let ts = transfer_syntax.unwrap_or(TransferSyntax::ImplicitVrLittleEndian);
        // Encapsulated syntaxes still carry an explicit-VR-LE dataset; only the
        // pixel data is compressed. Without this, an RLE file would be parsed
        // as implicit VR and every tag after the meta group would be garbage.
        cursor.dataset_ts = ts.dataset_encoding();
        cursor.encapsulated = ts.is_encapsulated();

        let mut builder = crate::series::SliceBuilder::new(ts);
        while cursor.has_more() {
            let el = cursor.next_dataset_element()?;
            builder.absorb(el)?;
        }
        builder.build()
    }

    /// The single-slice view of a parsed frame list.
    fn single(slices: &[crate::series::DicomSlice]) -> Result<crate::series::DicomSlice> {
        match slices.len() {
            1 => Ok(slices[0].clone()),
            n => Err(DicomError::InconsistentSeries(format!(
                "multi-frame object with {n} frames; use parse_bytes_all / parse_file_all"
            ))),
        }
    }
}

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
    /// Whether pixel data arrives encapsulated. Decides whether an
    /// undefined-length PixelData is fragment data to collect or a sequence
    /// to resolve.
    encapsulated: bool,
    /// The dataset's encoding (explicit or implicit VR LE), used for every
    /// element below the file meta group.
    dataset_ts: TransferSyntax,
}

/// Item (FFFE,E000).
const ITEM: Tag = (0xFFFE, 0xE000);
/// Item delimitation item (FFFE,E00D).
const ITEM_DELIM: Tag = (0xFFFE, 0xE00D);
/// Sequence delimitation item (FFFE,E0DD).
const SEQ_DELIM: Tag = (0xFFFE, 0xE0DD);

impl<'a> Cursor<'a> {
    fn has_more(&self) -> bool {
        self.pos < self.data.len()
    }

    /// The tag at the cursor without consuming it - the delimiter check
    /// for the bounded readers below.
    fn peek_tag(&self) -> Result<Tag> {
        if self.pos + 4 > self.data.len() {
            return Err(DicomError::UnexpectedEof {
                offset: self.pos,
                while_reading: "peek",
            });
        }
        let g = u16::from_le_bytes([self.data[self.pos], self.data[self.pos + 1]]);
        let e = u16::from_le_bytes([self.data[self.pos + 2], self.data[self.pos + 3]]);
        Ok((g, e))
    }

    fn skip(&mut self, n: usize, what: &'static str) -> Result<()> {
        if self.pos + n > self.data.len() {
            return Err(DicomError::UnexpectedEof {
                offset: self.pos,
                while_reading: what,
            });
        }
        self.pos += n;
        Ok(())
    }

    fn read_exact(&mut self, n: usize, what: &'static str) -> Result<&'a [u8]> {
        if self.pos + n > self.data.len() {
            return Err(DicomError::UnexpectedEof {
                offset: self.pos,
                while_reading: what,
            });
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    fn u16_le(&mut self, what: &'static str) -> Result<u16> {
        let b = self.read_exact(2, what)?;
        Ok(u16::from_le_bytes([b[0], b[1]]))
    }

    fn u32_le(&mut self, what: &'static str) -> Result<u32> {
        let b = self.read_exact(4, what)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    /// Reads the VR + length header for an explicit element (tag already
    /// consumed) and returns (vr, len). `len == usize::MAX` encodes
    /// undefined length.
    fn read_explicit_header(&mut self, tag: Tag) -> Result<(Vr, usize)> {
        let vr_b = self.read_exact(2, "VR")?;
        let mut vr_arr = [0u8; 2];
        vr_arr.copy_from_slice(vr_b);
        let vr = Vr::from_ascii(vr_arr).ok_or_else(|| DicomError::UnsupportedVr {
            tag,
            vr: String::from_utf8_lossy(&vr_arr).into_owned(),
        })?;
        let len = if vr.uses_long_length() {
            // Explicit long form: VR(2) + reserved(2) + length(4).
            let _reserved = self.read_exact(2, "reserved")?;
            let l = self.u32_le("long length")?;
            if l == 0xFFFF_FFFF {
                usize::MAX
            } else {
                l as usize
            }
        } else {
            self.u16_le("short length")? as usize
        };
        Ok((vr, len))
    }

    /// Reads one element in explicit VR LE (the file meta group's
    /// encoding, and the dataset encoding of every syntax but implicit).
    fn next_explicit(&mut self) -> Result<DicomElement> {
        let g = self.u16_le("tag group")?;
        let e = self.u16_le("tag element")?;
        let tag = (g, e);
        let (vr, len) = self.read_explicit_header(tag)?;
        self.finish_element(tag, vr, len)
    }

    /// Reads one element in the dataset's transfer syntax.
    fn next_dataset_element(&mut self) -> Result<DicomElement> {
        match self.dataset_ts {
            TransferSyntax::ExplicitVrLittleEndian => self.next_explicit(),
            TransferSyntax::ImplicitVrLittleEndian => self.next_implicit(),
            _ => unreachable!("dataset_ts narrows to explicit or implicit VR LE"),
        }
    }

    /// Reads one element in implicit VR LE.
    fn next_implicit(&mut self) -> Result<DicomElement> {
        let g = self.u16_le("tag group")?;
        let e = self.u16_le("tag element")?;
        let tag = (g, e);
        let len = self.u32_le("length")?;
        let vr = implicit_vr(tag);
        let len = if len == 0xFFFF_FFFF {
            usize::MAX
        } else {
            len as usize
        };
        self.finish_element(tag, vr, len)
    }

    /// Reads an element's value (and, for a sequence, its items), handling
    /// native, encapsulated and both sequence length forms.
    ///
    /// An undefined-length value on a compressed syntax is the encapsulated
    /// pixel data: a run of (FFFE,E000) item fragments terminated by a
    /// sequence delimiter. Those fragments are returned so the caller can
    /// decode them. Any other undefined length is a sequence (PS3.5
    /// reserves undefined length for SQ in implicit VR): its items are
    /// resolved rather than skipped, because the multi-frame functional
    /// groups carry the per-frame geometry. A defined-length SQ is parsed
    /// from the bytes just as its undefined-length sibling would be, so
    /// the two encodings of the same sequence cannot drift apart.
    fn finish_element(&mut self, tag: Tag, vr: Vr, len: usize) -> Result<DicomElement> {
        if len == usize::MAX {
            if tag == tags::PIXEL_DATA && self.encapsulated {
                let (value, fragments, basic_offset_table) = self.read_fragments()?;
                return Ok(DicomElement {
                    tag,
                    vr,
                    value,
                    items: Vec::new(),
                    pixel_fragments: Some(PixelFragments {
                        fragments,
                        basic_offset_table,
                    }),
                });
            }
            let items = self.read_items()?;
            return Ok(DicomElement {
                tag,
                vr,
                value: Vec::new(),
                items,
                pixel_fragments: None,
            });
        }
        let value = self.read_exact(len, "element value")?.to_vec();
        let items = if vr == Vr::Sq {
            let mut sub = Cursor {
                data: &value,
                pos: 0,
                encapsulated: false,
                dataset_ts: self.dataset_ts,
            };
            sub.read_items()?
        } else {
            Vec::new()
        };
        Ok(DicomElement {
            tag,
            vr,
            value,
            items,
            pixel_fragments: None,
        })
    }

    /// Reads the elements up to the next item/sequence delimiter or the
    /// end of the enclosing defined-length scope. Delimiters are left for
    /// the caller (they carry scope, not data).
    fn read_elements(&mut self) -> Result<Vec<DicomElement>> {
        let mut out = Vec::new();
        while self.has_more() {
            let tag = self.peek_tag()?;
            if tag == ITEM || tag == ITEM_DELIM || tag == SEQ_DELIM {
                break;
            }
            out.push(self.next_dataset_element()?);
        }
        Ok(out)
    }

    /// Reads the items of one sequence until its sequence delimiter (or
    /// the end of the enclosing defined-length scope, for a
    /// defined-length SQ parsed through a sub-cursor).
    fn read_items(&mut self) -> Result<Vec<Vec<DicomElement>>> {
        let mut items = Vec::new();
        loop {
            if !self.has_more() {
                break;
            }
            let tag = self.peek_tag()?;
            if tag == SEQ_DELIM {
                self.skip(8, "sequence delimiter")?;
                break;
            }
            if tag != ITEM {
                return Err(DicomError::BadValue {
                    tag,
                    reason: "expected an item in a sequence".into(),
                });
            }
            self.skip(4, "item tag")?;
            let len = self.u32_le("item length")?;
            if len == 0xFFFF_FFFF {
                let content = self.read_elements()?;
                if self.peek_tag()? != ITEM_DELIM {
                    let t = self.peek_tag()?;
                    return Err(DicomError::BadValue {
                        tag: t,
                        reason: "undefined-length item without an item delimiter".into(),
                    });
                }
                self.skip(8, "item delimiter")?;
                items.push(content);
            } else {
                let bytes = self.read_exact(len as usize, "item content")?;
                let mut sub = Cursor {
                    data: bytes,
                    pos: 0,
                    encapsulated: false,
                    dataset_ts: self.dataset_ts,
                };
                items.push(sub.read_elements()?);
            }
        }
        Ok(items)
    }

    /// Collects encapsulated pixel data fragments up to the sequence delimiter.
    ///
    /// PS3.5 Annex A.4: the first item is the Basic Offset Table, whose
    /// value is one u32 byte offset per frame into the concatenated
    /// fragment stream. It is resolved alongside the fragments (a
    /// multi-frame object's frame boundaries live there); the
    /// concatenated stream is returned first, as the element's `value`.
    fn read_fragments(&mut self) -> Result<(Vec<u8>, Vec<Vec<u8>>, Vec<u32>)> {
        let mut fragments: Vec<Vec<u8>> = Vec::new();
        let mut basic_offset_table: Vec<u32> = Vec::new();
        let mut offset_table_seen = false;
        loop {
            if !self.has_more() {
                return Err(DicomError::UnexpectedEof {
                    offset: self.pos,
                    while_reading: "encapsulated pixel data",
                });
            }
            let g = self.u16_le("fragment tag group")?;
            let e = self.u16_le("fragment tag element")?;
            match (g, e) {
                ITEM => {
                    let len = self.u32_le("fragment length")?;
                    if len == 0xFFFF_FFFF {
                        // A fragment must have a defined length; an undefined
                        // one means the stream is malformed rather than
                        // multi-frame, and cannot be skipped safely.
                        return Err(DicomError::BadValue {
                            tag: tags::PIXEL_DATA,
                            reason: "encapsulated fragment has undefined length".into(),
                        });
                    }
                    let data = self.read_exact(len as usize, "fragment data")?.to_vec();
                    if !offset_table_seen {
                        offset_table_seen = true;
                        // The Basic Offset Table is a separate item whose
                        // value is u32 frame offsets, not pixel data.
                        if data.len() % 4 != 0 {
                            return Err(DicomError::BadValue {
                                tag: tags::PIXEL_DATA,
                                reason: "Basic Offset Table is not a multiple of four bytes".into(),
                            });
                        }
                        basic_offset_table = data
                            .chunks_exact(4)
                            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                            .collect();
                    } else {
                        fragments.push(data);
                    }
                }
                SEQ_DELIM => {
                    let _len = self.u32_le("sequence delimiter length")?;
                    return Ok((fragments.concat(), fragments, basic_offset_table));
                }
                ITEM_DELIM => {
                    let _len = self.u32_le("item delimiter length")?;
                }
                other => {
                    return Err(DicomError::BadValue {
                        tag: other,
                        reason: format!(
                            "unexpected ({:04x},{:04x}) in encapsulated pixel data",
                            other.0, other.1
                        ),
                    });
                }
            }
        }
    }
}

/// Serialises an element back to explicit-VR bytes (used by the synthetic
/// writer and round-trip tests).
pub fn encode_element_explicit(tag: Tag, vr: Vr, value: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(12 + value.len());
    out.extend_from_slice(&tag.0.to_le_bytes());
    out.extend_from_slice(&tag.1.to_le_bytes());
    out.extend_from_slice(vr.ascii());
    if vr.uses_long_length() {
        out.extend_from_slice(&[0, 0]); // reserved
        out.extend_from_slice(&(value.len() as u32).to_le_bytes());
    } else {
        out.extend_from_slice(&(value.len() as u16).to_le_bytes());
    }
    out.extend_from_slice(value);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tags::{COLUMNS, MODALITY, ROWS};

    fn explicit_us(tag: Tag, v: u16) -> Vec<u8> {
        encode_element_explicit(tag, Vr::Us, &v.to_le_bytes())
    }

    fn part10(dataset: &[u8]) -> Vec<u8> {
        let mut buf = vec![0u8; 128];
        buf.extend_from_slice(b"DICM");
        // Meta group declaring explicit VR LE (as real files carry).
        buf.extend_from_slice(&encode_element_explicit(
            crate::tags::TRANSFER_SYNTAX_UID,
            Vr::Ui,
            b"1.2.840.10008.1.2.1 ",
        ));
        buf.extend_from_slice(dataset);
        buf
    }

    #[test]
    fn parse_explicit_minimal() {
        let mut dataset = Vec::new();
        dataset.extend_from_slice(&explicit_us(ROWS, 4));
        dataset.extend_from_slice(&explicit_us(COLUMNS, 6));
        dataset.extend_from_slice(&encode_element_explicit(MODALITY, Vr::Cs, b"CT "));
        let buf = part10(&dataset);
        let slice = DicomParser::parse_bytes(&buf).expect("parses");
        assert_eq!(slice.rows, 4);
        assert_eq!(slice.columns, 6);
        assert!(slice.pixel_data.is_empty());
        assert_eq!(slice.modality.as_deref(), Some("CT"));
    }

    #[test]
    fn rejects_non_dicom() {
        let junk = b"not a dicom file at all....";
        assert!(matches!(
            DicomParser::parse_bytes(junk),
            Err(DicomError::NotDicom(_))
        ));
    }

    #[test]
    fn skips_undefined_length_sequence() {
        // Explicit dataset containing an undefined-length SQ we do not
        // care about, followed by a usable element.
        let mut sq = encode_element_explicit((0x0008, 0x1110), Vr::Sq, &[]);
        // (the SQ itself was encoded with undefined length by hand below)
        let _ = &mut sq;
        let mut dataset = Vec::new();
        dataset.extend_from_slice(&{
            let mut v = Vec::new();
            v.extend_from_slice(&0x0008u16.to_le_bytes());
            v.extend_from_slice(&0x1110u16.to_le_bytes());
            v.extend_from_slice(b"SQ");
            v.extend_from_slice(&[0, 0]);
            v.extend_from_slice(&0xFFFF_FFFFu32.to_le_bytes());
            // Item with defined length containing one US element
            v.extend_from_slice(&0xFFFEu16.to_le_bytes());
            v.extend_from_slice(&0xE000u16.to_le_bytes());
            v.extend_from_slice(&(explicit_us(ROWS, 1).len() as u32).to_le_bytes());
            v.extend_from_slice(&explicit_us(ROWS, 1));
            // Sequence delimiter
            v.extend_from_slice(&0xFFFEu16.to_le_bytes());
            v.extend_from_slice(&0xE0DDu16.to_le_bytes());
            v.extend_from_slice(&0u32.to_le_bytes());
            v
        });
        dataset.extend_from_slice(&explicit_us(ROWS, 9));
        let buf = part10(&dataset);
        let slice = DicomParser::parse_bytes(&buf).expect("parses");
        assert_eq!(slice.rows, 9);
    }
}
