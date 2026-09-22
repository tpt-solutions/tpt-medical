//! Byte-level DICOM dataset parser (uncompressed little-endian syntaxes).

use std::path::Path;

use crate::error::{DicomError, Result};
use crate::tags::{implicit_vr, Tag, TransferSyntax, Vr};
use crate::{tags, DICM_MAGIC};

/// One decoded dataset element: tag, VR, and raw value bytes.
#[derive(Debug, Clone)]
pub struct DicomElement {
    /// Tag (group, element).
    pub tag: Tag,
    /// Value representation.
    pub vr: Vr,
    /// Raw value bytes.
    pub value: Vec<u8>,
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
    pub fn parse_file(path: &Path) -> Result<crate::series::DicomSlice> {
        let bytes = std::fs::read(path)?;
        Self::parse_bytes(&bytes)
    }

    /// Parses one DICOM Part-10 file from memory.
    pub fn parse_bytes(bytes: &[u8]) -> Result<crate::series::DicomSlice> {
        if bytes.len() < 132 || &bytes[128..132] != DICM_MAGIC {
            return Err(DicomError::NotDicom(Path::new("<memory>").to_path_buf()));
        }
        let mut cursor = Cursor {
            data: bytes,
            pos: 132,
        };

        // File Meta Information is always explicit VR little endian
        // (PS3.10 §7.1). Its last element hands us the transfer syntax.
        let mut transfer_syntax = None;
        let mut pending_dataset_element = None;
        while cursor.has_more() {
            let el = cursor.next_explicit()?;
            if el.tag.0 != 0x0002 {
                pending_dataset_element = Some(el);
                break;
            }
            if el.tag == tags::TRANSFER_SYNTAX_UID {
                transfer_syntax = Some(
                    TransferSyntax::from_uid(&el.as_text())
                        .map_err(DicomError::UnknownTransferSyntax)?,
                );
            }
        }
        let ts = transfer_syntax.unwrap_or(TransferSyntax::ImplicitVrLittleEndian);

        let mut builder = crate::series::SliceBuilder::default();
        if let Some(el) = pending_dataset_element {
            builder.absorb(el)?;
        }
        while cursor.has_more() {
            let el = cursor.next_element(ts)?;
            builder.absorb(el)?;
        }
        builder.build()
    }
}

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
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

    /// Reads one element in explicit VR LE.
    fn next_explicit(&mut self) -> Result<DicomElement> {
        let g = self.u16_le("tag group")?;
        let e = self.u16_le("tag element")?;
        let tag = (g, e);
        let (vr, len) = self.read_explicit_header(tag)?;
        let value = self.read_value(tag, vr, len)?;
        Ok(DicomElement { tag, vr, value })
    }

    /// Reads one element in the dataset's transfer syntax.
    fn next_element(&mut self, ts: TransferSyntax) -> Result<DicomElement> {
        match ts {
            TransferSyntax::ExplicitVrLittleEndian => self.next_explicit(),
            TransferSyntax::ImplicitVrLittleEndian => {
                let g = self.u16_le("tag group")?;
                let e = self.u16_le("tag element")?;
                let tag = (g, e);
                let len = self.u32_le("length")?;
                let vr = implicit_vr(tag);
                let value = if len == 0xFFFF_FFFF {
                    self.skip_undefined(ts)?;
                    Vec::new()
                } else {
                    self.read_exact(len as usize, "element value")?.to_vec()
                };
                Ok(DicomElement { tag, vr, value })
            }
        }
    }

    fn read_value(&mut self, tag: Tag, _vr: Vr, len: usize) -> Result<Vec<u8>> {
        if len == usize::MAX {
            self.skip_undefined(TransferSyntax::ExplicitVrLittleEndian)?;
            return Ok(Vec::new());
        }
        let _ = tag;
        Ok(self.read_exact(len, "element value")?.to_vec())
    }

    /// Consumes an undefined-length construct (sequence of items or
    /// encapsulated pixel data) up to its matching delimiter, honouring
    /// nesting. Items are encoded in the dataset's transfer syntax.
    fn skip_undefined(&mut self, ts: TransferSyntax) -> Result<()> {
        let mut depth = 1usize;
        while depth > 0 && self.has_more() {
            let g = self.u16_le("item tag group")?;
            let e = self.u16_le("item tag element")?;
            match (g, e) {
                ITEM => {
                    let len = self.u32_le("item length")?;
                    if len == 0xFFFF_FFFF {
                        depth += 1;
                    } else {
                        self.skip(len as usize, "item content")?;
                    }
                }
                ITEM_DELIM | SEQ_DELIM => {
                    let _len = self.u32_le("delimiter length")?;
                    depth -= 1;
                }
                _ => {
                    // Regular element inside an item, encoded per `ts`.
                    match ts {
                        TransferSyntax::ExplicitVrLittleEndian => {
                            let (vr, len) = self.read_explicit_header((g, e))?;
                            let _ = vr;
                            if len == usize::MAX {
                                depth += 1;
                            } else {
                                self.skip(len, "nested element")?;
                            }
                        }
                        TransferSyntax::ImplicitVrLittleEndian => {
                            let len = self.u32_le("nested length")?;
                            self.skip(len as usize, "nested element")?;
                        }
                    }
                }
            }
        }
        Ok(())
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
