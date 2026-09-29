//! External anchoring helpers: attaching a trail's digest to an external
//! service (an RFC 3161 timestamp authority or a transparency log) so the
//! non-repudiation guarantee does not rest on keys this system controls.
//!
//! This crate forms the record and checks the *binding* — that the
//! recorded digest is the SHA-256 of the payload being claimed.
//! Cryptographic verification of the anchor **token** (the TSA's
//! SignedData or the log's inclusion proof) belongs to the anchor
//! service's own verifier, and this type keeps the token opaque on
//! purpose: hand-rolling ASN.1 CMS parsing here would be a second
//! hand-rolled crypto surface, not a feature.

/// The kind of external anchor a record carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnchorKind {
    /// An RFC 3161 timestamp token (a signed `TimeStampResp`) from a
    /// timestamp authority.
    Rfc3161,
    /// A transparency-log anchor (e.g. an SCT or inclusion proof).
    TransparencyLog,
}

impl AnchorKind {
    /// Stable key for reports.
    pub fn key(self) -> &'static str {
        match self {
            AnchorKind::Rfc3161 => "rfc3161",
            AnchorKind::TransparencyLog => "transparency_log",
        }
    }
}

/// A record that a payload's digest was submitted to an external anchor.
///
/// Construction goes through [`AnchorRecord::anchor`], which computes and
/// stores the payload's SHA-256 digest; [`AnchorRecord::covers`] recomputes
/// it, so a record whose payload was silently swapped is detected
/// mechanically.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnchorRecord {
    /// What kind of anchor service holds the payload digest.
    pub kind: AnchorKind,
    /// SHA-256 hex digest of the anchored payload.
    pub digest_hex: String,
    /// The anchor service's opaque response (TSA reply bytes, SCT, proof) —
    /// verified by the service's own tooling, not parsed here.
    pub token_base64: String,
}

impl AnchorRecord {
    fn base64(data: &[u8]) -> String {
        // Minimal standard-alphabet base64 (RFC 4648, with padding) —
        // the crate stays dependency-free.
        const ALPHABET: &[u8; 64] =
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
        for chunk in data.chunks(3) {
            let b = [
                chunk[0],
                chunk.get(1).copied().unwrap_or(0),
                chunk.get(2).copied().unwrap_or(0),
            ];
            let n = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
            out.push(ALPHABET[(n >> 18) as usize & 63] as char);
            out.push(ALPHABET[(n >> 12) as usize & 63] as char);
            out.push(if chunk.len() > 1 {
                ALPHABET[(n >> 6) as usize & 63] as char
            } else {
                '='
            });
            out.push(if chunk.len() > 2 {
                ALPHABET[n as usize & 63] as char
            } else {
                '='
            });
        }
        out
    }

    /// Records that `payload`'s digest was anchored: the caller submits the
    /// digest (or the payload) to the anchor service, obtains `token`, and
    /// stores both here.
    pub fn anchor(kind: AnchorKind, payload: &[u8], token: &[u8]) -> Self {
        Self {
            kind,
            digest_hex: crate::sha256_hex(payload),
            token_base64: Self::base64(token),
        }
    }

    /// True when the recorded digest is still the SHA-256 of `payload` —
    /// the mechanical binding check. The token's cryptographic validity is
    /// the anchor service's verifier's job (see the module docs).
    pub fn covers(&self, payload: &[u8]) -> bool {
        self.digest_hex == crate::sha256_hex(payload)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn anchor_record_binds_to_the_payload_digest() {
        let payload = b"audit-chain-head-digest-or-tag";
        let record = AnchorRecord::anchor(AnchorKind::Rfc3161, payload, b"TSA-RESPONSE-BYTES");
        assert_eq!(record.kind.key(), "rfc3161");
        assert_eq!(record.digest_hex, crate::sha256_hex(payload));
        assert!(record.covers(payload));
        assert!(!record.covers(b"tampered payload"));
        assert_eq!(record.kind.key(), AnchorKind::Rfc3161.key());
    }

    #[test]
    fn base64_matches_the_rfc4648_vectors() {
        // RFC 4648 test vectors.
        let enc = |s: &str| AnchorRecord::base64(s.as_bytes());
        assert_eq!(enc(""), "");
        assert_eq!(enc("f"), "Zg==");
        assert_eq!(enc("fo"), "Zm8=");
        assert_eq!(enc("foo"), "Zm9v");
        assert_eq!(enc("foob"), "Zm9vYg==");
        assert_eq!(enc("fooba"), "Zm9vYmE=");
        assert_eq!(enc("foobar"), "Zm9vYmFy");
    }

    #[test]
    fn transparency_log_kind_round_trips() {
        let record = AnchorRecord::anchor(AnchorKind::TransparencyLog, b"tag", b"sct-bytes");
        assert_eq!(record.kind.key(), "transparency_log");
        assert!(record.covers(b"tag"));
    }
}
