//! SHA-256 hash chains and HMAC-SHA256 signatures for tamper-evident
//! audit logs.
//!
//! Pure-Rust, dependency-free implementations of FIPS 180-4 SHA-256 and
//! RFC 2104 HMAC, plus the audit-chain construction used by
//! `tpt-med-fda`: each record commits to the previous record's digest,
//! so any retroactive edit breaks every subsequent link.
//!
//! # Examples
//!
//! ```
//! use tpt_med_audit::{hash_chain, sha256_hex, hmac_sha256, verify_chain};
//!
//! let entries = vec!["create:sim-1".to_string(), "modify:sim-1".to_string()];
//! let chain = hash_chain(b"run-1", &entries);
//! assert_eq!(chain.len(), 2);
//! assert!(verify_chain(b"run-1", &entries, &chain));
//! let key = [7u8; 32];
//! let tag = hmac_sha256(&key, b"export package");
//! assert_eq!(tag.len(), 32);
//! assert!(!sha256_hex(b"").is_empty());
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

pub mod anchor;
#[cfg(feature = "ed25519")]
pub mod ed25519;

pub use anchor::{AnchorKind, AnchorRecord};

/// SHA-256 digest of `data` (FIPS 180-4).
pub fn sha256(data: &[u8]) -> [u8; 32] {
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];

    // Padding: message ‖ 0x80 ‖ zeros ‖ bit-length(u64 BE)
    let bit_len = (data.len() as u64) * 8;
    let mut msg = data.to_vec();
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bit_len.to_be_bytes());

    let mut w = [0u32; 64];
    for block in msg.chunks_exact(64) {
        for (i, word) in w.iter_mut().take(16).enumerate() {
            *word = u32::from_be_bytes([
                block[4 * i],
                block[4 * i + 1],
                block[4 * i + 2],
                block[4 * i + 3],
            ]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }

        let (mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh) =
            (h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }

    let mut out = [0u8; 32];
    for (i, word) in h.iter().enumerate() {
        out[4 * i..4 * i + 4].copy_from_slice(&word.to_be_bytes());
    }
    out
}

/// SHA-256 as lowercase hex.
pub fn sha256_hex(data: &[u8]) -> String {
    hex(&sha256(data))
}

/// Hex-encodes bytes.
pub fn hex(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

/// HMAC-SHA256 per RFC 2104.
pub fn hmac_sha256(key: &[u8], message: &[u8]) -> [u8; 32] {
    const BLOCK: usize = 64;
    let mut k = [0u8; BLOCK];
    if key.len() > BLOCK {
        k[..32].copy_from_slice(&sha256(key));
    } else {
        k[..key.len()].copy_from_slice(key);
    }
    let ipad: Vec<u8> = k.iter().map(|b| b ^ 0x36).collect();
    let opad: Vec<u8> = k.iter().map(|b| b ^ 0x5c).collect();
    let mut inner = ipad;
    inner.extend_from_slice(message);
    let inner_hash = sha256(&inner);
    let mut outer = opad;
    outer.extend_from_slice(&inner_hash);
    sha256(&outer)
}

/// A hash chain over ordered entries: `digest[i] = H(digest[i−1] ‖ entry)`,
/// seeded with a chain-specific seed (e.g. a run id). Any retroactive edit
/// to entry `i` invalidates digest `i` and every later link.
pub fn hash_chain(seed: &[u8], entries: &[String]) -> Vec<String> {
    let mut prev = sha256(seed);
    let mut digests = Vec::with_capacity(entries.len());
    for entry in entries {
        let mut buf = prev.to_vec();
        buf.extend_from_slice(entry.as_bytes());
        prev = sha256(&buf);
        digests.push(hex(&prev));
    }
    digests
}

/// Verifies a chain against the entries it claims to cover.
pub fn verify_chain(seed: &[u8], entries: &[String], digests: &[String]) -> bool {
    if entries.len() != digests.len() {
        return false;
    }
    hash_chain(seed, entries)
        .iter()
        .zip(digests)
        .all(|(a, b)| a == b)
}

/// Per-index verdict for one chain link, from [`verify_chain_detailed`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinkStatus {
    /// The stored digest equals the recomputed digest. Note: because the
    /// chain commits to the previous link, every link from the first
    /// tampered entry onward is also reported [`LinkStatus::Broken`] — a
    /// `Valid` link after a `Broken` one is impossible unless the seed (or
    /// the whole tail) was replaced.
    Valid,
    /// The stored digest does not match the recomputation. The first
    /// `Broken` index is the tampered entry.
    Broken,
    /// No stored digest exists at this index (truncated log).
    Missing,
}

/// One link's forensic report: index, verdict, and the digests compared.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinkReport {
    /// Zero-based entry index.
    pub index: usize,
    /// Verdict for this link.
    pub status: LinkStatus,
    /// Recomputed digest (what the chain requires here).
    pub expected: String,
    /// Stored digest, when one exists at this index.
    pub stored: Option<String>,
}

/// Walks the chain and reports each link individually, so a forensic tool
/// can show *where* the log was broken rather than only that it fails.
///
/// Length mismatches are reported as [`LinkStatus::Missing`] links beyond
/// the shorter side; the first `Broken` index is the first tampered entry
/// (everything after it recomputes against a tampered prefix and therefore
/// also reports `Broken`, unless the tail was regenerated wholesale).
pub fn verify_chain_detailed(
    seed: &[u8],
    entries: &[String],
    digests: &[String],
) -> Vec<LinkReport> {
    let mut prev = sha256(seed);
    let mut reports = Vec::with_capacity(entries.len().max(digests.len()));
    for i in 0..entries.len().max(digests.len()) {
        let entry = entries.get(i);
        let mut buf = prev.to_vec();
        buf.extend_from_slice(entry.map_or(&[][..], |e| e.as_bytes()));
        prev = sha256(&buf);
        let expected = hex(&prev);
        let stored = digests.get(i).cloned();
        let status = match stored.as_deref() {
            None => LinkStatus::Missing,
            Some(s) if s == expected => LinkStatus::Valid,
            Some(_) => LinkStatus::Broken,
        };
        reports.push(LinkReport {
            index: i,
            status,
            expected,
            stored,
        });
    }
    reports
}

/// Constant-time equality for tag comparison (avoids timing oracles).
pub fn ct_eq(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut diff = 0u8;
    for (x, y) in a.iter().zip(b) {
        diff |= x ^ y;
    }
    diff == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fips_180_4_vectors() {
        assert_eq!(
            sha256_hex(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(
            sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
        // One million 'a' — exercises multi-block + length handling.
        let million = vec![b'a'; 1_000_000];
        assert_eq!(
            sha256_hex(&million),
            "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"
        );
    }

    #[test]
    fn rfc_4231_hmac_vector_case_1() {
        // RFC 4231 test case 1 (HMAC-SHA-256)
        let key = [0x0b_u8; 20];
        let tag = hmac_sha256(&key, b"Hi There");
        assert_eq!(
            hex(&tag),
            "b0344c61d8db38535ca8afceaf0bf12b881dc200c9833da726e9376c2e32cff7"
        );
    }

    #[test]
    fn rfc_4231_hmac_vector_long_key() {
        // RFC 4231 test case 6 — key larger than the block size (131 bytes
        // of 0xaa), exercising the hash-the-key-first path.
        let key = vec![0xaa_u8; 131];
        let tag = hmac_sha256(
            &key,
            b"Test Using Larger Than Block-Size Key - Hash Key First",
        );
        assert_eq!(
            hex(&tag),
            "60e431591ee0b67f0d8a26aacbf5b77f8e0bc6213728c5140546040f0ee37f54"
        );
    }

    #[test]
    fn chain_detects_tampering() {
        let entries = vec![
            "create:case-1".into(),
            "modify:case-1".into(),
            "simulate:case-1".into(),
        ];
        let seed = b"run-2026-09-20";
        let digests = hash_chain(seed, &entries);
        assert!(verify_chain(seed, &entries, &digests));

        // Retroactive edit of entry 0 breaks every later link.
        let mut forged = entries.clone();
        forged[0] = "create:case-2".into();
        assert!(!verify_chain(seed, &forged, &digests));

        // Truncation is detected.
        assert!(!verify_chain(seed, &entries[..2], &digests));

        // A different seed produces a different chain.
        assert_ne!(hash_chain(b"other", &entries), digests);
    }

    #[test]
    fn detailed_report_pins_the_first_broken_link() {
        let seed = b"run-1";
        let entries = vec!["a".into(), "b".into(), "c".into()];
        let digests = hash_chain(seed, &entries);
        assert!(verify_chain(seed, &entries, &digests));

        // Tamper with entry 1: link 0 stays valid, 1 and 2 break (the chain
        // commits to the previous link, so the tail cannot validate).
        let mut forged = entries.clone();
        forged[1] = "tampered".into();
        let report = verify_chain_detailed(seed, &forged, &digests);
        assert_eq!(report[0].status, LinkStatus::Valid);
        assert_eq!(report[1].status, LinkStatus::Broken);
        assert_eq!(report[2].status, LinkStatus::Broken);
        assert_eq!(
            report.iter().position(|r| r.status == LinkStatus::Broken),
            Some(1),
            "first broken index is the tampered entry"
        );
        // The report carries both sides of the comparison.
        assert_ne!(
            report[1].stored.as_deref(),
            Some(report[1].expected.as_str())
        );
    }

    #[test]
    fn detailed_report_flags_truncation_as_missing() {
        let seed = b"run-2";
        let entries = vec!["a".into(), "b".into(), "c".into()];
        let digests = hash_chain(seed, &entries);
        // Truncated digest log: links without a stored digest are Missing.
        let report = verify_chain_detailed(seed, &entries, &digests[..1]);
        assert_eq!(report[0].status, LinkStatus::Valid);
        assert_eq!(report[1].status, LinkStatus::Missing);
        assert_eq!(report[2].status, LinkStatus::Missing);
        assert!(report[1].stored.is_none());
        assert!(!report[1].expected.is_empty());
    }

    #[test]
    fn detailed_report_flags_entry_gap_as_broken() {
        let seed = b"run-4";
        let entries = vec!["a".into(), "b".into(), "c".into()];
        let digests = hash_chain(seed, &entries);
        // A missing *entry* under a present digest is a gap in the record
        // itself — the recomputation cannot match, so it is Broken.
        let report = verify_chain_detailed(seed, &entries[..1], &digests);
        assert_eq!(report[0].status, LinkStatus::Valid);
        assert_eq!(report[1].status, LinkStatus::Broken);
        assert_eq!(report[2].status, LinkStatus::Broken);
    }

    #[test]
    fn detailed_report_on_fully_consistent_chain() {
        let seed = b"run-3";
        let entries = vec!["x".into()];
        let digests = hash_chain(seed, &entries);
        let report = verify_chain_detailed(seed, &entries, &digests);
        assert!(report.iter().all(|r| r.status == LinkStatus::Valid));
        assert_eq!(report[0].stored.as_deref(), Some(digests[0].as_str()));
    }

    #[test]
    fn ct_eq_behaviour() {
        assert!(ct_eq(b"abc", b"abc"));
        assert!(!ct_eq(b"abc", b"abd"));
        assert!(!ct_eq(b"abc", b"ab"));
    }
}
