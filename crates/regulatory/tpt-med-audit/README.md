# tpt-med-audit

SHA-256 hash chains and HMAC-SHA256 signatures for tamper-evident audit logs —
pure-Rust, dependency-free, FIPS 180-4 and RFC 2104 conformant.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--audit-orange)](https://crates.io/crates/tpt-med-audit)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--audit-blue)](https://docs.rs/tpt-med-audit)

| | |
|---|---|
| **Layer** | `regulatory` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | none (`std` only) |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

An audit trail nobody can trust is worse than no audit trail, because it
creates the *appearance* of traceability. For a 21 CFR Part 11 record, "we
wrote it down" is not the claim; the claim is "you cannot have changed it
without breaking something detectable".

A **hash chain** gives you that. Each record commits to the previous record's
digest:

```text
digest[i] = SHA256(seed ‖ digest[i-1] ‖ canonical(entry[i]))
```

Any retroactive edit to entry *i* changes `digest[i]`, which changes
`digest[i+1]`, and so on to the end. **Truncation is caught too**, because
the chain length is part of what the verifier recomputes. Adding, removing,
reordering and editing are all detectable from a chain alone.

An **HMAC-SHA256 tag** over the payload then binds the chain to a specific
key, so an attacker who can rewrite the whole log still cannot produce a valid
signature without the operator secret.

This crate exists because that is a genuinely small amount of code with a
genuinely large regulatory consequence, and because a WASM-compilable
audit trail should not drag in a crypto crate with a heavy dependency tree.

## Features

- **`sha256` / `sha256_hex`** — FIPS 180-4 SHA-256, correct padding and
  big-endian length encoding.
- **`hmac_sha256`** — RFC 2104 HMAC, including the block-size padding of keys
  longer than 64 bytes and the zero-padding of shorter keys.
- **`hash_chain(seed, entries)`** — build the chain digests.
- **`verify_chain(seed, entries, digests)`** — recompute and compare, detecting
  edits, reordering, truncation and seed splicing.
- **`ct_eq`** — constant-time byte-slice comparison, so `verify_chain` does
  not leak digest prefixes through timing.
- **`hex`** — lowercase hex encoding for digests and tags.
- Zero dependencies, no `unsafe`, compiles to `wasm32-unknown-unknown`, and
  the same code path is covered by native `cargo test`.

## Conventions and Limits — read this

- **This crate is not an audit trail.** It provides the cryptographic
  primitives. The trail itself, the record schema, signatures with *meaning*,
  timestamps and export live in
  [`tpt-med-fda`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/regulatory/tpt-med-fda).
- **A symmetric key gives integrity and attribution, not non-repudiation**
  (RFC 0003). Anyone holding the key can produce a valid tag. When
  non-repudiation is required, anchor the tag externally — an HSM, a
  transparency log, or RFC 3161 timestamping.
- **A hash chain is what detects tampering, not a signature.** Edits,
  reordering, truncation and seed splicing all break the chain; the HMAC binds
  the chain to one key.
- `verify_chain` recomputes and compares; it does not tell you *which* entry
  broke. For a forensic tool you would walk the chain and report the first
  mismatch index.
- `sha256` is **not** suitable for password hashing; there is no key
  stretching here.
- The chain construction is not hardened against a timing side channel (only
  the final comparison is constant-time). Acceptable for this threat model —
  an attacker who can time your log verification can usually read the log —
  but stated rather than glossed over.
## Usage

```rust
use tpt_med_audit::{hash_chain, hex, hmac_sha256, sha256_hex, verify_chain};

fn main() {
    let entries = vec!["create:sim-1".to_string(), "modify:sim-1".to_string()];

    // Build the chain: each entry commits to the previous digest.
    let seed = b"run-1";
    let chain = hash_chain(seed, &entries);
    assert_eq!(chain.len(), 2);
    assert!(verify_chain(seed, &entries, &chain));

    // A retroactive edit breaks every subsequent link.
    let tampered = vec!["create:sim-1".to_string(), "delete:sim-1".to_string()];
    assert!(!verify_chain(seed, &tampered, &chain));

    // Truncation is caught too: the chain length is part of the check.
    assert!(!verify_chain(seed, &entries[..1], &chain));

    // So is splicing a different run's seed in.
    assert!(!verify_chain(b"run-2", &entries, &chain));

    // Sign the package under an operator secret.
    let key = [7u8; 32];
    let tag = hmac_sha256(&key, b"export package");
    assert_eq!(tag.len(), 32);
    assert_eq!(hex(&tag).len(), 64);

    assert!(!sha256_hex(b"").is_empty());
}
```

## API Overview

| Item | Purpose |
|---|---|
| `sha256(data: &[u8]) -> [u8; 32]` | FIPS 180-4 SHA-256 digest |
| `sha256_hex(data: &[u8]) -> String` | Lowercase hex digest |
| `hmac_sha256(key: &[u8], message: &[u8]) -> [u8; 32]` | RFC 2104 HMAC; handles keys shorter and longer than the 64-byte block |
| `hash_chain(seed: &[u8], entries: &[String]) -> Vec<String>` | `digest[i] = H(seed ‖ digest[i-1] ‖ entry[i])`, as hex |
| `verify_chain(seed: &[u8], entries: &[String], digests: &[String]) -> bool` | Recompute and compare; detects edits, reordering, truncation and seed splicing |
| `ct_eq(a: &[u8], b: &[u8]) -> bool` | Constant-time slice comparison |
| `hex(bytes: &[u8]) -> String` | Lowercase hex encoding |

## Verification

Verified against the published standard test vectors, not against
themselves — a hand-rolled hash that is only tested against its own output
proves nothing:

- **FIPS 180-4 SHA-256 vectors** — the canonical test messages
  (`""`, `"abc"`, the 448-bit message, the million-`'a'` message) and their
  exact expected digests.
- **Padding boundary cases** — messages of length 55, 56, 63, 64 and 65
  bytes, which straddle the 64-byte block and force the extra length block.
  These are exactly the lengths where a hand-written implementation breaks.
- **RFC 4231 HMAC-SHA256 test vectors** — multiple key and message lengths,
  including keys longer than the block size (which must be hashed down) and
  keys shorter than it (which must be zero-padded).
- **Chain properties** — verified for a correct chain, and for each of the
  four attack shapes separately: single-entry edit, reordering, truncation,
  and seed substitution.
- **Empty and single-entry chains** are exercised, since off-by-one errors in
  the `digest[i-1]` seeding only appear at the boundaries.
- `ct_eq` is asserted equal-length equality and unequal on first/last-byte
  differences, so it is not trivially returning `true`.

## Known Limitations

- **No key management.** Keys are passed in as byte slices. Generation,
  storage, rotation and revocation are the caller's problem, and in a
  regulated deployment they are the *hard* problem.
- **No asymmetric signatures.** A symmetric HMAC tag cannot demonstrate
  non-repudiation (see the conventions above).
- **`verify_chain` is all-or-nothing** — it does not report the first failing
  index, which a forensic workflow will want.
- **No timestamp authority.** The chain is internally consistent by
  construction; that an entry was created *when* it claims is a separate
  question requiring an external time source.
- **Not constant-time with respect to entry count or length**, only the final
  comparison. The chain construction itself is not hardened against a timing
  side channel, which is acceptable for the threat model (an attacker who can
  time your log verification can usually just read your log) but should be
  stated rather than glossed over.

## Related Crates

- [`tpt-med-fda`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/regulatory/tpt-med-fda) — builds the 21 CFR Part 11 trail on these primitives (RFC 0003).
- [`tpt-med-core`](https://github.com/tpt-solutions/tpt-medical/tree/master/crates/core/tpt-med-core) — the `AuditEvent`s that become chain entries.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md).
**Never weaken or replace a cryptographic primitive without an
[RFC](../../../rfcs) and a security review.** Changes here are
behaviour-preserving refactors only, and every new case must arrive with its
standard test vector.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body. The software controls here support a 21 CFR Part 11
process; a validated system additionally requires procedural controls — SOPs,
operator training, record retention — that no library can provide.

