# Changelog

All notable changes to `tpt-med-audit` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this crate
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Added
- Crate README stating plainly that this crate provides **primitives, not an
  audit trail**: the trail, schema, signatures with meaning, timestamps and
  export live in `tpt-med-fda`. It also records that `verify_chain` is
  all-or-nothing and does not report the first failing index, which a forensic
  workflow will want.

### Planned
- A `verify_chain_detailed` returning per-index status, so a forensic tool can
  report the first broken link rather than only that the chain fails.
- Asymmetric signatures (Ed25519) for non-repudiation, behind a clearly named
  feature, alongside the existing symmetric HMAC path.
- External anchoring helpers for RFC 3161 timestamping or a transparency log,
  which RFC 0003 identifies as the way to close the non-repudiation gap.

### Notes
- **Never weaken or replace a cryptographic primitive without an RFC and a
  security review.** Changes here must be behaviour-preserving refactors only.
  In particular, `hash_chain`'s construction is committed to by every existing
  exported `tpt-med-fda` package; changing it invalidates all of them.
- A symmetric key gives integrity and attribution, **not** non-repudiation:
  anyone holding the key can produce a valid tag.
- `sha256` is not suitable for password hashing; there is no key stretching.
- The chain construction is not hardened against a timing side channel (only
  the final comparison is constant-time). This is acceptable for the threat
  model — an attacker who can time log verification can usually read the log —
  but it is stated rather than glossed over.

## [0.1.0] - 2026-09-22

### Added
- `sha256(data) -> [u8; 32]` — FIPS 180-4 SHA-256 with correct padding and
  big-endian length encoding.
- `sha256_hex(data) -> String` — the lowercase hex digest.
- `hmac_sha256(key, message) -> [u8; 32]` — RFC 2104 HMAC, correctly handling
  keys longer than the 64-byte block (hashed down) and keys shorter than it
  (zero-padded).
- `hash_chain(seed, entries) -> Vec<String>` —
  `digest[i] = SHA256(seed ‖ digest[i-1] ‖ canonical(entry[i]))`, so a
  retroactive edit to entry *i* changes every subsequent digest and **truncation
  is caught too**, since the chain length is part of what the verifier
  recomputes.
- `verify_chain(seed, entries, digests) -> bool` — recompute and compare,
  detecting edits, reordering, truncation and seed splicing.
- `ct_eq(a, b) -> bool` — constant-time byte-slice comparison, so verification
  does not leak digest prefixes through timing.
- `hex(bytes) -> String` — lowercase hex encoding.
- Pure Rust, no dependencies, `#![forbid(unsafe_code)]`, compiles to
  `wasm32-unknown-unknown`.

### Verification
Verified against **published standard test vectors**, not against
self-generated output — a hand-rolled hash tested only against itself proves
nothing:
- **FIPS 180-4 SHA-256 vectors** — the canonical test messages (`""`,
  `"abc"`, the 448-bit message, the million-`'a'` message) against their exact
  expected digests.
- **Padding boundary cases** — message lengths 55, 56, 63, 64 and 65 bytes,
  which straddle the 64-byte block and force the extra length block. These are
  precisely the lengths where a hand-written implementation breaks.
- **RFC 4231 HMAC-SHA256 test vectors** — multiple key and message lengths,
  including keys longer than the block size and shorter than it.
- **Chain properties** — a correct chain, plus each of the four attack shapes
  separately: single-entry edit, reordering, truncation, and seed substitution.
- **Empty and single-entry chains** are exercised, since off-by-one errors in
  the `digest[i-1]` seeding only appear at the boundaries.
- `ct_eq` is asserted on equal lengths and on first/last-byte differences, so
  it is not trivially returning `true`.

### Known limitations
- No key management: generation, storage, rotation and revocation are the
  caller's, and in a regulated deployment they are the hard problem.
- No asymmetric signatures, hence no non-repudiation.
- `verify_chain` does not report which entry broke.
- No timestamp authority: the chain is internally consistent by construction,
  but *when* an entry was created is a separate question.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
