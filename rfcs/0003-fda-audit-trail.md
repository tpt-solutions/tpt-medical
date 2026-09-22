# RFC 0003: FDA Audit Trail

- **Status:** Accepted
- **Started:** 2026-09-20
- **Crates:** `tpt-med-audit`, `tpt-med-fda`

## Summary

Design the 21 CFR Part 11 audit subsystem: hash-chained event records,
electronic signatures with meaning, and a canonical export package, built on
in-repo cryptographic primitives.

## Motivation

Simulation outputs feeding a 510(k)/De Novo submission must be reproducible
*and* attributable. Regulated deployments need to demonstrate: who did what,
when, why, and that the record was not altered after the fact. Existing Rust
options either pull heavyweight ring/openssl dependencies (awkward for WASM
and overkill for the threat model) or leave chaining policy to each
application — which is exactly where subtle, unauditable mistakes happen.

## Design

### Primitives (`tpt-med-audit`)

- **SHA-256** (FIPS 180-4), pure Rust, locked by the standard's test
  vectors (empty/`abc`/multi-block/one-million-`a`).
- **HMAC-SHA256** (RFC 2104), locked by RFC 4231 vectors including the
  >block-size key case.
- Constant-time tag comparison.
- Hash chain: `digest[i] = H(run_id ‖ digest[i−1] ‖ canonical(entry))` —
  retroactive edits break every later link; a different run id yields a
  different chain, so logs from different runs cannot be spliced.

### Trail semantics (`tpt-med-fda`)

- Every entry: sequence number, UTC timestamp (ISO-8601; civil-date
  conversion in-crate), actor token, object class/token, action,
  reason-for-change (required for Modify per §11.10(e)).
- **Actor tokens, not identities**: identity mapping is institution policy;
  the trail stores non-identifying tokens so the log itself is not a PHI
  store (see `tpt-med-core`).
- Electronic signatures (§11.50): signer, timestamp, meaning
  (author/reviewer/approver/responsible-party); re-signing after new
  entries is recorded explicitly.
- Export: canonical JSON (hand-rolled serializer with strict escaping —
  parsed back in tests) plus a detached HMAC-SHA256 tag under an operator
  key; tag binds the full payload including the digest chain.

### Verification story

`verify_integrity()` re-walks the chain; `export_tag` re-derives the HMAC.
Tampering with any entry, the digest index, or the exported bytes is
detectable. Tests demonstrate retroactive-edit, truncation, seed-splice, and
tag-mismatch detection.

### Explicit limits (documented, not hidden)

- HMAC-SHA256 with symmetric operator keys provides **integrity and
  attribution of the record set**, not non-repudiation against a malicious
  key holder. Deployments needing non-repudiation anchor the export tag
  into an external trust root (HSM countersignature, RFC 3161 timestamp);
  the package format reserves a `countersignatures` field for this.
- Part 11 compliance is a property of the *validated system* (procedures,
  training, retention), not of a library.

## Alternatives considered

- **ed25519 signatures** (via `ed25519-dalek`, dual MIT/Apache): the
  upgrade path for asymmetric signing; deferred until the countersignature
  anchor lands since it forces a dependency into the default build.
- **Merkle tree instead of a chain**: better for partial-disclosure proofs;
  the chain is simpler and sufficient for append-only audit; revisit if
  third-party partial verification becomes a requirement.

## Unresolved questions

- On-the-wire format for third-party verification (JSONL vs CBOR).
- Retention/snapshot policy hooks (nightly sealed checkpoints?).
