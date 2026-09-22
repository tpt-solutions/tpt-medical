# Audit trails (21 CFR Part 11)

`tpt-med-fda` records every operator action as a chained entry
(`digest[i] = SHA256(run_id ‖ digest[i−1] ‖ entry)`), supports electronic
signatures with meaning (§11.50), and exports a canonical JSON package with
a detached HMAC-SHA256 tag.

- Actor tokens are non-identifying (`tpt-med-core` privacy model) — the log
  is not a PHI store.
- `verify_integrity()` detects retroactive edits, truncation, and
  seed-splicing; the tag binds the entire payload.
- Limits (rfcs/0003): symmetric keys give integrity + attribution, not
  non-repudiation; anchor the tag externally (HSM/RFC 3161) when required.

The `fda-package` example runs a full screening analysis under audit and
exports the signed package to `test-data/golden/regulatory/`.
