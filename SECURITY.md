# Security Policy

## Supported Versions

| Version | Supported |
|---|---|
| 0.1.x   | ✅ |

## Reporting a Vulnerability

**Do not open a public issue for security vulnerabilities.**

Report privately via GitHub Security Advisories ("Report a vulnerability" on
the repository's Security tab) or email `security@tpt.solutions`.

Include: affected crate(s), version, a minimal reproducer (synthetic data
only), and your assessment of impact — especially for anything that could

- corrupt or falsify an audit trail / electronic signature,
- silently produce incorrect simulation results,
- leak patient data (e.g., via WASM host bindings).

You will receive an acknowledgement within 3 business days and a status
update at least every 7 days until resolution. Coordinated disclosure:
we aim to release fixes within 90 days.

## Scope Notes

- Cryptographic primitives live in `tpt-med-audit` (SHA-256 hash chains,
  HMAC-SHA256 signatures). Treat any weakness in chain verification or
  signature verification as critical.
- The WASM boundary (`tpt-med-wasm`) must never exfiltrate patient data;
  any host function that transmits data is a vulnerability by design.
