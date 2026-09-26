# Contributing to tpt-medical

Thank you for helping build open, auditable medical simulation tooling.
Contributions of every size are welcome — bug reports, verification cases,
documentation corrections, and proposals for new constitutive models.

## How contributions work here

**This project does not accept external pull requests.** Work arrives as
GitHub issues, and changes land from the maintainers' own branches. This is a
deliberate governance choice, not a closed-door policy: it keeps the review
load bounded, and it means a change that touches verification evidence or
regulatory output always lands with a maintainer who owns the numerical
consequences.

What this means in practice:

| You want to… | Do this |
|---|---|
| Report a bug | Open an issue with the [bug report template](.github/ISSUE_TEMPLATE/bug_report.md) |
| Ask for a feature | Open an issue with the [feature request template](.github/ISSUE_TEMPLATE/feature_request.md) |
| Propose a significant design change | Open an issue with the [RFC template](.github/ISSUE_TEMPLATE/rfc.md) and discuss it there |
| Contribute a fix or feature | Open the issue first. A maintainer will pick it up; you can offer a patch in the issue thread |
| Improve the docs | Open an issue describing what's wrong or missing |

Please do **not** open a pull request from a fork. One will be closed with a
pointer to the issue tracker. If you have already written a patch, attaching it
to the issue (a diff in a comment or a gist) means the work is not lost and a
maintainer can land it.

## Why issues-only

Most of the risk in this codebase is not in compiling — it is in a numerical
result changing silently. A PR diff makes it easy to approve a change to a
constitutive model or a correlation constant without noticing the stress field
moved. Routing contributions through an issue forces the question "what is the
numerical consequence, and which verification test covers it?" to be answered
in prose *before* any code exists, which is exactly the discipline ASME V&V 40
asks for.

## The RFC process

New constitutive models, solver algorithms, and regulatory features require an
RFC (RFC 0001–0005 are the accepted precedent).

1. **Open an issue** using the [RFC template](.github/ISSUE_TEMPLATE/rfc.md).
   The issue carries the same sections a committed RFC does: summary,
   motivation, detailed design, drawbacks and alternatives, unresolved
   questions.
2. **Discuss in the issue thread.** A maintainer will engage; a significant
   design change needs two approvals before it is accepted.
3. **A maintainer commits the RFC file** to `rfcs/NNNN-short-name.md` and sets
   its status to `Accepted`, `Deferred` or `Rejected`.
4. **Implementation references the accepted RFC** in the PR description and the
   crate `CHANGELOG.md`.

An RFC issue that is accepted is not a commitment to implement it; it is a
commitment that the design has been thought through.

## Code standards

These apply to every change, whoever authors it.

- **No `unsafe`** without an RFC and a `# Safety` contract review. The whole
  workspace is currently `#![forbid(unsafe_code)]`.
- **No panics in library code paths** for invalid input; return typed errors.
- **No PHI in tests, fixtures, or issues.** `test-data/` contains only
  synthetic datasets. Issues containing real patient data are closed and
  reported.
- **Numerical constants must cite their source** — the design spec, a published
  correlation, or a standard — in a doc comment.
- **Numerical code ships verification** against an analytical solution or
  published reference data, not against a stored snapshot of its own previous
  output.
- **Every crate has a README and a CHANGELOG.** The
  [per-crate documentation contract](docs/book/src/crates.md) lists the
  required README sections, and `scripts/check-crate-docs.sh` enforces it in
  CI.
- **Golden values are change-controlled.** A change to any value under
  `test-data/golden/` must explain the numerical change, name the verification
  test that covers it, and bump the relevant crate version, because numerical
  behaviour is part of the API. `scripts/diff-golden.sh` prints the drift
  table for a PR.

## Local gate (for maintainers and reviewers)

Anyone can run the gate locally before opening an issue, so a report can
include a green run:

```console
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo deny check licenses
bash scripts/check-crate-docs.sh
```

## Contribution licensing (CLA-free via DCO)

This project is CLA-free and uses the
[Developer Certificate of Origin](https://developercertificate.org). By opening
an issue that proposes a change, you agree that your proposal is usable under
the dual licence `MIT OR Apache-2.0`, and you certify the DCO with the
`Signed-off-by:` trailer on any commit you author.

## Release Cadence

SemVer, 6-week cadence. Crates publish independently; see the root
[CHANGELOG.md](CHANGELOG.md) and each crate's own `CHANGELOG.md`.

## Reporting Bugs

Use the [bug report template](.github/ISSUE_TEMPLATE/bug_report.md). Include a
minimal reproducer using synthetic data only. If a bug affects a numerical
result, report the reference or analytical value the result is compared against
and the observed value — "the number changed" is not reproducible, but "peak
von Mises moved from 36.09 MPa to 41.2 MPa" is.

