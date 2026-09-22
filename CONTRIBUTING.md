# Contributing to tpt-medical

Thank you for helping build open, auditable medical simulation tooling.
Contributions of all sizes are welcome — bug fixes, constitutive models,
documentation, verification test cases.

## Development Workflow

1. Fork the repository.
2. Create a feature branch: `feature/hyperelastic-model`.
3. Write code **and tests**. Numerical code requires verification against an
   analytical solution or published reference data (see ASME V&V 40).
4. Run the local gate:

   ```console
   cargo fmt --all
   cargo clippy --workspace --all-targets -- -D warnings
   cargo test --workspace
   cargo deny check licenses
   ```

5. Submit a PR with a DCO sign-off (`git commit -s`).

## Contribution Licensing (CLA-free via DCO)

By contributing, you agree your contribution is dual-licensed
`MIT OR Apache-2.0`, and you certify the
[Developer Certificate of Origin](https://developercertificate.org) with each
commit via `Signed-off-by:`.

## Code Standards

- **No `unsafe`** without an RFC and a `# Safety` contract review.
- **No panics in library code paths** for invalid input; return typed errors.
- **No PHI in tests.** All test fixtures use synthetic data
  (`test-data/` contains only synthetic datasets).
- Numerical constants must cite their source (spec §reference, published
  correlation, or standard) in a doc comment.
- New constitutive models, solver algorithms, or regulatory features require
  an [RFC](rfcs) and two approvals.

## RFC Process

Open a PR adding `rfcs/NNNN-short-name.md` (see existing RFCs for format).
Discussion happens on the PR; the RFC lands as `Accepted`, `Rejected`, or
`Deferred`. Implementation PRs reference the accepted RFC.

## Release Cadence

SemVer, 6-week cadence. Crates publish independently; see CHANGELOG.md.

## Reporting Bugs

Use the [bug report template](.github/ISSUE_TEMPLATE/bug_report.md). Include
a minimal reproducer using synthetic data only — issues containing real
patient data will be closed and reported.
