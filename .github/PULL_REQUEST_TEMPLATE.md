<!-- Maintainer-authored pull request. This project does not accept external
     pull requests — see CONTRIBUTING.md. Every PR here must reference an
     issue or an accepted RFC. -->

## Summary

What does this PR change and why? **Reference the issue** (e.g. `#123`) or
an **accepted RFC** (`rfcs/NNNN-...`) in the first line. A PR without one will
be asked for it before review.

## Change type

- [ ] Bug fix (non-breaking)
- [ ] New feature (non-breaking)
- [ ] Breaking change (API/numerical behavior)
- [ ] Documentation only
- [ ] RFC artefact (`rfcs/`), landing an issue that was accepted
- [ ] Infrastructure / CI only

## Verification & validation

For changes affecting simulation results or regulatory output:

- [ ] Verification test added or updated against an analytical solution or
      published reference (ASME V&V 40) — not against a snapshot of this
      code's own previous output
- [ ] `scripts/diff-golden.sh <base-sha>` run and the drift table pasted into
      this description
- [ ] If a `test-data/golden/` value changed: the numerical change is explained
      above, the covering verification test is named, and the affected crate
      version is bumped (numerical behaviour is API)
- [ ] Crate `README.md` and `CHANGELOG.md` updated (per-crate documentation
      contract in `docs/book/src/crates.md`)

## Checklist

- [ ] `cargo fmt --all` applied
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` clean
- [ ] `cargo test --workspace` passes
- [ ] `bash scripts/check-crate-docs.sh` passes
- [ ] `cargo deny check licenses` passes
- [ ] No `unsafe` introduced (or justified per accepted RFC)
- [ ] No real patient data (PHI) in code, tests, fixtures, or this description
- [ ] Commits DCO-signed (`git commit -s`)

