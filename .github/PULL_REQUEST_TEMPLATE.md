<!-- Thank you for contributing! All contributions are DCO-signed and dual-licensed MIT OR Apache-2.0. -->

## Summary

What does this PR change and why? Reference the issue or RFC
(`rfcs/NNNN-...`) if applicable.

## Change type

- [ ] Bug fix (non-breaking)
- [ ] New feature (non-breaking)
- [ ] Breaking change (API/numerical behavior)
- [ ] Documentation only
- [ ] RFC proposal (`rfcs/`)

## Verification & validation

For changes affecting simulation results or regulatory output:

- [ ] Verification test added/updated against an analytical solution or
      published reference (ASME V&V 40)
- [ ] Golden reference datasets checked (`test-data/golden/`)
- [ ] Numerical changes documented in CHANGELOG.md

## Checklist

- [ ] `cargo fmt --all` applied
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` clean
- [ ] `cargo test --workspace` passes
- [ ] No `unsafe` introduced (or justified per RFC)
- [ ] No real patient data (PHI) in code, tests, or fixtures
- [ ] Commits DCO-signed (`git commit -s`)
