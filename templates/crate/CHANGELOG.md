# Changelog

All notable changes to `tpt-med-example` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Crates in the `tpt-medical` workspace version and release **independently** on
a six-week cadence. Workspace-level and cross-cutting changes are recorded in
the [workspace CHANGELOG](../../../CHANGELOG.md); this file records only what
changes for consumers of this crate.

## [Unreleased]

### Planned
- TODO: what is coming, and which of those changes are semver-breaking.

### Notes
- Record here which changes are **breaking even though they look additive** —
  adding a variant to a public enum is one, because it breaks downstream
  exhaustive matches. This is the note people need *before* they hit it.

## [0.1.0] - YYYY-MM-DD

### Added
- Initial release: what this crate actually does, in enough detail that
  someone can tell whether it is what they need.
- Cite the source of every published constant, here or in the code.

### Verification
- What is verified and against what reference. For numerical code, name the
  analytical solution, published correlation, or standard — and say explicitly
  that a regression-anchored value on a synthetic phantom is not a clinical
  claim.

### Known limitations
- Mirror the README's Known Limitations. A changelog that lists only
  additions tells a reader nothing about what will break.

[Unreleased]: https://github.com/tpt-solutions/tpt-medical/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/tpt-solutions/tpt-medical/releases/tag/v0.1.0
