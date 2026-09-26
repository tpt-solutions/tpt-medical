//! End-to-end example binaries (see `src/bin`).

/// Per-crate versions resolved from the workspace `Cargo.lock` at build time.
///
/// [`tpt_med_fda::ReproducibilityManifest`](../tpt_med_fda/struct.ReproducibilityManifest.html)
/// wants the *actual* resolved version of every crate that took part in a
/// run, not the version of the binary recording it. Every `tpt-med-*` crate
/// pins its version to the workspace version today
/// (`version.workspace = true`), so [`CRATE_VERSIONS`] is currently uniform —
/// but it is read from `Cargo.lock`, the one place Cargo itself records the
/// resolved version of every crate, so it stays correct the day a crate is
/// released independently and its version diverges from the rest. See
/// `build.rs`.
pub mod crate_versions {
    include!(concat!(env!("OUT_DIR"), "/crate_versions.rs"));

    /// Looks up the resolved version of `name` in `Cargo.lock`.
    ///
    /// Falls back to this binary's own `CARGO_PKG_VERSION` if `name` is not
    /// present (a typo, or a crate added after the last `cargo build`
    /// regenerated the table) — a stale-but-plausible version beats a build
    /// failure in a manifest that is meant to be attached even when a run
    /// later fails.
    pub fn version_of(name: &str) -> &'static str {
        CRATE_VERSIONS
            .iter()
            .find(|(n, _)| *n == name)
            .map(|(_, v)| *v)
            .unwrap_or(env!("CARGO_PKG_VERSION"))
    }
}
