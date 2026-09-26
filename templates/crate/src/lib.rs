//! TODO: one line, lowercase, describing what this crate does.
//!
//! The module doc is the first thing rendered on docs.rs and the first thing
//! read in review, so write it for someone who has never heard of this crate.
//!
//! # What belongs in a crate
//!
//! One responsibility, stated in the `Why` section of the README. If the
//! module doc needs the word "and" twice, the crate does two things.
//!
//! # Verification
//!
//! Numerical code in this workspace must be verified against an analytical
//! solution or published reference, **not** against a stored snapshot of its
//! own previous output. Say here what this crate is verified against.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// A documented public item. Every public item needs a doc comment; the
/// workspace lints make this a warning and CI denies warnings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Example {
    /// A field with units, because a bare `f64` in a medical codebase is a
    /// bug waiting to happen. Say the unit in the name or the doc.
    pub value: f64,
}

impl Example {
    /// Constructor. Prefer `const fn` for plain newtypes.
    pub const fn new(value: f64) -> Self {
        Self { value }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Replace with a real test. A test that only asserts a constructor
    /// round-trips is close to worthless; test the property the crate exists
    /// to guarantee.
    #[test]
    fn constructor_round_trips() {
        assert_eq!(Example::new(2.5).value, 2.5);
    }
}
