# tpt-med-example

TODO: one line, lowercase, describing what this crate does.

[![Crates.io](https://img.shields.io/badge/crates.io-tpt--med--example-orange)](https://crates.io/crates/tpt-med-example)
[![Docs.rs](https://img.shields.io/badge/docs.rs-tpt--med--example-blue)](https://docs.rs/tpt-med-example)

| | |
|---|---|
| **Layer** | TODO: `core` / `imaging` / `solid` / `fluid` / `devices` / `surgical` / `regulatory` |
| **Status** | Alpha, `0.1.0` |
| **License** | MIT OR Apache-2.0 |
| **MSRV** | 1.82 |
| **Dependencies** | TODO: list them, or say "none" deliberately |
| **Changelog** | [CHANGELOG.md](CHANGELOG.md) |

---

## Why

The problem this crate solves, and the design decision that resolves it. State
the alternatives and why they were rejected — a reader deciding whether to
depend on this crate needs to know what it is *not*.

If the honest answer is "this exists because crate X needed a Y", say that.

## Features

- What actually ships, as a list. Not aspirations.
- What it deliberately does not do belongs in Known Limitations, not here.

## Conventions

- **Units.** Name the canonical unit for every quantity, and the workspace
  convention it follows.
- **Signs and indexing** where a caller could reasonably get it wrong.
- Anything where the obvious usage is subtly incorrect.

## Usage

```rust
// A compiling, runnable example using the real API. The crate root doc
// comment is the other doctest; at least one should exist.
use tpt_med_example::Example;

fn main() {
    let e = Example::new(2.5);
    assert_eq!(e.value, 2.5);
}
```

## API Overview

| Item | Purpose |
|---|---|
| `Example` | TODO |
| `Example::new(f64)` | TODO |

## Verification

What is verified, against what reference, and what that check *proves*. Not
"we have tests".

For numerical code, name the analytical solution, published correlation, or
standard. If a value is regression-anchored to a synthetic phantom, say that
too, and name the golden dataset under `test-data/golden/`.

## Known Limitations

Required, and not optional. A README that only lists strengths is marketing,
not documentation.

- The thing a user is most likely to try that this crate does not support.
- Any approximation that is only valid in a stated regime.
- What it borrows from a dependency, and what that dependency's limits imply.

## Related Crates

- The crates a user is most likely to reach for instead, or alongside.

## Contributing

See the workspace [CONTRIBUTING.md](../../../CONTRIBUTING.md). This project
does not accept external pull requests — open an issue first. New
constitutive models, solver algorithms, and regulatory features require an
accepted RFC.

## License

Licensed under either of [MIT](../../../LICENSE-MIT) or
[Apache-2.0](../../../LICENSE-APACHE), at your option.

## Regulatory Disclaimer

Research and development use only. Not cleared or approved by the FDA or any
other regulatory body for clinical diagnostic or treatment use.