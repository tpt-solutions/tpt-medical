# RFC NNNN: Short kebab-case title

- **Status:** Draft
- **Started:** YYYY-MM-DD
- **Crates:** `tpt-med-xxx` (implementation), and downstream consumers

> This template is for an issue in this project. The project does not accept
> external pull requests: open an issue with the [RFC template](
> ../../.github/ISSUE_TEMPLATE/rfc.md) using these sections, discuss it there,
> and a maintainer commits the file to `rfcs/NNNN-short-name.md` once the
> design is agreed. See [CONTRIBUTING.md](../../CONTRIBUTING.md).

## Summary

One paragraph. What is being decided, in terms a reader who has not read the
rest of the document can follow.

## Motivation

Why is this needed? What workflow or standard does it serve?

Name the **question of interest**, and state its **model risk** and **model
influence** under ASME V&V 40. A design that cannot say what question it
answers, and how bad a wrong answer would be, is not ready to be decided.

## Detailed design

The design in detail: data models, algorithms, APIs, numerics.

- Public API surface, and what a caller can and cannot do with it.
- Numerical scheme, including the discretisation and the solver.
- Units, coordinate conventions, and sign conventions.
- Error handling: typed errors, and what happens on invalid input.
- What `#![forbid(unsafe_code)]` means for this design.

### Alternatives considered

What else could solve this, and why it was rejected. An RFC that does not
name a rejected alternative has usually not done the thinking.

### Drawbacks

What this design makes *harder*. Every design has a cost; an RFC that lists
none is not credible.

## Verification strategy

**Required. An RFC for numerical work without a verification reference will
not be accepted.**

- Code verification: what analytical solution or standard test case proves the
  implementation is doing what the equations say.
- Calculation verification: mesh/time convergence, solver quality, and the
  quantified result.
- Validation: the reference data, its provenance, the comparison metric, and
  the agreement level expected.
- Which existing golden dataset under `test-data/golden/` this changes, and the
  drift `scripts/diff-golden.sh` should show.
- What remains **un**verified, stated as plainly as what is verified.

## Unresolved questions

Open items to settle before acceptance, and who settles them.
