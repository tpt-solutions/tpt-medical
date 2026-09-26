# templates/

Starting points for new work, based on the milestones in
`examples/src/bin/*.rs`. These are **scaffolding to copy**, not live crates:
nothing here is a workspace member, and `scripts/check-crate-docs.sh`
deliberately skips this directory.

## What is here

| Template | Copy it when… |
|---|---|
| [`crate/`](crate) | Adding a new library crate to the workspace |
| [`example-binary/`](example-binary) | Adding a new end-to-end milestone binary |
| [`rfc.md`](rfc.md) | Proposing a new constitutive model, solver algorithm, or regulatory feature |

## Why copy rather than `cargo new`

`cargo generate` / `cargo-generate` would be the usual answer, and it is a
reasonable future direction. It is not used here yet because the thing that
actually needs templating is not a crate skeleton — it is the *house style*:
the crate README section contract, the CHANGELOG shape, the `#![forbid(unsafe_code)]`
attribute, the workspace lint table, and the rule that numerical code ships
verification against a published reference. Those are a few hundred lines of
convention, and a template that is not checked by CI is a template that rots.

The rule that keeps these honest: `scripts/check-crate-docs.sh` enforces the
same README and CHANGELOG contract on real crates that these templates
promise. If the templates drift from the contract, the real crates are the
ones that fail.

## Using a template

```console
# A new library crate
cp -r templates/crate crates/<layer>/<new-crate>
$EDITOR crates/<layer>/<new-crate>/Cargo.toml     # name, description, keywords
$EDITOR crates/<layer>/<new-crate>/src/lib.rs

# Register it, or the CI member check will fail:
#   1. add the path to [workspace] members in the root Cargo.toml
#   2. add it to [workspace.dependencies] with a path + version
#   3. run: bash scripts/check-crate-docs.sh
```

Renaming is on you: the templates use `tpt-med-example` as a placeholder, and
`keywords`/`categories` must be crate-specific for crates.io to be useful for
discovery (see the note in the root `Cargo.toml`).
