#!/usr/bin/env bash
#
# check-crate-docs.sh - enforce the per-crate documentation and crates.io
# metadata convention across every workspace member.
#
# Checks, for every member of the workspace:
#
#   1. README.md exists and carries the required section set.
#   2. CHANGELOG.md exists, follows Keep a Changelog, and has an [Unreleased]
#      section plus a dated release section.
#   3. `readme` is wired up in the manifest.
#   4. keywords: 1-5 entries, each <= 20 characters, [A-Za-z0-9_-] only.
#   5. categories: at most 5 entries, each a slug from the published crates.io
#      category registry.
#
# Members are discovered with `find` and then cross-checked against the root
# manifest with a substring test, rather than by parsing the `members` array.
# Hand-parsing TOML in awk is fragile; if a manifest is ever added outside the
# root members list this script fails loudly instead of silently skipping it.
#
# The valid-category list below is vendored on purpose: the check must not need
# network access, and it doubles as documentation of what is actually
# registered. Regenerate with `curl -s https://crates.io/api/v1/categories`.
# Note that `medical-science` is NOT a valid category.
#
# Usage: scripts/check-crate-docs.sh   (from the repository root)
# Exit code: 0 = all good, 1 = one or more failures.

set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$PWD"
ROOT_MANIFEST="Cargo.toml"

# --- vendored crates.io category registry ----------------------------------
VALID_CATEGORIES="accessibility
aerospace
algorithms
api-bindings
artificial-intelligence
asynchronous
authentication
automotive
caching
command-line-interface
command-line-utilities
compilers
compression
computer-vision
concurrency
config
cryptography
data-structures
database
database-implementations
date-and-time
development-tools
email
embedded
emulators
encoding
external-ffi-bindings
filesystem
finance
game-development
game-engines
games
graphics
gui
hardware-support
internationalization
localization
mathematics
memory-management
multimedia
network-programming
no-std
os
parser-implementations
parsing
rendering
rust-patterns
science
security
simulation
template-engine
text-editors
text-processing
value-formatting
virtualization
visualization
wasm
web-programming
web-serialization"

# README sections every crate must document.
REQUIRED_SECTIONS=(
  "## Why"
  "## Features"
  "## Conventions"
  "## Usage"
  "## API Overview"
  "## Verification"
  "## Known Limitations"
  "## Related Crates"
  "## Contributing"
  "## License"
)

failures=0
checked=0

fail() {
  printf '  ERROR  %s\n' "$1" >&2
  failures=$((failures + 1))
}

# Read a scalar string field from the [package] table of a manifest.
pkg_field() {
  awk -v key="$2" '
    /^\[package\]/ { in_pkg = 1; next }
    /^\[/          { in_pkg = 0 }
    in_pkg && index($0, key " ") == 1 {
      line = $0
      sub(/\r$/, "", line)
      sub("^" key "[ \t]*=[ \t]*", "", line)
      gsub(/^"|"$/, "", line)
      print line
      exit
    }
  ' "$1"
}
# Read a string array field from [package], one element per line.
pkg_array() {
  awk -v key="$2" '
    /^\[package\]/ { in_pkg = 1; next }
    /^\[/          { in_pkg = 0 }
    in_pkg && index($0, key " ") == 1 {
      line = $0
      sub(/\r$/, "", line)
      sub("^" key "[ \t]*=[ \t]*", "", line)
      gsub(/^\[/, "", line)
      gsub(/\]$/, "", line)
      gsub(/"/, "", line)
      n = split(line, parts, ",")
      for (i = 1; i <= n; i++) {
        gsub(/^[ \t]+|[ \t]+$/, "", parts[i])
        if (parts[i] != "") print parts[i]
      }
      exit
    }
  ' "$1"
}

printf 'Per-crate documentation check (%s)\n\n' "$ROOT"

# Discover every package manifest, excluding the root and any build output.
while IFS= read -r manifest; do
  dir="$(dirname "$manifest")"
  rel="${dir#./}"

  # The root manifest is the workspace, not a member.
  [ "$rel" = "." ] && continue

  # templates/ holds scaffolding to copy, not live crates: they are not
  # workspace members and are excluded from the scan on purpose.
  case "$rel" in
    templates/*) continue ;;
  esac

  # Cross-check membership: the path must be listed in the root members array.
  if ! grep -qF "\"$rel\"" "$ROOT_MANIFEST"; then
    fail "$rel: has a Cargo.toml but is not listed in the root [workspace] members"
    continue
  fi

  name="$(pkg_field "$manifest" name)"
  [ -z "$name" ] && name="(unnamed: $rel)"
  checked=$((checked + 1))
  printf '%s\n' "$name"

  # --- README ---------------------------------------------------------------
  readme="$dir/README.md"
  if [ ! -f "$readme" ]; then
    fail "$name: missing README.md"
  else
    lines="$(wc -l < "$readme" | tr -d ' ')"
    if [ "$lines" -lt 40 ]; then
      fail "$name: README.md is only $lines lines; expected a comprehensive crate README (>= 40)"
    fi
    for section in "${REQUIRED_SECTIONS[@]}"; do
      if ! grep -qF "$section" "$readme"; then
        fail "$name: README.md is missing the '$section' section"
      fi
    done
    if ! grep -qi 'Regulatory Disclaimer' "$readme"; then
      fail "$name: README.md has no regulatory disclaimer"
    fi

    # Sections must appear in the contract order, not merely all be present.
    # A scrambled README is usually a botched edit that also drops content, so
    # order is checked explicitly rather than left to review.
    order=""
    for section in "${REQUIRED_SECTIONS[@]}"; do
      line="$(grep -nF "$section" "$readme" | head -1 | cut -d: -f1)"
      if [ -n "$line" ]; then
        order="$order $line"
      fi
    done
    sorted="$(printf '%s\n' $order | sort -n | tr '\n' ' ')"
    actual="$(printf '%s ' $order)"
    if [ "$sorted" != "$actual" ]; then
      fail "$name: README.md sections are out of contract order (line numbers:$actual)"
    fi

    # A section heading with no body is how a botched edit announces itself.
    # Presence and order both pass on an empty section, so this is checked
    # separately: every required section needs at least one content line
    # before the next heading.
    for section in "${REQUIRED_SECTIONS[@]}"; do
      # A CRLF checkout would make every `$0 == sec` comparison fail, since
      # the line carries a trailing \r. The repo does not pin .md to LF (only
      # *.sh), so core.autocrlf=true working trees are the normal case on
      # Windows. Normalise the record before matching, and normalise $0 for
      # the heading test too -- `/^## /` matches regardless, but a heading
      # line is the one that must be compared exactly.
      body="$(awk -v sec="$section" '
        { sub(/\r$/, "") }
        $0 == sec { inside = 1; next }
        /^## /      { inside = 0 }
        inside && NF { print; found = 1; exit }
        END         { exit (found ? 0 : 1) }
      ' "$readme")" || fail "$name: README.md section '$section' has no content"
    done

    # Unbalanced code fences would render a broken example and hide the rest
    # of the block from review.
    fences="$(grep -c '^```' "$readme" || true)"
    if [ $((fences % 2)) -ne 0 ]; then
      fail "$name: README.md has an odd number of code fences ($fences)"
    fi
  fi

  # --- CHANGELOG ------------------------------------------------------------
  changelog="$dir/CHANGELOG.md"
  if [ ! -f "$changelog" ]; then
    fail "$name: missing CHANGELOG.md"
  else
    if ! grep -qF '## [Unreleased]' "$changelog"; then
      fail "$name: CHANGELOG.md has no '## [Unreleased]' section"
    fi
    if ! grep -qE '^## \[[0-9]+\.[0-9]+\.[0-9]+\] - [0-9]{4}-[0-9]{2}-[0-9]{2}' "$changelog"; then
      fail "$name: CHANGELOG.md has no dated release section (e.g. '## [0.1.0] - 2026-09-22')"
    fi
    if ! grep -qi 'Keep a Changelog' "$changelog"; then
      fail "$name: CHANGELOG.md does not reference the Keep a Changelog format"
    fi
  fi
  # --- manifest metadata ----------------------------------------------------
  declared_readme="$(pkg_field "$manifest" readme)"
  if [ "$declared_readme" != "README.md" ]; then
    fail "$name: manifest 'readme' is '${declared_readme:-<unset>}', expected 'README.md'"
  fi

  # keywords: 1-5 entries, each <= 20 chars, [A-Za-z0-9_-]
  keywords="$(pkg_array "$manifest" keywords)"
  if [ -z "$keywords" ]; then
    fail "$name: manifest declares no keywords"
  else
    kw_count=0
    while read -r kw; do
      [ -z "$kw" ] && continue
      kw_count=$((kw_count + 1))
      if [ "${#kw}" -gt 20 ]; then
        fail "$name: keyword '$kw' is ${#kw} chars; crates.io allows at most 20"
      fi
      if ! printf '%s' "$kw" | grep -qE '^[A-Za-z0-9_-]+$'; then
        fail "$name: keyword '$kw' has characters crates.io disallows ([A-Za-z0-9_-] only)"
      fi
    done <<< "$keywords"
    if [ "$kw_count" -gt 5 ]; then
      fail "$name: $kw_count keywords declared; crates.io allows at most 5"
    fi
  fi

  # categories: <= 5 entries, each a registered slug
  categories="$(pkg_array "$manifest" categories)"
  if [ -z "$categories" ]; then
    fail "$name: manifest declares no categories"
  else
    cat_count=0
    while read -r cat; do
      [ -z "$cat" ] && continue
      cat_count=$((cat_count + 1))
      if ! grep -qxF "$cat" <<< "$VALID_CATEGORIES"; then
        fail "$name: category '$cat' is not a registered crates.io category"
      fi
    done <<< "$categories"
    if [ "$cat_count" -gt 5 ]; then
      fail "$name: $cat_count categories declared; crates.io allows at most 5"
    fi
  fi
done <<EOF
$(find . -name Cargo.toml -not -path './target/*' -not -path '*/target/*' | sort)
EOF

printf '\nChecked %d workspace members.\n' "$checked"

if [ "$checked" -eq 0 ]; then
  printf 'ERROR: no workspace members were discovered\n' >&2
  exit 1
fi

if [ "$failures" -ne 0 ]; then
  printf '\n%d problem(s) found.\n' "$failures" >&2
  exit 1
fi

printf 'All per-crate documentation and metadata checks passed.\n'