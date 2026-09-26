#!/usr/bin/env bash
#
# diff-golden.sh - render a before/after numeric-drift table for the golden
# reference datasets, for pasting into a pull request.
#
# Usage:
#   scripts/diff-golden.sh [BASE_REF]
#
# BASE_REF defaults to the merge base of HEAD against origin/master, falling
# back to HEAD~1. The "after" side is the working tree.
#
# What it does:
#   * For every JSON under test-data/golden/, loads the version at BASE_REF via
#     `git show` and the working-tree version.
#   * Walks both, collecting every numeric leaf keyed by its JSON path.
#   * Prints base, head, absolute delta and percent drift.
#   * Marks each row against the `tolerance_percent` the golden file itself
#     declares for that metric, when there is one.
#
# Exit code:
#   0  no numeric drift beyond the declared tolerance
#   1  drift beyond tolerance, or a golden file was added/removed
#   2  bad usage, or python3 unavailable
#
# Requires python3 (present on GitHub-hosted runners).

set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$PWD"
GOLDEN="test-data/golden"

if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 is required (not found on PATH)." >&2
  exit 2
fi

BASE="${1:-}"
if [ -z "$BASE" ]; then
  if git rev-parse --verify origin/master >/dev/null 2>&1; then
    BASE="$(git merge-base HEAD origin/master 2>/dev/null || echo HEAD~1)"
  else
    BASE="HEAD~1"
  fi
fi

if ! git rev-parse --verify "$BASE" >/dev/null 2>&1; then
  echo "ERROR: base ref '$BASE' does not resolve in this repository." >&2
  exit 2
fi

echo "Golden dataset drift: $(git rev-parse --short "$BASE") -> working tree"
echo

mapfile -t FILES < <(find "$GOLDEN" -name '*.json' | sort)
if [ "${#FILES[@]}" -eq 0 ]; then
  echo "No golden JSON files found under $GOLDEN." >&2
  exit 2
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

# Materialise the baseline copy of every golden file.
for f in "${FILES[@]}"; do
  if git cat-file -e "$BASE:$f" 2>/dev/null; then
    mkdir -p "$tmpdir/base/$(dirname "$f")"
    git show "$BASE:$f" > "$tmpdir/base/$f"
  fi
done

python3 - "$BASE" "$tmpdir" "${FILES[@]}" <<'PYEOF'
import json, sys, os

base_ref, tmpdir = sys.argv[1], sys.argv[2]
files = sys.argv[3:]


def numbers(node, path="$", out=None):
    """Collect {json_path: number} for every numeric leaf, skipping bools."""
    if out is None:
        out = {}
    if isinstance(node, bool):
        return out
    if isinstance(node, (int, float)):
        out[path] = node
    elif isinstance(node, list):
        for i, v in enumerate(node):
            numbers(v, f"{path}[{i}]", out)
    elif isinstance(node, dict):
        for k, v in node.items():
            numbers(v, f"{path}.{k}" if path else k, out)
    return out


def tolerances(node, out=None):
    """Map a metric name to its declared tolerance_percent, where present."""
    if out is None:
        out = {}
    if isinstance(node, dict):
        if "tolerance_percent" in node and isinstance(node["tolerance_percent"], (int, float)):
            name = node.get("name") or node.get("field")
            if isinstance(name, str):
                out[name] = float(node["tolerance_percent"])
        for v in node.values():
            tolerances(v, out)
    elif isinstance(node, list):
        for v in node:
            tolerances(v, out)
    return out


def load(p):
    try:
        with open(p) as f:
            return json.load(f)
    except Exception as e:
        return {"__error__": str(e)}


rows, added, removed, parse_errors = [], [], [], []

for f in files:
    head_path = f
    base_path = os.path.join(tmpdir, "base", f)
    if not os.path.exists(base_path):
        added.append(f)
        continue
    head, base = load(head_path), load(base_path)
    if "__error__" in base or "__error__" in head:
        parse_errors.append(f)
        continue
    bn, hn = numbers(base), numbers(head)
    tol = tolerances(head) or tolerances(base)
    for k in sorted(set(bn) | set(hn)):
        if k not in bn:
            added.append(f"{f} :: {k}")
            continue
        if k not in hn:
            removed.append(f"{f} :: {k}")
            continue
        b, h = bn[k], hn[k]
        if b == h:
            continue
        delta = h - b
        pct = (abs(delta) / abs(b) * 100.0) if b else float("inf")
        rows.append((f, k, b, h, delta, pct, tol.get(k.rsplit(".", 1)[-1])))

def fmt(x):
    if isinstance(x, float) and x == int(x) and abs(x) < 1e15:
        return str(int(x))
    return f"{x:.6g}" if isinstance(x, float) else str(x)

print("| File | Field | Base | Head | Delta | Drift | Within declared tolerance |")
print("|---|---|---:|---:|---:|---:|:--:|")
for f, k, b, h, d, p, t in rows:
    tol_s = f"yes (±{t:g} %)" if t is not None and p <= t else (
        f"**no** (±{t:g} %)" if t is not None else "n/a")
    print(f"| `{os.path.basename(f)}` | `{k}` | {fmt(b)} | {fmt(h)} | {fmt(d)} | {p:.3g} % | {tol_s} |")

if not rows and not added and not removed:
    print()
    print("No numeric drift and no added or removed golden fields. "
          "This is what an unchanged numerical result looks like.")
print()
print(f"{len(rows)} drifted field(s), {len(added)} added, {len(removed)} removed, "
      f"{len(parse_errors)} unreadable.")
if parse_errors:
    print("unreadable (not valid JSON at one side): " + ", ".join(parse_errors))
for label, items in (("ADDED", added), ("REMOVED", removed)):
    for it in items:
        print(f"{label}: {it}")

violations = [r for r in rows if r[6] is not None and r[5] > r[6]]
sys.exit(1 if (violations or added or removed or parse_errors) else 0)
PYEOF