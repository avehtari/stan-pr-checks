#!/bin/bash
# Bundle: packs pr3408_fix.tar.gz for run_roihu_bench.sh and
# run_bench_lccdf_tail.sh (the speed tables of #3408 Part 3).
#
# Build "TREE_REF + the PR diff + the working-tree changes" and pack it
# with the header variant for the speed test:
#   tree/           the merged tree
#   variants/orig/  the version at ORIG_REF of every stan/ file that the
#                   working tree modifies relative to ORIG_REF
#   lib_overlay/    lib files of TREE_REF that CmdStan 2.40.0 lacks
# For the speed tables of the PR text: ORIG_REF=origin/develop.
#
#   WT        required; a git clone or worktree of the PR branch (read only)
#   TREE_REF  optional; default origin/develop
#   ORIG_REF  optional; default ad3ba6982e (the PR head of 2026-10-21 noon)
# Output: pr3408_fix.tar.gz in this directory.
set -euo pipefail

: "${WT:?set WT to a clone of the PR branch}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${TMPDIR:-/tmp}/sm_pr3408_merged"
G="git --no-optional-locks -C $WT"

TREE_REF="${TREE_REF:-origin/develop}"
BASE=$($G merge-base HEAD "$TREE_REF")
echo "== tree $TREE_REF $($G rev-parse --short "$TREE_REF")," \
  "merge base $($G rev-parse --short "$BASE")," \
  "PR head $($G rev-parse --short HEAD)"
rm -rf "$OUT" "$OUT.patch"
mkdir -p "$OUT/tree"
$G archive "$TREE_REF" -- . ':!lib' | tar x -C "$OUT/tree"
$G diff --binary "$BASE" -- . ':!lib' > "$OUT.patch"
(cd "$OUT/tree" && git apply --whitespace=nowarn "$OUT.patch")
$G ls-files --others --exclude-standard -- stan test | while read -r f; do
  mkdir -p "$OUT/tree/$(dirname "$f")"
  cp "$WT/$f" "$OUT/tree/$f"
  echo "   untracked copied: $f"
done

CHANGED=$( ($G diff --name-only "$BASE" -- . ':!lib';
            $G ls-files --others --exclude-standard -- stan test) | sort -u)
BOTH=$($G diff --name-only "$BASE" "$TREE_REF" -- $CHANGED)
for f in $BOTH; do
  echo "== changed on both sides: $f; merged tree minus working tree:"
  diff "$WT/$f" "$OUT/tree/$f" | sed 's/^/   /' || true
done
for f in $CHANGED; do
  case " $(echo $BOTH) " in *" $f "*) continue ;; esac
  cmp -s "$WT/$f" "$OUT/tree/$f" || { echo "MERGE MISMATCH: $f"; exit 1; }
done
echo "== the other changed files equal the working tree"

ORIG_REF="${ORIG_REF:-ad3ba6982e}"
# only modified files: a file that the working tree adds has no version at
# ORIG_REF, and the baseline headers do not include it
ORIG_FILES=$($G diff --name-only --diff-filter=M "$ORIG_REF" -- stan)
for f in $ORIG_FILES; do
  mkdir -p "$OUT/variants/orig/$(dirname "$f")"
  $G show "$ORIG_REF:$f" > "$OUT/variants/orig/$f"
done
printf '%s\n' $ORIG_FILES > "$OUT/variants/orig_files.txt"
echo "== swap baseline $ORIG_REF: $(echo $ORIG_FILES | wc -w) files"

# lib files of TREE_REF that the CmdStan 2.40.0 lib (v5.4.0) lacks
LIBFILES=$($G diff --name-only --diff-filter=AM v5.4.0 "$TREE_REF" -- lib)
mkdir -p "$OUT/lib_overlay/lib"
# an empty path list would archive the whole repository
[ -n "$LIBFILES" ] \
  && $G archive "$TREE_REF" -- $LIBFILES | tar x -C "$OUT/lib_overlay"
echo "== lib overlay: $(find "$OUT/lib_overlay" -type f | wc -l) files"

cp "$HERE/bench_lccdf_tail.cpp" "$HERE/bench_cdfs.cpp" "$OUT/"
tar -czf "$HERE/pr3408_fix.tar.gz" -C "$OUT" tree variants lib_overlay \
  bench_lccdf_tail.cpp bench_cdfs.cpp
ls -la "$HERE/pr3408_fix.tar.gz"
