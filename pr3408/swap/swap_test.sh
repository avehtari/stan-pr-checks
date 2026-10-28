#!/bin/bash
# Bundle: regression swap for #3408.  Runs the given unit tests with the PR
# headers, then with the base-ref headers, then with the PR headers again.
# Supports the swap table of #3408 ("the new tests fail without the changes").
#
# Usage:
#   STAN_MATH=<clone of the PR branch> swap/swap_test.sh <test file> [...]
#
#   <test file>  path relative to STAN_MATH, with or without .cpp, for
#                example test/unit/math/prim/fun/grad_reg_inc_gamma_accuracy_test.cpp
#
# Environment:
#   STAN_MATH  required; a git clone or worktree.  It is read, never written.
#   BASE_REF   optional; the base (default origin/develop)
#   TBB_LIB    optional; a directory with a built libtbb.so.  Copied into the
#              scratch tree, so that make does not build TBB again.
#   SCRATCH    optional; the scratch directory (default: a new mktemp -d).
#              It must not exist yet or must be empty.
#   JOBS       optional; make -j (default 4)
#   KEEP       optional; KEEP=1 keeps the scratch tree at the end
#
# Method:
#   1. Copy the working tree of STAN_MATH (without .git) to SCRATCH.
#   2. List the stan/ files that differ between BASE_REF and the working
#      tree.  Modified (M) and deleted (D) files get the BASE_REF version
#      in the base stage.  Added files (A, untracked) stay; the base
#      headers do not include them.  test/ files are not swapped.
#   3. Stage pr_1: build and run every test.  Stage base: the same with
#      the base headers.  Stage pr_2: the same with the PR headers again.
#   4. Check that the md5 of the swapped files in pr_2 equals pr_1, and
#      that STAN_MATH is unchanged.
# A test that fails in pr_1 or pr_2 is a defect.  A new test should fail
# in the base stage.
set -o pipefail

if [ $# -lt 1 ]; then
  sed -n '2,20p' "$0" >&2
  exit 2
fi
: "${STAN_MATH:?set STAN_MATH to a clone of the PR branch}"
STAN_MATH="$(cd "$STAN_MATH" && pwd)"
BASE_REF="${BASE_REF:-origin/develop}"
JOBS="${JOBS:-4}"
G="git --no-optional-locks -C $STAN_MATH"

$G rev-parse --verify --quiet "$BASE_REF^{commit}" > /dev/null \
  || { echo "SETUP FAILED: $BASE_REF is not a commit in $STAN_MATH"; exit 1; }

tree_md5() { (cd "$1" && find stan test make -type f -print0 | sort -z \
                | xargs -0 md5sum | md5sum | cut -c1-16); }
SRC_MD5_BEFORE=$(tree_md5 "$STAN_MATH")

if [ -z "$SCRATCH" ]; then
  SCRATCH=$(mktemp -d "${TMPDIR:-/tmp}/swap_test.XXXXXX")
fi
mkdir -p "$SCRATCH" || exit 1
if [ -n "$(ls -A "$SCRATCH")" ]; then
  echo "SETUP FAILED: $SCRATCH is not empty"
  exit 1
fi
case "$SCRATCH/" in
  "$STAN_MATH"/*) echo "SETUP FAILED: SCRATCH is inside STAN_MATH"; exit 1 ;;
esac
W="$SCRATCH/tree"
KEEPDIR="$SCRATCH/keep"
mkdir -p "$W" "$KEEPDIR/pr" "$KEEPDIR/base" || exit 1

echo "== STAN_MATH $STAN_MATH, HEAD $($G rev-parse --short HEAD)," \
     "BASE_REF $BASE_REF $($G rev-parse --short "$BASE_REF")"
echo "== scratch $SCRATCH"
echo "== g++: $(g++ --version | head -1)"
[ -f "$STAN_MATH/make/local" ] && { echo "== NOTE: make/local exists:"; cat "$STAN_MATH/make/local"; }

# ── copy the working tree ───────────────────────────────────────────────
(cd "$STAN_MATH" && tar cf - --exclude=./.git .) | (cd "$W" && tar xf -) \
  || { echo "SETUP FAILED: copy"; exit 1; }
if [ -n "$TBB_LIB" ]; then
  [ -e "$TBB_LIB/libtbb.so" ] || { echo "SETUP FAILED: no libtbb.so in $TBB_LIB"; exit 1; }
  rm -rf "$W/lib/tbb" && cp -a "$TBB_LIB" "$W/lib/tbb" || exit 1
fi

# ── the files to swap ───────────────────────────────────────────────────
SWAP=""
ADDED=""
while read -r status f; do
  case "$status" in
    M|D) SWAP="$SWAP $f" ;;
    A) ADDED="$ADDED $f" ;;
    *) echo "SETUP FAILED: status $status of $f"; exit 1 ;;
  esac
done < <($G diff --no-renames --name-status "$BASE_REF" -- stan)
ADDED="$ADDED $($G ls-files --others --exclude-standard -- stan | tr '\n' ' ')"
echo "== stan/ files swapped (M, D): $(echo $SWAP | wc -w)"
for f in $SWAP; do echo "   $f"; done
echo "== stan/ files added by the PR, not swapped: $(echo $ADDED | wc -w)"
for f in $SWAP; do
  mkdir -p "$KEEPDIR/pr/$(dirname "$f")" "$KEEPDIR/base/$(dirname "$f")"
  if [ -e "$W/$f" ]; then
    cp -p "$W/$f" "$KEEPDIR/pr/$f" || exit 1
  fi
  $G show "$BASE_REF:$f" > "$KEEPDIR/base/$f" || exit 1
done

TESTS=""
for t in "$@"; do
  t="${t%.cpp}"
  [ -f "$W/$t.cpp" ] || { echo "SETUP FAILED: missing $t.cpp"; exit 1; }
  TESTS="$TESTS $t"
done

swap_md5() { (cd "$W" && for f in $SWAP; do
                [ -e "$f" ] && md5sum "$f" || echo "absent $f"; done | md5sum | cut -c1-16); }

install() {
  side="$1"
  for f in $SWAP; do
    if [ -e "$KEEPDIR/$side/$f" ]; then
      cp -p "$KEEPDIR/$side/$f" "$W/$f" || exit 1
    else
      rm -f "$W/$f"
    fi
  done
}

run_stage() {
  tag="$1"
  echo ""
  echo "######## stage $tag ########"
  echo "md5 of the swapped files: $(swap_md5)"
  failed=0
  for t in $TESTS; do
    base=$(basename "$t")
    rm -f "$W/$t" "$W/$t.o" "$W/$t.d"
    if ! (cd "$W" && make -j"$JOBS" "$t") > "$SCRATCH/build_${tag}_$base.log" 2>&1; then
      echo "  BUILD FAILED  $base"
      grep -m 10 "error:" "$SCRATCH/build_${tag}_$base.log"
      failed=$((failed + 1))
      continue
    fi
    log="$SCRATCH/run_${tag}_$base.log"
    (cd "$W" && "./$t") > "$log" 2>&1
    status=$?
    nok=$(grep -c '^\[       OK \]' "$log")
    nfail=$(grep -c '^\[  FAILED  \].*(' "$log")
    nexp=$(grep -c 'Failure$' "$log")
    if [ $status -eq 0 ]; then
      echo "  PASS  $base  ($nok cases)"
    else
      failed=$((failed + 1))
      echo "  FAIL  $base  (exit $status; cases ok $nok, failed $nfail;" \
           "failed expectations $nexp)"
      grep '^\[  FAILED  \].*(' "$log" | head -20
    fi
  done
  echo "== $tag: suites not passing: $failed"
}

run_stage pr_1
MD5_PR=$(swap_md5)
install base
run_stage base
install pr
run_stage pr_2
MD5_PR2=$(swap_md5)

echo ""
if [ "$MD5_PR" = "$MD5_PR2" ]; then
  echo "== restore check: swapped files equal stage pr_1 ($MD5_PR)"
else
  echo "== RESTORE MISMATCH: pr_1 $MD5_PR, pr_2 $MD5_PR2"
fi
SRC_MD5_AFTER=$(tree_md5 "$STAN_MATH")
if [ "$SRC_MD5_BEFORE" = "$SRC_MD5_AFTER" ]; then
  echo "== STAN_MATH unchanged ($SRC_MD5_BEFORE)"
else
  echo "== WARNING: STAN_MATH changed during the run"
fi
if [ "${KEEP:-0}" = "1" ]; then
  echo "== scratch kept: $SCRATCH"
else
  rm -rf "$SCRATCH"
fi
