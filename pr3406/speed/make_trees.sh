#!/bin/bash
# Make the tree tarball for run_bench_roihu.sh: t_<commit>/stan for each
# commit and lib/sundials_6.1.1/include of the first commit (the CmdStan
# 2.40.0 lib/ has no arkode headers).  Read-only use of the clone.
# Usage: CLONE=<stan math clone> ./make_trees.sh <out.tar.gz> <commit>...
set -eo pipefail
: "${CLONE:?set CLONE to a Stan Math clone}"
[ $# -ge 2 ] || { echo "usage: CLONE=<clone> $0 <out.tar.gz> <commit>..." >&2; exit 2; }
OUT=$(realpath -m "$1")
shift
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
for c in "$@"; do
  mkdir -p "$W/t_$c"
  git -C "$CLONE" archive "$c" stan | tar -x -C "$W/t_$c"
done
git -C "$CLONE" archive "$1" lib/sundials_6.1.1/include | tar -x -C "$W"
tar -czf "$OUT" -C "$W" $(cd "$W" && ls -d t_* lib)
for c in "$@"; do
  echo "t_$c md5 of stan/: $(cd "$W/t_$c" && find stan -type f | sort | xargs md5sum | md5sum | cut -d' ' -f1)"
done
echo "wrote $OUT"
