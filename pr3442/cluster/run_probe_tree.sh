#!/bin/bash
#SBATCH --job-name=pr3442_acc
#SBATCH --output=pr3442_acc_%j.log
#SBATCH --time=02:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#
# SLURM job for the accuracy probes of stan-dev/math#3442: compile one probe
# against each Stan Math tree and run it once per tree.
#
# Environment (pass with sbatch --export):
#   PROBE          required; a .cpp in the submit directory
#   TREES          required; tarballs of a stan/ directory, separated by ':'
#                  (cluster/make_trees.sh).  Tree name = file name without
#                  .tar.gz; the output of tree T is <probe>_<T>.out.
#   STAN_MATH_LIB  required; lib/ of a Stan Math clone with a built lib/tbb
#   MODULES        optional; modules to load, separated by ':'
#   PRECS, PTS     optional; read by some probes (see their headers)
# The trees are built in parallel in node-local scratch; then each binary
# runs alone, in the order of TREES.  Outputs go to the submit directory.
# Do not use `set -u` or `set -e`.
set -o pipefail

[ -f /etc/profile ] && source /etc/profile
[ -n "$MODULES" ] && module load ${MODULES//:/ }

: "${PROBE:?set PROBE}"
: "${TREES:?set TREES}"
LIB="${STAN_MATH_LIB:?set STAN_MATH_LIB}"
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
TREES="${TREES//:/ }"
echo "== host: $(hostname)   date: $(date -Is)   probe: $PROBE   trees: $TREES"
echo "== env: PRECS=${PRECS:-unset} PTS=${PTS:-unset}"
echo "== compiler: $(g++ --version | head -1)"
echo "== flags: -std=c++17 -O2 -DNDEBUG -D_REENTRANT -pthread"
[ -f "$PROBE" ] || { echo "SETUP FAILED: missing $PROBE"; exit 1; }
for t in $TREES; do
  [ -f "$t" ] || { echo "SETUP FAILED: missing $t"; exit 1; }
done
md5sum "$PROBE" ./*.hpp $TREES 2>/dev/null

WORK=$(mktemp -d "${TMPDIR:-/tmp}/pr3442_probe_XXXXXX") || exit 1
trap 'rm -rf "$WORK"' EXIT

tree_name() { basename "$1" .tar.gz; }

build() {
  name="$(tree_name "$1")"
  mkdir -p "$WORK/$name" && tar xzf "$1" -C "$WORK/$name" || return 1
  g++ -std=c++17 -O2 -DNDEBUG -D_REENTRANT -pthread -I . \
    -I "$WORK/$name" -I "$LIB/eigen_5.0.1" -I "$LIB/boost_1.87.0" \
    -I "$LIB/sundials_6.1.1/include" -I "$LIB/tbb_2020.3/include" \
    "$PROBE" -o "$WORK/probe_$name" \
    -L "$LIB/tbb" -ltbb -Wl,-rpath,"$LIB/tbb" \
    > "build_${PROBE%.cpp}_${name}.log" 2>&1
}
for t in $TREES; do build "$t" & done
wait

for t in $TREES; do
  name="$(tree_name "$t")"
  out="${PROBE%.cpp}_${name}.out"
  if [ ! -x "$WORK/probe_$name" ]; then
    echo "BUILD FAILED  $name"
    tail -20 "build_${PROBE%.cpp}_${name}.log"
    continue
  fi
  /usr/bin/time -f "   $name wall %e s" "$WORK/probe_$name" > "$out" 2> "${out%.out}.err"
  echo "ran  $name: exit $?  rows $(grep -c . "$out")  throws $(grep -c THROW "$out")"
  tail -1 "${out%.out}.err"
done
echo "== done: $(date -Is)"
