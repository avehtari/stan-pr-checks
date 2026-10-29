#!/bin/bash
#SBATCH --job-name=pr3442_speed
#SBATCH --output=pr3442_speed_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#
# SLURM job for the speed benchmarks of stan-dev/math#3442: compile each
# benchmark against each Stan Math tree with -O3, then run the binaries one
# after the other on one exclusive node, REPS interleaved rounds
# (round 1: tree A, tree B; round 2: tree A, tree B; ...).
#
# Environment (pass with sbatch --export):
#   PROBES         required; .cpp files in the submit directory, separated by ':'
#   TREES          required; tarballs of a stan/ directory, separated by ':'
#   STAN_MATH_LIB  required; lib/ of a Stan Math clone with a built lib/tbb
#   MODULES        optional; modules to load, separated by ':'
#   REPS           optional; number of rounds (default 3)
# Output: <probe>_bench.out, each line prefixed with <tree> <round>; table
# with speed/table_bench.py.  Do not use `set -u` or `set -e`.
set -o pipefail

[ -f /etc/profile ] && source /etc/profile
[ -n "$MODULES" ] && module load ${MODULES//:/ }

: "${PROBES:?set PROBES}"
: "${TREES:?set TREES}"
LIB="${STAN_MATH_LIB:?set STAN_MATH_LIB}"
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
PROBES="${PROBES//:/ }"
TREES="${TREES//:/ }"
REPS="${REPS:-3}"
FLAGS="-std=c++17 -O3 -DNDEBUG -D_REENTRANT -pthread"
echo "== host: $(hostname)   date: $(date -Is)   reps: $REPS"
echo "== probes: $PROBES   trees: $TREES"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-)"
echo "== governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
echo "== compiler: $(g++ --version | head -1)"
echo "== flags: $FLAGS"
for f in $PROBES $TREES; do
  [ -f "$f" ] || { echo "SETUP FAILED: missing $f"; exit 1; }
done
md5sum $PROBES $TREES

WORK=$(mktemp -d "${TMPDIR:-/tmp}/pr3442_bench_XXXXXX") || exit 1
trap 'rm -rf "$WORK"' EXIT

tree_name() { basename "$1" .tar.gz; }

for t in $TREES; do
  name="$(tree_name "$t")"
  mkdir -p "$WORK/$name" && tar xzf "$t" -C "$WORK/$name" || exit 1
done
for p in $PROBES; do
  for t in $TREES; do
    name="$(tree_name "$t")"
    # shellcheck disable=SC2086
    g++ $FLAGS -I "$WORK/$name" -I "$LIB/eigen_5.0.1" -I "$LIB/boost_1.87.0" \
      -I "$LIB/sundials_6.1.1/include" -I "$LIB/tbb_2020.3/include" \
      "$p" -o "$WORK/${p%.cpp}_$name" \
      -L "$LIB/tbb" -ltbb -Wl,-rpath,"$LIB/tbb" \
      > "build_${p%.cpp}_${name}.log" 2>&1 &
  done
done
wait

for p in $PROBES; do
  : > "${p%.cpp}_bench.out"
done
for r in $(seq 1 "$REPS"); do
  for p in $PROBES; do
    for t in $TREES; do
      name="$(tree_name "$t")"
      bin="$WORK/${p%.cpp}_$name"
      if [ ! -x "$bin" ]; then
        echo "BUILD FAILED  $p $name"
        tail -20 "build_${p%.cpp}_${name}.log"
        continue
      fi
      "$bin" | sed "s/^/$name\t$r\t/" >> "${p%.cpp}_bench.out"
      echo "ran  $p $name round $r: exit $?  $(date +%T)"
    done
  done
done
echo "== done: $(date -Is)"
