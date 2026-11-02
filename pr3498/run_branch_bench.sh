#!/bin/bash
#SBATCH --job-name=nl_branch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:50:00
#SBATCH --output=nl_branch_%j.log

# normal_lpdf develop / branch comparison (stan-dev/math#3498).  The
# tarball holds the branch stan/ tree (make_inputs.sh); the develop
# function is the renamed copy dev_normal_lpdf.hpp (normal_lpdf_dev), so one
# binary per (tree, case, propto) instantiates only one of the two.
# Builds check_branch.cpp and the bench binaries, runs the check, then 5
# interleaved timing rounds (dev, branch per case) on one core.
# Output: CHECK, HEAP, RESULT r<round> <tree> <case> <propto> <N> <median_ns>
#   <min_ns> <reps>.
# Do not use `set -u`: /etc/profile reads unset variables.

source /etc/profile
module --ignore_cache load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR" || exit 1
echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)   cpus: ${SLURM_CPUS_PER_TASK:-?}"
echo "== governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"

TARBALL=${TARBALL:-br_3a66fcf6aa.tar.gz}
SRCS="cases.hpp dev_normal_lpdf.hpp check_branch.cpp bench_branch.cpp syntax_cases.cpp $TARBALL"
for f in $SRCS; do
  [ -f "$f" ] || { echo "SETUP FAILED: missing $f"; echo "EXIT=1"; exit 1; }
done
md5sum $SRCS

# lib/ of a Stan Math clone at develop dc66941969 with a built TBB
LIB=${STAN_MATH_LIB:?set STAN_MATH_LIB to the lib/ directory of a built Stan Math clone}
[ -f "$LIB/tbb/libtbb.so" ] || { echo "SETUP FAILED: no $LIB/tbb/libtbb.so"; echo "EXIT=1"; exit 1; }
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; echo "EXIT=1"; exit 1; }
W="$TMPDIR/nl"
rm -rf "$W" && mkdir -p "$W" || exit 1
tar -xzf "$TARBALL" -C "$W" || { echo "SETUP FAILED: tar"; echo "EXIT=1"; exit 1; }
cp cases.hpp dev_normal_lpdf.hpp check_branch.cpp bench_branch.cpp syntax_cases.cpp "$W/" || exit 1

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-I $W -I $W/tree -isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/sundials_6.1.1/include -isystem $LIB/tbb_2020.3/include"
LINK="-L$LIB/tbb -Wl,-rpath,$LIB/tbb -ltbb"

# case index (cases.hpp) and propto: A-G propto 1, H propto 0, A and C propto 0
CONFIGS="0:1 1:1 2:1 3:1 4:1 5:1 6:1 7:0 0:0 2:0"
# BENCH_CONFIGS=c:p+c:p+... replaces the list (sbatch --export splits at ',')
[ -n "$BENCH_CONFIGS" ] && CONFIGS=$(echo "$BENCH_CONFIGS" | tr '+' ' ')
# CHECK_ONLY=1: build and run the check (and the syntax cases) only
[ -n "$CHECK_ONLY" ] && CONFIGS=""
MAXJ=$(( ${SLURM_CPUS_PER_TASK:-1} < 6 ? ${SLURM_CPUS_PER_TASK:-1} : 6 ))
wait_slot() {
  while [ "$(jobs -r | wc -l)" -ge "$MAXJ" ]; do sleep 2; done
}
echo "== build (at most $MAXJ compiles at a time)"
t0=$(date +%s)
( g++ $FLAGS $INCS -o "$W/check" "$W/check_branch.cpp" $LINK > "$W/build_check.log" 2>&1 \
    || echo "BUILD FAILED check" ) &
wait_slot
( g++ $FLAGS -fsyntax-only $INCS "$W/syntax_cases.cpp" > "$W/build_syntax.log" 2>&1 \
    && echo "SYNTAX_CASES OK" || echo "SYNTAX_CASES FAILED" ) &
for cfg in $CONFIGS; do
  c=${cfg%%:*}; p=${cfg##*:}
  for t in 0 1; do
    wait_slot
    ( g++ $FLAGS -DCASE=$c -DTREE=$t -DPROPTO=$p $INCS -o "$W/b_${t}_${c}_${p}" \
        "$W/bench_branch.cpp" $LINK > "$W/build_${t}_${c}_${p}.log" 2>&1 \
        || echo "BUILD FAILED $t $c $p" ) &
  done
done
wait
echo "== build done in $(( $(date +%s) - t0 )) s"
for f in "$W"/build_*.log; do
  if grep -q "error" "$f"; then echo "--- $f"; grep -m20 "error" "$f"; fi
done

echo "== check"
if [ -x "$W/check" ]; then
  "$W/check"
  echo "CHECK_EXIT=$?"
else
  echo "CHECK_EXIT=BUILD_FAILED"
fi

CORE=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}' | cut -d, -f1 | cut -d- -f1)
echo "== timing on core $CORE"
for round in 1 2 3 4 5; do
  echo "== round $round $(date -Is)"
  for cfg in $CONFIGS; do
    c=${cfg%%:*}; p=${cfg##*:}
    for t in 0 1; do
      b="$W/b_${t}_${c}_${p}"
      if [ -x "$b" ]; then
        taskset -c "$CORE" "$b" | sed "s/^RESULT /RESULT r$round /"
      else
        echo "RESULT r$round MISSING $t $c $p"
      fi
    done
  done
done
echo "== end $(date -Is)"
echo "EXIT=0"
