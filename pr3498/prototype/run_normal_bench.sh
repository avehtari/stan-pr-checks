#!/bin/bash
#SBATCH --job-name=nv_bench
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=nv_bench_%j.log

# normal_lpdf cost parts (prototype for stan-dev/math#3498).  Builds
# check_normal.cpp and one bench binary per variant against
# the origin/develop tree in the tarball, runs the check, then 5
# interleaved timing rounds on one core.
# Output: CHECK, RESULT r<round> <variant> <N> <median_ns> <min_ns> <reps>.
# Do not use `set -u`: /etc/profile reads unset variables.

source /etc/profile
module --ignore_cache load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR" || exit 1
echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)   cpus: ${SLURM_CPUS_PER_TASK:-?}"
echo "== governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"

TARBALL=dev_dc66941969.tar.gz
SRCS="normal_variants.hpp bench_normal.cpp check_normal.cpp $TARBALL"
for f in $SRCS; do
  [ -f "$f" ] || { echo "SETUP FAILED: missing $f"; echo "EXIT=1"; exit 1; }
done
md5sum $SRCS

# lib/ of a Stan Math clone at develop dc66941969 with a built TBB
LIB=${STAN_MATH_LIB:?set STAN_MATH_LIB to the lib/ directory of a built Stan Math clone}
[ -f "$LIB/tbb/libtbb.so" ] || { echo "SETUP FAILED: no $LIB/tbb/libtbb.so"; echo "EXIT=1"; exit 1; }
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; echo "EXIT=1"; exit 1; }
W="$TMPDIR/nv"
rm -rf "$W" && mkdir -p "$W" || exit 1
tar -xzf "$TARBALL" -C "$W" || { echo "SETUP FAILED: tar"; echo "EXIT=1"; exit 1; }
cp normal_variants.hpp bench_normal.cpp check_normal.cpp "$W/" || exit 1

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-I $W/tree -isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/sundials_6.1.1/include -isystem $LIB/tbb_2020.3/include"
LINK="-L$LIB/tbb -Wl,-rpath,$LIB/tbb -ltbb"

VARIANTS="0 1 2 3 4 5 6 7"
echo "== build"
t0=$(date +%s)
( g++ $FLAGS $INCS -o "$W/check" "$W/check_normal.cpp" $LINK > "$W/build_check.log" 2>&1 \
    || echo "BUILD FAILED check" ) &
for v in $VARIANTS; do
  ( g++ $FLAGS -DNV=$v $INCS -o "$W/b_$v" "$W/bench_normal.cpp" $LINK \
      > "$W/build_$v.log" 2>&1 || echo "BUILD FAILED $v" ) &
  # at most SLURM_CPUS_PER_TASK compiles at a time
  while [ "$(jobs -r | wc -l)" -ge "${SLURM_CPUS_PER_TASK:-1}" ]; do sleep 2; done
done
wait
echo "== build done in $(( $(date +%s) - t0 )) s"
for f in "$W"/build_*.log; do
  if [ -s "$f" ]; then echo "--- $f"; head -20 "$f"; fi
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
  for v in $VARIANTS; do
    if [ -x "$W/b_$v" ]; then
      taskset -c "$CORE" "$W/b_$v" | sed "s/^RESULT /RESULT r$round /"
    else
      echo "RESULT r$round MISSING $v"
    fi
  done
done
echo "== end $(date -Is)"
echo "EXIT=0"
