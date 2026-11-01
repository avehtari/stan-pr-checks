#!/bin/bash
#SBATCH --job-name=beta_bench
#SBATCH --output=beta_bench_%j.log
#SBATCH --time=00:40:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
# Speed of beta(double, double) in three regions and of beta_lcdf (var,
# b = 1e9), develop against the branch (PLAN §3.6).  One binary per
# function and per tree (bench_beta.cpp, -O3), built against the trees of
# beta_trees.tar.gz.  ROUNDS interleaved rounds on one core (taskset).
#
# Deployed with this job: bench_beta.cpp, beta_trees.tar.gz.
# Environment: LIB, a Stan Math lib/ directory with a built tbb/.
# Do not use `set -u`: /etc/profile reads unset variables.
set -o pipefail

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB="${LIB:?set LIB to a Stan Math lib/ directory}"
TBB="$LIB/tbb"
ROUNDS="${ROUNDS:-3}"
for f in bench_beta.cpp beta_trees.tar.gz; do
  [ -f "$f" ] || { echo "SETUP FAILED: no $f"; exit 1; }
done
echo "== host: $(hostname)   date: $(date -Is)   g++: $(g++ --version | head -1)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-)"
echo "== governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null)"
md5sum bench_beta.cpp beta_trees.tar.gz

WORK="${TMPDIR:-/tmp}/beta_bench_$SLURM_JOB_ID"
mkdir -p "$WORK" || exit 1
tar -xzf beta_trees.tar.gz -C "$WORK" || { echo "SETUP FAILED: untar"; exit 1; }
md5sum "$WORK"/*/stan/math/prim/fun/beta.hpp

for t in develop branch; do
  for fn in 1 2; do
    g++ -std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS \
      -DBENCH_FN=$fn \
      -I "$WORK/$t" -I "$LIB/eigen_5.0.1" -I "$LIB/boost_1.87.0" \
      -I "$LIB/sundials_6.1.1/include" -I "$LIB/tbb_2020.3/include" \
      bench_beta.cpp -o "$WORK/bench${fn}_$t" -L "$TBB" -ltbb \
      -Wl,-rpath,"$TBB" > "$WORK/build${fn}_$t.log" 2>&1 &
  done
done
wait
for t in develop branch; do
  for fn in 1 2; do
    [ -x "$WORK/bench${fn}_$t" ] || { echo "BUILD FAILED ${fn} $t"; cat "$WORK/build${fn}_$t.log"; exit 1; }
  done
done

CORE=$(taskset -cp $$ | sed 's/.*: //; s/[,-].*//')
echo "== pinned to core $CORE"
for r in $(seq 1 "$ROUNDS"); do
  for fn in 1 2; do
    for t in develop branch; do
      taskset -c "$CORE" "$WORK/bench${fn}_$t" 2>/dev/null \
        | sed "s/^TIME /TIME round$r $t /"
    done
  done
done
rm -rf "$WORK"
echo "EXIT=0"
