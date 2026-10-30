#!/bin/bash
#SBATCH --job-name=bisect_bench
#SBATCH --time=00:45:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=bisect_bench_%j.log

# Which #3406 commit made the var lcdfs slower?  bench_pr_table.cpp in quick
# mode (std_normal_lcdf, normal_lcdf, lognormal_lcdf, exp_mod_normal_lcdf;
# z on [-3, 3]; N = 1000; double and Matrix<var>) on every non-formatting
# commit from a28c1c1637 (2026-09-18) to f4cd10ee4a, plus develop
# e4745c87ff.  Builds in parallel, then runs one binary at a time,
# 2 passes x 3 rounds.
# Do not use `set -u`: /etc/profile reads unset variables.

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB=/scratch/work/ave/smtest_lib_d12f488749
WORK=${TMPDIR:-/tmp}/bisect_bench_$SLURM_JOB_ID
TREES=${TREES:-"a28c1c1637 7fa96c9adf 9a8deffbef 644ef934bd 8cff6fd3f3 7e99d35073
201dcab952 f37ea7cb67 686fc67261 5f33a10d31 f4cd10ee4a"}
TREES=${TREES//:/ }
TARBALL=${TARBALL:-bisect_trees.tar.gz}
echo "TREES=$TREES TARBALL=$TARBALL"

hostname
grep -m1 "model name" /proc/cpuinfo
echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1

mkdir -p "$WORK/t_develop" || exit 1
tar -xzf "$TARBALL" -C "$WORK" || { echo "FATAL: untar"; exit 1; }
tar -xzf develop_e4745c87ff_stan.tar.gz -C "$WORK/t_develop" || { echo "FATAL: untar dev"; exit 1; }
md5sum "$WORK"/t_*/stan/math/prim/prob/exp_mod_normal_lcdf.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
LINK="-L$LIB/tbb -Wl,-rpath,$LIB/tbb -ltbb"
ALL="develop $TREES"
for t in $ALL; do
  g++ $FLAGS -I "$WORK/t_$t" $INCS -o "bench_$t" bench_pr_table.cpp $LINK \
    > "build_$t.log" 2>&1 &
  while [ "$(jobs -r | wc -l)" -ge 8 ]; do sleep 2; done
done
wait
for t in $ALL; do
  [ -x "bench_$t" ] || { echo "FATAL build $t"; cat "build_$t.log"; exit 1; }
done
echo "builds done $(date -Is)"

for pass in 1 2; do
  for t in $ALL; do
    "./bench_$t" "$t" 3 quick
  done
done

rm -rf "$WORK"
echo "EXIT=0"
