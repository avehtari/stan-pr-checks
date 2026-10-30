#!/bin/bash
#SBATCH --job-name=c3406_bench
#SBATCH --time=00:25:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=c3406_bench_%j.log

# Speed of std_normal_lcdf_value_grad: the #3406 head (pr) against local
# commit 5f33a10d31 (commit).  Stan default flags.  Binaries alternate,
# 3 runs x 3 rounds.  Do not use `set -u`: /etc/profile reads unset
# variables.

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB=/scratch/work/ave/smtest_lib_d12f488749
WORK=${TMPDIR:-/tmp}/c3406_bench_$SLURM_JOB_ID

echo "=== node ==="
hostname
grep -m1 "model name" /proc/cpuinfo
echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1

mkdir -p "$WORK/commit" || exit 1
tar -xzf trees.tar.gz -C "$WORK" pr3406/stan || { echo "FATAL: untar pr"; exit 1; }
tar -xzf commit_5f33a10d31.tar.gz -C "$WORK/commit" stan || { echo "FATAL: untar commit"; exit 1; }
md5sum "$WORK"/pr3406/stan/math/prim/fun/std_normal_lcdf_impl.hpp \
  "$WORK"/commit/stan/math/prim/fun/std_normal_lcdf_impl.hpp \
  "$WORK"/commit/stan/math/prim/fun/erfcx.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
g++ $FLAGS -I "$WORK/pr3406" $INCS -o bench_pr bench_dedup.cpp || { echo "FATAL build pr"; exit 1; }
g++ $FLAGS -I "$WORK/commit" $INCS -o bench_commit bench_dedup.cpp \
  || { echo "FATAL build commit"; exit 1; }

for run in 1 2 3; do
  ./bench_pr pr 3
  ./bench_commit commit 3
done

rm -rf "$WORK"
echo "EXIT=0"
