#!/bin/bash
#SBATCH --job-name=cdfs3_bench
#SBATCH --time=00:40:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=cdfs3_bench_%j.log

# Speed of std_normal_cdf, normal_cdf, lognormal_cdf, exp_mod_normal_cdf on
# three trees: develop e4745c87ff, #3406 head + bf5861ab9c (exp(lcdf)), and
# the direct erfc form (worktree, not committed).  Stan default flags.
# Binaries alternate, 3 runs x 3 rounds.
# Do not use `set -u`: /etc/profile reads unset variables.

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB=/scratch/work/ave/smtest_lib_d12f488749
WORK=${TMPDIR:-/tmp}/cdfs3_bench_$SLURM_JOB_ID

hostname
grep -m1 "model name" /proc/cpuinfo
echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1

mkdir -p "$WORK/dev" "$WORK/explcdf" "$WORK/direct" || exit 1
tar -xzf develop_e4745c87ff_stan.tar.gz -C "$WORK/dev" || { echo "FATAL: untar dev"; exit 1; }
tar -xzf commit_bf5861ab9c.tar.gz -C "$WORK/explcdf" stan || { echo "FATAL: untar explcdf"; exit 1; }
tar -xzf direct_cdfs_stan.tar.gz -C "$WORK/direct" || { echo "FATAL: untar direct"; exit 1; }
md5sum "$WORK"/*/stan/math/prim/prob/normal_cdf.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
LINK="-L$LIB/tbb -Wl,-rpath,$LIB/tbb -ltbb"
for t in dev explcdf direct; do
  g++ $FLAGS -I "$WORK/$t" $INCS -o "bench_$t" bench_cdfs.cpp $LINK \
    || { echo "FATAL build $t"; exit 1; }
done

for run in 1 2 3; do
  ./bench_dev develop 3
  ./bench_explcdf explcdf 3
  ./bench_direct direct 3
done

rm -rf "$WORK"
echo "EXIT=0"
