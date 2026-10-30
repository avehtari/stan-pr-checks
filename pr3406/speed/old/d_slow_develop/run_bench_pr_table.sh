#!/bin/bash
#SBATCH --job-name=prtab_bench
#SBATCH --time=00:45:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=prtab_bench_%j.log

# Reproduction of the #3406 PR-text timing table (bench_pr_table.cpp) on
# four trees: develop 5252d51d47 (base of the PR at creation), the PR at
# creation a28c1c1637, develop e4745c87ff, and the head f4cd10ee4a.
# Stan default flags.  Binaries alternate, 2 runs x 2 rounds.
# Do not use `set -u`: /etc/profile reads unset variables.

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB=/scratch/work/ave/smtest_lib_d12f488749
WORK=${TMPDIR:-/tmp}/prtab_bench_$SLURM_JOB_ID

hostname
grep -m1 "model name" /proc/cpuinfo
echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1

mkdir -p "$WORK/dev0916" "$WORK/pr0918" "$WORK/devnow" "$WORK/head" || exit 1
tar -xzf develop_5252d51d47_stan.tar.gz -C "$WORK/dev0916" || { echo "FATAL: untar"; exit 1; }
tar -xzf pr3406_a28c1c1637_stan.tar.gz -C "$WORK/pr0918" || { echo "FATAL: untar"; exit 1; }
tar -xzf develop_e4745c87ff_stan.tar.gz -C "$WORK/devnow" || { echo "FATAL: untar"; exit 1; }
tar -xzf commit_f4cd10ee4a.tar.gz -C "$WORK/head" stan || { echo "FATAL: untar"; exit 1; }
md5sum "$WORK"/*/stan/math/prim/fun/Phi.hpp "$WORK"/*/stan/math/prim/prob/normal_lcdf.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
LINK="-L$LIB/tbb -Wl,-rpath,$LIB/tbb -ltbb"
for t in dev0916 pr0918 devnow head; do
  g++ $FLAGS -I "$WORK/$t" $INCS -o "bench_$t" bench_pr_table.cpp $LINK \
    || { echo "FATAL build $t"; exit 1; }
done

for run in 1 2; do
  for t in dev0916 pr0918 devnow head; do
    "./bench_$t" "$t" 2
  done
done

rm -rf "$WORK"
echo "EXIT=0"
