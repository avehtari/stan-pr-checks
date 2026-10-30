#!/bin/bash
#SBATCH --job-name=phi_bench
#SBATCH --time=00:25:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=phi_bench_%j.log

# Speed of Phi (double, and var with gradient): develop e4745c87ff against
# the #3406 head plus local commit 5f33a10d31.  Stan default flags.
# Binaries alternate, 3 runs x 3 rounds.
# Do not use `set -u`: /etc/profile reads unset variables.

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB=/scratch/work/ave/smtest_lib_d12f488749
WORK=${TMPDIR:-/tmp}/phi_bench_$SLURM_JOB_ID

hostname
grep -m1 "model name" /proc/cpuinfo
echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1

mkdir -p "$WORK/dev" "$WORK/commit" || exit 1
tar -xzf develop_e4745c87ff_stan.tar.gz -C "$WORK/dev" || { echo "FATAL: untar dev"; exit 1; }
tar -xzf commit_5f33a10d31.tar.gz -C "$WORK/commit" stan || { echo "FATAL: untar commit"; exit 1; }
md5sum "$WORK/dev/stan/math/prim/fun/Phi.hpp" "$WORK/commit/stan/math/prim/fun/Phi.hpp"

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
LINK="-L$LIB/tbb -Wl,-rpath,$LIB/tbb -ltbb"
g++ $FLAGS -I "$WORK/dev" $INCS -o bench_dev bench_phi.cpp $LINK || { echo "FATAL build dev"; exit 1; }
g++ $FLAGS -I "$WORK/commit" $INCS -o bench_commit bench_phi.cpp $LINK \
  || { echo "FATAL build commit"; exit 1; }

for run in 1 2 3; do
  ./bench_dev develop 3
  ./bench_commit commit 3
done

rm -rf "$WORK"
echo "EXIT=0"
