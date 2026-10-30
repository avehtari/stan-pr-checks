#!/bin/bash
#SBATCH --job-name=phiv_bench
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=phiv_bench_%j.log

# Speed of the Phi(double) variants (phi_variants.hpp) against the #3406
# head plus local commit 5f33a10d31.  Stan default flags.  One binary, the
# variants interleaved per round, 3 runs.
# Do not use `set -u`: /etc/profile reads unset variables.

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB=/scratch/work/ave/smtest_lib_d12f488749
WORK=${TMPDIR:-/tmp}/phiv_bench_$SLURM_JOB_ID

hostname
grep -m1 "model name" /proc/cpuinfo
echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1

mkdir -p "$WORK" && tar -xzf commit_5f33a10d31.tar.gz -C "$WORK" stan \
  || { echo "FATAL: untar"; exit 1; }
md5sum "$WORK/stan/math/prim/fun/std_normal_lcdf_impl.hpp" phi_variants.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-I $WORK -I . -isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
g++ $FLAGS $INCS -o bench_phiv bench_phi_variants.cpp || { echo "FATAL build"; exit 1; }

for run in 1 2 3; do
  ./bench_phiv 3
done

rm -rf "$WORK"
echo "EXIT=0"
