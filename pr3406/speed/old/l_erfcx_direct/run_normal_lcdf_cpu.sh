#!/bin/bash
#SBATCH --job-name=nlcdfcpu
#SBATCH --time=00:25:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=normal_lcdf_cpu_%j.log

# Kernel-level speed of the normal_lcdf value and gradient: origin/develop
# against the erfcx identities.
#
# Both variants are in one binary, so they run on the same node in the same
# job. Three alternating repetitions expose drift. SLURM's cgroup already
# confines the job to its allocated cores, so no taskset is needed; taskset
# -c 0 fails because core 0 is usually not in the allowed set.
#
# Do not use `set -u`: /etc/profile reads unset variables.
#
# Prerequisite: the stable-normal tree at $SMDIR, deployed with
#   rsync -az --delete stan/ triton.aalto.fi:~/smnormal/stan/
#   rsync -az --delete lib/eigen_5.0.1/ triton.aalto.fi:~/smnormal/lib/eigen_5.0.1/
# Boost and TBB are borrowed from the cmdstan copy; their versions match.
# Eigen does not: cmdstan has 3.4.0 and Stan Math develop needs 5.0.1.

cd "$SLURM_SUBMIT_DIR" || exit 1

SMDIR=${SMDIR:-$HOME/smnormal}
CMDSTAN_LIB=$HOME/.cmdstan/cmdstan-2.40.0-rc0/stan/lib/stan_math/lib
BOOST=$CMDSTAN_LIB/boost_1.87.0
TBB=$CMDSTAN_LIB/tbb_2020.3/include
EIGEN=$SMDIR/lib/eigen_5.0.1

echo "=== node ==="
hostname
grep -m1 "model name" /proc/cpuinfo
echo "cores: $(nproc), governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"
echo "allowed CPUs: $(grep -m1 Cpus_allowed_list /proc/self/status)"

echo
echo "=== toolchain ==="
source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1

echo
echo "=== paths ==="
for d in "$SMDIR/stan" "$BOOST" "$TBB" "$EIGEN"; do
  [ -d "$d" ] || { echo "FATAL: missing $d"; exit 1; }
  echo "ok $d"
done

echo
echo "=== the erfcx header under test ==="
md5sum "$SMDIR/stan/math/prim/fun/erfcx.hpp"
grep -c "fma" "$SMDIR/stan/math/prim/fun/erfcx.hpp"

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-I $SMDIR -isystem $BOOST -isystem $EIGEN -isystem $TBB"

echo
echo "=== build ==="
g++ $FLAGS $INCS -o bench_normal_lcdf bench_normal_lcdf_kernels.cpp \
  || { echo "FATAL build"; exit 1; }
echo "build ok"

echo
for rep in 1 2 3; do
  echo "########## rep $rep ##########"
  ./bench_normal_lcdf 2000000
done

echo
echo "EXIT=0"
