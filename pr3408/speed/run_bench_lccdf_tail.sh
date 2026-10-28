#!/bin/bash
#SBATCH --job-name=p3408_bench
#SBATCH --time=00:40:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=p3408_bench_%j.log
#
# Bundle: example SLURM job for bench_lccdf_tail.cpp alone (Triton job
# 21111830; an intermediate comparison, not a PR table).  Paths are
# variables; set them with --export.
#
# Speed of gamma_lccdf value plus alpha derivative: baseline headers in
# variants/orig of pr3408_fix.tar.gz (make_tarball.sh) against the tree.
# Stan default flags (-O3). The binaries alternate, 4 runs x 1 round.
#
#   STAN_MATH_LIB  lib/ of a Stan Math clone with a built lib/tbb
#   MODULES        optional; modules to load (Triton: "triton/2025.1-gcc")
# Do not use `set -u`: /etc/profile reads unset variables.

cd "${SLURM_SUBMIT_DIR:-$(pwd)}" || exit 1
LIB="${STAN_MATH_LIB:?set STAN_MATH_LIB to a lib/ directory with a built tbb}"
WORK=${TMPDIR:-/tmp}/p3408_bench_${SLURM_JOB_ID:-$$}

hostname
grep -m1 "model name" /proc/cpuinfo
echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"
[ -f /etc/profile ] && source /etc/profile
[ -n "$MODULES" ] && module load $MODULES
g++ --version | head -1

mkdir -p "$WORK" || exit 1
tar -xzf pr3408_fix.tar.gz -C "$WORK" || { echo "FATAL: untar"; exit 1; }
cd "$WORK" || exit 1
md5sum tree/stan/math/prim/fun/log_gamma_q_dgamma.hpp \
  variants/orig/stan/math/prim/fun/log_gamma_q_dgamma.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS -DTBB_INTERFACE_NEW"
INCS="-isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
LINK="-L$LIB/tbb -ltbb -Wl,-rpath,$LIB/tbb"
for v in orig fix; do
  case $v in
    fix) PRE="" ;;
    *) PRE="-I variants/$v" ;;
  esac
  g++ $FLAGS $PRE -I tree $INCS -o "bench_$v" bench_lccdf_tail.cpp $LINK \
    || { echo "FATAL build $v"; exit 1; }
done

for run in 1 2 3 4; do
  ./bench_orig orig 1
  ./bench_fix fix 1
done

cd / && rm -rf "$WORK"
echo "EXIT=0"
