#!/bin/bash
#SBATCH --account=<project>
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:40:00
#SBATCH --job-name=p3408_bench
#SBATCH --output=p3408_bench_%j.log
#
# Bundle: example SLURM job of both speed tables of #3408 Part 3 (CSC
# Roihu job 2430495).  Paths are variables; set them with --export.
#
# Speed of the final #3408 tree against the baseline headers in
# variants/orig (made by make_tarball.sh with ORIG_REF=origin/develop:
# develop's version of every stan/ file that the PR modifies). bench_cdfs
# (var, value plus all gradients) and bench_lccdf_tail (gamma_lccdf, alpha
# derivative, var, fvar<double>, fvar<var>), 4 alternating runs on one
# whole node.
#
#   JOBDIR         directory with pr3408_fix.tar.gz (from make_tarball.sh)
#   STAN_MATH_LIB  lib/ of a Stan Math clone with a built lib/tbb
# Do not use `set -u` or `set -e`: /etc/profile reads unset variables.
set -o pipefail

export CSC_ENV_INIT_NON_INTERACTIVE=yes
[ -f /etc/profile ] && source /etc/profile
command -v module > /dev/null && module purge
unset LD_PRELOAD

JOBDIR="${JOBDIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
: "${STAN_MATH_LIB:?set STAN_MATH_LIB to a lib/ directory with a built tbb}"

echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)   cores: $(nproc)"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"

cd "$JOBDIR" || exit 1
md5sum pr3408_fix.tar.gz
W="${TMPDIR:-/tmp}/p3408b"
rm -rf "$W" && mkdir -p "$W" || exit 1
tar xzf pr3408_fix.tar.gz -C "$W" || { echo "SETUP FAILED: untar"; exit 1; }
cp -rp "$STAN_MATH_LIB" "$W/lib" || { echo "SETUP FAILED: copy lib"; exit 1; }
cp -r "$W/lib_overlay/lib/." "$W/lib/" || { echo "SETUP FAILED: overlay"; exit 1; }
cd "$W" || exit 1
echo "== baseline files: $(wc -l < variants/orig_files.txt)"

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS -DTBB_INTERFACE_NEW"
INCS="-isystem $W/lib/eigen_5.0.1 -isystem $W/lib/boost_1.87.0"
INCS="$INCS -isystem $W/lib/tbb_2020.3/include -isystem $W/lib/sundials_6.1.1/include"
LINK="-L$W/lib/tbb -ltbb -Wl,-rpath,$W/lib/tbb"
for b in bench_cdfs bench_lccdf_tail; do
  g++ $FLAGS -I variants/orig -I tree $INCS -o ${b}_orig $b.cpp $LINK \
    || { echo "FATAL build $b orig"; exit 1; }
  g++ $FLAGS -I tree $INCS -o ${b}_fix $b.cpp $LINK \
    || { echo "FATAL build $b fix"; exit 1; }
done
for run in 1 2 3 4; do
  taskset -c 2 ./bench_cdfs_orig orig
  taskset -c 2 ./bench_cdfs_fix fix
  taskset -c 2 ./bench_lccdf_tail_orig orig 1
  taskset -c 2 ./bench_lccdf_tail_fix fix 1
done

cd / && rm -rf "$W"
echo "EXIT=0"
