#!/bin/bash
#SBATCH --job-name=st_bench
#SBATCH --output=st_bench_%j.log
#SBATCH --time=00:40:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
# Speed of student_t_lpdf, develop against the branch, scalar and
# vectorized calls (bench_st_vec.cpp). Reads the tree
# <tree> (CPU headers are never swapped there);
# the develop binary takes develop's student_t_lpdf.hpp from an overlay
# directory searched first. The two binaries run alternately ROUNDS times
# on one shared core; report the median ratio per kernel.
#
# Deployed with this job: bench_st_vec.cpp, dev_student_t_lpdf.hpp.
# Do not use `set -u`: /etc/profile reads unset variables.
set -o pipefail

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

cd "$SLURM_SUBMIT_DIR" || exit 1
TREE="${TREE:?set TREE to the Stan Math tree}"
TBB="${TBB:?set TBB to a directory with libtbb.so}"
REPS="${REPS:-200}"
ROUNDS="${ROUNDS:-5}"

echo "== host: $(hostname)   date: $(date -Is)   g++: $(g++ --version | head -1)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-)"
echo "== shared node: yes; ratios are interleaved"
for f in bench_st_vec.cpp dev_student_t_lpdf.hpp; do
  [ -f "$f" ] || { echo "SETUP FAILED: no $f"; exit 1; }
done
mkdir -p overlay/stan/math/prim/prob
cp dev_student_t_lpdf.hpp overlay/stan/math/prim/prob/student_t_lpdf.hpp
echo "== md5 of the headers used:"
md5sum overlay/stan/math/prim/prob/student_t_lpdf.hpp \
       "$TREE/stan/math/prim/prob/student_t_lpdf.hpp" \
       "$TREE/stan/math/prim/fun/digamma_diff.hpp" \
       "$TREE/stan/math/prim/fun/lbeta.hpp"

build() {
  label="$1"
  shift
  g++ -std=c++17 -O2 -DNDEBUG -pthread "$@" \
      -I "$TREE" -I "$TREE/lib/eigen_5.0.1" -I "$TREE/lib/boost_1.87.0" \
      -I "$TREE/lib/sundials_6.1.1/include" -I "$TREE/lib/tbb_2020.3/include" \
      bench_st_vec.cpp -o "bench_$label" -L "$TBB" -ltbb -Wl,-rpath,"$TBB" \
      2>&1 | tail -20
  [ -x "bench_$label" ] || { echo "BUILD FAILED $label"; exit 1; }
}
build develop -I overlay
build branch

for r in $(seq 1 "$ROUNDS"); do
  ./bench_develop "$REPS" | grep '^TIME' | sed "s/^TIME /TIME round$r develop /"
  ./bench_branch "$REPS" | grep '^TIME' | sed "s/^TIME /TIME round$r branch /"
done
echo "EXIT=0"
