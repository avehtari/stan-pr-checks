#!/bin/bash
#SBATCH --job-name=st_comp
#SBATCH --output=st_comp_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
# Cost of the parts of the student_t_lpdf nu terms (bench_comp.cpp), to
# find the cause of the branch slowdown of job 20815014. Reads the tree
# <tree> (CPU headers are never swapped there);
# the develop binary takes develop's student_t_lpdf.hpp from an overlay
# directory searched first. The two binaries run alternately ROUNDS times
# on one shared core; report the median per kernel and range.
#
# Deployed with this job: bench_comp.cpp, dev_student_t_lpdf.hpp.
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
echo "== shared node: yes; rounds are interleaved"
for f in bench_comp.cpp dev_student_t_lpdf.hpp; do
  [ -f "$f" ] || { echo "SETUP FAILED: no $f"; exit 1; }
done
mkdir -p overlay_comp/stan/math/prim/prob
cp dev_student_t_lpdf.hpp overlay_comp/stan/math/prim/prob/student_t_lpdf.hpp
echo "== md5 of the headers used:"
md5sum overlay_comp/stan/math/prim/prob/student_t_lpdf.hpp \
       "$TREE/stan/math/prim/prob/student_t_lpdf.hpp" \
       "$TREE/stan/math/prim/fun/digamma_diff.hpp" \
       "$TREE/stan/math/prim/fun/lbeta.hpp" \
       "$TREE/stan/math/prim/fun/lgamma_stirling_diff.hpp" \
       "$TREE/stan/math/prim/fun/lgamma.hpp" \
       "$TREE/stan/math/prim/fun/digamma.hpp"

build() {
  label="$1"
  shift
  g++ -std=c++17 -O2 -DNDEBUG -pthread "$@" \
      -I "$TREE" -I "$TREE/lib/eigen_5.0.1" -I "$TREE/lib/boost_1.87.0" \
      -I "$TREE/lib/sundials_6.1.1/include" -I "$TREE/lib/tbb_2020.3/include" \
      bench_comp.cpp -o "comp_$label" -L "$TBB" -ltbb -Wl,-rpath,"$TBB" \
      2>&1 | tail -20
  [ -x "comp_$label" ] || { echo "BUILD FAILED $label"; exit 1; }
}
rm -f comp_develop comp_branch
build develop -I overlay_comp
build branch

for r in $(seq 1 "$ROUNDS"); do
  ./comp_develop "$REPS" | grep '^TIME' | sed "s/^TIME /TIME round$r develop /"
  ./comp_branch "$REPS" | grep '^TIME' | sed "s/^TIME /TIME round$r branch /"
done
echo "EXIT=0"
