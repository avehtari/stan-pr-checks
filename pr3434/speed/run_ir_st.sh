#!/bin/bash
#SBATCH --job-name=st_ir
#SBATCH --output=st_ir_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
# Instructions per student_t_lpdf call (ir_st.cpp) under valgrind
# cachegrind, develop against the branch, per nu range. Reads the tree
# <tree>; the develop binary takes develop's
# student_t_lpdf.hpp from an overlay directory searched first. Same compile
# flags as run_bench_st.sh.
#
# Deployed with this job: ir_st.cpp, dev_student_t_lpdf.hpp.
# Do not use `set -u`: /etc/profile reads unset variables.
set -o pipefail

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

cd "$SLURM_SUBMIT_DIR" || exit 1
TREE="${TREE:?set TREE to the Stan Math tree}"
TBB="${TBB:?set TBB to a directory with libtbb.so}"
CALLS="${CALLS:-20000}"

echo "== host: $(hostname)   date: $(date -Is)   g++: $(g++ --version | head -1)"
for f in ir_st.cpp dev_student_t_lpdf.hpp; do
  [ -f "$f" ] || { echo "SETUP FAILED: no $f"; exit 1; }
done
mkdir -p overlay_ir/stan/math/prim/prob
cp dev_student_t_lpdf.hpp overlay_ir/stan/math/prim/prob/student_t_lpdf.hpp
echo "== md5 of the headers used:"
md5sum overlay_ir/stan/math/prim/prob/student_t_lpdf.hpp \
       "$TREE/stan/math/prim/prob/student_t_lpdf.hpp"

build() {
  label="$1"
  shift
  g++ -std=c++17 -O2 -DNDEBUG -pthread "$@" \
      -I "$TREE" -I "$TREE/lib/eigen_5.0.1" -I "$TREE/lib/boost_1.87.0" \
      -I "$TREE/lib/sundials_6.1.1/include" -I "$TREE/lib/tbb_2020.3/include" \
      ir_st.cpp -o "ir_$label" -L "$TBB" -ltbb -Wl,-rpath,"$TBB" \
      2>&1 | tail -20
  [ -x "ir_$label" ] || { echo "BUILD FAILED $label"; exit 1; }
}
rm -f ir_develop ir_branch
build develop -I overlay_ir
build branch

irefs() {
  valgrind --tool=cachegrind --cache-sim=no --cachegrind-out-file=/dev/null \
    "$@" 2>&1 | awk '/I *refs:/ {gsub(",", "", $NF); print $NF}'
}

for range in A B C; do
  for kernel in value grad; do
    for label in develop branch; do
      i0=$(irefs "./ir_$label" "$range" "$kernel" 0)
      i1=$(irefs "./ir_$label" "$range" "$kernel" "$CALLS")
      per=$(awk -v a="$i0" -v b="$i1" -v n="$CALLS" 'BEGIN {printf "%.1f", (b - a) / n}')
      echo "IR $range $kernel $label $per"
    done
  done
done
echo "EXIT=0"
