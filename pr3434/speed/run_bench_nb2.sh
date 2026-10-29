#!/bin/bash
#SBATCH --job-name=nb2_bench
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=nb2_bench_%j.log
#SBATCH --exclude=pht45

# neg_binomial_2_lpmf before (#3434 head 0377416157) and after evaluating
# log(mu + phi) once in the phi autodiff path.  bench_nb2.cpp; builds in
# parallel, then 3 passes x 3 rounds, alternated.
# Do not use `set -u`: /etc/profile reads unset variables.

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB="${LIB:?set LIB to a Stan Math lib/ directory}"
WORK=${TMPDIR:-/tmp}/nb2_bench_$SLURM_JOB_ID
TREES="head new"

hostname
grep -m1 "model name" /proc/cpuinfo
source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1

mkdir -p "$WORK" && tar -xzf nb2_trees.tar.gz -C "$WORK" || { echo "FATAL: untar"; exit 1; }
md5sum "$WORK"/t_*/stan/math/prim/prob/neg_binomial_2_lpmf.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
LINK="-L$LIB/tbb -Wl,-rpath,$LIB/tbb -ltbb"
for t in $TREES; do
  g++ $FLAGS -I "$WORK/t_$t" $INCS -o "$WORK/bench_$t" bench_nb2.cpp $LINK \
    > "$WORK/build_$t.log" 2>&1 &
done
wait
for t in $TREES; do
  [ -x "$WORK/bench_$t" ] || { echo "FATAL build $t"; head -20 "$WORK/build_$t.log"; exit 1; }
done
echo "builds done $(date -Is)"
for pass in 1 2 3; do
  for t in $TREES; do
    "$WORK/bench_$t" "$t" 3
  done
done
rm -rf "$WORK"
echo "EXIT=0"
