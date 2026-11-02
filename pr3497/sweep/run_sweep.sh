#!/bin/bash
#SBATCH --job-name=beta_sweep
#SBATCH --output=beta_sweep_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
# Accuracy sweep of beta(a, b), develop against the branch (PLAN §3.5):
# builds probe_sweep.cpp against both trees of beta_trees.tar.gz
# (make_trees.sh), runs it, and scores both outputs with score_sweep.py.
#
# Deployed with this job: probe_sweep.cpp, grid.hpp, score_sweep.py,
# beta_trees.tar.gz.
# Environment: LIB, a Stan Math lib/ directory with a built tbb/.
# Do not use `set -u`: /etc/profile reads unset variables.
set -o pipefail

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB="${LIB:?set LIB to a Stan Math lib/ directory}"
TBB="$LIB/tbb"
for f in probe_sweep.cpp grid.hpp score_sweep.py beta_trees.tar.gz; do
  [ -f "$f" ] || { echo "SETUP FAILED: no $f"; exit 1; }
done
[ -f "$TBB/libtbb.so" ] || { echo "SETUP FAILED: no $TBB/libtbb.so"; exit 1; }
echo "== host: $(hostname)   date: $(date -Is)   g++: $(g++ --version | head -1)"
md5sum probe_sweep.cpp grid.hpp score_sweep.py beta_trees.tar.gz

WORK="${TMPDIR:-/tmp}/beta_sweep_$SLURM_JOB_ID"
mkdir -p "$WORK" || exit 1
tar -xzf beta_trees.tar.gz -C "$WORK" || { echo "SETUP FAILED: untar"; exit 1; }
md5sum "$WORK"/*/stan/math/prim/fun/beta.hpp

for t in develop branch; do
  g++ -std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS \
    -I "$WORK/$t" -I "$LIB/eigen_5.0.1" -I "$LIB/boost_1.87.0" \
    -I "$LIB/sundials_6.1.1/include" -I "$LIB/tbb_2020.3/include" \
    probe_sweep.cpp -o "$WORK/probe_$t" -L "$TBB" -ltbb \
    -Wl,-rpath,"$TBB" > "$WORK/build_$t.log" 2>&1 &
done
wait
for t in develop branch; do
  [ -x "$WORK/probe_$t" ] || { echo "BUILD FAILED $t"; cat "$WORK/build_$t.log"; exit 1; }
  "$WORK/probe_$t" > "sweep_$t.txt" || { echo "RUN FAILED $t"; exit 1; }
done
wc -l sweep_develop.txt sweep_branch.txt

module load scicomp-python-env
echo "== python: $(python3 --version)  mpmath: $(python3 -c 'import mpmath; print(mpmath.__version__)')"
python3 score_sweep.py sweep_develop.txt sweep_branch.txt "$SLURM_CPUS_PER_TASK" \
  > score_sweep.txt || { echo "SCORE FAILED"; exit 1; }
cat score_sweep.txt
rm -rf "$WORK"
echo "EXIT=0"
