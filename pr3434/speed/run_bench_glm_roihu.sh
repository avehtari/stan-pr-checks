#!/bin/bash
#SBATCH --account=project_2020226
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00
#SBATCH --job-name=glm_bench
#SBATCH --output=glm_bench_%j.log
# 2026-10-12: medium nodes rc4149, rc4150, rc4151 fail every job in the task
# prolog (mixtureis/ROIHU.md, "Traps").
#SBATCH --exclude=rc4149,rc4150,rc4151
#
# Roihu version of run_bench_glm.sh: speed of PR #3434 against develop,
# bench_glm.cpp (the negative binomial GLM) and bench_dd.cpp (the rows of
# the PR speed table), each built against three stan/ trees from
# glm_trees.tar.gz:
#   dev   develop d12f488749
#   pre   #3434 at 468ff26dbb (before the review fixes)
#   head  #3434 at 1c8afb741a (current)
# One whole medium node (node-exclusive). Builds run in parallel; then
# ROUNDS interleaved rounds (dev pre head), one binary at a time.
# Header libraries from the branch lib/ (glm_libs.tar.gz); the built TBB
# 2020.3 from CmdStan 2.40.0. Everything is unpacked in the node-local
# $TMPDIR. Flags as in dd/run_bench_dd.sh (-O2 -DNDEBUG).
#
# Deployed with this job: bench_glm.cpp, bench_dd.cpp, glm_trees.tar.gz,
# glm_libs.tar.gz.
# Do not use `set -u` or `set -e`: /etc/profile reads unset variables.
set -o pipefail

export CSC_ENV_INIT_NON_INTERACTIVE=yes   # ROIHU.md gotcha R1
source /etc/profile
module purge
unset LD_PRELOAD                          # ROIHU.md gotcha R2

cd "$SLURM_SUBMIT_DIR" || exit 1
TBB="${TBB:?set TBB to a directory with libtbb.so}"
TREES="${TREES:-dev pre head}"
TREES="${TREES//:/ }"
ROUNDS="${ROUNDS:-5}"
REPS_GLM="${REPS_GLM:-4000}"
REPS_DD="${REPS_DD:-200}"

echo "== host: $(hostname)   date: $(date -Is)   g++: $(g++ --version | head -1)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-)   cores: $(nproc)"
for f in bench_glm.cpp bench_dd.cpp glm_trees.tar.gz glm_libs.tar.gz; do
  [ -f "$f" ] || { echo "SETUP FAILED: no $f"; exit 1; }
done
md5sum bench_glm.cpp bench_dd.cpp glm_trees.tar.gz glm_libs.tar.gz
[ -f "$TBB/libtbb.so" ] || { echo "SETUP FAILED: no built TBB in $TBB"; exit 1; }
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; exit 1; }

WORK="$TMPDIR/glm_bench_$SLURM_JOB_ID"
LIB="$WORK/lib"
mkdir -p "$WORK" "$LIB" || exit 1
tar -xzf glm_trees.tar.gz -C "$WORK" || { echo "SETUP FAILED: untar trees"; exit 1; }
tar -xzf glm_libs.tar.gz -C "$LIB" || { echo "SETUP FAILED: untar libs"; exit 1; }
echo "== md5 of the changed headers:"
md5sum "$WORK"/t_*/stan/math/prim/prob/neg_binomial_2_log_glm_lpmf.hpp \
       "$WORK"/t_*/stan/math/prim/prob/neg_binomial_2_lpmf.hpp

FLAGS="-std=c++17 -O2 -DNDEBUG -pthread -D_REENTRANT"
INCS="-I $LIB/eigen_5.0.1 -I $LIB/boost_1.87.0"
INCS="$INCS -I $LIB/sundials_6.1.1/include -I $LIB/tbb_2020.3/include"
LINK="-L$TBB -ltbb -Wl,-rpath,$TBB"
# develop has no digamma_diff; bench_dd.cpp times it only for the branch
for t in $TREES; do
  DD=""
  [ "$t" != dev ] && DD="-DHAS_DIGAMMA_DIFF"
  for b in glm dd; do
    g++ $FLAGS $DD -I "$WORK/t_$t" $INCS "bench_$b.cpp" -o "$WORK/${b}_$t" \
      $LINK > "$WORK/build_${b}_$t.log" 2>&1 &
  done
done
wait
LABELS=""
for t in $TREES; do
  ok=1
  for b in glm dd; do
    if [ ! -x "$WORK/${b}_$t" ]; then
      echo "BUILD FAILED ${b}_$t"
      head -20 "$WORK/build_${b}_$t.log"
      ok=0
    fi
  done
  [ "$ok" -eq 1 ] && LABELS="$LABELS $t"
done
echo "== built:$LABELS   $(date -Is)"
[ "$LABELS" = " $TREES" ] || { echo "SETUP FAILED: not all trees built"; exit 1; }

for r in $(seq 1 "$ROUNDS"); do
  for t in $LABELS; do
    "$WORK/glm_$t" "$REPS_GLM" | sed "s/^/round$r $t /"
    "$WORK/dd_$t" "$REPS_DD" | grep '^TIME' | sed "s/^/round$r $t /"
  done
done
echo "== load: $(cat /proc/loadavg)"
rm -rf "$WORK"
echo "EXIT=0"
