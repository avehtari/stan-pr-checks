#!/bin/bash
#SBATCH --account=project_2020226
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:40:00
#SBATCH --job-name=prog_roihu
#SBATCH --output=prog_roihu_%j.log
# 2026-10-12: medium nodes rc4149, rc4150, rc4151 fail every job in the task
# prolog (mixtureis/ROIHU.md, "Traps").
#SBATCH --exclude=rc4149,rc4150,rc4151
#
# The #3406 progression in one run: develop e4745c87ff, the PR as opened
# a28c1c1637, the PR before our commits 686fc67261, and the PR now
# 9b38ecbcb0.  bench_pr_table.cpp full mode (10 functions; N = 100 to
# 100000; z on [-3, 3] and [-30, 30]; double and Matrix<var>).  Whole
# medium node, system g++, CmdStan 2.40.0 libraries, sundials headers from
# the job.  Builds in parallel, then 2 passes over the trees, 2 rounds each.
# Do not use `set -u` or `set -e`: /etc/profile reads unset variables.
export CSC_ENV_INIT_NON_INTERACTIVE=yes   # ROIHU.md gotcha R1
source /etc/profile
module purge
unset LD_PRELOAD                          # ROIHU.md gotcha R2
CMDSTAN_LIB=/projappl/project_2020226/cmdstan-2.40.0/stan/lib/stan_math/lib
TREES="develop a28c1c1637 686fc67261 9b38ecbcb0"

echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)   cores: $(nproc)"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"
cd "$SLURM_SUBMIT_DIR" || exit 1
for f in bench_pr_table.cpp progression_trees.tar.gz sundials_include_arkode.tar.gz; do
  [ -f "$f" ] || { echo "SETUP FAILED: missing $f"; exit 1; }
done
[ -f "$CMDSTAN_LIB/tbb/libtbb.so" ] || { echo "SETUP FAILED: no built TBB"; exit 1; }
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; exit 1; }

W="$TMPDIR/prog"
rm -rf "$W" && mkdir -p "$W/sun" || exit 1
tar -xzf progression_trees.tar.gz -C "$W" || exit 1
tar -xzf sundials_include_arkode.tar.gz -C "$W/sun" || exit 1
md5sum "$W"/t_*/stan/math/prim/fun/Phi.hpp "$W"/t_*/stan/math/prim/prob/exp_mod_normal_lcdf.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $CMDSTAN_LIB/eigen_5.0.1 -isystem $CMDSTAN_LIB/boost_1.87.0"
INCS="$INCS -isystem $CMDSTAN_LIB/tbb_2020.3/include -isystem $W/sun/include"
LINK="-L$CMDSTAN_LIB/tbb -Wl,-rpath,$CMDSTAN_LIB/tbb -ltbb"
for t in $TREES; do
  g++ $FLAGS -I "$W/t_$t" $INCS -o "$W/bench_$t" bench_pr_table.cpp $LINK \
    > "$W/build_$t.log" 2>&1 &
done
wait
for t in $TREES; do
  [ -x "$W/bench_$t" ] || { echo "BUILD FAILED $t"; head -20 "$W/build_$t.log"; exit 1; }
done
echo "== builds done $(date -Is)"

for pass in 1 2; do
  for t in $TREES; do
    "$W/bench_$t" "$t" 2
  done
done
rm -rf "$W"
echo "EXIT=0"
