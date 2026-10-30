#!/bin/bash
#SBATCH --account=project_2020226
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00
#SBATCH --job-name=prtab_roihu
#SBATCH --output=prtab_roihu_%j.log
# 2026-10-12: medium nodes rc4149, rc4150, rc4151 fail every job in the task
# prolog (mixtureis/ROIHU.md, "Traps").
#SBATCH --exclude=rc4149,rc4150,rc4151
#
# The #3406 PR-text timing table (bench_pr_table.cpp, full mode) on a whole
# Roihu medium node: develop 5252d51d47, PR at creation a28c1c1637, develop
# e4745c87ff, head f4cd10ee4a.  Same as run_bench_pr_table.sh on Triton,
# with the system g++ and the CmdStan 2.40.0 libraries; the sundials
# headers come with the job (ARKODE).  Built in the node-local $TMPDIR.
# Do not use `set -u` or `set -e`: /etc/profile reads unset variables.
export CSC_ENV_INIT_NON_INTERACTIVE=yes   # ROIHU.md gotcha R1
source /etc/profile
module purge
unset LD_PRELOAD                          # ROIHU.md gotcha R2
CMDSTAN_LIB=/projappl/project_2020226/cmdstan-2.40.0/stan/lib/stan_math/lib

echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)   cores: $(nproc)"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"
echo "== glibc: $(ldd --version | head -1)"
cd "$SLURM_SUBMIT_DIR" || exit 1
INPUTS="bench_pr_table.cpp develop_5252d51d47_stan.tar.gz pr3406_a28c1c1637_stan.tar.gz
develop_e4745c87ff_stan.tar.gz commit_f4cd10ee4a.tar.gz sundials_include_arkode.tar.gz"
for f in $INPUTS; do
  [ -f "$f" ] || { echo "SETUP FAILED: missing $f"; exit 1; }
done
[ -f "$CMDSTAN_LIB/tbb/libtbb.so" ] || { echo "SETUP FAILED: no built TBB"; exit 1; }
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; exit 1; }

W="$TMPDIR/prtab"
rm -rf "$W" && mkdir -p "$W/dev0916" "$W/pr0918" "$W/devnow" "$W/head" "$W/sun" || exit 1
tar -xzf develop_5252d51d47_stan.tar.gz -C "$W/dev0916" || exit 1
tar -xzf pr3406_a28c1c1637_stan.tar.gz -C "$W/pr0918" || exit 1
tar -xzf develop_e4745c87ff_stan.tar.gz -C "$W/devnow" || exit 1
tar -xzf commit_f4cd10ee4a.tar.gz -C "$W/head" stan || exit 1
tar -xzf sundials_include_arkode.tar.gz -C "$W/sun" || exit 1
md5sum "$W"/*/stan/math/prim/fun/Phi.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $CMDSTAN_LIB/eigen_5.0.1 -isystem $CMDSTAN_LIB/boost_1.87.0"
INCS="$INCS -isystem $CMDSTAN_LIB/tbb_2020.3/include -isystem $W/sun/include"
LINK="-L$CMDSTAN_LIB/tbb -Wl,-rpath,$CMDSTAN_LIB/tbb -ltbb"
for t in dev0916 pr0918 devnow head; do
  g++ $FLAGS -I "$W/$t" $INCS -o "$W/bench_$t" bench_pr_table.cpp $LINK \
    > "$W/build_$t.log" 2>&1 &
done
wait
for t in dev0916 pr0918 devnow head; do
  [ -x "$W/bench_$t" ] || { echo "BUILD FAILED $t"; head -20 "$W/build_$t.log"; exit 1; }
done
echo "== builds done $(date -Is)"

for run in 1 2; do
  for t in dev0916 pr0918 devnow head; do
    "$W/bench_$t" "$t" 2
  done
done
rm -rf "$W"
echo "EXIT=0"
