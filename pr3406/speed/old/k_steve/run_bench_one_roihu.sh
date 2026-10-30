#!/bin/bash
#SBATCH --account=project_2020226
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --job-name=s3406_one
#SBATCH --output=s3406_one_%j.log
# medium nodes rc4149, rc4150, rc4151 fail every job in the task prolog
#SBATCH --exclude=rc4149,rc4150,rc4151
#
# As run_bench_steve_roihu.sh, but one binary per function
# (bench_one.cpp -DBENCH_ONLY=<function>), so that the other functions in
# the translation unit cannot change GCC's inlining decisions.
# Output lines: RESULT <tree>_<flags> <function> ... (same format).
# Do not use `set -u` or `set -e`: /etc/profile reads unset variables.
export CSC_ENV_INIT_NON_INTERACTIVE=yes   # ROIHU.md gotcha R1
source /etc/profile
module purge
unset LD_PRELOAD                          # ROIHU.md gotcha R2
CMDSTAN_LIB=/projappl/project_2020226/cmdstan-2.40.0/stan/lib/stan_math/lib
# TREES, FUNCS (space- or '+'-separated) and TARBALL can be set with
# --export; the tarball holds t_<tree>/stan for each tree.
TREES=$(echo "${TREES:-e4745c87ff e8bb9a0b38 e1fea0651c cf177ef96b}" | tr '+' ' ')
SETS="default native"
FUNCS=$(echo "${FUNCS:-Phi std_normal_cdf normal_cdf std_normal_lcdf normal_lcdf lognormal_lcdf
lognormal_cdf exp_mod_normal_lcdf exp_mod_normal_lccdf exp_mod_normal_cdf
exp_mod_normal_lpdf skew_normal_lpdf}" | tr '+' ' ')
TARBALL=${TARBALL:-s3406_trees.tar.gz}
echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)   cores: $(nproc)"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"
cd "$SLURM_SUBMIT_DIR" || exit 1
for f in bench_one.cpp "$TARBALL" sundials_include_arkode.tar.gz; do
  [ -f "$f" ] || { echo "SETUP FAILED: missing $f"; exit 1; }
done
md5sum bench_one.cpp "$TARBALL"
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; exit 1; }
W="$TMPDIR/s3406one"
rm -rf "$W" && mkdir -p "$W/sun" || exit 1
tar -xzf "$TARBALL" -C "$W" || exit 1
tar -xzf sundials_include_arkode.tar.gz -C "$W/sun" || exit 1
BASE="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $CMDSTAN_LIB/eigen_5.0.1 -isystem $CMDSTAN_LIB/boost_1.87.0"
INCS="$INCS -isystem $CMDSTAN_LIB/tbb_2020.3/include -isystem $W/sun/include"
LINK="-L$CMDSTAN_LIB/tbb -Wl,-rpath,$CMDSTAN_LIB/tbb -ltbb"
set_flags() {
  case "$1" in
    default) echo "" ;;
    native) echo "-march=native -mtune=native" ;;
  esac
}
for t in $TREES; do
  for s in $SETS; do
    for fn in $FUNCS; do
      g++ $BASE $(set_flags $s) -DBENCH_ONLY=$fn -I "$W/t_$t" $INCS \
        -o "$W/b_${t}_${s}_$fn" bench_one.cpp $LINK > "$W/build_${t}_${s}_$fn.log" 2>&1 &
    done
    wait
  done
done
missing=0
for t in $TREES; do
  for s in $SETS; do
    for fn in $FUNCS; do
      [ -x "$W/b_${t}_${s}_$fn" ] || { echo "BUILD FAILED ${t}_${s}_$fn"; missing=$((missing + 1)); }
    done
  done
done
echo "== builds done $(date -Is), missing $missing"

for pass in 1 2; do
  for fn in $FUNCS; do
    for t in $TREES; do
      for s in $SETS; do
        [ -x "$W/b_${t}_${s}_$fn" ] || continue
        if [ "${PERF:-0}" = "1" ]; then
          # user-space instructions and cycles of the whole run
          perf stat -x, -e instructions:u,cycles:u -o "$W/perf.txt" \
            taskset -c 2 "$W/b_${t}_${s}_$fn" "${t}_$s" 3
          echo "PERF ${t}_$s $fn pass $pass $(grep -E 'instructions|cycles' "$W/perf.txt" | cut -d, -f1,3 | tr '\n' ' ')"
        else
          taskset -c 2 "$W/b_${t}_${s}_$fn" "${t}_$s" 3
        fi
      done
    done
  done
done
rm -rf "$W"
echo "EXIT=0"
