#!/bin/bash
#SBATCH --account=project_2020226
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --job-name=p3406_speed
#SBATCH --output=p3406_speed_%j.log
#SBATCH --exclude=rc4149,rc4150,rc4151
#
# Rerun of the PR-text timing tables of stan-dev/math#3406 on one whole
# Roihu medium node: bench_one.cpp, one binary per function and tree,
# PASSES passes of (develop, head) per function, 2 rounds per call.
# Inputs in the submit directory: bench_one.cpp, build.sh (the pr3406
# build.sh) and TARBALL (t_<commit>/stan per tree, lib/sundials_6.1.1).
# CMDSTAN_LIB (--export) is a lib/ directory with eigen_5.0.1,
# boost_1.87.0, tbb_2020.3 and a built tbb/ (CmdStan 2.40.0).
# Do not use `set -u`: /etc/profile reads unset variables.
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile
module purge
unset LD_PRELOAD
: "${CMDSTAN_LIB:?set CMDSTAN_LIB with --export}"
TREES=$(echo "${TREES:-dc66941969 16c1ca18bf}" | tr '+' ' ')
FUNCS=$(echo "${FUNCS:-normal_lcdf normal_lccdf normal_cdf std_normal_lcdf
std_normal_lccdf std_normal_cdf Phi lognormal_lcdf lognormal_lccdf lognormal_cdf
exp_mod_normal_lpdf exp_mod_normal_lcdf exp_mod_normal_lccdf exp_mod_normal_cdf
skew_normal_lpdf}" | tr '+' ' ')
PASSES=${PASSES:-3}
ROUNDS=${ROUNDS:-2}
TARBALL=${TARBALL:-pr3406_speed_trees.tar.gz}
echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)   cores: $(nproc)"
echo "== governor: $(cat /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor 2>/dev/null)"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"
echo "== glibc: $(ldd --version | head -1)"
cd "$SLURM_SUBMIT_DIR" || exit 1
for f in bench_one.cpp build.sh "$TARBALL"; do
  [ -f "$f" ] || { echo "SETUP FAILED: missing $f"; exit 1; }
done
md5sum bench_one.cpp build.sh "$TARBALL"
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; exit 1; }
W="$TMPDIR/p3406speed"
rm -rf "$W" && mkdir -p "$W" || exit 1
tar -xzf "$TARBALL" -C "$W" || exit 1
for d in eigen_5.0.1 boost_1.87.0 tbb_2020.3 tbb; do
  [ -e "$CMDSTAN_LIB/$d" ] || { echo "SETUP FAILED: no $CMDSTAN_LIB/$d"; exit 1; }
  ln -s "$CMDSTAN_LIB/$d" "$W/lib/$d"
done
for t in $TREES; do
  echo "== tree $t md5 of stan/: $(cd "$W/t_$t" && find stan -type f | sort | xargs md5sum | md5sum | cut -d' ' -f1)"
done
echo "== flags: -std=c++17 -O3 -pthread -DBOOST_DISABLE_ASSERTS (build.sh default, no -march)"
for t in $TREES; do
  for fn in $FUNCS; do
    STAN_MATH="$W/t_$t" STAN_MATH_LIB="$W/lib" EXTRA_CXXFLAGS="-DBENCH_ONLY=$fn" \
      bash ./build.sh bench_one.cpp "$W/b_${t}_$fn" > "$W/build_${t}_$fn.log" 2>&1 &
  done
done
wait
missing=0
for t in $TREES; do
  for fn in $FUNCS; do
    if [ ! -x "$W/b_${t}_$fn" ]; then
      echo "BUILD FAILED ${t}_$fn"
      grep -m 3 'error' "$W/build_${t}_$fn.log" | cut -c1-300
      missing=$((missing + 1))
    fi
  done
done
echo "== builds done $(date -Is), missing $missing"

for pass in $(seq 1 "$PASSES"); do
  echo "== pass $pass $(date -Is)"
  for fn in $FUNCS; do
    for t in $TREES; do
      [ -x "$W/b_${t}_$fn" ] || continue
      taskset -c 2 "$W/b_${t}_$fn" "${t}_p$pass" "$ROUNDS"
      rc=$?
      [ "$rc" -eq 0 ] || echo "RUN FAILED ${t}_$fn pass $pass exit $rc"
    done
  done
done
rm -rf "$W"
echo "== end $(date -Is)"
echo "EXIT=0"
