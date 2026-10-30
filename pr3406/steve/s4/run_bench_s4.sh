#!/bin/bash
#SBATCH --job-name=pr3406_s4
#SBATCH --partition=batch-csl
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=pr3406_s4_%j.out
#
# PR #3406, claim S4 (commit 2015557859): scalar exp_mod_normal_lpdf at
# z ~ 9.4.  One binary per tree, 4 interleaved passes (A B C A B C ...),
# each run pinned to one core.  Whole node (--exclusive).
# Do not use `set -u`: /etc/profile reads unset variables.
source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

LIB=/scratch/work/ave/smtest_lib_d12f488749
TREES="1f8b6ba94d 2015557859 e4745c87ff"
echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)   cores: $(nproc)"
echo "== governor: $(cat /sys/devices/system/cpu/cpu2/cpufreq/scaling_governor 2>/dev/null)"
echo "== MHz cpu2: $(grep MHz /proc/cpuinfo | sed -n 3p)"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"
cd "$SLURM_SUBMIT_DIR" || exit 1
for t in $TREES; do
  [ -f "stan_$t.tar.gz" ] || { echo "SETUP FAILED: missing stan_$t.tar.gz"; exit 1; }
done
md5sum bench_emg_scalar.cpp stan_*.tar.gz
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; exit 1; }
W="$TMPDIR/s4"
rm -rf "$W" && mkdir -p "$W" || exit 1
for t in $TREES; do
  tar -xzf "stan_$t.tar.gz" -C "$W" || exit 1
done
FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
INCS="-isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
LINK="-L$LIB/tbb -Wl,-rpath,$LIB/tbb -ltbb"
echo "== flags: $FLAGS"
for t in $TREES; do
  g++ $FLAGS -I "$W/t_$t" $INCS -o "$W/b_$t" bench_emg_scalar.cpp $LINK \
    > "$W/build_$t.log" 2>&1 &
done
wait
for t in $TREES; do
  [ -x "$W/b_$t" ] || { echo "BUILD FAILED $t"; cat "$W/build_$t.log"; exit 1; }
done
echo "== builds done $(date -Is)"
for pass in 1 2 3 4; do
  for t in $TREES; do
    echo "PASS $pass"
    taskset -c 2 "$W/b_$t" "$t"
  done
done
echo "== MHz cpu2 at end: $(grep MHz /proc/cpuinfo | sed -n 3p)"
rm -rf "$W"
echo "EXIT=0"
