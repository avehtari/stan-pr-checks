#!/bin/bash
#SBATCH --account=project_2020226
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:40:00
#SBATCH --job-name=pr3409_steve
#SBATCH --output=pr3409_steve_%j.log
# medium nodes rc4149, rc4150, rc4151 fail every job in the task prolog
#SBATCH --exclude=rc4149,rc4150,rc4151
#
# Speed and accuracy of the incomplete beta roots on three trees:
#   dev    develop e4745c87ff (the develop merged into the branch on 10-20)
#   ours   bd8f27a265 (stable-inc-beta merged with that develop)
#   steve  d2eeacd676 (the reverse-mode / one-sided refactor on top of ours)
# Order: extract and verify, build in parallel, speed rounds one binary at a
# time on one pinned core (dev ours steve, twice), then the accuracy probes
# in parallel. Results are copied back to $HOME/pr3409_steve.
# Do not use `set -u` or `set -e`: /etc/profile reads unset variables.
set -o pipefail

export CSC_ENV_INIT_NON_INTERACTIVE=yes   # ROIHU.md gotcha R1
source /etc/profile
module purge
unset LD_PRELOAD                          # ROIHU.md gotcha R2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

JOBDIR="$HOME/pr3409_steve"
CMDSTAN_LIB=/projappl/project_2020226/cmdstan-2.40.0/stan/lib/stan_math/lib
W="$TMPDIR/pr3409"

echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)   cores: $(nproc)"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"

rm -rf "$W" && mkdir -p "$W" || exit 1
cp -rp "$CMDSTAN_LIB" "$W/lib" || { echo "SETUP FAILED: copy lib"; exit 1; }
cd "$JOBDIR" || exit 1
# the CmdStan 2.40.0 release lacks sundials include/arkode, which develop has
tar xzf lib_overlay.tar.gz -C "$W" || { echo "SETUP FAILED: overlay"; exit 1; }
test -f "$W/lib/sundials_6.1.1/include/arkode/arkode.h" \
  || { echo "SETUP FAILED: arkode.h"; exit 1; }
cp probe_speed3.cpp probe_hess_beta.cpp probe_cmp_beta.cpp points_beta.txt \
  points_beta_sub.txt "$W/" || { echo "SETUP FAILED: copy sources"; exit 1; }

# ── extract and verify ──
declare -A EXPECT
EXPECT[dev]="99540d4e20cabfd85b9cbe36ea0c104c 112047e084536b025e89071b96048809 b739259381b5900e5bc38e7896ba131f"
EXPECT[ours]="29725873b4e7c486071c289a03a4bfb2 112047e084536b025e89071b96048809 b739259381b5900e5bc38e7896ba131f"
EXPECT[steve]="7acccade657b0c03b37597abe4346c97 daf5a4688e37b466c032efecb3dbe8d0 04013130773a989c63d5cbec2c244b1e"
FILES="stan/math/prim/fun/grad_reg_inc_beta.hpp stan/math/prim/prob/beta_lcdf.hpp stan/math/rev/fun/inc_beta.hpp"
for t in dev ours steve; do
  mkdir -p "$W/$t" || exit 1
  tar xzf "tree_$t.tar.gz" -C "$W/$t" || { echo "EXTRACT FAILED ($t)"; exit 1; }
  got=""
  for f in $FILES; do got="$got $(md5sum "$W/$t/$f" | cut -d' ' -f1)"; done
  got="${got# }"
  if [ "$got" != "${EXPECT[$t]}" ]; then
    echo "MD5 MISMATCH ($t): $got"
    exit 1
  fi
  echo "== tree $t verified"
done
cd "$W" || exit 1

# ── build ──
FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS -DTBB_INTERFACE_NEW"
INCS="-isystem $W/lib/eigen_5.0.1 -isystem $W/lib/boost_1.87.0"
INCS="$INCS -isystem $W/lib/tbb_2020.3/include -isystem $W/lib/sundials_6.1.1/include"
LINK="-L$W/lib/tbb -ltbb -Wl,-rpath,$W/lib/tbb"
echo "== flags: $FLAGS"
build() {
  tag="$1"; src="$2"; tree="$3"; extra="$4"
  bin="${src%.cpp}_${tag}"
  if g++ $FLAGS $extra -I "$W/$tree" $INCS "$src" -o "$bin" $LINK \
      > "build_${bin}.log" 2>&1; then
    echo "built $bin"
  else
    echo "BUILD FAILED $bin"
    grep -m 20 -E "error" "build_${bin}.log"
  fi
}
build dev probe_speed3.cpp dev "" &
build ours probe_speed3.cpp ours "" &
build steve probe_speed3.cpp steve "-DHAVE_GRAD_RIB" &
build ours probe_cmp_beta.cpp ours "" &
build steve probe_cmp_beta.cpp steve "" &
build ours probe_hess_beta.cpp ours "" &
build steve probe_hess_beta.cpp steve "" &
wait
echo "== builds done: $(date -Is)"

# ── speed, one binary at a time, pinned ──
for round in 1 2; do
  for t in dev ours steve; do
    echo ""
    echo "######## speed $t, round $round ########"
    if [ -x "probe_speed3_$t" ]; then
      timeout 900 taskset -c 2 "./probe_speed3_$t"
      echo "== exit $?"
    else
      echo "== missing binary"
    fi
  done
done
echo "== speed done: $(date -Is)"

# ── accuracy, in parallel ──
run() {
  bin="$1"; shift; out="$1"; shift
  if [ -x "$bin" ]; then
    # Roihu has no /usr/bin/time
    timeout 900 "./$bin" "$@" > "$out" 2>&1
    echo "$bin exit $? rows $(grep -vc '^#' "$out") md5 $(md5sum "$out" | cut -d' ' -f1)"
  else
    echo "$bin missing"
  fi
}
run probe_cmp_beta_ours cmp_double_ours.txt points_beta.txt double &
run probe_cmp_beta_steve cmp_double_steve.txt points_beta.txt double &
run probe_cmp_beta_ours cmp_all_ours.txt points_beta_sub.txt all &
run probe_cmp_beta_steve cmp_all_steve.txt points_beta_sub.txt all &
run probe_hess_beta_ours hess_ours.txt points_beta_sub.txt &
run probe_hess_beta_steve hess_steve.txt points_beta_sub.txt &
wait

cp cmp_double_*.txt cmp_all_*.txt hess_*.txt build_*.log "$JOBDIR/" \
  || echo "COPY BACK FAILED"
cd / && rm -rf "$W"
echo "== done: $(date -Is)"
echo "EXIT=0"
