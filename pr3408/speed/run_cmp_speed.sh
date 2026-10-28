#!/bin/bash
#SBATCH --job-name=cmp_speed
#SBATCH --output=cmp_speed_%j.log
#SBATCH --time=00:45:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#
# Bundle: example SLURM job of the root speed table of #3408 Part 1
# (Triton job 20356756).  Paths are variables; set them with --export.
#
# Axis (c). Build probe_cmp_speed.cpp against both trees and run both on
# the same exclusive node, back to back, twice each to show the run-to-run
# spread.
#
#   DEV            develop tree (stan/)
#   BRA            PR tree (stan/)
#   STAN_MATH_LIB  lib/ of a Stan Math clone with a built lib/tbb
#   MODULES        optional; modules to load (Triton: "triton/2025.1-gcc")
set -o pipefail

[ -f /etc/profile ] && source /etc/profile
[ -n "$MODULES" ] && module load $MODULES

: "${DEV:?set DEV to the develop tree}"
: "${BRA:?set BRA to the PR tree}"
LIB="${STAN_MATH_LIB:?set STAN_MATH_LIB to a lib/ directory with a built tbb}"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1

echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu:  $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-)"
echo "== governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"

build() {
  tag="$1"; root="$2"
  rm -f "probe_cmp_speed_${tag}"
  g++ -std=c++17 -O2 -DNDEBUG -pthread \
      -I "$root" -I "$LIB/eigen_5.0.1" -I "$LIB/boost_1.87.0" \
      -I "$LIB/sundials_6.1.1/include" -I "$LIB/tbb_2020.3/include" \
      probe_cmp_speed.cpp -o "probe_cmp_speed_${tag}" \
      -L "$LIB/tbb" -ltbb -Wl,-rpath,"$LIB/tbb" \
      > "cmp_speed_build_${tag}.log" 2>&1 || { echo "BUILD FAILED ($tag)"; tail -30 "cmp_speed_build_${tag}.log"; return 1; }
  echo "built probe_cmp_speed_${tag}"
}

build develop "$DEV" || exit 1
build branch "$BRA" || exit 1

for round in 1 2; do
  for tag in develop branch; do
    echo ""
    echo "######## $tag, round $round ########"
    ./probe_cmp_speed_${tag}
  done
done
echo "== done: $(date -Is)"
