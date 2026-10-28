#!/bin/bash
#SBATCH --job-name=pr3409_speed
#SBATCH --output=pr3409_speed_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#
# Speed tables of stan-dev/math#3409: build probe_cmp_beta_speed.cpp against
# the develop tree and the PR tree and run both on the same exclusive node,
# back to back, twice each to show the run-to-run spread.
#
#   TREES          tarball with dev/stan and pr/stan (make_trees.sh)
#   STAN_MATH_LIB  lib/ of a Stan Math clone with a built lib/tbb; the same
#                  lib/ serves both trees
#   MODULES        optional; modules to load, separated by ':'
#                  (Triton: "triton/2025.1-gcc:scicomp-r-env", which gives
#                  the conda-forge g++)
set -o pipefail

[ -f /etc/profile ] && source /etc/profile
[ -n "$MODULES" ] && module --ignore_cache load ${MODULES//:/ }

: "${TREES:?set TREES to the tarball from make_trees.sh}"
LIB="${STAN_MATH_LIB:?set STAN_MATH_LIB to a lib/ directory with a built tbb}"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1

echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu:  $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-)"
echo "== governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
echo "== compiler: $(g++ --version | head -1)"
rm -rf trees && mkdir trees && tar -xzf "$TREES" -C trees || exit 1
cat trees/commits.txt
for f in stan/math/prim/fun/grad_reg_inc_beta.hpp \
         stan/math/prim/fun/inc_beta_dda.hpp; do
  printf "== md5 dev %s  pr %s  %s\n" \
    "$(md5sum "trees/dev/$f" | cut -d' ' -f1)" \
    "$(md5sum "trees/pr/$f" | cut -d' ' -f1)" "$f"
done

build() {
  tag="$1"
  rm -f "probe_speed_${tag}"
  g++ -std=c++17 -O3 -DNDEBUG -pthread \
      -I "trees/$tag" -I "$LIB/eigen_5.0.1" -I "$LIB/boost_1.87.0" \
      -I "$LIB/sundials_6.1.1/include" -I "$LIB/tbb_2020.3/include" \
      probe_cmp_beta_speed.cpp -o "probe_speed_${tag}" \
      -L "$LIB/tbb" -ltbb -Wl,-rpath,"$LIB/tbb" \
      > "speed_build_${tag}.log" 2>&1 \
    || { echo "BUILD FAILED ($tag)"; tail -30 "speed_build_${tag}.log"; return 1; }
  echo "built probe_speed_${tag}"
}

build dev || exit 1
build pr || exit 1

for round in 1 2; do
  for tag in dev pr; do
    echo ""
    echo "######## $tag, round $round ########"
    ./probe_speed_${tag}
  done
done
echo "== done: $(date -Is)"
