#!/bin/bash
#SBATCH --job-name=s3088
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --tmp=100G
#SBATCH --time=01:30:00
#SBATCH --output=run_3088_%j.log
# stan-dev/stan#3088.  Builds CmdStan develop 3ab132c (Stan a8ec3d18f,
# Math dc66941969, stanc3 v2.40.0) once, then runs run_3088.R.
# STEPS (default 1:2:3), N_SEEDS, S2_MODELS (colon-separated) pass through
# the environment.

set -e -o pipefail
source /etc/profile
module --ignore_cache load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

HERE=<home>/stan_3088_work
BASE=<scratch>/stan3088
T=$BASE/cmdstan
WORK=$BASE/work
MATHLIB=<scratch>/smtest_lib_d12f488749
STEPS=${STEPS:-1:2:3}
echo "== host $(hostname), $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
echo "== $(g++ --version | head -1)"

if [ ! -x "$T/bin/stansummary" ]; then
  echo "== build CmdStan in $T"
  rm -rf "$T" && mkdir -p "$T" && cd "$T"
  tar xzf "$HERE/cmdstan.tar.gz"
  rm -rf stan && mkdir -p stan/lib
  tar xzf "$HERE/stan.tar.gz" -C stan
  tar xzf "$HERE/math.tar.gz" -C stan/lib
  ln -s "$MATHLIB" stan/lib/stan_math/lib
  mkdir -p bin && cp /home/ave/cmdstan-2.40.0-v3/bin/stanc bin/stanc
  make -j16 build > build.log 2>&1 || { tail -30 build.log; exit 1; }
  echo "== build done"
fi
cd "$T"
bin/stanc --version
md5sum stan/src/stan/services/util/initialize.hpp \
  stan/lib/stan_math/stan/math/prim/constraint/offset_multiplier_constrain.hpp

mkdir -p "$WORK"
cp "$HERE"/*.stan "$HERE/run_3088.R" "$WORK/"
export CMDSTAN=$T
cd "$WORK"
Rscript run_3088.R "$WORK" ${STEPS//:/ }
