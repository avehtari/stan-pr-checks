#!/bin/bash
#SBATCH --job-name=s3088pr
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=run_probe_%j.log
# Probe of initialize.hpp (develop and PR) for the stan#3088 PR claims.
# Stan tree: <scratch>/stantest_3088 (stan_math_deploy.sh), with
# the develop header in .develop/.

set -e -o pipefail
source /etc/profile
module --ignore_cache load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

HERE=<home>/stan_3088_work
export STAN=<scratch>/stantest_3088
export STANC=$STAN/bin/stanc
W=<scratch>/stan3088/probe
echo "== host $(hostname), $(g++ --version | head -1)"
cat "$STAN/.deploy_info"
md5sum "$STAN/src/stan/services/util/initialize.hpp" \
  "$STAN/.develop/src/stan/services/util/initialize.hpp"
rm -rf "$W" && mkdir -p "$W/probe"
cp "$HERE"/probe_init.cpp "$HERE"/build.sh "$HERE"/compare_outputs.py \
  "$HERE"/partial_init_dependent.stan "$W/probe/"
cd "$W/probe"
"$STANC" --version
./build.sh probe_develop "$STAN/.develop/src"
./build.sh probe_pr
../out/probe_develop > ../out/probe_develop.txt
../out/probe_pr > ../out/probe_pr.txt
python3 compare_outputs.py ../out/probe_develop.txt ../out/probe_pr.txt \
  | tee ../out/compare.txt
echo EXIT=0
