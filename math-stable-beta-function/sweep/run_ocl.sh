#!/bin/bash
#SBATCH --job-name=beta_ocl
#SBATCH --output=beta_ocl_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=gpu-debug
#SBATCH --gres=gpu:1
# The OpenCL device function beta() of the branch against the CPU value
# and mpmath (PLAN §3.5 item 4): builds probe_ocl.cpp against the branch
# tree of beta_trees.tar.gz on the GPU node, runs it, and scores it with
# score_ocl.py.
#
# Deployed with this job: probe_ocl.cpp, grid.hpp, score_ocl.py,
# beta_trees.tar.gz.
# Environment: LIB, a Stan Math lib/ directory with a built tbb/.
# Do not use `set -u`: /etc/profile reads unset variables.
set -o pipefail

source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB="${LIB:?set LIB to a Stan Math lib/ directory}"
TBB="$LIB/tbb"
for f in probe_ocl.cpp grid.hpp score_ocl.py beta_trees.tar.gz; do
  [ -f "$f" ] || { echo "SETUP FAILED: no $f"; exit 1; }
done
echo "== host: $(hostname)   date: $(date -Is)   g++: $(g++ --version | head -1)"
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
md5sum probe_ocl.cpp grid.hpp score_ocl.py beta_trees.tar.gz

WORK="${TMPDIR:-/tmp}/beta_ocl_$SLURM_JOB_ID"
mkdir -p "$WORK" || exit 1
tar -xzf beta_trees.tar.gz -C "$WORK" branch || { echo "SETUP FAILED: untar"; exit 1; }
md5sum "$WORK"/branch/stan/math/opencl/kernels/device_functions/beta.hpp \
       "$WORK"/branch/stan/math/opencl/kernel_generator/elt_function_cl.hpp

g++ -std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS \
  -DSTAN_OPENCL -DOPENCL_DEVICE_ID=0 -DOPENCL_PLATFORM_ID=0 \
  -DCL_HPP_TARGET_OPENCL_VERSION=120 -DCL_HPP_MINIMUM_OPENCL_VERSION=120 \
  -DCL_HPP_ENABLE_EXCEPTIONS -DINTEGRATED_OPENCL=0 -Wno-ignored-attributes \
  -I "$WORK/branch" -I "$LIB/opencl_v2025.07.22" -I "$LIB/eigen_5.0.1" \
  -I "$LIB/boost_1.87.0" -I "$LIB/sundials_6.1.1/include" \
  -I "$LIB/tbb_2020.3/include" \
  probe_ocl.cpp -o "$WORK/probe_ocl" -L "$TBB" -ltbb -Wl,-rpath,"$TBB" \
  /usr/lib64/libOpenCL.so.1 > "$WORK/build.log" 2>&1
[ -x "$WORK/probe_ocl" ] || { echo "BUILD FAILED"; tail -40 "$WORK/build.log"; exit 1; }
"$WORK/probe_ocl" > ocl_branch.txt || { echo "RUN FAILED"; exit 1; }
wc -l ocl_branch.txt

module load scicomp-python-env
python3 score_ocl.py < ocl_branch.txt > score_ocl.txt || { echo "SCORE FAILED"; exit 1; }
cat score_ocl.txt
rm -rf "$WORK"
echo "EXIT=0"
