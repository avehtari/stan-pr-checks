#!/bin/bash
#SBATCH --job-name=emg_ocl_bench
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --output=emg_ocl_bench_%j.log

# GPU time of exp_mod_normal lcdf / lccdf / cdf (bench_ocl_emg.cpp) on
# develop e4745c87ff, the #3406 head 9b38ecbcb0 and the rewrite (t_emg).
# N = 1e6.  Builds in parallel, then 2 passes of 3 rounds.
# Do not use `set -u`: /etc/profile reads unset variables.

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB=/scratch/work/ave/smtest_lib_d12f488749
WORK=${TMPDIR:-/tmp}/emg_ocl_bench_$SLURM_JOB_ID
TREES="develop 9b38ecbcb0 emg"

hostname
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1
LIBCL=/usr/lib64/libOpenCL.so.1
[ -e "$LIBCL" ] || { echo "FATAL: no $LIBCL"; exit 1; }

mkdir -p "$WORK" && tar -xzf emg_ocl_trees.tar.gz -C "$WORK" || { echo "FATAL: untar"; exit 1; }
md5sum "$WORK"/t_*/stan/math/opencl/prim/exp_mod_normal_lcdf.hpp

FLAGS="-std=c++17 -O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS"
OCL="-DSTAN_OPENCL -DOPENCL_DEVICE_ID=0 -DOPENCL_PLATFORM_ID=0 -DINTEGRATED_OPENCL=0"
OCL="$OCL -DCL_HPP_TARGET_OPENCL_VERSION=120 -DCL_HPP_MINIMUM_OPENCL_VERSION=120"
OCL="$OCL -DCL_HPP_ENABLE_EXCEPTIONS -Wno-ignored-attributes -isystem $LIB/opencl_v2025.07.22"
INCS="-isystem $LIB/eigen_5.0.1 -isystem $LIB/boost_1.87.0"
INCS="$INCS -isystem $LIB/tbb_2020.3/include -isystem $LIB/sundials_6.1.1/include"
LINK="-L$LIB/tbb -Wl,-rpath,$LIB/tbb -ltbb $LIBCL"
for t in $TREES; do
  g++ $FLAGS $OCL -I "$WORK/t_$t" $INCS -o "$WORK/bench_$t" bench_ocl_emg.cpp $LINK \
    > "$WORK/build_$t.log" 2>&1 &
done
wait
for t in $TREES; do
  [ -x "$WORK/bench_$t" ] || { echo "FATAL build $t"; head -20 "$WORK/build_$t.log"; exit 1; }
done
echo "builds done $(date -Is)"

for pass in 1 2; do
  for t in $TREES; do
    "$WORK/bench_$t" "$t" 3 1000000 || echo "RUN FAILED $t"
  done
done
rm -rf "$WORK"
echo "EXIT=0"
