#!/bin/bash
#SBATCH --job-name=ocl_dev_bench
#SBATCH --partition=gpu-debug
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=ocl_dev_bench_%j.log

# Device-function time (bench_ocl_device.cpp, in-register walk, baseline
# subtracted, event profiling) of Phi, std_normal_lcdf_impl and
# std_normal_lcdf_derivative on develop e4745c87ff, a28c1c1637 (as opened)
# and 9b38ecbcb0 (now).  Binaries alternate, 2 passes x 3 rounds.
# Do not use `set -u`: /etc/profile reads unset variables.

cd "$SLURM_SUBMIT_DIR" || exit 1
LIB=/scratch/work/ave/smtest_lib_d12f488749
WORK=${TMPDIR:-/tmp}/ocl_dev_bench_$SLURM_JOB_ID

hostname
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
LIBCL=/usr/lib64/libOpenCL.so.1
[ -e "$LIBCL" ] || { echo "FATAL: no $LIBCL"; exit 1; }

mkdir -p "$WORK" && tar -xzf progression_trees.tar.gz -C "$WORK" || { echo "FATAL: untar"; exit 1; }
md5sum "$WORK"/t_*/stan/math/opencl/kernels/device_functions/Phi.hpp

build() {  # tree, defines
  g++ -std=c++17 -O2 -DCL_TARGET_OPENCL_VERSION=120 $2 -I "$WORK/t_$1" \
    -isystem "$LIB/opencl_v2025.07.22" -o "bench_$1" bench_ocl_device.cpp "$LIBCL" \
    || { echo "FATAL build $1"; exit 1; }
}
build develop ""
build a28c1c1637 "-DHAVE_LCDF"
build 9b38ecbcb0 "-DHAVE_ERFCX -DHAVE_LCDF"

for pass in 1 2; do
  for t in develop a28c1c1637 9b38ecbcb0; do
    if [ "$pass" = 1 ]; then
      "./bench_$t" "$t" 3 acc || echo "RUN FAILED $t"
    else
      "./bench_$t" "$t" 3 || echo "RUN FAILED $t"
    fi
  done
done
rm -rf "$WORK"
echo "EXIT=0"
