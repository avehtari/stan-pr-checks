#!/bin/bash
#SBATCH --job-name=nlcdfocl
#SBATCH --partition=gpu-debug
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=normal_lcdf_ocl_%j.log

# GPU speed of the normal_lcdf value and gradient: origin/develop against the
# erfcx identities. Single job, single GPU, minutes of walltime.
#
# Do not use `set -u`: /etc/profile reads unset variables and the job would
# fail before it starts.

cd "$SLURM_SUBMIT_DIR" || exit 1

SMDIR=${SMDIR:-$HOME/smnormal}
CLINC=${CLINC:-$HOME/erfcx_ocl/include}

echo "=== node ==="
hostname
nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader 2>&1

echo
echo "=== toolchain ==="
source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ || { echo "FATAL: no g++"; exit 1; }
g++ --version | head -1

echo
echo "=== the erfcx device function under test ==="
DEVFN=$SMDIR/stan/math/opencl/kernels/device_functions/erfcx.hpp
[ -f "$DEVFN" ] || { echo "FATAL: missing $DEVFN"; exit 1; }
md5sum "$DEVFN"
grep -c "fma" "$DEVFN"

echo
echo "=== locate libOpenCL ==="
# Link the library file by absolute path. The GPU nodes ship only the runtime
# soname libOpenCL.so.1, with no unversioned development symlink, so -lOpenCL
# fails at link time although the driver is present.
LIBCL=""
for d in /usr/lib64 /usr/lib/x86_64-linux-gnu "${CUDA_HOME:-/nonexistent}/lib64"; do
  for f in "$d/libOpenCL.so" "$d/libOpenCL.so.1"; do
    if [ -e "$f" ]; then
      echo "found $f"
      LIBCL="$f"
      break 2
    fi
  done
done
[ -n "$LIBCL" ] || { echo "FATAL: no libOpenCL found"; exit 1; }

[ -d "$CLINC/CL" ] || { echo "FATAL: no CL headers at $CLINC"; exit 1; }

echo
echo "=== build ==="
g++ -std=c++17 -O2 -DSTAN_OPENCL -isystem "$CLINC" -I "$SMDIR" \
    -o bench_normal_lcdf_ocl bench_normal_lcdf_ocl.cpp "$LIBCL" \
  || { echo "FATAL: build failed"; exit 1; }
echo "build ok"

echo
echo "=== run ==="
./bench_normal_lcdf_ocl "${ITERS:-2000}"
rc=$?
echo
echo "EXIT=$rc"
exit $rc
