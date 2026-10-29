#!/bin/bash
#SBATCH --job-name=nb_ocl
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=nb_ocl_%j.log
# Build and run the OpenCL tests of the negative binomial family after the
# fixes for jachymb's review comment (2026-10-12) on a GPU node, in
# <tree> (OpenCL configuration). Copy of
# run_ocl.sh with the file list and the stages changed.
#
# Stages, in one job:
#   new   the suites on the tree as deployed (the fixes)
#   base  the 468ff26dbb versions of the 5 changed OpenCL files
#         (base_*.hpp in the job directory) installed; the new reference
#         rows must fail; the fixes are restored afterwards and the restore
#         is checked by md5
# SUITES  colon-separated test paths below test/unit/math/opencl, without
#         _test.cpp
#
# Do not use `set -u`: /etc/profile reads unset variables.
cd "$SLURM_SUBMIT_DIR" || exit 1
JOBDIR=$PWD
SM="${SM:?set SM to the Stan Math tree with OpenCL}"
SUITES=${SUITES:-rev/neg_binomial_2_log_glm_lpmf:rev/neg_binomial_2_lpmf:rev/neg_binomial_2_log_lpmf:rev/neg_binomial_lpmf}

# changed files: <path below stan/math>:<468ff26dbb copy in the job directory>
CHANGED="opencl/kernels/neg_binomial_2_log_glm_lpmf.hpp:base_kernels_neg_binomial_2_log_glm_lpmf.hpp
opencl/prim/neg_binomial_2_log_glm_lpmf.hpp:base_prim_neg_binomial_2_log_glm_lpmf.hpp
opencl/prim/neg_binomial_2_lpmf.hpp:base_prim_neg_binomial_2_lpmf.hpp
opencl/prim/neg_binomial_2_log_lpmf.hpp:base_prim_neg_binomial_2_log_lpmf.hpp
opencl/prim/neg_binomial_lpmf.hpp:base_prim_neg_binomial_lpmf.hpp"

echo "=== node ==="
hostname
nvidia-smi --query-gpu=name --format=csv,noheader 2>&1
source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ >/dev/null || { echo "FATAL: no g++"; exit 1; }
cd "$SM" || { echo "FATAL: no $SM"; exit 1; }

for s in ${SUITES//:/ }; do
  [ -f "test/unit/math/opencl/${s}_test.cpp" ] \
    || { echo "FATAL: no test/unit/math/opencl/${s}_test.cpp"; exit 1; }
done
for line in $CHANGED; do
  [ -f "$JOBDIR/${line##*:}" ] || { echo "FATAL: no $JOBDIR/${line##*:}"; exit 1; }
done

show_md5() {
  echo "--- md5 of the changed OpenCL files ---"
  for line in $CHANGED; do md5sum "stan/math/${line%%:*}"; done
}

BACKUP=""
restore() {
  if [ -n "$BACKUP" ]; then
    for line in $CHANGED; do
      p=${line%%:*}
      cp "$BACKUP/$(echo "$p" | md5sum | cut -c1-8)_$(basename "$p")" \
        "stan/math/$p"
    done
    rm -rf "$BACKUP"
    BACKUP=""
    echo "=== restored the fixes ==="
    show_md5
  fi
}
trap restore EXIT

LIBCL=""
for d in /usr/lib64 /usr/lib/x86_64-linux-gnu; do
  for f in "$d/libOpenCL.so" "$d/libOpenCL.so.1"; do
    [ -e "$f" ] && { LIBCL="$f"; break 2; }
  done
done
cat > make/local <<EOC
STAN_OPENCL=true
OPENCL_DEVICE_ID=0
OPENCL_PLATFORM_ID=0
LDLIBS_OPENCL=$LIBCL
EOC

run_stage() {
  stage=$1
  # Remove the binaries first: if a build fails, an old binary from another
  # stage must not run. make -k builds every suite it can.
  TARGETS=""
  for s in ${SUITES//:/ }; do
    b="test/unit/math/opencl/${s}_test"
    rm -f "$b" "$b.o" "$b.d"
    TARGETS="$TARGETS $b"
  done
  echo; echo "=== stage $stage: build $(date +%T) ==="
  show_md5
  make -k -j8 $TARGETS > "$JOBDIR/build_${stage}_${SLURM_JOB_ID}.txt" 2>&1
  echo "build exit $? $(date +%T)"
  grep -E "error" "$JOBDIR/build_${stage}_${SLURM_JOB_ID}.txt" | head -20
  echo; echo "=== stage $stage: run ==="
  for s in ${SUITES//:/ }; do
    exe="test/unit/math/opencl/${s}_test"
    out="$JOBDIR/out_${stage}_$(echo "$s" | tr / _)_${SLURM_JOB_ID}.txt"
    echo "##### $s"
    if [ ! -x "$exe" ]; then
      echo "SUMMARY $stage $s BUILD_FAILED"
      continue
    fi
    "$exe" --gtest_color=no > "$out" 2>&1
    rc=$?
    npass=$(grep -c '^\[       OK \]' "$out")
    nfail=$(grep -c '^\[  FAILED  \] .*\.' "$out")
    nfail=$((nfail / 2))
    nexp=$(grep -c 'Failure$' "$out")
    echo "SUMMARY $stage $s exit=$rc passed=$npass failed=$nfail failed_expects=$nexp"
    grep '^\[  FAILED  \] .*(' "$out" | head -10
    grep -E 'Failure|Which is|The difference|vector|scalar' "$out" | head -16
  done
}

run_stage new

BACKUP=$(mktemp -d)
for line in $CHANGED; do
  p=${line%%:*}
  cp "stan/math/$p" "$BACKUP/$(echo "$p" | md5sum | cut -c1-8)_$(basename "$p")"
done
for line in $CHANGED; do
  cp "$JOBDIR/${line##*:}" "stan/math/${line%%:*}"
done
echo; echo "=== installed the 468ff26dbb versions ==="
run_stage base

restore
trap - EXIT
rm -f make/local
echo "=== done $(date +%T) ==="
