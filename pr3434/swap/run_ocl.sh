#!/bin/bash
#SBATCH --job-name=bb_ocl
#SBATCH --partition=gpu-debug
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=bb_ocl_%j.log
# Build and run the OpenCL tests of the stan-dev/math#3153 port on a GPU
# node, in <tree> (OpenCL configuration).
#
# STAGE=new      run the suites on the tree as deployed (the port)
# STAGE=develop  install the develop versions of the 7 changed OpenCL files
#                (dev_*.hpp in the job directory), run the suites, and
#                restore the port; the restore is checked by md5
# SUITES         colon-separated test paths below test/unit/math/opencl,
#                without _test.cpp; the default depends on STAGE
#
# Do not use `set -u`: /etc/profile reads unset variables.
cd "$SLURM_SUBMIT_DIR" || exit 1
JOBDIR=$PWD
SM="${SM:?set SM to the Stan Math tree with OpenCL}"
STAGE=${STAGE:-new}
ALL="device_functions/digamma_diff:kernel_generator/elt_function_cl"
ALL="$ALL:rev/beta_binomial_lpmf:rev/neg_binomial_lpmf:rev/neg_binomial_2_lpmf"
ALL="$ALL:rev/neg_binomial_2_log_lpmf:rev/neg_binomial_2_log_glm_lpmf"
ALL="$ALL:rev/student_t_lpdf:rev/lbeta"
REV="rev/beta_binomial_lpmf:rev/neg_binomial_lpmf:rev/neg_binomial_2_lpmf"
REV="$REV:rev/neg_binomial_2_log_lpmf:rev/neg_binomial_2_log_glm_lpmf"
REV="$REV:rev/student_t_lpdf"
# the beta_binomial test without its scalar-shape block: develop does not
# compile with scalar alpha and beta, so this variant shows the vector rows
# at run time; it lives only in the GPU tree, copied from the job directory
VEC=rev/beta_binomial_lpmf_vec
ALL="$ALL:$VEC"
REV="$REV:$VEC"
if [ -z "$SUITES" ]; then
  if [ "$STAGE" = "develop" ]; then SUITES=$REV; else SUITES=$ALL; fi
fi

# changed files: <path below stan/math>:<develop copy in the job directory>
CHANGED="opencl/prim/beta_binomial_lpmf.hpp:dev_prim_beta_binomial_lpmf.hpp
opencl/prim/neg_binomial_lpmf.hpp:dev_prim_neg_binomial_lpmf.hpp
opencl/prim/neg_binomial_2_lpmf.hpp:dev_prim_neg_binomial_2_lpmf.hpp
opencl/prim/neg_binomial_2_log_lpmf.hpp:dev_prim_neg_binomial_2_log_lpmf.hpp
opencl/prim/neg_binomial_2_log_glm_lpmf.hpp:dev_prim_neg_binomial_2_log_glm_lpmf.hpp
opencl/kernels/neg_binomial_2_log_glm_lpmf.hpp:dev_kernels_neg_binomial_2_log_glm_lpmf.hpp
opencl/kernel_generator/elt_function_cl.hpp:dev_kernel_generator_elt_function_cl.hpp
opencl/prim/student_t_lpdf.hpp:dev_prim_student_t_lpdf.hpp"

echo "=== node, stage $STAGE ==="
hostname
nvidia-smi --query-gpu=name --format=csv,noheader 2>&1
source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10
command -v g++ >/dev/null || { echo "FATAL: no g++"; exit 1; }
cd "$SM" || { echo "FATAL: no $SM"; exit 1; }

case ":$SUITES:" in
  *":$VEC:"*)
    [ -f "$JOBDIR/beta_binomial_lpmf_vec_test.cpp" ] \
      || { echo "FATAL: no $JOBDIR/beta_binomial_lpmf_vec_test.cpp"; exit 1; }
    cp "$JOBDIR/beta_binomial_lpmf_vec_test.cpp" test/unit/math/opencl/rev/ ;;
esac

# every test source must exist, or every suite reports BUILD FAILED
for s in ${SUITES//:/ }; do
  [ -f "test/unit/math/opencl/${s}_test.cpp" ] \
    || { echo "FATAL: no test/unit/math/opencl/${s}_test.cpp"; exit 1; }
done

show_md5() {
  echo "--- md5 of the changed and new OpenCL files ---"
  for line in $CHANGED; do md5sum "stan/math/${line%%:*}"; done
  md5sum stan/math/opencl/kernels/device_functions/digamma_diff.hpp \
    stan/math/opencl/kernels/device_functions/log_beta_ratio.hpp
}

BACKUP=""
restore() {
  if [ -n "$BACKUP" ]; then
    for line in $CHANGED; do
      cp "$BACKUP/$(basename "${line%%:*}")_$(echo "${line%%:*}" | md5sum | cut -c1-8)" \
        "stan/math/${line%%:*}"
    done
    rm -rf "$BACKUP"
    BACKUP=""
    echo "=== restored the port ==="
    show_md5
  fi
}
trap restore EXIT

if [ "$STAGE" = "develop" ]; then
  BACKUP=$(mktemp -d)
  for line in $CHANGED; do
    p=${line%%:*}
    d=${line##*:}
    [ -f "$JOBDIR/$d" ] || { echo "FATAL: no $JOBDIR/$d"; exit 1; }
    cp "stan/math/$p" "$BACKUP/$(basename "$p")_$(echo "$p" | md5sum | cut -c1-8)"
  done
  for line in $CHANGED; do
    cp "$JOBDIR/${line##*:}" "stan/math/${line%%:*}"
  done
  echo "=== installed the develop versions ==="
fi
show_md5

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

# Remove the binaries first: if a build fails, an old binary from another
# stage must not run. make -k builds every suite it can.
TARGETS=""
for s in ${SUITES//:/ }; do
  b="test/unit/math/opencl/${s}_test"
  rm -f "$b" "$b.o" "$b.d"
  TARGETS="$TARGETS $b"
done
echo; echo "=== build $(date +%T) ==="
make -k -j8 $TARGETS > "$JOBDIR/build_${STAGE}_${SLURM_JOB_ID}.txt" 2>&1
echo "build exit $? $(date +%T)"
grep -E "error" "$JOBDIR/build_${STAGE}_${SLURM_JOB_ID}.txt" | head -20

echo; echo "=== run ==="
for s in ${SUITES//:/ }; do
  exe="test/unit/math/opencl/${s}_test"
  out="$JOBDIR/out_${STAGE}_$(echo "$s" | tr / _)_${SLURM_JOB_ID}.txt"
  echo "##### $s"
  if [ ! -x "$exe" ]; then
    echo "SUMMARY $s BUILD_FAILED"
    continue
  fi
  "$exe" --gtest_color=no > "$out" 2>&1
  rc=$?
  npass=$(grep -c '^\[       OK \]' "$out")
  nfail=$(grep -c '^\[  FAILED  \] .*\.' "$out")
  nfail=$((nfail / 2))
  nexp=$(grep -c 'Failure$' "$out")
  echo "SUMMARY $s exit=$rc passed=$npass failed=$nfail failed_expects=$nexp"
  grep '^\[  FAILED  \] .*(' "$out" | head -10
  grep -E '^\[  (PASSED|FAILED)  \]|max relative' "$out" | head -20
  grep -E 'Failure|Expected|Which is|The difference|vector|scalar' "$out" \
    | head -24
done

restore
trap - EXIT
rm -f make/local test/unit/math/opencl/rev/beta_binomial_lpmf_vec_test.cpp
echo; echo "EXIT=0 $(date +%T)"
