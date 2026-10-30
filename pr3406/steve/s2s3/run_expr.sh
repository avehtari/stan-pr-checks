#!/bin/bash
#SBATCH --job-name=pr3406_expr
#SBATCH --partition=batch-csl,batch-skl
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --array=0-11
#SBATCH --output=pr3406_expr_%A_%a.out
#
# PR #3406, claims S2 and S3 (commit messages 201dcab952, f37ea7cb67).
# Each array task: unpack one commit, generate the expression tests for the
# task's functions with the nightly stanc3, compile with g++ -O0 under
# /usr/bin/time -v (S2: peak RSS), and for S3 link and run the tests.
# Do not use `set -u`: /etc/profile reads unset variables.
source /etc/profile
module load triton/2025.1-gcc scicomp-r-env fftw/3.3.10

LIB=/scratch/work/ave/smtest_lib_d12f488749
EMG="exp_mod_normal_lcdf exp_mod_normal_lccdf exp_mod_normal_cdf"
S2_FUNCS=($EMG normal_lcdf lognormal_lcdf)
S3_FUNCS="std_normal_lcdf std_normal_lccdf std_normal_cdf"

i=${SLURM_ARRAY_TASK_ID}
if [ "$i" -lt 5 ]; then
  CLAIM=S2; COMMIT=26256d4cdd; FUNCS=${S2_FUNCS[$i]}
elif [ "$i" -lt 10 ]; then
  CLAIM=S2; COMMIT=201dcab952; FUNCS=${S2_FUNCS[$((i - 5))]}
elif [ "$i" -eq 10 ]; then
  CLAIM=S3; COMMIT=201dcab952; FUNCS=$S3_FUNCS
else
  CLAIM=S3; COMMIT=f37ea7cb67; FUNCS=$S3_FUNCS
fi
TAG="${CLAIM}_${COMMIT}_$(echo $FUNCS | tr ' ' '+')"
echo "== task $i claim $CLAIM commit $COMMIT funcs: $FUNCS"
echo "== host: $(hostname)   date: $(date -Is)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
echo "== mem: $(free -g | awk '/Mem:/{print $2}') GB   job mem limit 64G"
echo "== g++: $(command -v g++)  $(g++ --version | head -1)"
echo "== python3: $(command -v python3) $(python3 --version 2>&1)"
cd "$SLURM_SUBMIT_DIR" || exit 1
SRC="src_$COMMIT.tar.gz"
[ -f "$SRC" ] || { echo "SETUP FAILED: missing $SRC"; exit 1; }
md5sum "$SRC"
[ -d "$LIB/sundials_6.1.1/include/arkode" ] || { echo "SETUP FAILED: lib"; exit 1; }
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; exit 1; }
W="$TMPDIR/expr_$i"
rm -rf "$W" && mkdir -p "$W" || exit 1
tar -xzf "$SRC" -C "$W" || exit 1
ln -s "$LIB" "$W/lib"
printf 'CXX=g++\nO=0\n' > "$W/make/local"
cd "$W" || exit 1

# ── stanc3 (nightly, as make/tests downloads it) and test generation ──
make -s test/expressions/stanc > stanc_dl.log 2>&1
[ -x test/expressions/stanc ] || { echo "SETUP FAILED: stanc"; cat stanc_dl.log; exit 1; }
echo "== stanc: $(test/expressions/stanc --version 2>&1 | head -1)"
python3 - $FUNCS <<'EOF'
import sys
sys.path[:0] = ["test", "test/expressions"]
import generate_expression_tests as g
g.main(sys.argv[1:], 1)
EOF
T=test/expressions/tests0_test
[ -f $T.cpp ] || { echo "GENERATION FAILED"; exit 1; }
echo "== generated $T.cpp: $(wc -l < $T.cpp) lines, $(grep -c '^TEST(' $T.cpp) TESTs" \
  " (Prim $(grep -c '^TEST(ExpressionTestPrim' $T.cpp)," \
  "Rev $(grep -c '^TEST(ExpressionTestRev' $T.cpp)," \
  "Fwd $(grep -c '^TEST(ExpressionTestFwd' $T.cpp))"
cp $T.cpp "$SLURM_SUBMIT_DIR/gen_${TAG}.cpp"
echo "== compile command:"
make -n $T.o 2>/dev/null | grep -- "-c " | tail -1

# ── compile the object under /usr/bin/time -v ──
/usr/bin/time -v make $T.o > compile.log 2> time.log
rc=$?
tail -40 compile.log
echo "== compile exit $rc"
grep -E "Maximum resident set size|Elapsed \(wall|User time|Exit status" time.log
rss_kb=$(awk -F: '/Maximum resident set size/{gsub(/ /,"",$2); print $2}' time.log)
echo "RESULT_RSS $TAG maxrss_kb=$rss_kb maxrss_GB=$(awk -v k="$rss_kb" 'BEGIN{printf "%.2f", k/1048576}') exit=$rc"
[ -f $T.o ] || { echo "RESULT $TAG BUILD FAILED"; exit 0; }

# ── S3: link and run ──
if [ "$CLAIM" = S3 ]; then
  make $T > link.log 2>&1 || { tail -20 link.log; echo "RESULT $TAG LINK FAILED"; exit 0; }
  ./$T --gtest_color=no > run_all.log 2>&1
  rc_all=$?
  echo "== full run exit $rc_all"
  grep -E "^\[  (FAILED|PASSED) |^\[==========\]|tests? ran" run_all.log
  echo "== failures (with messages):"
  awk '/^\[ RUN      \]/{buf=$0; next} {buf=buf "\n" $0} /^\[  FAILED  \] ExpressionTest/{print buf; buf=""}' \
    run_all.log | head -200
  ./$T --gtest_color=no --gtest_filter='ExpressionTestPrim*' > run_prim.log 2>&1
  rc_prim=$?
  echo "RESULT $TAG prim_exit=$rc_prim all_exit=$rc_all" \
    "prim_failed=$(grep -c '^\[  FAILED  \] ExpressionTestPrim' run_prim.log)"
  cp run_all.log "$SLURM_SUBMIT_DIR/run_${TAG}.log"
fi
echo "EXIT=0"
