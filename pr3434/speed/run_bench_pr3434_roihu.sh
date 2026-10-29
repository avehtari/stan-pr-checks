#!/bin/bash
#SBATCH --account=project_2020226
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:40:00
#SBATCH --job-name=pr3434_speed
#SBATCH --output=pr3434_speed_%j.log
# medium nodes rc4149, rc4150, rc4151 fail every job in the task prolog
# (mixtureis/ROIHU.md, "Traps").
#SBATCH --exclude=rc4149,rc4150,rc4151
#
# Speed of stan-dev/math#3434 against its merge base, all rows of the PR
# text speed section. Two stan/ trees from pr3434_trees.tar.gz (made with
# git archive):
#   dev  develop 510d263078 (merge base)
#   pr   #3434 64407093ba (merge of 510d263078 into the branch)
# Binaries per tree: bench_glm.cpp (the GLM), bench_dd.cpp (the gradients
# that use digamma_diff, lbeta, NB2 with 1000 counts), bench_st.cpp
# (student_t_lpdf, own binary), ir_st.cpp (instructions per student_t_lpdf
# call). One whole medium node (node-exclusive). Builds run in parallel;
# then ROUNDS interleaved rounds (dev pr), one binary at a time; then the
# instruction counts (valgrind cachegrind, and perf stat if it works).
# Header libraries from the clone lib/ at 64407093ba (pr3434_libs.tar.gz);
# the built TBB 2020.3 from CmdStan 2.40.0. Everything is unpacked in the
# node-local $TMPDIR. Flags -O2 -DNDEBUG.
#
# Deployed with this job: bench_glm.cpp, bench_dd.cpp, bench_st.cpp,
# ir_st.cpp, pr3434_trees.tar.gz, pr3434_trees.md5, pr3434_libs.tar.gz.
# Do not use `set -u` or `set -e`: /etc/profile reads unset variables.
set -o pipefail

export CSC_ENV_INIT_NON_INTERACTIVE=yes   # ROIHU.md gotcha R1
source /etc/profile
module purge
unset LD_PRELOAD                          # ROIHU.md gotcha R2

cd "$SLURM_SUBMIT_DIR" || exit 1
TBB="${TBB:-/projappl/project_2020226/cmdstan-2.40.0/stan/lib/stan_math/lib/tbb}"
TREES="dev pr"
ROUNDS="${ROUNDS:-5}"
REPS_GLM="${REPS_GLM:-4000}"
REPS_DD="${REPS_DD:-200}"
REPS_ST="${REPS_ST:-200}"
CALLS_IR="${CALLS_IR:-20000}"
SRC="bench_glm.cpp bench_dd.cpp bench_st.cpp ir_st.cpp"
DATA="pr3434_trees.tar.gz pr3434_trees.md5 pr3434_libs.tar.gz"

echo "== host: $(hostname)   date: $(date -Is)   g++: $(g++ --version | head -1)"
echo "== cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-)   cores: $(nproc)"
for f in $SRC $DATA; do
  [ -f "$f" ] || { echo "SETUP FAILED: no $f"; exit 1; }
done
md5sum $SRC $DATA
[ -f "$TBB/libtbb.so" ] || { echo "SETUP FAILED: no built TBB in $TBB"; exit 1; }
[ -n "$TMPDIR" ] && [ -d "$TMPDIR" ] || { echo "SETUP FAILED: no TMPDIR"; exit 1; }

WORK="$TMPDIR/pr3434_speed_$SLURM_JOB_ID"
LIB="$WORK/lib"
mkdir -p "$WORK" "$LIB" || exit 1
tar -xzf pr3434_trees.tar.gz -C "$WORK" || { echo "SETUP FAILED: untar trees"; exit 1; }
tar -xzf pr3434_libs.tar.gz -C "$LIB" || { echo "SETUP FAILED: untar libs"; exit 1; }
# every file of both trees against the md5 list made from the local archives
n_md5=$(wc -l < pr3434_trees.md5)
if (cd "$WORK" && md5sum -c --quiet "$SLURM_SUBMIT_DIR/pr3434_trees.md5"); then
  echo "== tree md5: all $n_md5 files equal to the local archives"
else
  echo "SETUP FAILED: tree md5 mismatch"
  exit 1
fi
echo "== md5 of the changed headers:"
for h in prim/prob/neg_binomial_2_log_glm_lpmf.hpp prim/prob/neg_binomial_2_lpmf.hpp \
         prim/prob/student_t_lpdf.hpp rev/fun/lbeta.hpp; do
  md5sum "$WORK"/t_*/stan/math/$h
done

FLAGS="-std=c++17 -O2 -DNDEBUG -pthread -D_REENTRANT"
INCS="-I $LIB/eigen_5.0.1 -I $LIB/boost_1.87.0"
INCS="$INCS -I $LIB/sundials_6.1.1/include -I $LIB/tbb_2020.3/include"
LINK="-L$TBB -ltbb -Wl,-rpath,$TBB"
BINS="glm dd st ir"
src_of() {
  case "$1" in
    ir) echo ir_st.cpp ;;
    *) echo "bench_$1.cpp" ;;
  esac
}
# develop has no digamma_diff; bench_dd.cpp times it only for the PR
for t in $TREES; do
  DD=""
  [ "$t" != dev ] && DD="-DHAS_DIGAMMA_DIFF"
  for b in $BINS; do
    g++ $FLAGS $DD -I "$WORK/t_$t" $INCS "$(src_of $b)" -o "$WORK/${b}_$t" \
      $LINK > "$WORK/build_${b}_$t.log" 2>&1 &
  done
done
wait
ok=1
for t in $TREES; do
  for b in $BINS; do
    if [ ! -x "$WORK/${b}_$t" ]; then
      echo "BUILD FAILED ${b}_$t"
      head -20 "$WORK/build_${b}_$t.log"
      ok=0
    fi
  done
done
[ "$ok" -eq 1 ] || { echo "SETUP FAILED: not all binaries built"; exit 1; }
echo "== built: $TREES x $BINS   $(date -Is)"

for r in $(seq 1 "$ROUNDS"); do
  for t in $TREES; do
    "$WORK/glm_$t" "$REPS_GLM" | sed "s/^/round$r $t /"
    "$WORK/dd_$t" "$REPS_DD" | grep '^TIME' | sed "s/^/round$r $t /"
    "$WORK/st_$t" "$REPS_ST" | grep '^TIME' | sed "s/^/round$r $t /"
  done
  echo "== round $r done   $(date -Is)"
done
echo "== load: $(cat /proc/loadavg)"

# ── instructions per student_t_lpdf call: (I(CALLS) - I(0)) / CALLS ──
echo "== valgrind: $(command -v valgrind) $(valgrind --version 2>&1)"
irefs() {
  valgrind --tool=cachegrind --cache-sim=no --cachegrind-out-file=/dev/null \
    "$@" 2>&1 | awk '/I *refs:/ {gsub(",", "", $NF); print $NF}'
}
perf_ins() {
  perf stat -x, -e instructions:u "$@" 2>&1 >/dev/null \
    | awk -F, '/instructions/ {print $1}'
}
echo "== perf: $(command -v perf)   perf_event_paranoid: $(cat /proc/sys/kernel/perf_event_paranoid)"
for range in A B C; do
  for kernel in value grad; do
    for t in $TREES; do
      i0=$(irefs "$WORK/ir_$t" "$range" "$kernel" 0)
      i1=$(irefs "$WORK/ir_$t" "$range" "$kernel" "$CALLS_IR")
      per=$(awk -v a="$i0" -v b="$i1" -v n="$CALLS_IR" \
        'BEGIN {if (a == "" || b == "") print "NA"; else printf "%.1f", (b - a) / n}')
      echo "IR cachegrind $range $kernel $t $per"
      p0=$(perf_ins "$WORK/ir_$t" "$range" "$kernel" 0)
      p1=$(perf_ins "$WORK/ir_$t" "$range" "$kernel" "$CALLS_IR")
      per=$(awk -v a="$p0" -v b="$p1" -v n="$CALLS_IR" \
        'BEGIN {if (a !~ /^[0-9]+$/ || b !~ /^[0-9]+$/) print "NA"; else printf "%.1f", (b - a) / n}')
      echo "IR perf $range $kernel $t $per"
    done
  done
done
rm -rf "$WORK"
echo "EXIT=0"
