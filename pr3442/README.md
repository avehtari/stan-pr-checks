# Reproduction scripts for stan-dev/math#3442

This directory holds the scripts behind the numbers in the text of
stan-dev/math#3442 ("Fix the underflow and gradient errors of the Wiener
functions").  Each section below names one claim of the PR text, the
commands that measure it, the expected output, and the trees that
produced it.

Trees: "develop" in the PR text is the merge base `43bc5684a3`; the PR
head is `ad0b7460d4`.  The scripts call them `base3442` and `head3442`.

Most numbers of the PR text come from intermediate states of the branch,
not from commits: "comb5" (before the score form for small sw and st0 and
before the underflow guard of the 7-parameter lcdf and lccdf) and "comb9"
(before the underflow guard).  The develop numbers come from a develop
tree with the same Wiener and `hcubature` files as `43bc5684a3`.  The
reruns of 2026-10-28 (§"Reruns") measure both commits.  Their outputs are
in `results/`; the expected outputs below are from these files.  Where a
rerun number differs from the earlier number, the section gives both.

## Layout

| path | content |
|---|---|
| `build.sh` | compiles one probe against one Stan Math tree |
| `refs/` | mpmath references: 4- and 5-parameter (p4, p5, test points), 7-parameter second route |
| `p45/` | the 4- and 5-parameter probes and their analysis |
| `w7/` | the 7-parameter probes, the C++ reference `probe_ref7.cpp`, the analyses |
| `speed/` | the speed benchmarks and the table script |
| `cluster/` | SLURM jobs (accuracy, speed) and `make_trees.sh` |
| `data/` | the references that the analyses read (§"Data files") |
| `results/` | the probe outputs and job logs of the reruns (§"Reruns") |

All commands run from this directory (`pr3442/`).  New outputs go to
`out/`, which git does not track.  To read the rerun outputs instead of
new ones, give `results` in place of `out` as the directory argument of
the analysis scripts.

## Requirements

- Python 3 with mpmath.  Tested with mpmath 1.3.0.
- g++ with C++17.  The cluster runs used the modules
  `triton/2025.1-gcc` and `scicomp-r-env` of Triton (conda-forge g++).
- A clone of stan-dev/math with the PR branch and a built `lib/tbb`:

  ```sh
  git clone git@github.com:stan-dev/math.git stan_math
  cd stan_math
  git fetch origin pull/3442/head:fix-wiener-accuracy-speed
  git checkout ad0b7460d4
  make -f make/standalone math-libs   # builds lib/tbb
  ```

Set up the two trees once:

```sh
export PR=$PWD/../stan_math                # the clone (stan/ and lib/)
export STAN_MATH_LIB=$PR/lib TBB_LIB=$PR/lib/tbb
mkdir -p out
STAN_MATH=$PR cluster/make_trees.sh out    # out/base3442.tar.gz, out/head3442.tar.gz
for t in base3442 head3442; do
  mkdir -p out/$t && tar xzf out/$t.tar.gz -C out/$t
done
# probe <source.cpp>: build against both trees, run, write out/<probe>_<tree>.out
probe() {
  n=$(basename "$1" .cpp)
  for t in base3442 head3442; do
    STAN_MATH=out/$t ./build.sh "$1" out/${n}_$t.bin \
      && out/${n}_$t.bin > out/${n}_$t.out
  done
}
```

`build.sh` uses `-O3`.  The cluster jobs use `-O2` for accuracy and
`-O3` for speed (`cluster/run_probe_tree.sh`, `cluster/run_bench_tree.sh`).
The 7-parameter lcdf and lccdf probe runs for tens of minutes per tree;
run it on a cluster:

```sh
# in a cluster directory that holds the probe, w7/swst0_points.hpp,
# cluster/run_probe_tree.sh, and the two tarballs one level up
sbatch --export=ALL,PROBE=probe_cdfacc.cpp,PRECS=1e-8,\
TREES=../base3442.tar.gz:../head3442.tar.gz,STAN_MATH_LIB=<lib>,MODULES=<m1:m2> \
  run_probe_tree.sh
```

## Reruns

Submitted 2026-10-28 on Triton (Aalto), both trees in each job.  Common
environment: `TREES=../base3442.tar.gz:../head3442.tar.gz`,
`STAN_MATH_LIB` = the `lib/` of a Stan Math clone with a built `lib/tbb`,
`MODULES=triton/2025.1-gcc:scicomp-r-env:fftw/3.3.10`.

| tag | job | node | script | other environment |
|---|---|---|---|---|
| `pr3442_acc_w4` | 21272912 | pht71 | `run_probe_tree.sh` | `PROBE=probe_w4.cpp` |
| `pr3442_acc_shift` | 21272915 | pht78 | `run_probe_tree.sh` | `PROBE=probe_shift.cpp` |
| `pr3442_acc_cdfacc` | 21272919 | pht82 | `run_probe_tree.sh` | `PROBE=probe_cdfacc.cpp`, `PRECS=1e-8` |
| `pr3442_acc_dacc` | 21272929 | pht82 | `run_probe_tree.sh` | `PROBE=probe_dacc.cpp` |
| `pr3442_acc_swst0` | 21272930 | skl9 | `run_probe_tree.sh` | `PROBE=probe_swst0.cpp` (default `PRECS` 1e-4/1e-8) |
| `pr3442_speed` | 21272931 | skl19 | `run_bench_tree.sh` | `PROBES=probe_bench.cpp:probe_bench4.cpp`, `REPS=3` |

Tarballs: `cluster/make_trees.sh` gives md5 `a3afee38e0a025ab641c754b346921c5`
(`base3442.tar.gz`) and `587e038bdc836a0386a54cdc8196db56`
(`head3442.tar.gz`), the same as the tarballs of the jobs.

Status: all 6 jobs completed on 2026-10-29 (exit 0 for each tree, no
throws, empty build logs).  Compiler g++ 14.3.0 (conda-forge).  The md5
lines in the logs (`results/pr3442_acc_<job>.log`,
`results/pr3442_speed_21272931.log`) agree with the two tarballs, the
probe sources and `w7/swst0_points.hpp` of this directory.  The job
scripts on Triton have the md5 of `cluster/`.  Outputs:
`results/<probe>_<tree>.out`, `results/probe_bench{,4}_bench.out`.

## 4-parameter `wiener_lpdf` (not changed)

Claim: "its own series; not changed, because it is accurate (300 points
against mpmath: value error at most 3.6e-9, gradients 5.2e-8)".

```sh
probe p45/probe_w4.cpp
python3 p45/analyze_w4.py data/ref_w4_p4p5.out out base3442 head3442
```

Expected: row `OLD p4`, `v max` 3.6e-09 and `g max` 5.2e-08 in both trees
(the file is not changed).  The p4 set: 300 points, a 0.5-3, v -5..5,
w 0.2-0.8, y - t0 0.05-3.05, t0 = 0.1, sv = 0.  References:
`refs/ref_w4_mpmath.py` at 50 digits (`python3 refs/ref_w4_mpmath.py
out/probe_w4_head3442.out > out/ref_w4.out` writes them again; the
p5 lcdf and lccdf rows need quadrature and take long).  Measured on
develop and comb5; rerun at `43bc5684a3` and `ad0b7460d4` (job
21272912, pht71): 3.6e-09 and 5.2e-08 in both trees.

## 5-parameter `wiener_lpdf`

Claim: "Three faults on the path for small densities, with value errors
up to 1e-2 at ordinary points"; table row "5-parameter lpdf | 1.1e-2 ->
4.6e-12 | 1.8e-2 -> 2.1e-14".

Same commands as above; rows `LPDF p4` and `LPDF p5`.  Expected
(rerun, job 21272912, pht71; the earlier develop and comb5 numbers are
the same):

| row | base value | head value | base grad | head grad |
|---|---|---|---|---|
| `LPDF p4` (sv = 0) | 1.1e-02 | 4.6e-12 | 1.8e-02 | 1.8e-14 |
| `LPDF p5` (sv 0-1.5) | 3.7e-03 | 1.9e-11 | 1.2e-02 | 2.1e-14 |

The table row of the PR text takes the maximum over both sets for
develop and for the PR gradients, but the p4 value for the PR value: the
maximum over both sets is 1.9e-11, not 4.6e-12.

The faults themselves at single points (log density down to -592) and on
a 20 000-point scan:

```sh
probe p45/probe_w5grad.cpp
python3 refs/ref_w5grad_mpmath.py out/probe_w5grad_head3442.out
```

Expected on develop: value error up to 0.47 at the 10 points, gradients
2.8e-2 relative, scan line `value worst 0.058`; on the fixed code value
2.1e-8, gradients 3.3e-13, scan 2.3e-10.  Measured on the single-fix
branch of the 5-parameter faults (Triton jobs 20871396, 20871471).  Not
rerun.  The 5 points of the new test `rev/prob/wiener5_lpdf_test.cpp`
have their references from `refs/ref_w5grad_mpmath.py` (function
`ref`).

## 4-parameter lcdf: gradient signs, `logMill`

Claim: "The small-time series of the derivatives with respect to a and w
had sign errors: gradient errors up to 1.6e-4, and up to 3e-3 in the
7-parameter lcdf and lccdf."  Table row "4-parameter lcdf | 1.1e-8 ->
1.1e-8 | 1.6e-4 -> 9.5e-9".

`probe p45/probe_w4.cpp`, then `p45/analyze_w4.py` as above: row
`LCDF p4`.  Expected: value 1.1e-08 in both trees, gradients 1.6e-04
(base, at `w`) and 9.5e-09 (head; job 21272912, pht71; comb5 gave the
same numbers).
The 7-parameter part (2.9e-3 in the lcdf, D set, 2 dimensions; 1.3e-3 at
the p5 set) is in §"7-parameter lcdf and lccdf".  The `logMill` change
has no probe here; its check is the prim test
`prim/prob/wiener4_lcdf_unnorm_test.cpp`.

## 4-parameter lccdf: cancellation of P - F

Claim: "P - F loses all digits where F is close to P (-inf at 8 of 300
ordinary points)."  Table row "4-parameter lccdf | -inf -> 1.4e-14 |
-inf -> 2.9e-9".

`probe p45/probe_w4.cpp`, then `p45/analyze_w4.py`: row `LCCDF p4`.
Expected on the base: column `inf` 8 (points with -inf), `n>1e-6` 51
(the 8 and 43 more points wrong by more than 1e-6); on the head: value
1.4e-14, gradients 2.9e-09 (job 21272912, pht71; comb5 gave the same
numbers).

## #3329: 7-parameter lpdf below the smallest double

Claims: the log shift of the integrand; table row "7-parameter lpdf,
#3329 point and scan | -inf -> 4.8e-11 | d/dsw -1e8 -> 5.3e-10".

```sh
probe w7/probe_shift.cpp
python3 w7/compare_shift.py out base3442 head3442        # values
python3 w7/analyze_scan.py data/ref7_scan_score.out out base3442 head3442
```

The scan: 16 points y = 3e-6 .. 3e-5 (a = 0.1, t0 = 0, w = 1e-4,
v = -10, sv = 0.1, sw = 1e-8, st0 = 1e-3), precision 1e-8.  Value
references: the `R` rows of `probe_shift.cpp` (Boost Gauss-Kronrod over the
log-scaled 5-parameter density, sw -> 0).  d/dsw references:
`probe_ref7.cpp` with `REF_SCORE=1` (differentiation under the integral):

```sh
STAN_MATH=out/base3442 ./build.sh w7/probe_ref7.cpp out/probe_ref7.bin
REF_SCORE=1 REF_LIST=$(python3 -c 'print(":".join("%.17g/0.1/0/1e-4/-10/0.1/1e-8/1e-3"
  % (3e-6 * 10 ** (k / 15)) for k in range(16)))') out/probe_ref7.bin > out/ref_scan.out
```

Expected (rerun, job 21272915, pht78; earlier: develop and comb9, job
20878532, the same numbers): base -inf at s0-s5 (6 of 16) and
d/dsw = -1e8 (= -1/sw) at s6-s8; head value error at most 4.8e-11 (s0,
`compare_shift.py`, against the `R` rows), d/dy 3.9e-12, d/dsw at most
5.3e-10 (s15, `analyze_scan.py`); d/dst0 = -1000 (exact) in both trees.
Against the score-form references of `analyze_scan.py`, the head value
error is at most 3.8e-11 (s11) and d/dy 4.3e-12.  The three `P` rows
(y = 1e-6) show about 4e-10 in `compare_shift.py` even if exact, because
the `R` reference leaves out sw = 1e-8 (sw^2 / (24 w^2) = 4e-10); the
4.8e-11 of the PR text is the maximum over the 16 scan points.

## d/dsw with sw > 0 and st0 = 0 (extreme grid)

Claim: "With sw > 0 and st0 = 0, d/dsw used an absolute series tolerance
(relative errors up to 2e3 where the density is small)."

```sh
probe w7/probe_shift.cpp      # (same run as above)
python3 w7/compare_ref7.py data/ref7_E_D.out out base3442 head3442
python3 w7/compare_shift.py out base3442 head3442   # section E: finite counts
```

The extreme grid: 128 points, y - t0 in {1e-6, 1e-4, 10, 50},
a in {0.1, 8}, v in {-30, -10, 10, 30}, w in {1e-4, 0.5}, st0 in {0, 1e-3},
sw = w / 2 or (1 - w) / 2, sv = 0.1, precision 1e-6.  References:
`probe_ref7.cpp` (default mode; `data/ref7_E_D.out`, Triton job 20868511).
Expected (rerun, job 21272915, pht78; earlier develop and comb9 gave the
same numbers): base `gerr` 2.0e+03 at `sw` for e84, e85, e88, e89; head
at most 2.1e-09 there.  Section E: base finite at 32 of 128 points (value
-inf at 96), head finite at 108 (value -inf at 20), no throws.  The 20
head -inf points all have a = 8 and y - t0 <= 1e-4 (12 at 1e-6, 8 at
1e-4).  Head errors over the 108 finite points: gradients at most
6.4e-08 (e35, d/dw), except e3, e7, e11, e15 (a = 0.1, y - t0 = 1e-6,
w = 0.5, st0 = 1e-3), where the value error is 6.4e-04 in both trees and
the head d/dsw error 6.4e-04 (base d/dsw = -1/sw = -4).

## d/dsw integration limits (st0 > y - t0)

Claim: "The d/dsw integral over the st0 variable used the limits of the
sw variable: for (y - t0)/st0 < 0.004, d/dsw was -1/sw."

No probe in this directory.  The check is the new test in
`rev/prob/wiener_full_test.cpp` (two points with (y - t0)/st0 = 0.0033;
develop -5 and -10, references 47.08 and 72.97).  The references of the
d/dsw test points are in `data/ref7_test_dsw.out` (`probe_ref7.cpp` with
`REF_SCORE=1` and `REF_LIST` = the points in columns 3-10 of that file).

## d/dsw and d/dst0 at small sw and st0

Claims: "in develop the d/dsw error at sw = 1e-8 is up to 1.1"; "Above
1e-3 the end-point forms stay: their error is about 3e-10/sw at the
default precision"; table row "7-parameter lpdf, sw or st0 = 1e-8 |
1.1 -> 3.3e-9".

```sh
probe w7/probe_swst0.cpp      # precisions 1e-4 and 1e-8
python3 w7/analyze_swst0.py data/ref7_swst0_s0.out,data/ref7_swst0_s1.out,\
data/ref7_swst0_s2.out,data/ref7_swst0_s3.out out base3442 head3442
```

The set (`w7/swst0_points.hpp`): 12 base points x 4 series x 5 sizes =
240 points; sw or st0 in {1e-8, 1e-6, 1e-4, 1e-2, 0.1} (st0 0.5 d at the
largest size).  References: `probe_ref7.cpp` with `REF_SWST0=1
REF_SCORE=1` and `REF_RANGE` 0:60, 60:120, 120:180, 180:240 (Triton jobs
20878282-20878285).  Expected at precision 1e-4, column `g var`
(base and head: rerun, job 21272930, skl9; develop and comb9 gave the
same numbers; comb5 from the earlier runs):

| series, size | base | comb5 | head |
|---|---|---|---|
| `wb 1e-8` (d/dsw, st0 > 0) | 1.1e+00 | 2.9e-02 | 3.3e-09 |
| `wb 1e-6` | 1.1e-02 | 2.9e-04 | 3.0e-11 |
| `wb 1e-4` | 1.1e-04 | 2.9e-06 | 2.1e-11 |
| `wb 1e-2` | 1.1e-06 | 2.9e-08 | 2.9e-08 |
| `wb big` (sw = 0.1) | 4.4e-07 | 4.2e-08 | 4.2e-08 |
| `tb 1e-8` (d/dst0, sw > 0) | 5.3e-04 | 5.3e-04 | 5.8e-11 |
| `tb 1e-6` | 5.3e-06 | 5.3e-06 | 5.8e-11 |
| `tb 1e-4` | 5.3e-08 | - | 9.8e-11 |

comb5 uses the end-point form at every size; its errors give the
2.9e-10/sw of the end-point form.  At the head the end-point form remains
for sw >= 1e-3: 2.9e-08 at sw = 1e-2 (2.9e-10/sw).  At sw = 0.1 the error
is 4.2e-08, which is not 1/sw scaling but the integration tolerance.

## Density tolerance of the 7-parameter lpdf

Claim: "The density integral uses a relative tolerance of at most 0.9e-6:
the log density error at the default precision goes from 1.7e-6 to
1.0e-7."  Table row "7-parameter lpdf, ordinary | 1.7e-6 -> 1.0e-7 |
2.3e-5 -> 2.1e-5".

```sh
probe w7/probe_dacc.cpp       # precision 1e-4; PRECS=1e-4/1e-6/1e-8 for more
python3 w7/analyze_tol.py data/ref7_D_lpdf.out out base3442 head3442
```

The D set: 300 points, a 0.5-3, v -3..3, w 0.2-0.8, t0 0.1-0.4, sv 0-1.5,
sw up to 0.8 of its range (0 at every third point), st0 0-0.3 (0 at
every third point), y - t0 0.05-2.05.  References: `probe_ref7.cpp
REF_D=1` (`data/ref7_D_lpdf.out`, error estimate at most 1e-11).
Expected at precision 1e-4 (rerun, job 21272929, pht82): base value
1.7e-06, gradients 2.3e-05 (at st0); head value 1.0e-07, gradients
2.4e-05 (at sw, point d192).  The PR text has 2.1e-05 for the head
gradients; that number came from a tree with the `hcubature` fix and
this tolerance only (Triton jobs 20873627, 20873633).  comb5 gave 2.4e-05
(FINDINGS of 2026-10-10, comparison with #3435).  Both are below the
requested precision 1e-4.

## `hcubature`

Claims: "the loop condition used the error and value of the first box
...; and a split box stayed in the list with error 0 and could be split
and counted again."  "makes the 7-parameter lpdf 11 to 14 times faster
when st0 > 0."

```sh
STAN_MATH=out/base3442 ./build.sh w7/probe_hcub.cpp out/probe_hcub.bin
out/probe_hcub.bin
```

The probe runs a copy of the develop algorithm with counters and the
fixed copy on two D points (d17, d221) at precisions 1e-4 to 1e-10.
Expected on develop (Triton job 20868088): d17 without shift at 1e-10
`stan_rel=6.000e+00` (the integral is 7 times the reference),
`dead_reselect=6`; d221 with shift at 1e-8 and 1e-10 `stan_rel=4.900e+01`,
`dead_reselect=49`; the develop copy runs 199 splits (the evaluation
limit) wherever the first box fails, the fixed copy 1 split, and
`fixed_rel` is at most 2.2e-16.  The head was not measured with this
probe; its `stan_rel` should be near `fixed_rel`.  The unit test
`prim/functor/hcubature_test.cpp` checks the step integral (develop
99.5 instead of 0.5).  "11 to 14 times faster" is the inverse of the
run-time ratios 0.09 and 0.07 (§"Run-time ratios").

## 7-parameter lcdf and lccdf

Claims: the integration tolerance at most 0.9e-6; "without them the log
ccdf was off by up to 1.4 where st0 > y - t0"; table rows "7-parameter
lcdf, 1-2 integration dims | 3.1e-7 -> 3.1e-7 | 2.9e-3 -> 3.9e-7" and
"7-parameter lccdf, 1-2 integration dims | -inf, 1.3 -> 8.4e-9 | 6.3e2 ->
2.7e-8".

```sh
PRECS=1e-8 probe w7/probe_cdfacc.cpp     # D set; cluster recommended
R=data/ref7_D_cdf_fd_r0.out,data/ref7_D_cdf_fd_r1.out,data/ref7_D_cdf_fd_r2.out,data/ref7_D_cdf_fd_r3.out
python3 w7/analyze_evlimit.py $R out base3442 head3442
python3 w7/analyze_cdf_fd.py d ${R//,/ } -- out base3442 head3442
probe p45/probe_w4.cpp                   # p5 set, sv only (1 dim)
python3 w7/analyze_cdf_fd.py p5 data/ref7_p5_cdf_fd.out -- out base3442 head3442
```

"1 dim" is the p5 set (sv > 0, sw = st0 = 0); "2 dims" and "3 dims" are
the D points with sv and one or both of sw, st0.  References:
`probe_ref7.cpp` with `REF_D=1 REF_CDF_FD=1` (values by Gauss-Kronrod,
split at the kink of the st0 weight; gradients by central differences,
noise at most 1.4e-9; Triton jobs 20877289, 20877290, 20877298, 20877299)
and with `REF_P5=1` for p5.  Expected `analyze_evlimit.py` output
(rerun, job 21272919, pht82, 4 cores of a shared node, `-O2`):

```
fn     tree         dim n    v max    v p99    g max    g p99    ms/pass
LCDF   base3442     2   200  3.1e-07  1.0e-08  2.9e-03  2.0e-03  14220.7
LCDF   base3442     3   100  4.3e-05  4.3e-05  1.3e-03  1.3e-03  9386.0
LCDF   head3442     2   200  3.1e-07  1.0e-08  3.9e-07  4.6e-08  7686.6
LCDF   head3442     3   100  4.3e-05  4.3e-05  8.4e-04  8.4e-04  6602.0
LCCDF  base3442     2   200  1.3e+00  1.1e+00  6.3e+02  4.8e+00  16941.6
LCCDF  base3442     3   100  1.4e+00  1.4e+00  7.8e+00  7.8e+00  10927.4
LCCDF  head3442     2   200  8.4e-09  2.5e-09  2.7e-08  2.4e-08  12340.0
LCCDF  head3442     3   100  1.6e-05  1.6e-05  4.2e-04  4.2e-04  9901.7
```

The errors are those of the earlier develop and comb5 runs.  Only the
`ms/pass` column differs (earlier: develop 12149.4, 8102.6, 14434.6,
9406.5; comb5 6457.3, 5554.6, 10375.9, 8299.1).  The base lccdf has a
value error > 1e-2 or -inf at 7 of the 200 2-dimension points and at 10
of the 100 3-dimension points.

`analyze_cdf_fd.py p5` (job 21272912): base lcdf gradients 1.3e-03,
lccdf 4 points far (-inf) and gradients 4.9e+02; head lcdf 9.5e-08 /
9.6e-08, lccdf 3.6e-11 / 2.1e-09 (as develop and comb5).  The "1.4" of
the late-start claim is the 3-dimension row; the 2-dimension rows give
1.3.

Small ccdf values (claims "finite down to a ccdf of about 1e-308",
"log ccdf -39, -120 and -299", "Where the value integral underflows to 0,
the value is -inf and the gradient integrals are skipped"):

```sh
PTS=3.1/2/0.1/0.5/5/0.01/0/0:9/2/0.1/0.5/5/0.01/0/0:22.1/2/0.1/0.5/5/0.01/0/0:\
60/2/0.1/0.5/5/0.01/0/0:1.6e-5/0.1/0/1e-4/-10/0.1/0/0:1e-6/0.1/0/1e-4/-10/0.1/0/0 \
PRECS=1e-8 probe w7/probe_cdfacc.cpp
```

Expected (Triton jobs 20880197, 20880468, 20880476; develop and the
branch state with the guard): lccdf at l0, l1, l2 develop -inf, PR
-39.05, -120.0, -299.4; l3 (log ccdf about -820) -inf in both; lcdf at
l4 -319.0241 in both, at l5 -inf (develop: NaN partials; PR: partials 0).
References for l1 and l2: `data/ref7_underflow_cdf_fd.out` (its ids l0,
l1 are the points y = 9 and y = 22.1).  Not rerun with this probe;
`w7/probe_pr3442.cpp` measures the same lccdf points and the lcdf at
y - t0 = 1e-6 at both commits (§"Candidate test rows at the merge base
and the head", P4 and P5).

## Accuracy table

| row of the PR table | section |
|---|---|
| 4-parameter lcdf | §"4-parameter lcdf: gradient signs, `logMill`" |
| 4-parameter lccdf | §"4-parameter lccdf: cancellation of P - F" |
| 5-parameter lpdf | §"5-parameter `wiener_lpdf`" |
| 7-parameter lpdf, ordinary | §"Density tolerance of the 7-parameter lpdf" |
| 7-parameter lpdf, #3329 point and scan | §"#3329: 7-parameter lpdf below the smallest double" |
| 7-parameter lpdf, sw or st0 = 1e-8 | §"d/dsw and d/dst0 at small sw and st0" |
| 7-parameter lcdf, lccdf, 1-2 integration dims | §"7-parameter lcdf and lccdf" |

Errors: largest over the points; values absolute on the log scale;
gradients relative to max(|ref|, 1).

The rows with the rerun numbers (`43bc5684a3` -> `ad0b7460d4`; jobs of
§"Reruns"):

| function | value | gradients |
|---|---|---|
| 4-parameter lcdf | 1.1e-8 -> 1.1e-8 | 1.6e-4 -> 9.5e-9 |
| 4-parameter lccdf | -inf -> 1.4e-14 | -inf -> 2.9e-9 |
| 5-parameter lpdf (p4; p5) | 1.1e-2 -> 4.6e-12; 3.7e-3 -> 1.9e-11 | 1.8e-2 -> 1.8e-14; 1.2e-2 -> 2.1e-14 |
| 7-parameter lpdf, ordinary | 1.7e-6 -> 1.0e-7 | 2.3e-5 -> 2.4e-5 |
| 7-parameter lpdf, #3329 point and scan | -inf -> 4.8e-11 | d/dsw -1e8 -> 5.3e-10 |
| 7-parameter lpdf, sw or st0 = 1e-8 | | 1.1 -> 3.3e-9 |
| 7-parameter lcdf, 1-2 integration dims | 3.1e-7 -> 3.1e-7 | 2.9e-3 -> 3.9e-7 |
| 7-parameter lccdf, 1-2 integration dims | -inf, 1.3 -> 8.4e-9 | 6.3e2 -> 2.7e-8 |

Differences from the PR text: the 7-parameter lpdf "ordinary" head
gradients (2.4e-5, PR 2.1e-5; §"Density tolerance of the 7-parameter
lpdf"), and the 5-parameter lpdf head value over both sets (1.9e-11, PR
4.6e-12, the p4 value).

## References and test points

Claim: "References: mpmath (50 digits) for the 4- and 5-parameter
functions; for the 7-parameter functions, my own adaptive Gauss-Kronrod
integration of the 5-parameter density (not `hcubature`), with gradients
under the integral or by finite differences."

- 4- and 5-parameter: `refs/ref_w4_mpmath.py` (p4, p5) and
  `refs/ref_w5grad_mpmath.py` (test points), 50 digits, gradients by
  `mpmath.diff`.  Each is one route: `ref_w4_mpmath.py` uses the
  small-time series for t < 1, `ref_w5grad_mpmath.py` the large-time
  series with 400 terms.
- 7-parameter: `w7/probe_ref7.cpp`.  It calls the internal 5-parameter
  functions of the tree it is built against with a relative series
  tolerance; the files in `data/` were built against develop.
- Second route for the 7-parameter test rows: `refs/ref_w7_mpmath.py`
  (mpmath quadrature, central differences), one point per call:
  `python3 refs/ref_w7_mpmath.py LCCDF id 9 2 0.1 0.5 5 0.01 0 0`.

Point sets: p4, p5 (§"4-parameter `wiener_lpdf`"), D (§"Density
tolerance"), the #3329 scan (§"#3329"), the 240-point set (§"d/dsw and
d/dst0 at small sw and st0"), the extreme grid.

## Run-time ratios

Claim: "Run time with gradients, this PR / develop: 4-parameter lcdf
0.89, lccdf 0.62; 5-parameter lpdf 1.09-1.13; 7-parameter lpdf 0.07 (sw
and st0), 0.09 (st0 only), 1.05 (sw only); 7-parameter lcdf 0.04-0.67,
lccdf 0.06-0.87."

```sh
# on a cluster, from a directory with the two benchmarks, the job script and
# the tarballs one level up
sbatch --export=ALL,PROBES=probe_bench.cpp:probe_bench4.cpp,REPS=3,\
TREES=../base3442.tar.gz:../head3442.tar.gz,STAN_MATH_LIB=<lib>,MODULES=<m1:m2> \
  run_bench_tree.sh
python3 speed/table_bench.py probe_bench4_bench.out    # 4- and 5-parameter
python3 speed/table_bench.py probe_bench_bench.out     # 7-parameter lpdf
```

The job builds both benchmarks against both trees with `-O3 -DNDEBUG
-D_REENTRANT -pthread`, then runs them on one exclusive node in 3
interleaved rounds (base, head, base, head, ...).  Each (set, mode) runs
whole passes over its points for at least 1 s after one warm-up pass.
Each timed call adds `sum * 1e-300` to y, so that the loop depends on
the earlier results; this changes no bit of y.  Mode `var` is the value
with the gradient, `dbl` the value only.  The rows of the claim: 4- and
5-parameter `lcdf_p4`, `lccdf_p4`, `lpdf_p4`, `lpdf_p5`; 7-parameter
`lcdf_p5`, `lccdf_p5` (sw = st0 = 0); 7-parameter lpdf `ord_both`,
`ord_st0`, `ord_sw`.  The 7-parameter lcdf and lccdf ratios at the D set
are the `ms/pass` columns of `analyze_evlimit.py` (§"7-parameter lcdf and
lccdf"), not a benchmark.

Rerun: job 21272931 (`pr3442_speed`), node skl19, Intel Xeon Gold 6148
at 2.40 GHz, governor `performance`, exclusive node, g++ 14.3.0
(conda-forge), `-O3 -DNDEBUG -D_REENTRANT -pthread`, 3 interleaved
rounds.  Head / base ratio of the medians; "rounds" is the range of the
per-round ratios; "earlier" is the ratio of the earlier runs (PR text,
FINDINGS of 2026-10-10):

| row | set | var | rounds | dbl | rounds | earlier var / dbl |
|---|---|---|---|---|---|---|
| 4-p lpdf (old) | `old_p4` | 1.001 | 1.001-1.003 | 1.001 | 1.000-1.003 | 1.00 / 1.00 |
| 4-p lcdf | `lcdf_p4` | 0.900 | 0.896-0.902 | 0.834 | 0.833-0.835 | 0.89 / 0.82 |
| 4-p lccdf | `lccdf_p4` | 0.646 | 0.645-0.652 | 1.336 | 1.331-1.337 | 0.62 / 1.30 |
| 5-p lpdf, sv = 0 | `lpdf_p4` | 1.123 | 1.123-1.123 | 1.177 | 1.176-1.183 | 1.13 / 1.17 |
| 5-p lpdf, sv > 0 | `lpdf_p5` | 1.090 | 1.088-1.091 | 1.137 | 1.136-1.140 | 1.09 / 1.13 |
| 7-p lpdf, sw, st0 | `ord_both` | 0.073 | 0.073-0.073 | 0.104 | 0.104-0.104 | 0.070 / 0.10 |
| 7-p lpdf, st0 | `ord_st0` | 0.090 | 0.090-0.091 | 0.082 | 0.081-0.084 | 0.090 / 0.082 |
| 7-p lpdf, sw | `ord_sw` | 1.058 | 1.056-1.062 | 1.084 | 1.081-1.107 | 1.05 / 1.08 |
| 7-p lpdf, extreme | `ext` | 1.006 | 1.006-1.010 | 0.225 | 0.225-0.226 | - |
| 7-p lcdf, p5 | `lcdf_p5` | 0.042 | 0.042-0.042 | 0.036 | 0.036-0.037 | 0.042 / 0.036 |
| 7-p lccdf, p5 | `lccdf_p5` | 0.057 | 0.056-0.057 | 0.052 | 0.052-0.052 | 0.056 / 0.050 |

The spread of each tree over the 3 rounds is at most 3.2 % (`ord_st0`
dbl, head).  7-parameter lcdf and lccdf at the D set (job 21272919, one
pass per tree, base then head, `-O2`, 4 cores of a shared node; with
gradients): lcdf 0.54 (2 dims) and 0.70 (3 dims), lccdf 0.73 and 0.91
(earlier 0.52, 0.67, 0.70, 0.87, job 20877547).  With these numbers the
claim reads: 4-parameter lcdf 0.90, lccdf 0.65; 5-parameter lpdf
1.09-1.12; 7-parameter lpdf 0.07, 0.09, 1.06; 7-parameter lcdf
0.04-0.70, lccdf 0.06-0.91.

Earlier runs (Triton, one core of a shared node, modules as above):
4- and 5-parameter job 20877927 (median of 5 rounds, spread <= 7.5 %,
comb5); 7-parameter lpdf job 20874806 (Xeon E5-2680 v3 at 2.50 GHz,
median of 3 rounds, spread <= 5 %, a branch state before comb5).  The
earlier benchmarks had no `sum * 1e-300` term; their calls take
microseconds, thousands of cycles, so the compiler did not remove them.

## Head against base, per point

`w7/regress_head_vs_base.py` reads the rerun outputs in `results/` and the
references in `data/`, and lists each (point, quantity) where the head
error is above 4 times the base error and above 1e-8.

```
python3 w7/regress_head_vs_base.py
```

Result: 233 pairs, none in the 4- and 5-parameter sets, no value error
above 1.4e-8.  227 pairs are 7-parameter lpdf gradients at precision 1e-4
(largest: d192 d/dsw 8.8e-9 -> 2.35e-5): the base is accurate there
because its `hcubature` runs to the evaluation limit; the head stops at
the requested tolerance.  The other 6 pairs are below 6.4e-8.

## Older 4-parameter `wiener_lpdf` at small densities

Claim: "The older 4-parameter `wiener_lpdf` returns -inf at short times
where the log density is below about -708 (35 of 600 grid points); where
it is finite, it is accurate."  This PR does not change that file.

`p45/probe_oldw4.cpp` evaluates `wiener_lpdf(y, a, t0, w, v)` and the
5-parameter lpdf with sv = 0 on a grid of 600 points (t in {1e-3, 0.01,
0.05, 0.2, 1, 3, 10, 30}, a in {0.5, 1, 2, 4, 8}, v in {-20, -5, 0, 5, 20},
w in {0.1, 0.5, 0.9}, t0 = 0.1).  `refs/ref_w4_mpmath.py` scores the
LPDF rows (50 digits).

```
mkdir -p out
STAN_MATH=<tree> OPT=-O2 ./build.sh p45/probe_oldw4.cpp out/probe_oldw4
out/probe_oldw4 > out/probe_oldw4.out
python3 refs/ref_w4_mpmath.py out/probe_oldw4.out > out/ref_oldw4.out
python3 p45/analyze_oldw4.py out/probe_oldw4.out out/ref_oldw4.out
```

Expected output (head ad0b7460d4; Triton job 21277321, `csl46`, g++ 14.3,
`-O2`; `results/probe_oldw4_head3442.out`, `results/ref_oldw4.out`):

```
OLD  non-finite 35; finite: value 8.52e-08, gradients 9.17e-07, lowest log density -6600
     non-finite at (t, a): [((0.001, 2.0), 5), ((0.001, 4.0), 10), ((0.001, 8.0), 10), ((0.01, 8.0), 10)]
     true log density there: -714.6 to -2.605e+04
LPDF non-finite 0; finite: value 1.75e-09, gradients 2.76e-12, lowest log density -2.605e+04
```

The series lengths of the older function come from an absolute error
bound (1e-6) on the standardized density, but no finite point is wrong.
The 35 non-finite points are all below log density -708: the small-time
sum underflows on the natural scale.

## Not fixed

Claims: "the 7-parameter lcdf and lccdf integrate natural-scale values and
return -inf below a cdf or ccdf of about 1e-308"; "With sv, sw and st0
all non-zero, the 3-dimensional integrals reach the evaluation limit of
6000 (errors up to 4.3e-5, as in develop)."

The first: points l3 to l5 of §"7-parameter lcdf and lccdf".  The second:
rows `dim 3` of `analyze_evlimit.py` (4.3e-05 in both trees).  The
evidence that the limit causes it came from trees with a limit of 20 000,
60 000 and 200 000 (Triton jobs 20877303, 20877547: lcdf value 1.7e-6,
8.2e-8, 4.7e-9).  These trees are not commits; to repeat, change the
`6000` in the 7-parameter lcdf and lccdf and run `probe_cdfacc.cpp`.

## Diff size

Claim: "code 7 files, +542/−149 (most in `wiener4_lccdf_unnorm.hpp`
+227/−49 and `wiener_full_lpdf.hpp` +162/−35); tests 6 files, +454/−10
(4 new files, 198 lines; `hcubature_test.cpp` +47; `wiener_full_test.cpp`
+209/−10)."

```sh
git -C $PR diff --shortstat 43bc5684a3 ad0b7460d4 -- stan test
git -C $PR diff --numstat 43bc5684a3 ad0b7460d4
```

Checked on 2026-10-28: all numbers agree.

## Tests that fail without the change

Claim: "Each new test fails with the develop versions of the 7 headers."

Result of 2026-10-29 (a stan-test run on Triton, g++ 14.3, CPU): develop
`510d263078` + the branch diff + an uncommitted test change in
`rev/prob/wiener_full_test.cpp` (the new rows of §"Candidate test rows
at the merge base and the head").  Jobs 21275336 (suites and swap),
21275606, 21275387, 21275337 (checklist).

- 59 of 59 Wiener and `hcubature` suites pass.
- With the develop versions of all 7 changed headers, 14 of 16 new or
  changed test cases fail.  `logMill_continuous_at_zero` and
  `logMill_large_drift_finite` pass (regression checks).
- cpplint, `test-math-dependencies`, header checks, `test-headers`,
  doxygen 1.9.7 pass.

The swap test is not in this directory: it copies the 7 headers of
develop into the test tree and runs the new test cases.

## References of the new test rows

The new rows of `rev/prob/wiener_full_test.cpp` use the digits of route
2, `refs/ref_w7_mpmath.py` (mpmath 1.3.0 at 30 digits, central
differences; Triton job 21272844, node skl16, 2 h 24 min):

```sh
while read fn id p; do python3 refs/ref_w7_mpmath.py $fn $id $p; done \
  < data/ref7_testrows_points.txt > out/ref7_testrows_mpmath.out
```

Expected: `data/ref7_testrows_mpmath.out` (the job added the run time
per row as the last column).  Route 1 is `w7/probe_ref7.cpp`: the E
points from `data/ref7_E_D.out`, the D points from
`data/ref7_D_cdf_fd_r*.out`, and k3.1, k52 from
`data/ref7_testrows_k_gk.out` (ids l0, l1; job 21273566):

```sh
STAN_MATH=out/base3442 ./build.sh w7/probe_ref7.cpp out/probe_ref7.bin
REF_CDF_FD=1 REF_LIST=3.1/2/0.1/0.5/5/0.01/0/0:52/2/0.1/0.5/5/0.01/0/0 \
  out/probe_ref7.bin > out/ref7_testrows_k_gk.out
```

Largest difference over the value and the 8 gradients (relative to
max(|ref|, 1)):

| row | route 1 vs 2 | head vs route 2 |
|---|---|---|
| lpdf e65 | 1.0e-8 (d/dsw) | 4.6e-9 |
| lpdf e100 | 2.5e-9 | 1.7e-8 |
| lpdf e113 | 4.1e-10 | 3.0e-11 |
| lpdf e84 | 1.6e-9 | 1.8e-9 |
| lpdf e87 | 3.7e-12 | 2.2e-7 (precision 1e-4) |
| lcdf d65, d178 | 3.9e-10, 5.3e-10 | 3.2e-8, 6.3e-11 |
| lccdf d65, d178 | 1.4e-10, 1.8e-10 | 4.1e-9, 1.6e-11 |
| lccdf k3.1, k52 | 5.5e-10, 2.2e-9 | 4.9e-12, 1.6e-13 |

## Candidate test rows at the merge base and the head

```sh
probe w7/probe_pr3442.cpp
```

Expected: `results/probe_pr3442_{base,head}3442.out` (Triton job
21272749, node pht100, g++ 14.3, `-O2`; log
`results/probe_tree_21272749.log`).  The job ran the probe before its
header comment changed; the code is the same.  Errors: value absolute
on the log scale; gradients relative to max(|ref|, 1); the largest over
the value and the 8 gradients.  References: route 1 (§"References of
the new test rows").

| row | prec | base | head | note |
|---|---|---|---|---|
| lpdf e65, e81, e100, e113 | all | -inf, NaN gradients | <= 2.0e-8 | P1, issue regions |
| lpdf e84, e88 (st0 = 0) | all | d/dsw 2.0e3 | <= 2.5e-10 | P2 |
| lpdf e87, e91 (sw, st0) | all | d/dst0 5.7, 7.8 | <= 2.2e-7 | P2; base 19 ms, head 0.8 ms |
| lcdf d65 (sv, st0) | 1e-8 | d/dw 2.9e-3 | 3.2e-8 | P3 |
| lcdf d178 (sv, sw) | 1e-8 | 2.0e-3 | 5.3e-10 | P3 |
| lccdf d65 | 1e-8 | value 0.83 | 4.1e-9 | P3; missing mass |
| lccdf d178 | 1e-8 | 1.5e-3 | 1.8e-10 | P3 |

E points: `y = 0.2 + d`, `t0 = 0.2`, `sv = 0.1`, `sw = 0.5 min(w, 1 - w)`;
e65 `d = 10, a = 0.1, v = -30, w = 1e-4, st0 = 1e-3`; e81 `d = 10, a = 8,
v = -30, w = 1e-4, st0 = 1e-3`; e100 `d = 50, a = 0.1, v = -10, w = 1e-4,
st0 = 0`; e113 `d = 50, a = 8, v = -30, w = 1e-4, st0 = 1e-3`; e84 `d = 10,
a = 8, v = -10, w = 1e-4, st0 = 0`; e87 `d = 10, a = 8, v = -10, w = 0.5,
st0 = 1e-3`.

P4, lccdf with drift variability only (`a = 2, t0 = 0.1, w = 0.5, v = 5,
sv = 0.01`), default precision 1e-8:

| y | base | head log ccdf | head time |
|---|---|---|---|
| 3.1 | -inf | -39.052352532656435 | 2.3 ms |
| 9, 22.1 | -inf | -119.99567604, -299.40720933 (= test references) | 2 ms |
| 51 | -inf | -693.703262078 | 2.1 ms |
| 52 | -inf | -707.309810359 (ccdf 6.6e-308) | 2.3 ms |
| 54 | -inf | -734.515508043 (ccdf below `DBL_MIN`) | 4.5 ms |
| 60 | -inf | -inf, partials 0 | 31 us |

P5, the two calls below the smallest double:

| call | base | head |
|---|---|---|
| lccdf y = 60 (above) | -inf; partials NaN, +-inf; 25 ms | -inf; partials 0; 31 us |
| lcdf y - t0 = 1e-6 (a = 0.1, w = 1e-4, v = -10, sv = 0.1) | -inf; all partials NaN | -inf; partials 0 |

Details of P2 (default precision): at e84 the base gives d/dsw 12116.7
(reference 6.033) and d/dsv 1546.0 (reference 963.04).  At e87 the base
gives d/dst0 -98.15 (reference 20.728).  The head tree with the develop
`hcubature.hpp` (job 21275329) gives e87 d/dst0 20.728069060, correct;
so the e87 error does not come from the `hcubature` defects.

## Changed test row (WienR)

Claim: "WienR -33.8123, mpmath -33.7956003."

```sh
python3 -c 'import sys; sys.path.insert(0, "refs"); import ref_w4_mpmath as r
from mpmath import mpf, mp
x = [mpf("2.2"), mpf("2.4"), mpf(0), mpf("0.9"), mpf("-4.9"), mpf(0)]
print(mp.nstr(r.lccdf4(*x), 12), [mp.nstr(g, 10) for g in r.grads(r.lccdf4, *x, False)])'
```

Expected: `-33.7956002778` and the gradients y, a, t0, w, v
`-12.88997018 0.3560443568 12.88997018 2.102385235 11.40027803` (the
values in `rev/prob/wiener_full_test.cpp`).  Checked on 2026-10-28.

## Data files

| file | content | made with |
|---|---|---|
| `ref_w4_p4p5.out` | mpmath references, p4 and p5 | `refs/ref_w4_mpmath.py` |
| `ref7_D_lpdf.out` | 7-parameter lpdf, D set | `probe_ref7.cpp REF_D=1` |
| `ref7_E_D.out` | lpdf, extreme grid and 4 D points | `probe_ref7.cpp` (job 20868511) |
| `ref7_D_cdf_fd_r0..3.out` | lcdf, lccdf, D set | `REF_D=1 REF_CDF_FD=1 REF_RANGE=...` |
| `ref7_p5_cdf_fd.out` | lcdf, lccdf, p5 set | `REF_P5=1` |
| `ref7_swst0_s0..3.out` | lpdf, 240-point set | `REF_SWST0=1 REF_SCORE=1 REF_RANGE=...` |
| `ref7_scan_score.out` | lpdf, #3329 scan | `REF_SCORE=1 REF_LIST=...` |
| `ref7_test_dsw.out` | lpdf d/dsw, test points | `REF_SCORE=1 REF_LIST=...` |
| `ref7_underflow_cdf_fd.out` | lccdf at log ccdf -120, -299 | `REF_CDF_FD=1 REF_LIST=...` |
| `ref7_testrows_points.txt` | the new test rows, one per line | input of `refs/ref_w7_mpmath.py` |
| `ref7_testrows_mpmath.out` | the new test rows, route 2 | `refs/ref_w7_mpmath.py` (job 21272844) |
| `ref7_testrows_k_gk.out` | lccdf k3.1, k52, route 1 | `REF_CDF_FD=1 REF_LIST=...` (job 21273566) |

The `probe_ref7.cpp` files were built against the develop tree.  All
files in `data/` are below 200 KB.  The probe outputs of the earlier
tree states are not included.

`results/` holds the outputs of the reruns at `43bc5684a3` and
`ad0b7460d4`: `<probe>_<tree>.out` and `.err`, the build logs (empty: no
warnings), the job logs `pr3442_acc_<job>.log`,
`pr3442_speed_21272931.log` and `probe_tree_21272749.log`, and the
benchmark outputs `probe_bench{,4}_bench.out`.  The largest files are
`probe_w4_{base,head}3442.out` (540 KB each).
