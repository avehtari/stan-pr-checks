# Reproduction scripts for stan-dev/math#3434

The logs in `results/` name the PR tree `edffbaf58a`.  That was a
local merge commit, later amended to `17b0b961ec`.  It has the same
`stan/` tree (`b1db0a3085`) as the PR head `64407093ba`, so every
result here holds for `64407093ba`.

This directory holds the scripts behind the numbers in the text of
stan-dev/math#3434 ("Stable beta_binomial_lpmf at large shapes" and 23
other functions).  Each section below names one claim of the PR text, the
commands that measure it, the expected output, and the measured commits.

Last run: 2026-10-28, on a laptop (12 cores, Ubuntu 22.04, g++ 11.4.0,
`-O2`, Python 3.10.12, mpmath 1.3.0).  The speed benchmarks ran on CSC
Roihu (job 2625261, §"Speed"), the swap tests on Triton on 2026-10-29
(§"Swap tests").

## Commits

| name in this file | commit | what |
|---|---|---|
| develop | `510d263078` | merge base of the merged head (develop of 2026-10-28) |
| PR | `64407093ba` | PR head after the merge of develop (2026-10-28); the probes ran on its tree |
| PR head | `5005dcaa82` | `64407093ba` with unused includes removed and comments changed; the swap tests ran on its tree |
| old head | `c962588c6d` | PR head that stevieb-bot checked (merge base `d12f488749`) |
| S6 pair | `0377416157`, `1c8afb741a` | parent and commit of the `log(mu + phi)` change |

## Layout

| path | content |
|---|---|
| `build.sh` | compiles one probe against one Stan Math tree |
| `refs/` | mpmath references, two routes each; the check of the test constants |
| `sweep/` | the worst-error table (P1): points, probes, scorers |
| `dd/` | `digamma_diff` over 7 281 points |
| `issue1578/` | the point of #1578 |
| `small_phi/` | the negative binomial family at small `phi` |
| `nb2_bitwise/` | `neg_binomial_2_lpmf` before and after the `log(mu + phi)` change |
| `student_t_series/` | the `student_t_lpdf` series against mpmath |
| `propto/` | the `propto` probe of 14 densities; the GLM `propto` error |
| `bnb_grid/` | `beta_neg_binomial_cdf` on a grid of 2 352 points |
| `lbr_select/` | `log_beta_ratio` without the selection of its form |
| `speed/` | speed benchmarks and their job scripts |
| `review_rows/` | the probes behind the test rows added after the review |
| `swap/` | job scripts of the swap tests before the merge of develop |
| `results/` | the outputs of the runs of 2026-10-28 and 2026-10-29 |

All commands run from this directory (`pr3434/`).  New outputs go to
`out/`, which git does not track; compare them with `results/`.

## Requirements

- Python 3 with mpmath.  Tested with mpmath 1.3.0.
- g++ with C++17.  Tested with g++ 11.4.
- A clone of stan-dev/math with the PR branch, and `stan/` of each commit:

  ```sh
  git clone git@github.com:stan-dev/math.git stan_math
  cd stan_math
  git fetch origin pull/3434/head:stable-beta-binomial-lpmf
  make -f make/standalone math-libs   # builds lib/tbb
  for c in 510d263078 64407093ba c962588c6d 0377416157 1c8afb741a; do
    mkdir -p ../trees/$c
    git archive $c stan | tar -x -C ../trees/$c
  done
  git archive 64407093ba test/unit/math | tar -x -C ../trees/64407093ba
  cd -
  ```

Set these variables once:

```sh
export T=$PWD/../trees                      # the stan/ trees
export STAN_MATH_LIB=$PWD/../stan_math/lib  # lib/ for every tree
export TBB_LIB=$STAN_MATH_LIB/tbb           # a directory with libtbb.so
export OPT=-O2
mkdir -p out
build() { STAN_MATH=$T/$1 ./build.sh $2 out/$3; }
```

## P1: the worst-error table

Claim: the table "develop, worst error at large shape" against "this PR"
(13 rows), and "The shape is `1e15` for the first nine rows and at most
`1e15` for the last four".

```sh
for t in dev:510d263078 pr:64407093ba; do
  L=${t%%:*}; C=${t#*:}
  build $C sweep/probe_sib.cpp probe_sib_$L
  build $C sweep/probe_st.cpp probe_st_$L
  build $C sweep/probe_bb.cpp probe_bb_$L
  for f in sib sib_extra sib_slow sib_nb2; do
    out/probe_sib_$L sweep/points_$f.txt > out/${f}_$L.txt
  done
  out/probe_st_$L sweep/points_st.txt > out/st_$L.txt
  out/probe_bb_$L sweep/points_bb.txt > out/bb_$L.txt
done
python3 sweep/score_p1.py dev pr --regions > out/p1_table.md
```

Run time: 3 min of builds, 1.5 min of probes (the 3 points of
`points_sib_slow.txt` take 6 s each on develop), 3 s of scoring.

Points (hex doubles, so C++ and Python read the same values):

| file | points | functions | regions |
|---|---|---|---|
| `points_bb.txt` | 13 177 | `beta_binomial_lpmf` | 14 regions, shapes 1e-2 to 1e18, `N` 1 to 1e6 |
| `points_sib*.txt` | 1 608 | 18 functions of the families | shape decades 1e0 to 1e15 per configuration |
| `points_st.txt` | 3 544 | `student_t_lpdf`, both `multi_student_t`, `yule_simon_lpmf`, `lbeta` | regions B1 to B8 of `gen_points_st.py` |

The scorer reports, per function and tree, the worst value error
(absolute) and the worst log-scale gradient error `|x (d/dx stan - d/dx
ref)|` over every point of every region, the relative error at that
point, the worst relative errors, and the counts of throws, `NaN` and
`inf`.  `--regions` adds the worst errors per region.  Expected
(`results/p1_table.md`), the main columns:

| function | develop value | develop grad | PR value | PR grad |
|---|---|---|---|---|
| `beta_binomial_lpmf` | 2.2e+02 | 7.7e+03 | 7.1e-09 (rel 2.8e-16) | 3.5e-10 (rel 2.3e-14) |
| `beta_neg_binomial_lpmf` | 1.0e+01 | 1.3e+01 | 2.6e-02 (rel 4.9e-17) | 8.2e-02 (rel 2.8e-16) |
| `neg_binomial_lpmf` | 2.3e-13 | 3.0 | 2.3e-13 | 1.0e-14 |
| `neg_binomial_2_lpmf` | 5.6e-13 | 2.7e+01 | 5.6e-13 | 1.7e-14 |
| `neg_binomial_2_log_lpmf` | 9.6e-13 | 2.6e+01 | 9.6e-13 | 2.0e-13 |
| `neg_binomial_2_log_glm_lpmf` | 2.3e+01 | 7.1 | 4.5e-13 | 2.8e-14 |
| `dirichlet_multinomial_lpmf` | 1.0e-12 | 2.8 | 1.0e-12 | 8.1e-15 |
| `lkj_corr_lpdf` | 3.9 | 8.2 (rel 1.5e-14) | 1.5e-01 (rel 2.7e-16) | 1.4e-01 (rel 2.6e-16) |
| `lkj_corr_cholesky_lpdf` | 3.9 | 8.4 (rel 1.5e-14) | 2.4e-02 | 2.4e-02 |
| `yule_simon_lccdf` | 7.6 | 2.1 | 1.7e-13 | 8.2e-15 |
| `beta_binomial_cdf` | 2.1e-01 | 1.9e+01 | 4.1e-08 | 1.1e-06 |
| `beta_binomial_lcdf` | `NaN` (16) | 4.1e+04 | 7.7e-08 | 2.3e-06 |
| `beta_binomial_lccdf` | 3.7e-01 | 2.2e+01 | 1.8e-07 | 2.1e-06 |
| `beta_neg_binomial_cdf` | 3.3e-01 | 1.5e+01 | 2.0e-09 | 2.0e-06 |
| `beta_neg_binomial_lccdf` | 2.9e-01 | 1.3e+01 | 2.0e-09 | 3.4e-06 |
| `student_t_lpdf` | 1.6e-01 | 1.4e-01 | 1.3e-11 | 4.5e-12 |
| `multi_student_t_cholesky_lpdf` | 1.8e-03 | 8.9e-05 | 8.1e-15 | 4.5e-15 |
| `multi_student_t_lpdf` | 1.6e-03 | 2.2e-07 | 5.6e-15 | 6.3e-15 |
| `yule_simon_lpmf` | 3.9e-03 | 3.6e-03 | 6.9e-10 (rel 2.1e-16) | 1.7e-10 (rel 1.9e-16) |
| `lbeta` (`var`), `x d/dx` | | 6.3e-01 | | 1.0e-14 |

Where an absolute error of the PR is above 1e-12, the value or the
gradient is large (`N = 1e6`, shapes of 1e15, counts up to 1e6), and the
relative error is at the rounding level.  The `points_bb.txt` shapes go up
to 1e18; with shapes up to 1e15 only, develop's worst is 0.19 (value) and
7.6 (gradient) over 11 086 points.

References: `refs/refs_bb.txt` (route A of `refs/ref_bb.py`, 80 digits),
`refs/refs_sib.txt` (`refs/ref_sib.py`, routes A, B and C) and
`refs/refs_st.txt` (`refs/ref_st.py`, routes A and B).  The checks of the
routes:

```sh
python3 refs/check_bb_routes.py sweep/points_bb.txt refs/refs_bb.txt 10
python3 sweep/score_sib.py refs/refs_sib.txt out/sib_pr.txt | head -30
python3 refs/ref_st.py sweep/points_st.txt out/refs_st.txt 8
```

Expected: `check_bb_routes.py` at most 4.9e-45 in every region (the 45
stored digits; route C, exact sums, for `N <= 2000`);
`results/sib_certification.md` at most 2.9e-40 between routes A and B;
`ref_st.py` at most 3.4e-26 (value) and 2.1e-25 (gradient), relative
(absolute below 1e-40, where a gradient is exactly 0).  `ref_st.py`
takes 2 s on 8 processes.

## S1, P5: the model of #3153

Claim: "With this PR the value is within `2.3e-14` of the exact value and
the gradient within `3.7e-14`, for `lc` from 0 to 44", and the table of
`lc` 20, 28, 32, 42.  `beta_binomial_lpmf(57 | 117, s, s)`,
`s = exp(lc) / 2`, gradient `d/dlc = s (d/dalpha + d/dbeta)`.

The points are region `REP` of `sweep/points_bb.txt` (`lc` = 0, 0.25,
..., 44; 177 points).

```sh
python3 sweep/score_bb.py out/bb_dev.txt out/bb_pr.txt > out/score_bb.txt
grep -A16 '^REP:' out/score_bb.txt; grep 'REP bounds' out/score_bb.txt
```

Expected (`results/score_bb.txt`): on the PR, value error at most
4.30e-14 (at `lc = 25.5`) and `d/dlc` error at most 3.69e-14 (at
`lc = 7.25`).  The old head `c962588c6d` gives the same bounds.  The
table rows reproduce: develop diff +4.4701e-03 and `d/dlc` -1.4027e-01 at
`lc = 32`; PR -6.5937e-13 and 6.8570e-13 there.

## S2: the reference constants of the tests

Claim: "All references are mpmath at 80 to 140 digits, and two
independent routes agree".

```sh
python3 refs/check_test_refs.py $T/64407093ba > out/test_refs.txt
```

The script reads the test files and recomputes every constant by two
routes (see its header): `large_shapes_test.cpp` (92 rows), the prim and
rev `beta_binomial` tests (12 values; 8 rows with log-shape gradients;
6 `lc` gradients), the OpenCL `beta_binomial` test (7 rows), the
`digamma_diff` prim test (8 values) and OpenCL device test (19 rows), and
the rows of the three negative binomial OpenCL tests (21 rows).  Every
other long literal of the six OpenCL test files must equal a checked
constant.  Run time 1 s.  Expected (`results/test_refs.txt`): routes
agree to at most 3.4e-26 relative (the `student_t` family references are
at 50 digits, route B at 110), every constant is within 1.2e-16 of route
A, and 0 literals are unchecked.

`refs/original/` holds the scripts that chose the test points
(2026-10-07).  They read probe outputs of that time and do not run here.
`refs/make_small_phi_rows.py` and `refs/make_dd_refs.py` run as they are.

## S3: #1578

Claim: the example of #1578 is already correct on develop; "the error of
`alpha * d/dalpha` on develop reaches 3.0 at `alpha = 1e15`, and with
this PR it is `1e-14`".

```sh
python3 refs/ref_sib.py issue1578/points_1578.txt out/refs_1578.txt
for t in dev pr; do out/probe_sib_$t issue1578/points_1578.txt > out/i1578_$t.txt; done
python3 issue1578/score_1578.py out/refs_1578.txt out/i1578_dev.txt out/i1578_pr.txt
```

`refs/refs_1578.txt` holds the references of this run (same content).

Expected (`results/score_1578.txt`): at `n = 0`, `alpha = 1e10 + 1`,
`beta = 0.8` both trees give `d/dalpha = -0.8109302162163288` and
`d/dbeta = 6944444445.13889` (relative errors below 1.5e-16).  At `n = 3`,
`alpha = 1e15`, `beta = 0.8` the error of `alpha * d/dalpha` is 3.0 on
develop and 0.039 on the PR; `alpha * d/dalpha` is about -8.1e14 there,
so 0.039 is a relative error of 4.8e-17.  The value 1e-14 holds with the
mean `alpha / beta` fixed: rows `NB` of P1, PR worst 1.0e-14
(`mu3:1e8`), develop worst 3.0 (`mu3:1e15`).  `issue1578/probe_1578.cpp`
is the probe of 2026-10-08 (prints the values only).

## S4, P2: `beta_binomial_lcdf` on develop

Claim: "On develop this function returned `NaN` at 16 of 81 test points
and was wrong by up to 61."

The 81 points are the `BBLC` lines of `sweep/points_sib.txt`
(`n` 0, 1, 5, 29, 57; `N = 117`; tied shapes and shapes in ratio 1:3,
decades 1e0 to 1e15).  P1 counts them: develop 16 `NaN` and 4 throws
(`alpha = beta = 1`, both trees), PR 0 `NaN`.  The worst finite error of
develop is 74.4 (`tied:1e15`, `n = 1`), not 61; the develop output at
`510d263078` is bitwise equal to the one of 2026-10-07.

## S5: small `phi`

Claim: "At `n = 1000`, `phi = 1e-8`, the error is 7.9e-7 (5.0e-9 in
`neg_binomial_lpmf`).  At `n = 1`, `phi = 1e-17`, the value is `-inf`
... the errors at these points are at most 5.3e-13."

```sh
python3 refs/ref_sib.py small_phi/points_small_phi.txt out/refs_small_phi.txt
for t in dev pr; do
  out/probe_sib_$t small_phi/points_small_phi.txt > out/small_phi_$t.txt
  python3 sweep/score_sib.py out/refs_small_phi.txt out/small_phi_$t.txt
done
```

`refs/refs_small_phi.txt` holds the references of this run.

Points: `neg_binomial_lpmf`, `neg_binomial_2_lpmf`,
`neg_binomial_2_log_lpmf` and the GLM (one observation) at
`(n, phi) = (1000, 1e-8), (1, 1e-17), (100, 1e-6)`, mean `n`.  Expected
(`results/small_phi_*.md`): develop 7.9e-7 (NB2, NB2L), 5.0e-9 (NB),
`-inf` at `phi = 1e-17` (NB, NB2, NB2L); the GLM on develop at most
1.3e-12.  PR at most 5.3e-13 (NB2 at `phi = 1e-8`).

## S6: `log(mu + phi)` once, bitwise identical

Claim: "The results are bitwise identical."

```sh
for c in 0377416157 1c8afb741a; do
  build $c nb2_bitwise/nb2_bitwise.cpp nb2_bitwise_$c
  out/nb2_bitwise_$c > out/nb2_bitwise_$c.txt
  build $c nb2_bitwise/bench_nb2.cpp bench_nb2_$c
  out/bench_nb2_$c L 0 > out/bench_nb2_values_$c.txt
done
cmp out/nb2_bitwise_0377416157.txt out/nb2_bitwise_1c8afb741a.txt
cmp out/bench_nb2_values_0377416157.txt out/bench_nb2_values_1c8afb741a.txt
```

`nb2_bitwise.cpp` prints the value and every adjoint in hex for 2 000
instances (`phi` from 0.14 to 148 and from 1e-8 to 1e15), with `mu` var,
`phi` var, both, a scalar `phi`, `double`, and each instance alone (12 006
lines).  `bench_nb2.cpp` with 0 rounds prints only its `VALUE` lines (the
check of 2026-10-20).  Expected: both `cmp` calls print nothing.

## S7, P6: `digamma_diff` over 7 281 points

Claim: "Over 7 281 points, `x` from `1e-3` to `1e18` and `d` from `1e-6`
to `1e7`, the error is at most 4.2 ulp; the plain difference reaches
`6.5e16` ulp."

```sh
build 64407093ba dd/probe_dd.cpp probe_dd_pr
out/probe_dd_pr dd/points_dd.txt > out/dd_pr.txt
python3 dd/score_dd.py score out/dd_pr.txt
```

`dd/points_dd.txt` comes from `python3 dd/score_dd.py gen` (seed
20261007).  Route A is `mp.digamma(x + d) - mp.digamma(x)` with 80 digits
plus the digits that cancel; route B is `mp.diff` of the `loggamma`
difference (skipped where `x + d` overflows or `x > 1e100`).  Run time
11 s.  Expected (`results/score_dd.txt`): worst 4.24 ulp (`x = 7.723`,
`d = 10`), plain difference 6.52e16 ulp, routes within 4.1e-81.

## S8: the `student_t_lpdf` series

Claim: "Against mpmath, from `h = 20` to `1e15`, the errors are at most
0.9 and 2.4 eps."

```sh
build 64407093ba student_t_series/probe_series.cpp probe_series_pr
out/probe_series_pr > out/series_pr.txt
python3 student_t_series/score_series.py out/series_pr.txt
```

The probe calls `internal::student_t_nu_constant_fun` and
`internal::student_t_nu_digamma_fun` at 2 012 log-spaced `h = nu / 2`.
Expected (`results/score_series.txt`): constant at most 0.51 eps,
digamma difference at most 1.76 eps; routes agree to 7e-45.  The
measurement of 2026-10-07 (0.89 and 2.4 eps) was of the candidate code
in `speed/bench_comp.cpp`, whose file did not keep the accuracy part.

## S9: the `propto` probe, 456 values

Claim: "matches the `propto` and full values of develop at all 456
points on the current head".

```sh
for t in dev:510d263078 pr:64407093ba; do
  build ${t#*:} propto/propto_probe.cpp propto_${t%%:*}
  out/propto_${t%%:*} > out/propto_${t%%:*}.txt
done
python3 propto/compare_propto.py out/propto_dev.txt out/propto_pr.txt
```

Expected: `TOTAL mismatches: 0` (tolerance 1e-12 absolute or 1e-13
relative).

## S10: `beta_neg_binomial_cdf` on a grid

Claim: "`beta_neg_binomial_cdf` then returns 0 with `NaN` gradients (15 of
2352 points of a wide grid), where develop returns a wrong finite value
with finite gradients."

```sh
for t in dev:510d263078 pr:64407093ba; do
  STAN_MATH=$T/${t#*:} EXTRA_CXXFLAGS=-DBNB_NO_PARTS \
    ./build.sh bnb_grid/probe_bnb.cpp out/probe_bnb_${t%%:*}
done
mkdir -p out/bnb && split -n r/12 -d bnb_grid/points_bnb.txt out/bnb/chunk_
for t in dev pr; do
  for c in out/bnb/chunk_??; do out/probe_bnb_$t $c 2 1 > $c.$t.out & done
  wait; cat out/bnb/chunk_??.$t.out > out/bnb_$t.txt
done
python3 bnb_grid/count_bnb.py out/bnb_dev.txt out/bnb_pr.txt
```

The grid is `bnb_grid/gen_points_bnb.py` (alpha 0.25 to 1e4, `r` and
`beta` 0.5 to 1e6, `n` 0 to 20 000).  The third argument 1 runs the cdf
only; each call runs in a child process with a time limit of 2 s
(argument 2; the run of 2026-10-09 used 20 s and found calls of up to
10 s).  Run time 6 min for both trees on 12 cores.  Then run the points
whose outcome differs between the trees and where one tree hit the limit
again with 30 s (16 s):

```sh
python3 - <<'EOF'
import sys; sys.path.insert(0, "bnb_grid"); import count_bnb as C
d, p = C.load("out/bnb_dev.txt"), C.load("out/bnb_pr.txt")
with open("out/bnb_recheck_points.txt", "w") as f:
  for k in d:
    if d[k][0] != p[k][0] and "limit" in (d[k][0], p[k][0]):
      f.write("BNB grid " + " ".join(float(x).hex() for x in k) + "\n")
EOF
for t in dev pr; do out/probe_bnb_$t out/bnb_recheck_points.txt 30 1 > out/bnb_recheck_$t.txt; done
python3 bnb_grid/count_bnb.py out/bnb_recheck_dev.txt out/bnb_recheck_pr.txt
```

Expected (`results/count_bnb.txt`, `results/count_bnb_recheck.txt`):
develop finite and the PR 0 with `NaN` gradients at 14 points of the
first run plus 1 of the recheck, `(0, 30, 10, 300)`: 15 points.  Develop
returns values from -7.2e-11 to 5.4e-11 there.  806 points hit the 2 s
limit on both trees and are not decided by this run (on 2026-10-09, with
20 s, 696 of them threw and 110 finished on both trees).  To decide them,
run all points with the limit 20; that takes about 22 min on 12 cores.

## P3: `log_beta_ratio` without the selection

Claim: "Without the selection, points with `N = 1e6` and one small shape
lost up to 4e-11 where develop was correct."

```sh
lbr_select/make_stirling_tree.sh $T/64407093ba $T/64407093ba_stirling
build 64407093ba_stirling sweep/probe_bb.cpp probe_bb_stirling
out/probe_bb_stirling sweep/points_bb.txt > out/bb_stirling.txt
python3 lbr_select/score_select.py out/bb_dev.txt out/bb_pr.txt out/bb_stirling.txt
```

The patch replaces `if (size_stirling <= size_lbeta)` by `if (true ||
...)` in `prim/fun/log_rising_factorial_ratio.hpp`.  Expected
(`results/score_select.txt`): 131 points where the Stirling-only tree is
worse than 4 times develop and the PR and above 1e-13, all at `N = 1000`
and `N = 1e6`.  The worst absolute loss is 2.8e-10 (relative 3.7e-11, at
`N = 1e6`, both shapes 2.9e5); 87 of the 131 points have one shape below
1000, with up to 2.4e-10 absolute and 6.1e-11 relative.  The PR and
develop have the same errors at these points.

## P4: GLM `propto` error with data `phi`

Claim: "For the GLM with data `phi`, this brings back a rounding error of
about `eps * phi log(phi)` in the `propto` value."

```sh
for t in dev:510d263078 pr:64407093ba; do
  build ${t#*:} propto/glm_propto.cpp glm_propto_${t%%:*}
  out/glm_propto_${t%%:*} > out/glm_propto_${t%%:*}.txt
done
python3 propto/score_glm_propto.py out/glm_propto_dev.txt out/glm_propto_pr.txt
```

Five observations, `alpha` and `beta` var, `phi` data from 1e3 to 1e15.
Route A is the kept part `sum y theta - (y + phi) log(exp(theta) + phi)`;
route B is the full log pmf minus the dropped terms.  Expected
(`results/score_glm_propto.txt`): routes agree to 5e-91; the PR `propto`
error grows from 5.7e-12 (`phi = 1e3`) to 18 (`phi = 1e15`), below
`eps |propto value|` (7.8e-12 to 38), and develop's errors have the same
size.  The full value of the PR is within 3.5e-13.

## Speed

Run of 2026-10-28 on develop `510d263078` and the PR `64407093ba`: CSC
Roihu job 2625261, node rc4215 (one whole `medium` node,
node-exclusive), AMD EPYC 9965, g++ 11.5.0, `-std=c++17 -O2 -DNDEBUG
-pthread -D_REENTRANT`, 5 interleaved rounds (dev pr), median.  Header
libraries from the clone `lib/` at `64407093ba`; TBB 2020.3 from CmdStan
2.40.0.  The md5 of all 3 830 files of both trees matched the local
`git archive` output.

```sh
S=out/speed && mkdir -p $S/trees/t_dev $S/trees/t_pr
git -C ../stan_math archive 510d263078 stan | tar -x -C $S/trees/t_dev
git -C ../stan_math archive 64407093ba stan | tar -x -C $S/trees/t_pr
(cd $S/trees && find t_dev t_pr -type f | sort | xargs md5sum > ../pr3434_trees.md5)
tar -czf $S/pr3434_trees.tar.gz -C $S/trees t_dev t_pr
tar -czf $S/pr3434_libs.tar.gz -C ../stan_math/lib eigen_5.0.1 boost_1.87.0 \
  sundials_6.1.1/include tbb_2020.3/include
# from ~/proj/wip:
TRITON_HOST=vehtaria@roihu-cpu.csc.fi tools/triton_deploy_submit.sh submit \
  --dir pr3434_speed_1028 --job speed/run_bench_pr3434_roihu.sh \
  speed/bench_glm.cpp speed/bench_dd.cpp speed/bench_st.cpp speed/ir_st.cpp \
  $S/pr3434_trees.tar.gz $S/pr3434_trees.md5 $S/pr3434_libs.tar.gz
python3 speed/summarize.py pr3434_speed_<job>.log
```

`speed/bench_st.cpp` (new) times `student_t_lpdf` in its own binary, for
`nu` 2-20, 2-40 and 1e2-1e8, with a control kernel.  The job also counts
instructions per `student_t_lpdf` call with `ir_st.cpp` (valgrind
cachegrind and `perf stat -e instructions:u`; both work on Roihu and
agree).  Log and summary: `results/pr3434_speed_2625261.log`,
`results/pr3434_speed_2625261_summary.txt`.

The older runs (develop `d12f488749`):

`speed/` holds the benchmarks of the speed section of the PR text:

- `bench_glm.cpp`, `bench_dd.cpp`, `run_bench_glm_roihu.sh`,
  `summarize.py`: the time-ratio table.  Roihu jobs 2464563, 2466177 and
  2466306 (log files here), one whole `medium` node, AMD EPYC 9965,
  g++ 11.5, `-O2 -DNDEBUG`, 5 interleaved rounds, median.  The table of
  the PR text is job 2466306 (develop `d12f488749` against `c962588c6d`).
- `glm_compare.cpp`, `compare_glm.py`: the accuracy check of the GLM
  speed-up (2 240 lines, 112 cases).
- `bench_nb2.cpp`, `run_bench_nb2.sh`: `neg_binomial_2_lpmf` with
  `log(mu + phi)` once.
- `bench_comp.cpp`, `run_bench_comp.sh`, `run_bench_st.sh`,
  `bench_st_vec.cpp`, `ir_st.cpp`, `run_ir_st.sh`,
  `dev_student_t_lpdf.hpp`: `student_t_lpdf` time per part (Triton job
  20819212), time (jobs 20821707, 20821709) and instruction counts
  (valgrind cachegrind on the binaries of job 20821847).

The job scripts take the tree paths from `TREE`, `SM`, `LIB` and `TBB`.

## Swap tests

Run of 2026-10-29 on the tree of `5005dcaa82` against develop
`510d263078`, Triton, CPU builds with g++ 14.3.0 on batch nodes, OpenCL
built on a CPU node and run on a V100.  Output (the `RESULT` lines, the
swapped files and the develop failure lines): `results/swap_20261029.txt`.

Procedure: deploy develop plus the branch diff, run the suites, replace
the swapped files with their develop versions (all at once), build and
run the same suites again, restore the branch files, and run them a third
time.  Files that develop does not have stay in the tree.

- CPU (jobs 21277331, 21277408): the 24 changed non-OpenCL files under
  `stan/` that exist on develop (listed after `swapped:` in the output).
  `prim/fun/digamma_diff.hpp` and `prim/fun/log_rising_factorial_ratio.hpp`
  stay.  Suites: `prim/prob/beta_binomial_test`,
  `rev/prob/beta_binomial_lpmf_test`, `rev/prob/large_shapes_test`,
  `mix/prob/neg_binomial_2_log_glm_lpmf_test`,
  `rev/prob/neg_binomial_2_log_glm_lpmf_test`.  On develop 7 of the 9 new
  tests fail; `large_shapes_propto_dropped_terms` and
  `NegBinomial2LogGLM_vector_matches_instances` pass, by design.
- OpenCL (build job 21277332, GPU job 21277560): the 32 changed files
  under `stan/` that exist on develop, 24 CPU and 8 OpenCL files.  The 4
  new files stay.  On develop, every `opencl_large_shapes_reference` test
  and `opencl_n_outside_support` fail; `lp_propto.val()` is -8.21632 for
  `n = 6, N = 5`, and the `alpha` adjoints are not 0.
  `beta_binomial_lpmf_scalar_shapes_test` does not compile on develop
  (`lbeta(const double&, const double&)`), and the `digamma_diff` tests do
  not compile because develop has no `digamma_diff`.

The scripts of the runs before the merge of `510d263078` are in `swap/`:
`run_bb_swap.sh` (CPU, the 23 develop headers of that time; Triton jobs
20812673, 20814979, 20823134), `run_nb_fix.sh` (the negative binomial
review fixes; jobs 20916196, 20976802), `run_ocl.sh` and `run_nb_ocl.sh`
(OpenCL on a V100; jobs 20823167, 20916534, 20976804).

## Test rows added after the review

`review_rows/` chooses and checks the rows of the review commit
`64407093ba`: the `YSC` and `YSLC` rows and the `NB` row at the point of
#1578 in `large_shapes_test.cpp`, and the `fvar<var>` Hessian test of
`neg_binomial_2_log_glm_lpmf` in
`mix/prob/neg_binomial_2_log_glm_lpmf_test.cpp`.

```sh
for c in 510d263078 64407093ba; do
  STAN_MATH=<tree of $c> ./build.sh review_rows/probe_ys.cpp out/ys_$c
  out/ys_$c > out/review_rows_$c.txt
  STAN_MATH=<tree of $c> ./build.sh review_rows/probe_glm_hess.cpp out/hess_$c
  out/hess_$c
done
python3 review_rows/score_ys.py      # results/review_rows_ys.txt
python3 review_rows/rows_e1_e2.py    # results/review_rows_refs.txt
```

- `yule_simon_cdf` and `_lcdf` are `1 - c` and `log1p(-c)` with
  `c = n B(n, alpha + 1)`.  Two routes: `loggamma` and `digamma`, and the
  exact product `n! / prod (alpha + j)` with `sum 1 / (alpha + j)`; they
  agree to 1e-65 or better.  With the test tolerance
  `tol * max(1, |ref|)` no point fails on develop, so the test compares
  these two tags relative to the reference.  Develop's relative error of
  `alpha d/dalpha`: 5e-7 at `alpha = 1e8`, 4.0e-3 at 1e12, 1.0 at 1e15
  (`n = 1`).  PR: at most 2.3e-15 for `n <= 5`, and 1.3e-13 at `n = 20`,
  where `log c` is about -645 and one rounding of it is 1.4e-13.
- `neg_binomial_lpmf(0 | 1e10 + 1, 0.8)`: the closed form against
  `mp.diff` at 140 digits, agreement 4e-92.  Develop and the PR give the
  same value and `alpha` partial; the `beta` partials differ by 1 ulp.
- `results/review_rows_hessian.txt`: the first two lines are develop
  `510d263078`, the last two the PR (`d^2/dtheta^2` at `y = 1`,
  `phi = 0.07`; reference `-(y + phi) phi mu / (mu + phi)^2`, checked
  against `mp.diff`, agreement 4e-97).  Develop: 0 at `theta = 35` and a
  relative error of 1.4e-3 at 27; PR: at most 5.6e-16.
