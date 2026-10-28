# Reproduction scripts for stan-dev/math#3408

This repository holds the scripts behind the numbers in the text of
stan-dev/math#3408 ("Incomplete gamma gradient roots and log-space cdfs").
Each section below names one claim of the PR text, the commands that
measure it, the expected output, and the commits that were measured.

Last run: 2026-10-28, on a laptop (12 cores, Ubuntu 22.04, g++ 11.4.0,
Python 3.10.12, mpmath 1.3.0).

## Layout

| directory | content |
|---|---|
| `build.sh` | compiles one probe against one Stan Math tree |
| `refs/` | mpmath reference generators (two routes each) |
| `roots/` | the gradient roots: 1 800-point sweep, 49 990 points, the forward-mode hang |
| `cdfs/` | the 18 cdf functions: survey, 324 reference rows, `gamma_lccdf` modes, grid |
| `speed/` | the speed benchmarks and example SLURM jobs (not run again) |
| `swap/` | `swap_test.sh`: run tests with the PR headers and with the base headers |

All commands below run from this directory (`pr3408/`).  Outputs go to
`out/`.

## Requirements

- Python 3 with mpmath.  Tested with mpmath 1.3.0 and 1.4.1.
- g++ with C++17.  Tested with g++ 11.4, 11.5 and 14.3.
- R with ggplot2, only for the heat map of the grid (`cdfs/plot_grid.R`).
- A clone of stan-dev/math with the PR branch:

  ```sh
  git clone git@github.com:stan-dev/math.git stan_math
  cd stan_math
  git fetch origin pull/3408/head:stable-inc-gamma
  git checkout stable-inc-gamma
  make -f make/standalone math-libs   # builds lib/tbb and SUNDIALS
  ```

  `make -f make/standalone math-libs` builds `lib/tbb/libtbb.so`.  The
  probes link to it.

- A develop tree.  The PR text calls the merge base `5252d51d47`
  "develop".  The probes need its `stan/` directory only:

  ```sh
  mkdir -p ../dev_5252d51d47
  git archive 5252d51d47 stan | tar x -C ../dev_5252d51d47
  ```

Set these variables once:

```sh
export PR=$PWD/../stan_math                 # the PR tree (stan/ and lib/)
export DEV=$PWD/../dev_5252d51d47           # stan/ of the merge base
export STAN_MATH_LIB=$PR/lib                # lib/ for both trees
export TBB_LIB=$PR/lib/tbb                  # a directory with libtbb.so
mkdir -p out
```

`build.sh` takes `STAN_MATH` (the tree), and optionally `STAN_MATH_LIB`,
`TBB_LIB`, `CXX`, `OPT` (default `-O3`) and `EXTRA_CXXFLAGS`.  It compiles
with `-std=c++17 -O3 -pthread -DBOOST_DISABLE_ASSERTS` and links
`-ltbb` with an rpath.  Example:

```sh
STAN_MATH=$DEV ./build.sh roots/sweep_roots.cpp out/sweep_roots_dev
STAN_MATH=$PR  ./build.sh roots/sweep_roots.cpp out/sweep_roots_pr
```

A probe that includes `stan/math/mix.hpp` takes about 1 min to compile.
Every run below takes less than 2 min on the laptop.

## Commits

| name | commit | what |
|---|---|---|
| develop | `5252d51d47` | merge base of the PR; "develop" in the PR text |
| PR head | `22d781f290` | the PR head when the text was posted (2026-10-21) |
| merged head | `00217e3f54` | the PR branch after a merge of develop `510d263078` |
| branch of Part 1 | `1a9df72569` | the root change alone (3 files), used for the root tables |

The re-runs of 2026-10-28 used develop and the merged head.  In the merged
head, the `stan/` files of the PR are those of `22d781f290`; the merge
changed only include lists.

## Claim 1: sweep of 1 800 points in nine regions

PR text, Part 1, "Accuracy": the table of the worst relative error per
region, `develop` against this PR.

```sh
python3 roots/gen_points.py 200 out/points.txt     # md5 84ceeb7f76a3125821c80bebb9c20118
STAN_MATH=$DEV ./build.sh roots/sweep_roots.cpp out/sweep_roots_dev
STAN_MATH=$PR  ./build.sh roots/sweep_roots.cpp out/sweep_roots_pr
out/sweep_roots_dev out/points.txt > out/sweep_dev.txt
out/sweep_roots_pr  out/points.txt > out/sweep_pr.txt
python3 roots/score_sweep.py out/sweep_dev.txt 8
python3 roots/score_sweep.py out/sweep_pr.txt 8
```

Columns: `up` is `grad_reg_inc_gamma`, `lo` is `grad_reg_lower_inc_gamma`,
`ei` is Eigen `igamma_der_a`.  On the PR tree the three columns are equal,
because both roots call Eigen.  Expected worst relative errors (develop
run; the `ei` column is the PR):

| region | up worst | lo worst | ei worst (PR) | up non-finite |
|---|---|---|---|---|
| R1 series, a<1, z<8 | 1.49e-03 | 1.51e-10 | 8.51e-15 | 0 |
| R2 series, 1<=a<=10 | 7.67e-03 | 9.60e-03 | 2.13e-15 | 0 |
| R3 series, 10<a<=170 | 2.64e+29 | 5.80e-01 | 1.21e-13 | 0 |
| R4 asymptotic, a<=z<=2a | 5.05e-01 | 1.14e+07 | 1.67e-13 | 0 |
| R5 asymptotic, z>=3a | 3.74e-04 | 2.07e+165 | 5.36e-12 | 0 |
| R6 seam z=8 | 1.23e-02 | 9.85e-12 | 4.35e-15 | 0 |
| R7 integer a | 1.91e+15 | 2.55e+01 | 3.26e-14 | 0 |
| R8 a>171.6 | inf | 1.14e+01 | 3.81e-13 | 176 |
| R9 z>=3a, a>=12 | 9.55e-07 | 1.87e+127 | 2.50e-11 | 0 |

Measured on: develop `5252d51d47` (Triton job 20335653, 2026-09-18) and,
again, develop and the merged head (2026-10-28).  The develop output of
2026-10-28 is byte-identical to the output of job 20335653.

Differences from the PR text, all in the PR column, all from the
reference: R3 1.2e-13 (PR text 1.9e-13), R5 5.4e-12 (3.8e-12), R9
2.5e-11 (4.3e-12).  The C++ output is identical; `refs/ref_mpmath.py` was
changed after the PR table was made (2026-09-20, the tail nodes of the
quadrature route).  At the worst R5 and R9 points, `mp.diff` (route A)
agrees with Eigen to 1.5e-14 and 3.4e-14, and with `ref_mpmath.dQda` to
5.3e-12 and 2.5e-11.  Against `mp.diff` at all 400 points of R5 and R9,
the PR roots have a worst relative error of 4.5e-14 (R5) and 9.4e-14
(R9).  So the values above for R5 and R9 are limits of the reference, not
of the roots.

## Claim 2: 49 990 points, non-finite counts and worst errors

PR text: "`develop` returns NaN or inf at 6 379 of them for
`grad_reg_inc_gamma`; this PR is finite at every one"; worst relative
error outside `a ≈ z ≈ 1e4` 3e-11; near `z = a` 3e-11 for `a < 1 000`,
2.9e-7 for `a < 10 000`, 5.2e-6 at `a = z = 1e4`.

```sh
python3 roots/gen_points_cmp.py 5000 out/points_cmp.txt   # md5 55cecb4a421ef015a9902fc132b5fa01
STAN_MATH=$DEV ./build.sh roots/probe_cmp_roots.cpp out/probe_cmp_roots_dev
STAN_MATH=$PR  ./build.sh roots/probe_cmp_roots.cpp out/probe_cmp_roots_pr
out/probe_cmp_roots_dev out/points_cmp.txt double > out/cmp_double_dev.txt
out/probe_cmp_roots_pr  out/points_cmp.txt double > out/cmp_double_pr.txt
python3 roots/score_cmp_roots.py double out/cmp_double_dev.txt out/cmp_double_pr.txt 8
python3 roots/bin_near_diag.py out/cmp_double_pr.txt
```

The references come from `roots/refs_cache_49990.txt` (a, z, dQ/da in
hex; md5 `a1d4f45450b33f394a25dfb1746c32b5`).  Triton job 20356673
computed it with `refs/ref_mpmath.py`.  To compute it again, delete the
file; `score_cmp_roots.py` then needs about 10 min on 6 cores.

Expected output of `score_cmp_roots.py` (`grad_reg_inc_gamma`; `dev nf` and
`bra nf` count non-finite results at points with a normal reference):

| region | dev worst | dev nf | PR median | PR worst | PR nf |
|---|---|---|---|---|---|
| C1 whole box | inf | 183 | 1.27e-15 | 2.84e-11 | 0 |
| C2 a<1 | 5.48e-03 | 0 | 5.29e-16 | 2.84e-11 | 0 |
| C3 z<<a | inf | 109 | 2.30e-15 | 1.13e-12 | 0 |
| C4 z near a | inf | 2 089 | 1.99e-14 | 2.93e-07 | 0 |
| C5 z>=3a | inf | 18 | 8.53e-16 | 2.50e-11 | 0 |
| C6 seam z=8 | 1.49e-01 | 0 | 8.62e-16 | 2.73e-13 | 0 |
| C7 integer a | inf | 2 827 | 6.49e-14 | 3.13e-11 | 0 |
| C8 a near 171.62 | inf | 406 | 5.92e-14 | 2.93e-11 | 0 |
| C9 a*log z near 709 | inf | 873 | 7.70e-14 | 2.52e-11 | 0 |
| C10 corners | inf | 874 | 2.25e-15 | 5.17e-06 | 0 |

The output is identical to that of Triton job 20356673 (develop
`5252d51d47` against `1a9df72569`).  Expected output of
`bin_near_diag.py`:

| bin | n | worst | median |
|---|---|---|---|
| C4 [1, 10) | 1 217 | 3.19e-14 | 7.18e-16 |
| C4 [10, 100) | 1 290 | 1.21e-13 | 6.02e-15 |
| C4 [100, 1000) | 1 251 | 3.16e-11 | 1.09e-13 |
| C4 [1000, 10000) | 1 083 | 2.93e-07 | 1.83e-12 |
| C10 a, z > 9 000 | 556 | 5.17e-06 | 1.63e-11 |

Notes:

- The `dev nf` column sums to **7 379**, not 6 379 as in the PR text.
  Over all 49 990 points, including those with a subnormal or zero
  reference, develop returns a non-finite `grad_reg_inc_gamma` at 10 812
  points and a non-finite `grad_reg_lower_inc_gamma` at 6 475.  The PR
  returns no non-finite value.
- The PR worst values of about 3e-11 in C1, C2, C5, C7, C8 and C9 are
  limits of the cached reference.  At the worst point of each of these
  regions, the PR is within 1.2e-15 to 9.6e-14 of `mp.diff`, and the cache
  is 2.5e-11 to 3.1e-11 from `mp.diff`.  The C4 and C10 values (2.9e-7 and
  5.2e-6) are real: the cache agrees with `mp.diff` to 1e-18 there.
- The 3.16e-11 of the bin C4 [100, 1000) is a limit of the cache too.  At
  its worst point, `(a, z) = (927.83, 1612.09)`, the PR is within 6.0e-14
  of `mp.diff`, the cache is 3.2e-11 from `mp.diff`, and route C (the
  series `ref_mpmath.dPda_C`) agrees with `mp.diff` to 4e-17.
- At the worst C10 point, `(a, z) = (10152.96, 10152.78)`, the three
  routes agree to 7e-17: the cache (route B), `mp.diff` (route A) and
  route C.  The PR error there, 5.17e-6, is real.

## Claim 3: spread between autodiff orders

PR text: "roots, autodiff-order spread: inf in 7 of 10 regions (develop),
2.2e-14 (this PR)".

```sh
out/probe_cmp_roots_dev roots/points_capped.txt all > out/cmp_all_dev.txt
out/probe_cmp_roots_pr  roots/points_capped.txt all > out/cmp_all_pr.txt
python3 roots/score_cmp_roots.py all out/cmp_all_dev.txt out/cmp_all_pr.txt 8
```

`roots/points_capped.txt` holds 2 000 points of `points_cmp.txt`, 200 per
region, all with `z < 500`.  Develop does not terminate in forward mode
for larger `z` (claim 8).  The script that drew this subset was not kept;
the file is the input of Triton job 20355783.  `score_cmp_roots.py` reads
the references of these points from `refs_cache_49990.txt`.

Expected output, last table ("spread between the four autodiff
configurations"):

| region | dev up | PR up | dev lo | PR lo |
|---|---|---|---|---|
| C1 whole box | inf | 4.45e-15 | 0 | 4.45e-15 |
| C2 a<1 | 0 | 1.54e-15 | 0 | 1.54e-15 |
| C3 z<<a | inf | 4.21e-16 | 0 | 4.21e-16 |
| C4 z near a | inf | 1.96e-14 | 0 | 1.96e-14 |
| C5 z>=3a | 0 | 2.10e-15 | 0 | 2.10e-15 |
| C6 seam z=8 | 0 | 1.25e-15 | 0 | 1.25e-15 |
| C7 integer a | inf | 1.53e-14 | 0 | 1.53e-14 |
| C8 a near 171.62 | inf | 9.99e-15 | 0 | 9.99e-15 |
| C9 a*log z near 709 | inf | 1.83e-14 | 0 | 1.83e-14 |
| C10 corners | inf | 2.18e-14 | 0 | 2.18e-14 |

Measured on: Triton jobs 20355783 and 20356673 (develop `5252d51d47`,
branch `1a9df72569`), and again on develop and the merged head
(identical output).

## Claim 4: the 18 cdfs at 324 points

PR text, summary table and Part 3 "Accuracy": develop has a non-finite
value or gradient at 199 of 324 test points, the PR at 0; worst gradient
error where finite 1.6e+93 (develop) and 2.7e-10 (PR); worst value error
2e-13 and gradient error 8e-13, except the `nu` partial of
`scaled_inv_chi_square` (2.7e-10); "the values change only where they
were `-inf` or lost digits".

Two probes.  `ref324_probe.cpp` evaluates the 324 rows of the reference
test of the PR head (`cdfs/incomplete_gamma_cdfs_reference_324.hpp`, a
copy of `test/unit/math/rev/prob/incomplete_gamma_cdfs_reference.hpp` at
`22d781f290`).  `survey_probe.cpp` evaluates 17 of the functions (not
`poisson_cdf`) at the same 18 points (306 cases).

```sh
STAN_MATH=$DEV ./build.sh cdfs/ref324_probe.cpp out/ref324_dev
STAN_MATH=$PR  ./build.sh cdfs/ref324_probe.cpp out/ref324_pr
out/ref324_dev > out/ref324_dev.txt
out/ref324_pr  > out/ref324_pr.txt
python3 cdfs/ref324_score.py develop=out/ref324_dev.txt PR=out/ref324_pr.txt

STAN_MATH=$DEV ./build.sh cdfs/survey_probe.cpp out/survey_dev
STAN_MATH=$PR  ./build.sh cdfs/survey_probe.cpp out/survey_pr
out/survey_dev cdfs/survey_cases.txt > out/survey_dev.txt
out/survey_pr  cdfs/survey_cases.txt > out/survey_pr.txt
cd cdfs && python3 survey_score.py ../out/survey_pr.txt && cd ..
cd cdfs && python3 survey_score.py ../out/survey_pr.txt --detail | grep 2.7e-10; cd ..
```

Expected output of `ref324_score.py` (summary lines):

```text
== develop: 324 rows
   throws 0; rows with a non-finite value or gradient 199; total 199 of 324
   worst finite gradient error 1.64e+93  gamma_lcdf(300.0 | 150.0, 2.0) [upper_tail] ...
   values off by more than 1e-12: 66
== PR: 324 rows
   throws 0; rows with a non-finite value or gradient 0; total 0 of 324
   worst finite value error    2.43e-13  scaled_inv_chi_square_lccdf(...) [bulk_big]
   worst finite gradient error 2.74e-10  scaled_inv_chi_square_lcdf(1.0 | 800.0, 1.0) ...
   values off by more than 1e-12: 0
   develop value non-finite: 28; finite and off by more than 1e-12: 38; within 1e-12 ...: 0
```

The last line checks the claim on the values.  A value changes from
develop to the PR at 66 rows.  At 28 of them develop returns `-inf`.  At
the other 38 develop is off by more than 1e-12: by 1 (it returns 0 for a
log probability of -1e-32 to -1e-114), 1.5e-4 (subnormal probability),
4.2e-8 or 1.1e-9.  At no row does a correct develop value change.

Expected output of `survey_score.py` on the PR: every cell at most
2e-13 (value) and 8e-13 (gradient), except `scaled_inv_chi_square_*`
gradients: 2e-12 (bulk) and 3e-10 (bulk_big).  The `--detail` line gives
`2.7e-10` for the second gradient component (`nu`) at `a = z = 400`.
On develop the survey shows `inf` in most tail cells.

Measured on: the survey on the working tree of `7d3ea2f5c6` (2026-10-21,
NOTE §10), and again on develop and the merged head (2026-10-28).  The
324-row probe is new (2026-10-28); the source of the 199 and 1.6e+93 of
the PR text was not recorded.  Both numbers are reproduced exactly.

## Claim 5: `gamma_lccdf` in all autodiff modes

PR text: "For `gamma_lccdf` in all autodiff modes, first derivatives are
within 5.6e-12 and second derivatives within 1.9e-10 (7e-9 near
`y = alpha + 1` at `alpha = 5`)."

```sh
STAN_MATH=$DEV ./build.sh cdfs/probe_cpp.cpp out/probe_cpp_dev
STAN_MATH=$PR  ./build.sh cdfs/probe_cpp.cpp out/probe_cpp_pr
out/probe_cpp_dev > out/cpp_dev.txt 2> out/cpp_dev_ffd2.txt
out/probe_cpp_pr  > out/cpp_pr.txt  2> out/cpp_pr_ffd2.txt
python3 cdfs/score_cpp.py dev=out/cpp_dev.txt pr=out/cpp_pr.txt
python3 cdfs/score_cpp_ffd2.py dev=out/cpp_dev_ffd2.txt pr=out/cpp_pr_ffd2.txt
```

`score_cpp.py` scores `d/dalpha` in `var`, `fvar<double>`, `fvar<var>`,
`fvar<fvar<double>>`, and the `fvar<var>` second derivative (`d2:`).
`score_cpp_ffd2.py` scores the `fvar<fvar<double>>` second derivative.
Expected, PR tree: first derivatives at most 5.6e-12 (at `alpha = 1e4`,
`y = 10 100`); second derivatives at most 1.9e-10 (same point) except
7.2e-9 and 6.9e-9 at `alpha = 5`, `y = 6.5` and `7.5`; both second
derivative modes give the same errors.  Develop: non-finite first
derivatives at 9 of 17 points; the `fvar<fvar<double>>` second
derivative is 0 at 16 points and NaN at one.

Measured on: the working tree of `8d376d5763` (2026-10-21, NOTE §8), and
again on develop and the merged head (2026-10-28).

## Claim 6: the grid of 266 cells

PR text: "On the grid of the plot by @mjcarter95, all 266 cells are
finite and the worst relative error of the shape derivative is 3.4e-12."

```sh
STAN_MATH=$PR  ./build.sh cdfs/grid_probe.cpp out/grid_probe_pr
STAN_MATH=$DEV ./build.sh cdfs/grid_probe.cpp out/grid_probe_dev
python3 cdfs/plot_grid.py out/grid_probe_pr out/grid_pr.csv
python3 cdfs/plot_grid.py out/grid_probe_dev out/grid_dev.csv
Rscript cdfs/plot_grid.R out/grid_pr.csv out/grid_pr.png "gamma_lccdf d/dalpha, PR"
```

Expected: PR `worst 3.35e-12 non-finite 0 cells 266`; develop `worst
2.92e-01 non-finite 137 cells 266`.  `out/grid_pr.csv` is identical to
`cdfs/gamma_lccdf_dalpha_grid.csv`, the CSV of the figure
`cdfs/gamma_lccdf_dalpha_fixed.png` in the PR thread.  Each run takes
about 50 s (mpmath, two routes per cell).

Measured on: the working tree of `7d3ea2f5c6` (2026-10-21), and again on
develop and the merged head (2026-10-28).

## Claim 7: the points of the review

The points that the review of 2026-10-27 asked about: the series branch
of the old root at `(170, 100)`, `(165, 80)`, `(160, 120)`, `(150, 150)`;
the 1e-6 default at `(0.5, 0.25)` and `(8, 4)` with an exact `g`; the
`gamma_lccdf` shape derivative at `(150, 1200)`, `(10, 1000)`,
`(200, 1600)`; the `fvar` example of #2006; `gamma_lcdf(20000 | 800, 0.1)`.

```sh
STAN_MATH=$DEV ./build.sh cdfs/probe_review.cpp out/probe_review_dev
STAN_MATH=$PR  ./build.sh cdfs/probe_review.cpp out/probe_review_pr
out/probe_review_dev
out/probe_review_pr
python3 refs/ref_review.py
```

Expected:

| quantity | point | reference | develop | PR |
|---|---|---|---|---|
| dQ/da | (170, 100) | 6.6120002024305e-11 | -inf | 6.6120002024311e-11 |
| dQ/da | (165, 80) | 4.7880951404265e-17 | -inf | 4.7880951404264e-17 |
| dQ/da | (160, 120) | — | inf | 8.742371113095304e-05 |
| dQ/da | (150, 150) | 0.0325916061570337 | 0.0100767351829512 | 0.0325916061570312 |
| dQ/da | (0.5, 0.25) | 0.797947316783223 | 0.797947390469200 | 0.797947316783223 |
| dQ/da | (8, 4) | 0.0421220900217987 | 0.0421182882820318 | 0.0421220900217987 |
| lccdf d/dalpha | (150, 1200) | 2.0837288947176 | 2 | 2.0837288947176424 |
| lccdf d/dalpha | (10, 1000) | 4.6570107370619 | NaN | 4.6570107370619187 |
| lccdf d/dalpha | (200, 1600) | 2.0826567480374 | NaN | 2.0826567480373623 |
| #2006 tangent | (3, 36) | -11.500539782310905 | -6704987257582.83 | -11.500539782310895 |
| gamma_lcdf | (800, 2000) | -1.48339094805287e-205 | 0 | -1.48339094805311e-205 |

`ref_review.py` prints both routes; every pair agrees to better than
1e-55.  The `gamma_lccdf` references come from `refs/gen_refs.py`.

Measured on: develop, `b172d616`, `ad3ba698`, `22d781f290` and
`00217e3f54` (2026-10-28).

## Claim 8: forward-mode `gamma_q` does not terminate on develop

PR text: "a call at `a = 0.168, z = 4 748` did not return in 15 minutes".

```sh
STAN_MATH=$DEV ./build.sh roots/probe_diag_slow.cpp out/diag_dev
STAN_MATH=$PR  ./build.sh roots/probe_diag_slow.cpp out/diag_pr
timeout 60 out/diag_dev out/points_cmp.txt 60 50 > out/diag_dev.txt 2> out/diag_dev.err
echo "exit $?"          # 124: the timeout stopped it
tail -n 1 out/diag_dev.err
timeout 60 out/diag_pr out/points_cmp.txt 60 50 > out/diag_pr.txt 2> out/diag_pr.err
echo "exit $?"          # 0
```

The probe writes a marker to stderr before each call.  Expected on
develop: exit 124 after 60 s; the last marker is
`START n=7 upper_fd a=0.16833942110303479 z=4747.6532029474602`
(`grad_reg_inc_gamma` with `fvar<double>`).  Expected on the PR: exit 0
at once, `# examined 60 points`, and `DONE  n=7 upper_fd 0.0 ms ok` in
`out/diag_pr.err`.

Measured on: Triton job 20355549 (develop, 15 min timeout), and again on
develop and the merged head with a 60 s timeout (2026-10-28).

## Claim 9: the new tests fail without the changes

PR text, "Testing": the swap table.  `swap/swap_test.sh` copies the
working tree of a clone to a scratch directory, and runs the given tests
three times: with the PR headers, with the `BASE_REF` version of every
`stan/` file that the PR modifies or deletes, and with the PR headers
again.  It does not write to the clone.

```sh
STAN_MATH=$PR BASE_REF=origin/develop TBB_LIB=$TBB_LIB JOBS=8 \
  swap/swap_test.sh test/unit/math/prim/fun/grad_reg_inc_gamma_accuracy_test.cpp
```

Expected (merged head, 2026-10-28, 25 files swapped):

```text
######## stage pr_1 ########
  PASS  grad_reg_inc_gamma_accuracy_test  (6 cases)
######## stage base ########
  FAIL  grad_reg_inc_gamma_accuracy_test  (exit 1; cases ok 0, failed 6; failed expectations 21)
######## stage pr_2 ########
  PASS  grad_reg_inc_gamma_accuracy_test  (6 cases)
== restore check: swapped files equal stage pr_1 (...)
== STAN_MATH unchanged (...)
```

The table in the PR text came from swaps of the headers of one part at a
time (Triton job 20336069; CSC Roihu jobs 2428840, 2430049 and 2437062).
`swap_test.sh` swaps all files of the PR together.  To swap one part,
use a `BASE_REF` that holds the other parts (for example the commit
before the part).  The suites of the table:

- `test/unit/math/prim/fun/grad_reg_inc_gamma_accuracy_test.cpp`,
  `test/unit/math/mix/fun/grad_reg_inc_gamma_test.cpp`,
  `test/unit/math/mix/fun/gamma_q_test.cpp`
- `test/unit/math/rev/prob/gamma_lccdf_test.cpp`,
  `test/unit/math/mix/prob/gamma_lccdf_test.cpp`
- `test/unit/math/rev/prob/incomplete_gamma_cdfs_reference_test.cpp`
- `test/unit/math/mix/fun/gamma_p_test.cpp`,
  `test/unit/math/prim/prob/poisson_ccdf_log_test.cpp`

The mix suites take several minutes each to compile.  Run them on a
build machine.

## The references of the tests

| test file of the PR | generator | check |
|---|---|---|
| `grad_reg_inc_gamma_accuracy_test.cpp` | `refs/ref_mpmath.py` (`dQda`) | route A against B |
| `gamma_lccdf_tail_reference.hpp` | `refs/gen_refs.py` | output equals the 13 rows |
| `incomplete_gamma_cdfs_reference.hpp` | `refs/gen_cdf_refs.py OUT` | output equals the header |
| survey cases | `refs/survey_cases.py` (run in `cdfs/`) | equals `cdfs/survey_*.txt` |

On 2026-10-28 `gen_cdf_refs.py` wrote a header identical to the header in
the merged head (334 rows: the 324 rows of `22d781f290` and 10 rows of the
root table), in 70 s.  `gen_refs.py` and `survey_cases.py` reproduced
their saved outputs byte for byte.  `refs/probe_lccdf_tail.py` and
`refs/probe_lower_series.py` print the Python prototypes of the log-space
continued fraction and lower series.

`python3 refs/ref_mpmath.py` prints the agreement of route A (`mp.diff`)
and route B (quadrature).  Two known limits:

- At `(12500, 12000)` route B agrees with route A to 4 digits.  The
  substitution `t = z s^(1/a)` hides the density peak at large `a` (see
  the comment above `dQda`).  The public `dQda` uses the series there,
  which agrees with `mp.diff`.
- The incomplete beta section stops in the d/dz table with a mpmath
  `hypsum()` error at `(12500, 12500, 0.5)`.  The beta roots are not part
  of #3408.

## Speed (not run again)

The speed numbers of the PR text were measured once on cluster nodes.
The scripts are in `speed/`.  The SLURM jobs are examples; set the
variables at their top.  `speed/make_tarball.sh` packs the tree and the
develop headers for the two jobs of Part 3.

| PR table | machine | compiler, flags | trees | job |
|---|---|---|---|---|
| Part 1, roots | Triton pe63, Xeon E5-2680 v3 2.50 GHz | g++ (1), `-O2 -DNDEBUG` | (2) | 20356756 |
| Part 3, cdfs | Roihu rc5155, AMD EPYC 9965 | g++ 11.5.0, `-O3` (3) | (4) | 2430495 |
| Part 3, gamma_lccdf | Roihu rc5155, AMD EPYC 9965 | g++ 11.5.0, `-O3` (3) | (4) | 2430495 |

1. GCC of the Triton modules `triton/2025.1-gcc scicomp-r-env`.  The job
   log does not record the version; on 2026-10-28 the modules give g++
   14.3.0.  Exclusive node, `performance` governor, two rounds, spread
   under 3 %.  Script: `speed/run_cmp_speed.sh`, `speed/probe_cmp_speed.cpp`.
2. develop `5252d51d47` against `1a9df72569` (the root change: 3 files).
3. `-O3 -pthread -D_REENTRANT -DBOOST_DISABLE_ASSERTS -DTBB_INTERFACE_NEW`,
   pinned with `taskset -c 2`, 4 alternating runs, medians.  Script:
   `speed/run_roihu_bench.sh` with `speed/bench_cdfs.cpp` and
   `speed/bench_lccdf_tail.cpp`.
4. The PR working tree of 2026-10-21 19:31 EEST (commit `7d3ea2f5c6`,
   before `22d781f290`, which changed only `rev/fun/gamma_p.hpp` and
   `poisson_lccdf`) applied to develop, against develop's version of the
   60 `stan/` files that differ (`make_tarball.sh` with
   `ORIG_REF=origin/develop`).  The log does not print the develop commit;
   the test runs of the same day used develop `e50e79f188`.

`speed/run_bench_lccdf_tail.sh` is Triton job 21111830, an intermediate
comparison of two `gamma_lccdf` variants.  It is not a table of the PR
text.

## Differences from the PR text

| claim | PR text | now | cause |
|---|---|---|---|
| 2: develop non-finite | 6 379 | 7 379 | the sum of the per-region counts is 7 379 |
| 1: PR worst R3 | 1.9e-13 | 1.2e-13 | reference changed after the table |
| 1: PR worst R5 | 3.8e-12 | 5.4e-12 | reference error; 4.5e-14 against `mp.diff` |
| 1: PR worst R9 | 4.3e-12 | 2.5e-11 | reference error; 9.4e-14 against `mp.diff` |
| 2: PR worst outside z = a | 3e-11 | 3.1e-11 | same value; reference error, < 1e-13 by `mp.diff` |
