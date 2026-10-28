# Reproduction scripts for stan-dev/math#3409

This directory holds the scripts behind the numbers in the text of
stan-dev/math#3409 ("More stable incomplete beta gradient roots", part 1
of 3).  Each section below names one claim of the PR text, the commands
that measure it, the expected output, and the measured commits.

Last run: 2026-10-28, on a laptop (12 cores, Ubuntu 22.04, g++ 11.4.0,
Python 3.10.12, mpmath 1.3.0), except the speed job (Triton, see
§"Speed").

## Layout

| path | content |
|---|---|
| `build.sh` | compiles one probe against one Stan Math tree |
| `refs/` | mpmath references, two routes each; the test constants |
| `sweep/` | the 18 000-point accuracy sweep in nine regions |
| `points/` | the defect examples, the throws per function, the scalar types |
| `tfp/` | TensorFlow Probability at integer shapes |
| `ddz/` | `inc_beta_ddz` with the overflow policy against Boost's default |
| `speed/` | the speed benchmark and its SLURM job |
| `swap/` | `swap_test.sh`: the tests with the PR headers and with the base headers |
| `results/` | the outputs of the run of 2026-10-28 that the PR text quotes |

All commands below run from this directory (`pr3409/`).  New outputs go
to `out/`, which git does not track; compare them with `results/`.

## Requirements

- Python 3 with mpmath.  Tested with mpmath 1.3.0.
- g++ with C++17.  Tested with g++ 11.4 and 14.3.
- A clone of stan-dev/math with the PR branch:

  ```sh
  git clone git@github.com:stan-dev/math.git stan_math
  cd stan_math
  git fetch origin pull/3409/head:stable-inc-beta
  git checkout stable-inc-beta
  make -f make/standalone math-libs   # builds lib/tbb
  ```

- A develop tree.  The PR text calls the merge base `579253050c`
  "develop".  The probes need its `stan/` directory only:

  ```sh
  mkdir -p ../dev_579253050c
  git -C ../stan_math archive 579253050c stan | tar x -C ../dev_579253050c
  ```

Set these variables once:

```sh
export PR=$PWD/../stan_math              # the PR tree (stan/ and lib/)
export DEV=$PWD/../dev_579253050c        # stan/ of the merge base
export STAN_MATH_LIB=$PR/lib             # lib/ for both trees
export TBB_LIB=$PR/lib/tbb               # a directory with libtbb.so
mkdir -p out
```

## References for the unit tests

Claim: "The reference is mpmath at 80 digits, computed by two independent
routes."

```sh
python3 refs/gen_test_refs.py 8 > out/test_refs.txt
```

Route A is `mp.diff` of `I` at 110 digits; route B is quadrature of the
density at 80 digits (`refs/ref_mpmath.py`).  Both work on the side of
the symmetry relation where `I <= 1/2`.  Expected: every row shows `inf`
digits (the routes agree to the 60 digits of the comparison) at all 21
points and at the 3 lcdf gradients, and the values equal the constants
in the test files.  Run time about 30 s on 8 cores.

## Accuracy over 18 000 points

Claims: the accuracy table of the PR text, the develop counts, the worst
errors of the PR.

```sh
for t in dev:$DEV pr:$PR; do
  STAN_MATH=${t#*:} EXTRA_CXXFLAGS=-DPROBE_PRIM_ONLY \
    ./build.sh sweep/probe_cmp_beta.cpp out/cmp_dbl_${t%%:*}
  out/cmp_dbl_${t%%:*} sweep/points_beta.txt double > out/cmp_double_${t%%:*}.txt
done
python3 sweep/table_pr.py out/cmp_double_dev.txt out/cmp_double_pr.txt \
  sweep/refs_beta_18000.txt
```

`sweep/points_beta.txt` comes from `python3 sweep/sweep_beta.py gen 2000
sweep/points_beta.txt` (seed 20260920).  `sweep/refs_beta_18000.txt` holds
`I`, `dI/da` and `dI/db` at each point, from `ref_mpmath.py`; `python3
sweep/sweep_beta.py score ...` computes it (about 12 min on 20 cores).
The probe calls each root on its own, so a throw marks only its own
column.  Expected output, `579253050c` against `2764cf584d`:

```
region                         n |  A nf A >1e-2  B nf B >1e-2 | PR nf    worst     rel*  sub
R1 a,b<1                    2000 |     0       0     0       0 |     0  5.1e-13  1.2e-14    0
R2 1<=a,b<=20 z<0.75        2000 |     0       0     0       0 |     0  3.4e-14  1.9e-14    0
R3 seam z>0.75 a<500        2000 |     0      82     3     125 |     0  1.0e-13  4.3e-13    0
R4 seam b>a b>500 z<=0.75   1985 |     1      87   550     933 |     0  1.0e-12  1.8e-12    1
R5 b>>a z small             2000 |     0     408     2      12 |     0  1.9e-13  1.7e-12    0
R6 large a,b z uniform      1523 |    46      36   932     235 |     0  7.2e-13  3.8e-12   11
R7 large a,b z near mode    2000 |     0       0  1325      10 |     0  1.3e-13  3.3e-12    0
R8 negbin n<=200            1965 |     0       4     0      26 |     0  5.3e-13  2.4e-13    2
R9 negbin 200<n<=5000       1999 |    29      22   292    1203 |     0  4.2e-13  6.3e-13    0
total                            |    76     639  3104    2544 |     0
```

Path A is `inc_beta_dda` / `inc_beta_ddb`, path B the 9-argument
`grad_reg_inc_beta`.  A point counts once if either partial fails.  Run
time: develop 4 s, PR 0.1 s.

Claim: at the worst point of each region, the two reference routes agree
to 60 digits, so the code causes the error.

```sh
python3 sweep/check_worst.py out/cmp_double_pr.txt sweep/refs_beta_18000.txt
```

Expected (`results/check_worst.txt`): `inf` in the `routes` column in every
region, and the same relative error against the new reference as
against the cached one (3.78e-12 in R6, 3.28e-12 in R7, 1.84e-12 in R4).

## Defect examples

Claims: the values on develop at the defect points of the PR text, and
the throws.

```sh
for t in dev:$DEV pr:$PR; do
  for p in probe_points probe_throws; do
    STAN_MATH=${t#*:} ./build.sh points/$p.cpp out/${p}_${t%%:*}
  done
  out/probe_points_${t%%:*} > out/points_${t%%:*}.txt
  out/probe_throws_${t%%:*} sweep/points_beta.txt > out/throws_${t%%:*}.txt
done
```

Expected on develop (`results/points_dev.txt`): the 9-argument root is NaN
at `(600, 600, 0.5)` and 0.70535 at `(352, 590, 0.758)`;
`inc_beta_dda(1.145, 6786, 1.8e-4)` = -0.20271;
`inc_beta_dda(446, 5, 0.76)` = -1.78e-15; `neg_binomial_lcdf(132 | 313,
0.068)` d/dalpha = 0.  At the PR head (`results/points_pr.txt`) all of these
equal the references.  Both trees give NaN for d/dy of `beta_lcdf(0.5 |
600, 600)` and for d/dbeta of `neg_binomial_lcdf(132 | 313, 0.068)` (the
value code of the two lcdfs, not the root).

`results/throws_dev.txt` ends with `throws: dda 98 ddb 98 root9 24 beta 0`;
`results/throws_pr.txt` with all counts 0.

## Scalar types

Claim: the shape gradients of `inc_beta` in `var`, `fvar<double>` and
`fvar<fvar<double>>` against the `double` root.

```sh
for t in dev:$DEV pr:$PR; do
  STAN_MATH=${t#*:} ./build.sh points/probe_types.cpp out/probe_types_${t%%:*}
  out/probe_types_${t%%:*} sweep/points_beta_sub.txt > out/types_${t%%:*}.txt
done
```

Expected: on develop every cell is `0.0e+00/  0/ 0`.  At the PR head the
`var` columns and the two-shape `fvar` columns are 0; the one-shape `fvar`
and `fvar<fvar>` columns differ by up to 2.5e-13 relative (R5, `a` only),
because they call the one-sided root.

## TensorFlow Probability at integer shapes

Claim: "TFP's own d/db is off by 4.2e-4 at (2, 3, 0.25) and by 6.7 % at
(10, 1, 0.5).  Its d/da is off by 5.6 % at (1, 100, 0.01)."

```sh
python3 -m pip install --target /tmp/tfp_lib --no-deps \
  tensorflow-probability==0.25.0 absl-py decorator cloudpickle gast dm-tree
PYTHONPATH=/tmp/tfp_lib python3 tfp/tfp_betainc.py > out/tfp_betainc.txt
```

Expected: relative errors 4.2e-04, 6.7e-02 and 5.6e-02 at the integer
points, and below 5e-9 when both shapes move by 1e-7.

## `inc_beta_ddz` and Boost's `ibeta_derivative`

Claim: "bit-identical to the default call at 180 360 grid points, and
`+inf` at the 3 240 points where the default throws."

```sh
STAN_MATH=$PR ./build.sh ddz/check_ddz.cpp out/check_ddz
out/check_ddz
```

Expected: `points 180360, bit-identical 180360, default threw 3240` and
`inc_beta_ddz(1.05, 0.35, 1) = inf, (0.5, 2, 0) = inf`.

## Regression check

Claim: the new tests fail with the develop headers and pass with the PR.

```sh
T=test/unit/math
STAN_MATH=$PR BASE_REF=579253050c TBB_LIB=$TBB_LIB swap/swap_test.sh \
  $T/prim/fun/grad_reg_inc_beta_test $T/prim/fun/inc_beta_dda_test \
  $T/prim/fun/inc_beta_ddb_test $T/prim/fun/inc_beta_ddz_test \
  $T/rev/prob/beta_cdf_log_test $T/rev/prob/neg_binomial_cdf_log_test \
  $T/prim/fun/grad_reg_inc_beta_one_sided_test \
  $T/prim/fun/inc_beta_dda_vectorized_test \
  $T/prim/fun/inc_beta_ddz_vectorized_test
```

Expected in the base stage (`results/swap.txt`): in the first six files, 15
of 19 cases fail (6 of 8, 3 of 3, 3 of 3, 1 of 2, 1 of 1, 1 of 2).  The
four that pass are `grad_reg_inc_beta.2`, `grad_reg_inc_beta.moderate`,
`inc_beta_ddz` and `NegBinomial_derivatives_lcdf`, at points where
develop is right.  The last three files do not compile, because they
test the new overloads.  All 25 cases pass in the stages `pr_1` and
`pr_2`.  Triton jobs 21268803 and 21270239 (g++ 14.3) gave the same
result before the last comment edits.

## Speed

Claims: the function-level and distribution-level speed tables.

```sh
WT=$PR speed/make_trees.sh            # packs speed/trees.tar.gz
sbatch --export=ALL,TREES=trees.tar.gz,STAN_MATH_LIB=<lib>,MODULES=<m1:m2> \
  speed/run_cmp_beta_speed.sh         # run from the directory of the files
```

The job builds `speed/probe_cmp_beta_speed.cpp` against both trees with
`-O3` and runs each twice on one exclusive node.  The run of the PR text
is `results/pr3409_speed_21271376.log`: Triton job 21271376, Xeon Gold 6230
at 2.10 GHz, `performance` governor, conda-forge g++ 14.3.0, trees
`579253050c` and `2764cf584d` with the working tree of this revision.
The two rounds agree within 1 %.  The tables of the PR text give the mean
of the two rounds of net ns.

`speed/probe_speed3.cpp`, `run_pr3409_steve_roihu.sh` and
`table_speed.py` are the benchmark behind the speed comment of
2026-10-22 on the PR (original port against the reverse-mode code of
#3489).  They compare three trees: develop `e4745c87ff`, the original
port `bd8f27a265`, and `d2eeacd676`.  Roihu jobs 2443389 and 2443471,
AMD EPYC 9965, g++ 11.5, `-O3`.  I did not run them again.
