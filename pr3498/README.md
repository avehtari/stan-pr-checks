# Checks for stan-dev/math#3498 (faster vectorized normal_lpdf)

This directory holds the scripts behind the numbers in the PR "Make
vectorized normal_lpdf faster without heap arrays" (stan-dev/math#3498,
branch `faster-normal-lpdf`).  Sections 1 to 7 quote the claims of the
PR text and give the commands and the expected output.  Section 8
answers a review question about `to_ref`.

Last run: 2026-11-01 and 2026-11-02.  All jobs ran on Triton (Aalto),
with g++ 14.3.0 (conda-forge) and `-O3`; the timings on Xeon Gold 6248
nodes with the `performance` governor, the OpenCL test on a Tesla V100.
In the logs, `<scratch>`, `<home>` and `<workstation>` replace the
machine paths.

## Commits

| name in this file | commit | what |
|---|---|---|
| develop | `dc66941969` | merge base of the PR (develop of 2026-10-29) |
| PR | `421386dec5` | PR head (8 commits) |
| 3a66 | `3a66fcf6aa` | first commit: the change of `normal_lpdf.hpp` |
| 620a | `620ab901f5` | fourth commit |
| 760f | `760f6a1f65` | sixth commit: PR head before the review |
| dde6 | `dde60be194` | seventh commit, after the review: `sigma` is checked first |

From `3a66fcf6aa` to `760f6a1f65`, the commits change `stan/` only in
comment lines, so the timings and the check (made with `3a66fcf6aa` and
`620ab901f5`) apply to `760f6a1f65`.  `dde60be194` removes the code that
kept the develop error order when `sigma` is invalid; `421386dec5`
changes a comment and adds tests.  The timings were not made again
after `dde60be194`.  Section 7 gives the check of `dde60be194`.

## Layout

| path | content |
|---|---|
| `make_inputs.sh` | makes `dev_normal_lpdf.hpp` and `br_<commit>.tar.gz` from a clone |
| `dev_normal_lpdf.hpp` | the develop `normal_lpdf`, renamed `normal_lpdf_dev` |
| `cases.hpp` | 25 type cases (argument kinds of `y`, `mu`, `sigma`) |
| `check_branch.cpp` | branch against develop: values, gradients, exception texts, heap allocations |
| `bench_branch.cpp` | time of one gradient evaluation; one binary per (version, case, propto) |
| `syntax_cases.cpp` | compile-only instantiation of many more argument types |
| `run_branch_bench.sh` | SLURM job: builds and runs the check, then 5 interleaved timing rounds |
| `summarize.py` | the ratio table from a job log |
| `mpmath_reference.py` | the constants of the mpmath test in `mix/prob/normal_test.cpp` |
| `to_ref_probe.cpp` | does develop's `to_ref` of a value view copy the values? (section 8) |
| `prototype/` | the cost measurement before the change (variants `f1` to `f4`) |
| `results/jobs/` | logs of the check and timing jobs |
| `results/tests/<tag>/` | logs of the unit-test, checklist and OpenCL jobs |

## Requirements

- A clone of stan-dev/math with the PR branch, and the `lib/` of a clone
  at `dc66941969` with a built TBB (`make -f make/standalone math-libs`).
- g++ with C++17; the job script loads the Triton modules
  `triton/2025.1-gcc scicomp-r-env`.  Elsewhere, remove the `module`
  line.
- Python 3 (tested with 3.10.12); mpmath (tested with 1.3.0) for
  `mpmath_reference.py`.

```sh
CLONE=<clone> ./make_inputs.sh dc66941969 3a66fcf6aa
STAN_MATH_LIB=<lib> TARBALL=br_3a66fcf6aa.tar.gz sbatch run_branch_bench.sh
python3 summarize.py nl_branch_<job>.log
```

`BENCH_CONFIGS=c:p+c:p+...` (index into `cases.hpp`, propto) replaces the
default list, cases A-H; `CHECK_ONLY=1` builds and runs only the check.

## 1. Speed

> In reverse mode with vector arguments, with this PR `normal_lpdf` takes
> 0.17 to 0.99 of the `develop` time in the cases that I measured (see
> Speed).

> In case A, `develop` takes 4.6 ns per element at N = 100 and 11 to 15
> ns from N = 1e4. The branch takes 1.9 to 4.0 ns. The largest ratio of a
> single round in the table is 0.89 (case C, N = 1e3).

> The largest ratio that I measured is for a data vector `sigma` with a
> scalar data `y` and an autodiff `mu` vector. It is 0.96 to 0.99 at N <=
> 1e3, and 0.62 to 0.72 from N = 1e4.

One evaluation creates the autodiff arguments, calls `normal_lpdf`, runs
`grad()` and `recover_memory()`.  Each binary times N = 1e2, 1e3, 1e4,
1e5; per N it takes the median of 7 batches of about 20 ms.  The job
runs 5 rounds, develop and branch alternating; a cell is the median over
the rounds.

| job | node | cases | log |
|---|---|---|---|
| 21346884 | csl38 | A-H (the PR table without Y) | `results/jobs/nl_branch_21346884.log` |
| 21347580 | csl7 | A, I, J, K, N, T, U, W | `results/jobs/nl_branch_21347580.log` |
| 21348124 | csl10 | T, Y | `results/jobs/nl_branch_21348124.log` |

```sh
python3 summarize.py results/jobs/nl_branch_21346884.log
```

Expected: the A-H rows of the PR table; case A develop 4.64 ns per
element at N = 100, 11.06 and 15.35 at 1e4 and 1e5; branch 2.37, 2.05,
1.89, 3.96.  The largest `max` in any A-H row is 0.89 (C, N = 1000).
Row Y of the PR table and case T (`y` data scalar, `mu` `Matrix<var>`,
`sigma` data vector) come from job 21348124: T 0.98, 0.99, 0.63, 0.72.
Job 21347580 gives T 0.96, 0.99, 0.62, 0.65.

## 2. The cost before the change (prototype)

> On `develop`, two parts of a vectorized call cost 56 to 76 % of its
> time in case A below (measured with a prototype):

> In the prototype, removing `y_scaled_sq` or `scaled_diff` saved 1.7 to
> 4.4 ns per element from N = 1e4.

`prototype/` holds the variants `f1` to `f4` of case A
(`normal_variants.hpp`), job 21345403 on `csl35`
(`prototype/nv_bench_21345403.log`).  Its job needs the tarball
`dev_dc66941969.tar.gz`:
`git -C <clone> archive --prefix=tree/ dc66941969 stan | gzip > prototype/dev_dc66941969.tar.gz`.

Median ns per element over the 5 rounds:

| N | dev | f1 | f2 | f3 | f4 |
|---|---|---|---|---|---|
| 100 | 5.00 | 3.33 | 2.65 | 2.42 | 2.08 |
| 1e3 | 4.31 | 2.88 | 2.24 | 2.07 | 1.91 |
| 1e4 | 8.72 | 7.48 | 4.26 | 2.19 | 2.06 |
| 1e5 | 10.83 | 9.99 | 8.32 | 3.92 | 3.43 |

- `f4` has neither the checks nor the heap arrays: (dev - f4) / dev is
  58 %, 56 %, 76 % and 68 %.
- `f1` to `f2` removes `y_scaled_sq`: 3.22 ns at 1e4, 1.67 at 1e5.
  `f2` to `f3` removes `scaled_diff`: 2.07 at 1e4, 4.40 at 1e5.

## 3. Behavior against develop

> 1512 of the 1692 comparisons give the same exception text, or values
> and gradients with an error below 1e-12 times max(1, |reference|). The
> other 180 comparisons have an invalid `sigma` together with an invalid
> `y` or `mu`; there the branch gives the `sigma` exception text. The
> largest error is 2.4e-15 for the values and 1.4e-15 for the gradients.
> A `malloc` counter finds no allocation of 8 N bytes or more on the
> branch in any of the 25 type cases; `develop` makes up to 7 such
> allocations per evaluation.

`check_branch.cpp` calls the branch `normal_lpdf` and `normal_lpdf_dev`
in one binary, N = 1000, for each type case of `cases.hpp` and 16 to 19
input cases (`make_inputs()` in the file), with `propto` true and false,
in reverse mode and with `fvar<double>` tangents (not for varmat).  A
`CHECK` line gives the outcome and two errors: `rel` scales by
max(1, |ref|, cond) (cond: the sum of the absolute terms for a scalar
argument), `rel_plain` scales the gradients by max(1, |ref|).  `HEAP`
lines count the `malloc` calls of at least 8 N bytes during one warm
gradient evaluation (interposer on `__libc_malloc`; `HEAP_SELFTEST`
must be 1).

Expected at `620ab901f5` (`results/jobs/nl_branch_21348124.log`):
`CHECK_SUMMARY ok 1692 fail 0`; largest `rel` 2.4e-15 (a value), largest `rel_plain` 1.4e-15;
`HEAP ... branch 0` in all 50 lines; largest develop count 7 (case I).
At `dde60be194` (`results/jobs/nl_branch_21361765.log`):
`CHECK_SUMMARY ok 1512 fail 180`, the same largest errors, see
section 7.

```sh
awk '$1=="CHECK" && $NF=="value" {if ($8+0>r) r=$8+0; if ($10+0>p) p=$10+0}
     END {print r, p}' results/jobs/nl_branch_21348124.log
```

## 4. Unit tests and the mutants

> The new tests also pass with the `develop` header, because each call
> in them has at most one invalid argument. To show that the tests cover
> the new code, I removed each of the two checks after the sum of
> squares (one per code path) in turn. At `760f6a1f65`, without the
> first, 8 assertions fail; without the second, 5 assertions fail.

The jobs of this section ran at `760f6a1f65`, CPU builds on batch
nodes.  At the PR head, `results/tests/nlu_t2/` (job 21362182) and
`results/tests/nlf_t3/` (job 21362690) give the same suites; see
section 5.  The mutants were not run again.

- `results/tests/nlf_a/` (job 21348417): `prim/prob/normal_test` (8
  tests), `mix/prob/normal_test` (5), `rev/prob/normal_log_test` (1),
  `rev/prob/normal_id_glm_lpdf_test` (26) pass with the PR header
  (stage `branch`) and with the develop header (stage `develop`).
- `results/tests/nlf_m1/` (job 21348418): mutant M1 removes the check
  after the sum in the branch where the partials of `mu` or `y` hold the
  scaled residuals.  Stage `branch` is the mutant: 8 assertions of
  `AgradRev.mathMixScalFun_normal_lpdf_messages` fail; prim passes.
  Stage `develop` is the PR header here: all pass.
- `results/tests/nlf_m2/` (job 21348419): mutant M2 removes the other
  check after the sum.  2 assertions of
  `ProbDistributionsNormal.lpdfAllDataChecks` and 3 of
  `AgradRev.mathMixScalFun_normal_lpdf_messages` fail.

```sh
grep -c "Failure" results/tests/nlf_m1/smtest_21348418_branch_mix_prob_normal_test.run
grep -c "Failure" results/tests/nlf_m2/smtest_21348419_branch_{prim,mix}_prob_normal_test.run
```

## 5. Other suites and the checklist

> These tests pass: `prim/prob/normal_test`, `mix/prob/normal_test`,
> `rev/prob/normal_log_test`, `rev/prob/normal_id_glm_lpdf_test` and the
> generated `test/prob/normal` lpdf tests.

> `opencl/rev/normal_lpdf_test` passes on a V100.

> cpplint, `make test-math-dependencies`, the header check, `make
> test-headers` and clang-format 10 pass. `make doxygen` with doxygen
> 1.9.7 gives no warnings.

> I ran the unit tests, the generated tests, the OpenCL test and the
> checklist items at `421386dec5`, with the change on `develop` at
> `dc66941969`, with g++ 14.3.0.

At the PR head.  Jobs 21362182 and 21362474 ran on the working tree on
top of `dde60be194`; its three changed files had the md5 of
`421386dec5`.

| check | commit | job | log |
|---|---|---|---|
| prim, mix `normal_test`, branch and develop header; `test/prob/normal` `-f normal_0`, 5346 checks | `421386dec5` | 21362182 | `results/tests/nlu_t2/` |
| `rev/prob/normal_log_test` (1), `normal_id_glm_lpdf_test` (26) | `421386dec5` | 21362690 | `results/tests/nlf_t3/` |
| `opencl/rev/normal_lpdf_test`, 6 tests, V100 | `421386dec5` | 21362691, 21362694 | `results/tests/nl_d3/` |
| cpplint, deps, header, test-headers | `421386dec5` | 21362474 | `results/tests/nlu_c2/` |
| doxygen, 0 warnings | `421386dec5` | 21362696 | `results/tests/nlf_c3/` |

Earlier runs:

| check | commit | job | log |
|---|---|---|---|
| `test/prob/normal`, `-f normal_0`: 7 programs, 5346 checks | `f83aee2c10` | 21346941 | `results/tests/nlpdf_t1/` |
| `opencl/rev/normal_lpdf_test`, 6 tests | `620ab901f5` | 21347595, 21347611 | `results/tests/nl_d2/` |
| cpplint, deps, header, test-headers, doxygen | `620ab901f5` | 21347651 | `results/tests/nl_c1/` |
| cpplint, header, doxygen | `760f6a1f65` | 21348529 | `results/tests/nlf_c/` |

The `test/prob` run of job 21346941 used the PR code (`stan/` as at the
head except comments).  In that job the mix suite did not compile: the
test file of `f83aee2c10` had an error, fixed in the next commit.
clang-format 10 has no log; `stan_clang_format.sh check` matched the
changed files before each commit.

## 6. The mpmath references of the mix test

`test/unit/math/mix/prob/normal_test.cpp` compares values and gradients
with constants from mpmath (40 digits) at the double inputs of the test,
with tolerance 1e-14 max(1, |ref|).

```sh
python3 mpmath_reference.py
```

Expected: `scalar_sigma value -6.3968830540223898`, `vector_sigma value
-5.3337600932202835`, and the `d_mu`, `d_sigma` lists of the test file.

## 7. Commit `dde60be194`: `normal_lpdf` checks `sigma` first

> Values and gradients agree with `develop` up to rounding. Exception
> types are the same as on `develop`. Exception texts are also the same,
> except when `sigma` is invalid together with `y` or `mu`: then the
> text is about `sigma`.

`check_positive(sigma)` now runs before the arithmetic, also when `y` or
`mu` is invalid.  The checks of `y` and `mu` stay where they were: after
the sum of squares, and before the size-zero and all-data `propto`
returns.  Only an input with an invalid `sigma` and another invalid
argument gives an exception text different from develop.

- `results/jobs/nl_branch_21361765.log` (job 21361765, `CHECK_ONLY=1`,
  the uncommitted working tree that became `dde60be194`): 1692
  comparisons, 1512 the same as develop.  The other 180 are the input
  cases `y_nan_sigma_neg` and `mu_inf_sigma_neg` (90 each): develop
  gives the `y` or `mu` message, the branch the `sigma` message.  Both
  messages come from `check_not_nan`, `check_finite` or
  `check_positive`, which all throw `std::domain_error`.
  `HEAP ... branch 0` in all 50 lines.
- The unit tests and the checklist at the PR head `421386dec5`: section 5.
- The mutant counts of section 4 (8 and 5) are from `760f6a1f65`; they
  were not run again.

```sh
grep "^CHECK .* FAIL" results/jobs/nl_branch_21361765.log | awk '{print $3}' | sort | uniq -c
```

## 8. `to_ref` of a value view (review question)

The body of `to_ref(T&& a)` is `return std::forward<T>(a);`, but its
return type is `ref_type_t<T&&>`.  `ref_type_t` is `plain_type_t<T>`
when `Eigen::Ref` cannot refer to the expression (`stan/math/prim/meta/ref_type.hpp`).  The
value view of a `Matrix<var>` reads the values through the `vari`
pointers, so `Eigen::Ref` cannot refer to it.  The return statement then
evaluates the view into a new array.  Develop's `normal_lpdf` calls
`to_ref(as_value_column_array_or_scalar(...))` for each argument (for
`y`: `y_ref`).

`to_ref_probe.cpp` (build line in the file; ran with g++ 11.4.0 on the
workstation; `to_ref.hpp`, `ref_type.hpp` and
`as_value_column_array_or_scalar.hpp` as at `dc66941969`), N = 1000:

```
Matrix<var>  to_ref returns a plain type: 1   allocations >= 8000 bytes: 1
VectorXd     to_ref returns a plain type: 0   allocations >= 8000 bytes: 0
```
