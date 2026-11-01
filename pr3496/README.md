# Checks for stan-dev/math#3496 (discrete cdfs, part of issue #56)

This directory holds the scripts behind the numbers in the PR "Fix
discrete cdfs when one element of a vector has probability 0 or 1"
(stan-dev/math#3496, part of issue #56).  Until 2026-10-31 this directory
was `discrete_cdfs_edge_elements/`; the PR body links that name at commit
`e7fceb3`.  Sections 1 to 8 quote the claims of the PR text and give the
commands and the expected output.

An edge element is an element with probability 0 or 1.  A swap test
runs the new test with the develop versions of the changed headers.

Last run: 2026-10-31.  The probes ran on a laptop (Ubuntu 22.04, g++
11.4.0, `-O1`, Python 3.10.12, mpmath 1.3.0).  The unit tests, the swap
test, the generated tests and the checklist ran on Triton batch nodes
(g++ 14.3.0), and the OpenCL tests on a V100.  `results/versions.txt`
records the versions of the laptop run.  In the Triton logs, the
placeholders `<workstation>`, `<scratch>` and `<home>` replace the
machine paths.

## Commits

| name in this file | commit | what |
|---|---|---|
| develop | `dc66941969` | merge base of the PR (develop of 2026-10-29) |
| PR | `1315143b74` | PR head, tree `1741248cd7` |
| #3433 | `709a1a4d76` | head of stan-dev/math#3433, stacked on #3490, #3489 and #3409; for section 4 only |

Most Triton jobs in `results/` ran before the squash of three follow-up
commits, at the commits that their logs name.  The last ones ran at
`c33769a28f` (jobs 21321024, 21321025, 21321126) and at `1315143b74`
(job 21321997).  Both commits have the tree `1741248cd7`.

## Layout

| path | content |
|---|---|
| `build.sh` | compiles one probe against one Stan Math tree |
| `make_trees.sh` | exports `stan/` and `test/unit/` of develop, the PR and #3433 into `out/trees/` |
| `run_probes.sh` | runs sections 1 to 5 and writes their outputs to `out/` |
| `edges/` | the probe behind the table of the PR text; the UBSan probe |
| `refs/` | mpmath references of the quoted numbers, two routes each |
| `large_n/` | `neg_binomial_2_lcdf` at large finite counts |
| `overlap/` | the headers that the open PRs also change |
| `swap/` | the swap test for a local clone |
| `results/` | the outputs of 2026-10-31; `results/triton*/` are the Triton job logs |

Run all commands from this directory.  New outputs go to `out/`, which
git does not track; compare them with `results/`.

## Requirements

- Python 3 with mpmath (tested with 1.3.0).
- g++ with C++17 (tested with 11.4 and 14.3).
- A clone of stan-dev/math with the PR branch and a built TBB
  (`make -f make/standalone math-libs`).  `make_trees.sh` fetches the
  head of #3433 if the clone does not have it.

```sh
CLONE=<clone> ./make_trees.sh
CLONE=<clone> STAN_MATH_LIB=<clone>/lib TBB_LIB=<clone>/lib/tbb ./run_probes.sh
```

## 1. The table: edge-element defects on develop

> The functions in the table below had one or both of two defects. One
> element with probability 1 made the whole result 1 (`cdf`) or 0
> (`lcdf`, `lccdf`), and the function ignored the other elements. ...
> An element with probability 0 gave the correct `-inf`. The gradient
> was not 0 ...

`edges/probe_edges.cpp` evaluates each cdf, lcdf and lccdf of the 8
families with one ordinary and one edge element.  It uses both orders,
reverse mode and scalar parameters.  An edge element with probability 1
must give the value and the gradient of the ordinary element alone.  An
edge element with probability 0 must give 0 or `-inf` and a zero
gradient.  The tolerance is 1e-12 relative.

Expected (`results/edges_develop.txt`, `results/edges_pr.txt`):

- The output has 92 rows.  On the PR, 0 rows are `BAD`.
- On develop, 39 rows are `BAD`.  Each one belongs to a row of the PR
  table.
- With an edge element of probability 0, all functions except
  `poisson_binomial_lccdf` are `BAD` only in the order
  `{ordinary, edge}`: there the function had already set the partials of
  the ordinary element.
- `poisson_binomial_lccdf` is `BAD` in both orders.  With 3 trials per
  element the value is finite (-37.5353) instead of `-inf`.  With 1 trial per element
  (`poisson_binomial_lccdf_1`) it is `-inf` with a nonzero gradient.
- `poisson_binomial_cdf` and `_lcdf` at `y = N` are not `BAD` on
  develop.  They differ from the ordinary element alone by 1.1e-16 in
  the value and the gradient.  The PR text names them below the table.

## 1b. The functions that already followed the rule

The rule: a first pass returns 0 or `-inf` before it sets any partial,
and the main loop skips each element with probability 1.

> The PR does not change the discrete cdf functions that already
> followed the rule: `bernoulli_cdf` and `_lcdf`, `binomial_cdf` and
> `_lcdf`, `beta_binomial_cdf` and `_lcdf`, and `discrete_range`.

In `results/edges_develop.txt`, all rows of `bernoulli_cdf`,
`bernoulli_lcdf`, `binomial_cdf`, `binomial_lcdf`, `beta_binomial_cdf`
and `beta_binomial_lcdf` are `OK`.  The `discrete_range` functions have
only integer arguments, so they have no gradient.  The
`CLAIM discrete_range` line shows that
`discrete_range_lccdf({2, 0} | 1, 5)` and `discrete_range_lcdf({2, 9} |
1, 5)` equal the value of the element 2 alone (-0.510826 and -0.916291).

#3408 uses the same rule in the `poisson` cdf functions.  In
the head of #3408, `git show pr/3408:stan/math/prim/prob/poisson_lccdf.hpp`
shows the first pass (lines 63 to 72) and the `continue` for `n < 0`
(lines 76 to 78).

## 2. The quoted numbers

> Example: `neg_binomial_lccdf({3, -1} | 2.5, 1.5)` returned 0, but the
> correct value is `neg_binomial_lccdf(3 | 2.5, 1.5) = -2.04`.

> For `neg_binomial_2_lcdf({3, INT_MAX} | 4, 2)` on `develop`, the `phi`
> derivative is 19.41; the correct value, the `phi` derivative of
> `neg_binomial_2_lcdf(3 | 4, 2)`, is -0.0604.

> For `theta = (0.5, 0.2, 0.7)` and `y = 3` it returned -36.74, which is
> `log(2^-53)`.

The `CLAIM` lines of `results/edges_develop.txt` and
`results/edges_pr.txt` give the Stan values.

| quantity | develop | PR |
|---|---|---|
| `neg_binomial_lccdf({3, -1} \| 2.5, 1.5)` | 0 | -2.04043 |
| `d/dphi neg_binomial_2_lcdf({3, INT_MAX} \| 4, 2)` | 19.4057 | -0.0604443 |
| `poisson_binomial_lccdf(3 \| 0.5, 0.2, 0.7)` | -36.736800569677101, equal to `log(2^-53)` | `-inf` |

`refs/ref_claims.py` computes the references (`results/ref_claims.txt`).
The two routes of each number agree to 5e-41 or better:

- `neg_binomial_lccdf(3 | 2.5, 1.5)` = -2.0404283800026 (pmf sum; upper
  incomplete beta);
- `neg_binomial_2_lcdf(3 | 4, 2)` = -0.617864120139397, `d/dmu`
  -0.203562340966921, `d/dphi` -0.060444349736812 (pmf sum; lower
  incomplete beta);
- `poisson_binomial(0.5, 0.2, 0.7)`: `P(Y = 0..3)` = 0.12, 0.43, 0.38,
  0.07 (recursion over the trials; enumeration), so `P(Y <= 3) = 1`.
  With 3 trials, `P(Y > 3) = 0`, and the lccdf at `y = 3` is `-inf`.

> The existing test `prim/prob/poisson_binomial_ccdf_log_test` expected
> that value. This PR changes the expectation to `-inf`, as the other
> `y = N` cases in that file expect.

`test/unit/math/prim/prob/poisson_binomial_ccdf_log_test.cpp` expects
`-inf` at `y = N` in lines 11 to 14 and 25 to 28 (`N = 0` and `N = 1`).
The PR changes the expectation in line 41 (`N = 3`) to `-inf`.

## 3. The int overflow in neg_binomial_lcdf

> With `n = INT_MAX`, `n + 1` overflows `int` (UBSan reports it on
> `develop`); this PR uses `n + 1.0`.

`edges/ubsan_nb_lcdf.cpp`, built with `-fsanitize=undefined`
(`results/ubsan_nb_lcdf.txt`): on develop, `runtime error: signed
integer overflow: 2147483647 + 1` at `neg_binomial_lcdf.hpp:77`.  On the
PR, `neg_binomial_lcdf({3, INT_MAX} | 2.5, 1.5) = -0.139231`.

## 4. Large finite counts in neg_binomial_2_lcdf

> Large finite counts still go through `beta_lcdf`, and this PR does not
> change that path. There the `phi` derivative is also wrong. At
> `mu = 4`, `phi = 2`, it is 11.8 for `n = 10^6` and 18.7 for
> `n = 10^9`. The correct derivative is below `1e-176000`. At the head of
> #3433, which is stacked on #3409, #3489 and #3490, the derivative
> becomes 0.

`neg_binomial_2_lcdf(n | mu, phi)` is `beta_lcdf(phi / (mu + phi) | phi,
n + 1)`.  The PR removes only the `n = INT_MAX` elements before the call.
`large_n/probe_nb2_lcdf_large_n.cpp` evaluates the value and both
derivatives at 9 points.  `large_n/score_large_n.py` compares them with
two mpmath routes.  Route A is the log of the lower incomplete beta.
Where `n > 100 mu`, route A is `log1p(-S)`, with `S` the direct sum of
the pmf over `k > n`.  Route B is `log1p(-I)`, with `I` the upper
incomplete beta.  For the values at `phi = 2`, the routes are
independent in all rows.  They agree to 4e-51 or better.

Expected output: `results/nb2_large_n_score_{develop,pr,pr3433}.txt`.
The cells give relative errors.  Where the true value is near 0, a cell
gives the computed value.
The develop and PR outputs are identical: the PR changes only the path
for `n = INT_MAX`.

| n | mu | PR (= develop) | head of #3433 (with #3409, #3489, #3490) |
|---|---|---|---|
| 10 | 4 | all within 2.4e-14 | within 2.9e-15 |
| 1e3 | 4 | `d/dphi` -1.9e-13; true 2.9e-172 | correct |
| 1e6, 1e9, 2^31 - 2 | 4 | `d/dphi` 11.8, 18.7, 19.5; true below 1e-176000 | 0 |
| 1e3, 1e6 | = n | within 3.1e-9 | within 8.1e-9 |
| 1e9, 2^31 - 2 | = n | value 6.9e-8, `d/dmu` 3.0e-6, `d/dphi` 1.4e-5 | `d/dphi` 3.5e-5, the others the same |

At `n = mu = 1e9`, the head of #3433 does not fix the error.  I did not
check which PR of the stack makes the tail derivative 0.

## 5. The merge-order note

> Six of them change 12 of the 18 headers in this PR. #3490 changes the
> four `neg_binomial` cdf headers, and #3433 changes `binomial_lccdf` and
> `neg_binomial_lccdf`. #3434 changes `beta_binomial_lccdf` and the
> `beta_neg_binomial` and `yule_simon` cdf headers. #3437 and #3445 also
> change `beta_binomial_lccdf`, and #3441 also changes the three
> `beta_neg_binomial` headers. #3408, #3409 and #3489 change none of these
> headers. #3490 and #3433 are stacked on #3409 and #3489.

Several of these PRs are stacked, so `overlap/overlap.sh` compares each
PR with its own base branch (`results/overlap.txt`):

| PR | base branch | headers of this PR that it changes |
|---|---|---|
| #3408 | `develop` | 0 |
| #3409 | `develop` | 0 |
| #3489 | `stable-inc-beta` (#3409) | 0 |
| #3490 | `stable-inc-beta-rev` (#3489) | 4: `neg_binomial_2_cdf`, `neg_binomial_cdf`, `neg_binomial_lccdf`, `neg_binomial_lcdf` |
| #3433 | `stable-inc-beta-cdfs` (#3490) | 2: `binomial_lccdf`, `neg_binomial_lccdf` |
| #3434 | `develop` | 7: `beta_binomial_lccdf`, the three `beta_neg_binomial` and the three `yule_simon` headers |
| #3437 | `stable-beta-binomial-lpmf` (#3434) | 1: `beta_binomial_lccdf` |
| #3441 | `stable-beta-binomial-lpmf` (#3434) | 3: the three `beta_neg_binomial` headers |
| #3445 | `fix-hypergeometric-3F2-termination` | 1: `beta_binomial_lccdf` |

12 headers in total.  The same lists come from the GitHub file list of
each PR (`gh api repos/stan-dev/math/pulls/<N>/files`).

## 6. Swap test

> All 11 test cases pass with the change and fail with the `develop`
> headers (503 failed expectations). For this test, I replaced all 17
> changed `prim` headers with their `develop` versions at once.

Triton job 21321024 (node csl17, g++ 14.3.0) ran the new test three
times: with the PR headers, with the develop versions of the 17 changed
`prim` headers, and with the PR headers again.  Logs:
`results/triton_21321024/`.

| stage | cases | failed |
|---|---|---|
| PR | 11 | 0 |
| develop headers | 11 | 11 |
| PR again | 11 | 0 |

Failed expectations per test case with the develop headers (count of
`: Failure` lines in the develop `.run` log): `bernoulli` 8, `binomial`
10, `beta_binomial` 16, `neg_binomial` 36, `neg_binomial_2` 32,
`beta_neg_binomial` 54, `yule_simon` 24, `poisson_binomial` 15,
`vector_parameters` 184, `poisson_binomial_theta_per_element` 80,
`neg_binomial_2_lcdf_expect_ad` 44; 503 in total.

To run the swap test in a local clone at the PR head:

```sh
TREE=<clone at the PR head> ./swap/swap.sh
```

The script stops if the clone has changes in `stan/` or a tree other
than `1741248cd7`.  It restores the PR headers on exit, and it exits
with 0 only if the three stages give pass, fail, pass.

## 7. OpenCL

> The existing OpenCL test `opencl_matches_cpu_small_n_negative` of
> `bernoulli_lccdf` ... passes on a V100. With the `develop` OpenCL file
> and the fixed CPU file, it fails. Two new OpenCL tests ... compare
> OpenCL with OpenCL and also pass on `develop`.

Triton jobs 21321025 (build) and 21321044 (V100):
`test/unit/math/opencl/rev/bernoulli_lccdf_test` passes 9 of 9 cases.
With the develop `stan/math/opencl/prim/bernoulli_lccdf.hpp`,
`opencl_matches_cpu_small_n_negative` fails, and the two new tests
`opencl_broadcast_n_negative` and `opencl_broadcast_theta_n_negative`
pass.  Logs: `results/triton_ocl_21321025/`.  These OpenCL builds also
compile the body of the OpenCL header, which the CPU header check skips
(`#ifdef STAN_OPENCL`).

## 8. Other test runs and the checklist

> The new test and 35 existing unit suites of the 8 families (`prim`,
> `rev`, `mix`) pass. The generated `test/prob` tests of the 7 families
> that have them (all except `poisson_binomial`) pass.

Logs: `results/triton_jobs/`.  Jobs 21319265, 21320438, 21320579,
21320777 and 21320800 ran these 35 existing suites with
`make <suite>` in a tree of develop plus the PR diff:

- `mix/prob`: `bernoulli_cdf`, `bernoulli_lcdf`, `bernoulli_lccdf`,
  `neg_binomial`;
- `prim/prob`: `bernoulli`, `beta_binomial`, `beta_neg_binomial`,
  `binomial`, `binomial_lpmf`, `neg_binomial_2`, `neg_binomial`,
  `yule_simon`, `poisson_binomial`, and the `cdf_log` and `ccdf_log`
  tests of `bernoulli`, `beta_binomial`, `binomial`, `neg_binomial_2`,
  `neg_binomial` and `poisson_binomial`;
- `prim/fun`: `poisson_binomial_log_probs`;
- `rev/prob`: `binomial_lpmf`, `neg_binomial_2_cdf`, `neg_binomial_2`,
  `neg_binomial_cdf`, `neg_binomial`, and the `cdf_log` and `ccdf_log`
  tests of `neg_binomial_2` and `neg_binomial`.

The `FAIL` and `BUILD_FAILED` lines in these logs come from commits
before the squash:

- the new test at `9518476257` (job 21319265, before the
  `neg_binomial_2_lcdf` commit);
- `poisson_binomial_ccdf_log_test` at `56872c02d1` (job 21320579, before
  the changed test expectation);
- one link failure of `prim/prob/neg_binomial_test` in job 21320579.
  Its log does not show the message, and job 21320777 built and passed
  the same suite in the same tree.

Every suite passed in a later job.

Generated tests: `./runTests.py -j6 test/prob/<family>` for `bernoulli`,
`binomial`, `beta_binomial`, `neg_binomial`, `neg_binomial_2`,
`beta_neg_binomial` and `yule_simon`, every program of each family (job
21319265).  After job 21319265, the PR commits changed only
`neg_binomial`, `neg_binomial_2` and `poisson_binomial` code, and
comments in the CPU and OpenCL `bernoulli_lccdf`.  Jobs 21320438 and 21320579 ran
`test/prob/neg_binomial_2` and `test/prob/neg_binomial` again.

Checklist at `1315143b74` (job 21321997): `make test-math-dependencies`,
a compile check of each of the 18 changed headers, `make test-headers`,
and `make doxygen` pass.  `results/doxygen_setup.txt` gives the doxygen
version (1.9.7) and the `WARN_AS_ERROR = FAIL_ON_WARNINGS` line of
`doxygen/doxygen.cfg`.  cpplint passed at `c33769a28f`, tree
`1741248cd7` (job 21321126).  `results/clang_format.txt`: clang-format
10 reports no changes in the 21 changed files at `1315143b74`.
