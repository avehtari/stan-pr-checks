# beta(a, b) without cancellation at large arguments

This directory holds the scripts behind the numbers of stan-dev/math#3497,
"Compute beta(a, b) without cancellation at large arguments" (branch
`stable-beta-function`).  The PR text links this directory under its
first name, `math-stable-beta-function/`, at commit `23abe61`.  Sections 1 to 8 quote the claims of the PR text
and give the commands, the expected output and the commits.

## Commits and trees

| name | commit | what |
|---|---|---|
| develop | `dc66941969` | stan-dev/math develop of 2026-10-29, the base of the branch |
| branch | `b40012df74` | head of `stable-beta-function` at the first review (3 commits) |
| head | `47159a62ca` | head after the review of 2026-11-03 (6 commits); `std::isfinite` and OpenCL `isfinite` replace `< INFTY` and `< INFINITY` (same results for all arguments), `rev/prob/beta_lcdf_test` tolerance 1e-7 |
| a64 | `a64c3679fa` | second commit of the branch; its `stan/` differs from `b40012df74` only in a comment |
| pr3433 | `709a1a4d76` | head of stan-dev/math#3433 (stacked on #3490, #3489 and #3409); section 6 only |

The Triton jobs of sections 3, 5 and 7 ran at `a64c3679fa`.  The unit
tests of section 8 ran at both commits (table in section 8).

```sh
CLONE=<stan-dev/math clone> ./make_trees.sh develop:dc66941969 \
    branch:b40012df74 pr3433:709a1a4d76
STAN_MATH_LIB=<clone>/lib TBB_LIB=<dir with libtbb.so> ./run_probes.sh
```

`run_probes.sh` writes to `out/` the files whose expected versions are in
`results/local/` and `results/background/`.  It builds the local probes
with g++ 11.4.0 and `-O1`, and it scores them with Python 3.10.12 and
mpmath 1.3.0 (`results/local/versions.txt`).  The Triton jobs used g++
14.3.0, `-O3` and mpmath 1.3.0.  Their logs are in `results/jobs/` and
`results/tests/`, with the machine paths replaced by `<workstation>`,
`<scratch>` and `<home>`.

## 1. The defect and the table of the summary

> `stan::math::beta(a, b)` computes `exp(lgamma(a) + lgamma(b) -
> lgamma(a + b))`. For large `b`, `lgamma(b)` and `lgamma(a + b)` are both
> about `b log(b)`, so the relative error of `beta(a, b)` is about
> `b log(b) eps`. ...
>
> Relative error of `beta(a, b)` against mpmath (probes built with g++
> 11.4.0, `-O1`): (table)

`probe_table.cpp` prints `beta(a, b)` at the six points of the table.
`refs/score_table.py` scores the develop and branch outputs with
`reference()` of `refs/gen_beta_reference.py`.  `reference()` uses two
mpmath routes, and the closed form `1 / (b (b + 1))` for `a = 2`.  Its
mpmath precision grows with the arguments.

Expected: `results/local/score_table.txt`.

## 2. Edge inputs

> On `develop`, `lgamma(1.7e308)` overflows, and the formula gives
> `inf - inf`. My edge-case probe has 21 pairs with one argument 1.7e308
> and a finite positive other argument, and `develop` gives NaN for all
> of them. With this PR, these values are finite, or 0 or `inf` where the
> true value underflows or overflows. No other result at 0, `inf`,
> negative or NaN arguments changes.

`probe_edges.cpp` evaluates `beta(a, b)` at the 289 pairs of 17 inputs:
NaN, 0, -0.5, -1, -2.5, `inf`, 4.9e-324, 1e-300, 1e-3, 0.5, 1, 2, 9.99,
10, 1e3, 1e9, 1.7e308.  For each pair, it also prints `exp(lbeta(a, b))`
and the value that the rule of the branch selects (the rule column).  `refs/compare_edges.py`
compares develop and branch.

Expected in `results/local/edges_compare.txt`:

- 229 pairs identical;
- 21 class changes (the classes are NaN, `inf`, 0 and finite), all NaN on
  develop, all with one argument 1.7e308;
- 39 other finite value changes, all with at least one argument at least
  10;
- the branch equals its rule column at all 289 pairs.

`results/local/edges_develop.txt` and `edges_branch.txt` hold the values.

## 3. Bitwise identity below 10, and the accuracy sweep

> When both arguments are below 10, `lbeta` returns `lgamma(x) +
> lgamma(y) - lgamma(x + y)` with `x` the smaller argument. ... The
> accuracy sweep below confirms this at all 25 600 points of the region.

> `a` and `b` each take 495 values from 1e-3 to 2^31 (40 per decade):
> 245 025 points. ... (table) ... At 136 points, the error of this PR is
> above 4 times the `develop` error and above 1e-15. At these points the
> error of this PR is at most 1.4e-13, and at most 1.45 times
> `(1 + |log B|) eps`, with `B = beta(a, b)`. ...

`sweep/grid.hpp` defines the grid.  `sweep/probe_sweep.cpp` prints
`beta(a, b)` on the grid and on the edge inputs of section 2.
`sweep/score_sweep.py` scores develop and branch against two mpmath
routes at 50 digits.  The job script `sweep/run_sweep.sh` builds the probe
against both trees of `out/beta_trees.tar.gz`.  To make this file, run
`make_trees.sh develop:dc66941969 branch:a64c3679fa`, then `tar -C
out/trees -czf out/beta_trees.tar.gz develop branch`.  In the jobs of
sections 3, 5 and 7, the tree `branch` holds `a64c3679fa`, not
`b40012df74`.

```sh
triton_deploy_submit.sh submit --dir beta_fn_checks2 --job sweep/run_sweep.sh \
    --export LIB=<Stan Math lib/ with a built tbb/> sweep/probe_sweep.cpp \
    sweep/grid.hpp sweep/score_sweep.py out/beta_trees.tar.gz
```

Job 21346661 (Triton `pht91`, g++ 14.3.0, `-O3`, mpmath 1.3.0, 32
processes).  Expected: `results/jobs/score_sweep.txt` and
`results/jobs/beta_sweep_21346661.log`.  In `score_sweep.txt`:

- the `small`, `mixed` and `large` rows (both below 10, one below 10,
  both at least 10) give the table of the PR text;
- `bitwise_diff` counts the points that differ from develop (0 in
  `small`); the PR column "bitwise equal" is `points - bitwise_diff`;
- the `underflow` row has the 88 186 points with a reference below
  2.2e-308 (errors in units of 4.9e-324);
- `nonfinite_d/b` is 0 in every row;
- the `worse` column sums to 136;
- the `WORSE` lines list these points: develop errors are at most
  2.7e-14, branch errors are at most 1.4e-13, and `branch_scaled` is at
  most 1.45;
- the two routes agree to 1.77e-40.

## 4. The callers and the beta_lcdf x derivative

> Example: I evaluated the `x` derivative of `beta_lcdf(x | 2, b)` at
> `b x` in {0.5, 1, 2, 8}. On `develop`, its relative error is up to
> 2.0e-6 at `b = 2^30` and 3.0e-6 at `b = 2^31`. With this PR, it is at
> most 2.4e-15.

`probe_caller.cpp` evaluates `beta_lcdf(x | 2, b)` and its `x` derivative
in `var` for `b` in {2^20, 2^30, 2^31} and `b x` in {0.5, 1, 2, 8}.  At
these points `x`, `1 - x` and `b x` are exact.  `refs/score_caller.py`
scores them against mpmath at 50 digits (`mp.betainc`, and the binomial
sum for `a = 2`; the routes agree to 5.6e-51).

Expected: `results/local/score_caller_develop.txt` (`d_dx` 1.82e-10 at
`2^20`, 2.02e-6 at `2^30`, 2.98e-6 at `2^31`) and
`score_caller_branch.txt` (worst 2.38e-15).  The list of callers in the
PR text comes from `git grep` of `beta(` in `stan/` at `dc66941969`.

## 5. OpenCL

> OpenCL (V100): the probe used the 18 reference rows of the new tests
> that have no argument 1.7e308. At all 18 rows, the device result equals
> the CPU result. At 1 000 grid points per region, the worst device
> relative errors are 5.9e-15, 3.4e-14 and 1.7e-13 (both below 10, one
> below 10, both at least 10). At the same points, the worst CPU relative
> errors are 5.9e-15, 3.4e-14 and 1.5e-13. In the region with both
> arguments at least 10, 793 of the 1 000 points underflow; the worst
> errors above are for the other 207 points.

`sweep/probe_ocl.cpp` computes `beta` on the device (`matrix_cl`) and on
the CPU.  It uses the 18 rows of the `a64c3679fa` version of
`beta_reference.hpp` (without the two 1.7e308 rows) and 1 000 grid points
of each region.  `sweep/score_ocl.py` gives the ulp distance between
device and CPU and the relative error of both against mpmath.  The job
script is `sweep/run_ocl.sh` (Triton partition `gpu-debug`).

Job 21346663 (Triton `dgx10`, Tesla V100-SXM2-16GB).  Expected:
`results/jobs/score_ocl.txt` and `results/jobs/beta_ocl_21346663.log`.
The device and the CPU differ by up to 31, 238 and 939 ulp in the three
regions; both have the relative errors above.

## 6. Out of scope

> The `b` derivative of `beta_lcdf` at the head of #3433 still has a
> relative error of 7.3e-6 at `a = 0.5`, `b = 1e9`. ...
> Boost `ibeta(2, 1e9, 2e-9)` has a relative value error of 3.7e-8. ...
> When `1 - x` is rounded, this term has a relative error of 5.5e-8 at
> `b = 1e9`. ...

These numbers come from the probes of an earlier issue draft (not
posted), in `background/`:

- `background/probe_beta_lcdf.cpp` and `background/score_beta_lcdf.py`:
  `beta_lcdf` and its three derivatives, and Boost `ibeta`, at 63 points
  (`a` in {0.5, 2, 10}, `b` from 1e3 to 2^31, `b x` in {0.5, 2, 8}),
  against two mpmath routes.  Expected:
  `results/background/score_develop.txt` and `score_pr3433.txt`.
  `score_pr3433.txt` is built against `out/trees/pr3433`
  (`709a1a4d76`).  Its worst `d/db` is 7.3e-6 at `a = 0.5`, `b = 1e9`,
  `b x = 0.5`.  Boost `ibeta` has its worst error, 3.7e-8, at `a = 2`,
  `b = 1e9`, `b x = 2`.
- `background/probe_dx_parts.cpp`: the parts of the `x` derivative on
  develop.  Expected: `results/background/dx_parts_score.txt` (`pow(1 - x,
  b - 1)`: 5.5e-8 at `b = 1e9`).  The scoring command:

  ```sh
  STAN_MATH=out/trees/develop ./build.sh background/probe_dx_parts.cpp out/probe_dx_parts
  ./out/probe_dx_parts > out/dx_parts.txt
  python3 -c "
  import mpmath as mp
  mp.mp.dps = 40
  for l in open('out/dx_parts.txt'):
      a, b, x, be, pw, el = [mp.mpf(v) for v in l.split()]
      rb = mp.beta(a, b)
      rp = (1 - x) ** (b - 1)
      print(mp.nstr(a, 3), mp.nstr(b, 4), mp.nstr(abs(be - rb) / rb, 3),
            mp.nstr(abs(pw - rp) / rp, 3), mp.nstr(abs(el - rb) / rb, 3))
  "
  ```

- `background/probe_beta_lccdf_dx.cpp` and `score_lccdf_dx.py`: the `x`
  derivative of `beta_lccdf` on develop and at #3433
  (`results/background/score_lccdf_dx_*.txt`).

`run_probes.sh` reproduces the `results/background/` files of the first
and third probes byte for byte.  It does not run the second probe; the
command above does.  The `yule_simon` item of the PR text
quotes `stan/math/prim/prob/yule_simon_{cdf,lcdf,lccdf}.hpp` at
`dc66941969` (lines 73, 73, 72).

## 7. Speed

> I timed `beta(double, double)` on 4 096 log-uniform input pairs per
> region, and `beta_lcdf` with `var` arguments at `b = 1e9` (Triton job
> 21346662). Each function has its own binary (`-O3`), and the job runs 3
> interleaved rounds on one core of a Xeon Gold 6248. ... (table)

`speed/bench_beta.cpp` times `beta` on 4 096 log-uniform input pairs per
region (`-DBENCH_FN=1`).  It also times `beta_lcdf` with `var` arguments
and its gradient at `b = 1e9` (`-DBENCH_FN=2`).  Each time is the minimum
over 15 blocks.  The `bump = acc * 1e-300` line stops the compiler from
hoisting the loop body, and it makes each call depend on the previous
result.  `speed/run_bench.sh` builds one binary per function and per
tree, and it runs 3 interleaved rounds pinned to one core (`taskset`).

Job 21346662 (Triton `csl3`, Intel Xeon Gold 6248 at 2.50 GHz, governor
`performance`, g++ 14.3.0).  Expected: `results/jobs/beta_bench_21346662.log`
(`TIME round<r> <tree> <function> <ns per call>`).  The ratios of the PR
text are branch / develop of the same round.  The develop column is the
middle value of the 3 rounds, rounded.

## 8. Tests

> With the `develop` `prim/fun/beta.hpp`, 4 of the 5 new CPU test cases
> fail (28 failed expectations). ... With the `develop` OpenCL device
> function and `elt_function_cl.hpp`, all 3 new OpenCL test cases fail.
> ... 33 existing unit suites ... (185 programs). ...

Reference rows: `refs/gen_beta_reference.py` prints
`test/unit/math/prim/fun/beta_reference.hpp`.  Its routes agree to
6.2e-67 (`results/local/beta_reference_routes.txt`).
`refs/gen_caller_reference.py` prints the rows of
`rev/prob/beta_lcdf_test.cpp` (`results/local/caller_reference.txt`).

The unit tests ran with the author's `stan_math_test.sh` jobs, with the
change on develop `dc66941969`.  "Swap" means: the job replaces the named
files by their develop versions, runs the suites, restores the files and
runs them again.  Logs in `results/tests/<tag>/`:

| tag | job | commit | what |
|---|---|---|---|
| `betafn_t3` | 21348406 | `b40012df74` | `prim/fun/beta_test`, `rev/prob/beta_lcdf_test`, `mix/fun/beta_test`; swap of `prim/fun/beta.hpp` |
| `betafn_o3` | 21348407, GPU 21348415 | `b40012df74` | 4 OpenCL suites on a V100; swap of the OpenCL device function and `elt_function_cl.hpp` |
| `betafn_c3` | 21348478 | `b40012df74` | cpplint, deps, 2 headers, doxygen 1.9.7 |
| `betafn_t2` | 21346673 | `a64c3679fa` | 34 suites (33 existing and `rev/prob/beta_lcdf_test`) |
| `betafn_p2` | 21346950 | `a64c3679fa` | generated `test/prob` with `-f cdf`, 6 families |
| `betafn_c2` | 21348143 | `a64c3679fa` | cpplint, deps, 3 headers, test-headers, doxygen 1.9.7 |
| `betafn_t4` | 21386527 | `372cebb637` | the 34 suites of `betafn_t2`; swap of `prim/fun/beta.hpp` on the 3 suites of `betafn_t3` |
| `betafn_o4` | 21386540, GPU 21386589 | `372cebb637` | the 4 OpenCL suites of `betafn_o3`, same swap |
| `betafn_c4` | 21386827 | `372cebb637` | cpplint, deps, 2 headers, test-headers, doxygen 1.9.7 |

`47159a62ca` differs from `372cebb637` only in a test comment.
`results/ci/windows_x64_rev_tests_b40012df74.txt` holds the failure lines
of the GitHub Actions job `rev tests (windows-latest)` at `b40012df74`
(run 41405748336): with tolerance 1e-12, `rev/prob/beta_lcdf_test` failed
with relative errors of 1.7e-8 (value at `b = 2^30`), 2.2e-8 and 3.7e-9
(derivative at `b = 2^30` and `2^31`).  On develop, `betafn_t4` gives
relative derivative errors of 2.0e-6 and 3.0e-6 (1e-7 tolerance, 2
`Failure` lines in `*_develop_rev_prob_beta_lcdf_test.run`).

- The 28 failed expectations are the `Failure` lines of the develop stage
  in `betafn_t3`: `beta_large_arguments` 20, `beta_huge_argument_underflow`
  2, `beta_int_arguments` 4, and the `beta_lcdf` test 2.  Count them with
  `grep -c Failure results/tests/betafn_t3/*_develop_*.run`.
- `betafn_o3`: `smtest_ocl_run_21348415.log` has `FAILED_CASE` lines for
  `prim_matches_reference`, `rev_value_matches_reference` and
  `edge_values_match_cpu` on the develop stage.  The last one fails
  because the swap keeps the branch CPU `beta()`, which is finite at
  1.7e308, while the develop device function gives NaN.
- `betafn_p2`: this directory has
  `smtest_21346950_prob_<family>_summary.txt` (programs built, gtest
  totals, failures) instead of the per-family compile logs; 21 + 14 + 9 +
  12 + 12 + 117 = 185 programs.
- For the cpplint logs of `betafn_c2` and `betafn_c3` and the
  test-headers log of `betafn_c2`, this directory keeps only the last 5
  lines (`*_tail.txt`); the `RESULT` lines of the main job logs give the
  status.
- clang-format 10: `results/local/clang_format.txt`
  (`stan_clang_format.sh branch origin/develop` in the clone at
  `b40012df74`, on the workstation).
