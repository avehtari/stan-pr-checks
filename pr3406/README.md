# stan-dev/math#3406: scripts behind the numbers

The pull request moves the normal family of log CDFs onto one
`std_normal_lcdf` kernel built on the `erfcx` helpers, and computes `Phi`,
the normal cdfs and the `exp_mod_normal` lcdf, lccdf and cdf from `erfcx`.
This directory holds the scripts behind the numbers in the PR description
and in the comments, and the probes made for the fact check of
2026-10-28 (stevieb-bot, commit `150f54df16`).

## Requirements

- Python 3 with mpmath 1.4 (most scripts set 50 to 700 digits).
- g++ with C++17.  Local probes: g++ 11.4.  Triton tests: g++ 14.3.
  Roihu timings: g++ 11.5.
- A Stan Math tree for each commit.  `build.sh` compiles one source against
  `STAN_MATH`; `STAN_MATH_LIB` and `TBB_LIB` point to a `lib/` and a built
  TBB when the tree is a `git archive` of `stan/` only.  Example:

  ```bash
  git -C <clone> archive dc66941969 stan | tar -x -C trees/dc66941969
  STAN_MATH=trees/dc66941969 STAN_MATH_LIB=<clone>/lib \
    TBB_LIB=<clone>/lib/tbb ./build.sh ranges/probe_test_points.cpp out/tp_dev
  ```

## Commits

| name | commit | note |
|---|---|---|
| develop (review) | `e4745c87ff` | merge base that the fact check used |
| develop (now) | `dc66941969` | merged into the branch on 2026-10-30 |
| PR as opened | `a28c1c1637` | andrjohns, 2026-09-18 |
| reviewed head | `150f54df16` | |
| merged head | `16c1ca18bf` | `150f54df16` + `dc66941969` |
| head with the new tests | `45c69aba9a` | `16c1ca18bf` + tests and comments |
| first merge | `1e46958b98` | not pushed; differs from `16c1ca18bf` only in `prim/prob/neg_binomial_2_lpmf.hpp`; the accuracy READMEs call it "head" |

Develop `dc66941969` and `e4745c87ff` have the same versions of the
headers that the accuracy scripts call, so the develop numbers hold for
both.

The job scripts in `speed/old/` and `steve/` are the scripts as they ran on
Triton and Roihu, with the cluster paths of that time (for example the
`lib/` directory in `/scratch/work/ave/`).  Set those paths for another
machine.

## Where each fact is backed

Facts in the order of the fact check.  "§" names a section of the README
file in the second column.

| fact | file § section |
|---|---|
| original cdfs 1.4 to 1.9 times develop's time | `README_speed.md` §2a |
| `9a8deffbef` slowed the `var` lcdfs | `README_speed.md` §2b |
| OpenCL `exp_mod_normal` timings | `README_speed.md` §2j |
| #2881 | unit tests (below) |
| valid input ranges | `README_ranges.md` §Valid input ranges |
| develop gradient error 7.6e-06 to 5.6e-05 | `README_normal.md` §N1 |
| `Phi` 7.7e5 ulp on [-5, -0.67] | `README_normal.md` §N3 |
| previous `exp_mod_normal` form fails | `README_emg.md` §M1 |
| deep-tail rows before `1f8b6ba94d` | `README_emg.md` §M4 |
| NaN at the median for `sigma` below 1 / DBL_MAX | unit test `rev/prob/normal_lcdf_small_sigma_test.cpp` |
| compile check of the three cdfs | unit test (below) |
| 256 propto checks (`9a8deffbef`) | `README_steve.md` §S1 |
| description timing tables | `README_speed.md` §1 |
| 601-point kernel table | `README_normal.md` §N1 |
| `erfcx` direct against `std_normal_erfcx` | `README_speed.md` §2l |
| wiener4 log Mills ratio | `README_ranges.md` §wiener4 |
| accuracy before and after `5f33a10d31` | `README_normal.md` §N2 |
| speed of `5f33a10d31` | `README_speed.md` §2f |
| `Phi` as opened, speed | `README_speed.md` §2e |
| `Phi` worst ulp, three versions | `README_normal.md` §N3 |
| three cdfs, accuracy | `README_normal.md` §N4 |
| bitwise identical results | `README_emg.md` §B1 |
| ten-function table on the EPYC 9965 | `README_speed.md` §2c |
| "the description measured a slow develop" | `README_speed.md` §2d |
| `erfcx` form of `Phi` more accurate | `README_normal.md` §N3 |
| `skew_normal_lpdf` speed (`927099efb8`) | `README_speed.md` §2g |
| OpenCL `Phi` and `std_normal_lcdf` timings | `README_speed.md` §2h |
| `exp_mod_normal` 510 to 900 ulp | `README_emg.md` §M2 |
| `exp_mod_normal` speed (`a17bba2e9e`) | `README_speed.md` §2i |
| OpenCL host comparison 5e-14 | `README_emg.md` §M3 |
| `e1fea0651c`, `cf177ef96b` speed | `README_speed.md` §2k |
| about 2.7 x 2^-52 | `README_normal.md` §N5 |
| 105-point grid | `README_emg.md` §M4 |
| compile memory 33 GB to 9 GB | `README_steve.md` §S2 |
| expression tests failed (`f37ea7cb67`) | `README_steve.md` §S3 |
| `exp_mod_normal_lpdf` at z ~ 9.4 | `README_steve.md` §S4 |
| varmat helper defects | `README_ranges.md` §varmat helpers |
| sign branch of `cf177ef96b` | `README_steve.md` §S5 |
| reference constants of the tests | `README_normal.md` §N6, `README_emg.md` §M5 |

## Unit tests added after the fact check

In the branch, not in this directory.  Swap test: the develop versions of
the 32 changed headers that exist on develop, copied in at once (Triton,
g++ 14.3, jobs 21301439 and 21302434).

| test | develop |
|---|---|
| `mix/prob/normal_ccdf_log_test` `normal_lccdf_infinite_endpoints` | fails, `sigma.adj()` NaN |
| `mix/prob/normal_cdf_log_test` `normal_lcdf_infinite_endpoints` | fails, `sigma.adj()` NaN |
| `prim/fun/Phi_test` `Phi(-4.988500030550561)` | fails, relative error 8.5e-11 |
| `mix/prob/normal_cdfs_argument_types_test` | passes (instantiation and value test) |
| `rev/prob/normal_family_range_ends_test` (4 tests) | all fail (`-inf`, NaN, 0) |
| `rev/fun/as_array_or_scalar_array_test` | does not compile |
| `rev/meta/ref_type_arena_matrix_test` | fails (2 `is_reference` checks) |
| `prim/fun/log1p_test` `log1p_double_containers` | passes (behavior of both versions) |

The other tests that the PR adds (the `exp_mod_normal` reference rows, the
branch accuracy points, the new `mix` and OpenCL tests) were checked on
develop `e4745c87ff` by the fact check itself.

Runs at the head (Triton, g++ 14.3, CPU builds on batch nodes):

| job | content | result |
|---|---|---|
| 21301439 | 67 CPU suites at `1e46958b98` | 66 pass; `rev/prob/neg_binomial_2_test` crashes (see "The merge of develop") |
| 21301978 | 5 `neg_binomial_2` suites and the 20 generated `test/prob/neg_binomial_2` programs at `16c1ca18bf` | pass (12 960 cases) |
| 21302434 | the 3 new `rev` suites and 5 changed suites, with the swap | pass; swap as in the table above |
| 21301623, 21301922 | 31 OpenCL suites at `1e46958b98` (V100) | 30 pass; `neg_binomial_2_lpmf_test` stops |
| 21302081 | `opencl/rev/neg_binomial_2_lpmf_test` at `16c1ca18bf` (V100) | pass |
| 21301869, 21302476 | cpplint, test-math-dependencies, header check of 52 headers, doxygen 1.9.7 | pass |
| 21303290, 21303811 | the 10 new and changed CPU suites at `45c69aba9a`, with the swap | pass; swap as in the table above |
| 21303288, 21303395 | 5 OpenCL `exp_mod_normal` suites and `elt_function_cl_test` at `45c69aba9a` (V100) | pass |
| 21303636 | cpplint, test-math-dependencies, header check, doxygen 1.9.7 at `45c69aba9a` | pass |

Points and references of the range tests: `ranges/probe_test_points.cpp`,
`ranges/ref_test_points.py` (two routes for log Phi at 700 digits).
Probe of the infinite endpoints and the `Phi` point:
`tests_probe/probe_e1_e2.cpp`.

## The merge of develop

`merge/probe_nb2.cpp` shows a heap-use-after-free in `neg_binomial_2_lpmf`
with develop's version of the function and the branch's vectorized
`log1p` (run it with `EXTRA_CXXFLAGS="-g -fsanitize=address"`).  The merge
commit keeps the branch's version of the function.

## Outputs

`out/` is not tracked.  Each README names the output files of its
scripts and quotes the numbers.
