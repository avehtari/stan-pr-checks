# Speed claims of stan-dev/math#3406

Section 1 is a rerun of the timing tables in the PR description.  Section
2 lists the older speed claims of the review comments, with copies of
their scripts, job scripts and result files in `speed/old/`.  Section 2
was not run again.

Requirements: g++ with C++17, a Stan Math tree (`STAN_MATH`) and a `lib/`
with Eigen 5.0.1, Boost 1.87.0, TBB 2020.3 (built `libtbb.so`) and SUNDIALS
6.1.1 headers (`STAN_MATH_LIB`), as for `build.sh`.  Do not time on a
laptop: use one whole node of a cluster.

## 1. Rerun of the description's timing tables (2026-10-30)

Claim (PR description):

> While performance is either unchanged or significantly improved across inputs:

The description gives develop / branch ns per element for 15 functions,
four argument types and N = 100 to 100000.  It does not give the benchmark,
the machine, the compiler or the inputs.

### Files

| file | content |
|---|---|
| `speed/bench_one.cpp` | the benchmark; one function per binary (`-DBENCH_ONLY=<fn>`) |
| `speed/run_bench_roihu.sh` | job: builds 15 x 2 binaries, then 3 passes of (develop, head) per function |
| `speed/make_trees.sh` | makes the tree tarball with `git archive` |
| `speed/summarize.py` | tables, min and max per function, spread, value check |
| `speed/p3406_speed_2662883.log` | raw output of the job (`RESULT` and `CHECK` lines) |
| `speed/bench_pr_table.cpp` | the older multi-function benchmark (§2 b, c, d, g, i) |

`bench_one.cpp` extends `old/k_steve/bench_one.cpp`: the four argument
types of the description, N = 100 to 100000, the three lccdfs, and a
`CHECK` line per case.  Without `-DBENCH_ONLY` it does not compile: GCC's
inlining limits apply to the whole translation unit, so one binary for
many functions can change the time of a function that did not change.

Argument types (the four tables of the description):

| label | y | parameters | timed |
|---|---|---|---|
| `double` | `Eigen::VectorXd` | `double` | value |
| `var` | `Eigen::Matrix<var, -1, 1>` | `var` | value and `grad()` |
| `varmat` | `var_value<Eigen::VectorXd>` | `var` | value and `grad()` |
| `fvar` | `Eigen::Matrix<fvar<double>, -1, 1>`, tangent 1 | `fvar<double>`, tangent 0 | value and tangent |

Inputs: `z` uniform on [-3, 3] (`bulk`) or [-30, 30] (`wide`), seed
`20261020 + N`.  `Phi` and the `std_normal_*` functions take `z`; the
others take `y = mu + sigma z` with `mu = 0.3`, `sigma = 1.7`,
`lambda = 1`; the lognormal functions take `exp(mu + sigma z)`.  `Phi` is
timed as `sum(Phi(y))`.  The `var` and `varmat` times include the
construction of `y` on the autodiff stack and `recover_memory()`.  About
2e6 elements per `double` measurement and 1e6 per other measurement.  The
loop adds `acc * 1e-300` to the inputs, so that the compiler cannot hoist
the call (`STAN_MATH_DISTRIBUTIONS.md` §4); non-finite results do not
enter `acc`.

### Commands

```bash
CLONE=<stan math clone> speed/make_trees.sh out/speed/pr3406_speed_trees.tar.gz \
  dc66941969 16c1ca18bf
# deploy speed/bench_one.cpp, build.sh, the tarball and speed/run_bench_roihu.sh
# to one directory, then:
sbatch --export=ALL,CMDSTAN_LIB=<cmdstan-2.40.0>/stan/lib/stan_math/lib run_bench_roihu.sh
python3 speed/summarize.py speed/p3406_speed_2662883.log bulk
python3 speed/summarize.py speed/p3406_speed_2662883.log wide
```

One binary by hand:

```bash
STAN_MATH=<tree> EXTRA_CXXFLAGS=-DBENCH_ONLY=normal_lcdf ./build.sh \
  speed/bench_one.cpp out/speed/b_normal_lcdf
taskset -c 2 out/speed/b_normal_lcdf <label> 2
```

### Setup of the run

- Trees: develop `dc66941969` (2026-10-29) and the PR head `16c1ca18bf`
  (develop merged into the branch, so the diff is the PR only), from
  `git archive <commit> stan`.  md5 of all files under `stan/`:
  `c278326d628606c523f04f052eca7095` and
  `e8870beb6a210b5490aaf45c4b89237d` (printed by the job and by
  `make_trees.sh`).
- Job: CSC Roihu 2662883, 2026-10-30, one whole `medium` node `rc4281`
  (node exclusive), AMD EPYC 9965, governor `performance`, pinned to core
  2.  System g++ 11.5.0, glibc 2.34.
- Flags: `-std=c++17 -O3 -pthread -DBOOST_DISABLE_ASSERTS` (`build.sh`,
  Stan's default `-O3`, no `-march`, no `-DNDEBUG`).  Eigen, Boost and TBB
  from CmdStan 2.40.0.
- 3 passes x 2 rounds per binary; in each pass develop then head per
  function.  The tables give the minimum over the 6 rounds.

### Results, `bulk` (z on [-3, 3])

develop / head ns per element (speedup = develop / head, as in the
description):

| function | type | N = 100 | N = 1000 | N = 10000 | N = 100000 |
|---|---|---|---|---|---|
| normal_lcdf | double | 16.0 / 10.7 (1.50x) | 15.6 / 10.9 (1.43x) | 18.6 / 12.0 (1.55x) | 26.2 / 16.0 (1.64x) |
| normal_lcdf | var | 40.9 / 15.8 (2.60x) | 41.2 / 15.3 (2.70x) | 48.7 / 18.6 (2.62x) | 52.0 / 24.8 (2.10x) |
| normal_lcdf | varmat | 39.5 / 13.3 (2.96x) | 39.5 / 12.9 (3.07x) | 47.0 / 13.7 (3.44x) | 50.3 / 22.1 (2.27x) |
| normal_lcdf | fvar | 39.9 / 13.8 (2.89x) | 39.8 / 13.5 (2.94x) | 49.9 / 21.5 (2.32x) | 55.0 / 28.4 (1.94x) |
| normal_lccdf | double | 16.3 / 11.1 (1.47x) | 15.7 / 11.0 (1.43x) | 18.5 / 11.9 (1.55x) | 26.3 / 15.9 (1.65x) |
| normal_lccdf | var | 45.0 / 16.3 (2.77x) | 44.6 / 15.4 (2.90x) | 51.9 / 18.7 (2.78x) | 55.5 / 24.8 (2.24x) |
| normal_lccdf | varmat | 40.8 / 13.8 (2.95x) | 40.3 / 13.0 (3.09x) | 47.2 / 13.8 (3.42x) | 50.5 / 22.0 (2.29x) |
| normal_lccdf | fvar | 41.1 / 14.6 (2.82x) | 40.4 / 13.8 (2.93x) | 53.7 / 21.7 (2.47x) | 55.8 / 28.4 (1.97x) |
| normal_cdf | double | 10.8 / 7.5 (1.43x) | 10.4 / 7.4 (1.41x) | 10.9 / 7.4 (1.48x) | 18.0 / 12.7 (1.42x) |
| normal_cdf | var | 19.7 / 13.6 (1.45x) | 19.4 / 13.1 (1.49x) | 23.0 / 15.1 (1.52x) | 27.0 / 22.4 (1.20x) |
| normal_cdf | varmat | 18.1 / 11.0 (1.64x) | 17.6 / 10.5 (1.68x) | 20.7 / 10.3 (2.01x) | 24.2 / 15.9 (1.53x) |
| normal_cdf | fvar | 18.2 / 11.8 (1.53x) | 17.7 / 11.3 (1.56x) | 23.5 / 18.2 (1.29x) | 28.6 / 26.1 (1.10x) |
| std_normal_lcdf | double | 14.7 / 10.6 (1.38x) | 14.4 / 10.7 (1.35x) | 16.7 / 11.8 (1.41x) | 22.9 / 15.9 (1.44x) |
| std_normal_lcdf | var | 37.6 / 14.3 (2.62x) | 37.8 / 14.2 (2.65x) | 44.5 / 15.2 (2.92x) | 47.6 / 23.7 (2.01x) |
| std_normal_lcdf | varmat | 36.0 / 11.7 (3.08x) | 36.1 / 11.6 (3.10x) | 42.5 / 12.7 (3.35x) | 45.6 / 17.0 (2.68x) |
| std_normal_lcdf | fvar | 36.8 / 12.6 (2.93x) | 36.7 / 12.4 (2.95x) | 46.0 / 19.4 (2.38x) | 51.0 / 26.1 (1.95x) |
| std_normal_lccdf | double | 15.2 / 11.1 (1.37x) | 14.6 / 10.8 (1.34x) | 16.6 / 12.0 (1.39x) | 23.2 / 15.8 (1.47x) |
| std_normal_lccdf | var | 41.7 / 15.1 (2.77x) | 41.4 / 14.5 (2.85x) | 48.2 / 15.4 (3.13x) | 51.7 / 23.8 (2.17x) |
| std_normal_lccdf | varmat | 37.2 / 12.4 (3.00x) | 36.6 / 12.0 (3.06x) | 43.1 / 12.9 (3.33x) | 46.4 / 17.2 (2.70x) |
| std_normal_lccdf | fvar | 37.9 / 13.0 (2.90x) | 37.1 / 12.5 (2.97x) | 49.7 / 19.3 (2.57x) | 51.6 / 26.1 (1.98x) |
| std_normal_cdf | double | 9.7 / 7.3 (1.33x) | 9.5 / 7.0 (1.36x) | 9.6 / 7.0 (1.38x) | 15.1 / 12.3 (1.23x) |
| std_normal_cdf | var | 17.3 / 11.5 (1.50x) | 16.9 / 11.1 (1.53x) | 19.4 / 13.3 (1.46x) | 22.4 / 20.6 (1.09x) |
| std_normal_cdf | varmat | 15.3 / 9.4 (1.64x) | 14.9 / 9.0 (1.65x) | 17.4 / 8.9 (1.95x) | 20.4 / 14.5 (1.41x) |
| std_normal_cdf | fvar | 15.8 / 9.6 (1.64x) | 15.1 / 9.2 (1.63x) | 20.2 / 15.0 (1.34x) | 24.8 / 23.1 (1.07x) |
| Phi | double | 9.4 / 6.7 (1.40x) | 8.9 / 6.5 (1.38x) | 9.5 / 6.6 (1.43x) | 14.9 / 11.7 (1.27x) |
| Phi | var | 20.2 / 16.8 (1.20x) | 19.7 / 16.5 (1.20x) | 20.4 / 17.9 (1.14x) | 26.6 / 23.2 (1.14x) |
| Phi | varmat | 12.6 / 9.9 (1.27x) | 12.0 / 9.5 (1.27x) | 12.7 / 10.1 (1.26x) | 18.2 / 14.8 (1.23x) |
| Phi | fvar | 14.5 / 10.9 (1.34x) | 14.1 / 11.4 (1.24x) | 17.3 / 14.9 (1.16x) | 20.3 / 17.4 (1.17x) |
| lognormal_lcdf | double | 22.2 / 14.7 (1.50x) | 21.7 / 14.7 (1.47x) | 28.5 / 15.8 (1.81x) | 33.4 / 20.1 (1.66x) |
| lognormal_lcdf | var | 22.7 / 20.3 (1.12x) | 21.8 / 19.6 (1.11x) | 24.4 / 22.9 (1.07x) | 32.9 / 29.2 (1.13x) |
| lognormal_lcdf | varmat | 20.3 / 17.8 (1.14x) | 19.4 / 17.2 (1.13x) | 19.9 / 18.0 (1.10x) | 26.9 / 26.6 (1.01x) |
| lognormal_lcdf | fvar | 21.1 / 18.2 (1.16x) | 20.1 / 17.7 (1.14x) | 25.3 / 25.7 (0.98x) | 33.7 / 31.7 (1.06x) |
| lognormal_lccdf | double | 21.7 / 15.1 (1.44x) | 21.6 / 14.9 (1.45x) | 29.0 / 15.8 (1.83x) | 33.2 / 20.2 (1.65x) |
| lognormal_lccdf | var | 22.8 / 20.6 (1.10x) | 21.8 / 19.7 (1.11x) | 24.3 / 23.0 (1.06x) | 32.9 / 29.1 (1.13x) |
| lognormal_lccdf | varmat | 20.1 / 18.5 (1.09x) | 19.2 / 17.4 (1.10x) | 19.6 / 18.2 (1.07x) | 26.4 / 26.5 (1.00x) |
| lognormal_lccdf | fvar | 20.9 / 19.2 (1.09x) | 19.9 / 18.1 (1.10x) | 25.0 / 26.2 (0.95x) | 33.1 / 32.0 (1.03x) |
| lognormal_cdf | double | 17.7 / 11.4 (1.55x) | 17.7 / 11.1 (1.60x) | 22.4 / 11.1 (2.03x) | 29.9 / 16.4 (1.82x) |
| lognormal_cdf | var | 19.6 / 17.9 (1.10x) | 18.7 / 17.0 (1.10x) | 21.3 / 19.2 (1.11x) | 30.1 / 26.4 (1.14x) |
| lognormal_cdf | varmat | 17.0 / 15.8 (1.08x) | 16.1 / 14.8 (1.08x) | 16.4 / 14.6 (1.12x) | 23.0 / 23.8 (0.96x) |
| lognormal_cdf | fvar | 17.6 / 15.9 (1.11x) | 16.7 / 14.9 (1.12x) | 21.9 / 22.0 (1.00x) | 30.2 / 28.9 (1.04x) |
| exp_mod_normal_lpdf | double | 13.8 / 10.1 (1.37x) | 13.1 / 9.8 (1.33x) | 13.3 / 10.5 (1.27x) | 19.8 / 12.7 (1.55x) |
| exp_mod_normal_lpdf | var | 20.5 / 15.7 (1.31x) | 19.4 / 14.6 (1.33x) | 24.6 / 17.3 (1.42x) | 28.6 / 21.8 (1.31x) |
| exp_mod_normal_lpdf | varmat | 17.9 / 12.2 (1.48x) | 17.0 / 12.0 (1.41x) | 19.2 / 12.5 (1.54x) | 25.9 / 15.4 (1.69x) |
| exp_mod_normal_lpdf | fvar | 18.8 / 13.7 (1.37x) | 17.8 / 12.7 (1.40x) | 25.0 / 19.1 (1.31x) | 32.3 / 24.2 (1.33x) |
| exp_mod_normal_lcdf | double | 25.4 / 18.4 (1.38x) | 24.2 / 18.8 (1.29x) | 30.6 / 24.2 (1.27x) | 39.0 / 30.4 (1.28x) |
| exp_mod_normal_lcdf | var | 37.0 / 29.6 (1.25x) | 33.4 / 29.7 (1.12x) | 49.7 / 35.9 (1.38x) | 55.8 / 39.7 (1.41x) |
| exp_mod_normal_lcdf | varmat | 34.2 / 26.9 (1.27x) | 30.8 / 27.0 (1.14x) | 44.1 / 33.2 (1.33x) | 48.8 / 36.8 (1.33x) |
| exp_mod_normal_lcdf | fvar | 35.3 / 27.5 (1.29x) | 31.3 / 27.5 (1.14x) | 49.7 / 36.2 (1.37x) | 56.1 / 41.8 (1.34x) |
| exp_mod_normal_lccdf | double | 25.4 / 21.2 (1.20x) | 24.3 / 21.2 (1.14x) | 30.7 / 27.0 (1.13x) | 39.0 / 31.6 (1.23x) |
| exp_mod_normal_lccdf | var | 37.0 / 32.4 (1.14x) | 33.4 / 33.1 (1.01x) | 49.9 / 38.2 (1.30x) | 55.8 / 40.9 (1.36x) |
| exp_mod_normal_lccdf | varmat | 33.9 / 29.8 (1.14x) | 30.7 / 30.4 (1.01x) | 44.4 / 35.8 (1.24x) | 49.0 / 38.1 (1.29x) |
| exp_mod_normal_lccdf | fvar | 35.4 / 30.4 (1.17x) | 31.3 / 30.9 (1.01x) | 49.6 / 38.9 (1.28x) | 55.8 / 43.2 (1.29x) |
| exp_mod_normal_cdf | double | 21.9 / 16.4 (1.34x) | 20.8 / 16.6 (1.25x) | 27.0 / 20.3 (1.33x) | 35.5 / 28.7 (1.24x) |
| exp_mod_normal_cdf | var | 33.8 / 25.3 (1.33x) | 29.8 / 25.4 (1.17x) | 46.3 / 30.6 (1.51x) | 52.3 / 34.8 (1.50x) |
| exp_mod_normal_cdf | varmat | 31.5 / 22.9 (1.38x) | 27.8 / 22.6 (1.23x) | 41.5 / 27.9 (1.49x) | 46.3 / 31.8 (1.46x) |
| exp_mod_normal_cdf | fvar | 32.6 / 23.5 (1.39x) | 28.1 / 23.6 (1.19x) | 46.6 / 35.9 (1.30x) | 52.7 / 41.7 (1.26x) |
| skew_normal_lpdf | double | 13.5 / 10.6 (1.28x) | 13.1 / 10.7 (1.22x) | 14.9 / 11.9 (1.25x) | 19.9 / 15.8 (1.26x) |
| skew_normal_lpdf | var | 21.2 / 15.6 (1.36x) | 20.2 / 15.1 (1.33x) | 24.4 / 18.3 (1.33x) | 31.3 / 24.8 (1.27x) |
| skew_normal_lpdf | varmat | 18.7 / 13.4 (1.39x) | 17.8 / 13.1 (1.36x) | 19.7 / 14.0 (1.41x) | 28.7 / 22.4 (1.29x) |
| skew_normal_lpdf | fvar | 19.3 / 13.8 (1.40x) | 18.3 / 13.4 (1.37x) | 26.3 / 21.2 (1.24x) | 32.9 / 28.0 (1.17x) |

Smallest and largest speedup per function, `bulk` / `wide`:

| function | bulk min | bulk max | wide min | wide max |
|---|---|---|---|---|
| normal_lcdf | 1.43 double 1000 | 3.44 varmat 10000 | 2.65 fvar 1e5 | 5.89 varmat 10000 |
| normal_lccdf | 1.43 double 1000 | 3.42 varmat 10000 | 2.68 fvar 1e5 | 5.84 varmat 10000 |
| normal_cdf | 1.10 fvar 1e5 | 2.01 varmat 10000 | 1.13 fvar 1e5 | 2.02 varmat 10000 |
| std_normal_lcdf | 1.35 double 1000 | 3.35 varmat 10000 | 2.74 fvar 1e5 | 6.01 varmat 10000 |
| std_normal_lccdf | 1.34 double 1000 | 3.33 varmat 10000 | 2.76 double 1e5 | 5.87 varmat 10000 |
| std_normal_cdf | 1.07 fvar 1e5 | 1.95 varmat 10000 | 1.08 fvar 1e5 | 1.93 varmat 10000 |
| Phi | 1.14 var 10000 | 1.43 double 10000 | 1.04 fvar 10000 | 1.17 fvar 100 |
| lognormal_lcdf | 0.98 fvar 10000 | 1.81 double 10000 | 0.98 varmat 1e5 | 2.10 double 10000 |
| lognormal_lccdf | 0.95 fvar 10000 | 1.83 double 10000 | 0.95 varmat 1e5 | 2.05 double 10000 |
| lognormal_cdf | 0.96 varmat 1e5 | 2.03 double 10000 | 0.91 varmat 1e5 | 2.01 double 10000 |
| exp_mod_normal_lpdf | 1.27 double 10000 | 1.69 varmat 1e5 | 1.17 var 1e5 | 1.51 varmat 10000 |
| exp_mod_normal_lcdf | 1.12 var 1000 | 1.41 var 1e5 | (*) | (*) |
| exp_mod_normal_lccdf | 1.01 var 1000 | 1.36 var 1e5 | 1.09 var 1000 | 1.63 double 10000 |
| exp_mod_normal_cdf | 1.17 var 1000 | 1.51 var 10000 | 1.02 var 1000 | 1.51 var 1e5 |
| skew_normal_lpdf | 1.17 fvar 1e5 | 1.41 varmat 10000 | 1.14 fvar 1e5 | 1.76 fvar 10000 |

(*) For `wide` inputs develop's `exp_mod_normal_lcdf` returns NaN (value
and derivatives) for every N: it is non-finite for `z <= -10`.  The wide
row of this function times a NaN path on develop and is not a valid
comparison (summarize.py prints 1.00 to 1.33).

Notes:

- Spread between passes (largest per-pass minimum / smallest - 1):
  median 0.4 % (bulk) and 0.5 % (wide) over 480 cells each; largest 6.3 %
  (develop `exp_mod_normal_lpdf`, `double`, N = 10000, bulk); all other
  bulk cells at most 3.6 %, wide cells at most 4.3 %.  So differences above
  about 5 % are real.
- Cycle check: the shortest time, 6.5 ns for `Phi(double)`, is about 24
  cycles at 3.7 GHz; no call was removed by the compiler.
- `CHECK` lines: develop and head values agree to 3e-14 relative or
  better, except `exp_mod_normal_lcdf` (1.1e-13) and `exp_mod_normal_cdf`
  (2.6e-11 at a product near 1e-101, where develop is inaccurate).  The
  wide-range cdfs underflow to 0 in both trees.
- Compared with the description: the direction is the same (head faster
  in 235 of 240 bulk cells and 236 of 240 wide cells), but the speedups
  are smaller.  Description
  min / max, this run (bulk) min / max: `normal_lcdf` 1.66 / 3.56 against
  1.43 / 3.44; `Phi` 1.17 / 1.65 against 1.14 / 1.43; `lognormal_lcdf`
  1.37 / 2.60 against 0.98 / 1.81; `exp_mod_normal_lpdf` 1.78 / 2.36
  against 1.27 / 1.69; `skew_normal_lpdf` 1.58 / 2.49 against 1.17 / 1.41.
  The head is slower than develop only for the lognormal functions with
  `varmat` or `fvar` at N >= 10000: 0.95 to 1.00 (bulk), 0.91 to 1.00
  (wide; 0.91 is `lognormal_cdf`, `varmat`, N = 100000).  The
  description's develop times are 1.7 to 3.1 times the EPYC times (§2 d).

## 2. Older speed claims (copied, not run again)

The files are copies of the originals in the workstation directories
`normal_lcdf_work/{audit,dedup,bench,steve_1022,rescued_tmp}/`.  A
`deployed*/` subdirectory holds the copies from the job directory on the
cluster, where they differ from the local copy or where only the cluster
had them.  The job scripts are unchanged and contain the cluster paths of
the time.  The tree tarballs are in `out/speed/old_trees/` (not tracked);
`t_<name>/stan` per tree.  Unless a section says otherwise: g++ 14.3.0
(Triton, `scicomp-r-env`), flags `-O3 -DBOOST_DISABLE_ASSERTS`, one core,
binaries alternated, `RESULT` lines in the logs.

### a. Original PR cdfs 1.4 to 1.9 times develop's time

> The original PR also computes these three as `exp(lcdf)`, and with `double` arguments they take 1.4 to 1.9 times as long as on develop.

([comment 6276777866](https://github.com/stan-dev/math/pull/3406#issuecomment-6276777866))

- Files `speed/old/a_cdfs/`: `bench_cdfs.cpp`, `run_bench_cdfs.sh`
  (develop against head), `run_bench_cdfs3.sh` (develop, `exp(lcdf)`,
  direct), `cdfs_bench_21064786_results.txt`,
  `cdfs3_bench_21064914_results.txt`, the two job logs.
- Trees: develop `e4745c87ff` (`develop_e4745c87ff_stan.tar.gz`); `exp(lcdf)`
  = head `686fc67261` + `5f33a10d31` + `bf5861ab9c`
  (`commit_bf5861ab9c.tar.gz`); direct = + `0d4ac36e34`
  (`direct_cdfs_stan.tar.gz`).
- Machine: Triton `pe81`, Xeon E5-2680 v3; jobs 21064786 and 21064914;
  3 runs x 3 rounds.
- The quoted 1.4 to 1.9 are the `double`, K = 1 rows of job 21064914
  (1.45 to 1.92).  Job 21064786 gives 1.39 to 1.97 for the same rows.  With
  K = 64 the ratios are 0.85 to 1.27.  Spread up to 17.7 % (21064786).

### b. `9a8deffbef` caused the `var` lcdf slowdown

> One commit is responsible: 9a8deffbef ("Refactor normal log-CDF dispatch and density expressions") replaced evaluated temporaries with lazy expressions.

([comment 6285594457](https://github.com/stan-dev/math/pull/3406#issuecomment-6285594457))

- Files `speed/old/b_bisect/`: `deployed_bisect/` (`run_bench_bisect.sh`,
  `bench_pr_table.cpp` as run, `bisect_bench_21065499.log`),
  `deployed_fixcheck/` (the same for job 21065590),
  `bisect_bench_21065499_{results,summary}.txt`,
  `fixcheck_bench_{results,summary}.txt`, and the later local
  `run_bench_bisect.sh`.
- Command: `bench_pr_table <tree> 3 quick` (N = 1000, bulk, `double` and
  `Matrix<var>`), 2 passes x 3 rounds.
- Trees: job 21065499 every non-formatting commit `a28c1c1637` ...
  `f4cd10ee4a` (`bisect_trees.tar.gz`) plus develop `e4745c87ff`; job
  21065590 `a28c1c1637`, `686fc67261`, `f4cd10ee4a`, fix `c30ee6f9ef`
  (`fix_trees.tar.gz`).
- Machines: Triton `pe33` (21065499) and `pe2` (21065590), Xeon E5-2680 v3.

### c. Ten-function table on the EPYC 9965

> With the commits from my comments, the PR is now faster than develop for every function in the table except `exp_mod_normal_lcdf`, apart from `exp_mod_normal_lpdf` and `skew_normal_lpdf` with `var`, which are at develop's time.

([comment 6285594457](https://github.com/stan-dev/math/pull/3406#issuecomment-6285594457))

- Files `speed/old/c_progression/`: `run_bench_progression_roihu.sh`,
  `prog_roihu_2379049_{results,summary}.txt`, `deployed/` (the job script
  and `bench_pr_table.cpp` as run, log 2379049, log 2379005 of the aborted
  first attempt).
- Command: `bench_pr_table <tree> 2` (full mode: 10 functions, N = 100 to
  100000, bulk and wide), 2 passes x 2 rounds.
- Trees: develop `e4745c87ff`, `a28c1c1637` (as opened), `686fc67261`,
  `9b38ecbcb0` (`progression_trees.tar.gz`).
- Machine: Roihu job 2379049, whole `medium` node `rc5233`, AMD EPYC 9965,
  system g++ 11.5.0.  Spread at most 5.5 % at N = 1000.

### d. "The description measured a slow develop"

> So the description measured a slow develop, not a faster branch.

([comment 6285594457](https://github.com/stan-dev/math/pull/3406#issuecomment-6285594457))

- Files `speed/old/d_slow_develop/`: `run_bench_pr_table.sh` (Triton),
  `run_bench_pr_table_roihu.sh`, `prtab_bench_21065374_{results,summary}.txt`,
  `prtab_roihu_2374342_results.txt`, `prtab_roihu_vs_triton_summary.txt`,
  `deployed_triton/` and `deployed_roihu/` (`bench_pr_table.cpp` as run,
  logs).
- Trees: develop `5252d51d47` (base of the PR at creation), `a28c1c1637`,
  develop `e4745c87ff`, head `f4cd10ee4a` (`develop_5252d51d47_stan.tar.gz`,
  `pr3406_a28c1c1637_stan.tar.gz`, `develop_e4745c87ff_stan.tar.gz`,
  `commit_f4cd10ee4a.tar.gz`).
- Machines: Triton job 21065374, `pe81`, Xeon E5-2680 v3 (2 runs x 2
  rounds); Roihu job 2374342, `rc6171`, AMD EPYC 9965, g++ 11.5.0.
- The claim is an inference: absolute ns of the description against
  these runs (develop 1.7 to 3.1 times the EPYC time, the branch 0.80 to
  1.06 times).  The description's platform is not known.

### e. `Phi` as opened: 4.4x, 1.9x, 1.4x

> it takes 4.4 times as long for `|x| < 0.66`, 1.9 times on [-3, 3] and on [-37, -5], and 1.4 times on [5, 8].

([comment 6276090897](https://github.com/stan-dev/math/pull/3406#issuecomment-6276090897))

- Files `speed/old/e_phi/`: `deployed_21064693/` (`bench_phi_variants.cpp`,
  `phi_variants.hpp`, `run_bench_phi_variants.sh` as run; the local copies
  in `e_phi/` are later versions with more variants),
  `phiv_bench_21064693_results.txt` and log; `bench_phi.cpp`,
  `run_bench_phi.sh`, `phi_bench_21064615_results.txt` and log.
- Trees: develop `e4745c87ff` and head `686fc67261` + `5f33a10d31`
  (`commit_5f33a10d31.tar.gz`).  `bench_phi_variants` builds one binary
  on the second tree; the variants (`dev`, `pr`, `d1`, `d2`) are copies of
  the formulas in `phi_variants.hpp`.  `bench_phi` builds one binary per
  tree.
- Machine: Triton `pe81`, Xeon E5-2680 v3; 3 runs x 3 rounds.
- The `double` ratios of the comment (4.4, 1.9, 1.9, 1.4) are those of job
  21064693 (`pr` / `dev`: 4.39, 1.88, 1.91, 1.42).  Job 21064615
  (`bench_phi.cpp`, two trees) gives 3.99, 1.89, 1.87, 1.44.  The
  comment's `Phi(var)` range "1.2 to 1.5" is from job 21064615 (1.18 to
  1.54).

### f. `5f33a10d31` speed 0.66 to 1.02

> Time relative to before (below 1 is faster):

([comment 6274784156](https://github.com/stan-dev/math/pull/3406#issuecomment-6274784156))

- Files `speed/old/f_dedup/`: `bench_dedup.cpp`, `run_bench_commit.sh`,
  `c3406_bench_21064232_results.txt` and log (local and deployed copies
  are identical).
- Trees: head `686fc67261` (`pr3406/stan` in `dedup_trees.tar.gz`) and
  `686fc67261` + `5f33a10d31` (`commit_5f33a10d31.tar.gz`).
- Machine: Triton job 21064232, `pe81`, Xeon E5-2680 v3; 3 runs x 3
  rounds; the benchmark times `internal::std_normal_lcdf_value_grad`.

### g. `skew_normal_lpdf`: 0.85 with `var`, 7 % with `double`

> That made `double` calls 7 % slower, because there `z` is read only twice and the allocation costs more than the second read.

([comment 6291485931](https://github.com/stan-dev/math/pull/3406#issuecomment-6291485931))

- Files `speed/old/g_skew_normal/`: `lpdf_triton_results.txt`,
  `lpdf_triton2_results.txt`, `lpdf_triton_summary.txt`;
  `deployed_21069776/` and `deployed_21069800/` (`run_bench_bisect.sh`,
  `bench_pr_table.cpp` as run, logs).  The submit script was found in the
  Triton job directories `~/lpdf_triton*/`: it is the bisection script of
  §b with `TREES=a28c1c1637:9b38ecbcb0:fix`.
- Trees: develop `e4745c87ff`, `a28c1c1637`, `9b38ecbcb0`, fix
  (`lpdf_trees.tar.gz`: first version, all temporaries evaluated;
  `lpdf_trees2.tar.gz`: `927099efb8`, md5 of `skew_normal_lpdf.hpp`
  checked).
- Machine: Triton `csl3`, Xeon Gold 6248; jobs 21069776 (first version,
  the 7 %) and 21069800 (`927099efb8`, the 0.85); `quick` mode, 2 passes x 3
  rounds.
- The posted 0.84 / 0.96 / 0.85 (`var`) are medians or means of the 6
  rounds.  The minimum gives 0.85 / 0.98 / 0.87: develop was 6 % faster in
  pass 2 than in pass 1 (44.5 against 42.0 ns).  The 7 % is 0.85 -> 0.91
  of develop's time (`double`, `lpdf_triton_summary.txt`).

### h. OpenCL device timings of `Phi` and `std_normal_lcdf` on a V100

> I also timed the OpenCL device functions on a Tesla V100, with each work item walking through the input range in registers and a baseline kernel subtracted.

([comment 6291485931](https://github.com/stan-dev/math/pull/3406#issuecomment-6291485931))

- Files `speed/old/h_ocl_device/`: `bench_ocl_device.cpp`,
  `run_bench_ocl_device.sh` (job 21070484, identical on the cluster),
  `ocl_dev_bench2_{results,summary}.txt`;
  `deployed_21070459/` (first version), `ocl_dev_bench_{results,summary}.txt`
  and log 21070459.  `ocl_dev_bench2_acc.txt` (the per-point device values)
  and the log of job 21070484 are 14 MB each and are not in the
  repository.
- Trees: develop `e4745c87ff`, `a28c1c1637`, `9b38ecbcb0`
  (`progression_trees.tar.gz`).
- Machine: Triton `dgx1` (21070484) and `dgx20` (21070459), Tesla
  V100-SXM2-16GB; 2 passes x 3 rounds, median of 6, ps per evaluation.

### i. `a17bba2e9e` makes the `exp_mod_normal` cdfs faster (EPYC)

> The change also makes the functions faster.

([comment 6291485931](https://github.com/stan-dev/math/pull/3406#issuecomment-6291485931))

- Files `speed/old/i_emg_cpu/`: `emg_roihu{2,3,4}.txt` (from the
  workstation `/tmp`), `deployed_roihu/` (job script, logs 2385189,
  2386087, 2386276, 2386768), `emg_bench_triton_{results,summary}.txt`,
  `deployed_triton/` (logs 21071563, 21071704, 21071862, 21164745).
- Command: `bench_pr_table <tree> 3 emg` (N = 1000, bulk and wide).
  `speed/bench_pr_table.cpp` is the latest version (identical with the
  Roihu and Triton job directories); the exact 2026-10-20 version of the
  Roihu runs is not kept (the file was overwritten on 2026-10-23).
- Trees (`emg_trees.tar.gz`): develop `e4745c87ff`, `9b38ecbcb0`
  ("before"), `emg` = `a17bba2e9e` ("after"; md5 of
  `exp_mod_normal_lcdf.hpp` checked), `emg1` (an earlier version).
- Machine: the posted table is job 2386768 (`emg_roihu4.txt`), Roihu
  `rc5223`, AMD EPYC 9965, g++ 11.5.0; it matches to 0.01.  "Xeon Gold 6230
  gives the same pattern": Triton job 21071563
  (`emg_bench_triton_results.txt`).

### j. OpenCL `exp_mod_normal` timings (V100, N = 1e6)

> The `var` path is slower than before because each of the four wrappers runs the whole device function, and the device compiler does not merge the four copies.

([comment 6294609254](https://github.com/stan-dev/math/pull/3406#issuecomment-6294609254))

- Files `speed/old/j_ocl_emg/`: `bench_ocl_emg.cpp`, `run_bench_ocl_emg.sh`,
  `emg_ocl_bench.txt` (= log 21079254), log 21079254.
- Trees (`emg_ocl_trees.tar.gz`): develop `e4745c87ff`, `9b38ecbcb0`,
  `emg` = the OpenCL port `39886e0eb8` (md5 of
  `opencl/prim/exp_mod_normal_lcdf.hpp` checked).
- Machine: Triton job 21079254, `dgx3`, Tesla V100-SXM2-32GB; 2 passes x 3
  rounds; spread 10 to 60 %.
- Not measured: whether the device compiler merges the four copies.  The
  claim rests on the time difference only.

### k. `e1fea0651c` 1.2 to 1.3x and `cf177ef96b` 1.2 to 1.4x

> `exp_mod_normal_cdf` is 1.2 to 1.3 times faster after e1fea0651c, and after cf177ef96b `std_normal_lcdf`, `normal_lcdf`, `lognormal_lcdf`, `exp_mod_normal_lpdf` and `skew_normal_lpdf` are 1.2 to 1.4 times faster for z on [-30, 30].

([comment 6337624690](https://github.com/stan-dev/math/pull/3406#issuecomment-6337624690))

- Files `speed/old/k_steve/`: `bench_one.cpp`, `run_bench_one_roihu.sh`,
  `s3406_one_2444780.log`, `s3406_one_2444780_results.txt`,
  `summarize_bench.py`.
- Trees (`s3406_trees.tar.gz`): `e4745c87ff`, `e8bb9a0b38`, `e1fea0651c`,
  `cf177ef96b`.
- Machine: Roihu job 2444780, `rc4227`, AMD EPYC 9965, g++ 11.5.0;
  default and `-march=native` flags; 2 passes x 3 rounds; spread at most
  12 %.  Command: `python3 summarize_bench.py s3406_one_2444780_results.txt`.

### l. September comment: `erfcx` direct against the internal `std_normal_erfcx`

> Comparing your internal `std_normal_erfcx` path against calling `erfcx` directly, value and gradient computed together:

([comment 5732263602](https://github.com/stan-dev/math/pull/3406#issuecomment-5732263602))

- Files `speed/old/l_erfcx_direct/`: `bench_normal_lcdf_kernels.cpp`,
  `bench_normal_lcdf_ocl.cpp`, `run_normal_lcdf_cpu.sh`,
  `run_normal_lcdf_ocl.sh`, and all 7 logs of the Triton job directory
  `~/normal_lcdf_bench/`.  The logs were not on the workstation; the
  Triton copies were found on 2026-10-30.  Sources and job scripts are
  identical on both sides.
- The posted table is section "PR 3406 vs the erfcx identity" of job
  20331463 (CPU, column `ratio` of the throughput loop: 1.16, 1.10, 1.06,
  1.11) and job 20331487 (GPU: 1.35, 2.21, 1.03, 1.93).
- The benchmark holds a copy of the PR kernel of 2026-09-18 (namespace
  `pr3406`) and calls `erfcx` from the `erfcx` PR tree of that day.  It is
  a kernel benchmark without autodiff.  The kernel was removed in
  `5f33a10d31`, so it cannot run against the head.
- Machines: 20331463 Triton `csl38`, Xeon Gold 6248; 20331487 `dgx1`,
  Tesla V100-SXM2-16GB.
- The comment gives only the z < 0 rows.  The same section has two more
  rows where the direct call is slower: [0.1, 2.9] 0.58 (CPU) and 0.72
  (GPU), [2.9, 6] 0.57 and 0.61.
