# Checks for stan-dev/stan#3457 (partial inits, part of stan#3088)

This directory holds the scripts and logs behind the numbers in the PR
"Draw missing parameters of a partial init on the unconstrained scale"
(stan-dev/stan#3457, branch `fix-partial-init-dependent-transforms`,
addresses part of stan-dev/stan#3088) and in the documentation PR
stan-dev/docs#1004.  Until 2026-11-01 this directory was
`stan_partial_init/`; both PR texts link that name at commit `02b93f1`.
Sections 1 to 7 quote the claims of the PR text and give the commands and
the expected output.

Terms: a *supplied* parameter has a user initial value; a *missing*
parameter has none.

Last run: 2026-11-01.  The tests, the probe and the CmdStan runs ran on
Triton batch nodes (Intel Xeon, g++ 14.3.0); lint ran on a laptop.  In the logs, the
placeholders `<workstation>`, `<scratch>` and `<home>` replace the machine
paths.

## Commits

| name in this file | commit | what |
|---|---|---|
| develop | Stan `a8ec3d18f`, Stan Math `dc66941969` | merge base of the PR (develop of 2026-10-30) |
| PR | Stan `78b495834` | PR head |
| CmdStan develop | `3ab132c` | CmdStan with Stan `a8ec3d18f`, stanc3 v2.40.0 (section 1 and `cmdstan/` step 2) |

The unit-test job 21321656 (job name `smtest`) ran on the uncommitted
edits that became commit `3b0d6734b`; its log records
`uncommitted_edits=yes`.  The md5 values that the deploy reported for the
three files equal those of `3b0d6734b`.  The amend to `78b495834` changed
only the doc comment of `partial_init_unconstrained`: every changed line
is a comment line (`git diff 3b0d6734b 78b495834`).  The probe job
21342688 ran at `78b495834`.  Both jobs used the `bin/stanc` of the Stan
test tree, v2.39.0-15 (`results/run_probe_21342688.log`); the deploy
script copies it from CmdStan 2.40.0-rc0, whose stanc3 has that version.

The CmdStan runs used a build of CmdStan `3ab132c` with Stan `a8ec3d18f`
and Stan Math `dc66941969` from `git archive` tarballs
(`triton/run_3088.sh`).  The logs record CmdStan version 2.40.0 and the
md5 `49a48975` of `initialize.hpp`, which is the develop version.

## Layout

| path | content |
|---|---|
| `probe/` | `probe_init.cpp`, built with the develop version and with the PR version of `initialize.hpp`; `build.sh`; `compare_outputs.py`; the test model |
| `swap/swap_test.sh` | the swap test for a local Stan clone |
| `cmdstan/` | `run_3088.R` and the models of the CmdStan runs (develop only) |
| `grad/grad_check.cpp` | gradients of the offset/multiplier form at small scales (for the #3088 comment and the docs PR) |
| `triton/` | the Triton job scripts, with placeholders for the paths |
| `results/` | outputs of 2026-11-01 and the Triton logs |

Run commands from this directory.  New outputs go to `out/`, which git
does not track; compare them with `results/`.

## Requirements

- A Stan clone at the PR branch with `lib/stan_math` checked out and its
  TBB built (`make -f lib/stan_math/make/standalone math-libs` in the
  clone).  The sundials headers must include `arkode/`.
- A `stanc` binary.
- g++ with C++17, Python 3 (no packages).
- For `cmdstan/`: R with cmdstanr and posterior, and a CmdStan build of
  develop.

```sh
git -C <clone> show a8ec3d18f:src/stan/services/util/initialize.hpp \
  > <dir>/stan/services/util/initialize.hpp       # mkdir -p the path first
STAN=<clone> STANC=<stanc> probe/build.sh probe_develop <dir>
STAN=<clone> STANC=<stanc> probe/build.sh probe_pr
out/probe_develop > out/probe_develop.txt
out/probe_pr > out/probe_pr.txt
python3 probe/compare_outputs.py out/probe_develop.txt out/probe_pr.txt
STAN=<clone> swap/swap_test.sh
```

The probe uses the test model of the PR
(`probe/partial_init_dependent.stan`, a copy of
`src/test/test-models/good/services/partial_init_dependent.stan`).  For
each case it initializes with 20 seeds (1 seed for `init_radius = 0`; 3
seeds for `lb_10`) and prints per seed:

- `max_abs_missing`: the largest `|params[n] - u_random[n]|` over the
  missing entries, where `u_random` is a copy of the first random draw;
- `max_rel_supplied`: the largest relative difference between a supplied
  value and its value after `write_array`;
- `max_abs_unconstrained`: the largest `|params[n]|`;
- `rng_next`: the next draw of the RNG after `initialize`.

For the runs without partial inits (`none_r2`, `none_r0`, `full_r2`) it
also prints all parameters in hexadecimal for the bitwise comparison.

## Probe results (job 21342688, `results/compare.txt`)

| case | supplied | develop: max_abs_missing | PR: max_abs_missing | max_rel_supplied (both) |
|---|---|---|---|---|
| `sigma_1e-20` | `sigma = 1e-20` | 7.67e20 | 0 | 7.52e-16 |
| `mu_5` | `mu = 5` | 28.5 | 0 | 0 |
| `mu_5_r0` | `mu = 5`, `init_radius = 0` | 5 | 0 | 0 |
| `x_3` | `x = 3` | 6.66e-16 | 0 | 1.48e-16 |
| `a_1_c_4` | `a = 1`, `c = 4` | 11.3 | 0 | 0 |
| `b_2_d_-3` | `b = 2`, `d = -3` | 7.42 | 0 | 1.48e-16 |
| `lb_10` | `lb = 10` | throws, 100 rejections per seed | same | |

```
runs without partial inits: 41 of 41 parameter vectors bitwise equal
next RNG draw after initialize: 142 of 142 equal
```

## 1. The reprex of #3088 on develop

> With `real<lower=0> sigma; real<multiplier=sigma> x;` and `sigma =
> 1e-20` supplied (the reprex in #3088), `x` starts at unconstrained
> values from 1e19 to 6e20, and `lp__` at the initial point is -9e37 to
> -1.5e41. In all 20 chains (20 seeds, CmdStan develop), `sigma` did not
> change in 1000 draws.
> With `offset=mu` and `mu = 5` supplied, the unconstrained `x` goes up
> to 28 instead of staying in (-2, 2).

`cmdstan/run_3088.R`, step 1, on CmdStan develop (job 21320589, log
`results/triton_cmdstan/run_3088_21320589.log`).  The step calls the
`diagnose` method of each model binary with the partial init and reads
the unconstrained values and `lp__` of the initial point:

```
r1_m[sigma=1e-20]            init ok 20/20; rejected inits per seed: min 0 max 0
   unconstrained[2]: |min|  1.36e+19  |max|  5.53e+20
   lp at init: min -1.53e+41  max -9.32e+37
r1_om[mu=5]                  init ok 20/20; rejected inits per seed: min 0 max 0
   unconstrained[3]: |min|     0.694  |max|        28
r1_m       stuck (sd_sigma == 0) 20/20; lp_mean range -1.37e+41 -9.58e+36; stepsize med  2.72e-17
r1_manual  stuck (sd_sigma == 0)  0/20; lp_mean range     -1.79     -1.55; stepsize med     0.578
```

The two `stuck` lines are the sampling step 1b of the same job.  "Stuck"
means that `sigma` has standard deviation 0 over the 1000 draws.  The
probe shows the same defect with the test model of the PR (case
`sigma_1e-20`: unconstrained values up to 7.67e20 on develop, at most
46.1 = |log(1e-20)| with the PR).

## 2. Develop already handles the reverse dependency

> Develop already handles the reverse dependency. A supplied parameter
> whose transform uses missing parameters gets the correct unconstrained
> value, because develop frees it with the same random constrained values
> that the missing parameters start from.

Probe case `x_3` on develop: the missing entries differ from the random
draws by at most 6.66e-16 (relative 2.84e-15), and `x` comes back within
1.48e-16.

## 3. Missing and supplied parameters with the PR

> Missing parameters start exactly at the random unconstrained draws (at
> 0 with `init_radius = 0`), also when their offset or multiplier uses a
> supplied parameter. Supplied parameters keep their values on the
> constrained scale, up to the rounding of `transform_inits` followed by
> `write_array`, as on develop.

Probe, PR build: `max_abs_missing` is 0 in every case; `mu_5_r0` gives 0
(develop: 5).  `max_rel_supplied` is at most 7.52e-16 in both builds.

## 4. Runs without partial inits and the RNG stream

> Runs without supplied values or with all values supplied do not change,
> and the number of random draws per attempt is the same as on develop, so
> the RNG stream does not change.

Probe: the 41 parameter vectors of `none_r2` (20 seeds), `none_r0` and
`full_r2` (20 seeds) are bitwise equal between the builds.  The next RNG
draw after `initialize` is equal in all 142 runs except the 3 `lb_10`
runs, which throw.

## 5. The fixed-point iteration needs more than one pass

> More than one pass is needed when a supplied parameter depends on a
> missing parameter that depends on a supplied one.

Test `chain_provided_b_d` and probe case `b_2_d_-3`: `b` (supplied)
depends on `a` (missing), `c` (missing) on `b`, and `d` (supplied) on `c`.
Develop gives missing entries up to 7.42 from their random draws; the PR gives 0.

## 6. The bound case is not fixed

> With `real a; real<lower=a> y;` and `a = 10` supplied, the random `y`
> (with the default R = 2) is below 2 + e² ≈ 9.4, every attempt is rejected, and initialization
> fails.

The random `y` is `a_r + exp(u)` with `a_r` and `u` in (-2, 2), so it is
below 2 + e² ≈ 9.389.  Probe case `lb_10` (the model's `lb` and `y`):
both builds throw after 100 rejections per seed.

## 7. The unit tests and the swap test

> | `initialize.hpp` | `initialize_partial_test` | `initialize_test` |
> |---|---|---|
> | this PR | 7/7 pass | 16/16 pass |
> | develop | 1/7 pass | 16/16 pass |

Job 21321656 (`results/triton_tests/smtest_21321656.log` and the `.run`
files): the job ran each suite with the PR header, with the develop
header, and with the PR header again.  On develop,
`initialize_partial_test` fails `small_multiplier`, `provided_offset`,
`radius_zero`, `provided_dependent_parameter`, `chain_provided_a_c` and
`chain_provided_b_d`.  `provided_dependent_parameter` fails because its
check of the missing entries uses exact equality, and develop differs in
the last bits (section 2).  With the PR header, `fail_init_test` and the
33 suites in `src/test/unit/services/sample` also pass (one `.run` file
each).  `swap/swap_test.sh` runs the new tests with the develop header in a local
Stan clone.

Lint (`results/lint.txt`): clang-format 14 and cpplint 2.0.2 with the
filters of `make/cpplint` report nothing on the two C++ files.

## Other data in this directory

`cmdstan/` step 2 and `grad/` give the numbers of the #3088 comment and
of the documentation PR (stuck chains of the offset/multiplier form with
default inits; gradient errors at small scales).  Their logs are in
`results/triton_cmdstan/` (jobs 21320542, 21320544, 21321652) and
`results/grad_check_out.txt`.  `grad_check.cpp` builds like the probe,
against a Stan Math tree, with `-I <math> -I <math>/lib/eigen_5.0.1
-I <math>/lib/boost_1.87.0 -D_REENTRANT` and the TBB library.
