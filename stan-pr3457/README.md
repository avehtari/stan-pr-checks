# Checks for stan-dev/stan#3457 (partial inits, part of stan#3088)

This directory holds the scripts and logs behind the numbers in the PR
"Draw missing parameters of a partial init on the unconstrained scale"
(stan-dev/stan#3457, branch `fix-partial-init-dependent-transforms`,
addresses part of stan-dev/stan#3088) and in the documentation PR
stan-dev/docs#1004.  Until 2026-11-01 this directory was
`stan_partial_init/`; the docs PR text and the first #3457 text link
that name at commit `02b93f1`.  The 2026-11-02 #3457 text links
`stan-pr3457/`.
Sections 1 to 7 quote the claims of the PR text and give the commands and
the expected output.

Terms: a *supplied* parameter has a user initial value; a *missing*
parameter has none.

Last run: 2026-11-03 (the two `initialize` suites, lint and the header
test at the PR head `2324dfd29`; the probe and the other unit tests at
`0f517b597`); 2026-11-02 (PR v3, kept in `*_2026-11-02/`); 2026-11-01
(CmdStan runs, sections 1, 2 and 8).
The tests, the probe and the CmdStan runs ran on
Triton batch nodes (Intel Xeon, g++ 14.3.0).  Lint ran on a laptop on
2026-11-01; from 2026-11-02 cpplint ran on the Triton login node and
clang-format on a laptop.  In the logs, the
placeholders `<workstation>`, `<scratch>` and `<home>` replace the machine
paths.

## Commits

| name in this file | commit | what |
|---|---|---|
| develop | Stan `a8ec3d18f`, Stan Math `dc66941969` | merge base of the PR (develop of 2026-10-30) |
| PR v1 | Stan `78b495834` | first PR head (2026-11-01 results) |
| PR v2 | Stan `7a63912c4` | flat vectors and `unconstrain_array` instead of `get_dims` |
| PR v3 | Stan `a5786800d` | one decision per tuple slot (2026-11-02 results; the jobs ran on the uncommitted files, same md5) |
| PR v4 | Stan `0f517b597` | rejects initial values with only some slots of a tuple (2026-11-03 results) |
| PR | Stan `2324dfd29` | third review: clearer names, shorter comments, local RNG, no size fallback (`triton_tests_2026-11-03_final/`) |
| (squashed) | `eb235255b`, `d25d70673`, `93f83b502`, `5d935d206` | unpushed commits of the third review, squashed into `2324dfd29` with later comment edits; their logs are in `triton_tests_2026-11-03_review3/` |
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
| `results/` | outputs of 2026-11-01 and the Triton logs (`triton_cmdstan_pr/`: section 8); 2026-11-02 (PR v3) in `probe_2026-11-02/` and `triton_tests_2026-11-02/`; 2026-11-03 (PR v4) in `probe_2026-11-03/` and `triton_tests_2026-11-03/`; PR in `triton_tests_2026-11-03_final/`; the squashed commits in `triton_tests_2026-11-03_review3/` |
| `random_var_context_tuples/` | probe of the develop defect of section 9 |

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

## Probe results (job 21365693, `results/probe_2026-11-03/compare_output.txt`)

PR at `0f517b597`, built with `tools/stan_probe.sh` of the workspace
(develop headers vs the PR headers), stanc v2.39.0-15:

| case | supplied | develop: max_abs_missing | PR: max_abs_missing | max_rel_supplied (PR) |
|---|---|---|---|---|
| `sigma_1e-20` | `sigma = 1e-20` | 1.05e21 | 0 | 7.52e-16 |
| `mu_5` | `mu = 5` | 28.5 | 0 | 0 |
| `mu_5_r0` | `mu = 5`, `init_radius = 0` | 5 | 0 | 0 |
| `x_3` | `x = 3` | 8.7 | 0 | 1.48e-16 |
| `a_1_c_4` | `a = 1`, `c = 4` | 11.3 | 0 | 0 |
| `b_2_d_-3` | `b = 2`, `d = -3` | 9.96 | 0 | 1.48e-16 |
| `a_c_e` | `a = 1`, `c = 4`, `e = -2` | 11.3 | 0 | 4.44e-16 |
| `z_2` | `z = 2` (`offset=at[2].1`) | 8.7 | 0 | 4.44e-16 |
| `lb_10` | `lb = 10` | throws, 100 rejections per seed | same | |

```
runs without partial inits: 41 of 41 parameter vectors bitwise equal
next RNG draw after initialize: 182 of 182 equal
```

The PR v3 run of 2026-11-02 (job 21362578, `results/probe_2026-11-02/`)
had two more cases with a partly supplied tuple (`t1_sigma`,
`w21_sigma`; PR v3: 0, develop: 1.05e21) and 222 runs.  The PR now rejects
these inits (section 7).

On develop, every case also includes the missing array of tuples `at`.
`random_var_context` puts its values into the wrong slots (section 9),
so on develop `at`, and `z` through `offset=at[2].1`, move away from
their draws in every case with `init_radius = 2`.  (With
`init_radius = 0` all draws are 0, and nothing moves.)  This gives the
8.7 of `x_3`.

The 2026-11-01 run at PR v1 (job 21342688, `results/compare.txt`) used
the smaller model without `e`, `at`, `z`, `tv` and `w`.  Section 2 uses
it.

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
`sigma_1e-20`: unconstrained values up to 1.05e21 on develop in job
21365693 and 7.67e20 in job 21342688, at most 46.1 = |log(1e-20)| with
the PR).

## 2. Develop already handles the reverse dependency

> Develop already handles the reverse dependency. A supplied parameter
> whose transform uses missing parameters gets the correct unconstrained
> value, because develop frees it with the same random constrained values
> that the missing parameters start from.

Probe case `x_3` on develop in the 2026-11-01 run (job 21342688, model
without `at`): the missing entries differ from the random draws by at
most 6.66e-16 (relative 2.84e-15), and `x` comes back within 1.48e-16.

## 3. Missing and supplied parameters with the PR

> Missing parameters start exactly at the random unconstrained draws (at
> 0 with `init_radius = 0`), also when their offset or multiplier uses a
> supplied parameter. Supplied parameters keep their values on the
> constrained scale, up to the rounding of `transform_inits` followed by
> `write_array`, as on develop.

Probe (job 21365693), PR build: `max_abs_missing` is 0 in every case;
`mu_5_r0` gives 0 (develop: 5).  `max_rel_supplied` is at most 7.52e-16 in both
builds.

## 4. Runs without partial inits and the RNG stream

> Runs without supplied values or with all values supplied do not change,
> and the number of random draws per attempt is the same as on develop, so
> the RNG stream does not change.

Probe (job 21365693): the 41 parameter vectors of `none_r2` (20 seeds),
`none_r0` and `full_r2` (20 seeds) are bitwise equal between the builds.
The next RNG draw after `initialize` is equal in all 182 runs except the
3 `lb_10` runs, which throw.

## 5. The fixed-point iteration needs more than one pass

> More than one pass is needed when a supplied parameter depends on a
> missing parameter that depends on a supplied one.

Test `chain_provided_a_c_e` and probe case `a_c_e`: `a`, `c` and `e`
are supplied, `b` and `d` are missing, and each has `offset=` the
previous one.  The first pass corrects `c`, the second pass corrects
`e`, and a third pass confirms.  Develop gives missing entries up to 11.3
from their random draws; the PR gives 0.  Test `chain_provided_b_d` and
probe case `b_2_d_-3` need one correction pass (develop: 9.96).

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
> | this PR | 14/14 pass | 16/16 pass |
> | develop | 1/14 pass | 16/16 pass |

Job 21364915 (`results/triton_tests_2026-11-03/smtest_21364915.log` and
the `.run` files): the job ran 36 suites with the PR header, and the two
`initialize` suites also with the develop header and with the PR header
again.  On develop, all `initialize_partial_test` cases fail except
`bound_on_provided_value_is_rejected`.  `partial_tuple_is_rejected`
fails because develop accepts an incomplete tuple.  `provided_dependent_parameter`
fails because its check of the missing entries uses exact equality and
develop differs in the last bits (section 2), and because `at` and `z`
move on develop (section 9).  With the PR header,
`fail_init_test` and the 33 suites in `src/test/unit/services/sample`
also pass (one `.run` file each).  The job ran on the uncommitted files
of `0f517b597` (same md5 as in `results/lint_2026-11-03.txt`).

Earlier versions: job 21361639 (tag `pi_t5`) shows that
`offset_from_array_of_tuples` fails with PR v1, and job 21362559 (tag
`pi_t10`) shows that the four `partially_provided_*` tests fail with PR
v2 and the other 12 pass.  Job 21362615 (`results/triton_tests_2026-11-02/`)
is the PR v3 run with 16 tests.  The 2026-11-01 run of PR v1 (job 21321656,
7 tests) is in `results/triton_tests/`.  `swap/swap_test.sh` runs the
new tests with the develop header in a local Stan clone.

Lint (`results/lint.txt`, 2026-11-01): clang-format 14 and cpplint 2.0.2
with the filters of `make/cpplint` report nothing on the two C++ files.
2026-11-02 (`results/lint_2026-11-02.txt`, with the md5 sums of the
files): cpplint 2.0.2 on Triton (same filters) reports nothing on
`initialize.hpp` and `initialize_partial_test.cpp`; clang-format 14 on
the laptop reports nothing on these and `initialize_test.cpp`.
2026-11-03 (`results/lint_2026-11-03.txt`): cpplint (Triton) and
clang-format 14 (laptop) report nothing on the two files of `0f517b597`.

The third-review commit `2324dfd29` changes `partial_init_unconstrained`
without a change of behavior.  Job 21382891
(`results/triton_tests_2026-11-03_final/`) ran `initialize_partial_test`
and `initialize_test` with the PR and the develop header: 14/14 and
16/16 with the PR, 1/14 with develop.  Job 21382947 ran cpplint and the
header test of `initialize.hpp`.  `results/lint_2026-11-03_final.txt`
has the md5 of `initialize.hpp` and the clang-format output.  `fail_init_test`, the 33
`services/sample` suites and the probe last ran at `0f517b597`.  The
logs of the squashed commits are in `results/triton_tests_2026-11-03_review3/`
and `results/lint_2026-11-03_review3.txt`.

## 8. The reprex of #3088 with the PR, end to end

Not quoted in the PR text (added 2026-11-01).  `cmdstan/run_3088.R`,
step 1, with CmdStan `3ab132c` built with Stan `78b495834` (the PR head)
and Stan Math `dc66941969`, through `tools/cmdstan_run.sh` of the
workspace (job 21345228, log `results/triton_cmdstan_pr/csrun_21345228.log`).
The same 20 seeds as section 1:

| case | develop (section 1) | PR |
|---|---|---|
| `multiplier=sigma`, `sigma = 1e-20`: largest unconstrained `x` | 5.53e20 | 1.82 |
| same: `lp__` at init | -9.32e37 to -1.53e41 | -47.7 to -46.1 |
| same: stuck chains in sampling | 20/20 | 0/20 |
| `offset=mu`, `mu = 5`: largest unconstrained `x` | 28 | 1.8 |
| `lower=a`, `a = 10`: initialization | fails, 100 rejections per seed | same |

With the PR, the built-in `multiplier=sigma` model and the manual model
of the first post (`r1_manual`) have the same smallest and largest
unconstrained values and `lp__` at init over the 20 seeds (the log
prints only these summaries).

## 9. random_var_context and arrays of tuples (develop defect)

`random_var_context_tuples/` holds the probe of stan-dev/stan#3459
about `random_var_context` (found in WardBrian's review of #3457).  The
model
`rvc_tuples.stan` has a missing `array[2] tuple(real<lower=10>,
real<upper=0>) at`.  `probe_rvc_tuples.cpp` prints the `write_array`
order and `vals_r`, and runs `initialize` without inits and with the
partial init `a = 0.5`.  Job 21362046 (tag `rvc_p1`), built with
`tools/stan_probe.sh` against develop `a8ec3d18f` and against the
working-tree `initialize.hpp` that was later committed as PR v2
`7a63912c4` (the job log records head `78b49583` with uncommitted
edits; the agent checked md5 `fa1043d7` on Triton before the job, the
md5 of the file in `7a63912c4`):

| run | develop | #3457 |
|---|---|---|
| `vals_r("at.1")` | at[1].1, at[1].2 (wrong) | same |
| no inits, radius 2 | OK | OK |
| `a = 0.5`, radius 0 | fails, 1 rejection | same |
| `a = 0.5`, radius 2, seeds 1-3 | fails, 100 rejections | same |

The outputs are in `random_var_context_tuples/results/`
(`out_base.txt` = develop, `out_branch.txt` = #3457; identical).

## Other data in this directory

`cmdstan/` step 2 and `grad/` give the numbers of the #3088 comment and
of the documentation PR (stuck chains of the offset/multiplier form with
default inits; gradient errors at small scales).  Their logs are in
`results/triton_cmdstan/` (jobs 21320542, 21320544, 21321652) and
`results/grad_check_out.txt`.  `grad_check.cpp` builds like the probe,
against a Stan Math tree, with `-I <math> -I <math>/lib/eigen_5.0.1
-I <math>/lib/boost_1.87.0 -D_REENTRANT` and the TBB library.
