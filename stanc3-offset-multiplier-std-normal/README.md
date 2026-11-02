# stanc3: std_normal for a parameter with offset/multiplier (prototype)

Scripts and results behind the stanc3 feature request "Generate
std_normal for a parameter with offset/multiplier and a matching normal
statement" (from stan-dev/stan#3088).  The prototype is not a PR; the
patch is in `stanc3_prototype.patch`.

## The prototype

`stanc3_prototype.patch` applies to stanc3 `master` `af10bc294`.  It adds
the pass `src/analysis_and_optimization/Offset_multiplier_std_normal.ml`
and the setting `offset_multiplier_std_normal` in `Optimize.ml`/`.mli`.
The pass is off at `--O0` and `--O1`.  It is on at `--Oexperimental` and
with `-foffset-multiplier-std-normal`.  It runs first in
`optimization_suite`, on `log_prob` and `reverse_mode_log_prob` only.
`optimize_ad_levels` also starts from the top-level `FnReadParam`
declarations of `log_prob`, so that `x_raw__` keeps an autodiff type.

The pass applies to a parameter `x` with `<offset=m, multiplier=s>` and a
statement `x ~ normal(m, s)` or `target += normal_lpdf(x | m, s)`.  The
location and scale must be structurally equal to the offset and
multiplier.  The generated `log_prob` is then

    x_raw__ = in__.read<...>(dims);                 // no transform
    x = offset_multiplier_constrain(x_raw__, m, s);
    x_log_scale__ = L;                               // L = sum(log(s)) over the elements of x
    jacobian += x_log_scale__;
    target += std_normal_lpdf<propto__>(x_raw__) - x_log_scale__;   // each matching statement

`x` is `real`, `vector`, `row_vector` or `array[] real`.  `write_array`,
`transform_inits` and `generate_quantities` do not change.  The patch has
no tests.  The only `dune runtest` differences are the help-text lines
that name the new pass.

## Requirements and versions

- stanc3 `af10bc294` with the patch, built with
  `dune build src/stanc/stanc.exe --profile static` (OCaml 5.5.1).  Copy
  `_build/default/src/stanc/stanc.exe` as `stanc_om_static`.
- CmdStan `3ab132c` with Stan `78b495834` (stan-dev/stan#3457; the
  default-init path is the same as `develop` `a8ec3d18f`) and Stan Math
  `dc66941969`.  The jobs did not use its own `bin/stanc` (a stanc3
  nightly, `10361e1`) for the compared models: a wrapper ran
  `stanc_om_static` (§"How to run").
- R with cmdstanr and posterior.

## How to run

Put `run_om_std_normal.R`, `om_probe_main.hpp`, `models/*.stan` and
`stanc_om_static` into one directory `DIR`, set `CMDSTAN` and run

    N_CORES=32 Rscript run_om_std_normal.R DIR [step ...]   # steps 1 2 3 4

"base" is `stanc_om_static` without the flag.  "om" is `stanc_om_static`
with `-foffset-multiplier-std-normal`.  CmdStan's make runs stanc as
`$(WINE) bin/stanc ...`; the script sets `WINE` to a wrapper that runs
`stanc_om_static` instead.

cmdstanr builds its model methods with `bin/stanc` directly, so the script
does not use these methods.  Steps 1 and 2 use the `log_prob` method of the
executables (`sig_figs=18`, `propto` true).  Step 4 links the probe
`om_probe_main.hpp` instead of CmdStan's `main.o`
(`make CMDSTAN_MAIN_O=...`).  The results below are from Triton job
21347618 (32 cores, 32 min), log `results/csrun_21347618.log`.

## Step 4: same target as the current code (`models/om_nodata.stan`)

The model has 43 unconstrained parameters.  The pass rewrites 13 of them.
It does not rewrite the indexed statement `q[i] ~ ...`, the tuple
component, the scale `2 * sigma` or the matrix.  The cases include a loop,
an `if`, a `profile`, truncation, a repeated statement and a user function
as the multiplier.  They also include data, int and vector multipliers.
The script compares each variant with base at `--O0` at 20 random points.

| variant | propto 0: max abs lp diff | propto 1: lp diff | max rel grad diff |
|---|---|---|---|
| om `--O0` | 5.7e-14 | -0.6527 | 8.4e-15 |
| om `--O1` | 5.7e-14 | -0.6527 | 8.4e-15 |
| om `--Oexperimental` | 5.7e-14 | -0.6527 | 8.4e-15 |
| base `--Oexperimental -fno-offset-multiplier-std-normal` | 2.8e-14 | 2.8e-14 | 8.5e-15 |

With `propto` true, the rewrite keeps `-L` for the data multipliers of
`g`, `g2` and `h`, and `normal_lupdf` drops it.  The expected shift is
`-(4 log 0.7 + 3 log 2) = -0.6527` at every point.

## Step 1: `models/om_cases.stan`, `log_prob` method

At 20 random points, om against base, the max rel grad diff is 1.1e-14.
The max abs lp diff is 1.1, from the data multiplier of `f` (expected
constant `-3 log 0.7 = 1.07`).

## Step 2: gradients of `j_om` against `j_ncp`

`j_om` has `vector<offset=mu, multiplier=sigma>[1000] re` and
`re ~ normal(mu, sigma)`; `j_ncp` has `raw ~ std_normal()`.  Absolute
gradient differences against `j_ncp`, `jacobian` true, `raw` uniform on
(-2, 2):

| mu | sigma | base: mu | base: raw (max) | om: mu | om: raw (max) |
|---|---|---|---|---|---|
| 0 | 1e-14 | 2.7 | 6.7e-16 | 1.8e-13 | 0 |
| 0 | 1e-30 | 9e+15 | 4.4e-16 | 2e-13 | 0 |
| 1 | 1e-14 | 6.3 | 0.011 | 2e-11 | 0 |
| 1 | 1e-18 | 1.2e-11 | 2 | 1.2e-11 | 0 |
| 1000 | 1e-14 | 1.1e-08 | 2 | 1.1e-08 | 0 |

`results/csrun_21347618.log` has all 48 rows (sigma 1e-2 to 1e-30).

## Step 3: default inits

20 seeds per model and variant, 1000 warmup iterations, 1000 sampling
iterations (200 for `d_om`).  A run is stuck if the median of the scale
parameter (`sigma`, or `s` in `k_om`) after warmup is below 1e-6.  The
divergence counts are sums over the 20 seeds.

| model | base stuck | om stuck | base divergent (warmup) | om divergent (warmup) |
|---|---|---|---|---|
| `j_om` (normal, 1000 groups) | 2 | 0 | 382 | 409 |
| `k_om` (prior only, n = 200) | 5 | 0 | 3444 | 155 |
| `d_om` (binomial, I = 10 000) | 7 | 0 | 3321 | 477 |

The base counts are the same as the CmdStan `develop` counts in
`../stan-pr3457/` (same seeds).
