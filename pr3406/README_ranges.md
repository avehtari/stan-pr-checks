# stan-dev/math#3406: ranges, wiener4 log Mills ratio, varmat helpers

This file holds three sections for the main `README.md`.  All commands
run from `pr3406/`.

Trees: `DEV` is the merge base `e4745c87ff` (the "develop" of the
review).  `HEAD` is `1e46958b98`, the PR head `150f54df16` merged with
develop `dc66941969`.  Each tree needs `stan/` and `lib/`.

Last run: 2026-10-30, laptop (12 cores, Ubuntu 22.04, g++ 11.4.0,
Python 3.10.12, mpmath 1.3.0).

```sh
export DEV=<tree of e4745c87ff> HEAD=<tree of 1e46958b98>
export TBB_LIB=<a lib/tbb with libtbb.so>
mkdir -p out/ranges out/wiener4 out/varmat
```

## Valid input ranges

Claim: "The valid ranges of inputs (before over/underflow) across all
functions are expanded" (PR description, 12-row table).

```sh
STAN_MATH=$DEV  OPT=-O2 ./build.sh ranges/probe_ranges.cpp out/ranges/probe_dev
STAN_MATH=$HEAD OPT=-O2 ./build.sh ranges/probe_ranges.cpp out/ranges/probe_head
python3 -u ranges/probe_ranges.py --dev out/ranges/probe_dev \
  --head out/ranges/probe_head > out/ranges/ranges.txt     # about 15 min
```

Parameters: `mu = 0`, `sigma = 1`; `exp_mod_normal` `lambda = 1`;
`skew_normal` `xi = 0`, `omega = 1`, `alpha = 1`.  The axis is `z = y`,
except for `lognormal`, where it is `z = log(y)`.  The gradient is
`d/dy` (a `var` `y`, the other arguments `double`); for `lognormal` it is
also `d/dy`, not `d/dz`.  The value comes from the `double` call.

`ranges/probe_ranges.cpp` is a line server: it reads `<function> <y>`
and prints the value and `d/dy` in hex.  `ranges/probe_ranges.py`
drives both binaries.  `ranges/refs_normal.py` holds the references.

Method, per function, tree, value or gradient, and side:

1. Scan from the center outward: 40 points per decade for
   `|z|` in [0.01, 1000], then 4 per decade up to `DBL_MAX`
   (`lognormal`: `z` steps of 0.05 to 40, then 1, to the ends of `y`).
2. At the first failing scan point, bisect on the bit pattern of the
   double between it and the last passing point (relative width 1e-7).
3. Print the end to 3 significant digits, and the value and reference
   at the first failing point.  `L` marks a representable limit: the
   scan reached `DBL_MAX` (or the smallest `y`), or the exact result is
   above `DBL_MAX`.  If points beyond the first failure pass again, the
   output says "valid again".

A point passes when the result is finite, has a relative error of at
most `tol` with an absolute floor of `2^-1074`, and is not saturated.
Saturated means equal to a limit value of the function (0 for a log
cdf; 0 or 1 for a cdf; 0 for a gradient) where the exact value is
different.

| criterion | tol | saturated result |
|---|---|---|
| T (the table's meaning) | 1e-2 | fails |
| A (accurate) | 1e-8 | fails |
| C (correctly rounded) | 1e-8 | passes when within tol |

References: mpmath at 40 digits (more for the `exp_mod_normal` lower
tail), two routes that differ in `erfcx(x)`.  Route A is
`erfc(x) exp(x^2)`, or the asymptotic series above `x = 1e6`.  Route B
is `Gamma(1/2, x^2) exp(x^2) / sqrt(pi)`, or the Lentz continued
fraction above `x = 25`.  The output prints the route agreement at each
test point.  The largest disagreement on the scan grid was 7e-41
(normal), 5e-50 (`exp_mod_normal`).

Known limit of the scorer: above about `z = 9.5`, `1 - Phi(z)` is below
1e-40, so the 40-digit reference of a cdf is exactly 1.  The cdf rows
then show "valid again" beyond the saturation end.  This does not
change the ends.

Expected output: the three tables at the end of `out/ranges/ranges.txt`,
then the test candidates (a point inside the head range where develop
fails, with both trees and both routes), then the points that the
review names (-1.5e154, -1e154, -1e150, -1.5e308, -2e5, 745, 746).

Criterion T, `*` = differs from the PR cell:

| function | dev value | head value | dev grad | head grad |
|---|---|---|---|---|
| normal_lcdf, std_ | [-9.46e30, 38.5]* | [-1.9e154 L, 38.5]* | [-1.27e308, 38.6] | [-1.8e308 L, 38.6] |
| normal_lccdf, std_ | [-38.5, 9.46e30]* | [-38.5, 1.9e154 L]* | [-38.6, 1.27e308] | [-38.6, 1.8e308 L] |
| normal_cdf, std_ | [-37.5, 8.24] | [-38.5, 8.29] | [-37.5, 38.6] | [-38.4, 38.6] |
| Phi | [-37.5, 8.24] | [-38.5, 8.29] | [-38.6, 38.6] | [-38.6, 38.6] |
| lognormal_lcdf | [-38.5, 7.83]* | [-744 L, 38.5] | [-38.4, 37.6] | [-703 L, 37.6] |
| lognormal_lccdf | [-7.8, 38.5]* | [-38.5, 710 L] | [-38.5, 38.4] | [-38.5, 710 L] |
| lognormal_cdf | [-38.5, 8.29] | [-38.5, 8.29] | [-27.1, 37.6] | [-38.4, 37.6] |
| exp_mod_normal_lpdf | [-37.5, 1.8e308 L] | [-1.9e154 L, 1.8e308 L]* | [-37.4, 1.8e308 L] | [-1.8e308 L, 1.8e308 L] |
| exp_mod_normal_lcdf | [-7.3, 33.9]* | [-1.9e154 L, 746]* | [-6.59, 746]* | [-1.27e308, 746]* |
| exp_mod_normal_lccdf | [-6.64, 746]* | [-38.4, 1.27e308] | [-6.91, 746]* | [-38.5, 1.8e308 L] |
| exp_mod_normal_cdf | [-6.64, 37.9]* | [-38.4, 37.9] | [-6.91, 746]* | [-38.3, 746] |
| skew_normal_lpdf | [-38.5, 1.34e154]* | [-1.34e154 L, 1.34e154]* | [-38.4, 1.8e308 L] | [-8.99e307 L, 1.8e308 L]* |

Criterion A differs from T in these cells: the develop gradient of
`normal_lcdf` and `normal_lccdf` fails 1e-8 already at `|z|` = 0.002
(worst 5.35e-05 at `z` = 2.11, the former Taylor branches); develop
`lognormal_lcdf` value [-38.2, 5.72]; develop `exp_mod_normal` ends move
to about -5.

## wiener4 log Mills ratio

Claim (PR comment): `wiener4_lcdf_unnorm` line 128,
`std_normal_lcdf(-x) - std_normal_lpdf(x)`, loses accuracy without
bound, also with correctly rounded inputs; the form
`log(sqrt(pi/2) * erfcx(x / sqrt(2)))` does not.

```sh
for t in dev head; do
  [ $t = dev ] && S=$DEV || S=$HEAD
  STAN_MATH=$S OPT=-O2 ./build.sh wiener4/probe_mills.cpp out/wiener4/probe_mills_$t
  out/wiener4/probe_mills_$t > out/wiener4/mills_$t.txt
done
python3 wiener4/score_mills.py out/wiener4/mills_dev.txt out/wiener4/mills_head.txt
```

The probe calls the tree's `std_normal_lcdf`, `std_normal_lpdf` and
`erfcx`.  "CR inputs" subtracts the correctly rounded `lcdf(-x)` and
`lpdf(x)` in double.  ulp = relative error / 2^-53; the reference is
mpmath at 60 digits, two routes (`log ncdf - log npdf`, and `erfcx` from
the incomplete gamma), agreement 1e-55 or better.

Expected (head; develop differs only where its lcdf is less accurate):

| x | difference | CR inputs | erfcx form |
|---|---|---|---|
| 10 | 8.0 | 8.0 | 0.65 |
| 30 | 87.0 | 87.0 | 0.03 |
| 100 | 340 | 340 | 1.26 |
| 300 | 8 601 | 8 601 | 0.21 |
| 1e4 | 6.0e6 | 1.3e6 | 0.17 |
| worst, log grid 1 to 1e4 | 6.0e6 | 4.0e6 | 3.41 |

The comment's 6 132 ulp at `x = 300` does not reproduce: this script
gives 8 601 ulp at `x = 300`, with the tree's functions and with
correctly rounded inputs.

## varmat helpers

Claim: "`as_array_or_scalar()` now returns an array unchanged instead of
re-wrapping it; `ref_type_if` keeps an rvalue `arena_matrix` by value
rather than referencing a dead temporary."

```sh
DEV=$DEV HEAD=$HEAD varmat/run_varmat.sh        # 8 ASan builds, ~5 min
varmat/reach.sh dev=$DEV head=$HEAD opened=<tree of a28c1c1637>
```

`varmat/probe_varmat.cpp` has four cases, one per `-DCASE_x`.
`run_varmat.sh` builds each with
`-fsanitize=address -fno-omit-frame-pointer -g -O1` and runs it.
`varmat/reach.sh` copies `stan/` of each tree, puts back the develop
helpers with a `static_assert` that fires on the defect, and compiles
`varmat/reach_probe.cpp` (15 functions; `y` as `var_value<VectorXd>`
and `Matrix<var>`; `mu`, `sigma` as `var` or `var_value<VectorXd>`).

Expected:

| case | develop | head |
|---|---|---|
| `as_array_or_scalar` twice on `var_value<VectorXd>` | compile error | values and adjoints correct |
| `as_array_or_scalar(var_value<ArrayXd>)` | compile error | values and adjoints correct |
| `ref_type_t<arena_matrix<VectorXd>&&>` | `const` reference | value |
| `ref_type_t<arena_matrix<Matrix<var>>&&>` | `const` reference | value |
| `const auto& r = to_ref(<rvalue arena_matrix>)`, read `r` | ASan stack-use-after-scope | no ASan report |
| `reach.sh` (dev, head, a28c1c1637) | compiles | compiles |

`reach.sh` compiles on all three trees, so no `static_assert` fires.
None of the probed functions passes an rvalue `arena_matrix` to
`ref_type_if` or an array `var_value` to `as_array_or_scalar`.  This
holds on develop, at the head with the develop helpers, and at
`a28c1c1637` (the PR as opened).  The PR's two varmat test files
(`mix/prob/normal_ccdf_log_test.cpp`, `std_normal_ccdf_log_test.cpp`)
also compile against the instrumented head.  A grep of develop finds
no `to_ref` or `to_ref_if` call on `to_arena(...)`, `arena_matrix(...)`
or `std::move(...)`.
