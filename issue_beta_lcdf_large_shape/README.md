# beta_lcdf and beta() at large shapes (issue draft)

This directory holds the scripts behind the numbers of an issue draft
about `stan::math::beta(a, b)` and `beta_lcdf` at large shapes (`b` from
1e3 to 2^31).  The sections quote the claims of the draft and give the
commands and the expected output.

Last run: 2026-10-31 on a laptop (`results/versions.txt`: g++ 11.4.0,
`-O1`, Python 3.10.12, mpmath 1.3.0).

## Commits

| name | commit | what |
|---|---|---|
| develop | `dc66941969` | stan-dev/math develop of 2026-10-29 |
| #3433 | `709a1a4d76` | head of stan-dev/math#3433, which is stacked on #3490, #3489 and #3409; "#3433" below means this head |

Export `stan/` of both commits, for example with
`../pr3496/make_trees.sh`, which writes
`out/trees/develop` and `out/trees/pr3433` in that directory.  Set
`STAN_MATH_LIB` to the `lib/` of a stan-dev/math clone and `TBB_LIB` to a
directory with a built `libtbb.so`.

## 1. beta_lcdf at large b

> On `develop` (`dc66941969`), the relative error of the `x` derivative
> of `beta_lcdf` is the error of `beta(a, b)`. ...

`probe_beta_lcdf.cpp` evaluates `beta_lcdf(x | a, b)` and its three
derivatives, and Boost `ibeta(a, b, x)`, for `a` in {0.5, 2, 10}, `b` in
{1e3, 1e5, 1e6, 1e7, 1e8, 1e9, 2^31} and `x = t / b` with `t` in {0.5, 2,
8} (63 points).  `score_beta_lcdf.py` compares them with two mpmath
routes: the log of the lower incomplete beta, and `log1p(-I)` with `I`
the upper incomplete beta.  The routes agree to 2e-32 or better.

```sh
T=../pr3496/out/trees
for t in develop pr3433; do
  STAN_MATH=$T/$t OPT=-O1 ./build.sh probe_beta_lcdf.cpp out/probe_$t
  ./out/probe_$t > out/beta_lcdf_$t.txt
  python3 score_beta_lcdf.py < out/beta_lcdf_$t.txt > out/score_$t.txt
done
```

Expected: `results/score_develop.txt`, `results/score_pr3433.txt`
(relative errors).  The largest errors over the 63 points:

| quantity | develop | #3433 | point of the largest develop error (a, b, b x) |
|---|---|---|---|
| value | 7.1e-8 | 7.1e-8 | 2, 1e9, 2 |
| `d/dx` | 3.0e-6 | 3.0e-6 | 2, 2^31, 0.5 |
| `d/da` | 6.8e-3 | 4.2e-7 | 0.5, 2^31, 8 |
| `d/db` | 1.2e-2 | 7.3e-6 | 0.5, 2^31, 8 |
| Boost `ibeta` value | 3.7e-8 | 3.7e-8 | 2, 1e9, 2 |

The `d/dx` column grows like `b log(b) eps`: about 3e-10 at `b = 1e6`, 1.5e-6 at
`b = 1e9` (`a = 2`).

## 2. The cause of the x-derivative error: beta(a, b)

> `stan::math::beta(a, b)` computes `exp(lgamma(a) + lgamma(b) -
> lgamma(a + b))`. ... `exp(lbeta(a, b))` does not have this problem ...

The `x` derivative on develop is `pow(1 - x, b - 1) * pow(x, a - 1) / P /
beta(a, b)`.  `probe_dx_parts.cpp` prints `beta(a, b)`, `pow(1 - x, b -
1)` and `exp(lbeta(a, b))` at `x = 2 / b`; the inline script below
scores them against mpmath.

```sh
STAN_MATH=$T/develop OPT=-O1 ./build.sh probe_dx_parts.cpp out/probe_dx_parts
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

Expected (`results/dx_parts_score.txt`): the error of `beta(a, b)` is the
error of the `x` derivative (2.98e-6 at `a = 2`, `b = 2^31`).
`exp(lbeta(a, b))` is within 3e-14.  `pow(1 - x, b - 1)` has a relative
error of 5.5e-8 at `b = 1e9`, where `1 - x` is rounded; at `b = 2^31`, `x = 2^-30` and
`1 - x` is exact.

## 3. The shape derivatives with #3433

> At the head of #3433, with the new `grad_reg_inc_beta` of #3409, it is
> 7.3e-6. The remaining error is probably the plain difference
> `digamma(a + b) - digamma(b)` in `grad_reg_inc_beta`. ...

The #3433 `grad_reg_inc_beta` does not use `beta(a, b)`; it uses
`lbeta`.  It multiplies by `digamma(a + b) - digamma(b)`
(`stan/math/prim/fun/grad_reg_inc_beta.hpp`, lines 132 and 137 at
`709a1a4d76`).  At `a = 2`, `b = 1e9` that difference is about 2e-9,
while `digamma(b)` is about 20.7, so the plain difference has a relative
error of about `20.7 eps / 2e-9 = 2.3e-6`.  The probe gives 1.5e-6 for
`d/db` there.  The scripts do not test this cause.

## 4. The x derivative with inc_beta_ddz

> At the head of #3433, `beta_lccdf` computes the `x` derivative with
> `inc_beta_ddz`. At 36 of the points below, its relative error is below
> 1.4e-15 for `a = 0.5` and `a = 10`. On `develop`, it is up to 3.0e-6.
> For `a = 2` and `b x` of 2 or more it reaches 5.4e-8 (at `b = 1e9`),
> because `beta_lccdf` divides by `1 - I`, and Boost `ibeta` has a value
> error there (see below).

`probe_beta_lccdf_dx.cpp` evaluates the `x` derivative of `beta_lccdf` at
36 of the points of section 1.  `score_lccdf_dx.py` compares it with two
mpmath routes: the closed form with the upper incomplete beta, and
`mp.diff` of `log1p(-I)` with the lower incomplete beta.  The routes
agree to 6e-33 or better.

```sh
for t in develop pr3433; do
  STAN_MATH=$T/$t OPT=-O1 ./build.sh probe_beta_lccdf_dx.cpp out/probe_lccdf_dx_$t
  ./out/probe_lccdf_dx_$t > out/beta_lccdf_dx_$t.txt
  python3 score_lccdf_dx.py < out/beta_lccdf_dx_$t.txt > out/score_lccdf_dx_$t.txt
done
```

Expected: `results/score_lccdf_dx_develop.txt` (up to 3.0e-6, the error
of `beta(a, b)`) and `results/score_lccdf_dx_pr3433.txt`: within 1.4e-15 for
`a = 0.5` and `a = 10`.  For `a = 2` and `b x` = 2 and 8: 1.6e-15 and
6.9e-15 at `b = 1e3`, 5.4e-11 and 8.0e-12 at `b = 1e6`, 5.4e-8 and
4.2e-9 at `b = 1e9`, 2.4e-16 and 5.6e-14 at `b = 2^31`.  The
function divides by `1 - I`, so the value error of Boost `ibeta` (section
1) enters at these points.  At `a = 2`, `b = 1e9`, `b x = 2`, `I` is about
`1 - 3 exp(-2) = 0.594`.  The relative `ibeta` error 3.7e-8 then gives
`3.7e-8 * 0.594 / 0.406 = 5.4e-8` in `1 - I`.
