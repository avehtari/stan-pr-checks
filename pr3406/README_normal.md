# std_normal_lcdf, Phi and the normal cdfs (claims N1 to N6)

This fragment covers the claims about `std_normal_lcdf`, `Phi`,
`std_normal_cdf`, `normal_cdf` and `lognormal_cdf` in the comments on
stan-dev/math#3406.  All commands run from `pr3406/`.  Outputs go to
`out/normal/`.

Last run: 2026-10-30, laptop, g++ 11.4.0 `-O3`, Python 3.10.12,
mpmath 1.3.0.

## Trees

Each tree holds `stan/` (and `test/` where a script reads the tests) of
one commit, plus `lib/` of a clone:

```sh
G=../stan_math              # clone of stan-dev/math with the PR branch
for c in e4745c87ff a28c1c1637 686fc67261 5f33a10d31 bf5861ab9c \
         0d4ac36e34 9b38ecbcb0 cf177ef96b 1f8b6ba94d 1e46958b98; do
  mkdir -p ../trees/$c
  git -C $G archive $c stan test/unit | tar -x -C ../trees/$c
  ln -s $(cd $G && pwd)/lib ../trees/$c/lib
done
export TBB_LIB=$G/lib/tbb OPT=-O3
mkdir -p out/normal
```

| commit | role |
|---|---|
| `e4745c87ff` | merge base, "develop" |
| `a28c1c1637` | the PR as opened (2026-09-18) |
| `686fc67261` | PR head before `5f33a10d31` |
| `5f33a10d31` | erfcx helpers reused in `std_normal_lcdf` |
| `bf5861ab9c` | `Phi` as `0.5 erfc` |
| `0d4ac36e34` | the three cdfs from `erfc` |
| `9b38ecbcb0` | `Phi` and the cdfs from `erfcx` |
| `cf177ef96b`, `1f8b6ba94d` | commits of 2026-10-21 and 2026-10-22 |
| `1e46958b98` | head: PR head `150f54df16` merged with develop |

A helper for the builds below:

```sh
b() { STAN_MATH=../trees/$1 ./build.sh $2 out/normal/$3_$1; }
```

## N1. Kernel accuracy and develop's gradient error (comment of 2026-09-18)

Claims:

> Worst ulp over 601 points per range: [table of nine ranges]

> For comparison, develop's worst in-range relative error is 7.6e-06 to
> 5.6e-05, which is 3.4e+10 to 2.5e+11 ulp, across every branch except the
> far Mills tail.

Commands:

```sh
for c in e4745c87ff a28c1c1637 1e46958b98; do
  b $c std_normal_lcdf/probe_lcdf_grid.cpp lcdf_grid
  out/normal/lcdf_grid_$c > out/normal/lcdf_grid_$c.txt
done
python3 std_normal_lcdf/score_lcdf_grid.py out/normal/lcdf_grid_{e4745c87ff,a28c1c1637,1e46958b98}.txt \
  > out/normal/score_lcdf_grid.txt
```

The probe calls the public `std_normal_lcdf` (value with `double`,
gradient with `var`) and `normal_lcdf(y | 0, 1)` (gradient with `var`).
Grid: 2001 points with both ends in each of the nine ranges of the
comment, and in each of the twelve gradient branches of develop
(`d_*`, given in `s = y / sqrt(2)` in the develop source); plus the 15
points of `expect_branch_accuracy`.  Reference: three routes, `log(erfc(-y
/ sqrt 2) / 2)`, `log(ncdf(y))` (mpmath calls `erfc` for `y < 0`, so this
route is not independent) and quadrature of the density on the small
side.  The C route agrees with A to 3.7e-37 or better at every point.
Run time about 3 min on 12 cores.

Units: `ulp` = relative error / 2^-53; `ulp52` = relative error / 2^-52.
The 2026-09-18 comment uses 2^-52 for the develop figures
(5.62e-05 / 2.22e-16 = 2.5e+11).

Expected output, kernel table (value ulp / gradient ulp):

| range | comment | `a28c1c1637` | head |
|---|---|---|---|
| [-300, -40] | 1.69 / 0.85 | 3.12 / 1.32 | 2.65 / 5.12 |
| [-40, -10] | 1.51 / 0.87 | 2.7 / 1.4 | 2.56 / 4.64 |
| [-10, -5.66] | 1.60 / 0.85 | 2.7 / 1.43 | 2.14 / 4.69 |
| [-5.8, -5.5] | 1.54 / 3.12 | 2.96 / 4.55 | 2.39 / 6.3 |
| [-5.66, -1] | 1.71 / 5.92 | 2.33 / 7.15 | 2.1 / 6.91 |
| [-0.7, -0.6] | 2.37 / 2.75 | 3.32 / 4.41 | 3.28 / 5.12 |
| [-0.6, 0.6] | 1.86 / 2.36 | 2.61 / 3.02 | 6 / 2.79 |
| [0.6, 4] | 7.25 / 4.81 | 8.81 / 6.14 | 7.75 / 6.14 |
| [4, 10] | 30.37 / 31.53 | 33.4 / 33.2 | 33.8 / 33.2 |

The comment's table came from a Python transcription of the kernel (the
comment says so), on 601 points.  The C++ kernel on 2001 points gives
higher values in every row.  The head changes the lower-tail gradient
(1.3 to 5.1 ulp, from `cf177ef96b`, see N5) and the central value (from
`5f33a10d31`, see N2).

Expected output, develop (`e4745c87ff`), worst gradient relative error per
branch:

| branch (s) | `std_normal_lcdf` | `normal_lcdf(y \| 0, 1)` |
|---|---|---|
| A&S, s > 2.9 | 5.91e-05 @ 4.10341 | 5.91e-05 @ 4.10341 |
| Taylor 2.7 | 3.26e-05 @ 4.10122 | 3.26e-05 @ 4.10122 |
| Taylor 2.3 | 2.78e-05 @ 3.53525 | 2.8e-05 @ 3.53553 |
| Taylor 1.85 | 2.27e-05 @ 2.96942 | 6.09e-05 @ 2.12132 |
| Taylor 1.15 | 6.04e-05 @ 2.12083 | 6.09e-05 @ 2.12132 |
| Taylor 0.45 | 7.56e-06 @ 1.13088 | 7.61e-06 @ 1.13137 |
| A&S, no fit | 4.6e-05 @ -2.96829 | 4.63e-05 @ -2.96985 |
| A&S, fit -2.1 | 3.29e-05 @ -5.51416 | 4.63e-05 @ -2.96985 |
| A&S, fit -3.9 | 2.76e-05 @ -9.89949 | 1.49e-05 @ -6.0744 |
| A&S, fit -7 | 3.65e-05 @ -12.3885 | 3.65e-05 @ -12.3885 |
| A&S, fit -17 | 2.36e-05 @ -26.0442 | 2.36e-05 @ -26.0442 |
| Mills, s < -29 | 1.31e-11 @ -41.0122 | 6.16e-06 @ -41.0122 |

- The range over the branches except Mills is 7.56e-06 to 6.04e-05
  (`std_normal_lcdf`) and 7.61e-06 to 6.09e-05 (`normal_lcdf`).  In
  ulp52: 3.4e+10 to 2.72e+11.  The upper end of the comment, 5.6e-05, is
  low: the true maximum is at the end `s = 1.5` of the Taylor-1.15
  branch.
- At the test point `y = 2.121320343559643` (`1.5 sqrt 2`),
  `normal_lcdf` gives 6.09e-05 (the bot's number) and `std_normal_lcdf`
  gives 1.25e-05.  `normal_lcdf` scales by `1 / (sigma sqrt 2)`, which
  gives `s <= 1.5` (Taylor 1.15); `std_normal_lcdf` multiplies by
  `INV_SQRT_TWO`, which gives `s > 1.5` (Taylor 1.85).
- In the head, the worst gradient error is 3.69e-15 (33.2 ulp, at
  `y = 8.281`); below `y = 4` it is at most 7.72e-16 (6.96 ulp).  All 15
  test points are at most 5.76e-16.

## N2. Before and after `5f33a10d31`

Claims (comment on `5f33a10d31`): the table "2001 points per range, value
/ gradient", the central range 6.43 and 8.77 ulp, and "9.6 ulp on a
denser grid of 20001 points".

Commands:

```sh
for c in 686fc67261 5f33a10d31 1e46958b98; do
  b $c std_normal_lcdf/probe_dedup.cpp dedup
  out/normal/dedup_$c > out/normal/dedup_$c.txt
done
python3 std_normal_lcdf/score_dedup.py out/normal/dedup_686fc67261.txt out/normal/dedup_5f33a10d31.txt
python3 std_normal_lcdf/score_dedup.py out/normal/dedup_5f33a10d31.txt out/normal/dedup_1e46958b98.txt
for c in 5f33a10d31 1e46958b98; do
  b $c std_normal_lcdf/probe_central_order.cpp central_order
  out/normal/central_order_$c > out/normal/central_order_$c.txt
  python3 std_normal_lcdf/score_central_order.py out/normal/central_order_$c.txt
  out/normal/dedup_$c 20001 > out/normal/dedup20001_$c.txt
  python3 std_normal_lcdf/score_central_order.py out/normal/dedup20001_$c.txt --lo -0.66 --hi 0.66
done
```

The probe outputs of `686fc67261` and `5f33a10d31` are bitwise equal to
the outputs of 2026-10-19 (`out_pr.txt`, `out_commit.txt`).

Expected output (value / gradient ulp; "same" counts bitwise equal points):

| range | posted before | posted after | rerun `686fc67261` | rerun `5f33a10d31` | head |
|---|---|---|---|---|---|
| [-300, -40] | 2.73 / 1.33 | 2.73 / 1.33 | 2.73 / 1.33 | 2.73 / 1.33 | 2.61 / 5.51 |
| [-40, -5.66] | 2.93 / 1.38 | 2.93 / 1.38 | 2.93 / 1.38 | 2.93 / 1.38 | 2.61 / 5.31 |
| [-5.66, -0.66] | 2.46 / 7.01 | 2.46 / 6.50 | 2.46 / 7.01 | 2.46 / 6.50 | 2.46 / 6.50 |
| [-0.66, 0.66] | 2.93 / 4.01 | 6.43 / 2.74 | 2.93 / 4.01 | 6.43 / 2.74 | 6.43 / 2.74 |
| [0.62, 0.70] | 4.58 / 3.66 | 8.77 / 3.14 | 4.58 / 3.66 | 8.77 / 3.14 | 8.77 / 3.14 |
| [0.66, 5.66] | 11.43 / 10.56 | 11.98 / 10.56 | 11.43 / 10.56 | 11.98 / 10.56 | 11.98 / 10.56 |
| [5.66, 37] | 512.84 / 512.48 | 509.62 / 512.48 | 512.84 / 512.48 | 509.62 / 512.48 | 511.18 / 512.48 |

Denser grids, worst value ulp for `z >= 0`:

| grid | `5f33a10d31` | head |
|---|---|---|
| `probe_central_order`, 20001 points, 3 summation orders | 9.56 @ 0.646366292 (all 3) | 9.56 |
| kernel, `probe_dedup 20001`, points in [-0.66, 0.66] | 9.83 @ 0.646838267 | 9.83 |

- The 9.6 ulp of the comment is the `probe_central_order` result (a copy
  of the formula with `erfcx_small`, not the kernel).  The kernel itself
  reaches 9.83 ulp on the denser kernel grid.
- Routes `erfc` and `ncdf` agree to 5.3e-51 on these grids.
- The head differs from `5f33a10d31` in the lower tails (gradient 1.33
  to 5.51 ulp) and above 5.5; see N5 for the cause.

## N3. `Phi` over 47001 points

Claims (comments of 2026-10-19 and 2026-10-20): the `Phi` table on
[-38.5, 8.5] for develop, the PR before, `bf5861ab9c` and `9b38ecbcb0`;
"its `0.5 * (1 + erf(x / sqrt(2)))` branch on [-5, -0.67] loses up to
7.7e5 ulp".

Commands:

```sh
for c in e4745c87ff a28c1c1637 bf5861ab9c 9b38ecbcb0 1e46958b98; do
  b $c phi/probe_phi_tree.cpp phi_tree
  out/normal/phi_tree_$c > out/normal/phi_tree_$c.txt
done
O=out/normal
paste -d' ' $O/phi_tree_e4745c87ff.txt <(cut -d' ' -f2 $O/phi_tree_a28c1c1637.txt) \
  <(cut -d' ' -f2 $O/phi_tree_bf5861ab9c.txt) <(cut -d' ' -f2 $O/phi_tree_9b38ecbcb0.txt) \
  <(cut -d' ' -f2 $O/phi_tree_1e46958b98.txt) > $O/phi_trees_all.txt
N=develop_e4745c87ff,pr_a28c1c1637,bf5861ab9c,9b38ecbcb0,head_1e46958b98
python3 phi/score_phi.py $O/phi_trees_all.txt --names $N --table --normal-only
python3 phi/score_phi.py $O/phi_trees_all.txt --names $N --worst-in -5 -0.67
# the variants of phi_variants.hpp, and their bitwise check against the trees
b 5f33a10d31 phi/probe_phi_variants.cpp phi_variants
out/normal/phi_variants_5f33a10d31 > $O/phi_variants_5f33a10d31.txt
python3 phi/score_phi.py $O/phi_variants_5f33a10d31.txt --names dev,pr,d1,d2,e,e2 --table --normal-only
python3 phi/compare_variants.py $O/phi_variants_5f33a10d31.txt \
  e4745c87ff=$O/phi_tree_e4745c87ff.txt a28c1c1637=$O/phi_tree_a28c1c1637.txt \
  bf5861ab9c=$O/phi_tree_bf5861ab9c.txt 9b38ecbcb0=$O/phi_tree_9b38ecbcb0.txt \
  1e46958b98=$O/phi_tree_1e46958b98.txt
```

`--normal-only` skips points where `Phi(x) < 2^-1022`, as the posted
tables did.  Without it the row [-38.5, -30) is 9.01e+15 in every column
(subnormal results).  Reference: `erfc`; quadrature of the density agrees
to 1.9e-48 at every worst point.  The variant output is bitwise equal to
`out_phi_variants_e2.txt` of 2026-10-20.

Expected output, worst ulp:

| x | develop | posted "PR before" | `pr` variant | `a28c1c1637` | `bf5861ab9c` | `9b38ecbcb0` | head |
|---|---|---|---|---|---|---|---|
| [-38.5, -30) | 9.01e+15 (cut off) | 1822 | 1.82e+03 | 1.82e+03 | 1.7e+03 | 515 | 514 |
| [-30, -20) | 1.22e+03 | 922 | 922 | 922 | 1.22e+03 | 258 | 258 |
| [-20, -10) | 470 | 440 | 440 | 440 | 470 | 130 | 130 |
| [-10, -5) | 118 | 93 | 92.8 | 92.8 | 118 | 34 | 35.4 |
| [-5, -1) | 7.66e+05 | 19.1 | 19.1 | 18.8 | 30.5 | 12.5 | 12.5 |
| [-1, -0.6) | 3.21 | 4.6 (row [-1, 0)) | 4.57 | 3.83 | 3.05 | 3.9 | 3.9 |
| [-0.6, 0) | 1.97 | | 3.38 | 3.66 | 1.59 | 2.26 | 2.26 |
| [0, 1) | 1.3 | 2.6 | 2.64 | 1.76 | 1.3 | 1.41 | 1.41 |
| [1, 5) | 1.18 | 1.3 (row [1, 8.5)) | 1.32 | 1.07 | 1.18 | 1.32 | 1.32 |
| [5, 8.5) | 0.749 | | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |

- Develop on [-5, -0.67]: 7.66e+05 ulp at `x = -4.988500030550561`.
- The posted "PR before" column is the `pr` variant: `exp` of the kernel
  of `5f33a10d31`, the PR head at that time, not of `a28c1c1637`.
- Bitwise check (points of 47001 that differ): `dev` = `e4745c87ff` (0),
  `d1` = `bf5861ab9c` (0), `e2` = `9b38ecbcb0` (0).  `pr` differs from
  `a28c1c1637` at 1518 points.  The head differs from `e2` at 14126
  points; its worst errors stay within 1.4 ulp of `9b38ecbcb0`.

## N4. The three cdfs (comment on `0d4ac36e34`)

Claim: "Worst relative error in ulp against mpmath (100 digits; partials
by `mpmath.diff`), develop / original PR / after 0d4ac36e34" for `z` from
-38 to 8.

Commands:

```sh
for c in e4745c87ff a28c1c1637 5f33a10d31 0d4ac36e34 1e46958b98; do
  b $c cdfs/probe_cdfs.cpp cdfs
  out/normal/cdfs_$c > out/normal/cdfs_$c.txt
done
python3 cdfs/score_cdfs.py out/normal/cdfs_{e4745c87ff,a28c1c1637,5f33a10d31,0d4ac36e34,1e46958b98}.txt
```

Bitwise: `e4745c87ff` = `out_cdfs_dev.txt`, `5f33a10d31` =
`out_cdfs_new.txt`, `0d4ac36e34` = `out_cdfs_direct.txt` (all of
2026-10-19).  `a28c1c1637` differs from `out_cdfs_new.txt`.  `d/dy`
by `mpmath.diff` agrees with the analytic density to 8.2e-91.  Run time
about 3 min.

Expected output, worst ulp, columns develop / `a28c1c1637` / `5f33a10d31`
/ `0d4ac36e34` / head.  Every posted number equals the rerun of
`e4745c87ff` / `5f33a10d31` / `0d4ac36e34`:

| function | z | value | `d/dy` |
|---|---|---|---|
| `std_normal_cdf` | [-38, -30) | 1.5e+03 / 9.5e+02 / 9.5e+02 / 1.5e+03 / 3.2 | 1e+03 / 9.5e+02 / 9.5e+02 / 1.7 / 2.4 |
| `std_normal_cdf` | [-10, -3) | 3.5e+05 / 51 / 51 / 89 / 3.9 | 65 / 51 / 51 / 1.5 / 2.5 |
| `std_normal_cdf` | [-3, 3) | 60 / 6.3 / 6.3 / 9.8 / 4 | 6.3 / 3.3 / 3.3 / 3.2 / 1.8 |
| `normal_cdf` | [-30, -10) | 1.7e+03 / 1.2e+03 / 1.2e+03 / 1.1e+03 / 9.9e+02 | 1.7e+03 / 1.2e+03 / 1.2e+03 / 9.9e+02 / 9.9e+02 |
| `normal_cdf` | [-10, -3) | 3.5e+05 / 1.2e+02 / 1.2e+02 / 1.6e+02 / 89 | 63 / 1.2e+02 / 1.2e+02 / 87 / 87 |
| `normal_cdf` | [-3, 3) | 59 / 7.5 / 5 / 11 / 5.9 | 10 / 7.6 / 4.6 / 5.9 / 4.6 |
| `lognormal_cdf` | [-30, -10) | 1.7e+03 / 1.2e+03 / 1.2e+03 / 1.1e+03 / 9.9e+02 | 9e+15 / 1.2e+03 / 1.2e+03 / 9.9e+02 / 9.9e+02 |
| `lognormal_cdf` | [-10, -3) | 62 / 1.2e+02 / 1.2e+02 / 1.5e+02 / 88 | 64 / 1.1e+02 / 1.1e+02 / 85 / 85 |
| `lognormal_cdf` | [-3, 3) | 12 / 7.9 / 4.7 / 10 / 4.9 | 9.2 / 6.8 / 3.2 / 4.9 / 4.9 |

- The posted "original PR" column is `5f33a10d31` (the PR head at that
  time).  `a28c1c1637` differs in [-3, 3) only.
- 9e+15 for develop's `lognormal_cdf` `d/dy` is the "returns 0" of the
  posted table.
- Not in the posted table: `lognormal_cdf` `d/dy` on [-38, -30) is
  2.8e+07 / 2.8e+07 / 2.9e+06 / 3.1e+07 (PR trees and head); `d/dsigma`
  on [-3, 3) is 1.9e+15 in all trees (zero crossing, as posted).
  No tree gives a non-finite result for the three cdfs.

## N5. erfcx and std_normal_lcdf within about 2.7 x 2^-52

Claim (comment after `1f8b6ba94d`): "Accuracy against mpmath stays within
about 2.7 x 2^-52 for `erfcx` and `std_normal_lcdf`."

Commands:

```sh
for c in cf177ef96b 1f8b6ba94d 1e46958b98; do
  b $c std_normal_lcdf/probe_steve.cpp steve
  out/normal/steve_$c > out/normal/out_probe_$c.txt
done
(cd out/normal && python3 ../../std_normal_lcdf/score_steve.py \
  out_probe_cf177ef96b.txt out_probe_1f8b6ba94d.txt out_probe_1e46958b98.txt)
```

Unit: relative error / 2^-52.  The `cf177ef96b` output is bitwise equal
to the output of 2026-10-22.  Expected output (worst / median / count
non-finite), equal in all three trees:

```
erfcx                        1.32 /   0.313 / nf 0
slope -5.66<z<0              2.74 /   0.617 / nf 0
slope z<=-5.66               2.08 /   0.709 / nf 0
value -5.66<z<0             0.722 /   0.323 / nf 0
value z<=-5.66               1.08 /   0.215 / nf 146
value z>=0                4.5e+15 /      20 / nf 0
slope z>=0                4.5e+15 /      20.6 / nf 0
```

- The bound holds for `erfcx` on [4, 1e300] and for `std_normal_lcdf` at
  `z < 0`: worst 2.74.
- It does not hold for `z >= 0` (worst 4.5e+15, median 20): there the
  result rounds to 0 or loses digits where `log Phi(z)` is near 0.  The
  2026-10-22 note excluded these rows.  N1 gives the `z >= 0` errors per
  range (up to 33.8 ulp on [4, 10]).
- The 146 non-finite values are at `|z| > 1.3e154`, where `-z^2/2`
  overflows.

Before `cf177ef96b` (tree `e8bb9a0b38`, same commands, 2026-10-30), for
comparison:

```
section                             e8bb9a0b38.txt              1e46958b98.txt
erfcx                       0.922 /   0.253 / nf 0       1.32 /   0.313 / nf 0
slope -5.66<z<0              2.29 /   0.671 / nf 0       2.74 /   0.617 / nf 0
slope z<=-5.66              0.619 / 4.9e-46 / nf 0       2.08 /   0.709 / nf 0
value -5.66<z<0             0.722 /   0.323 / nf 0      0.722 /   0.323 / nf 0
value z<=-5.66             1.19 /   0.221 / nf 146     1.08 /   0.215 / nf 146
```

`cf177ef96b` changed the public `erfcx` for `x >= 4` (one division over a
common denominator): worst error 0.92 to 1.32 (units of 2^-52) on [4,
1e300].

## N6. Reference constants of the tests

Constants: `std_normal_lcdf_impl_test.cpp` (`scalar_tail_kernel`, 20
rows, value and slope), `Phi_test.cpp` (`TEST(MathFunctions, Phi)`), and
the 15 rows of `expect_branch_accuracy`.

Commands:

```sh
python3 refs_normal/gen_refs.py ../trees/1e46958b98
python3 refs_normal/gen_refs.py ../stan_math     # working tree with the new Phi row
```

Route A: `erfc` at 60 digits.  Route B: quadrature of the density on the
small side at 45 digits (asymptotic series of the Mills ratio for
`|z| > 1e6`).  Columns: the constant, route A, `|B/A - 1|`,
`|constant/A - 1|`.  Run time about 10 s.

Expected output:

- Routes agree to 3.2e-46 or better at all 60 rows (59 rows at
  `1e46958b98`).
- Every constant is within 8.5e-17 relative of route A (the reference
  rounded to double), except the rows `z = 1e8` and `z = 1e150` (value
  and slope).  Their constant is 0: the true values (about
  `-exp(-5e15)` and smaller) are below the smallest double.
- `Phi(-4.988500030550561)`: constant 3.0424950967132873e-07, route A
  3.0424950967132873e-07, routes 1.3e-46, constant 3.0e-17.  This row is
  in the working tree of the clone, not in `1e46958b98`.
- `slope(-1e150)`: the constant 1e+150 equals route A, and route B (the
  asymptotic series) gives a relative difference of 0.0e+00.
