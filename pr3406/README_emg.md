# exp_mod_normal claims and bitwise claims of stan-dev/math#3406

This file covers the `exp_mod_normal` lcdf / lccdf / cdf claims (M1 to M5)
and the bitwise-identity claims (B1) of the comments on stan-dev/math#3406.
Each section gives the quoted claim, the commands, the expected output, the
commits, and the result of the rerun.

Last run: 2026-10-30, on a laptop (12 cores, Ubuntu 22.04, g++ 11.4.0,
`-O3`, Python 3.10.12, mpmath 1.3.0).

## Layout

| path | content |
|---|---|
| `emg/probe_emg_phi.cpp` | lcdf / lccdf, value and 4 partials, 170-point grid |
| `emg/probe_emg_dense.cpp` | the same on the dense grid (1250 points) |
| `emg/score_emg_grid.py` | mpmath scorer (uses `score_emg_phi.py`) |
| `emg/check_ref_rows.cpp` | the 13 + 6 test rows with the test comparison |
| `emg/probe_deep_grid.cpp` | non-finite counts on the 105-point deep grid |
| `emg/run_probes.sh` | builds and runs the four probes per tree |
| `emg/ocl/` | device code compiled on the host against the CPU kernel |
| `emg/orig/` | the original probes and outputs of 2026-10-19 to 10-22 |
| `refs_emg/` | the two reference generators, the test header, a comparison |
| `bitwise/` | `compare.cpp` and `run_compare.sh` |

All commands run from `pr3406/`.  Outputs go to `out/emg/` and
`out/bitwise/`.

## Trees

The probes need one Stan Math tree (`stan/`, `test/unit`, `lib/`) per
commit in one directory, `TREES_DIR`:

```sh
export TREES_DIR=$PWD/../trees           # any directory
export TBB_LIB=<a Stan Math clone>/lib/tbb
for c in e4745c87ff 686fc67261 927099efb8 a17bba2e9e 39886e0eb8 cf177ef96b \
         d081fff3a4 1f8b6ba94d 1e46958b98 0d4ac36e34 f4cd10ee4a c30ee6f9ef \
         9b38ecbcb0; do
  mkdir -p $TREES_DIR/$c
  git -C <clone> archive $c stan test/unit | tar -x -C $TREES_DIR/$c
  ln -s <clone>/lib $TREES_DIR/$c/lib
done
```

| commit | role |
|---|---|
| `e4745c87ff` | merge base, "develop" in the review |
| `686fc67261` | PR head of the 2026-10-19 audit |
| `927099efb8` | last commit with the `log_diff_exp` form |
| `a17bba2e9e` | the erfcx-difference form |
| `39886e0eb8` | the OpenCL port of a17bba2e9e |
| `cf177ef96b`, `d081fff3a4` | before the deep-tail fix |
| `1f8b6ba94d` | the deep-tail fix |
| `1e46958b98` | PR head `150f54df16` merged with develop |

Build and run all emg probes (about 6 min with `JOBS=6`):

```sh
JOBS=6 emg/run_probes.sh 1e46958b98 927099efb8 e4745c87ff a17bba2e9e \
  d081fff3a4 1f8b6ba94d cf177ef96b 686fc67261
```

## M1. The previous form fails in the left tail and for small v

Claim: "The previous form, `log F = log_diff_exp(log Phi(z), b)` with
`b = v^2/2 - v z + log Phi(z - v)` and `v = lambda * sigma`, has two
failures at this PR's head (grid of `z` in [-1e4, 100] and `v` in [1e-8,
100], against mpmath)", with this table:

| function | region | value | dy | dsigma | dlambda |
|---|---|---|---|---|---|
| lcdf | left tail (`z <= -30`) | 2.9e5 ulp | 3.4e7 | 3.4e7 | 1.5e13 |
| lcdf | small `v` (`<= 1e-4`) | non-finite (10 points) | | | |
| lccdf | small `v` | 1.4e12 | 1.3e12 | 1.3e12 | 6.0e2 |

The bot asks: "say which reference rows fail with the form before
a17bba2e9e".

### M1a. The grid script behind the table

Grid: 17 `z` from -1e4 to 100 and 10 `v` from 1e-8 to 100, `mu = 0`,
`sigma = 1` (170 points, 340 output lines with both functions).  The
reference is mpmath at 100 digits, the partials by `mpmath.diff`; the
second route for `d/dy` is the analytic derivative.

```sh
python3 emg/score_emg_grid.py out/emg/ref_cache_grid.json \
  out/emg/probe_emg_phi_{927099efb8,686fc67261,e4745c87ff,a17bba2e9e,1f8b6ba94d,1e46958b98}.txt \
  > out/emg/score_grid.txt
```

The reference is cached in the json file; the first run computes it on
all cores.
Expected output for `927099efb8`:

```
lcdf: 170 points; non-finite components 10 at 2 points
regime                    value       dy      dmu   dsigma  dlambda
rest                    9.2e+02  3.0e+03  3.0e+03  3.4e+03  5.5e+04
small v (<= 1e-4)           inf      inf      inf      inf      inf
left tail (z <= -30)    2.9e+05  3.4e+07  3.4e+07  3.4e+07  1.5e+13
large v (>= 10)         7.5e+02  7.5e+02  7.5e+02  4.4e+03  1.1e+04
lccdf: 170 points; non-finite components 0 at 0 points
small v (<= 1e-4)       1.4e+12  1.3e+12  1.3e+12  1.3e+12  6.0e+02
```

Rerun result:

- `927099efb8` and `686fc67261` reproduce every number of the table.
  The output of 2026-10-19 (`emg/orig/out_pr.txt`) gives the same table.
- "non-finite (10 points)" is 10 non-finite components at 2 points,
  `(z, v) = (-1000, 1e-8)` and `(-1e4, 1e-4)`, all 5 components each.
- `e4745c87ff` (develop): lcdf 350 non-finite components at 70 points,
  lccdf 145 at 29 points.
- `1e46958b98`: no non-finite result; see M2.

### M1b. The reference rows with the form before a17bba2e9e

`emg/check_ref_rows.cpp` runs the 13 rows of
`exp_mod_normal_lcdf_reference` and the 6 rows of
`exp_mod_normal_lcdf_deep_tail_reference` (header copied to `refs_emg/`
from `1e46958b98`) with the checks of `exp_mod_normal_lcdf_test.cpp`:
for each row and each of lcdf, lccdf and cdf, the `double` call and the
`var` value, `dy`, `dmu`, `dsigma`, `dlambda`, with `|expected - actual| <=
1e-12 |expected|`.  It also runs the two vectorized-sum tests.  Output:
`out/emg/check_ref_rows_<commit>.txt`.

Result for `927099efb8` (rows numbered from 0; `.` pass, `F` outside the
tolerance, `N` non-finite):

| row | z, v | lcdf | lccdf | cdf |
|---|---|---|---|---|
| 1 | -1, 1e-6 | FFFFFF (7.0e-11) | FFFFF. (1.1e-9) | FF.... (1.1e-9) |
| 2 | -30, 1e-6 | FFFFFF (2.7e-9) | FFFFF. (4.7e-8) | FFFFF. (1.3e-6) |
| 3 | -100, 1 | .....F (2.1e-11) | pass | pass |
| 6 | 3, 1e-3 | ....F. (2.9e-12) | ....F. (1.5e-12) | ....F. (2.9e-12) |

- The relative error of the value (or of the failing partial) is in
  brackets.  The other 9 rows pass for all three functions.
- Vectorized sums: lcdf value, `dmu`, `dsigma` fail, and `dy` at 2 rows and
  `dlambda` at 3 rows; lccdf `dy` at 2 rows.
- 33 row checks plus 10 vectorized checks = 43 failed checks, equal to the
  "43 failed expectations" of the swap test in the PR text.  No non-finite
  value in the 13 rows.
- Deep-tail rows on `927099efb8`: lcdf rows 0 to 5 fail every check
  (value `-inf`, partials NaN); cdf rows 0 to 5 have NaN partials; lccdf
  passes.  60 failed checks, all non-finite.
- `686fc67261` gives the same lists.  `e4745c87ff` (develop) fails 179
  checks, 134 of them non-finite; it adds lcdf rows 4, 9, 12, lccdf rows 4,
  10, 11, 12 and cdf rows 3, 4, 12.
- `a17bba2e9e`, `d081fff3a4`: the 13 rows and the vectorized sums pass;
  only deep-tail rows fail (M4).  `1f8b6ba94d`, `1e46958b98`: everything
  passes, `SUMMARY failing row checks 0, non-finite results 0`.

## M2. Accuracy after a17bba2e9e

Claim: "Worst error over the value and the four partials after the change,
no non-finite result: 510 ulp (lcdf) and 480 ulp (lccdf) on the grid
above, and 640 and 900 ulp on a denser grid of 2500 points with `|z|` in
[1e-4, 100] and `v` in [1e-10, 100].  One exception: `dsigma` near
`z = 1/v` at large `v` (3.9e3 ulp at `(0.0316, 31.6)`)".

The dense grid has 50 `z` (`±10^k`, `k` from -4 to 2 in steps of 0.25)
and 25 `v` (`10^k`, `k` from -10 to 2 in steps of 0.5): 1250 `(z, v)`
points per function.  "2500 points" counts the lcdf and the lccdf
together (2500 output lines).

```sh
python3 emg/score_emg_grid.py out/emg/ref_cache_dense.json \
  out/emg/probe_emg_dense_{927099efb8,a17bba2e9e,1f8b6ba94d,1e46958b98}.txt \
  > out/emg/score_dense.txt
```

The first run computes the reference on all cores (a few minutes on 12
cores); later runs read the cache.

| grid | function | a17bba2e9e | 1f8b6ba94d | 1e46958b98 | claim |
|---|---|---|---|---|---|
| 170 | lcdf | 512 (100, 3) | 512 | 512 | 510 |
| 170 | lccdf | 479 (100, 100) | 479 | 479 | 480 |
| dense, no dsigma | lcdf | 640 (56.2341, 10) | 640 | 640 | 640 |
| dense, no dsigma | lccdf | 902 (56.2341, 100) | 902 | 902 | 900 |
| dense, dsigma | lcdf | 3896 (0.0316228, 31.6228) | 3896 | 3896 | 3.9e3 |
| dense, dsigma | lccdf | 3898 (0.0316228, 31.6228) | 3898 | 3898 | 3.9e3 |

- Non-finite results: 0 on both grids at the three commits.
- The second-largest `dsigma` error, 2.0e3 at `(0.01, 100)`, is also on
  `z = 1/v`.
- `1f8b6ba94d` changes only the left-tail regime (grid lcdf `dy` 66 to
  150 ulp; dense lcdf `dy` 49 to 94 ulp); the worst values do not change.
- For comparison, `927099efb8` on the dense grid: lcdf 1.2e15 ulp, lccdf
  2.3e14 ulp.
- Route check: `d/dy` from `mpmath.diff` and the analytic `d/dy` agree to
  2.9e-85 (grid) and 6.9e-89 (dense).

## M3. OpenCL device code against the CPU code

Claim: "Compiled on the host, the device code agrees with the CPU code to
5e-14 relative on the 2500-point grid (1.2e-12 at the `dsigma` zero)".

```sh
emg/ocl/run_cmp.sh 39886e0eb8 1e46958b98
```

`gen_dev_src.cpp` prints the device source (`erfcx` and
`exp_mod_normal_lcdf` device functions) of the tree; `cmp_host.cpp`
compiles it as C++ in a namespace and compares
`exp_mod_normal_lcdf_part` with `internal::exp_mod_normal_lcdf_kernel` for
the value, `d_z`, `d_v` and `d_s`.  Two grids: `dense` (the 1250 points of
M2) and `fine` (243 `z` x 121 `v` = 29403 points, the grid of the rescued
`cmp.cpp`).

The rescued files in `emg/ocl/rescued/` are incomplete: `dev_as_c.c` is a
generated file of one tree with a C prelude, and `cmp.cpp` has no build
script.  `gen_dev_src.cpp`, `cmp_host.cpp` and `run_cmp.sh` replace them.
The device source of `39886e0eb8` equals the rescued `dev.cl` (whitespace
ignored).

Worst relative difference over value, `d_z`, `d_v` (and `d_s`):

| grid | function | 39886e0eb8 | 1e46958b98 |
|---|---|---|---|
| dense | lcdf | 2.04e-14 (1.20e-12) | 2.45e-14 (1.20e-12) |
| dense | lccdf | 1.36e-14 (1.20e-12) | 2.44e-14 (1.20e-12) |
| fine | lcdf | 4.87e-14 (1.20e-12) | 4.87e-14 (1.20e-12) |
| fine | lccdf | 4.84e-14 (1.20e-12) | 4.84e-14 (1.20e-12) |

- The `d_s` maximum is at `(0.01, 100)`, on `z = 1/v`, in all four runs.
- No finite / non-finite mismatch.
- The 5e-14 of the claim matches the fine grid.  On the 2500-line dense
  grid of the claim the maximum is 2.45e-14.

## M4. Non-finite results on the 105-point deep grid

Claim (1f8b6ba94d): "On a grid of 105 points from z = -1e300 to 10 and v =
1e-300 to 1e10, non-finite results drop from 38 to 0 for the cdf, from 1 to
0 for the lccdf, and from 38 to 28 for the lcdf, where all 28 are at |z| >=
1e155 [...] (develop has 89, 89 and 56)."  The bot asks: "which deep-tail
rows failed with the kernel before 1f8b6ba94d".

Grid: 15 `z` x 7 `v`, `var` arguments.  A point counts if a partial is
non-finite; for the cdf also if the value is non-finite; for the lcdf and
lccdf if the value is NaN (`-inf` is permitted).  Source:
`emg/probe_deep_grid.cpp` (sections 4 and 5 of `emg/orig/probe_steve.cpp`).
Output: `out/emg/probe_deep_grid_<commit>.txt`.

| commit | cdf | lcdf | lccdf |
|---|---|---|---|
| e4745c87ff (develop) | 89 | 89 | 56 |
| 927099efb8 | 75 | 75 | 1 |
| a17bba2e9e | 38 | 38 | 1 |
| cf177ef96b | 38 | 38 | 1 |
| d081fff3a4 | 38 | 38 | 1 |
| 1f8b6ba94d | 0 | 28 | 0 |
| 1e46958b98 | 0 | 28 | 0 |

- All counts of the claim are reproduced.  The 28 lcdf points at
  `1f8b6ba94d` are the 7 `v` at each of `z = -1e155, -1e160, -1e200,
  -1e300`, with value `-inf` and `dsigma = inf`.
- The original outputs (`emg/orig/out_probe_cf177ef96b.txt`,
  `out_probe_e4745c87ff.txt`) give the same counts.

Deep-tail rows of the test (M1b program) on `d081fff3a4`, the parent of
`1f8b6ba94d`:

| row | z, v | lcdf | lccdf | cdf |
|---|---|---|---|---|
| 0 | -1e5, 1e-300 | ..NNNN (partials inf) | pass | ..NNNN (NaN) |
| 1 | -1e20, 1e-300 | ..NNNN (partials inf) | pass | ..NNNN (NaN) |
| 2 | -1e100, 1e-300 | NNNNNN (value -inf) | pass | ..NNNN (NaN) |
| 3 | -1e150, 1e-10 | ..NNNN (partials inf) | pass | ..NNNN (NaN) |
| 4 | -1e154, 0.1 | ..NNNN (partials inf) | pass | ..NNNN (NaN) |
| 5 | -1e154, 1 | pass | pass | pass |

42 failed checks, all non-finite.  `cf177ef96b` and `a17bba2e9e` give the
same list.  At `1f8b6ba94d` all rows pass.

## M5. Reference generators

Claim: the 13 rows are mpmath at 80 digits; the 6 deep-tail rows are
mpmath at 800 digits from two routes.

```sh
python3 refs_emg/ref_emg_lcdf_tests.py > out/emg/ref_rows_13.txt
python3 refs_emg/ref_emg_lcdf_deep_tail.py > out/emg/ref_rows_deep.txt
python3 refs_emg/compare_header.py refs_emg/exp_mod_normal_lcdf_reference.hpp \
  out/emg/ref_rows_13.txt out/emg/ref_rows_deep.txt
```

Run time: 3 s and a few minutes.  The agreement of the routes goes to
stderr (`out/emg/ref_rows_13.err`, `out/emg/ref_rows_deep.err`).

- `ref_emg_lcdf_tests.py` works at 80 digits (`mp.mp.dps = 80`), as the
  docstring says.  It had one route; route B is new: the values from
  `F = g (E(u) - E(t))` and `1 - F = g (E(-u) + E(t))` with
  `E(x) = exp(x^2) erfc(x)`, the partials from the closed form, at 120
  digits.  Smallest agreement: 73.3 digits (row `z = -30, v = 1e-6`).
- `ref_emg_lcdf_deep_tail.py` compared the routes only for `|z| <= 1e20`.
  Now it runs the direct route at all 6 rows and prints the agreement per
  row: from 7.49e-490 to 6.29e-205 (the largest at `z = -1e100`, where the
  direct route cancels about 400 digits).
- Comparison with the header of `1e46958b98`:
  `reference (13 rows): ... 0 of 182 numbers differ` and
  `deep tail (6 rows): ... 0 of 84 numbers differ`.  The rescued row files
  (`refs_emg/emg_ref_rows.txt`, `refs_emg/deep_rows.txt`) also equal the
  header.

## B1. Bitwise-identical results

Claims: f4cd10ee4a and c30ee6f9ef "change only when Eigen evaluates an
expression, so all results are bitwise identical"; the same for
927099efb8 (`skew_normal_lpdf`).

```sh
bitwise/run_compare.sh > out/bitwise/summary.txt
```

`bitwise/compare.cpp` prints in hex the value and all partials of
`lognormal_cdf`, `lognormal_lcdf`, `lognormal_lccdf`, `normal_lcdf`,
`normal_lccdf`, `exp_mod_normal_lcdf`, `exp_mod_normal_lccdf`,
`exp_mod_normal_cdf` and `skew_normal_lpdf` (`propto` false and true), on
11 points from `z = -40` to `40`, for these argument types: `double`
scalars, `std::vector<double>`, Eigen vectors of `double`, `var` scalars,
`double` y with `var` parameters, `std::vector<var>`, Eigen vectors of
`var`, `var_value<VectorXd>`, `fvar<double>` scalars (one tangent per
argument) and `std::vector<fvar<double>>`.  500 lines and 2870 numbers per
build, g++ 11.4.0 `-O3`, x86-64 without `-march`.

| parent vs commit | lines | differing lines |
|---|---|---|
| 0d4ac36e34 vs f4cd10ee4a | 500 | 0 |
| f4cd10ee4a vs c30ee6f9ef | 500 | 4 |
| 9b38ecbcb0 vs 927099efb8 | 500 | 0 |

The 4 lines of c30ee6f9ef hold 10 differing numbers, all in
`exp_mod_normal_lcdf` / `_lccdf` with vector arguments:

```
exp_mod_normal_lcdf std_vector_var field 1 of 35 -0x1.0c526dbaedb63p+11 -0x1.0c526dbaedb62p+11
exp_mod_normal_lcdf eigen_var field 11 of 25 0x1.266a31209398ap-3 0x1.266a31209398bp-3
exp_mod_normal_lccdf var_value_vector field 8 of 25 -0x1.1977b0b31b02fp-10 -0x1.1977b0b31b038p-10
exp_mod_normal_lccdf var_value_vector field 19 of 25 -0x1.c7c6eb921278cp-9 -0x1.c7c6eb9212799p-9
```

- The differences are 1 to 13 ulp.  So the claim holds for f4cd10ee4a and
  927099efb8, and not for the `exp_mod_normal` part of c30ee6f9ef.
- The `exp_mod_normal` code of c30ee6f9ef is the `log_diff_exp` form,
  which a17bba2e9e replaced.  The normal and lognormal outputs of
  c30ee6f9ef are bitwise identical.
- The cause was not examined.  A probable cause is that an evaluated
  array and a lazy expression go through different Eigen packet paths
  (`exp`, `log`, sums).
