# Reproduction scripts for stan-dev/math#3438

This directory holds the scripts behind the numbers in the text and the
comments of stan-dev/math#3438 ("Return the final Newton iterate as the
Laplace latent mean, and fix the solver-3 log determinant").  Each
section below names one claim, the commands that measure it, the
expected output, and the measured commits.

Last run: 2026-10-28, on a laptop (Ubuntu 22.04, g++ 11.4.0, `-O1`,
Python 3.10.12, mpmath 1.3.0).  The R script ran on 2026-10-09 with
CmdStan 2.40.0 and CmdStan `develop` (see §"Examples A, B and C").

## Layout

| path | content |
|---|---|
| `build.sh` | compiles one C++ program against one Stan Math tree |
| `pr_laplace_mean_lu_logdet.R` | examples A, B and C of the PR text (cmdstanr) |
| `probe_saddle.cpp` | example C final iterate; example D (odd permutation) |
| `issue3469/example.cpp` | the example program of stan-dev/math#3469, unchanged |
| `issue3469/scan.cpp` | the 105-point grid of the same example |
| `refs/ref_issue3469.py` | mpmath: value at θ = 0 and Laplace value, #3469 example |
| `refs/ref_example_d.py` | mpmath: mode and Laplace value of example D |
| `results/` | the outputs of the run of 2026-10-28 |

All commands run from this directory (`pr3438/`).  New outputs go to
`out/`, which git does not track; compare them with `results/`.

## Requirements

- Python 3 with mpmath.
- g++ with C++17.
- A clone of stan-dev/math with the PR branch, and the `stan/` directory
  of the merge base `e4745c87ff`, which the texts call "develop":

  ```sh
  git clone git@github.com:stan-dev/math.git stan_math
  cd stan_math
  git fetch origin pull/3438/head:laplace-mean-lu-logdet
  git checkout laplace-mean-lu-logdet
  make -f make/standalone math-libs   # builds lib/tbb
  mkdir -p ../dev_e4745c87ff
  git archive e4745c87ff stan | tar x -C ../dev_e4745c87ff
  cd -
  ```

Set these variables once:

```sh
export PR=$PWD/../stan_math              # the PR tree (stan/ and lib/)
export DEV=$PWD/../dev_e4745c87ff        # stan/ of the merge base
export STAN_MATH_LIB=$PR/lib             # lib/ for both trees
export TBB_LIB=$PR/lib/tbb               # a directory with libtbb.so
mkdir -p out
build() {  # build <source> <name>: one binary per tree
  STAN_MATH=$DEV OPT=-O1 ./build.sh "$1" "out/$2_dev"
  STAN_MATH=$PR OPT=-O1 ./build.sh "$1" "out/$2_pr"
}
```

The PR changes no file in `lib/`, so both trees can use the same `lib/`.

## Examples A, B and C (the table of the PR text)

Claim: the table of the PR text, and "The posterior sd is 5.6 to 9.4
times smaller than the prior sd."

```sh
Rscript pr_laplace_mean_lu_logdet.R                         # CmdStan default
CMDSTAN_PATH=<CmdStan develop> STAN_MATH=$PR/ Rscript pr_laplace_mean_lu_logdet.R
```

The script needs cmdstanr, posterior, numDeriv and Matrix.  It is the
script in the PR text, with generic paths in the usage lines.  The
`develop` column is from CmdStan 2.40.0; the "this PR" column is from
CmdStan `develop` with `STAN_MATH` set to the PR tree.

## Example C: the final iterate from (0.01, 3)

Claim: "From the start `(0.01, 3)`, the Newton iteration of solver 3 ends
at `(1.8e-6, 5.9e-9)`, near the saddle point.  At that iterate,
`K^-1 + W` has the eigenvalues -0.099 and 101.3."

```sh
build probe_saddle.cpp probe_saddle
out/probe_saddle_dev > out/probe_saddle_dev.txt
out/probe_saddle_pr > out/probe_saddle_pr.txt
diff out/probe_saddle_dev.txt results/probe_saddle_e4745c87ff.txt
diff out/probe_saddle_pr.txt results/probe_saddle_3c5e80ef4c.txt
```

On the merge base, `internal::laplace_marginal_density_est` returns the
final iterate, because the old log determinant does not throw.  Expected
(merge base, block "C, start (0.01, 3)"):

```text
theta (1.80065e-06, 5.92318e-09), |grad log p(theta | y)| 1.78e-07
W (-0.43133, 100), eigenvalues of K^-1 + W (-0.0990923, 101.334)
det(B) -30.1244, permutation sign +1, LU diagonal (-0.725321, 41.5325)
laplace_marginal_tol: -nan
```

`K^-1 + W` has one negative eigenvalue, so the final iterate is the
saddle point, 1.8e-6 from `(0, 0)`.  The PR binary (3c5e80ef4c) throws
the error `det(I + Sigma * W) is not positive at the final iterate, ...`
for this start.  Later PR commits shorten this message to
`det(I + Sigma * W) <= 0 (saddle point or singular)`.

## Example D: an odd permutation with positive pivots

Claim: "It returned a finite value when the permutation was odd and all
pivots were positive (example D)."

Same binaries as above, blocks "D".  Example D: `K = [1 -1.5; -1.5 4]`,
two Cauchy observations at `±0.3 sqrt(3)` with `s = 0.3`, one normal
observation `z = 0` with `sz = 0.1`.  At `θ[1] = 0` the Cauchy curvature
is at its minimum, `W[1] = -1 / (2 s^2)`, and `(0, 0)` (the default
initial guess) is a saddle point.

| start | final iterate | det(B) | e4745c87ff | 3c5e80ef4c |
|---|---|---|---|---|
| (0, 0) | (0, 0) | −576.8 | −4.449184855 | throws |
| (0.01, 0.01) | (−0.3504, 0.0030) | 1608 | NaN | −4.847223761 |
| (−0.3, 0) | (−0.3504, 0.0030) | 1608 | NaN | −4.847223761 |

At `(0, 0)` the LU factorization swaps the rows (permutation sign −1)
and both pivots are positive (8.33, 69.2), so the old
`sum(log(diag(U)))` is finite.  At the maximum, one pivot is negative.
The log posterior is symmetric in θ, so it has two maxima of equal
height, `±(−0.3504, 0.0030)`.  From `(0.01, 0.01)` the iteration passes
the saddle point and ends at the maximum with `θ[1] < 0`.  The unit test
therefore starts at `(−0.3, 0)`, inside one basin.  The last four blocks
of the probe output use the settings of the unit test (tolerance 1e-12,
1000 steps): from `(−0.3, 0)` the errors against the mpmath reference are
1.8e-15 in the value and 6.3e-16 in the mean.  The PR head after
3c5e80ef4c has the shorter message `det(I + Sigma * W) <= 0 (saddle
point or singular)`.

The reference:

```sh
python3 refs/ref_example_d.py > out/ref_example_d.txt
diff out/ref_example_d.txt results/ref_example_d.txt
```

It computes the mode by Newton and by `mp.findroot` (50 digits), and the
Laplace value −4.8472237608447798 at the mode `(−0.35038406125098763,
0.0029862277947527355)`.  The unit test
`laplace_marginal_tol.solver_3_odd_permutation_saddle_point_throws`
uses these values.

## The #3469 example: 8 points

Claim (PR comment of 2026-10-18): "This PR prints the same values,
except at μ = 1, α = 1e9, where develop gives +inf and this PR throws."

```sh
build issue3469/example.cpp example
out/example_dev; out/example_pr
```

## The #3469 example: 105-point grid

Claim (PR comment of 2026-10-18): "develop gives +inf at 87 points, this
PR throws at the same 87 points, and the other 18 values are identical."

```sh
build issue3469/scan.cpp scan
out/scan_dev > out/scan_dev.txt
out/scan_pr > out/scan_pr.txt
paste out/scan_dev.txt out/scan_pr.txt | awk '
  {n++; if ($3 == "inf") inf++; if ($6 == "throws") thr++
   if ($3 == "inf" && $6 == "throws") both++
   if ($3 != "inf" && $6 != "throws") {other++; if ($3 == $6) same++}}
  END {print n, "points;", inf, "inf;", thr, "throws;", both, "both;",
       other, "other;", same, "identical"}'
```

Expected: `105 points; 87 inf; 87 throws; 87 both; 18 other; 18
identical`.  The grid is μ ∈ {−5, 0, 1, 2, 5} and log10 α = 7 to 12 in
steps of 0.25.  The outputs are byte-identical to `results/scan_*.txt`
(first run 2026-10-18, rerun 2026-10-28).

## The #3469 example: values at θ = 0

Claim (PR comment of 2026-10-18): "the μ = −30, α ≥ 1e5 rows still
return the value at θ = 0."

```sh
python3 refs/ref_issue3469.py > out/ref_issue3469.txt
diff out/ref_issue3469.txt results/ref_issue3469.txt
```

80 digits.  The value at θ = 0, `ll(0) − log det B(0) / 2`, comes from
`mp.det` and from the eigenvalues of `e^μ K` (closed form).  The Laplace
value comes from the mode by damped Newton and by `mp.findroot`.  The
last two columns print the differences between the routes.  Expected
rows for μ = −30:

```text
 -30 1.0e+5 | -429.759638992 | -22.1750873716 | ...
 -30 1.0e+7 | -431.690693525 | -26.7802575096 | ...
 -30 1.0e+9 | -436.285290904 | -31.3854276956 | ...
```

The merge base and the PR print −429.759639, −431.6906935 and
−436.2852909 for these rows: the value at θ = 0 to all printed digits.

## Not in this directory

The GP-model point of #3469 (ρ = 1, μ = 1, α = 1e8).  The issue author
found it with his model; the model and the data are not posted.
