# PR #3406: checks of five claims in SteveBronder's commit messages

Last updated: 2026-10-30

The fact check of 2026-10-28 (stevieb-bot) marked these claims "Not
Verified" or did not check them.  The commits are in the PR branch.
Scripts: `steve/`.  Fetched outputs: `out/steve/`.

Machine for S2 and S3: Triton `batch-csl`, Intel Xeon Gold 6248 (2.5 GHz),
g++ 14.3.0 (conda-forge, module `scicomp-r-env`), Python 3.10, stanc3
nightly `10361e1`, `lib/` of develop `d12f488749` (it has the ARKODE
headers).
Machine for S4: Triton `pht114` (`batch-skl`), Intel Xeon Gold 6230
(2.1 GHz), whole node (`--exclusive`), governor `performance`, g++ 14.3.0.

| claim | commit | result |
|---|---|---|
| S1 | 9a8deffbef | not reproduced: no test at the commit gives 256 |
| S2 | 201dcab952 | reproduced: 33.1 -> 9.5 GB (EMG), -37 % (normal, lognormal) |
| S3 | f37ea7cb67 | reproduced: 6 Prim tests fail at the parent, 0 at the commit |
| S4 | 2015557859 | partly: 21 % faster (double), still 29 % slower than develop |
| S5 | cf177ef96b | the `x = abs(z)` ternary was removed; `if (z < 0)` at L114 stays |

---

## S1. 256 proportional-density checks (commit 9a8deffbef)

> Validation: targeted value and nested-autodiff tests pass, along with 256
> proportional-density checks.

Commands:

```bash
git show 9a8deffbef -- test/
git ls-tree -r --name-only 9a8deffbef | grep -i propto     # nothing
# generated test/prob instances at 9a8deffbef (generator compiled locally)
cd steve/s1 && git -C $CLONE archive 9a8deffbef test/prob | tar x
g++ -O0 -I$CLONE/lib/boost_1.87.0 test/prob/generate_tests.cpp -o gen
./gen test/prob/exp_mod_normal/exp_mod_normal_test.hpp 100
./gen test/prob/skew_normal/skew_normal_test.hpp 100
```

What the commit changes in `test/`:

- `test/unit/math/prim/prob/std_normal_lcdf_impl_test.cpp`: `scalar_tail_kernel`
  and `vectorized_tail_kernel` (adapted to the new value-only return type;
  the vectorized test now also runs on the input replicated 8 times) and a
  new `value_only_storage` test.
- `test/unit/math/mix/prob/std_normal_lcdf_impl_test.cpp`: one line (return type).

None of these tests calls a density with `propto`.  The tree at 9a8deffbef
has no `expect_propto*` helper.  `expect_propto_vectorized.hpp` and the
propto tests for `exp_mod_normal_lpdf` and `skew_normal_lpdf` came later,
in 8cff6fd3f3 (2026-10-13).

The only propto checks at 9a8deffbef for the densities that the commit
changes (`exp_mod_normal_lpdf`, `skew_normal_lpdf`) are the `Propto` tests
of the generated `test/prob` framework.  `test_propto()` returns at once
for all-double arguments and for any vector argument.  Counts from the
generator output (default build, no `STAN_PROB_TEST_ALL`):

| distribution | `Propto` instances | non-trivial | valid points | `EXPECT_NEAR` |
|---|---|---|---|---|
| `exp_mod_normal_lpdf` | 1296 + 1215 + 1215 | 15 v + 15 ffv + 15 vv = 45 | 3 | 135 |
| `skew_normal_lpdf` | 1296 + 1215 + 1215 | 45 | 4 | 180 |

- Sum of non-trivial instances: 90.  Sum of `EXPECT_NEAR` calls: 315.
  Var-only (`v`) instances: 30, `EXPECT_NEAR` calls: 105.
- No count from the tree gives 256.  256 = 4^4 is the number of type
  combinations of a 4-argument density with 4 types per argument (for
  example double, var, `std::vector<double>`, `std::vector<var>`).  That is
  a possible origin, but no test at 9a8deffbef enumerates it.

Result: **not reproduced**.  The 256 checks are not in the commit or in the
tree at the commit.  No cluster job.

---

## S2. Compile memory of the generated expression tests (commit 201dcab952)

> Peak g++ -O0 memory compiling the generated expression tests drops from 33
> GB to 9 GB for each exp_mod_normal lcdf/lccdf/cdf and by about 35% for the
> normal and lognormal ones

Job: `steve/s2s3/run_expr.sh`, array job **21302464** (tasks 0-9),
64 GB per task, 2 cores.  Each task unpacks one commit into `$TMPDIR`,
links the shared `lib`, writes `make/local` = `CXX=g++`, `O=0`, downloads
the nightly stanc3 with `make test/expressions/stanc`, and generates the
test for one function with `generate_expression_tests.main([fn], 1)` (the
call that `runTests.py test/expressions --only-functions fn` makes).  It
then runs `/usr/bin/time -v make test/expressions/tests0_test.o`.  Compile
command (from `make -n`):

```text
g++ -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes
  -Wno-class-memaccess -isystem lib/tbb_2020.3/include -O0 -I . -isystem
  lib/eigen_5.0.1 -isystem lib/boost_1.87.0 -isystem lib/sundials_6.1.1/include
  ... -DBOOST_DISABLE_ASSERTS -DGTEST_HAS_PTHREAD=0 -c ... tests0_test.cpp
```

CI (`Jenkinsfile`, stage `expr`) uses `CXX=clang++ -Werror`, `O=0` and a
64 GiB pod.  This check uses g++ -O0, as the claim says.

Peak resident set size ("Maximum resident set size"), one compile each:

| function | TESTs | 26256d4cdd (parent) | 201dcab952 | change |
|---|---|---|---|---|
| `exp_mod_normal_lcdf` | 720 | 33.14 GB (28:38) | 9.53 GB (3:48) | -71 % |
| `exp_mod_normal_lccdf` | 720 | 33.15 GB (29:42) | 9.53 GB (3:46) | -71 % |
| `exp_mod_normal_cdf` | 720 | 32.97 GB (37:55) | 9.54 GB (3:10) | -71 % |
| `normal_lcdf` | 168 | 3.53 GB (1:06) | 2.24 GB (0:30) | -37 % |
| `lognormal_lcdf` | 168 | 3.68 GB (0:56) | 2.30 GB (0:38) | -37 % |

- GB = 2^30 bytes (kB / 1048576).  Times in parentheses are wall times.
- All compiles exited 0.
- The generated files at the two commits hold the same tests in a
  different order (Python set order); after the removal of digits and a
  sort, the files are identical.
- The `exp_mod_normal_cdf` cell at 26256d4cdd was fetched later
  (`out/steve/expr/pr3406_expr_21302464_2.out`).

Result: **reproduced** for all 5 cells.

---

## S3. Prim expression tests of std_normal lcdf/lccdf/cdf (commit f37ea7cb67)

> This failed the Prim expression tests for std_normal_lcdf, std_normal_lccdf
> and std_normal_cdf.

Same job, tasks 10 (201dcab952 = f37ea7cb67^) and 11 (f37ea7cb67): one
generated file with `std_normal_lcdf std_normal_lccdf std_normal_cdf`
(18 TESTs: 6 Prim, 6 Rev, 6 Fwd), built with `make tests0_test` (g++ -O0)
and run in full.

| commit | result |
|---|---|
| 201dcab952 (parent) | 12 of 18 pass; 6 fail, all `ExpressionTestPrim` |
| f37ea7cb67 | 18 of 18 pass |

Failing tests at the parent: `ExpressionTestPrim.std_normal_lccdf3`,
`std_normal_lcdf4`, `std_normal_cdf6`, `std_normal_lccdf7`,
`std_normal_lcdf8`, `std_normal_cdf11`.  The message is the same in each:

```text
test/expressions/tests0_test.cpp:17: Failure
Expected: (matrix0_expr1_counter) <= (int5), actual: 2 vs 1
```

That is the "evaluated at most once" check: the double vector expression
argument was evaluated twice.  The Rev and Fwd tests pass at both commits.

Result: **reproduced**.  Logs: `out/steve/expr/pr3406_expr_21302464_{10,11}.out`,
`run_S3_*.log`.

---

## S4. Scalar exp_mod_normal_lpdf at z ~ 9.4 (commit 2015557859)

> This fixes the scalar exp_mod_normal_lpdf benchmark (z ~ 9.4), where
> develop's erfc returns early.

The parent of 2015557859 is **1f8b6ba94d** ("Keep exp_mod_normal cdf
gradients finite in the deep lower tail", by avehtari), not e1fea0651c.

In `exp_mod_normal_lpdf` at 2015557859, z = (y - mu - lambda sigma^2) /
sigma is the argument of `std_normal_lcdf_value_grad`.  Develop computes
`log(erfc(-z / sqrt(2)))`; at z = 9.4 the erfc argument is -6.6, where
glibc erfc returns 2 early.

Benchmark: `steve/s4/bench_emg_scalar.cpp`; job `steve/s4/run_bench_s4.sh`,
job **21302510** (node `pht114`, whole node, governor `performance`; an
earlier copy on `batch-csl`, job 21302477, was cancelled before it started).

- mu = 0, sigma = 1, lambda = 1; a block of 1000 inputs with z uniform on
  [9.3, 9.5] (controls: [1.5, 2.5] and [-0.1, 0.1]).
- double: 4 double arguments, value; 4 x 10^6 calls per number.
- var: 4 var arguments, value + `grad()` + `recover_memory()`; 10^6 calls.
- Throughput loop with the accumulator bump; one binary per tree; `-O3`;
  4 interleaved passes (A B C A B C ...), `taskset -c 2`.
- The value sums agree across the three trees to the last digit (z2 `sum_dy`
  differs in the 16th digit on develop).

ns per call, mean of 4 passes (pass-to-pass spread < 1 %, except develop:
z0 double 3.7 %, z2 var 1.9 %, z9.4 var 1.6 %):

| z | type | 1f8b6ba94d | 2015557859 | e4745c87ff | 2015/1f8b | 2015/dev |
|---|---|---|---|---|---|---|
| 9.4 | double | 33.32 | 26.47 | 20.59 | 0.79 | 1.29 |
| 9.4 | var | 118.6 | 106.9 | 105.2 | 0.90 | 1.02 |
| 2 | double | 42.67 | 42.68 | 56.05 | 1.00 | 0.76 |
| 2 | var | 131.7 | 131.5 | 157.5 | 1.00 | 0.83 |
| 0 | double | 36.84 | 35.87 | 34.21 | 0.97 | 1.05 |
| 0 | var | 115.9 | 117.6 | 128.8 | 1.01 | 0.91 |

- Cycle check: 20.6 ns at 2.1 GHz is 43 cycles per call; no call was
  optimized away.
- The commit makes the double call at z ~ 9.4 21 % faster and the var call
  10 % faster.  The double call is still 29 % slower than develop; the var
  call is 2 % slower.

Result: **partly reproduced**.  The commit improves the case, but it does
not remove the gap to develop for doubles on this machine.

---

## S5. Branch on the sign of z (commit cf177ef96b)

> The scalar function no longer branches on the sign of z for doubles.

Command: `git show cf177ef96b -- stan/math/prim/fun/std_normal_lcdf_impl.hpp`;
`git show 150f54df16:stan/math/prim/fun/std_normal_lcdf_impl.hpp`.

For doubles, cf177ef96b replaces the ternary `x = (z < 0 ? -z : z) *
INV_SQRT_TWO` with `x = std::fabs(z) * INV_SQRT_TWO` (L92-97 at 150f54df16),
and it compiles out the `z <= -4 * SQRT_TWO` lower-tail branch (L75-88, now
`if constexpr (!std::is_arithmetic_v<T>)`).  The branch `if (z < 0)` at L114
remains for every type: it selects the lower-tail value `log(erfcx(x))`
(L115) or the upper-tail value `log1p(-tail)` (L123-126), so the commit
message is correct only for the computation of `x`.

No cluster job.

## Tree tarballs

The job scripts expect tree tarballs next to them; they are not in the
repository.  Recreate them from a clone of the PR branch:

```sh
# S2, S3: steve/s2s3/src_<commit>.tar.gz for 26256d4cdd, 201dcab952, f37ea7cb67
git -C <clone> archive --format=tar.gz -o $PWD/steve/s2s3/src_<commit>.tar.gz \
  <commit> make makefile runTests.py stan test
# S4: steve/s4/stan_<commit>.tar.gz for e4745c87ff, 1f8b6ba94d, 2015557859
git -C <clone> archive --format=tar.gz --prefix=t_<commit>/ \
  -o $PWD/steve/s4/stan_<commit>.tar.gz <commit> stan
```
