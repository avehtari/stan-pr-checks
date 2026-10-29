# Checks for Stan pull requests

This repository holds the scripts behind the numbers in my pull
requests to stan-dev/math and other Stan repositories.  The scripts
are reference generators, accuracy sweeps, probes, speed benchmarks,
and swap tests (the new tests run against the headers of the base
branch, to show that they fail without the change).
A pull request links to its directory at a fixed commit, for example
`https://github.com/<owner>/<repo>/tree/<commit>/pr3408`.

## Layout

One directory per pull request, named `pr<number>`.  Each directory is
self-contained: it has its own `README.md`, `build.sh` and data files,
and its commands run from that directory.  A directory does not import
from another one; copy a shared script instead.

| directory | pull request | topic |
|---|---|---|
| `pr3408/` | stan-dev/math#3408 | incomplete gamma gradient roots, 18 cdfs in log space |
| `pr3409/` | stan-dev/math#3409 | incomplete beta gradient roots (part 1 of 3) |
| `pr3434/` | stan-dev/math#3434 | lgamma, lbeta, digamma differences at large shapes |
| `pr3438/` | stan-dev/math#3438 | Laplace latent mean, solver-3 log determinant |
| `pr3442/` | stan-dev/math#3442 | Wiener functions: underflow, gradients, `hcubature` |

## Conventions for a PR directory

- `README.md`: the requirements, then one section per claim of the PR
  text.  Each section gives the quoted claim, the commands, the expected
  output and the commits that produced it.  For a speed claim it also
  gives the machine, compiler, flags and job.
- `build.sh`: compiles one C++ source against a Stan Math tree given by
  `STAN_MATH`.  No absolute paths in any file.
- `refs/`: mpmath reference generators.  Each reference comes from two
  independent mpmath computations, and the script prints their
  agreement.
- A file header of two or three lines names the purpose and the claim.
- Outputs go to `out/`, which is not tracked.
- After a merge or a new commit of the pull request, run the accuracy
  scripts again and update the expected outputs.

## Requirements

Python 3 with mpmath, g++ with C++17, and a clone of the repository of
the pull request.  Each `README.md` lists the versions it was tested
with.
