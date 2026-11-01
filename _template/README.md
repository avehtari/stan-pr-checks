# Template for a PR directory

Copy this directory to `pr<N>/` (or a working name until the PR has a
number), then replace every `<...>` and delete this section.

```sh
cp -r _template pr<N>
cd pr<N>
CLONE=<clone> ./make_trees.sh develop:<merge base> pr:<PR head>
STAN_MATH_LIB=<clone>/lib TBB_LIB=<clone>/lib/tbb ./run_probes.sh
cp out/*.txt results/        # only after a deliberate rerun
<wip>/tools/triton_fetch_logs.sh --tag <tag> --dest results/triton_<job> --scrub
```

Checks before the PR text is posted (`STEVIEB_BOT_REPLIES.md` §2 step 7):

1. `claim-check` on the PR text and on this README, with the job tags,
   `results/` and the clone.
2. Every `>` quote below matches the PR text word for word ("..." marks
   a cut).
3. `<wip>/tools/pr_overlap.sh --pr <N> --paths stan/` for the merge-order
   note.

---

# Checks for stan-dev/math#<N> (<topic>)

This directory holds the scripts behind the numbers in the PR "<title>"
(stan-dev/math#<N>).  Sections 1 to <k> quote the claims of the PR text
and give the commands and the expected output.

Last run: <date>.  The probes ran on <machine> (`results/versions.txt`).
The unit tests and the swap test ran on <cluster> (<compiler>).  In the
cluster logs, `<workstation>`, `<scratch>` and `<home>` replace the
machine paths.

## Commits

| name in this file | commit | what |
|---|---|---|
| develop | `<sha>` | merge base of the PR (develop of <date>) |
| PR | `<sha>` | PR head, tree `<tree sha>` |

## Layout

| path | content |
|---|---|
| `build.sh` | compiles one probe against one Stan Math tree |
| `make_trees.sh` | exports `stan/` and `test/unit/` of the commits into `out/trees/` |
| `run_probes.sh` | runs the probes and scorers; writes `out/` |
| `probes/` | the C++ probes |
| `refs/` | mpmath references, two routes each |
| `results/` | the outputs of <date> and the cluster job logs |

Run all commands from this directory.  New outputs go to `out/`, which
git does not track; compare them with `results/`.

## Requirements

- Python 3 with mpmath (tested with <version>).
- g++ with C++17 (tested with <version>).
- A clone of stan-dev/math with the PR branch, and a built TBB
  (`make -f make/standalone math-libs`).

## 1. <claim>

> <quote of the PR text, word for word>

<What the probe computes, at which points, against which references.>

```sh
<commands>
```

Expected (`results/<file>`): <the numbers of the quote, with the point of
each worst value and the relative error next to an absolute one>.

## <k>. Swap test

> <quote>

<Job id, node, compiler; the stages; the failed expectations per test
case, counted from the `.run` log of the develop stage.>
