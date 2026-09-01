# Parallel workstream handouts

Eight pieces of work, written to be run in parallel in separate git worktrees.
Each handout is self-contained: it states what is wrong, the evidence, the
constraints that must not be broken, and how to verify.

**Most of this has shipped.** Read the State column before picking anything
up. Three rows still hold live work, and one of them is a decision rather than
code: C (decide), E (upstream), F (backlog). H's remaining item is blocked
rather than available -- see below before picking it up.

Everything here was investigated against `main` at **v1.27.0** on
Windows 11, Python 3.12.7, with pandas 2.3.2 / anndata 0.12.2 /
spatialdata 0.7.3 / zarr 3.1.3. `main` is **v3.0.0** now, and three of those
four are below the floors `pyproject.toml` declares today -- anndata
`>=0.13.2`, spatialdata `>=0.8.0`, zarr `>=3.1.6`. Re-measure before trusting
a figure below: the anndata and zarr moves in particular changed how the store
gets written, which is the subject of half these handouts.

| # | Handout | Branch | Worktree | State |
|---|---|---|---|---|
| A | [pandas-3-string-dtypes.md](pandas-3-string-dtypes.md) | `fix/pandas3-string-dtypes` | `../Thyra-pandas3` | **SHIPPED** -- issue #117 closed, anndata floor raised, coercion deleted |
| B | [optimize-chunks.md](optimize-chunks.md) | `fix/optimize-chunks-sparse` | `../Thyra-chunks` | **RESOLVED** -- flag survives as a hidden deprecated no-op |
| C | [write-amplification.md](write-amplification.md) | `perf/write-amplification` | `../Thyra-perf` | **MEASURED, decision owed** -- **[findings](findings-write-amplification.md)** |
| D | [toolchain-hygiene.md](toolchain-hygiene.md) | `chore/toolchain-hygiene` | `../Thyra-toolchain` | **RESOLVED** -- `pre-commit run --all-files` passes clean |
| E | [upstream-lazy-table-pr.md](upstream-lazy-table-pr.md) | *upstream* | *scverse/spatialdata* | **OPEN upstream** -- scverse/spatialdata#1055 still unmerged (2026-08-03) |
| F | [interpolation-resampling.md](interpolation-resampling.md) | *not created* | `../Thyra-interp` | **backlog** -- a capability gap, not a defect |
| G | [resampling-scils-alignment.md](resampling-scils-alignment.md) | `fix/resampling-scils-alignment` | `../Thyra-resampling` | **DONE** -- v2.1.0-v2.2.3, item 5 left |
| H | [loose-ends-after-scils-alignment.md](loose-ends-after-scils-alignment.md) | one per item | `../Thyra-loose` | items 3 and 4 shipped; item 2 is **BLOCKED**; item 1 is not Thyra's |

Handout H is what was found *around* G and deliberately not folded into it.
Three of its four items are settled, so do not read it front-to-back:

- **Item 1** -- Ousia importing Thyra against a spatialdata/anndata/zarr
  cluster that violates Thyra's declared pins, silently, because the
  dependency is a bare editable path dep with no version constraint -- is a
  defect in *Ousia's* dependency declaration. Nothing in this repository
  fixes it, and no work here should wait on it.
- **Item 2**, `main` having no branch protection, is the only item left and is
  **BLOCKED** as of 2026-08-03: the release-actor bypass the item assumed
  would exist does not. The branch-protection API does still answer
  "Branch not protected", so the defect is real -- it is the fix that has no
  route. Read "The trap" in the item before starting. The concurrency half
  shipped separately in #144.
- **Item 3**, `release.yml`'s dead "No release needed" branch, shipped in
  `11c808f`.
- **Item 4**, taking spectrum representation explicitly, shipped in `5c7aa29`
  as `--spectrum-type`.

Handout G is the one to read before touching resampling. It records what the
SCiLS Lab manual actually specifies -- SCiLS is the baseline Thyra's
resampling was designed against -- and settles two questions that earlier
material left open: intensity is a count rather than a density, and the
default paths measure exact. It also corrects F's description of
`tic_preserving`.

Handout E is work in `scverse/spatialdata`, not here. It is listed because
it is the critical path for Ousia reading Thyra output lazily, and because
its findings change how B and C should be judged.

Handout F is backlog, not a defect: a capability `main` does not have, whose
only prior art was about to be deleted with a stale branch. It is written up
so the decision can be made deliberately rather than by attrition.

## The lazy-reading picture

Worth reading before B or C. spatialdata PR #1055 (yours) adds
`read_zarr(..., lazy=True)`. It was tested against stores from all three
Thyra write paths and **works today**:

```
in_memory      X=Array  obs=Dataset2D  block_equal=True
streaming_pcs  X=Array  obs=Dataset2D  block_equal=True
streaming_coo  X=Array  obs=Dataset2D  block_equal=True
```

Two consequences:

- **Thyra's on-disk format is already right for lazy reading.** Nothing
  needs to change here to support it. What must not happen is a change that
  breaks it, which is why B and C both carry a "assert `read_lazy` still
  works" constraint.
- **Wide `var` axes are the design target, not a defect.** The PR's own
  benchmark is 100,000 pixels x 100,000 m/z bins. That materially weakens
  handout C's lead 2 (default bin counts upsampling 45x): lazy loading is
  the answer to wide axes, so shrinking the default is less obviously
  right than it looked. Handout C says so.

## Ordering and overlap

Historical, and no longer binding. **A has shipped**, so the constraint it
imposed on everything else is gone, and B and D are resolved too.

What is left of this section: C is free to start from current `main`. But
what C needs now is a **decision** on
[the findings](findings-write-amplification.md) -- specifically whether the
default bin count stays as it is on the `tic_preserving` path -- not more
measurement. Read the findings before writing any code for it.

## Worktree setup

Each handout repeats this, but for convenience:

```bash
git worktree add -b <branch> ../<Thyra-dir> main
```

Two things to know about working in a worktree here:

1. **The venv is shared.** `thyra` is installed into the Poetry venv as an
   editable install pointing at the *main* checkout, so
   `poetry run python` from a worktree still imports the main checkout's
   code. To exercise worktree code, put the worktree first on the path:

   ```bash
   PYTHONPATH=$(pwd) poetry run pytest
   ```

   This bit me while verifying handout A: two runs reported identical
   results because both imported `main`.

2. **Line endings.** Fixed by handout D. Commits no longer fail the
   `mixed-line-ending` hook on the first attempt, and
   `pre-commit run --all-files` is safe to run as a pre-push check: it passes
   and leaves `git status` empty. Worktrees created before that landed still
   hold CRLF files; `git add --renormalize .` settles them without producing
   a diff.

   One thing worktrees still hide: bandit silently skips everything under
   `.github/`, because its default exclusions contain the bare string `.git`
   and a worktree's `.git` is a file rather than a directory. Verify bandit
   findings from an ordinary checkout.

## Shared conventions

- Conventional commits (`feat`, `fix`, `docs`, `refactor`, `chore`, `perf`,
  `build`, `style`, `test`).
- Commit messages describe the change, never the tooling used to make it. No
  authorship or co-authorship trailers for anything but a person.
- No emojis anywhere: code, comments, docstrings, commits, docs.
- Before every commit:
  ```bash
  poetry run black . && poetry run isort . && poetry run flake8 && poetry run pytest
  ```
- `poetry run mkdocs build --strict` must stay clean.
- Two pre-existing `no-any-return` mypy errors exist in
  `thyra/tools/make_example_data.py:62` and
  `thyra/metadata/extractors/waters_extractor.py:130`. They are not yours
  unless you are on handout D.

## Branch cleanup — done 2026-07-30

This section used to list two dead branches. The sweep it asked for has
happened and went further: 56 stale local branches and 5 remote ones are
gone, leaving `main`, `gh-pages`, and the three live lane branches.

Each was classified before deletion, not deleted by age:

| Basis | Count |
|---|---|
| ancestor of `main` (`git branch -d`, git verifies the merge itself) | 22 |
| squash/rebase-merged — every commit's patch already upstream (`git cherry` all `-`) | 14 |
| superseded, verified individually | 20 |

SHAs are recorded in `.git/deleted-branches-2026-07-30.txt`; restore any with
`git branch <name> <sha>`.

Two findings came out of the sweep and are the reason it was worth doing
branch by branch rather than in bulk:

- `feature/lazy-loading-support` carried the **only** test of the
  `encoding-type` / `encoding-version` attrs the PCS path hand-writes —
  `git grep encoding-type tests/` on `main` was empty. Its implementation
  was genuinely superseded (`main` distinguishes `"string"` for 0-d scalars
  in `uns` from `"string-array"` for arrays, which the branch's binary
  `_set_encoding_attrs(is_string=...)` helper could not express), but the
  test was not. Ported to `main` in #108 *before* the branch was deleted.
- `feature/conservative-interpolation-implementation` held a capability
  `main` still lacks. Archived as tag `archive/conservative-interpolation`
  and written up as handout F.

Everything else was already on `main` in equal or better form.
