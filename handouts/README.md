# Parallel workstream handouts

Four independent pieces of work, written to be run in parallel in separate
git worktrees. Each handout is self-contained: it states what is wrong, the
evidence, the constraints that must not be broken, and how to verify.

Everything here was investigated against `main` at **v1.27.0** on
Windows 11, Python 3.12.7, with pandas 2.3.2 / anndata 0.12.2 /
spatialdata 0.7.3 / zarr 3.1.3.

| # | Handout | Branch | Worktree | Priority |
|---|---|---|---|---|
| A | [pandas-3-string-dtypes.md](pandas-3-string-dtypes.md) | `fix/pandas3-string-dtypes` | `../Thyra-pandas3` | **first** |
| B | [optimize-chunks.md](optimize-chunks.md) | `fix/optimize-chunks-sparse` | `../Thyra-chunks` | independent |
| C | [write-amplification.md](write-amplification.md) | `perf/write-amplification` | `../Thyra-perf` | after A -- **[findings](findings-write-amplification.md)** |
| D | [toolchain-hygiene.md](toolchain-hygiene.md) | `chore/toolchain-hygiene` | `../Thyra-toolchain` | independent |
| E | [upstream-lazy-table-pr.md](upstream-lazy-table-pr.md) | *upstream* | *scverse/spatialdata* | independent |
| F | [interpolation-resampling.md](interpolation-resampling.md) | *not created* | `../Thyra-interp` | backlog |

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

**A must land first.** It is the only one fixing a defect users can hit
today, and it touches `_save_output` in
`thyra/converters/spatialdata/base_spatialdata_converter.py`.

**C also touches that file** (`_process_single_spectrum`, the write
orchestration). C should rebase on `main` once A has merged rather than
racing it. If C starts before A lands, expect one small conflict in
`_save_output`.

**B and D share no files with anything else** and can run start-to-finish
independently.

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
- Never mention Claude or AI authorship in commit messages.
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
