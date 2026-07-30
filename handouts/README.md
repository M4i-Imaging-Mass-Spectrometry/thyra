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
| C | [write-amplification.md](write-amplification.md) | `perf/write-amplification` | `../Thyra-perf` | after A |
| D | [toolchain-hygiene.md](toolchain-hygiene.md) | `chore/toolchain-hygiene` | `../Thyra-toolchain` | independent |

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

2. **Line endings.** Every commit that touches a file will fail the
   `mixed-line-ending` hook once, rewrite the file to LF, and need a second
   `git add` + `git commit`. This is expected until handout D lands. Do not
   run `pre-commit run --all-files` as a check: it rewrites ~137 unrelated
   files.

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

## Branch cleanup (do this once, anywhere)

Two remote branches are dead and should be deleted so they stop looking
like work in progress:

- `origin/fix/streaming-pcs-image-write` — already merged. `git cherry -v
  main origin/fix/streaming-pcs-image-write` reports `-`, and
  `tests/unit/converters/test_streaming_pcs_roundtrip.py` is on `main`.
- `origin/feature/lazy-loading-support` — superseded, 83 commits behind.
  `main` already writes every `encoding-type` / `encoding-version` attr the
  branch adds, and does it **better**: `main` distinguishes `"string"` for
  0-d scalars in `uns` from `"string-array"` for arrays, which the branch's
  binary `_set_encoding_attrs(is_string=...)` helper cannot express.
  Verified that `anndata.experimental.read_lazy()` already works on `main`
  output (see handout C for the transcript).

`origin/fix/streaming-coo-arrow-string-write` is **live and wanted** — it is
the subject of handout A. Do not delete it.

The stale worktree at `../Thyra-coo-fix` sits on that branch 19 commits
behind `main`. Handout A supersedes it; remove it with
`git worktree remove ../Thyra-coo-fix` once A has landed.
