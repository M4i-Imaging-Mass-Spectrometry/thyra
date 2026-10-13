# Contributing to Thyra

Thank you for helping. This page holds the rules you need before your first pull
request. The full guide, with code style, test conventions and the release process,
is [docs/contributing.md](docs/contributing.md).

To report a security problem, do not open an issue. See [SECURITY.md](SECURITY.md).

## Set up

You need Python 3.12, 3.13 or 3.14, [uv](https://docs.astral.sh/uv/) and Git.

```bash
git clone https://github.com/<your-username>/thyra.git
cd thyra
uv sync --locked --group test --group dev
uv run pre-commit install
```

The dependencies live in [PEP 735](https://peps.python.org/pep-0735/) groups:
`test` (pytest), `dev` (the linters) and `docs` (the documentation site). A plain
`uv sync` installs only `dev`, so pytest would be missing.

## Check your change

Run the tests for the part you changed, one folder at a time while you work:

```bash
uv run pytest tests/unit/readers
uv run pytest tests/unit/path/to/test_file.py -k name
```

CI runs the whole unit suite on Linux and Windows. Tests marked `integration` need
real instrument data and are skipped by default.

Run the same lint check CI runs. It can block a release, so it must pass:

```bash
git add -A
uv run pre-commit run --all-files
```

`--all-files` skips untracked files, which is why `git add -A` comes first.

If you changed `docs/` or `mkdocs.yml`, run the strict build CI runs, and the
plain-language test:

```bash
uv sync --locked --group test --group dev --group docs
uv run mkdocs build --strict
uv run pytest tests/unit/test_docs_plain_language.py
```

## Title your pull request as a conventional commit

Pull requests are squash-merged, and the pull request title becomes the commit title
on `main`. The release tool reads it to decide the next version. A check fails the
pull request when the title does not follow the format:

```text
type(optional scope): short description
```

| type | effect |
|---|---|
| `feat` | new feature, next minor version |
| `fix`, `perf` | bug fix or speed-up, next patch version |
| `docs`, `test`, `ci`, `build`, `refactor`, `style`, `chore` | no release |

For example: `fix(resampling): refuse a bin count with a width`.

Add `!` after the type (`feat!: ...`) only for a change that breaks existing use. It
makes the next release a new major version, so ask a maintainer first.

Your commit messages become the body of that commit, and the release tool reads them
too. A paragraph that starts with a type, such as `fix: ...`, counts as one more
change, and a line starting `BREAKING CHANGE:` makes a major version. Write them in
the same format or as plain sentences. The maintainer who merges may rewrite the
body.

## Rules you cannot guess

- **No file paths or personal names in tracked files.** That covers code, comments,
  tests, docs and commit messages: no drive letters, no home folders, no lab shares,
  no colleagues' names. Pre-commit catches Windows drive-letter and home-folder
  paths, but not every form, so check your diff.
- **Name datasets, never path them.** Write "the PHI ToF-SIMS demo set", not where it
  sits on your disk.
- **No emojis** in code, comments, docs or commits.
- **User docs are plain language, and a test enforces it.** Pages under Home, "Get
  started" and "Guides" have a 25-word sentence limit and a word budget. Put the
  detail on a technical page and the reasons in
  [design decisions](docs/design-decisions.md). See "Writing the documentation" in
  the [full guide](docs/contributing.md#writing-the-documentation).
- **Some changes need a design-decision entry.** If your change alters what a
  converted store contains, changes a default, or retires an option, add an entry to
  [docs/design-decisions.md](docs/design-decisions.md). Take the next free D-number
  from `main`.
- **Those changes also need a release-notes line.** The generated notes list commit
  titles only. Write one or two plain sentences for users in the pull request
  description, under "Release notes".
- **Merging publishes nothing.** Releases are batched and cut on a schedule. See
  "Releases are batched" in the [full guide](docs/contributing.md#releases-are-batched-not-per-merge).

## Code of conduct

By taking part you agree to the [Code of Conduct](docs/code-of-conduct.md).
