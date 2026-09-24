# Contributing to Thyra

Thank you for your interest in contributing to Thyra! This document provides guidelines for contributing to this project.

## Table of Contents

- [Development Environment Setup](#development-environment-setup)
- [Code Style Guidelines](#code-style-guidelines)
- [Testing Requirements](#testing-requirements)
- [Pull Request Process](#pull-request-process)
- [Issue Reporting Guidelines](#issue-reporting-guidelines)
- [Communication Channels](#communication-channels)
- [Writing the documentation](#writing-the-documentation)

## Development Environment Setup

### Prerequisites

- Python 3.12 or 3.13
- [uv](https://docs.astral.sh/uv/) for dependency management
- Git

### Setup Steps

1. **Fork and Clone the Repository**
   ```bash
   git clone https://github.com/<your-username>/thyra.git
   cd thyra
   ```

2. **Install Dependencies**
   ```bash
   uv sync
   ```

3. **Install Pre-commit Hooks**
   ```bash
   uv run pre-commit install
   ```

   CI runs `pre-commit run --all-files` on every pull request, so skipping this
   step does not skip the checks. It only moves the failure from your machine
   to the PR.

4. **Verify Installation**
   ```bash
   uv run pytest -m "not integration"
   ```

## Code Style Guidelines

We use automated tools to maintain consistent code style:

### Formatting
- **Black** for code formatting
- **isort** for import sorting
- Line length: 88 characters (Black default)

### Linting
- **flake8** for code linting
- **bandit** for security checks

### Type Hints
- Use type hints for all public functions
- Follow PEP 484 conventions
- `thyra/py.typed` is the PEP 561 marker, so those annotations are part of the
  published API rather than an internal convenience: a type checker in a
  downstream project reads them instead of resolving every thyra symbol to
  `Any`. The marker and the `Typing :: Typed` classifier in `pyproject.toml`
  stand or fall together, and the `clean-venv-install` CI job asserts the
  installed distribution carries it.
- `disallow_untyped_defs` and `disallow_incomplete_defs` are still off under
  `[tool.mypy]`, so a function left unannotated is not caught here -- it is
  published as implicit `Any` inside a package that advertises itself as
  typed. That is what makes the first bullet a rule rather than a preference.

### Running Code Quality Checks

One command reproduces the CI lint job exactly, because CI runs this command:

```bash
uv run pre-commit run --all-files
```

It applies black, isort, flake8, mypy, bandit and pydocstyle with the settings
in `.flake8` and `pyproject.toml`, plus the file-hygiene and lab-share-path
hooks. **`--all-files` skips untracked files**, so `git add -A` first or a new
file you have just written is not checked.

The individual tools are still there if you want to run one on its own:

```bash
# Format code
uv run black .
uv run isort .

# Run linting
uv run flake8

# Run security checks
uv run bandit -r thyra/
```

## Testing Requirements

### Test Types

1. **Unit Tests** - Fast tests for individual functions
   ```bash
   uv run pytest -m "not integration"
   ```

2. **Integration Tests** - End-to-end workflow tests
   ```bash
   uv run pytest -m "integration"
   ```

### Test Coverage

- Aim for >80% code coverage for new code
- Run tests with coverage:
  ```bash
  uv run pytest --cov=thyra --cov-report=html
  ```

### Writing Tests

- Place unit tests in `tests/unit/`
- Place integration tests in `tests/integration/`
- Use descriptive test names: `test_should_convert_imzml_when_valid_file_provided`
- Do not hand-write a `unit` or `integration` marker. `tests/conftest.py`
  stamps one on every collected test from the directory it lives in, so the
  directory you chose above is the whole of the decision
- Capture Thyra's log records with the `thyra_logs` fixture, not `caplog`.
  `setup_logging` sets `propagate = False` on the `thyra` logger, so once
  any test has invoked the CLI caplog's root handler stops seeing Thyra
  records -- and an assertion against an empty capture passes rather than
  fails. `thyra_logs` attaches to the named logger and carries `.text`,
  `.messages` and the records themselves
- Set command-line arguments with `monkeypatch.setattr(sys, "argv", [...])`,
  never by assigning `sys.argv`. An assignment is never undone, so the last
  CLI test decides what every test after it sees

Both rules are enforced by `tests/unit/test_log_capture_convention.py`,
which reads the test sources rather than running them.

## Pull Request Process

### Before Submitting

1. **Create a Feature Branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make Your Changes**
   - Follow code style guidelines
   - Add tests for new functionality
   - Update documentation if needed

3. **Run Quality Checks**
   ```bash
   git add -A                            # --all-files skips untracked files
   uv run pre-commit run --all-files
   uv run pytest
   ```

   The last two are what the CI `lint` and `test` jobs run. Running them here
   is the only way to see a red check before a reviewer does.

4. **Commit Your Changes**
   - Use clear, descriptive commit messages
   - Follow conventional commit format when possible:
     ```
     feat: add dry-run mode for conversion preview
     fix: resolve memory leak in large dataset processing
     docs: update installation instructions
     ```

### Submitting the Pull Request

1. **Push Your Branch**
   ```bash
   git push origin feature/your-feature-name
   ```

2. **Open a Pull Request**
   - Use the provided PR template
   - Fill out all sections completely
   - Link related issues

3. **Respond to Review Feedback**
   - Address all reviewer comments
   - Make requested changes promptly
   - Ask questions if feedback is unclear

## Issue Reporting Guidelines

### Before Opening an Issue

1. **Check Existing Issues** - Search for similar issues first
2. **Reproduce the Problem** - Ensure you can consistently reproduce the issue
3. **Gather Information** - Collect relevant system info, error messages, and sample data (if shareable)

### Issue Types

Use the appropriate issue template:

- **Bug Report** - For reporting software defects
- **Feature Request** - For suggesting new functionality
- **Question** - For general questions about usage

### Information to Include

**For Bug Reports:**
- Operating system and version
- Python version
- Thyra version
- Input file format and size (if relevant)
- Complete error message and stack trace
- Steps to reproduce

**For Feature Requests:**
- Clear description of the desired functionality
- Use case and motivation
- Proposed implementation approach (if you have ideas)

## Communication Channels

### Getting Help

- **GitHub Issues** - For bug reports and feature requests
- **GitHub Discussions** - For general questions and community discussions

## Writing the documentation

The documentation source is in the `docs/` folder. Many readers are
scientists who have never used a terminal, so the site has layers, and every
topic should have all of them:

1. **Beginner pages**, under "Get started" and "Guides". They say what to do
   and what the reader gets, in plain words.
2. **Advanced boxes** inside those pages, for optional depth. They are always
   collapsed and their title starts with "Advanced:", so a beginner sees one
   line and knows it is optional:

    ```markdown
    ??? advanced "Advanced: How Thyra chooses the shared m/z axis"
        The explanation, indented four spaces.
    ```

3. **Technical reference** pages, for readers who want every option, every
   format and every detail.

### Rules for beginner pages

- Keep sentences to 15 to 20 words on average, and split any over 25.
- Give each paragraph one idea, in three to five sentences at most.
- Address the reader as "you", in the present tense and the active voice.
- Put one action in each numbered step. Show the command, then what the
  reader will see.
- Define a term the first time it appears, or link it to the
  [Glossary](glossary.md).
- Leave out version history ("used to", "since 3.x"), issue numbers and
  measurements. History belongs in the changelog. Reasons and measurements
  belong in [Design decisions](design-decisions.md), linked in one line.
- Cut before you add. A change that makes a beginner page longer needs a
  reason.

When a code change needs a documentation change, write the shortest
sentence the beginner page needs. Put the detail on the technical page and
the reasons in Design decisions.

`docs/includes/abbreviations.md` lists terms that show their meaning on
hover, on every page. Add a term there when pages use it without explaining
it, but only if a mass spectrometrist might not know it: computing words, or
terms from another part of the field. Every reader knows m/z and TIC, so
underlining them would only add clutter.

Test every code example, and preview the site before opening a pull request:

```bash
uv sync --group docs
uv run mkdocs serve
```

## Development Workflow

### Typical Contribution Flow

1. **Choose an Issue**
   - Look for issues labeled `good-first-issue` if you're new
   - Comment on the issue to indicate you're working on it

2. **Develop Your Solution**
   - Follow the development environment setup
   - Make small, focused commits
   - Write tests as you go

3. **Test Thoroughly**
   - Run the full test suite
   - Test with different input formats if relevant
   - Verify performance impact for large datasets

4. **Document Your Changes**
   - Update docstrings for modified functions
   - Update user documentation if needed
   - Add changelog entry for significant changes

## Code of Conduct

Please note that this project is governed by our [Code of Conduct](code-of-conduct.md). By participating, you agree to abide by its terms.

## Versioning Policy

Thyra follows [Semantic Versioning](https://semver.org/) (SemVer):

### Version Format: `MAJOR.MINOR.PATCH`

- **MAJOR** - Incremented for incompatible API changes
- **MINOR** - Incremented for backwards-compatible functionality additions
- **PATCH** - Incremented for backwards-compatible bug fixes

### Release Process

1. **Automated Versioning** - We use `python-semantic-release` for automated version bumping
2. **Commit Message Format** - Use conventional commits to trigger appropriate version bumps:
   ```
   feat: add new converter format (triggers MINOR)
   fix: resolve memory leak (triggers PATCH)
   feat!: redesign API structure (triggers MAJOR)
   ```

3. **Breaking Changes** - Always include `!` in commit type or `BREAKING CHANGE:` in footer
4. **Changelog** - Automatically generated from commit messages into
   `CHANGELOG.md`. Never edit it by hand; the [Changelog](changelog.md) page
   includes that file verbatim.

### Releases are batched, not per-merge

**Merging a pull request does not publish a release.** Releases are cut on a
schedule or on demand, and each one covers every commit merged since the
previous tag.

`semantic-release version` reads all commits since the last tag and applies the
highest bump among them, so a batch of eight `fix:` merges becomes one patch
version whose changelog section lists all eight. Releasing per merge would have
turned the same eight into eight versions and eight PyPI uploads.

`.github/workflows/release.yml` therefore has no `push` trigger. It runs:

- **on a cron**, Mondays at 06:00 UTC, releasing whatever has accumulated;
- **on demand**, whenever you want a release sooner:

  ```bash
  gh workflow run release.yml
  ```

- **publish-only**, to re-upload the current version to PyPI without cutting a
  new one, via the `publish_only` input on the Run workflow button.

A run with nothing releasable since the last tag - only `chore:`, `ci:`,
`docs:`, `refactor:`, `style:`, `test:` or `build:` commits - reports
`No release` and exits green. Only `feat:` (minor) and `fix:`/`perf:` (patch)
move the version.

Nothing is lost by waiting: unreleased commits sit on `main` and the next run
picks them up. What you should *not* do is merge a fix and then expect
`pip install thyra` to have it minutes later - check the tags, or trigger a
release yourself.

### Grouping issues into release batches

Open issues are grouped into **milestones**, one per planned batch. A milestone
is the unit of release: work through its issues on one branch, open one pull
request that closes all of them, merge, then cut a release.

Milestones rather than draft pull requests, because a draft PR needs a branch
with commits on it - six placeholder branches would each run CI, go stale
against a moving `main`, and say nothing a milestone does not. Open the pull
request when you start writing the code, not when you plan the batch.

Group by the code the issues touch, not by how they were found. Issues that
edit the same reader or the same function belong in one batch: fixing them in
separate pull requests means solving the same merge conflict once per pull
request.

### Development Versions

- **Alpha/Beta** releases may be created for testing: `1.2.0-alpha.1`
- **Release candidates** before major releases: `2.0.0-rc.1`

## Questions?

If you have questions about contributing that aren't covered here, please:

1. Check the existing [GitHub Discussions](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/discussions)
2. Open a new discussion if your question hasn't been asked
3. Tag maintainers if you need urgent clarification

Thank you for contributing to Thyra!
