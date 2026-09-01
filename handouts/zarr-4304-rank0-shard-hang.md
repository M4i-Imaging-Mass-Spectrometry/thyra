# Handout J -- zarr rank-0 shard hang (zarr-python#4304)

Upstream work in `zarr-developers/zarr-python`, not in this repository. Found while designing
handout I.

Status as of 2026-09-01:

- Filed: https://github.com/zarr-developers/zarr-python/issues/4304 (open)
- Fix not written. No fork of zarr-python exists on this account yet.
- Reproduced first-hand on zarr 3.2.1 and 3.3.0; the loop is unchanged on `main`.
- NOT a duplicate of zarr-python#3751, which was an `IndexError` on an explicit `shards=()`
  write, fixed by #3793.

Why it matters here: it is the reason handout I's design must never put `shards` into
anndata's `dataset_kwargs`. Every SpatialData table carries rank-0 scalars in
`uns/spatialdata_attrs`, so a `shards` key reaching them would hang the write rather than
raise. Thyra is not affected today, because anndata's scalar writers do not route through the
sharding path.

# PART 1 -- THE ISSUE AS FILED

**Zarr version:** 3.3.0 (also on 3.2.1, and the loop is unchanged on `main`)

**Numcodecs version:** 0.16.5

**Python Version:** 3.13.3

**Operating System:** Windows

**Installation:** pip into a fresh venv

---

## Description

`shards="auto"` on a 0-d array works fine and gives you `shards=()`. Set
`array.target_shard_size_bytes` (the key added in #3568) and the same call hangs forever. No error,
no warning, nothing printed, it just spins in pure Python until you kill it.

So a config key that's only meant to affect shard size ends up deciding whether the call returns at
all.

`_guess_num_chunks_per_axis_shard` in `src/zarr/core/chunk_grids.py` only gets called when
`shards="auto"` and the key is set. Without the key the caller just uses 2, so you never reach it.
For a 0-d array both loop conditions are constant:

```python
bytes_per_chunk = math.prod(chunk_shape) * item_size   # math.prod(()) == 1, so == item_size
if max_bytes < bytes_per_chunk:                        # won't fire for any budget >= item_size
    return 1
num_axes = len(chunk_shape)                            # 0
chunks_per_shard = 1
while (bytes_per_chunk * ((chunks_per_shard + 1) ** num_axes)) <= max_bytes and all(
    c * (chunks_per_shard + 1) <= a for c, a in zip(chunk_shape, array_shape, strict=True)
):
    chunks_per_shard += 1
```

`(chunks_per_shard + 1) ** 0` is 1, so the first check is just `item_size <= max_bytes` and never
changes. `zip((), ())` is empty so `all()` is True. Neither depends on the counter, so it loops
forever.

Should return 1. A 0-d array has no axes to grow a shard along, which is what the no-budget path
already works out on its own when it gives you `shards=()`.

Fix is right after `num_axes = len(chunk_shape)`:

```python
if num_axes == 0:
    return 1
```

## Steps to reproduce

```python
# /// script
# requires-python = ">=3.11"
# dependencies = ["zarr==3.3.0"]
# ///
"""`shards="auto"` on a 0-d array never terminates when `array.target_shard_size_bytes` is set.

Each case runs in a subprocess with a 15 s timeout, so this script always exits.
"""

import subprocess
import sys
import warnings

CASES = {"none": None, "128MiB": 128 * 1024 * 1024}

if len(sys.argv) > 1:  # child: run one case
    import zarr

    warnings.simplefilter("ignore")  # "shard shape inference is experimental"
    budget = CASES[sys.argv[1]]
    with zarr.config.set({"array.target_shard_size_bytes": budget} if budget else {}):
        arr = zarr.create_array(store={}, shape=(), dtype="int64", shards="auto")
    print(f"returned, shards={arr.shards}")
else:  # parent: run both, bounded
    import zarr

    print(f"zarr {zarr.__version__}\n")
    for name in CASES:
        print(f"target_shard_size_bytes = {name}:")
        try:
            out = subprocess.run(
                [sys.executable, __file__, name], capture_output=True, text=True, timeout=15
            )
            print(f"    {out.stdout.strip() or out.stderr.strip().splitlines()[-1]}")
        except subprocess.TimeoutExpired:
            print("    DID NOT TERMINATE after 15 s")
```

Gives:

```
zarr 3.3.0

target_shard_size_bytes = none:
    returned, shards=()
target_shard_size_bytes = 128MiB:
    DID NOT TERMINATE after 15 s
```

Calling it directly shows the same thing, and that the early return still works:

```python
from zarr.core.chunk_grids import _guess_num_chunks_per_axis_shard

_guess_num_chunks_per_axis_shard(chunk_shape=(), item_size=8, max_bytes=4, array_shape=())
# 1, early return fires

_guess_num_chunks_per_axis_shard(chunk_shape=(), item_size=8, max_bytes=128 * 1024 * 1024, array_shape=())
# never returns
```

## Notes

Not a dup of #3751 ("Cannot write a 0 dimensional array with sharding", closed). That one was an
`IndexError` out of `_ShardIndex.get_chunk_slices_vectorized` on `arr[()] = 42.0` with an explicit
`shards=()`, fixed by #3793. This is array creation under `shards="auto"`, and it hangs instead of
raising.

Where I ran into it: anndata writes zarr v3 arrays with `shards="auto"` and reads the same config
key, so anything that sets a shard budget and then writes a 0-d array ends up here. anndata dodges
it at the moment because its scalar writers don't go through the sharding path, so nobody's actually
broken by this today. You can still get there from plain zarr, as above.

#4218 is refactoring `chunk_grids.py` but doesn't touch this loop.

# PART 2 -- THE HANDOFF PROMPT

`````
Fix zarr-developers/zarr-python issue #4304 and open a PR.
https://github.com/zarr-developers/zarr-python/issues/4304

## The bug

`_guess_num_chunks_per_axis_shard` in `src/zarr/core/chunk_grids.py` never terminates for a
0-dimensional array. It is only reached when `shards="auto"` and the zarr config key
`array.target_shard_size_bytes` is set (otherwise the caller uses the constant 2 and never calls it).

```python
bytes_per_chunk = math.prod(chunk_shape) * item_size   # math.prod(()) == 1, so == item_size
if max_bytes < bytes_per_chunk:                        # won't fire for any budget >= item_size
    return 1
num_axes = len(chunk_shape)                            # 0 for a 0-d array
chunks_per_shard = 1
while (bytes_per_chunk * ((chunks_per_shard + 1) ** num_axes)) <= max_bytes and all(
    c * (chunks_per_shard + 1) <= a for c, a in zip(chunk_shape, array_shape, strict=True)
):
    chunks_per_shard += 1
```

`(chunks_per_shard + 1) ** 0` is 1, so the first condition is the constant `item_size <= max_bytes`.
`zip((), ())` is empty so `all()` is vacuously True. Neither depends on the counter, so it spins
forever with no error and no output.

Verified on zarr 3.2.1 and 3.3.0, and the loop is unchanged on `main`.

## The fix

Add a guard immediately after `num_axes = len(chunk_shape)`:

```python
if num_axes == 0:
    return 1
```

Returning 1 is correct and safe: the caller then does
`for a_shape, c_shape in zip(array_shape, chunk_shape_flat, strict=True)`, which is empty at rank 0,
so `_shards_out` stays `()` whatever this function returns. The result is `shards=()`, which is
exactly what the no-budget path already produces today. The fix makes the budget and no-budget paths
agree.

Do not refactor the loop, do not "fix" other rank-0 edge cases you notice, do not touch anything
else in the file. Two lines plus tests plus a changelog fragment. That is the whole PR.

## Steps

1. `gh repo fork zarr-developers/zarr-python --clone` (there is no existing fork on this account and
   no local clone). Work on a branch off `main`, e.g. `fix/4304-rank0-shard-hang`.
2. Before changing anything, reproduce it yourself so you know the fix works. Run this in a
   subprocess with a timeout so nothing runs away:

```python
import subprocess, sys, warnings
CASES = {"none": None, "128MiB": 128 * 1024 * 1024}
if len(sys.argv) > 1:
    import zarr
    warnings.simplefilter("ignore")
    budget = CASES[sys.argv[1]]
    with zarr.config.set({"array.target_shard_size_bytes": budget} if budget else {}):
        arr = zarr.create_array(store={}, shape=(), dtype="int64", shards="auto")
    print(f"returned, shards={arr.shards}")
else:
    for name in CASES:
        try:
            out = subprocess.run([sys.executable, __file__, name], capture_output=True,
                                 text=True, timeout=15)
            print(name, "->", out.stdout.strip())
        except subprocess.TimeoutExpired:
            print(name, "-> DID NOT TERMINATE")
```

   Expected before the fix: `none -> returned, shards=()` and `128MiB -> DID NOT TERMINATE`.
   Expected after: both return `shards=()`.

3. Apply the guard.
4. Add tests to `tests/test_chunk_grids.py`. Two of them:

```python
def test_guess_num_chunks_per_axis_shard_0d():
    # regression test for https://github.com/zarr-developers/zarr-python/issues/4304
    assert (
        _guess_num_chunks_per_axis_shard(
            chunk_shape=(), item_size=8, max_bytes=128 * 1024 * 1024, array_shape=()
        )
        == 1
    )


def test_create_0d_array_auto_shards_with_target_shard_size():
    with zarr.config.set({"array.target_shard_size_bytes": 128 * 1024 * 1024}):
        arr = zarr.create_array(store={}, shape=(), dtype="int64", shards="auto")
    assert arr.shards == ()
```

   Match the file's existing import and naming style rather than copying these verbatim. The second
   test is the user-facing one and matters more than the first.

   Note on test design: the repo has no `pytest-timeout`. It sets `faulthandler_timeout = 600` with
   `faulthandler_exit_on_timeout = true`, so if someone reintroduces the bug these tests stall for
   10 minutes and then get killed with a traceback pointing at the loop, rather than failing fast.
   That is acceptable, do not build a threading harness around it. Mention it in the PR description
   so a reviewer is not surprised.

5. Add a towncrier changelog fragment: `changes/<PR-number>.bugfix.md`, one or two sentences. The
   config is `[tool.towncrier]` in `pyproject.toml` with `directory = 'changes'` and
   `issue_format` pointing at `/pull/{issue}`, so name it after the PR number once you have it.
   Some existing fragments use the issue number instead, so either is defensible, prefer the PR
   number. Look at a couple of existing `changes/*.bugfix.md` files for tone.

6. Run the tests. The repo uses hatch (`[tool.hatch.envs.test]`, matrix over Python 3.12 / 3.13 /
   3.14, `minimal` and `optional` dep sets). Running `pytest tests/test_chunk_grids.py` in a venv
   with the package installed is enough for local confidence; CI is the referee.

7. Lint: ruff 0.16.3, `line-length = 100`, flake8-annotations is on but `tests/**` is exempt from
   ANN001 and ANN201. Run the repo's pre-commit if it is configured.

8. Open the PR against `zarr-developers/zarr-python:main`. Link it to #4304 with a closing keyword.
   In the description: state the one line cause, the fix, and that budget and no-budget now agree at
   rank 0. Add one sentence offering to move the guard if maintainers would rather it live in the
   caller (`resolve_outer_and_inner_chunks`) instead of inside the function. Keep the whole thing
   short.

## Style

No emojis anywhere, including commit messages and the PR body. No em dashes. Plain, short sentences.
Do not add Claude or any AI tool as co-author or committer, and do not add "Assisted-by" or
"Generated with" trailers. Commit as the repo's configured git user only.

## Definition of done

- Guard added, both tests pass, and the reproducer from step 2 returns `shards=()` for both cases.
- Changelog fragment present.
- Lint clean.
- PR open against upstream, linked to #4304, CI green.
- Nothing else in the diff. Expected shape: 1 source file, 1 test file, 1 changelog file.

## Known context

- PR #4218 is currently refactoring `chunk_grids.py` heavily but does not touch this loop. If it
  lands first, expect a trivial rebase. Do not try to coordinate with it.
- Issue #3751 ("Cannot write a 0 dimensional array with sharding") looks similar and is NOT this
  bug. It was an `IndexError` from `_ShardIndex.get_chunk_slices_vectorized` on `arr[()] = 42.0`
  with an explicit `shards=()`, fixed by #3793. Do not conflate them, and do not claim this PR fixes
  #3751.
`````
