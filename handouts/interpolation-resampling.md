# F. PCHIP interpolation: designed, half-built, never shipped

**Branch:** `feat/pchip-resampling` (not created)
**Worktree:** `../Thyra-interp`
**Priority:** backlog. Nothing is broken; this is a capability gap.

This handout exists because the work behind it was about to be lost. A
2026-07-30 sweep deleted 38 stale local branches whose content was already on
`main`, and 18 more that had genuinely unique commits. Of those 18, exactly
one held something `main` still does not have. Everything else was
superseded, reimplemented, or abandoned.

The code is preserved at tag **`archive/conservative-interpolation`**
(`106d6ff`). Tags are refs, so it will survive `git gc`; the branch it came
from is gone.

---

## The gap

`main` offers three resampling methods and no interpolation:

```python
class ResamplingMethod(Enum):
    NONE = "none"
    NEAREST_NEIGHBOR = "nearest_neighbor"
    TIC_PRESERVING = "tic_preserving"
```

`thyra/resampling/strategies/` holds `nearest_neighbor.py` and
`tic_preserving.py`. There is no interpolating strategy, and
`git grep -i pchip thyra/` is empty.

That is a real limitation for profile-mode data. Nearest-neighbour snaps each
peak to a bin and distorts peak shape; TIC-preserving redistributes intensity
so the total ion count survives rebinning, which is the right call for
quantitative work but still does not reconstruct the profile between samples.
Shape-preserving cubic interpolation is the third option, and it is the one
mMass-style pipelines reach for when the m/z grid changes.

Whether that matters enough to build is the open question. It is a judgement
about the science, not about the code, which is why this is a handout rather
than a ticket.

## What already shipped, so you do not rebuild it

The archived branch carries `INTERPOLATION_PLAN.md`, a two-phase plan.
**Phase 1 shipped**, reimplemented and improved, and is on `main` today:

| Plan item | Status on `main` |
|---|---|
| `BaseMassAxisGenerator` abstract class | `thyra/resampling/mass_axis/base_generator.py` |
| `LinearMassAxisGenerator` | `mass_axis/linear_generator.py` |
| Physics-based TOF / Orbitrap / FTICR generators | `mass_axis/linear_tof_generator.py`, `reflector_tof_generator.py`, `orbitrap_generator.py`, `fticr_generator.py` |
| Auto-select generator from instrument metadata | `thyra/resampling/instrument_detectors.py` + `decision_tree.py` |

So the axis half of the plan is done, and done better than the archived
branch did it -- the branch only had a linear generator.

**Phase 2 did not ship.** The archived branch has a working first cut:

```
msiconvert/interpolators/base_interpolator.py    172 lines
msiconvert/interpolators/pchip_interpolator.py   125 lines
tests/unit/interpolators/test_pchip_interpolator.py    281 lines
tests/unit/interpolators/test_integration_phase2.py    296 lines
```

## Read it, do not merge it

The archived code is from **August 2025** and cannot be applied:

- It targets the `msiconvert/` package, before the rename to `thyra/`.
- It predates `thyra/resampling/` entirely, so its own mass-axis generators
  duplicate what `main` now has and its `BaseInterpolator` does not implement
  the `strategies/base.py` interface a resampling strategy has to satisfy
  today.
- Its `poetry.lock` is a year stale.

Treat it as a design document with tests attached, not as a patch. The
useful parts are the plan's edge-case list (single points, duplicate m/z
values, NaN handling) and the test cases, which encode what "correct"
looked like to whoever wrote them.

## What to do, if it is worth doing

1. Decide the science question first: does profile-mode conversion need
   shape-preserving interpolation, given `TIC_PRESERVING` already exists?
   If the answer is no, delete the tag and this handout.
2. If yes, add `PCHIP = "pchip"` to `ResamplingMethod` and implement it as a
   normal strategy in `thyra/resampling/strategies/pchip.py` against
   `strategies/base.py`. Do not resurrect the `interpolators/` package
   layout -- `main`'s resampling module is the architecture now.
3. Port the archived tests to the current interfaces rather than writing
   fresh ones. They are the most valuable thing in the tag.
4. Wire it into `ResamplingDecisionTree` only if there is a rule for when it
   should be chosen automatically. Otherwise leave it opt-in.

## Constraints

- Do not add an `interpolators/` package. `thyra/resampling/` is the home for
  anything that changes the mass axis, and splitting that across two package
  trees is what made the 2025 attempt unmergeable.
- `scipy` is already a dependency (`>=1.7.0`), so `scipy.interpolate.PchipInterpolator`
  needs no new dependency. Do not add one.
- Interpolation must stay opt-in. `NONE` is the default for a reason: the
  cheapest conversion is the one that does not touch the data.

## Verification

```bash
uv run pytest -q
uv run black . && uv run isort . && uv run flake8
```

A new strategy needs a round-trip test proving that interpolating onto the
*same* axis is the identity, and a TIC-drift test quantifying how much total
ion count the method loses relative to `TIC_PRESERVING`. Without the second
one there is no way to advise users which to pick.
