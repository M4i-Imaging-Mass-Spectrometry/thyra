# tests/unit/readers/test_mass_axis_streaming.py
"""Tests for the streaming processed-mode common mass axis build.

``_extract_continuous_mass_axis`` used to collect every spectrum's m/z array
into a list, concatenate it, and hand the result to ``np.unique`` -- three
live copies of the whole m/z payload. It now streams: values go into a
scratch buffer that is sorted, deduplicated and merged into a running axis
whenever it fills.

The contract is that the result stays **bit-identical** to what
``np.unique(np.concatenate(all_mzs))`` returned, so nearly everything here
compares against that reference rather than against hand-written expectations.
``np.testing.assert_allclose``, which the older reader test uses, is not
sufficient: it cannot see a dtype change, a NaN difference, or a signed-zero
difference.

The most important thing in this file is that the whole matrix runs against a
monkeypatched batch size. At the real batch of 2^20 values every test input
fits in one batch, the merge path is never entered, and the tests are
worthless as a guard on the part most likely to break. Both prototypes of this
algorithm shipped with a crash that only a tiny batch size exposes: one on an
all-NaN batch, one on an empty accumulator.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from thyra.readers.imzml import imzml_reader as mod
from thyra.readers.imzml.imzml_reader import (
    ImzMLReader,
    _dedupe_sorted,
    _filter_absent,
    _merge_disjoint,
)

# Batch sizes small enough to force many folds, plus the real one.
BATCH_SIZES = [1, 2, 3, 7, 1 << 20]


class _StubParser:
    """Minimal stand-in for ImzMLParser over a list of m/z runs.

    Returns read-only arrays, as pyimzml does -- it hands back a
    ``np.frombuffer`` view over a ``bytes`` object. Anything that tried to sort
    a spectrum in place would raise here rather than silently corrupt the
    caller's data.
    """

    def __init__(self, runs, raise_on=()):
        self._runs = [np.asarray(r) for r in runs]
        self._raise_on = set(raise_on)
        self.coordinates = [(i + 1, 1, 1) for i in range(len(self._runs))]

    def getspectrum(self, idx):
        if idx in self._raise_on:
            raise RuntimeError(f"boom at {idx}")
        mzs = self._runs[idx]
        view = np.frombuffer(mzs.tobytes(), dtype=mzs.dtype)
        return view, np.ones(view.size, dtype=np.float32)


def _reference(runs):
    """What the previous implementation returned."""
    non_empty = [np.asarray(r) for r in runs if np.asarray(r).size]
    return np.unique(np.concatenate(non_empty))


def _build(runs, max_len=None, raise_on=()):
    reader = SimpleNamespace(max_mass_axis_length=max_len)
    return ImzMLReader._extract_continuous_mass_axis(
        reader, _StubParser(runs, raise_on=raise_on)
    )


def _assert_bit_identical(got, ref):
    """Compare exactly, including dtype, NaN payloads and signed zeros."""
    assert got.dtype == ref.dtype, f"dtype {got.dtype} != {ref.dtype}"
    assert got.shape == ref.shape, f"shape {got.shape} != {ref.shape}"
    raw = np.uint64 if got.dtype == np.float64 else np.uint32
    assert np.array_equal(got.view(raw), ref.view(raw))


@pytest.fixture(params=BATCH_SIZES)
def batch(request, monkeypatch):
    """Run each test at several batch sizes, including tiny ones."""
    monkeypatch.setattr(mod, "_MASS_AXIS_BATCH_VALUES", request.param)
    monkeypatch.setattr(mod, "_MASS_AXIS_PROBE_BLOCK", max(1, request.param))
    return request.param


class TestMatchesReference:
    """Bit-identity against np.unique(np.concatenate(...))."""

    @pytest.mark.parametrize(
        "runs",
        [
            pytest.param([[1.0, 2.0, 3.0]], id="single-spectrum"),
            pytest.param([[5.0]], id="single-peak"),
            pytest.param([[1.0], [2.0], [3.0], [4.0]], id="one-peak-each"),
            pytest.param([[1.0, 2.0], [1.0, 2.0], [1.0, 2.0]], id="all-identical"),
            pytest.param([[1.0, 1.0, 1.0, 2.0, 2.0, 3.0]], id="dupes-within"),
            pytest.param([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]], id="dupes-across"),
            pytest.param([[1.0, 2.0], [10.0, 11.0]], id="disjoint"),
            pytest.param(
                [list(np.arange(50.0, 100.0)), list(np.arange(0.0, 50.0))],
                id="second-entirely-below-first",
            ),
            pytest.param(
                [[float(20 - i)] for i in range(20)], id="strictly-descending"
            ),
            pytest.param([[3.0, 1.0, 2.0]], id="unsorted-within-spectrum"),
        ],
    )
    def test_shapes(self, runs, batch):
        _assert_bit_identical(_build(runs), _reference(runs))

    def test_many_spectra_shared_grid(self, batch):
        grid = np.arange(100.0, 200.0, 0.5)
        rng = np.random.default_rng(0)
        runs = [np.sort(rng.choice(grid, size=20, replace=False)) for _ in range(40)]
        _assert_bit_identical(_build(runs), _reference(runs))

    def test_many_spectra_all_distinct(self, batch):
        rng = np.random.default_rng(1)
        runs = [np.sort(rng.uniform(100.0, 200.0, 25)) for _ in range(40)]
        _assert_bit_identical(_build(runs), _reference(runs))

    def test_spectrum_larger_than_batch(self, batch):
        # Forces the reallocation branch when one spectrum exceeds the buffer.
        runs = [np.arange(0.0, 100.0), [500.0], np.arange(200.0, 400.0)]
        _assert_bit_identical(_build(runs), _reference(runs))


class TestEmptyAndDegenerate:
    """Empty spectra are skipped; nothing at all is an error."""

    @pytest.mark.parametrize(
        "runs",
        [
            pytest.param([[], [1.0, 2.0], [3.0]], id="leading-empty"),
            pytest.param([[1.0, 2.0], [3.0], []], id="trailing-empty"),
            pytest.param([[1.0], [], [2.0], [], [3.0]], id="interleaved-empty"),
        ],
    )
    def test_empty_spectra_are_skipped(self, runs, batch):
        _assert_bit_identical(_build(runs), _reference(runs))

    def test_all_spectra_empty_raises(self, batch):
        with pytest.raises(ValueError, match="No spectra found"):
            _build([[], [], []])

    def test_zero_spectra_raises(self, batch):
        with pytest.raises(ValueError, match="No spectra found"):
            _build([])


class TestFloatEdgeCases:
    """NaN, infinities and values that must not be collapsed."""

    @pytest.mark.parametrize(
        "runs",
        [
            pytest.param([[1.0, np.nan]], id="one-nan"),
            pytest.param([[np.nan]], id="only-nan"),
            pytest.param([[np.nan], [np.nan], [np.nan]], id="all-nan-every-batch"),
            pytest.param([[1.0, np.nan, np.nan, np.nan]], id="several-nan"),
            pytest.param([[1.0, np.nan], [2.0], [np.nan, 3.0]], id="nan-interleaved"),
        ],
    )
    def test_nan_matches_reference(self, runs, batch):
        got = _build(runs)
        _assert_bit_identical(got, _reference(runs))
        # numpy collapses NaNs and sorts them last.
        assert np.count_nonzero(np.isnan(got)) == 1
        assert np.isnan(got[-1])

    def test_all_nan_batch_does_not_crash(self, batch):
        """Regression: an all-NaN batch broke an early version of the dedupe.

        Without the ``first_nan > 0`` guard this raises
        "operands could not be broadcast together".
        """
        runs = [[np.nan, np.nan], [np.nan], [1.0]]
        _assert_bit_identical(_build(runs), _reference(runs))

    def test_infinities(self, batch):
        runs = [[-np.inf, 1.0], [np.inf], [0.0, np.inf]]
        _assert_bit_identical(_build(runs), _reference(runs))

    def test_one_ulp_apart_is_not_collapsed(self, batch):
        a = 700.0
        b = np.nextafter(a, np.inf)
        got = _build([[a], [b]])
        assert got.size == 2
        _assert_bit_identical(got, _reference([[a], [b]]))

    def test_subnormals(self, batch):
        tiny = np.finfo(np.float64).tiny
        runs = [[tiny, tiny / 2.0], [0.0]]
        _assert_bit_identical(_build(runs), _reference(runs))


class TestDtype:
    """The axis dtype follows the file's mzPrecision."""

    def test_float32_stays_float32(self, batch):
        runs = [np.array([3.0, 1.0], dtype=np.float32), np.array([2.0], np.float32)]
        got = _build(runs)
        assert got.dtype == np.float32
        _assert_bit_identical(got, _reference(runs))

    def test_float64_stays_float64(self, batch):
        runs = [np.array([3.0, 1.0], dtype=np.float64)]
        assert _build(runs).dtype == np.float64


class TestMaxMassAxisLength:
    """The guard that stops a raw axis growing without bound."""

    def test_raises_when_exceeded(self, batch):
        runs = [[float(i)] for i in range(100)]
        with pytest.raises(ValueError, match="exceeded"):
            _build(runs, max_len=10)

    def test_message_names_the_limit_and_the_count(self, batch):
        runs = [[float(i)] for i in range(100)]
        with pytest.raises(ValueError) as exc:
            _build(runs, max_len=10)
        text = str(exc.value)
        assert "10" in text
        assert "spectra" in text
        # It must say what to do about it.
        assert "resampling" in text

    def test_not_swallowed_by_the_per_spectrum_handler(self, batch):
        """Regression: the guard must escape the read's try/except.

        The previous code wrapped the whole loop body, so a ValueError raised
        during the merge would have been logged as a warning and ignored.
        """
        runs = [[float(i)] for i in range(100)]
        with pytest.raises(ValueError):
            _build(runs, max_len=5)

    def test_generous_limit_does_not_raise(self, batch):
        runs = [[float(i)] for i in range(50)]
        assert _build(runs, max_len=10_000).size == 50

    def test_explicit_none_is_unlimited(self, batch):
        runs = [[float(i)] for i in range(50)]
        assert _build(runs, max_len=None).size == 50


class TestDefaultMaxMassAxisLength:
    """The reader's default cap, and how a caller overrides it.

    SCiLS Lab limits its own common mass axis to 10 million bins (2026b User
    Guide, p.76). Thyra adopts that number: it is the tool the resampling was
    designed against, and previously the cap defaulted to None only because
    no defensible figure was available.
    """

    def _reader(self, **options):
        """Construct a reader without touching the filesystem.

        ``ImzMLReader.__init__`` is lazy -- it stores paths and defers the
        parser -- so it is safe to build one for a path that does not exist.
        """
        return ImzMLReader(Path("does-not-exist.imzML"), **options)

    def test_default_is_ten_million(self):
        assert mod.DEFAULT_MAX_MASS_AXIS_LENGTH == 10_000_000
        assert self._reader().max_mass_axis_length == 10_000_000

    def test_caller_can_lower_it(self):
        assert self._reader(max_mass_axis_length=5_000).max_mass_axis_length == 5_000

    def test_caller_can_raise_it(self):
        reader = self._reader(max_mass_axis_length=50_000_000)
        assert reader.max_mass_axis_length == 50_000_000

    def test_explicit_none_restores_unlimited(self):
        """Absent means "use the default"; None means "no cap at all"."""
        assert self._reader(max_mass_axis_length=None).max_mass_axis_length is None


class TestPerSpectrumErrors:
    """A failing spectrum is warned about, not fatal."""

    def test_failing_spectra_are_skipped(self, batch, caplog):
        runs = [[1.0], [2.0], [3.0], [4.0]]
        got = _build(runs, raise_on=(1, 2))
        assert np.array_equal(got, np.array([1.0, 4.0]))

    def test_all_spectra_failing_raises(self, batch):
        with pytest.raises(ValueError, match="No spectra found"):
            _build([[1.0], [2.0]], raise_on=(0, 1))


class TestFuzz:
    """Randomised comparison against the reference.

    This is what found the crashes in both prototypes of this algorithm.
    """

    @pytest.mark.parametrize("seed", range(60))
    def test_random_runs_match_reference(self, seed, monkeypatch):
        rng = np.random.default_rng(seed)
        monkeypatch.setattr(mod, "_MASS_AXIS_BATCH_VALUES", int(rng.integers(1, 16)))
        monkeypatch.setattr(mod, "_MASS_AXIS_PROBE_BLOCK", int(rng.integers(1, 16)))

        n_runs = int(rng.integers(0, 12))
        pool = rng.uniform(0.0, 50.0, 25)
        runs = []
        for _ in range(n_runs):
            size = int(rng.integers(0, 10))
            vals = rng.choice(pool, size=size, replace=True).astype(np.float64)
            if size and rng.random() < 0.15:
                vals[rng.integers(0, size)] = np.nan
            if size and rng.random() < 0.08:
                vals[rng.integers(0, size)] = np.inf
            runs.append(np.sort(vals))

        if not any(r.size for r in runs):
            with pytest.raises(ValueError):
                _build(runs)
            return

        _assert_bit_identical(_build(runs), _reference(runs))


class TestHelpers:
    """The pure helpers, targeted directly."""

    @pytest.mark.parametrize(
        "values",
        [
            [],
            [1.0],
            [1.0, 1.0],
            [1.0, 2.0, 2.0, 3.0],
            [np.nan],
            [np.nan, np.nan],
            [1.0, np.nan, np.nan],
        ],
    )
    def test_dedupe_sorted_matches_np_unique(self, values):
        a = np.array(values, dtype=np.float64)
        got = _dedupe_sorted(a)
        ref = np.unique(a)
        assert np.array_equal(got, ref, equal_nan=True)

    def test_dedupe_sorted_does_not_mutate_input(self):
        a = np.array([1.0, 1.0, 2.0])
        before = a.copy()
        _dedupe_sorted(a)
        assert np.array_equal(a, before)

    def test_filter_absent_empty_accumulator(self):
        """Regression: an empty accumulator used to raise IndexError."""
        s = np.array([1.0, 2.0])
        assert np.array_equal(_filter_absent(np.array([]), s), s)

    def test_filter_absent_removes_present_values(self):
        acc = np.array([1.0, 3.0, 5.0])
        s = np.array([2.0, 3.0, 6.0])
        assert np.array_equal(_filter_absent(acc, s), np.array([2.0, 6.0]))

    def test_filter_absent_all_present(self):
        acc = np.array([1.0, 2.0, 3.0])
        assert _filter_absent(acc, np.array([1.0, 3.0])).size == 0

    def test_filter_absent_nan_already_present(self):
        """A NaN in the accumulator must not be re-added every fold."""
        acc = np.array([1.0, np.nan])
        assert _filter_absent(acc, np.array([np.nan])).size == 0

    def test_merge_disjoint_unions_and_sorts(self):
        out = _merge_disjoint([np.array([1.0, 4.0])], [np.array([2.0, 3.0])])
        assert np.array_equal(out, np.array([1.0, 2.0, 3.0, 4.0]))

    def test_merge_disjoint_releases_its_inputs(self):
        """The boxes are emptied so the caller's reference cannot pin memory."""
        box_a, box_b = [np.array([1.0])], [np.array([2.0])]
        _merge_disjoint(box_a, box_b)
        assert box_a[0] is None
        assert box_b[0] is None
