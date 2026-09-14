"""Five readers, one raw mass-axis builder (issue #294).

Each of these built "the union of every spectrum's m/z" its own way, and
three of the five did it badly:

- **Waters** and **solariX** appended every spectrum's array to a list and
  then paid ``np.concatenate`` + ``np.unique`` -- three payload-sized
  allocations, with a transient of O(*total peaks*) rather than of the
  distinct values it was computing.
- **mzPeak** folded with ``np.union1d`` per row group, re-copying the whole
  axis every group.
- **timsTOF** accumulated into a ``set[float]``: O(unique) in elements but
  ~10x that in bytes, and it reported its own cost as ``len(unique) * 8``.

What matters is that the *answer* is unchanged. Every one of these asserts
the reader returns exactly what ``np.unique`` over its inputs returns, which
is the contract every downstream consumer was already written against.
"""

from __future__ import annotations

import numpy as np
import pytest

from thyra.core.mass_axis import MassAxisAccumulator
from thyra.errors import ConversionRefused


class TestTheFoldReturnsWhatNpUniqueReturned:
    """The property the whole refactor rests on."""

    @pytest.mark.parametrize(
        "batches",
        [
            pytest.param([[3.0, 1.0, 2.0]], id="one_spectrum"),
            pytest.param([[1.0, 2.0], [2.0, 3.0], [1.0, 3.0]], id="overlapping"),
            pytest.param([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], id="disjoint"),
            pytest.param([[1.0, 1.0, 1.0]], id="all_duplicates"),
            pytest.param([[], [1.0], []], id="empty_spectra_interleaved"),
            pytest.param([[1.0, np.nan], [np.nan, 2.0]], id="nan"),
            pytest.param([[-1.0, 0.0, 1.0]], id="negative_and_zero"),
        ],
    )
    def test_it_matches(self, batches):
        arrays = [np.asarray(b, dtype=np.float64) for b in batches]
        accumulator = MassAxisAccumulator()
        for a in arrays:
            accumulator.add(a)

        expected = np.unique(np.concatenate([a for a in arrays if a.size]))
        np.testing.assert_array_equal(accumulator.finish(), expected)

    def test_many_small_batches_agree_with_one_big_one(self):
        """The fold path, exercised: enough values to force real merges."""
        rng = np.random.default_rng(294)
        values = rng.integers(0, 5_000, size=40_000).astype(np.float64)

        accumulator = MassAxisAccumulator()
        for chunk in np.array_split(values, 500):
            accumulator.add(chunk)

        np.testing.assert_array_equal(accumulator.finish(), np.unique(values))


class TestTheNaNDivergenceIsGone:
    """``set`` dedupes NaN by identity and ``sorted()`` on NaN is
    meaningless, so the timsTOF builder returned an UNSORTED axis with the
    NaNs duplicated. Every other builder returned ``np.unique``'s answer.

    It never shipped -- the converter refuses a non-finite axis -- but the
    five builders disagreed about it, which is the thing being removed.
    """

    def test_the_old_set_shape_really_did_diverge(self):
        """The control, so the claim is not taken on trust."""
        u: set = set()
        for a in ([1.0, np.nan], [np.nan, 2.0]):
            u.update(np.asarray(a, dtype=np.float64))
        old = np.array(sorted(u))

        assert old.size == 4  # two NaNs survived
        assert not np.array_equal(old, np.unique(old))  # and it is unsorted

    def test_the_fold_collapses_them_like_np_unique(self):
        accumulator = MassAxisAccumulator()
        for a in ([1.0, np.nan], [np.nan, 2.0]):
            accumulator.add(np.asarray(a, dtype=np.float64))

        axis = accumulator.finish()
        assert int(np.count_nonzero(np.isnan(axis))) == 1
        np.testing.assert_array_equal(axis[:-1], [1.0, 2.0])


class TestAnEmptySourceIsStillRefused:
    """No caller passes a loop variable to ``finish()`` any more, so an
    empty source reaches the refusal instead of an ``UnboundLocalError``."""

    def test_no_spectra_at_all(self):
        with pytest.raises(ConversionRefused, match="No spectra found"):
            MassAxisAccumulator().finish()

    def test_every_spectrum_empty(self):
        accumulator = MassAxisAccumulator()
        for _ in range(3):
            accumulator.add(np.asarray([], dtype=np.float64))
        with pytest.raises(ConversionRefused, match="No spectra found"):
            accumulator.finish()
        assert accumulator.n_seen == 3

    def test_finish_takes_no_argument(self):
        """The signature is the fix. A caller cannot pass a loop variable
        that may never have been bound, because there is nowhere to put it."""
        import inspect

        params = list(inspect.signature(MassAxisAccumulator.finish).parameters)
        assert params == ["self"]


class TestTheCapStillWorks:
    def test_it_refuses_once_the_axis_outgrows_the_cap(self):
        accumulator = MassAxisAccumulator(total_spectra=100, max_length=3)
        with pytest.raises(ConversionRefused, match="exceeded 3 unique"):
            for i in range(10):
                accumulator.add(np.asarray([float(i), i + 0.5], dtype=np.float64))
            accumulator.finish()

    def test_the_message_counts_spectra_without_a_caller_tracking_them(self):
        accumulator = MassAxisAccumulator(total_spectra=100, max_length=3)
        with pytest.raises(ConversionRefused) as excinfo:
            for i in range(10):
                accumulator.add(np.asarray([float(i), i + 0.5], dtype=np.float64))
            accumulator.finish()
        assert "of 100 spectra" in str(excinfo.value)

    def test_an_unknown_total_claims_no_denominator(self):
        """solariX walks a cursor and would need a second COUNT(*) to say."""
        accumulator = MassAxisAccumulator(max_length=3)
        with pytest.raises(ConversionRefused) as excinfo:
            for i in range(10):
                accumulator.add(np.asarray([float(i), i + 0.5], dtype=np.float64))
            accumulator.finish()
        message = str(excinfo.value)
        assert "spectra" in message and " of " not in message.split("spectra")[0][-20:]

    def test_no_cap_means_no_refusal(self):
        accumulator = MassAxisAccumulator()
        for i in range(50):
            accumulator.add(np.asarray([float(i)], dtype=np.float64))
        assert accumulator.finish().size == 50


class TestEveryReaderUsesIt:
    """A completeness audit: a reader that builds a union axis must fold.

    Discovered from the source rather than listed, for the reason the
    drop-tally audit is: a hardcoded list can only confirm the builders
    that were converted. `np.concatenate`/`np.union1d` over a per-spectrum
    accumulation is the shape being removed.
    """

    READERS = [
        "thyra/readers/waters/waters_reader.py",
        "thyra/readers/bruker/solarix/solarix_reader.py",
        "thyra/readers/mzpeak/mzpeak_reader.py",
        "thyra/readers/bruker/timstof/timstof_reader.py",
        "thyra/readers/imzml/imzml_reader.py",
    ]

    @pytest.mark.parametrize("path", READERS, ids=lambda p: p.split("/")[-1])
    def test_it_imports_the_shared_builder(self, path):
        from pathlib import Path

        source = Path(path).read_text(encoding="utf-8")
        assert "MassAxisAccumulator" in source, f"{path} builds its own axis"

    @pytest.mark.parametrize("path", READERS, ids=lambda p: p.split("/")[-1])
    def test_it_no_longer_collects_then_uniques(self, path):
        """Checked against the AST, and only inside the axis builders.

        Two narrowings, both learned by getting it wrong. Every one of
        these files now *documents* the shape it used to have, so a
        substring search hits the prose explaining why it is gone. And
        ``np.unique(np.concatenate(...))`` is a perfectly ordinary thing
        to write about something that is not a mass axis -- imzML uses it
        to sample row indices -- so the whole-module search flagged that
        too. Only calls inside the functions that build an axis count.
        """
        import ast
        from pathlib import Path

        builders = {
            "get_common_mass_axis",
            "build_raw_mass_axis",
            "_extract_continuous_mass_axis",
        }

        def called(node):
            f = node.func
            if isinstance(f, ast.Attribute):
                return f.attr
            return f.id if isinstance(f, ast.Name) else ""

        tree = ast.parse(Path(path).read_text(encoding="utf-8"))
        offenders = []
        for fn in ast.walk(tree):
            if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if fn.name not in builders:
                continue
            for node in ast.walk(fn):
                if not isinstance(node, ast.Call):
                    continue
                name = called(node)
                if name == "union1d":
                    offenders.append(f"{fn.name}: np.union1d")
                elif name == "concatenate":
                    # Any concatenate inside a builder means it is holding
                    # the spectra to join them. Not just the nested
                    # ``np.unique(np.concatenate(...))`` spelling: Waters
                    # split those across two statements, so matching only
                    # the nested form let the shape through.
                    offenders.append(f"{fn.name}: np.concatenate")

        assert not offenders, f"{path} still builds an axis with {offenders}"
