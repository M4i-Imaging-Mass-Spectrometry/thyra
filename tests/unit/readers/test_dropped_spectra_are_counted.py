"""A dropped spectrum is counted, not narrated once per spectrum (#308).

Every reader frame loop announced each dropped spectrum with its own
``logger.warning`` and counted none of them. On a 17,423-pixel
acquisition whose source is systematically unreadable -- a truncated
``.tdf_bin``, a region whose ``MaldiFrameInfo`` rows are missing -- that
is 17,423 warning lines, and twice over, because the converter reads the
source in two passes.

The one number a user needs, *how much of my data did not make it*, was
in none of them. The run then either succeeded with a short store or hit
the empty-conversion refusal (#242), neither of which names the drops as
the cause.

So: name the first one with its cause, log the rest at DEBUG, and report
the total once when the pass ends -- the shape the sibling tables
already use (``msms_table``'s ``n_dropped``).
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import pytest

from thyra.core.drop_tally import DropTally
from thyra.errors import ConversionRefused

_TALLY_LOGGER = "tests.drop_tally"


class TestTheTally:
    def _tally(self):
        return DropTally(logging.getLogger(_TALLY_LOGGER), "frames with no coordinates")

    def test_nothing_dropped_says_nothing(self, caplog):
        tally = self._tally()
        with caplog.at_level(logging.DEBUG, logger=_TALLY_LOGGER):
            tally.summarise(10)
        assert caplog.records == []

    def test_the_first_drop_warns_with_its_cause(self, caplog):
        tally = self._tally()
        with caplog.at_level(logging.DEBUG, logger=_TALLY_LOGGER):
            tally.drop("frame 7", OSError("truncated"))
        warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1
        assert "frame 7" in warnings[0].getMessage()
        assert "truncated" in warnings[0].getMessage()

    def test_the_rest_go_to_debug(self, caplog):
        tally = self._tally()
        with caplog.at_level(logging.DEBUG, logger=_TALLY_LOGGER):
            for frame_id in range(1, 6):
                tally.drop(f"frame {frame_id}")
        warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
        debugs = [r for r in caplog.records if r.levelno == logging.DEBUG]
        assert len(warnings) == 1
        assert len(debugs) == 4

    def test_the_total_is_reported_once_with_the_fraction(self, caplog):
        tally = self._tally()
        with caplog.at_level(logging.DEBUG, logger=_TALLY_LOGGER):
            for frame_id in range(1, 4):
                tally.drop(f"frame {frame_id}")
            tally.summarise(5)
        summaries = [
            r.getMessage() for r in caplog.records if "Dropped 3 of 5" in r.getMessage()
        ]
        assert len(summaries) == 1
        assert "frames with no coordinates" in summaries[0]

    def test_an_unknown_total_still_reports_the_count(self, caplog):
        tally = self._tally()
        with caplog.at_level(logging.DEBUG, logger=_TALLY_LOGGER):
            tally.drop("frame 1")
            tally.summarise(None)
        assert any("Dropped 1 frames" in r.getMessage() for r in caplog.records)


class TestTheFrameLoopUsesIt:
    """``_iter_frames`` is the shared loop behind four public iterators."""

    def _reader(self, coords_for):
        """A ``BrukerReader`` reduced to what ``_iter_frames`` reads."""
        from types import SimpleNamespace

        from thyra.readers.bruker.timstof.timstof_reader import BrukerReader

        stub = SimpleNamespace(
            _region_frame_ids=None,
            progress_callback=None,
            _get_maldi_frame_ids=lambda: [1, 2, 3, 4, 5],
            _get_coordinate_offsets=lambda: None,
            _get_frame_coordinates_cached=lambda frame_id, offsets: coords_for(
                frame_id
            ),
        )
        return lambda: BrukerReader._iter_frames(stub)

    def test_three_missing_of_five_warn_once_and_summarise_once(self, caplog):
        iter_frames = self._reader(
            lambda frame_id: None if frame_id in (2, 3, 4) else (frame_id, 0, 0)
        )
        module = "thyra.readers.bruker.timstof.timstof_reader"
        with caplog.at_level(logging.DEBUG, logger=module):
            kept = list(iter_frames())

        assert [frame_id for frame_id, _ in kept] == [1, 5]
        messages = [r.getMessage() for r in caplog.records if r.name == module]
        warnings = [
            r.getMessage()
            for r in caplog.records
            if r.name == module and r.levelno == logging.WARNING
        ]
        # One naming the first dropped frame, one carrying the total.
        assert len(warnings) == 2
        assert "frame 2" in warnings[0]
        assert "Dropped 3 of 5" in warnings[1]
        # Frames 3 and 4 are still findable, at DEBUG.
        assert any("frame 3" in m for m in messages)
        assert any("frame 4" in m for m in messages)

    def test_a_clean_pass_says_nothing(self, caplog):
        iter_frames = self._reader(lambda frame_id: (frame_id, 0, 0))
        module = "thyra.readers.bruker.timstof.timstof_reader"
        with caplog.at_level(logging.WARNING, logger=module):
            assert len(list(iter_frames())) == 5
        assert [r for r in caplog.records if r.name == module] == []


class TestARefusalIsNotDemotedToADrop:
    """``ConversionRefused`` subclasses ``ValueError``, so the broad
    ``except Exception`` in these loops caught the deliberate
    mobility/m-z length-mismatch refusal and logged it as one bad
    spectrum among thousands.

    Driven, not read: asserting that the source text contains
    ``except ConversionRefused:`` passes on a clause that is unreachable,
    on one that re-raises the wrong thing, and on a comment that happens
    to contain the string.
    """

    def test_the_imzml_mobility_loop_lets_it_out(self, monkeypatch, tmp_path):
        """The mismatch refusal is raised inside the loop's own ``try``."""
        from thyra.readers.imzml import imzml_reader

        reader = imzml_reader.ImzMLReader.__new__(imzml_reader.ImzMLReader)
        refusal = ConversionRefused("3 mobility values for 5 m/z values")

        def explode(self, *args, **kwargs):
            raise refusal

        monkeypatch.setattr(
            imzml_reader.ImzMLReader, "_ensure_parser_initialized", lambda self: None
        )
        monkeypatch.setattr(
            imzml_reader.ImzMLReader, "_mobility", object(), raising=False
        )
        monkeypatch.setattr(
            imzml_reader.ImzMLReader,
            "_get_spectrum_coordinates",
            explode,
            raising=False,
        )

        parser = SimpleNamespace(coordinates=[(1, 1, 1), (2, 1, 1)])
        monkeypatch.setattr(imzml_reader.ImzMLReader, "parser", parser, raising=False)
        monkeypatch.setattr(
            imzml_reader.ImzMLReader, "is_continuous", False, raising=False
        )
        monkeypatch.setattr(
            imzml_reader.ImzMLReader, "_quiet_mode", True, raising=False
        )

        with pytest.raises(ConversionRefused, match="mobility values"):
            list(reader.iter_mobility_spectra())

    def test_the_imzml_spectrum_loop_lets_it_out(self, monkeypatch):
        """The same clause on the summed route."""
        from thyra.readers.imzml import imzml_reader

        reader = imzml_reader.ImzMLReader.__new__(imzml_reader.ImzMLReader)
        monkeypatch.setattr(
            imzml_reader.ImzMLReader,
            "_get_spectrum_coordinates",
            lambda self, parser, idx: (_ for _ in ()).throw(
                ConversionRefused("refused")
            ),
            raising=False,
        )
        tally = DropTally(logging.getLogger("tests.refusal"), "spectra")
        pbar = SimpleNamespace(update=lambda n: None)

        with pytest.raises(ConversionRefused, match="refused"):
            reader._process_single_spectrum(object(), 0, pbar, tally)
        # And it was not counted as a bad spectrum on the way out.
        assert tally.n_dropped == 0

    def test_an_ordinary_failure_is_still_counted(self, monkeypatch):
        """The broad clause still does its job for everything else."""
        from thyra.readers.imzml import imzml_reader

        reader = imzml_reader.ImzMLReader.__new__(imzml_reader.ImzMLReader)
        monkeypatch.setattr(
            imzml_reader.ImzMLReader,
            "_get_spectrum_coordinates",
            lambda self, parser, idx: (_ for _ in ()).throw(OSError("truncated")),
            raising=False,
        )
        tally = DropTally(logging.getLogger("tests.refusal"), "spectra")
        pbar = SimpleNamespace(update=lambda n: None)

        assert reader._process_single_spectrum(object(), 0, pbar, tally) is None
        assert tally.n_dropped == 1

    def test_conversion_refused_really_is_a_value_error(self):
        """The premise. If this ever stops holding, the guards above can
        be reconsidered -- but never by reading a test, only the raise."""
        assert issubclass(ConversionRefused, ValueError)


def _reader_functions():
    """Every method of the two reader classes, as ``(label, node, source)``.

    Walked out of the class body rather than filtered by ``hasattr``:
    ``hasattr`` also matches a module-level function or another class's
    method that happens to share a name, and mis-labels it as the
    reader's.
    """
    import ast
    import importlib
    import inspect

    out = []
    for module_name, owner_name in (
        ("thyra.readers.bruker.timstof.timstof_reader", "BrukerReader"),
        ("thyra.readers.imzml.imzml_reader", "ImzMLReader"),
    ):
        module = importlib.import_module(module_name)
        tree = ast.parse(inspect.getsource(module))
        owner = next(
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.ClassDef) and node.name == owner_name
        )
        for node in owner.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                out.append((f"{owner_name}.{node.name}", node, ast.unparse(node)))
    return out


def _swallowing_generators():
    """Generators that catch an exception and carry on to the next spectrum.

    This is the exact shape the issue is about: one warning per dropped
    spectrum, a total nowhere. Discovered rather than listed, because a
    hardcoded list can only confirm the loops that were fixed -- two TDF
    mobility loops were missed on the first pass for precisely that
    reason, and a listed parametrisation could never have noticed.

    A handler that re-raises is not swallowing, and a function that does
    not yield is not handing the caller the rest of the run:
    ``_isolation_windows`` continues past a ``sqlite3.OperationalError``
    but is probing which schema variant the file has, not discarding
    data. Both exclusions were measured against the tree, not assumed.

    **It cannot see a loop whose handler lives in a helper it calls** --
    which is the shape of imzML's ``_iter_spectra_single`` and
    ``_iter_spectra_batch``, where the drop happens inside
    ``_process_single_spectrum``. Those are caught by the second rule
    below instead, because they open a tally. A *new* delegating loop
    that drops with no tally at all would be invisible to both, so this
    is an aid to review rather than a proof of completeness.
    """
    import ast

    found = []
    for label, node, source in _reader_functions():
        if not any(isinstance(n, (ast.Yield, ast.YieldFrom)) for n in ast.walk(node)):
            continue
        for handler in (n for n in ast.walk(node) if isinstance(n, ast.ExceptHandler)):
            body = list(ast.walk(handler))
            if any(isinstance(s, ast.Raise) for s in body):
                continue  # re-raises: a decision, not a drop
            found.append((label, source))
            break
    return found


@pytest.mark.parametrize(
    "label, source",
    _swallowing_generators(),
    # Label only. Letting pytest build an id from the source puts whole
    # function bodies into the test id, and into every CI log line.
    ids=[label for label, _ in _swallowing_generators()],
)
def test_a_loop_that_swallows_counts_what_it_swallowed(label, source):
    assert "tally.drop(" in source, f"{label} drops spectra without counting them"


@pytest.mark.parametrize(
    "label, source",
    [(lbl, src) for lbl, _, src in _reader_functions() if "DropTally(" in src],
    ids=[lbl for lbl, _, src in _reader_functions() if "DropTally(" in src],
)
def test_whoever_opens_a_tally_closes_it(label, source):
    """Counting without reporting is the half of the defect that is easy
    to reintroduce: the loop looks fixed and still says nothing."""
    assert "tally.summarise(" in source, f"{label} counts drops but never reports them"
