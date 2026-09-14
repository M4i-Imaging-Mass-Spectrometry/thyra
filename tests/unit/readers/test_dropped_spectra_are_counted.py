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
    ``except Exception`` in the imzML loops caught the deliberate
    mobility/m-z length-mismatch refusal and logged it as one bad
    spectrum among thousands."""

    def test_the_imzml_spectrum_loop_re_raises_it(self):
        import inspect

        from thyra.readers.imzml import imzml_reader

        source = inspect.getsource(imzml_reader.ImzMLReader._process_single_spectrum)
        assert "except ConversionRefused:" in source
        assert source.index("except ConversionRefused:") < source.index(
            "except Exception"
        )

    def test_the_imzml_mobility_loop_re_raises_it(self):
        import inspect

        from thyra.readers.imzml import imzml_reader

        source = inspect.getsource(imzml_reader.ImzMLReader.iter_mobility_spectra)
        assert "except ConversionRefused:" in source
        assert source.index("except ConversionRefused:") < source.index(
            "except Exception"
        )

    def test_conversion_refused_really_is_a_value_error(self):
        """The premise. If this ever stops holding, the guards above can
        be reconsidered -- but never by reading a test, only the raise."""
        assert issubclass(ConversionRefused, ValueError)


@pytest.mark.parametrize(
    "module_name, attribute",
    [
        ("thyra.readers.bruker.timstof.timstof_reader", "_iter_frames"),
        ("thyra.readers.bruker.timstof.timstof_reader", "_iter_spectra_raw"),
        ("thyra.readers.bruker.timstof.timstof_reader", "iter_frame_scans"),
        ("thyra.readers.imzml.imzml_reader", "_iter_spectra_single"),
        ("thyra.readers.imzml.imzml_reader", "_iter_spectra_batch"),
        ("thyra.readers.imzml.imzml_reader", "iter_mobility_spectra"),
    ],
)
def test_every_drop_loop_summarises(module_name, attribute):
    """Each loop that can drop a spectrum reports its total once."""
    import importlib
    import inspect

    module = importlib.import_module(module_name)
    owner = module.BrukerReader if "timstof" in module_name else module.ImzMLReader
    source = inspect.getsource(getattr(owner, attribute))
    assert "tally.summarise(" in source, f"{attribute} counts drops but never says so"
