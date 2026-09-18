# tests/unit/converters/test_resample_override_contradiction.py
"""An explicit ``--resample-method`` that overrules the detector says so.

The detector has a verdict for every source, and until #246 only the
``auto`` path ever asked for it. An explicit method was applied with
nothing checked and nothing said, which is how ``tic_preserving`` on a
Bruker TDF -- for which the detector chooses nearest-neighbour --
interpolated across the gaps of a sparse centroid list and filled the
axis. Measured on a 713-frame PASEF acquisition:

======================  ==================  ===============
                          nearest_neighbor   tic_preserving
======================  ==================  ===============
stored non-zeros                   302,106      423,386,757
table                              7.6 MB           583 MB
peak RSS                            0.5 GB           5.9 GB
======================  ==================  ===============

Per-pixel TIC is identical either way, so the TIC identity cannot see it.
What breaks is the siblings, quietly: the heatmap marginal against the
stored mean spectrum came to rel 68, with nothing recorded anywhere.

Same bug class as #168 on PHI ToF-SIMS -- which was fixed *by* adding a
detector, and that is exactly why a detector is not enough on its own. A
detector only ever steers ``auto``.

Nothing stored changes here (design decision D15): the warning names
``--resample-gap-tolerance``, which already exists and, measured on the
same acquisition, brings those 423,386,757 entries down to 4,801,946.
"""

from __future__ import annotations

import logging

import pytest

from thyra.resampling.axis_planner import AxisPlanner, normalize_resampling_config
from thyra.resampling.types import ResamplingMethod

_LOGGER = "thyra.resampling.axis_planner"


#: A sparse centroid source the detector reads as nearest-neighbour. The
#: PHI format flag is the cheapest unambiguous way to get that verdict
#: without a vendor file; ``PhiToFSIMSDetector`` matches on it alone.
_SPARSE_SOURCE = {
    "source_path": "acquisition.raw",
    "essential_metadata": {
        "spectrum_type": "profile spectrum",
        "dimensions": (4, 4, 1),
        "mass_range": (1.0, 100.0),
        "source_path": "acquisition.raw",
        "total_peaks": 64,
        "n_spectra": 16,
    },
    "format_specific": {"format": "PHI SmartSoft-TOF raw"},
}


def _run(config, metadata=None):
    """The planner's own setup: the method is chosen when it is built."""
    return AxisPlanner.from_metadata(
        metadata if metadata is not None else _SPARSE_SOURCE, config
    )


def _config(**overrides):
    settings = {"method": "nearest_neighbor", "axis_type": "auto"}
    settings.update(overrides)
    return normalize_resampling_config(settings)


def _contradictions(records):
    return [r.getMessage() for r in records if "detector chose" in r.getMessage()]


class TestTheDetectorIsConsultedOnAnOverride:
    def test_it_warns_when_the_override_contradicts(self, thyra_logs):
        with thyra_logs(_LOGGER, logging.WARNING) as records:
            planner = _run(_config(method="tic_preserving"))

        messages = _contradictions(records)
        assert len(messages) == 1, [r.getMessage() for r in records]
        assert "TIC_PRESERVING" in messages[0]
        assert "NEAREST_NEIGHBOR" in messages[0]
        # The override still applies: this reports, it does not overrule.
        assert planner.method is ResamplingMethod.TIC_PRESERVING

    def test_it_names_the_flag_that_fixes_it(self, thyra_logs):
        """A warning with no way out is just noise."""
        with thyra_logs(_LOGGER, logging.WARNING) as records:
            _run(_config(method="tic_preserving"))

        assert "--resample-gap-tolerance" in _contradictions(records)[0]

    def test_it_says_so_when_the_tolerance_is_already_set(self, thyra_logs):
        """Do not tell someone to pass a flag they have already passed."""
        with thyra_logs(_LOGGER, logging.WARNING) as records:
            _run(_config(method="tic_preserving", gap_tolerance_da=0.01))

        message = _contradictions(records)[0]
        assert "0.01" in message
        assert "discarded rather than interpolated" in message

    def test_it_is_silent_when_the_override_agrees(self, thyra_logs):
        """No noise on an explicit flag that matches the detector."""
        with thyra_logs(_LOGGER, logging.WARNING) as records:
            planner = _run(_config(method="nearest_neighbor"))

        assert _contradictions(records) == []
        assert planner.method is ResamplingMethod.NEAREST_NEIGHBOR

    def test_it_is_silent_on_auto(self, thyra_logs):
        """``auto`` *is* the detector; it cannot contradict itself."""
        with thyra_logs(_LOGGER, logging.WARNING) as records:
            planner = _run(_config(method="auto"))

        assert _contradictions(records) == []
        assert planner.method is ResamplingMethod.NEAREST_NEIGHBOR


class TestTheCheckNeverBreaksAConversion:
    def test_a_source_the_detector_cannot_classify_still_converts(self, thyra_logs):
        """Detection is advisory here, so a failure is a debug line.

        The method the caller asked for has to be applied either way: a
        conversion that refused because an advisory check raised would be
        a worse bug than the one being fixed.
        """

        class _Exploding(AxisPlanner):
            def detection_metadata(self):
                raise RuntimeError("no metadata")

        with thyra_logs(_LOGGER, logging.WARNING) as records:
            planner = _Exploding(None, _config(method="tic_preserving"))

        assert planner.method is ResamplingMethod.TIC_PRESERVING
        assert _contradictions(records) == []

    def test_the_rest_of_the_setup_still_runs(self):
        """The warning sits in the middle of the planner's own setup.

        The tolerance it names reaches the operator off the config the
        planner was built from, so the method and that config are what
        the warning must leave intact.
        """
        planner = _run(_config(method="tic_preserving", gap_tolerance_da=0.25))

        assert planner.method is ResamplingMethod.TIC_PRESERVING
        assert planner.config.gap_tolerance_da == pytest.approx(0.25)
