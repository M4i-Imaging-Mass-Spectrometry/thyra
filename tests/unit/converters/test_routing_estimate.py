# tests/unit/converters/test_routing_estimate.py
"""The size estimate must be honest, and routing must not depend on it.

``_estimate_output_size_gb`` scores a dataset as
``n_pixels * n_mz_bins * 4``. The resampling branch already resolved the
real bin count (issue #87), but the raw-axis branch still guessed it from
``(max_mass - min_mass) / 0.01`` -- a 10 mDa spacing the data need not have.

That guess is wrong in both directions. A continuous file carrying 4,000
points over 250-1200 m/z was scored as though it had 95,000, inflating it
24x; a processed file whose spectra share no m/z values was scored far too
low, which routed the very largest datasets to the method that holds the
most in memory.

**The routing half of that is now history.** ``_should_use_pcs`` no longer
consults the estimate at all: PCS is faster and lighter at every size
measured, so ``"auto"`` picks it unconditionally and
``PCS_SIZE_THRESHOLD_GB`` is gone. The estimate survives as a log line, and
these tests survive with it -- an inaccurate number in a support log is a
smaller problem than an inaccurate route, but it is still a problem, and the
24x inflation is the kind of thing that gets re-derived if nobody wrote down
that the fallback is unreliable.

``convert()`` runs ``_initialize_conversion()`` before reaching the
estimate, so the axis is already built and there is nothing to guess. These
tests pin that the built axis is preferred, that the old heuristics still
apply when it is not available -- which is the case when the estimator is
called directly, as the older tests in ``test_streaming_converter.py`` do --
and that no combination of sizes changes the route.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pytest

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.converters.spatialdata.streaming_converter import (
    SPATIALDATA_AVAILABLE,
    StreamingSpatialDataConverter,
)


class _Meta:
    def __init__(self, dimensions, mass_range):
        self.dimensions = dimensions
        self.mass_range = mass_range
        self.n_spectra = dimensions[0] * dimensions[1] * dimensions[2]
        self.pixel_size = (20.0, 20.0)
        self.estimated_memory_gb = 1.0
        self.coordinate_bounds = (0, dimensions[0], 0, dimensions[1])
        self.is_3d = dimensions[2] > 1
        self.has_mass_axis = True
        self.source_path = "mock"


class _Reader:
    def __init__(self, dimensions, mass_range):
        self._meta = _Meta(dimensions, mass_range)

    def get_essential_metadata(self):
        return self._meta


def _converter(
    dimensions: Tuple[int, int, int],
    mass_range: Tuple[float, float] = (250.0, 1200.0),
    axis: Optional[np.ndarray] = None,
    resampling_config=None,
):
    """Build a converter, optionally with the mass axis already resolved."""
    with tempfile.TemporaryDirectory() as tmpdir:
        conv = StreamingSpatialDataConverter(
            reader=_Reader(dimensions, mass_range),
            output_path=Path(tmpdir) / "out.zarr",
            use_csc="auto",
            resampling_config=resampling_config,
        )
    conv._common_mass_axis = axis
    return conv


class TestPrefersTheBuiltAxis:
    """When the axis exists, its length is what counts."""

    def test_uses_real_axis_length(self):
        n_pixels = 10_000
        axis = np.linspace(250.0, 1200.0, 4_000)
        conv = _converter((100, 100, 1), axis=axis)

        expected = n_pixels * 4_000 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)

    def test_real_axis_beats_the_raw_heuristic(self):
        """The heuristic would say 95,000 bins; the axis says 4,000."""
        axis = np.linspace(250.0, 1200.0, 4_000)
        conv = _converter((100, 100, 1), axis=axis)
        with_axis = conv._estimate_output_size_gb()

        conv._common_mass_axis = None
        without_axis = conv._estimate_output_size_gb()

        # (1200 - 250) / 0.01 = 95,000, i.e. 23.75x the real count.
        assert without_axis == pytest.approx(with_axis * 95_000 / 4_000, rel=1e-6)

    def test_real_axis_beats_the_resampling_plan_too(self):
        """A built axis is authoritative even when resampling is configured."""
        axis = np.linspace(250.0, 1200.0, 1_000)
        conv = _converter(
            (100, 100, 1),
            axis=axis,
            resampling_config={"method": "nearest_neighbor", "target_bins": 500_000},
        )
        expected = 10_000 * 1_000 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)


class TestFallbacksStillApply:
    """Without a built axis, the previous behaviour is unchanged."""

    def test_raw_heuristic_when_no_axis_and_no_resampling(self):
        conv = _converter((100, 100, 1), mass_range=(250.0, 1200.0), axis=None)
        expected = 10_000 * 95_000 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)

    def test_resampling_plan_when_no_axis(self):
        conv = _converter(
            (50, 50, 1),
            mass_range=(100.0, 1000.0),
            axis=None,
            resampling_config={"method": "nearest_neighbor", "target_bins": 5_000},
        )
        expected = 2_500 * 5_000 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)


class TestRoutingIsSizeIndependent:
    """``"auto"`` means PCS, whatever the estimate says.

    Each case below picked a *different* route under the old threshold, so
    together they pin that the size no longer reaches the decision.
    """

    def test_tiny_dataset_still_routes_to_pcs(self):
        """A 400-byte dataset. Went to COO under the threshold."""
        conv = _converter((10, 10, 1), axis=np.linspace(250.0, 1200.0, 10))

        assert conv._estimate_output_size_gb() < 1.0
        assert conv._should_use_pcs() is True

    def test_narrow_real_axis_routes_to_pcs(self):
        """Comfortably under the old 30 GB gate, so this went to COO."""
        conv = _converter((1000, 600, 1), axis=np.linspace(250.0, 1200.0, 4_000))

        assert conv._estimate_output_size_gb() < 30.0
        assert conv._should_use_pcs() is True

    def test_wide_real_axis_routes_to_pcs(self):
        """Over the old gate, so this one already went to PCS. Still does."""
        conv = _converter((250, 200, 1), axis=np.linspace(250.0, 1200.0, 200_000))

        assert conv._estimate_output_size_gb() > 30.0
        assert conv._should_use_pcs() is True

    def test_route_does_not_move_when_the_estimate_does(self):
        """The strongest form: swing the estimate 24x, route must not budge.

        Dropping the built axis makes the estimator fall back to the 10 mDa
        heuristic, which inflates this dataset from ~9 GB to ~223 GB -- a
        swing that straddled the old threshold and flipped the route. If a
        size gate is ever reintroduced, this is what fails.
        """
        conv = _converter((1000, 600, 1), axis=np.linspace(250.0, 1200.0, 4_000))

        with_axis = conv._estimate_output_size_gb()
        assert conv._should_use_pcs() is True

        conv._common_mass_axis = None
        without_axis = conv._estimate_output_size_gb()

        assert without_axis > with_axis * 20  # 95,000 / 4,000 = 23.75x
        assert conv._should_use_pcs() is True

    def test_explicit_use_csc_still_wins(self):
        """COO stays reachable: it is an escape hatch, not dead code."""
        conv = _converter((10, 10, 1), axis=np.linspace(250.0, 1200.0, 10))

        conv._use_csc = True
        assert conv._should_use_pcs() is True
        conv._use_csc = False
        assert conv._should_use_pcs() is False

    def test_the_threshold_constant_is_gone(self):
        """A leftover constant would read as a live gate. It is not one.

        Left behind, the next person tunes it and nothing happens.
        """
        assert not hasattr(StreamingSpatialDataConverter, "PCS_SIZE_THRESHOLD_GB")


class TestDegenerate:
    """Empty and tiny axes must not raise."""

    def test_empty_axis_is_zero_sized(self):
        conv = _converter((10, 10, 1), axis=np.array([]))
        assert conv._estimate_output_size_gb() == 0.0
        # Zero is a size like any other now; it does not divert the route.
        assert conv._should_use_pcs() is True

    def test_single_bin_axis(self):
        conv = _converter((10, 10, 1), axis=np.array([500.0]))
        expected = 100 * 1 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)


@pytest.mark.skipif(
    not SPATIALDATA_AVAILABLE,
    reason="SpatialData dependencies not available",
)
class TestTheRouteIsActuallyTaken:
    """End-to-end: which branch ``convert()`` reaches, not what a predicate says.

    ``_should_use_pcs`` returning True is not the same as the conversion
    taking the PCS branch -- ``convert()`` has to consult it and act on the
    answer. Asserting the predicate alone would keep passing if that wiring
    were bypassed, so this drives a real conversion and records which of the
    two entry points ran.
    """

    @staticmethod
    def _route_taken(**kwargs) -> str:
        reader = MockMSIReader(
            MockMSIConfig(n_x=6, n_y=4, n_mz_bins=400, peaks_per_spectrum=(20, 40))
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            converter = StreamingSpatialDataConverter(
                reader=reader,
                output_path=Path(tmpdir) / "out.zarr",
                dataset_id="route",
                pixel_size_um=10.0,
                **kwargs,
            )

            calls: list[str] = []

            def _spy(name):
                original = getattr(converter, name)

                def wrapper(*args, **kwargs_):
                    calls.append(name)
                    return original(*args, **kwargs_)

                return wrapper

            for entry_point in ("_convert_to_csc_no_cache", "_stream_build_coo"):
                setattr(converter, entry_point, _spy(entry_point))

            assert converter.convert() is True, "conversion fixture failed"

        if "_convert_to_csc_no_cache" in calls:
            return "pcs"
        if "_stream_build_coo" in calls:
            return "coo"
        raise AssertionError(f"neither route ran: {calls}")

    def test_default_conversion_takes_pcs(self):
        """No ``use_csc`` at all -- the case every in-repo caller hits.

        Nothing in ``thyra/`` passes ``use_csc``, so this is what the CLI and
        the public API did, and before this change it was COO.
        """
        assert self._route_taken() == "pcs"

    def test_explicit_true_takes_pcs(self):
        assert self._route_taken(use_csc=True) == "pcs"

    def test_explicit_false_still_takes_coo(self):
        """COO must stay reachable, or the tests covering it stop meaning anything."""
        assert self._route_taken(use_csc=False) == "coo"
