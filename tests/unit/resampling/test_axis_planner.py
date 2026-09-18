# tests/unit/resampling/test_axis_planner.py
"""The axis decision, asked without a conversion.

Until #352 every one of these questions -- which law, how many bins, how
wide, and what the store's provenance would record -- could only be asked
by standing up a converter with a reader and a writable output path, and
the answers lived on sixteen private attributes assigned in three places.
The planner takes the same two inputs the conversion has (what the
detectors are told about the source, and what the caller asked for) and
answers them on their own.

The numbers here are the defaults ``docs/design-decisions.md`` records:
5 mDa at m/z 1000 on every axis but ``linear_tof``, which takes 17 mDa at
m/z 300 to sit close to the axis SCiLS Lab lays for FlexImaging data; a
``tof`` axis sized in bins per measured peak width instead (D20); and a
linearisation that must come back usable, or every placement onto the
axis silently falls back to a binary search (D21).
"""

from __future__ import annotations

import numpy as np
import pytest

from thyra.errors import ConversionRefused
from thyra.resampling import common_axis
from thyra.resampling.axis_planner import AxisPlanner
from thyra.resampling.mass_axis.tof_generator import (
    DEFAULT_BINS_PER_FWHM,
    PHI_NANOTOF_LAW,
)
from thyra.resampling.types import AxisType, ResamplingConfig, ResamplingMethod

#: Matched by ``RapiflexDetector`` on the format flag alone: a constant
#: axis, and the cheapest source that reaches ``TIC_PRESERVING``.
_RAPIFLEX = {"format_specific": {"format": "Rapiflex"}}

#: Matched by ``PhiToFSIMSDetector``, which declares the two-term width
#: law fitted on 311 peaks across twelve nanoTOF acquisitions (D20).
_PHI = {"format_specific": {"format": "PHI SmartSoft-TOF raw"}}

_RANGE = (100.0, 1100.0)


def _planner(metadata=None, **config) -> AxisPlanner:
    """A planner over a 100-1100 m/z source, and the caller's request."""
    told = {"essential_metadata": {"mass_range": _RANGE}}
    told.update(metadata or {})
    return AxisPlanner.from_metadata(told, ResamplingConfig(**config))


class TestTheDecisionNeedsNoConverter:
    """A metadata dict and a config are the whole input."""

    def test_the_detector_decides_the_law_and_the_method(self):
        planner = _planner(_RAPIFLEX)

        plan = planner.resolve()

        assert plan.axis_type is AxisType.CONSTANT
        assert planner.method is ResamplingMethod.TIC_PRESERVING
        assert planner.method_was_auto
        assert (plan.min_mz, plan.max_mz) == _RANGE

    def test_a_source_nothing_matches_falls_to_the_default_detector(self):
        planner = _planner()

        assert planner.resolve().axis_type is AxisType.CONSTANT
        assert planner.method is ResamplingMethod.NEAREST_NEIGHBOR

    def test_five_milli_dalton_bins_at_a_thousand_are_the_default(self):
        """1,000 m/z of range in 5 mDa bins is 200,000 of them."""
        plan = _planner(_RAPIFLEX).resolve()

        assert plan.target_bins == 200_000

    def test_linear_tof_asks_for_seventeen_milli_dalton_bins_at_three_hundred(self):
        """The SCiLS Lab convention, and the one axis type with its own."""
        planner = _planner(axis_type=AxisType.LINEAR_TOF)

        assert planner.reference_params(AxisType.LINEAR_TOF) == (0.017, 300.0)

    def test_the_callers_width_beats_every_default(self):
        planner = _planner(_RAPIFLEX, mass_width_da=0.05, reference_mz=500.0)

        assert planner.reference_params(AxisType.CONSTANT) == (0.05, 500.0)
        assert planner.resolve().target_bins == 20_000

    def test_a_declared_width_law_reaches_the_plan(self):
        """D20: a ``tof`` axis is sized in bins per measured peak width."""
        planner = _planner(_PHI)
        plan = planner.resolve()

        assert plan.axis_type is AxisType.TOF
        width, reference_mz = planner.reference_params(AxisType.TOF)
        assert reference_mz == 1000.0
        assert width > 0.0
        assert plan.target_bins > 0


class TestWhatTheStrategyIsBuiltFrom:
    """``settle`` returns the triple ``build_strategy`` takes, and the plan."""

    def _settle(self, **config):
        planner = _planner(_RAPIFLEX, **config)
        return planner, planner.settle(lambda n_bins: None)

    def test_the_axis_is_the_one_the_plan_described(self):
        planner, settled = self._settle(target_bins=4_000, mass_width_da=None)

        assert len(settled.axis) == planner.resolve().target_bins
        assert settled.axis[0] == pytest.approx(_RANGE[0])
        assert settled.axis[-1] == pytest.approx(_RANGE[1])
        assert np.all(np.diff(settled.axis) > 0)

    def test_the_declared_range_travels_beside_the_axis(self):
        """Not the axis's own span: a physics axis reports bin centres."""
        _, settled = self._settle(target_bins=4_000)

        assert settled.axis_range == _RANGE

    def test_the_linearisation_comes_back_usable(self):
        """D21. Without it every placement onto the axis is a search."""
        _, settled = self._settle(target_bins=4_000)

        assert settled.linearisation is not None

    def test_the_provenance_records_what_was_actually_done(self):
        """With "auto" settings the request says nothing about any of it."""
        _, settled = self._settle(target_bins=4_000)

        assert settled.provenance == {
            "method": ResamplingMethod.TIC_PRESERVING,
            "axis_type": AxisType.CONSTANT,
            "target_bins": 4_000,
            "min_mz": 100.0,
            "max_mz": 1100.0,
            "mass_width_da": 0.005,
            "reference_mz": 1000.0,
        }

    def test_a_tof_axis_records_its_law_too(self):
        planner = _planner(_PHI, target_bins=4_000)
        settled = planner.settle(lambda n_bins: None)

        assert settled.provenance["tof_a"] == pytest.approx(PHI_NANOTOF_LAW[0])
        assert settled.provenance["tof_b"] == pytest.approx(PHI_NANOTOF_LAW[1])
        assert settled.provenance["bins_per_fwhm"] == DEFAULT_BINS_PER_FWHM


class TestTheMemoryBudgetIsAskedFirst:
    """Issue #251: a bin count that cannot fit is refused unbuilt."""

    def test_nothing_is_laid_before_the_budget_has_answered(self, monkeypatch):
        """300M bins were refused in 2.3 s, at 7.44 GB of peak RSS."""

        def _no_axis(*args, **kwargs):
            raise AssertionError("the axis was built before the budget answered")

        monkeypatch.setattr(
            common_axis.CommonAxisBuilder, "build_uniform_axis", _no_axis
        )
        planner = _planner(_RAPIFLEX, target_bins=4_000)

        def _refuse(n_bins):
            raise ConversionRefused(f"{n_bins:,} bins is too many")

        with pytest.raises(ConversionRefused, match="4,000 bins"):
            planner.settle(_refuse)

    def test_the_budget_is_asked_the_count_the_plan_resolved(self):
        asked = []
        planner = _planner(_RAPIFLEX, target_bins=4_000)

        planner.settle(asked.append)

        assert asked == [4_000]


class _Reader:
    """A source with an axis of its own and nothing else."""

    has_shared_mass_axis = True

    def __init__(self, axis=(100.0, 200.0, 300.0)):
        self._axis = np.asarray(axis, dtype=np.float64)

    def get_common_mass_axis(self):
        return self._axis


class TestTheSourcesOwnAxis:
    """``--no-resample``: deciding not to lay one is still a decision."""

    def test_the_axis_is_the_readers_own(self):
        settled = AxisPlanner(_Reader()).settle(lambda n_bins: None)

        assert list(settled.axis) == [100.0, 200.0, 300.0]

    def test_nothing_about_a_raw_axis_is_a_resampling(self):
        """No generator laid it, so it has no coordinate and no operator."""
        settled = AxisPlanner(_Reader()).settle(lambda n_bins: None)

        assert settled.linearisation is None
        assert settled.method is None
        assert settled.provenance is None
        assert settled.axis_range is None

    def test_the_budget_is_asked_once_it_is_known(self):
        asked = []

        AxisPlanner(_Reader()).settle(asked.append)

        assert asked == [3]

    def test_an_empty_axis_is_refused(self):
        with pytest.raises(ConversionRefused, match="Common mass axis is empty"):
            AxisPlanner(_Reader(axis=())).settle(lambda n_bins: None)
