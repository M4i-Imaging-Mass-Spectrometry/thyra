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

import logging

import numpy as np
import pytest

from thyra.errors import ConversionRefused
from thyra.resampling import common_axis
from thyra.resampling.axis_planner import AxisPlanner, normalize_resampling_config
from thyra.resampling.mass_axis.tof_generator import (
    DEFAULT_BINS_PER_FWHM,
    PHI_NANOTOF_LAW,
    TOFAxisGenerator,
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
            # The width the 4,000 points have, not the 5 mDa default (D56).
            "mass_width_da": pytest.approx(1000.0 / 3_999),
            "reference_mz": 1000.0,
        }

    def test_a_tof_axis_records_its_law_too(self):
        planner = _planner(_PHI, target_bins=4_000)
        settled = planner.settle(lambda n_bins: None)

        assert settled.provenance["tof_a"] == pytest.approx(PHI_NANOTOF_LAW[0])
        assert settled.provenance["tof_b"] == pytest.approx(PHI_NANOTOF_LAW[1])
        # The count fixed the axis, so k is the one 4,000 bins give (D56).
        law = TOFAxisGenerator(*PHI_NANOTOF_LAW)
        k = 4_000 / (1e3 * (law.cumulative(1100.0) - law.cumulative(100.0)))
        assert settled.provenance["bins_per_fwhm"] == pytest.approx(k)
        assert settled.provenance["bins_per_fwhm"] != DEFAULT_BINS_PER_FWHM


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


class TestARequestIsHonouredWholeOrRefused:
    """Issue #445: two API requests used to be half-dropped without a word.

    A width with a bin count kept the width; half a TOF law took both
    halves from the detected instrument. Either way the store held an axis
    the caller had not asked for. The check sits in
    ``normalize_resampling_config``, which every conversion goes through
    for both the dict and the dataclass.
    """

    @pytest.mark.parametrize(
        "config",
        [
            {"target_bins": 4_000, "width_at_mz": 0.01},
            ResamplingConfig(target_bins=4_000, mass_width_da=0.01),
        ],
        ids=["dict", "dataclass"],
    )
    def test_a_width_with_a_bin_count_is_refused(self, config):
        with pytest.raises(
            ConversionRefused, match="both 'target_bins' and 'width_at_mz'"
        ):
            normalize_resampling_config(config)

    @pytest.mark.parametrize("given, missing", [("tof_a", "tof_b"), ("tof_b", "tof_a")])
    def test_half_a_law_is_refused_from_a_dict(self, given, missing):
        with pytest.raises(ConversionRefused, match=f"'{given}' without '{missing}'"):
            normalize_resampling_config({"axis_type": "tof", given: 0.0185})

    def test_half_a_law_is_refused_from_the_dataclass(self):
        with pytest.raises(ConversionRefused, match="'tof_b' without 'tof_a'"):
            normalize_resampling_config(
                ResamplingConfig(axis_type=AxisType.TOF, tof_b=PHI_NANOTOF_LAW[1])
            )

    @pytest.mark.parametrize(
        "config",
        [
            {"target_bins": 4_000},
            {"width_at_mz": 0.01},
            {"axis_type": "tof"},
            {"axis_type": "tof", "tof_a": 0.0185, "tof_b": 9.1e-6},
        ],
    )
    def test_a_whole_request_still_passes(self, config):
        normalize_resampling_config(config)


class TestANonFiniteRangeIsRefused:
    """D55: an infinite end used to reach the bin count and overflow."""

    @pytest.mark.parametrize(
        "mass_range", [(100.0, np.inf), (-np.inf, 1100.0), (100.0, np.nan)]
    )
    def test_from_the_source(self, mass_range):
        told = {"essential_metadata": {"mass_range": mass_range}}
        planner = AxisPlanner.from_metadata(told, ResamplingConfig())

        with pytest.raises(ConversionRefused, match="not finite") as refused:
            planner.resolve()
        assert "source's own mass range" in str(refused.value)

    def test_from_the_caller(self):
        planner = _planner(max_mz=float("inf"))

        with pytest.raises(ConversionRefused, match="not finite") as refused:
            planner.resolve()
        assert "--resample-max-mz" in str(refused.value)


class TestALowerEndAtZeroIsRefusedOnAnInverseAxis:
    """A bin count that divides by the lower end cannot start at 0 (D63)."""

    @pytest.mark.parametrize(
        "axis_type", [AxisType.REFLECTOR_TOF, AxisType.ORBITRAP, AxisType.FTICR]
    )
    def test_from_the_source(self, axis_type):
        told = {"essential_metadata": {"mass_range": (0.0, 1100.0)}}
        planner = AxisPlanner.from_metadata(told, ResamplingConfig(axis_type=axis_type))

        with pytest.raises(ConversionRefused, match="needs a lower end above 0") as r:
            planner.resolve()
        assert "source's own mass range" in str(r.value)

    @pytest.mark.parametrize(
        "axis_type", [AxisType.REFLECTOR_TOF, AxisType.ORBITRAP, AxisType.FTICR]
    )
    def test_from_the_caller(self, axis_type):
        planner = _planner(axis_type=axis_type, min_mz=0.0)

        with pytest.raises(ConversionRefused, match="needs a lower end above 0") as r:
            planner.resolve()
        assert "--resample-min-mz" in str(r.value)

    @pytest.mark.parametrize("axis_type", [AxisType.CONSTANT, AxisType.LINEAR_TOF])
    def test_a_law_that_lays_a_range_from_zero_is_left_alone(self, axis_type):
        told = {"essential_metadata": {"mass_range": (0.0, 1100.0)}}
        planner = AxisPlanner.from_metadata(told, ResamplingConfig(axis_type=axis_type))

        assert planner.resolve().min_mz == 0.0

    def test_a_tof_law_without_a_linear_term_is_refused(self):
        told = {"essential_metadata": {"mass_range": (0.0, 1100.0)}}
        planner = AxisPlanner.from_metadata(
            told, ResamplingConfig(axis_type=AxisType.TOF, tof_a=0.0, tof_b=1e-5)
        )

        with pytest.raises(ConversionRefused, match="needs a lower end above 0"):
            planner.resolve()

    def test_a_tof_law_with_a_linear_term_lays_from_zero(self):
        told = {"essential_metadata": {"mass_range": (0.0, 1100.0)}}
        planner = AxisPlanner.from_metadata(
            told, ResamplingConfig(axis_type=AxisType.TOF, tof_a=0.01, tof_b=1e-8)
        )

        plan = planner.resolve()
        assert plan.min_mz == 0.0
        assert plan.target_bins >= 2

    def test_a_non_finite_range_is_still_refused_as_not_finite(self):
        told = {"essential_metadata": {"mass_range": (0.0, np.inf)}}
        planner = AxisPlanner.from_metadata(
            told, ResamplingConfig(axis_type=AxisType.REFLECTOR_TOF)
        )

        with pytest.raises(ConversionRefused, match="not finite"):
            planner.resolve()


class TestBinsPerFwhmIsHonouredOrRefused:
    """C7 (D56): no setting is dropped, or recorded as applied when it was not."""

    @pytest.mark.parametrize(
        "config",
        [
            {"axis_type": "tof", "width_at_mz": 0.01, "bins_per_fwhm": 2.0},
            ResamplingConfig(
                axis_type=AxisType.TOF, mass_width_da=0.01, bins_per_fwhm=2.0
            ),
        ],
    )
    def test_beside_a_width(self, config):
        with pytest.raises(ConversionRefused, match="'bins_per_fwhm' and 'width_at"):
            normalize_resampling_config(config)

    @pytest.mark.parametrize(
        "config",
        [
            {"axis_type": "tof", "target_bins": 5000, "bins_per_fwhm": 2.0},
            ResamplingConfig(
                axis_type=AxisType.TOF, target_bins=5000, bins_per_fwhm=2.0
            ),
        ],
    )
    def test_beside_a_bin_count(self, config):
        with pytest.raises(ConversionRefused, match="'bins_per_fwhm' and 'target"):
            normalize_resampling_config(config)

    def test_on_a_named_axis_that_is_not_tof(self):
        planner = _planner(axis_type=AxisType.CONSTANT, bins_per_fwhm=2.0)

        with pytest.raises(ConversionRefused, match="resolved to constant"):
            planner.resolve()

    def test_on_a_detected_axis_that_is_not_tof(self):
        """Only the resolved axis type says whether it applies."""
        planner = _planner(_RAPIFLEX, bins_per_fwhm=2.0)

        with pytest.raises(ConversionRefused, match="sizes a tof axis"):
            planner.resolve()

    def test_alone_on_a_tof_axis_it_sizes_the_axis(self):
        planner = _planner(_PHI, bins_per_fwhm=2.0)
        settled = planner.settle(lambda n_bins: None)

        assert settled.provenance["bins_per_fwhm"] == 2.0


class TestTheRecordedWidthIsTheAxisWidth:
    """C8 (D56): the step records the width the axis has."""

    def _settle(self, mass_range, **config):
        told = {"essential_metadata": {"mass_range": mass_range}}
        planner = AxisPlanner.from_metadata(told, ResamplingConfig(**config))
        return planner.settle(lambda n_bins: None)

    def test_a_width_the_floor_raised(self):
        """1 Da over 500-510 is 10 bins; the floor lays 100 points 0.101 Da apart."""
        settled = self._settle(
            (500.0, 510.0), axis_type=AxisType.CONSTANT, mass_width_da=1.0
        )

        assert settled.provenance["target_bins"] == 100
        recorded = settled.provenance["mass_width_da"]
        np.testing.assert_allclose(np.diff(settled.axis), recorded, rtol=1e-9)
        assert recorded == pytest.approx(10.0 / 99)

    def test_a_width_that_set_the_count_is_recorded_as_given(self):
        settled = self._settle(_RANGE, axis_type=AxisType.CONSTANT, mass_width_da=0.05)

        assert settled.provenance["mass_width_da"] == 0.05

    def test_a_bin_count_on_a_constant_axis(self):
        settled = self._settle(_RANGE, axis_type=AxisType.CONSTANT, target_bins=100)

        recorded = settled.provenance["mass_width_da"]
        np.testing.assert_allclose(np.diff(settled.axis), recorded, rtol=1e-9)
        assert recorded == pytest.approx(1000.0 / 99)

    def test_a_floored_tof_axis_records_the_k_it_has(self):
        law = (0.5, 1e-8)
        settled = self._settle(
            (500.0, 510.0),
            axis_type=AxisType.TOF,
            tof_a=law[0],
            tof_b=law[1],
            bins_per_fwhm=0.01,
        )

        generator = TOFAxisGenerator(*law)
        k = 100 / (1e3 * (generator.cumulative(510.0) - generator.cumulative(500.0)))
        assert settled.provenance["target_bins"] == 100
        assert settled.provenance["bins_per_fwhm"] == pytest.approx(k)
        assert settled.provenance["mass_width_da"] == pytest.approx(
            generator.bin_width_at(1000.0, k)
        )


class TestAnEmptySourceRangeNamesTheWayForward:
    """A5 (D56)."""

    def test_every_peak_at_one_mz(self):
        told = {"essential_metadata": {"mass_range": (500.0, 500.0)}}
        planner = AxisPlanner.from_metadata(told, ResamplingConfig())

        with pytest.raises(ConversionRefused, match="is empty") as refused:
            planner.resolve()
        assert "--no-resample" in str(refused.value)
        assert "--resample-min-mz" in str(refused.value)

    def test_typed_bounds_are_not_told_to_drop_resampling(self):
        planner = _planner(min_mz=600.0, max_mz=600.0)

        with pytest.raises(ConversionRefused, match="is empty") as refused:
            planner.resolve()
        assert "--no-resample" not in str(refused.value)

    def test_a_non_finite_source_range_is_not_offered_the_raw_axis(self):
        """A non-finite m/z is refused on the raw axis too (D55)."""
        told = {"essential_metadata": {"mass_range": (100.0, np.inf)}}
        planner = AxisPlanner.from_metadata(told, ResamplingConfig())

        with pytest.raises(ConversionRefused, match="not finite") as refused:
            planner.resolve()
        assert "--no-resample" not in str(refused.value)


class TestAWidenedSourceRangeIsReportedWhereItSetsTheAxis:
    """K2 (D56): the warning follows what the axis does."""

    _WATERS = {
        "essential_metadata": {"mass_range": (99.9878, 1200.2158)},
        "acquisition_params": {
            "function_1_acq_range_start": 100.0,
            "function_1_acq_range_end": 1200.0,
        },
    }

    def _warnings(self, thyra_logs, metadata=None, **config):
        planner = AxisPlanner.from_metadata(
            metadata or self._WATERS, ResamplingConfig(**config)
        )
        with thyra_logs("thyra.resampling", logging.WARNING) as records:
            planner.resolve()
        return records.text

    def test_the_source_range_sets_the_axis(self, thyra_logs):
        text = self._warnings(thyra_logs, axis_type=AxisType.CONSTANT)

        assert "widened down to 99.9878 and up to 1200.2158" in text

    def test_typed_bounds_set_the_axis(self, thyra_logs):
        text = self._warnings(
            thyra_logs, axis_type=AxisType.CONSTANT, min_mz=700.0, max_mz=900.0
        )

        assert "widened" not in text

    def test_one_typed_bound_leaves_the_other_side(self, thyra_logs):
        text = self._warnings(thyra_logs, axis_type=AxisType.CONSTANT, min_mz=700.0)

        assert "widened up to 1200.2158" in text
        assert "down to" not in text

    def test_a_range_inside_the_acquisition_range_is_silent(self, thyra_logs):
        inside = dict(self._WATERS, essential_metadata={"mass_range": (100.0, 1200.0)})

        text = self._warnings(thyra_logs, inside, axis_type=AxisType.CONSTANT)

        assert "widened" not in text


class TestEveryLawRecordsTheWidthItLays:
    """C8 (D56) on the physics laws: the recorded width is the local step."""

    @pytest.mark.parametrize(
        "axis_type",
        [
            AxisType.LINEAR_TOF,
            AxisType.REFLECTOR_TOF,
            AxisType.ORBITRAP,
            AxisType.FTICR,
        ],
    )
    def test_a_bin_count(self, axis_type):
        told = {"essential_metadata": {"mass_range": _RANGE}}
        config = ResamplingConfig(axis_type=axis_type, target_bins=20_000)
        settled = AxisPlanner.from_metadata(told, config).settle(lambda n: None)

        reference_mz = settled.provenance["reference_mz"]
        i = int(np.searchsorted(settled.axis, reference_mz))
        local = float(settled.axis[i] - settled.axis[i - 1])
        assert settled.provenance["mass_width_da"] == pytest.approx(local, rel=1e-3)


class TestTheFloorIsReadFromTheCount:
    """The laid step can come out wider than asked; the floor still applied."""

    def test_a_width_just_too_wide_for_a_hundred_bins(self, thyra_logs):
        told = {"essential_metadata": {"mass_range": (500.0, 510.0)}}
        config = ResamplingConfig(axis_type=AxisType.CONSTANT, mass_width_da=0.1005)
        planner = AxisPlanner.from_metadata(told, config)

        with thyra_logs("thyra.resampling", logging.INFO) as records:
            settled = planner.settle(lambda n: None)

        assert settled.provenance["target_bins"] == 100
        assert settled.provenance["mass_width_da"] == pytest.approx(10.0 / 99)
        assert "floor of 100 bins" in records.text
        assert "finer than asked" not in records.text

    def test_an_inverted_source_range_is_not_told_to_set_bounds_around_it(self):
        told = {"essential_metadata": {"mass_range": (600.0, 500.0)}}
        planner = AxisPlanner.from_metadata(told, ResamplingConfig())

        with pytest.raises(ConversionRefused, match="is empty") as refused:
            planner.resolve()
        assert "around its m/z" not in str(refused.value)


class TestANonFiniteSettingIsRefused:
    """NaN or infinity means nothing as a width, bound or law (#503, D62).

    It used to reach the processing record, where JSON cannot hold it. The
    CLI refuses one already; the API now does too, in either spelling.
    """

    @pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
    @pytest.mark.parametrize(
        "key",
        ["gap_tolerance_da", "min_mz", "max_mz", "reference_mz", "width_at_mz"],
    )
    def test_in_a_dict(self, key, value):
        with pytest.raises(ConversionRefused, match="finite number"):
            normalize_resampling_config({key: value})

    @pytest.mark.parametrize(
        "field", ["gap_tolerance_da", "mass_width_da", "bins_per_fwhm"]
    )
    def test_in_a_dataclass(self, field):
        with pytest.raises(ConversionRefused, match=field):
            normalize_resampling_config(ResamplingConfig(**{field: float("nan")}))

    def test_finite_settings_pass(self):
        config = normalize_resampling_config(
            {"gap_tolerance_da": 0.01, "min_mz": 100, "max_mz": 1000.0}
        )
        assert config.gap_tolerance_da == 0.01
