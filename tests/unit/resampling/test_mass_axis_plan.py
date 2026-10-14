# tests/unit/resampling/test_mass_axis_plan.py
"""One mass axis for several sources, decided without a file (issue #449, D43).

The planner core takes what the detector chain is told about each source, so
every rule is exercised here with a metadata dict and no reader. What the
plan promises is the last class: the axis it describes is, bit for bit, the
axis every source's own conversion builds from the plan's settings.
"""

from __future__ import annotations

import json
import re

import numpy as np
import pytest

from thyra.errors import ConversionRefused
from thyra.mass_axis_plan import _plan_from_metadata, _Source
from thyra.metadata.schema import mass_axis_sha256
from thyra.resampling.axis_planner import AxisPlanner, normalize_resampling_config
from thyra.resampling.decision_tree import ResamplingDecisionTree
from thyra.resampling.mass_axis.tof_generator import (
    DEFAULT_BINS_PER_FWHM,
    MRT_TOF_LAW,
    TIMSTOF_TOF_LAW,
)
from thyra.resampling.plan import MassAxisPlan
from thyra.resampling.types import AxisType, ResamplingConfig, ResamplingMethod

_RAPIFLEX = {"format_specific": {"format": "Rapiflex"}}
_TIMS = {"GlobalMetadata": {"InstrumentName": "timsTOF fleX"}}
_MRT_CENTROID = {
    "format_specific": {"format": "Waters MassLynx raw", "is_mrt": True},
    "essential_metadata": {"spectrum_type": "centroid spectrum"},
}
_MRT_PROFILE = {
    "format_specific": {"format": "Waters MassLynx raw", "is_mrt": True},
    "essential_metadata": {"spectrum_type": "profile spectrum"},
}
_CENTROID = {"essential_metadata": {"spectrum_type": "centroid spectrum"}}


def _src(name, base, mass_range, **kwargs) -> _Source:
    """A source over ``mass_range`` with the detector facts of ``base``."""
    metadata = {k: dict(v) if isinstance(v, dict) else v for k, v in base.items()}
    essential = dict(metadata.get("essential_metadata", {}))
    essential["mass_range"] = mass_range
    metadata["essential_metadata"] = essential
    return _Source(name, metadata, **kwargs)


def _two(**overrides) -> MassAxisPlan:
    return _plan_from_metadata(
        [
            _src("a", _RAPIFLEX, (100.0, 1100.0)),
            _src("b", _RAPIFLEX, (200.0, 900.0)),
        ],
        **overrides,
    )


class TestRange:
    def test_the_range_is_the_intersection(self):
        plan = _two()

        assert (plan.min_mz, plan.max_mz) == (200.0, 900.0)
        assert plan.origins["min_mz"] == "intersection of the sources' ranges"
        assert plan.origins["max_mz"] == "intersection of the sources' ranges"

    def test_a_forced_bound_wins_and_the_other_is_still_intersected(self):
        plan = _two(min_mz=150.0)

        assert (plan.min_mz, plan.max_mz) == (150.0, 900.0)
        assert plan.origins["min_mz"] == "forced"
        assert plan.origins["max_mz"] == "intersection of the sources' ranges"

    def test_a_source_the_range_does_not_fit_is_named(self):
        plan = _two(min_mz=150.0)

        (warning,) = [w for w in plan.warnings if w.startswith("b:")]
        assert "[200, 900]" in warning and "always zero" in warning

    def test_ranges_with_nothing_in_common_are_refused(self):
        sources = [
            _src("a", _RAPIFLEX, (100.0, 300.0)),
            _src("b", _RAPIFLEX, (500.0, 900.0)),
        ]

        with pytest.raises(ConversionRefused, match="no m/z in common") as err:
            _plan_from_metadata(sources)

        assert "a:" in str(err.value) and "b:" in str(err.value)

    def test_a_source_with_no_stated_range_is_refused_unless_both_bounds_given(self):
        sources = [
            _src("a", _RAPIFLEX, (100.0, 1100.0)),
            _src("b", _RAPIFLEX, (0.0, 0.0), range_known=False),
        ]

        with pytest.raises(ConversionRefused, match="b"):
            _plan_from_metadata(sources)
        with pytest.raises(ConversionRefused, match="b"):
            _plan_from_metadata(sources, min_mz=150.0)

        plan = _plan_from_metadata(sources, min_mz=150.0, max_mz=800.0)
        assert (plan.min_mz, plan.max_mz) == (150.0, 800.0)

    def test_an_imzml_source_that_sets_a_bound_is_flagged(self):
        sources = [
            _src("a.imzML", _RAPIFLEX, (150.0, 900.0), source_format="imzml"),
            _src("b", _RAPIFLEX, (100.0, 1100.0), source_format="bruker"),
        ]

        plan = _plan_from_metadata(sources)

        (warning,) = [w for w in plan.warnings if "acquisition window" in w]
        assert warning.startswith("a.imzML:")
        assert "minimum and maximum" in warning

    def test_a_header_range_is_planned_as_read_and_its_origin_says_so(self):
        """The text step is applied when the header is read (D51); the
        planner takes the moved range as it stands and says where it came from."""
        sources = [
            _src(
                "a.imzML",
                _RAPIFLEX,
                (299.9, 600.1),
                source_format="imzml",
                range_is_header=True,
            )
        ]

        plan = _plan_from_metadata(sources)

        assert (plan.min_mz, plan.max_mz) == (299.9, 600.1)
        assert "header text moved out" in plan.origins["min_mz"]
        assert "header text moved out" in plan.origins["max_mz"]

    def test_a_source_that_is_not_a_header_range_says_nothing_of_one(self):
        plan = _plan_from_metadata([_src("a", _RAPIFLEX, (300.0, 600.0))])

        assert (plan.min_mz, plan.max_mz) == (300.0, 600.0)
        assert "header" not in plan.origins["min_mz"]

    def test_a_forced_bound_is_used_as_given(self):
        sources = [
            _src("a.imzML", _RAPIFLEX, (300.0, 600.0), range_is_header=True),
        ]

        plan = _plan_from_metadata(sources, min_mz=300.0, max_mz=600.0)

        assert (plan.min_mz, plan.max_mz) == (300.0, 600.0)
        assert plan.origins["min_mz"] == "forced"

    def test_the_imzml_edge_warning_fires_on_a_header_range(self):
        sources = [
            _src(
                "a.imzML",
                _RAPIFLEX,
                (150.0, 900.0),
                source_format="imzml",
                range_is_header=True,
            ),
            _src("b", _RAPIFLEX, (100.0, 1100.0), source_format="bruker"),
        ]

        plan = _plan_from_metadata(sources)

        (warning,) = [w for w in plan.warnings if "acquisition window" in w]
        assert warning.startswith("a.imzML:")
        assert "minimum and maximum" in warning

    def test_a_forced_min_above_the_planned_max_says_so(self):
        with pytest.raises(ConversionRefused) as err:
            _two(min_mz=950.0)

        message = str(err.value)
        assert "forced min_mz 950" in message and "planned max_mz 900" in message
        assert "Pass a lower min_mz, or pass max_mz too" in message

    def test_a_forced_max_below_the_planned_min_says_so(self):
        with pytest.raises(ConversionRefused) as err:
            _two(max_mz=150.0)

        message = str(err.value)
        assert "forced max_mz 150" in message and "planned min_mz 200" in message
        assert "Pass a higher max_mz, or pass min_mz too" in message

    def test_a_source_that_is_not_imzml_is_not_flagged_for_it(self):
        sources = [
            _src("a", _RAPIFLEX, (150.0, 900.0), source_format="phi"),
            _src("b", _RAPIFLEX, (100.0, 1100.0), source_format="imzml"),
        ]

        plan = _plan_from_metadata(sources)

        assert not [w for w in plan.warnings if "acquisition window" in w]

    def test_a_forced_bound_is_not_blamed_on_the_imzml_source(self):
        sources = [
            _src("a.imzML", _RAPIFLEX, (150.0, 900.0), source_format="imzml"),
            _src("b", _RAPIFLEX, (100.0, 1100.0)),
        ]

        plan = _plan_from_metadata(sources, min_mz=150.0, max_mz=900.0)

        assert not [w for w in plan.warnings if "acquisition window" in w]

    def test_a_forced_empty_range_is_refused(self):
        with pytest.raises(ConversionRefused, match="empty"):
            _two(min_mz=500.0, max_mz=400.0)

    def test_no_source_is_refused(self):
        with pytest.raises(ConversionRefused, match="at least one"):
            _plan_from_metadata([])

    def test_one_source_pins_its_own_axis(self):
        plan = _plan_from_metadata([_src("a", _RAPIFLEX, (100.0, 1100.0))])

        assert (plan.min_mz, plan.max_mz) == (100.0, 1100.0)
        assert plan.sources == ("a",)


class TestTheMemoryGuard:
    def test_an_axis_too_wide_for_memory_is_refused_unbuilt(self, monkeypatch):
        """The plan lays its axis through the guard a conversion uses (#473)."""
        from thyra.resampling import memory_guard

        # About 100,000 bins against 1 MB "free": refused by the guard, and
        # small enough that a plan laid without the guard fails fast instead
        # of allocating.
        monkeypatch.setattr(memory_guard, "available_memory_gb", lambda: 0.001)
        source = ("a", {"essential_metadata": {"mass_range": (100.0, 200.0)}})

        with pytest.raises(ConversionRefused, match=r"bins") as err:
            _plan_from_metadata([source], axis_type="constant", width_at_mz=1e-3)

        message = str(err.value)
        assert re.search(r"100,00\d bins", message), message
        assert "width_at_mz" in message and "--resample-bins" not in message


class TestAxisType:
    def _mixed(self):
        return [
            _src("a", _RAPIFLEX, (100.0, 1100.0)),
            _src("b", _CENTROID, (100.0, 1100.0)),
        ]

    def test_detectors_that_disagree_are_refused_naming_each(self):
        with pytest.raises(ConversionRefused, match="axis_type") as err:
            _plan_from_metadata(self._mixed())

        message = str(err.value)
        assert "a: constant" in message and "b: reflector_tof" in message

    def test_a_forced_type_is_accepted(self):
        plan = _plan_from_metadata(self._mixed(), axis_type="reflector_tof")

        assert plan.axis_type is AxisType.REFLECTOR_TOF
        assert plan.origins["axis_type"] == "forced"

    def test_agreement_is_recorded(self):
        assert _two().origins["axis_type"] == "all sources agree"

    def test_an_unknown_name_is_refused_like_the_resampling_config(self):
        with pytest.raises(ConversionRefused, match="axis_type"):
            _two(axis_type="sideways")

    def test_a_law_without_a_tof_axis_is_refused(self):
        with pytest.raises(ConversionRefused, match="tof axis"):
            _two(tof_a=0.02, tof_b=1e-5)


class TestTofLaw:
    def _mrt(self, **kwargs):
        sources = [
            _src("a", _MRT_CENTROID, (100.0, 1100.0)),
            _src("b", _MRT_CENTROID, (200.0, 900.0)),
        ]
        return _plan_from_metadata(sources, **kwargs)

    def test_an_agreed_law_is_used(self):
        plan = self._mrt()

        assert plan.axis_type is AxisType.TOF
        assert (plan.tof_a, plan.tof_b) == MRT_TOF_LAW
        assert plan.bins_per_fwhm == DEFAULT_BINS_PER_FWHM
        assert plan.width_at_mz is None
        assert "agree" in plan.origins["tof_law"]
        assert plan.origins["bins_per_fwhm"] == "default"

    def test_different_laws_are_refused(self):
        sources = [
            _src("a", _MRT_CENTROID, (100.0, 1100.0)),
            _src("b", _TIMS, (100.0, 1100.0)),
        ]

        with pytest.raises(ConversionRefused, match="different TOF width laws"):
            _plan_from_metadata(sources, axis_type="tof")

    def test_a_forced_pair_wins_over_different_laws(self):
        sources = [
            _src("a", _MRT_CENTROID, (100.0, 1100.0)),
            _src("b", _TIMS, (100.0, 1100.0)),
        ]

        plan = _plan_from_metadata(sources, axis_type="tof", tof_a=0.03, tof_b=2e-5)

        assert (plan.tof_a, plan.tof_b) == (0.03, 2e-5)
        assert plan.origins["tof_law"] == "forced"

    @pytest.mark.parametrize("half", [{"tof_a": 0.02}, {"tof_b": 1e-5}])
    def test_half_a_pair_is_refused(self, half):
        with pytest.raises(ConversionRefused, match="tof_a and tof_b"):
            self._mrt(**half)

    def test_no_declared_law_is_refused(self):
        with pytest.raises(ConversionRefused, match="needs tof_a and tof_b"):
            _two(axis_type="tof")

    def test_a_declaring_source_sets_the_law_for_one_that_does_not(self):
        sources = [
            _src("a", _TIMS, (100.0, 1100.0)),
            _src("b", _RAPIFLEX, (100.0, 1100.0)),
        ]

        plan = _plan_from_metadata(sources, axis_type="tof")

        assert (plan.tof_a, plan.tof_b) == TIMSTOF_TOF_LAW

    def test_a_forced_width_gives_the_bins_per_fwhm(self):
        plan = self._mrt(width_at_mz=0.004, reference_mz=800.0)

        assert plan.bins_per_fwhm != DEFAULT_BINS_PER_FWHM
        assert "forced width" in plan.origins["bins_per_fwhm"]
        assert plan.width_at_mz is None
        assert "800" in plan.origins["width_at_mz"]

    def test_a_width_and_bins_per_fwhm_each_fix_the_count(self):
        with pytest.raises(ConversionRefused, match="bins_per_fwhm"):
            self._mrt(width_at_mz=0.004, bins_per_fwhm=4.0)

    def test_a_tof_config_holds_the_law_and_never_a_width(self):
        config = self._mrt().to_resampling_config()

        assert config["axis_type"] == "tof"
        assert {"tof_a", "tof_b", "bins_per_fwhm"} <= set(config)
        assert "width_at_mz" not in config and "target_bins" not in config


class TestWidth:
    def test_the_default_width_is_the_axis_types_own(self):
        plan = _two()

        assert plan.width_at_mz == 0.005
        assert plan.reference_mz == 1000.0

    def test_linear_tof_defaults_to_its_own_reference(self):
        plan = _two(axis_type="linear_tof")

        assert (plan.width_at_mz, plan.reference_mz) == (0.017, 300.0)

    def test_a_forced_width_is_used_and_recorded(self):
        plan = _two(width_at_mz=0.02, reference_mz=500.0)

        assert (plan.width_at_mz, plan.reference_mz) == (0.02, 500.0)
        assert plan.origins["width_at_mz"] == "forced"
        assert plan.origins["reference_mz"] == "forced"

    def test_the_widest_is_picked_after_scaling_to_one_m_z(self, monkeypatch):
        """Declared at different m/z, the order of the raw numbers is not the order.

        x: 0.010 Da at m/z 100 is 0.100 Da at m/z 1000 on a reflector axis
        (width grows with m/z), y: 0.050 Da at 1000. A maximum over the raw
        numbers picks y.
        """
        declared = {"x": (0.010, 100.0), "y": (0.050, 1000.0)}

        def fake(self, metadata):
            return declared[metadata["who"]]

        monkeypatch.setattr(ResamplingDecisionTree, "select_reference_width", fake)
        sources = [
            _src("x", {"who": "x"}, (100.0, 1100.0)),
            _src("y", {"who": "y"}, (100.0, 1100.0)),
        ]

        plan = _plan_from_metadata(sources, axis_type="reflector_tof")

        assert plan.width_at_mz == pytest.approx(0.100)
        assert plan.origins["width_at_mz"].startswith("widest declared width")
        assert "x: 0.1 Da" in plan.origins["width_at_mz"]
        assert "y: 0.05 Da" in plan.origins["width_at_mz"]
        (warning,) = [w for w in plan.warnings if "bin width" in w]
        assert warning.startswith("y:") and "loses resolution" in warning

    @pytest.mark.parametrize(
        "axis, exponent",
        [
            ("constant", 0.0),
            ("linear_tof", 0.5),
            ("reflector_tof", 1.0),
            ("orbitrap", 1.5),
            ("fticr", 2.0),
        ],
    )
    def test_a_width_scales_with_the_law_of_the_axis(self, monkeypatch, axis, exponent):
        monkeypatch.setattr(
            ResamplingDecisionTree,
            "select_reference_width",
            lambda self, metadata: (0.010, 500.0),
        )
        plan = _plan_from_metadata(
            [_src("x", {}, (100.0, 1100.0))], axis_type=axis, reference_mz=1000.0
        )

        assert plan.width_at_mz == pytest.approx(0.010 * 2.0**exponent)

    def test_a_narrower_source_is_named(self):
        sources = [
            _src("fine", _RAPIFLEX, (100.0, 1100.0)),
            _src("other", _CENTROID, (100.0, 1100.0)),
        ]

        plan = _plan_from_metadata(sources, axis_type="constant", width_at_mz=0.05)

        by_name = {w.split(":")[0]: w for w in plan.warnings if "bin width" in w}
        assert "loses resolution" in by_name["fine"]
        # tic_preserving rebins bins wider than its samples (D50), so a
        # coarse plan no longer costs it peaks and draws no extra caution.
        assert "tic_preserving" not in by_name["fine"]

    def test_a_forced_method_is_the_plans_method(self):
        sources = [_src("fine", _RAPIFLEX, (100.0, 1100.0))]

        plan = _plan_from_metadata(sources, width_at_mz=0.05, method="nearest_neighbor")

        assert plan.method is ResamplingMethod.NEAREST_NEIGHBOR

    def test_a_forced_width_at_or_below_every_sources_own_warns_nobody(self):
        plan = _two(width_at_mz=0.005)

        assert not [w for w in plan.warnings if "bin width" in w]


class TestSpectrumTypes:
    def test_mixed_types_are_warned_about(self):
        sources = [
            _src("a", _MRT_CENTROID, (100.0, 1100.0)),
            _src("b", _MRT_PROFILE, (100.0, 1100.0)),
        ]

        plan = _plan_from_metadata(sources, axis_type="tof")

        (warning,) = [w for w in plan.warnings if "spectrum types" in w]
        assert "a: centroid spectrum" in warning and "b: profile spectrum" in warning

    def test_the_same_type_is_quiet(self):
        assert not [w for w in _two().warnings if "spectrum types" in w]


class TestThePlanItself:
    def test_the_config_never_holds_bins_and_a_width(self):
        config = _two().to_resampling_config()

        assert "width_at_mz" in config and "target_bins" not in config
        assert config["mass_axis_sha256"] == _two().sha256
        assert "method" not in config

    def test_a_forced_method_travels(self):
        assert _two(method="tic_preserving").to_resampling_config()["method"] == (
            "tic_preserving"
        )

    def test_the_config_normalises_and_carries_the_fingerprint(self):
        plan = _two()

        config = normalize_resampling_config(plan.to_resampling_config())

        assert config.mass_axis_sha256 == plan.sha256
        assert config.target_bins is None

    def test_a_plan_normalises_to_the_same_config(self):
        plan = _two()

        assert normalize_resampling_config(plan) == normalize_resampling_config(
            plan.to_resampling_config()
        )

    @pytest.mark.parametrize("bad", ["abc", "A" * 64, "g" * 64, "a" * 63, 7])
    def test_a_fingerprint_that_is_not_lower_case_hex_is_refused(self, bad):
        with pytest.raises(ConversionRefused, match="mass_axis_sha256"):
            normalize_resampling_config({"mass_axis_sha256": bad})

    def test_a_config_object_with_a_bad_fingerprint_is_refused(self):
        with pytest.raises(ConversionRefused, match="mass_axis_sha256"):
            normalize_resampling_config(ResamplingConfig(mass_axis_sha256="abc"))

    def test_the_dict_round_trips(self):
        plan = _two(method="tic_preserving")

        again = MassAxisPlan.from_dict(json.loads(json.dumps(plan.to_dict())))

        assert again == plan

    def test_a_tof_plan_round_trips(self):
        plan = _plan_from_metadata(
            [_src("a", _MRT_CENTROID, (100.0, 1100.0))], width_at_mz=0.004
        )

        assert MassAxisPlan.from_dict(plan.to_dict()) == plan

    def test_save_and_load_round_trip_as_utf8(self, tmp_path):
        plan = _two()
        path = tmp_path / "axis.json"

        plan.save(path)

        assert MassAxisPlan.load(path) == plan
        assert json.loads(path.read_text(encoding="utf-8"))["format"] == (
            "thyra-mass-axis-plan"
        )

    def test_a_plan_file_in_utf8_with_a_mark_or_utf16_loads(self, tmp_path):
        plan = _two()
        text = json.dumps(plan.to_dict(), ensure_ascii=False)
        encoded = {
            "utf-8-sig": text.encode("utf-8-sig"),
            "utf-16-le": text.encode("utf-16"),
            "utf-16-be": b"\xfe\xff" + text.encode("utf-16-be"),
        }
        for name, raw in encoded.items():
            path = tmp_path / f"{name}.json"
            path.write_bytes(raw)

            assert MassAxisPlan.load(path) == plan, name

    def test_a_plan_file_in_another_encoding_is_refused_in_one_line(self, tmp_path):
        data = _two().to_dict()
        data["sources"] = ["caf\u00e9"]
        path = tmp_path / "latin.json"
        path.write_bytes(json.dumps(data, ensure_ascii=False).encode("latin-1"))

        with pytest.raises(ConversionRefused, match="not UTF-8, UTF-16 or UTF-32") as e:
            MassAxisPlan.load(path)

        assert "\n" not in str(e.value) and "latin.json" in str(e.value)

    @pytest.mark.parametrize("version", [True, 1.0, "1", 2, None])
    def test_a_version_that_is_not_the_whole_number_one_is_refused(self, version):
        data = _two().to_dict()
        data["version"] = version

        with pytest.raises(ConversionRefused, match="version"):
            MassAxisPlan.from_dict(data)

    def test_a_key_no_plan_has_is_refused_by_name(self):
        data = _two().to_dict()
        data["target_bins"] = 5000

        with pytest.raises(ConversionRefused, match="'target_bins'"):
            MassAxisPlan.from_dict(data)

    def test_an_edited_bin_count_is_refused(self):
        plan = _two()
        data = plan.to_dict()
        data["n_bins"] = plan.n_bins + 1

        with pytest.raises(ConversionRefused, match="does not match its settings") as e:
            MassAxisPlan.from_dict(data)

        assert f"{plan.n_bins:,}" in str(e.value)

    @pytest.mark.parametrize("key", ["tof_a", "tof_b", "bins_per_fwhm"])
    def test_a_tof_plan_without_its_law_is_refused_naming_the_key(self, key):
        data = _two(axis_type="tof", tof_a=0.03, tof_b=2e-5).to_dict()
        data[key] = None

        with pytest.raises(ConversionRefused, match=repr(key)) as e:
            MassAxisPlan.from_dict(data)

        assert "--" not in str(e.value)

    def test_a_plan_without_its_width_is_refused_naming_the_key(self):
        data = _two().to_dict()
        data["width_at_mz"] = None

        with pytest.raises(ConversionRefused, match="'width_at_mz'"):
            MassAxisPlan.from_dict(data)

    def test_loading_a_plan_logs_no_axis_calculation(self, thyra_logs):
        data = _two().to_dict()

        with thyra_logs("thyra.resampling.axis_planner") as records:
            MassAxisPlan.from_dict(data)

        assert not [r for r in records if "bins" in r.getMessage()]

    def test_an_inverted_range_is_refused_naming_the_keys(self):
        data = _two().to_dict()
        data["min_mz"], data["max_mz"] = data["max_mz"], data["min_mz"]

        with pytest.raises(ConversionRefused, match="'min_mz'.*below its 'max_mz'"):
            MassAxisPlan.from_dict(data)

    def test_a_dict_that_is_not_a_plan_is_refused(self):
        data = _two().to_dict()
        data["format"] = "something-else"

        with pytest.raises(ConversionRefused, match="Not a mass axis plan"):
            MassAxisPlan.from_dict(data)
        del data["format"]
        with pytest.raises(ConversionRefused, match="Not a mass axis plan"):
            MassAxisPlan.from_dict(data)

    def test_a_file_that_is_not_json_is_refused(self, tmp_path):
        path = tmp_path / "axis.json"
        path.write_text("not json", encoding="utf-8")

        with pytest.raises(ConversionRefused, match="axis.json"):
            MassAxisPlan.load(path)

    def test_a_missing_file_is_refused(self, tmp_path):
        with pytest.raises(ConversionRefused, match="nope.json"):
            MassAxisPlan.load(tmp_path / "nope.json")

    def test_a_plan_with_a_bad_fingerprint_is_refused_on_load(self):
        data = _two().to_dict()
        data["sha256"] = "x"

        with pytest.raises(ConversionRefused, match="SHA-256"):
            MassAxisPlan.from_dict(data)

    @pytest.mark.parametrize(
        "key, value",
        [
            ("min_mz", "abc"),
            ("min_mz", [1]),
            ("min_mz", 0),
            ("max_mz", float("inf")),
            ("reference_mz", True),
            ("width_at_mz", -0.1),
            ("n_bins", "x"),
            ("n_bins", 1.5),
            ("n_bins", 1),
            ("method", "none"),
            ("method", ["tic_preserving"]),
            ("origins", [1, 2]),
            ("warnings", "one"),
            ("sources", [1]),
        ],
    )
    def test_a_hand_edited_field_is_refused_not_raised(self, key, value):
        """A plan file is JSON a person may edit; a bad value is a refusal."""
        data = _two().to_dict()
        data[key] = value

        with pytest.raises(ConversionRefused, match="mass axis plan"):
            MassAxisPlan.from_dict(data)

    def test_a_zero_tof_term_is_a_valid_law(self):
        plan = _plan_from_metadata(
            [_src("a", _MRT_CENTROID, (100.0, 1100.0))],
            axis_type="tof",
            tof_a=0.0185,
            tof_b=0.0,
        )

        assert MassAxisPlan.from_dict(plan.to_dict()) == plan


_AXES = [
    ("constant", _RAPIFLEX, {}),
    ("linear_tof", _RAPIFLEX, {}),
    ("reflector_tof", _CENTROID, {}),
    ("orbitrap", _RAPIFLEX, {}),
    ("fticr", _RAPIFLEX, {"width_at_mz": 0.01}),
    ("tof", _MRT_CENTROID, {"tof_a": 0.03, "tof_b": 2e-5}),
]


class TestTheAxisIsWhatEverySourceBuilds:
    @pytest.mark.parametrize("axis_type, base, extra", _AXES, ids=[a[0] for a in _AXES])
    def test_bit_for_bit(self, axis_type, base, extra):
        sources = [
            _src("a", base, (100.0, 1100.0)),
            _src("b", base, (200.0, 900.0)),
            _src("c", _TIMS, (50.0, 1500.0)),
        ]
        plan = _plan_from_metadata(sources, axis_type=axis_type, **extra)
        config = normalize_resampling_config(plan.to_resampling_config())

        axes = [
            AxisPlanner.from_metadata(s.metadata, config)
            .settle(lambda n_bins: None)
            .axis
            for s in sources
        ]

        first = axes[0]
        assert len(first) == plan.n_bins
        assert mass_axis_sha256(first) == plan.sha256
        for other in axes[1:]:
            assert np.array_equal(first, other)
