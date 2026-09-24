"""The ``uns`` block built from arguments, with no converter anywhere.

That is the point of the split (issue #351): the block used to be
fourteen converter methods reading twelve converter attributes, so the
only way to ask what a given input produces was to build a converter and
run a conversion. Here the assembler is handed a stub reader and an
:class:`UnsContext` and asked directly, which is also the shape a second
write path would use if one ever existed again.
"""

import json
import logging
from types import SimpleNamespace
from typing import Any, Dict, List, Optional

import numpy as np
import pytest

from thyra.converters.spatialdata.uns_assembler import UnsAssembler, UnsContext
from thyra.core.base_converter import PixelSizeSource
from thyra.metadata.schema import MSI_METADATA_UNS_KEY
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata

_REGIONS = [{"region_number": 1, "n_spectra": 9}]


def _essential() -> EssentialMetadata:
    return EssentialMetadata(
        dimensions=(3, 3, 1),
        coordinate_bounds=(0.0, 2.0, 0.0, 2.0),
        mass_range=(100.0, 1000.0),
        pixel_size=(10.0, 10.0),
        n_spectra=9,
        total_peaks=90,
        source_path="stub.imzML",
        spectrum_type="centroid spectrum",
    )


def _comprehensive() -> ComprehensiveMetadata:
    return ComprehensiveMetadata(
        essential=_essential(),
        format_specific={"files": ["stub.imzML", "stub.ibd"], "offsets": [1, 2, 3]},
        acquisition_params={"polarity": "positive"},
        instrument_info={"model": "stub"},
        raw_metadata={"cvParams": [{"name": "MS1 spectrum"}]},
    )


class _StubReader:
    """Everything the assembler asks a reader for, and nothing else."""

    file_type = "imzml"

    def __init__(
        self,
        *,
        comprehensive: Optional[ComprehensiveMetadata] = None,
        fragmentation: Any = None,
        applied_calibration: Optional[Dict[str, Any]] = None,
    ) -> None:
        self._comprehensive = (
            _comprehensive() if comprehensive is None else comprehensive
        )
        self._fragmentation = fragmentation
        self._applied_calibration = applied_calibration
        self.fragmentation_reads = 0

    def get_comprehensive_metadata(self) -> ComprehensiveMetadata:
        return self._comprehensive

    def get_applied_mz_calibration(self) -> Optional[Dict[str, Any]]:
        return self._applied_calibration

    @property
    def has_ion_mobility(self) -> bool:
        return False

    def get_mobility_axis(self) -> Any:
        return None

    @property
    def has_fragmentation(self) -> bool:
        return self._fragmentation is not None

    def get_fragmentation(self) -> Any:
        self.fragmentation_reads += 1
        return self._fragmentation


class _StubSchedule:
    """The five things the assembler reads off a fragmentation schedule."""

    is_msms = True

    def __init__(self, n_windows: int) -> None:
        self.windows = tuple(
            SimpleNamespace(target=300.0 + 10 * i) for i in range(n_windows)
        )

    @property
    def n_precursors(self) -> int:
        return len({w.target for w in self.windows})

    @property
    def merges_precursors(self) -> bool:
        return self.n_precursors > 1

    def to_uns(self) -> Dict[str, Any]:
        return {"ms_level": 2, "targets": [w.target for w in self.windows]}

    def to_extractor_report(self) -> Dict[str, Any]:
        return {"present": True, "ms_level": 2}


def _assembler(reader: Optional[_StubReader] = None) -> UnsAssembler:
    return UnsAssembler(
        reader if reader is not None else _StubReader(),
        pixel_size_xy=lambda: (10.0, 20.0),
        pixel_size_detection_info={"source_format": "imzml"},
        resampling_config=None,
    )


def _context(**overrides: Any) -> UnsContext:
    defaults: Dict[str, Any] = {
        "mobility_table_key": None,
        "msms_table_key": None,
        "mobility_grid": None,
        "resolved_resampling_plan": None,
        "region_info": _REGIONS,
        "pixel_size_source": PixelSizeSource.AUTO_DETECTED,
        "mobility_heatmap": lambda: None,
    }
    defaults.update(overrides)
    return UnsContext(**defaults)


class TestTheBlockFromArgumentsAlone:
    def test_essential_metadata_carries_the_dimensions_and_the_version(self):
        uns = _assembler().build(_context())

        essential = uns["essential_metadata"]
        # Tuples become lists: zarr cannot store a tuple.
        assert essential["dimensions"] == [3, 3, 1]
        assert essential["mass_range"] == [100.0, 1000.0]
        assert essential["source_path"] == "stub.imzML"
        assert essential["spectrum_type"] == "centroid spectrum"
        assert essential["thyra_version"]

    def test_a_string_list_is_json_while_a_numeric_one_stays_a_list(self):
        uns = _assembler().build(_context())

        assert json.loads(uns["format_specific"]["files"]) == [
            "stub.imzML",
            "stub.ibd",
        ]
        assert uns["format_specific"]["offsets"] == [1, 2, 3]
        assert json.loads(uns["raw_metadata"]["cvParams"]) == [{"name": "MS1 spectrum"}]

    def test_regions_is_the_json_of_the_context_summary(self):
        uns = _assembler().build(_context())

        assert json.loads(uns["regions"]) == _REGIONS

    def test_no_region_summary_omits_the_key_rather_than_writing_an_empty_one(self):
        uns = _assembler().build(_context(region_info=None))

        assert "regions" not in uns

    def test_the_msi_metadata_block_names_the_siblings_the_context_names(self):
        uns = _assembler().build(
            _context(mobility_table_key="stub_z0_mobility", msms_table_key=None)
        )

        analysis = uns[MSI_METADATA_UNS_KEY]["ms_analysis"]
        assert analysis["ion_mobility"]["resolved_table"] == "stub_z0_mobility"
        assert "resolved_table" not in analysis.get("fragmentation", {})

    def test_the_pixel_size_is_the_pair_the_accessor_returns(self):
        uns = _assembler().build(_context())

        pixel_size = uns[MSI_METADATA_UNS_KEY]["ms_analysis"]["pixel_size_um"]
        assert (pixel_size["x"], pixel_size["y"]) == (10.0, 20.0)
        assert uns[MSI_METADATA_UNS_KEY]["provenance"]["source_format"] == "imzml"

    def test_the_key_order_is_the_store_order(self):
        """The order the sections are composed in is the store's key order."""
        reader = _StubReader(fragmentation=_StubSchedule(2))
        uns = _assembler(reader).build(
            _context(mobility_heatmap=lambda: {"counts": np.zeros((2, 2))})
        )

        assert list(uns) == [
            "essential_metadata",
            "format_specific",
            "acquisition_params",
            "instrument_info",
            "raw_metadata",
            "regions",
            MSI_METADATA_UNS_KEY,
            "mobility_heatmap",
            "msms_schedule",
        ]


class TestTheProcessingStep:
    def test_the_resolved_plan_wins_over_the_requested_config(self):
        config = SimpleNamespace(method="auto", target_bins=None, min_mz=None)
        assembler = UnsAssembler(
            _StubReader(),
            pixel_size_xy=lambda: (10.0, 10.0),
            pixel_size_detection_info=None,
            resampling_config=config,
        )

        steps = assembler._processing_provenance(
            _context(
                resolved_resampling_plan={
                    "method": "tic_preserving",
                    "target_bins": 4096,
                    "min_mz": None,
                }
            )
        )

        assert [s.name for s in steps] == ["conversion", "mass axis resampling"]
        parameters = steps[1].parameters
        assert parameters["method"] == "tic_preserving"
        assert parameters["target_bins"] == 4096
        # Unset on both sides: dropped, not serialised as a null.
        assert "min_mz" not in parameters


class TestTheCalibrationStep:
    """What the reader applied is a step, bound to MS:1001485."""

    def test_a_reader_that_applies_none_records_no_step(self):
        steps = _assembler()._processing_provenance(_context())

        assert [s.name for s in steps] == ["conversion"]

    def test_the_readers_parameters_become_the_step(self):
        reader = _StubReader(applied_calibration={"use_recalibrated_state": True})

        steps = _assembler(reader)._processing_provenance(_context())

        assert [s.name for s in steps] == ["conversion", "m/z calibration"]
        step = steps[1]
        assert step.action_term is not None
        assert step.action_term.accession == "MS:1001485"
        assert step.action_term.name == "m/z calibration"
        assert step.parameters == {"use_recalibrated_state": True}
        assert step.software.name == "thyra"

    def test_it_comes_before_resampling_because_it_happens_first(self):
        assembler = UnsAssembler(
            _StubReader(applied_calibration={"calibration": "appended"}),
            pixel_size_xy=lambda: (10.0, 10.0),
            pixel_size_detection_info=None,
            resampling_config=SimpleNamespace(method="nearest_neighbor"),
        )

        steps = assembler._processing_provenance(_context())

        assert [s.name for s in steps] == [
            "conversion",
            "m/z calibration",
            "mass axis resampling",
        ]

    def test_the_step_reaches_the_block(self):
        reader = _StubReader(applied_calibration={"calibration": "header"})

        block = _assembler(reader).build(_context())[MSI_METADATA_UNS_KEY]

        steps = json.loads(block["processing"])
        assert steps[1] == {
            "name": "m/z calibration",
            "action_term": {"accession": "MS:1001485", "name": "m/z calibration"},
            "software": steps[0]["software"],
            "parameters": {"calibration": "header"},
        }


class TestTheMobilityHeatmap:
    def test_a_block_is_recorded_and_the_callable_runs_once_per_build(self):
        calls: List[int] = []

        def heatmap() -> Dict[str, Any]:
            calls.append(1)
            return {"counts": np.arange(4.0)}

        uns = _assembler().build(_context(mobility_heatmap=heatmap))

        assert uns["mobility_heatmap"]["counts"].tolist() == [0.0, 1.0, 2.0, 3.0]
        assert len(calls) == 1

    def test_nothing_is_recorded_when_the_callable_says_there_is_none(self):
        uns = _assembler().build(_context(mobility_heatmap=lambda: None))

        assert "mobility_heatmap" not in uns


class TestTheFragmentationSchedule:
    def test_it_is_read_once_per_conversion_however_many_tables_are_built(self):
        reader = _StubReader(fragmentation=_StubSchedule(3))
        assembler = _assembler(reader)

        assembler.build(_context())
        assembler.build(_context())

        assert reader.fragmentation_reads == 1

    def test_the_chimera_warning_is_said_once(self, thyra_logs):
        reader = _StubReader(fragmentation=_StubSchedule(3))
        assembler = _assembler(reader)

        with thyra_logs("thyra.converters.spatialdata", logging.WARNING) as records:
            assembler.build(_context())
            assembler.build(_context())

        merged = [r for r in records if "isolates 3 precursors" in r.message]
        assert len(merged) == 1
        assert "msms_schedule" in merged[0].message

    def test_a_single_precursor_is_not_warned_about(self, thyra_logs):
        assembler = _assembler(_StubReader(fragmentation=_StubSchedule(1)))

        with thyra_logs("thyra.converters.spatialdata", logging.WARNING) as records:
            assembler.build(_context())

        assert not [r for r in records if "precursors per pixel" in r.message]

    def test_the_schedule_block_names_the_resolved_table(self):
        reader = _StubReader(fragmentation=_StubSchedule(2))

        uns = _assembler(reader).build(_context(msms_table_key="stub_z0_msms"))

        assert uns["msms_schedule"]["resolved_table"] == "stub_z0_msms"
        # A purely numeric list is kept as a list; only the rest is JSON.
        assert uns["msms_schedule"]["targets"] == [300.0, 310.0]


class TestUnnamingADeclinedSibling:
    @staticmethod
    def _table(key: str) -> SimpleNamespace:
        return SimpleNamespace(
            uns={
                "mobility_axis": {"resolved_table": key, "n_channels": 4},
                MSI_METADATA_UNS_KEY: {
                    "ms_analysis": {
                        "ion_mobility": {
                            "present": True,
                            "resolved_table": key,
                            "grid": {"n_bins": 4},
                        }
                    }
                },
            }
        )

    def test_every_table_named_loses_the_pointer(self):
        summed = self._table("stub_z0_mobility")
        sibling = self._table("stub_z0_mobility")
        tables = {"stub_z0": summed, "stub_z0_msms": sibling}

        _assembler().unname_declined_siblings(
            tables, ("stub_z0", "stub_z0_msms"), ["ion_mobility"]
        )

        for table in (summed, sibling):
            assert "resolved_table" not in table.uns["mobility_axis"]
            block = table.uns[MSI_METADATA_UNS_KEY]["ms_analysis"]["ion_mobility"]
            # forget_resolved_table drops the grid with the pointer, and
            # leaves the acquisition fact alone.
            assert "resolved_table" not in block
            assert "grid" not in block
            assert block["present"] is True

        # The plain-name block keeps everything else it carried.
        assert summed.uns["mobility_axis"]["n_channels"] == 4

    def test_a_none_key_and_a_table_without_uns_are_skipped(self):
        tables = {"stub_z0": SimpleNamespace(), "stub_z0_mobility": self._table("x")}

        _assembler().unname_declined_siblings(
            tables, (None, "stub_z0", "absent", "stub_z0_mobility"), ["ion_mobility"]
        )

        block = tables["stub_z0_mobility"].uns["mobility_axis"]
        assert "resolved_table" not in block


class TestTheReaderBoundary:
    def test_a_reader_that_cannot_say_costs_the_whole_block_and_is_logged(
        self, thyra_logs
    ):
        class _Refusing(_StubReader):
            def get_comprehensive_metadata(self):
                raise RuntimeError("the vendor library said no")

        with thyra_logs("thyra.converters.spatialdata", logging.WARNING) as records:
            uns = _assembler(_Refusing()).build(_context())

        assert uns == {}
        assert any(
            "Could not read metadata for uns provenance" in r.getMessage()
            for r in records
        )


@pytest.mark.parametrize("section", ["acquisition_params", "instrument_info"])
def test_a_section_the_reader_left_empty_is_omitted(section):
    comprehensive = _comprehensive()
    setattr(comprehensive, section, {})

    uns = _assembler(_StubReader(comprehensive=comprehensive)).build(_context())

    assert section not in uns
