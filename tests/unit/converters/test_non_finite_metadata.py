"""A NaN or an infinity in the vendor metadata costs one field, not the block.

The models refuse such a number (#503, D62). The trap in that rule: the
uns assembler catches the model's error, so a value that reached a model
used to write the store with no ``msi_metadata`` block at all, while the
conversion reported success. Every test here converts, or builds the
block the converter writes, and asserts the block is still there.
"""

import dataclasses
import json
import logging
import math
from pathlib import Path
from types import SimpleNamespace

import pytest

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.convert import _determine_pixel_size
from thyra.converters.spatialdata import SpatialDataConverter
from thyra.converters.spatialdata.uns_assembler import UnsAssembler, UnsContext
from thyra.core.base_converter import PixelSizeSource
from thyra.errors import ConversionRefused
from thyra.metadata.schema import MSI_METADATA_UNS_KEY
from thyra.metadata.schema.store_io import read_msi_metadata_blocks

NON_FINITE = [float("nan"), float("inf"), float("-inf")]
NON_FINITE_IDS = ["nan", "inf", "-inf"]


def _config() -> MockMSIConfig:
    return MockMSIConfig(
        n_x=3, n_y=2, n_z=1, n_mz_bins=64, peaks_per_spectrum=(3, 5), seed=5
    )


class _NonFiniteMobilityReader(MockMSIReader):
    """A source whose mobility range is NaN, as a Bruker text field can be."""

    def get_comprehensive_metadata(self):
        comprehensive = super().get_comprehensive_metadata()
        format_specific = dict(comprehensive.format_specific or {})
        format_specific["ion_mobility"] = {
            "present": True,
            "one_over_k0_range": [float("nan"), 1.99],
            "num_scans_max": float("inf"),
        }
        return dataclasses.replace(comprehensive, format_specific=format_specific)


class _NonFinitePitchReader(MockMSIReader):
    """A source that states a pitch of NaN by infinity."""

    def get_essential_metadata(self):
        essential = super().get_essential_metadata()
        object.__setattr__(essential, "pixel_size", (float("nan"), float("inf")))
        return essential


def _convert(reader, output: Path, **kwargs) -> dict:
    converter = SpatialDataConverter(
        reader=reader,
        output_path=output,
        dataset_id="m",
        include_optical=False,
        **kwargs,
    )
    assert converter.convert() is True
    blocks = read_msi_metadata_blocks(output)
    assert set(blocks) == {"m_z0"}, "the store lost its msi_metadata block"
    return blocks["m_z0"]


class TestTheBlockIsStillWritten:
    def test_a_non_finite_mobility_range_costs_the_range_only(self, tmp_path):
        block = _convert(
            _NonFiniteMobilityReader(_config()),
            tmp_path / "out.zarr",
            pixel_size_um=10.0,
        )

        mobility = block["ms_analysis"]["ion_mobility"]
        assert mobility["present"] is True
        assert "range_lower" not in mobility and "num_scans" not in mobility

    def test_a_non_finite_source_pitch_leaves_the_default(self, tmp_path, thyra_logs):
        """No pitch given and none usable: the store says it used the default."""
        with thyra_logs("thyra.core.base_converter", logging.WARNING) as logs:
            block = _convert(_NonFinitePitchReader(_config()), tmp_path / "out.zarr")

        assert block["ms_analysis"]["pixel_size_um"] == {"x": 1.0, "y": 1.0}
        assert block["provenance"]["pixel_size_source"] == "default"
        assert any("not a finite positive" in r.getMessage() for r in logs)


class TestAPitchIsRefused:
    """#503 path 5: ``<= 0`` alone let NaN and infinity through."""

    @pytest.mark.parametrize("value", NON_FINITE, ids=NON_FINITE_IDS)
    def test_by_the_converter(self, value, tmp_path):
        with pytest.raises(ConversionRefused, match="finite positive"):
            SpatialDataConverter(
                reader=MockMSIReader(_config()),
                output_path=tmp_path / "out.zarr",
                dataset_id="m",
                pixel_size_um=value,
            )

    @pytest.mark.parametrize("value", NON_FINITE, ids=NON_FINITE_IDS)
    def test_a_detected_y_by_the_converter(self, value, tmp_path):
        with pytest.raises(ConversionRefused, match="detected_y_um"):
            SpatialDataConverter(
                reader=MockMSIReader(_config()),
                output_path=tmp_path / "out.zarr",
                dataset_id="m",
                pixel_size_um=10.0,
                pixel_size_source=PixelSizeSource.AUTO_DETECTED,
                pixel_size_detection_info={"detected_y_um": value},
            )

    @pytest.mark.parametrize(
        "pitch",
        [(float("nan"), 10.0), (10.0, float("inf")), (0.0, 10.0)],
        ids=["nan-x", "inf-y", "zero"],
    )
    def test_when_detected_by_convert(self, pitch):
        """Refused like a missing pitch, so the user can state the right one."""
        reader = _NonFinitePitchReader(_config())
        reader.get_essential_metadata = lambda: SimpleNamespace(  # type: ignore
            pixel_size=pitch
        )

        with pytest.raises(ConversionRefused, match="--pixel-size"):
            _determine_pixel_size(reader, None, "mock")


class TestTheProcessingRecord:
    """#503 path 6: a value a caller slips past the config check is left out."""

    def test_a_non_finite_config_value_leaves_the_block(self):
        config = SimpleNamespace(
            method="tic_preserving", gap_tolerance_da=float("nan"), min_mz=100.0
        )
        assembler = UnsAssembler(
            MockMSIReader(_config()),
            pixel_size_xy=lambda: (10.0, 10.0),
            pixel_size_detection_info=None,
            resampling_config=config,
        )
        context = UnsContext(
            mobility_table_key=None,
            msms_table_key=None,
            mobility_grid=None,
            resolved_resampling_plan={"tof_a": float("inf")},
            region_info=None,
            mass_axis_fingerprint=None,
            pixel_size_source=PixelSizeSource.USER_PROVIDED,
            mobility_heatmap=lambda: None,
        )

        block = assembler.build(context)[MSI_METADATA_UNS_KEY]

        steps = json.loads(block["processing"])
        resampling = steps[-1]["parameters"]
        assert resampling == {"method": "tic_preserving", "min_mz": 100.0}
        assert all(
            math.isfinite(v) for v in resampling.values() if isinstance(v, float)
        )
