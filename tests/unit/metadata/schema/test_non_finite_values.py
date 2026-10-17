"""A NaN or an infinity never reaches the ``msi_metadata`` block (#503, D62).

JSON has no NaN or infinity, so the published schema never allowed them,
but the models took them and a conversion stored them. The builder now
leaves each such vendor value out, the models refuse one that slips past,
and ``thyra validate`` names one an older store still holds.
"""

import json
import logging
from types import SimpleNamespace

import numpy as np
import pytest
from pydantic import ValidationError

from thyra.core.mobility import MobilityAxis
from thyra.metadata.document import _pixel_size
from thyra.metadata.schema import build_msi_metadata
from thyra.metadata.schema.builder import _mobility_grid, finite_parameters
from thyra.metadata.schema.models import (
    MSI_METADATA_SCHEMA_VERSION,
    IonMobility,
    IsolationWindow,
    MobilityGrid,
    MSIMetadata,
    PixelSizeUm,
    ProcessingStep,
    SoftwareRef,
)
from thyra.metadata.schema.validate import validate_document
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata

NON_FINITE = [float("nan"), float("inf"), float("-inf")]
NON_FINITE_IDS = ["nan", "inf", "-inf"]


def _comprehensive(ion_mobility):
    essential = EssentialMetadata(
        dimensions=(3, 2, 1),
        coordinate_bounds=(0.0, 2.0, 0.0, 1.0),
        mass_range=(50.0, 1000.0),
        pixel_size=(20.0, 20.0),
        n_spectra=6,
        total_peaks=492,
        source_path="synthetic_tims.d",
    )
    return ComprehensiveMetadata(
        essential=essential,
        format_specific={"ion_mobility": ion_mobility},
        acquisition_params={},
        instrument_info={},
        raw_metadata={},
    )


def _mobility_block(ion_mobility):
    meta = build_msi_metadata(_comprehensive(ion_mobility), pixel_size_um=(20.0, 20.0))
    return meta.to_uns_dict()["ms_analysis"]["ion_mobility"]


class TestTheModelsRefuseThem:
    """The rule a value that slips past the builder meets (#503 path 7)."""

    @pytest.mark.parametrize("value", NON_FINITE, ids=NON_FINITE_IDS)
    def test_an_unbounded_field(self, value):
        with pytest.raises(ValidationError, match="finite number"):
            IonMobility(present=True, range_lower=value, range_upper=1.5)

    def test_a_bounded_field_refuses_infinity(self):
        """``gt=0`` refused NaN already; +inf passed it."""
        with pytest.raises(ValidationError, match="finite number"):
            PixelSizeUm(x=float("inf"), y=10.0)

    @pytest.mark.parametrize("value", NON_FINITE, ids=NON_FINITE_IDS)
    def test_an_isolation_window(self, value):
        with pytest.raises(ValidationError):
            IsolationWindow(target=300.0, collision_energy=value)

    @pytest.mark.parametrize("value", NON_FINITE, ids=NON_FINITE_IDS)
    def test_a_processing_parameter(self, value):
        """The union of types gives no float branch a pass either."""
        with pytest.raises(ValidationError):
            ProcessingStep(
                name="conversion",
                software=SoftwareRef(name="thyra", version="0"),
                parameters={"intensity_threshold": value},
            )

    def test_the_root_model_inherits_the_rule(self):
        """``MSIMetadata`` sets its own config; pydantic merges the base's."""
        assert MSIMetadata.model_config.get("allow_inf_nan") is False

    def test_the_packed_strings_are_strict_json(self):
        """A value past validation raises rather than writing ``NaN``."""
        step = ProcessingStep.model_construct(
            name="conversion",
            software=SoftwareRef(name="thyra", version="0"),
            parameters={"intensity_threshold": float("nan")},
        )
        meta = build_msi_metadata(None, pixel_size_um=(10.0, 10.0))
        meta = meta.model_copy(update={"processing": [step]})
        with pytest.raises(ValueError):
            meta.to_uns_dict()


class TestTheMobilityRange:
    """#503 path 1: the Bruker GlobalMetadata range, the Waters drift times."""

    @pytest.mark.parametrize("value", NON_FINITE, ids=NON_FINITE_IDS)
    @pytest.mark.parametrize("edge", [0, 1], ids=["lower", "upper"])
    def test_a_non_finite_edge_leaves_the_range_out(self, value, edge):
        mobility_range = [0.9, 1.99]
        mobility_range[edge] = value

        block = _mobility_block({"present": True, "one_over_k0_range": mobility_range})

        assert block["present"] is True
        assert "range_lower" not in block and "range_upper" not in block

    def test_a_finite_range_is_kept(self):
        block = _mobility_block({"present": True, "one_over_k0_range": [0.9, 1.99]})

        assert (block["range_lower"], block["range_upper"]) == (0.9, 1.99)

    @pytest.mark.parametrize("value", NON_FINITE, ids=NON_FINITE_IDS)
    def test_a_non_finite_scan_count_is_left_out(self, value):
        """``int(inf)`` raised OverflowError and stopped the conversion."""
        block = _mobility_block({"present": True, "num_scans_max": value})

        assert "num_scans" not in block

    def test_the_waters_report_spans_the_finite_drift_times(self):
        """``np.min`` passed a NaN drift time into the reported range."""
        axis = MobilityAxis(
            kind_accession="MS:1002476",
            kind_name="ion mobility drift time",
            unit_accession="UO:0000028",
            unit_name="millisecond",
            values=np.array([np.nan, 0.5, 7.25, np.inf]),
        )

        assert axis.to_extractor_report()["range"] == [0.5, 7.25]

    def test_drift_times_with_no_finite_value_report_no_range(self):
        axis = MobilityAxis(
            kind_accession="MS:1002476",
            kind_name="ion mobility drift time",
            unit_accession="UO:0000028",
            unit_name="millisecond",
            values=np.array([np.nan, np.nan]),
        )

        assert "range" not in axis.to_extractor_report()

    @pytest.mark.parametrize("value", NON_FINITE, ids=NON_FINITE_IDS)
    def test_a_grid_with_a_non_finite_edge_is_dropped(self, value):
        assert (
            _mobility_grid(
                {"law": "linear", "lower": value, "upper": 1.5, "n_channels": 10}
            )
            is None
        )
        assert isinstance(
            _mobility_grid(
                {"law": "linear", "lower": 0.5, "upper": 1.5, "n_channels": 10}
            ),
            MobilityGrid,
        )


class TestProcessingParameters:
    """#503 path 6: a parameter that is not finite is left out of its step."""

    @pytest.mark.parametrize(
        "value",
        [*NON_FINITE, np.float32("nan"), np.float64("inf")],
        ids=[*NON_FINITE_IDS, "np32-nan", "np64-inf"],
    )
    def test_a_non_finite_number_is_left_out(self, value):
        assert finite_parameters({"gap_tolerance_da": value, "method": "x"}) == {
            "method": "x"
        }

    def test_everything_else_is_kept_as_given(self):
        given = {"a": 1, "b": 2.5, "c": "tof", "d": True, "e": np.float32(1.5)}

        assert finite_parameters(given) == given


class TestThePixelSizeOfADocument:
    """#503 path 5 on the ``thyra metadata`` route: no pitch, not no document."""

    @pytest.mark.parametrize(
        "pitch",
        [(float("nan"), 10.0), (10.0, float("inf")), (0.0, 10.0)],
        ids=["nan-x", "inf-y", "zero"],
    )
    def test_an_unusable_pitch_states_none(self, pitch, thyra_logs):
        with thyra_logs("thyra.metadata.document", logging.WARNING) as logs:
            assert _pixel_size(SimpleNamespace(pixel_size=pitch)) is None
        assert len(logs) == 1

    def test_a_usable_pitch_is_kept(self):
        assert _pixel_size(SimpleNamespace(pixel_size=(20, 25.0))) == (20.0, 25.0)


def _stored_block_with(**ms_analysis):
    """A block as an older Thyra stored it, NaN included."""
    return {
        "schema_version": MSI_METADATA_SCHEMA_VERSION,
        "ms_analysis": {"pixel_size_um": {"x": 10.0, "y": 10.0}, **ms_analysis},
        "processing": json.dumps(
            [
                {
                    "name": "conversion",
                    "software": {"name": "thyra", "version": "0"},
                    "parameters": {"intensity_threshold": float("nan")},
                }
            ]
        ),
        "provenance": {"thyra_version": "0"},
    }


class TestValidateOnAnOlderStore:
    """The verdict on a store written before the rule (D62): an error."""

    def test_each_value_is_one_error_that_says_what_to_do(self):
        block = _stored_block_with(
            ion_mobility={
                "present": True,
                "range_lower": float("inf"),
                "range_upper": 1.5,
            }
        )
        assert "NaN" in block["processing"]  # how json.dumps wrote it then

        meta, issues = validate_document(block)

        assert meta is None
        assert [(i.severity, i.location) for i in issues] == [
            ("error", "ms_analysis.ion_mobility.range_lower"),
            ("error", "processing.0.parameters.intensity_threshold"),
        ]
        assert all("converting the source again" in i.message for i in issues)

    def test_a_finite_block_still_validates(self):
        block = _stored_block_with()
        block["processing"] = block["processing"].replace("NaN", "50.0")

        meta, issues = validate_document(block)

        assert meta is not None and not issues
        assert meta.processing[0].parameters["intensity_threshold"] == 50.0

    def test_an_unknown_key_holding_nan_is_still_unknown(self):
        block = _stored_block_with(bogus=float("nan"))

        _, issues = validate_document(block)

        assert any(
            i.location == "ms_analysis.bogus" and "not permitted" in i.message
            for i in issues
        )


def test_an_infinite_ms_level_leaves_the_fragmentation_block_out():
    """``int(inf)`` raised OverflowError past the builder's own catch."""
    meta = build_msi_metadata(
        None, pixel_size_um=(1.0, 1.0), fragmentation={"ms_level": float("inf")}
    )

    assert meta.ms_analysis.fragmentation is None
