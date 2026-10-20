"""What ``convert_msi`` does with arguments nobody meant to pass.

The CLI validates its own options with click, so these gaps were only ever
reachable through the Python API -- which is how Ousia and every notebook
call Thyra. Three shapes were found in the 2026-09-08 sweep (issue #261):

- a no-op argument accepted in silence when misspelled (``streaming``);
- an argument accepted and then used, turning every file into a refusal
  that named the caller's nonsense value back at them
  (``max_mass_axis_length``);
- a removed argument whose refusal reaches the caller as ``False`` rather
  than as an exception (``sparse_format``), which design decision D10
  described the other way round.

All three end the same way at the front door, which is the point: one
``ERROR`` line saying what to do, and ``False``. That is what every other
refusal does since issue #234, and what the CLI turns into exit 1.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest

from thyra.convert import _validate_streaming, convert_msi, sample_value_problem
from thyra.core.mass_axis import (
    validate_max_mass_axis_length as _validate_max_mass_axis_length,
)
from thyra.errors import ConversionRefused


class TestStreaming:
    """``streaming`` selects nothing, and still has to say when it is wrong.

    Every conversion streams since v3.23 (design decision D11), so the
    argument is kept only so existing calls keep running. That is not a
    licence to accept anything: a caller who writes ``streaming="flase"``
    believes they changed something.
    """

    @pytest.mark.parametrize("value", [True, False, "auto"])
    def test_the_three_real_values_are_accepted(self, value):
        assert _validate_streaming(value) is True

    @pytest.mark.parametrize("value", ["yes", "true", "false", None, 1, 2.0, []])
    def test_anything_else_is_refused(self, value, thyra_logs):
        with thyra_logs("thyra.convert", logging.ERROR) as records:
            assert _validate_streaming(value) is False
        said = " ".join(r.getMessage() for r in records)
        assert "streaming must be True, False or 'auto'" in said

    def test_zero_is_refused_rather_than_read_as_false(self, thyra_logs):
        """The sharpest case, and the reason this check exists.

        ``0 == False`` but ``0 is not False``, so ``streaming=0`` -- which a
        caller writes meaning ``False`` -- slipped past the ``is False``
        test that exists precisely to answer that request, and converted in
        silence.
        """
        with thyra_logs("thyra.convert", logging.ERROR) as records:
            assert _validate_streaming(0) is False
        assert "got 0" in " ".join(r.getMessage() for r in records)

    def test_it_is_refused_before_the_source_is_opened(self, tmp_path, thyra_logs):
        """A path that does not exist must not be the complaint.

        The argument check runs in ``_validate_input_parameters``, ahead of
        every path check and every read, so a typo fails while the caller is
        still looking at their own arguments rather than after a metadata
        scan.
        """
        with thyra_logs("thyra.convert", logging.ERROR) as records:
            assert (
                convert_msi(
                    tmp_path / "no-such-file.imzML",
                    tmp_path / "out.zarr",
                    streaming="yes",
                )
                is False
            )
        said = " ".join(r.getMessage() for r in records)
        assert "streaming must be" in said
        assert "does not exist" not in said


class TestMaxMassAxisLength:
    """The raw-axis cap is a count, so it is a positive integer or nothing."""

    @pytest.mark.parametrize("value", [None, 1, 10_000_000])
    def test_a_count_or_no_limit_is_accepted(self, value):
        assert _validate_max_mass_axis_length(value) == value

    @pytest.mark.parametrize("value", [-1, 0, 2.5, "ten", [10]])
    def test_anything_else_is_refused(self, value):
        with pytest.raises(ConversionRefused, match="positive integer or None"):
            _validate_max_mass_axis_length(value)

    @pytest.mark.parametrize("value", [True, False])
    def test_a_bool_is_refused_rather_than_read_as_a_count(self, value):
        """``bool`` is an ``int`` subclass; ``True`` would mean a cap of one."""
        with pytest.raises(ConversionRefused, match="positive integer or None"):
            _validate_max_mass_axis_length(value)

    def test_the_message_does_not_blame_the_file(self):
        """It used to.

        A cap of ``-1`` was accepted and then compared against a growing
        axis, so every file was refused for "exceeding -1 unique m/z
        values" -- a sentence about the data, for a mistake in the call.
        """
        with pytest.raises(ConversionRefused) as excinfo:
            _validate_max_mass_axis_length(-1)
        assert "exceeded" not in str(excinfo.value)
        assert "max_mass_axis_length must be" in str(excinfo.value)


class TestSparseFormatReachesTheCaller:
    """``sparse_format`` was removed in v3.22, and the removal still speaks.

    Design decision D10 said passing the keyword "raises". The converter
    class does raise (pinned in
    ``tests/unit/converters/test_spatialdata_converter.py``), but nothing
    called ``convert_msi`` and looked: its ``except ConversionRefused``
    turns every refusal into one logged ``ERROR`` and ``False``, which is
    the whole point of that handler (issue #234).

    So the front door reports this the same way it reports a Waters raster
    whose stage never moved, and D10's sentence was the thing that was
    wrong, not the code. Pinned here so the two cannot drift apart again:
    what matters is that the removal is *not silent*, and that the message
    is the one the caller needs.
    """

    FIXTURE = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "fixtures"
        / "iontof_sparse.imzML"
    )

    def test_it_is_reported_and_returns_false(self, tmp_path, thyra_logs):
        with thyra_logs("thyra.convert", logging.ERROR) as records:
            result = convert_msi(
                self.FIXTURE,
                tmp_path / "out.zarr",
                dataset_id="pinned",
                sparse_format="csr",
            )

        assert result is False
        said = " ".join(r.getMessage() for r in records)
        assert "sparse_format was removed" in said
        assert "tocsr()" in said, "the message must say what to call instead"

    def test_nothing_is_written(self, tmp_path):
        out = tmp_path / "out.zarr"
        assert (
            convert_msi(self.FIXTURE, out, dataset_id="pinned", sparse_format="csc")
            is False
        )
        assert not out.exists(), "a refused conversion must leave no store"


class TestSampleValues:
    """``sample_id`` and ``sample_name`` are the user's labels (issue #447).

    They are written into the store's metadata, so a bad one is refused
    with the other arguments, before the source is looked at.
    """

    @pytest.mark.parametrize("argument", ["sample_id", "sample_name"])
    @pytest.mark.parametrize(
        "value, said",
        [
            ("", "empty or blank"),
            ("   ", "empty or blank"),
            (7, "must be a string"),
            ("a\nb", "control characters"),
            ("a\x00b", "control characters"),
            ("x" * 257, "at most 256"),
        ],
    )
    def test_a_bad_value_is_refused_before_the_source_is_opened(
        self, tmp_path, thyra_logs, argument, value, said
    ):
        with thyra_logs("thyra.convert", logging.ERROR) as records:
            assert (
                convert_msi(
                    tmp_path / "no-such-file.imzML",
                    tmp_path / "out.zarr",
                    **{argument: value},
                )
                is False
            )
        message = " ".join(r.getMessage() for r in records)
        assert argument in message
        assert said in message
        assert "does not exist" not in message

    @pytest.mark.parametrize(
        "value", [None, "S1", "liver section 3", "x" * 256, "caf\u00e9 \u00b5m"]
    )
    def test_an_ordinary_label_passes(self, value):
        assert sample_value_problem("sample_id", value) is None


_FIXTURE = (
    Path(__file__).resolve().parents[1] / "data" / "fixtures" / "iontof_sparse.imzML"
)


@pytest.fixture
def no_reader(monkeypatch):
    """Fail the test if convert_msi gets as far as building a reader."""
    import thyra.convert

    def refuse(*args, **kwargs):
        raise AssertionError("the reader was built")

    monkeypatch.setattr(thyra.convert, "_create_reader", refuse)


def _errors(records):
    return [r.getMessage() for r in records if r.levelno >= logging.ERROR]


class TestKeywordsConvertMsiDoesNotUse:
    """Q6 of the 4.3.3 stress report, issue #519.

    A keyword the converter did not use was ignored but recorded in the
    store as applied; a reader option with a typo changed nothing. Both
    are refused now, in one line, before the source is opened.
    """

    @pytest.mark.parametrize(
        "kwargs, named, said",
        [
            ({"spectrum_type": "centroid"}, ["spectrum_type"], "reader_options"),
            (
                {"spectrum_type": "centroid", "max_mass_axis_length": 50},
                ["spectrum_type", "max_mass_axis_length"],
                "reader_options",
            ),
            ({"log_file": Path("run.log")}, ["log_file"], "takes no keyword"),
            (
                {"pixel_size_y_um": 99.0, "z_spacing_source": "measured"},
                ["pixel_size_y_um", "z_spacing_source"],
                "takes no keyword",
            ),
            (
                {"reader_options": {"intensity_treshold": 150.0}},
                ["intensity_treshold"],
                "no reader takes",
            ),
        ],
    )
    def test_refused_in_one_line_before_the_reader_is_built(
        self, tmp_path, thyra_logs, no_reader, kwargs, named, said
    ):
        out = tmp_path / "out.zarr"
        with thyra_logs("thyra", logging.ERROR) as records:
            assert convert_msi(_FIXTURE, out, pixel_size_um=10.0, **kwargs) is False
        errors = _errors(records)
        assert len(errors) == 1, errors
        assert all(name in errors[0] for name in named)
        assert said in errors[0]
        assert not out.exists()

    def test_the_documented_keywords_pass_the_check(self):
        from thyra.convert import CONVERTER_KEYWORDS, keyword_problems
        from thyra.reader_options import READER_OPTION_FORMATS

        assert (
            keyword_problems(
                {name: None for name in CONVERTER_KEYWORDS},
                {key: None for key in READER_OPTION_FORMATS},
            )
            == []
        )

    def test_an_option_the_format_ignores_is_named_once(self, tmp_path, thyra_logs):
        with thyra_logs("thyra", logging.WARNING) as records:
            assert convert_msi(
                _FIXTURE,
                tmp_path / "out.zarr",
                pixel_size_um=10.0,
                reader_options={"tdf_spectrum": "vendor_centroid"},
            )
        ignored = [r.getMessage() for r in records if "is ignored" in r.getMessage()]
        assert len(ignored) == 1
        assert "tdf_spectrum" in ignored[0] and "imzML" in ignored[0]


class TestTheListsAreTheDocstring:
    """The allowlists are the names ``convert_msi`` documents (#519)."""

    @staticmethod
    def _section(start: str, end: str) -> str:
        doc = convert_msi.__doc__ or ""
        return doc[doc.index(start) : doc.index(end)]

    def test_every_reader_key_is_documented_and_nothing_else(self):
        import re

        from thyra.reader_options import READER_OPTION_FORMATS

        section = self._section("reader_options:", "**kwargs:")
        documented = re.findall(r"^\s+- (\w+):", section, flags=re.MULTILINE)
        assert documented == list(READER_OPTION_FORMATS)

    def test_every_converter_keyword_is_documented(self):
        from thyra.convert import CONVERTER_KEYWORDS

        section = self._section("**kwargs:", "include_optical:")
        for name in CONVERTER_KEYWORDS:
            assert f"``{name}``" in section, name
        assert "max_mass_axis_length" not in section
        assert "spectrum_type" not in section


class TestNumbersAndFlags:
    """Numpy scalars are numbers; a bool is not; handle_3d is a bool (#519)."""

    @pytest.mark.parametrize(
        "value", [25, 25.0, np.float64(25), np.float32(25), np.int64(25)]
    )
    def test_any_real_number_is_a_pixel_size(self, value):
        from thyra.convert import _validate_numeric_parameters

        assert _validate_numeric_parameters(value, value) is True

    @pytest.mark.parametrize("value", [True, np.bool_(True), "25"])
    def test_a_bool_or_a_string_is_not(self, value, thyra_logs):
        from thyra.convert import _validate_numeric_parameters

        with thyra_logs("thyra.convert", logging.ERROR) as records:
            assert _validate_numeric_parameters(None, value) is False
        assert "Z spacing must be a number" in records[0].getMessage()

    @pytest.mark.parametrize("value", ["no", 1, 0, None])
    def test_handle_3d_takes_only_a_bool(self, tmp_path, thyra_logs, value):
        with thyra_logs("thyra.convert", logging.ERROR) as records:
            assert convert_msi(_FIXTURE, tmp_path / "o.zarr", handle_3d=value) is False
        assert "handle_3d must be True or False" in records[0].getMessage()

    def test_a_numpy_z_spacing_reaches_the_converter(self, tmp_path):
        from unittest.mock import MagicMock

        from thyra.converters.spatialdata.converter import SpatialDataConverter

        converter = SpatialDataConverter(
            MagicMock(),
            tmp_path / "o.zarr",
            pixel_size_um=10.0,
            handle_3d=True,
            z_spacing_um=np.float32(5),
        )
        assert converter._z_spacing_um_arg == 5.0

    def test_a_bool_z_spacing_is_refused_by_the_converter(self, tmp_path):
        from unittest.mock import MagicMock

        from thyra.converters.spatialdata.converter import SpatialDataConverter

        with pytest.raises(ConversionRefused, match="z_spacing_um"):
            SpatialDataConverter(
                MagicMock(), tmp_path / "o.zarr", pixel_size_um=10.0, z_spacing_um=True
            )

    def test_an_int_pixel_size_stores_float_columns(self, tmp_path):
        import zarr

        out = tmp_path / "o.zarr"
        assert convert_msi(
            _FIXTURE,
            out,
            pixel_size_um=20,
            handle_3d=np.bool_(False),
            resampling_config={"method": "nearest_neighbor", "target_bins": 50},
        )
        root = zarr.open_group(str(out), mode="r")
        obs = root["tables/msi_dataset_z0/obs"]
        assert obs["spatial_x"].dtype == np.float64
        assert obs["spatial_y"].dtype == np.float64
        assert root.attrs["conversion_options"]["handle_3d"] is False


class TestOutputPath:
    """The CLI's output-path checks reach the API (#519)."""

    def test_a_file_in_the_way_is_refused_before_the_reader_is_built(
        self, tmp_path, thyra_logs, no_reader
    ):
        blocker = tmp_path / "afile.txt"
        blocker.write_text("x")
        with thyra_logs("thyra", logging.ERROR) as records:
            assert (
                convert_msi(_FIXTURE, blocker / "sub" / "o.zarr", pixel_size_um=10.0)
                is False
            )
        errors = _errors(records)
        assert len(errors) == 1
        assert "is a file, not a directory" in errors[0]


class TestResamplingValues:
    """The CLI's value checks, for a dict and a ResamplingConfig alike (#519)."""

    @pytest.mark.parametrize(
        "config, said",
        [
            ({"width_at_mz": -0.5}, "'width_at_mz'"),
            ({"width_at_mz": 0.0}, "'width_at_mz'"),
            ({"min_mz": -50.0}, "'min_mz'"),
            ({"min_mz": 500.0, "max_mz": 400.0}, "below"),
            ({"gap_tolerance_da": 0.0}, "'gap_tolerance_da'"),
            ({"gap_tolerance_da": -1.0}, "'gap_tolerance_da'"),
            ({"reference_mz": -1.0}, "'reference_mz'"),
            ({"target_bins": "400"}, "whole number"),
            ({"target_bins": 400.5}, "whole number"),
            ({"target_bins": 1}, "at least 2"),
            ({"target_bins": -400}, "at least 2"),
            ({"target_bins": True}, "whole number"),
            ({"tof_a": -1.0, "tof_b": 1.0}, "'tof_a'"),
            ({"tof_a": 0.0, "tof_b": 0.0}, "no width"),
        ],
    )
    @pytest.mark.parametrize("as_dataclass", [False, True])
    def test_refused(self, config, said, as_dataclass):
        from thyra.resampling.axis_planner import normalize_resampling_config
        from thyra.resampling.types import ResamplingConfig

        if as_dataclass:
            fields = dict(config)
            if "width_at_mz" in fields:
                fields["mass_width_da"] = fields.pop("width_at_mz")
            config = ResamplingConfig(**fields)
        with pytest.raises(ConversionRefused, match=said):
            normalize_resampling_config(config)

    def test_a_string_method_in_the_dataclass_is_resolved(self):
        from thyra.resampling.axis_planner import normalize_resampling_config
        from thyra.resampling.types import ResamplingConfig, ResamplingMethod

        config = normalize_resampling_config(
            ResamplingConfig(method="tic_preserving", target_bins=np.int64(400))
        )
        assert config.method is ResamplingMethod.TIC_PRESERVING
        assert type(config.target_bins) is int

    def test_refused_before_the_reader_is_built(self, tmp_path, thyra_logs, no_reader):
        with thyra_logs("thyra", logging.ERROR) as records:
            assert (
                convert_msi(
                    _FIXTURE,
                    tmp_path / "o.zarr",
                    pixel_size_um=10.0,
                    resampling_config={"width_at_mz": -0.5},
                )
                is False
            )
        assert len(_errors(records)) == 1

    def test_a_gap_tolerance_under_nearest_neighbor_is_named_as_ignored(
        self, tmp_path, thyra_logs
    ):
        with thyra_logs("thyra", logging.INFO) as records:
            assert convert_msi(
                _FIXTURE,
                tmp_path / "o.zarr",
                pixel_size_um=10.0,
                resampling_config={
                    "method": "nearest_neighbor",
                    "target_bins": 50,
                    "gap_tolerance_da": 0.01,
                },
            )
        said = [r.getMessage() for r in records]
        assert any("gap-tolerance" in m and "is ignored" in m for m in said)
        assert not any("are zeroed" in m for m in said)
