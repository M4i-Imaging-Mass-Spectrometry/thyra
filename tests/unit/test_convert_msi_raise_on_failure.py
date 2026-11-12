"""``convert_msi(..., raise_on_failure=True)`` raises the reason it failed with (D87).

``convert_msi`` returns ``bool``, and a program that uses Thyra as a library
could not tell its user why a ``False`` came back: the reason went to the
log. The opt-in parameter raises it instead. The default is pinned here as
firmly as the new behaviour, because 50 tests and an unknown number of
callers assert ``is True`` / ``is False``.

What these tests hold in place:

* the default returns the same value and logs the same ERROR lines;
* a refusal raises ``ConversionRefused`` and anything else raises
  ``ConversionFailed``, both built from the text of the log line;
* the exception is new, with no cause and no context (issue #249: a
  retained exception pins the reader and its scratch memmaps);
* the raise comes after the cleanup, so the partial store is already set
  aside when the caller sees it;
* nothing recorded for one call leaks into the next.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest

from thyra import convert as convert_module
from thyra.convert import CONVERTER_KEYWORDS, convert_msi
from thyra.converters.spatialdata.base_spatialdata_converter import (
    BaseSpatialDataConverter,
)
from thyra.errors import ConversionFailed, ConversionRefused

_FIXTURES = Path(__file__).resolve().parent.parent / "data" / "fixtures"
_DECLARED = _FIXTURES / "pixel_size_unit_declared.imzML"
_NO_PIXEL_SIZE = _FIXTURES / "pixel_size_area_old_name_no_extent.imzML"


def _errors(records):
    return [r.getMessage() for r in records if r.levelno >= logging.ERROR]


def _convert(out: Path, **kwargs):
    return convert_msi(_DECLARED, out, resampling_config={}, **kwargs)


class TestTheExceptionTypes:
    def test_a_refusal_is_a_failure_and_a_value_error(self):
        refusal = ConversionRefused("no")
        assert isinstance(refusal, ConversionFailed)
        assert isinstance(refusal, ValueError)

    def test_a_failure_is_not_a_refusal_or_a_value_error(self):
        failure = ConversionFailed("no")
        assert not isinstance(failure, ConversionRefused)
        assert not isinstance(failure, ValueError)

    def test_both_are_public(self):
        import thyra

        assert thyra.ConversionFailed is ConversionFailed
        assert thyra.ConversionRefused is ConversionRefused


class TestTheDefaultIsUnchanged:
    def test_an_unknown_keyword_returns_false_and_logs_one_error(
        self, tmp_path, thyra_logs
    ):
        out = tmp_path / "out.zarr"
        with thyra_logs("thyra", logging.ERROR) as records:
            assert convert_msi(_DECLARED, out, pixel_size_um=10.0, bogus=1) is False
        assert _errors(records) == [
            "Nothing was converted. convert_msi takes no keyword 'bogus'; its "
            "converter keywords are " + ", ".join(CONVERTER_KEYWORDS) + "."
        ]

    def test_raising_logs_the_same_lines(self, tmp_path, thyra_logs):
        out = tmp_path / "out.zarr"
        with thyra_logs("thyra", logging.ERROR) as records:
            convert_msi(_DECLARED, out, pixel_size_um=10.0, bogus=1)
        by_default = _errors(records)
        with thyra_logs("thyra", logging.ERROR) as records:
            with pytest.raises(ConversionRefused):
                convert_msi(
                    _DECLARED, out, pixel_size_um=10.0, raise_on_failure=True, bogus=1
                )
        assert _errors(records) == by_default

    def test_the_flag_is_not_a_converter_keyword(self):
        assert "raise_on_failure" not in CONVERTER_KEYWORDS

    def test_the_flag_is_not_forwarded_to_the_converter(self, tmp_path, monkeypatch):
        seen = {}
        real = convert_module._create_converter

        def spy(*args, **kwargs):
            seen.update(kwargs)
            return real(*args, **kwargs)

        monkeypatch.setattr(convert_module, "_create_converter", spy)
        assert _convert(tmp_path / "out.zarr", raise_on_failure=True) is True
        assert "raise_on_failure" not in seen


class TestArgumentsAndPaths:
    def test_an_unknown_keyword_names_it(self, tmp_path):
        with pytest.raises(ConversionRefused, match="'bogus'"):
            _convert(tmp_path / "out.zarr", raise_on_failure=True, bogus=1)

    def test_a_bad_resampling_config_gives_the_planners_text(
        self, tmp_path, thyra_logs
    ):
        out = tmp_path / "out.zarr"
        bad = {"width_at_mz": -0.5}
        with thyra_logs("thyra", logging.ERROR) as records:
            assert convert_msi(_DECLARED, out, resampling_config=bad) is False
        planner_text = _errors(records)
        assert len(planner_text) == 1
        with pytest.raises(ConversionRefused) as raised:
            convert_msi(_DECLARED, out, resampling_config=bad, raise_on_failure=True)
        assert str(raised.value) == planner_text[0]

    def test_a_missing_input_is_refused(self, tmp_path):
        with pytest.raises(ConversionRefused, match="Input path does not exist"):
            convert_msi(
                tmp_path / "nothing.imzML",
                tmp_path / "out.zarr",
                raise_on_failure=True,
            )

    def test_an_existing_output_is_refused_and_left_alone(self, tmp_path):
        out = tmp_path / "out.zarr"
        out.mkdir()
        with pytest.raises(ConversionRefused, match="already exists"):
            _convert(out, raise_on_failure=True)
        assert out.is_dir()

    @pytest.mark.parametrize(
        "kwargs, said",
        [
            ({"dataset_id": "a b"}, "Dataset ID"),
            ({"format_type": ""}, "Format type"),
            ({"pixel_size_um": -1.0}, "Pixel size"),
            ({"handle_3d": "no"}, "handle_3d must be True or False"),
            ({"streaming": "yes"}, "streaming must be"),
            ({"sample_id": 5}, "sample_id"),
        ],
    )
    def test_a_bad_argument_is_refused(self, tmp_path, kwargs, said):
        with pytest.raises(ConversionRefused, match=said):
            convert_msi(
                _DECLARED, tmp_path / "out.zarr", raise_on_failure=True, **kwargs
            )

    @pytest.mark.parametrize("path_kwarg", ["input_path", "output_path"])
    def test_an_empty_path_is_refused(self, tmp_path, path_kwarg):
        kwargs = {"input_path": _DECLARED, "output_path": tmp_path / "out.zarr"}
        kwargs[path_kwarg] = ""
        with pytest.raises(ConversionRefused, match="must be a valid string"):
            convert_msi(raise_on_failure=True, **kwargs)

    def test_an_output_under_a_file_raises_the_logged_reason(
        self, tmp_path, thyra_logs
    ):
        blocker = tmp_path / "not_a_folder"
        blocker.write_text("x")
        out = blocker / "out.zarr"
        with thyra_logs("thyra", logging.ERROR) as records:
            assert _convert(out) is False
        logged = _errors(records)
        assert len(logged) == 1
        with pytest.raises(ConversionRefused) as raised:
            _convert(out, raise_on_failure=True)
        assert str(raised.value) == logged[0]

    def test_a_folder_of_several_runs_is_refused(self, tmp_path, monkeypatch):
        def two_runs(path):
            raise ConversionRefused("Two runs: a.d and b.d")

        monkeypatch.setattr(convert_module, "resolve_run", two_runs)
        with pytest.raises(ConversionRefused, match="Two runs"):
            _convert(tmp_path / "out.zarr", raise_on_failure=True)

    @pytest.mark.parametrize("value", ["yes", 1, None])
    def test_a_non_bool_flag_is_refused(self, tmp_path, value):
        out = tmp_path / "out.zarr"
        if value:
            with pytest.raises(ConversionRefused, match="raise_on_failure must be"):
                _convert(out, raise_on_failure=value)
        else:
            assert _convert(out, raise_on_failure=value) is False
        assert not out.exists()

    @pytest.mark.parametrize("argument", ["handle_3d", "streaming"])
    def test_a_value_whose_repr_raises_is_still_a_refusal(self, tmp_path, argument):
        class Unprintable:
            def __repr__(self):
                raise RuntimeError("no repr")

        out = tmp_path / "out.zarr"
        assert _convert(out, **{argument: Unprintable()}) is False
        with pytest.raises(ConversionRefused, match="Unprintable"):
            _convert(out, raise_on_failure=True, **{argument: Unprintable()})

    def test_a_numpy_bool_flag_is_accepted(self, tmp_path):
        assert _convert(tmp_path / "out.zarr", raise_on_failure=np.True_) is True


class TestARefusalFromTheReader:
    def test_a_source_with_no_pixel_size_raises_the_readers_text(
        self, tmp_path, thyra_logs
    ):
        out = tmp_path / "out.zarr"
        with thyra_logs("thyra", logging.ERROR) as records:
            assert convert_msi(_NO_PIXEL_SIZE, out) is False
        logged = _errors(records)
        assert len(logged) == 1 and "Pixel size not found" in logged[0]
        with pytest.raises(ConversionRefused) as raised:
            convert_msi(_NO_PIXEL_SIZE, out, raise_on_failure=True)
        assert str(raised.value) == logged[0]
        assert not out.exists()


class TestAnUnexpectedFailure:
    def _fail_with(self, monkeypatch, exc):
        def boom(converter):
            raise exc

        monkeypatch.setattr(convert_module, "_perform_conversion", boom)

    def test_it_raises_a_failure_that_is_not_a_refusal(self, tmp_path, monkeypatch):
        self._fail_with(monkeypatch, RuntimeError("boom"))
        with pytest.raises(ConversionFailed) as raised:
            _convert(tmp_path / "out.zarr", raise_on_failure=True)
        assert not isinstance(raised.value, ConversionRefused)
        assert not isinstance(raised.value, ValueError)
        assert "boom" in str(raised.value)

    def test_it_is_new_and_carries_no_original(self, tmp_path, monkeypatch):
        self._fail_with(monkeypatch, RuntimeError("boom"))
        with pytest.raises(ConversionFailed) as raised:
            _convert(tmp_path / "out.zarr", raise_on_failure=True)
        assert raised.value.__cause__ is None
        assert raised.value.__context__ is None
        assert raised.value.__traceback__ is not None
        # Nothing below the raising frame: no frame of the failed run.
        names = []
        tb = raised.value.__traceback__
        while tb is not None:
            names.append(tb.tb_frame.f_code.co_name)
            tb = tb.tb_next
        assert "boom" not in names

    def test_a_refusal_in_the_conversion_stays_a_refusal(self, tmp_path, monkeypatch):
        original = ConversionRefused("not this source")
        self._fail_with(monkeypatch, original)
        with pytest.raises(ConversionRefused, match="not this source") as raised:
            _convert(tmp_path / "out.zarr", raise_on_failure=True)
        # A new exception from the text, not the one the run raised: that
        # one's traceback holds the frames of the failed run.
        assert raised.value is not original
        assert raised.value.__cause__ is None
        assert raised.value.__context__ is None

    def test_an_interrupt_before_the_source_is_read(self, tmp_path, monkeypatch):
        def interrupted(*args, **kwargs):
            raise KeyboardInterrupt()

        monkeypatch.setattr(convert_module, "_create_reader", interrupted)
        assert _convert(tmp_path / "out.zarr") is False
        with pytest.raises(KeyboardInterrupt):
            _convert(tmp_path / "out2.zarr", raise_on_failure=True)

    def test_an_empty_message_names_the_type(self, tmp_path, monkeypatch):
        self._fail_with(monkeypatch, MemoryError())
        with pytest.raises(ConversionFailed, match="MemoryError"):
            _convert(tmp_path / "out.zarr", raise_on_failure=True)

    def test_an_interrupt_raises_a_new_interrupt(self, tmp_path, monkeypatch):
        interrupt = KeyboardInterrupt()
        self._fail_with(monkeypatch, interrupt)
        with pytest.raises(KeyboardInterrupt) as raised:
            _convert(tmp_path / "out.zarr", raise_on_failure=True)
        assert raised.value is not interrupt
        assert raised.value.__context__ is None

    def test_the_default_still_returns_false(self, tmp_path, monkeypatch):
        self._fail_with(monkeypatch, RuntimeError("boom"))
        assert _convert(tmp_path / "out.zarr") is False
        self._fail_with(monkeypatch, KeyboardInterrupt())
        assert _convert(tmp_path / "out2.zarr") is False

    def test_the_cleanup_runs_before_the_raise(self, tmp_path, monkeypatch):
        out = tmp_path / "out.zarr"

        def write_then_fail(converter):
            out.mkdir()
            (out / "zarr.json").write_text("{}")
            raise RuntimeError("boom")

        monkeypatch.setattr(convert_module, "_perform_conversion", write_then_fail)
        with pytest.raises(ConversionFailed):
            _convert(out, raise_on_failure=True)
        assert not out.exists()
        assert (tmp_path / "out.zarr.failed").is_dir()


class TestTheConverterSaysWhy:
    """Five reasons used to stay inside the converter; each reaches the caller."""

    def test_a_converter_that_returns_false_with_a_reason(self, tmp_path, monkeypatch):
        class Stub:
            failure_reason = "the stub could not finish"
            failure_kind = "error"
            store_complete = False

            def convert(self):
                return False

        monkeypatch.setattr(convert_module, "_create_converter", lambda *a, **k: Stub())
        with pytest.raises(ConversionFailed, match="the stub could not finish") as e:
            _convert(tmp_path / "out.zarr", raise_on_failure=True)
        assert not isinstance(e.value, ConversionRefused)

    def test_a_converter_that_returns_false_with_no_reason(self, tmp_path, monkeypatch):
        class Stub:
            store_complete = False

            def convert(self):
                return False

        monkeypatch.setattr(convert_module, "_create_converter", lambda *a, **k: Stub())
        with pytest.raises(ConversionFailed, match="the log has the details"):
            _convert(tmp_path / "out.zarr", raise_on_failure=True)

    def test_the_converters_refusal(self, tmp_path, monkeypatch):
        def refuse(self, state):
            raise ConversionRefused("this axis has no width law")

        monkeypatch.setattr(BaseSpatialDataConverter, "_save_output", refuse)
        with pytest.raises(ConversionRefused, match="no width law"):
            _convert(tmp_path / "out.zarr", raise_on_failure=True)

    def test_the_converters_interrupt(self, tmp_path, monkeypatch):
        def interrupted(self, state):
            raise KeyboardInterrupt

        monkeypatch.setattr(BaseSpatialDataConverter, "_save_output", interrupted)
        with pytest.raises(KeyboardInterrupt):
            _convert(tmp_path / "out.zarr", raise_on_failure=True)

    def test_the_converters_unexpected_error(self, tmp_path, monkeypatch):
        def broken(self, state):
            raise RuntimeError("converter broke")

        monkeypatch.setattr(BaseSpatialDataConverter, "_save_output", broken)
        with pytest.raises(ConversionFailed) as raised:
            _convert(tmp_path / "out.zarr", raise_on_failure=True)
        assert not isinstance(raised.value, ConversionRefused)
        assert str(raised.value) == "Error during conversion: converter broke"
        assert raised.value.__cause__ is None and raised.value.__context__ is None

    def test_a_failed_save_quarantines_the_partial_store_first(
        self, tmp_path, monkeypatch
    ):
        import spatialdata

        def write_part_then_fail(self, path, *args, **kwargs):
            Path(path).mkdir(parents=True)
            (Path(path) / "zarr.json").write_text("{}")
            raise OSError(28, "No space left on device")

        monkeypatch.setattr(spatialdata.SpatialData, "write", write_part_then_fail)
        out = tmp_path / "out.zarr"
        with pytest.raises(ConversionFailed) as raised:
            _convert(out, raise_on_failure=True)
        assert str(raised.value).startswith("Error saving SpatialData:")
        assert "No space left" in str(raised.value)
        assert not isinstance(raised.value, ConversionRefused)
        assert raised.value.__cause__ is None and raised.value.__context__ is None
        assert not out.exists()
        assert sorted(p.name for p in tmp_path.iterdir()) == ["out.zarr.failed"]

    def test_a_failed_rename_names_the_temporary_folder(self, tmp_path, monkeypatch):
        from thyra.converters.spatialdata import base_spatialdata_converter as module

        def refuse_rename(rename):
            raise PermissionError(13, "Access is denied")

        monkeypatch.setattr(module, "retry_on_os_error", refuse_rename)
        out = tmp_path / "out.zarr"
        with pytest.raises(ConversionFailed) as raised:
            _convert(out, raise_on_failure=True)
        assert "The store is complete, but it could not be renamed" in str(raised.value)
        assert "Access is denied" in str(raised.value)
        assert not out.exists()
        # The finished store is still there, under the name the message gives.
        kept = [p for p in tmp_path.iterdir() if ".incomplete-" in p.name]
        assert len(kept) == 1 and kept[0].name in str(raised.value)

    def test_the_default_still_returns_false_for_a_failed_save(
        self, tmp_path, monkeypatch
    ):
        import spatialdata

        def fail(self, *args, **kwargs):
            raise OSError(28, "No space left on device")

        monkeypatch.setattr(spatialdata.SpatialData, "write", fail)
        assert _convert(tmp_path / "out.zarr") is False


class TestNothingLeaksBetweenCalls:
    def test_a_failure_then_a_success(self, tmp_path):
        with pytest.raises(ConversionRefused):
            _convert(tmp_path / "a.zarr", raise_on_failure=True, bogus=1)
        assert _convert(tmp_path / "b.zarr", raise_on_failure=True) is True

    def test_the_default_after_a_raising_call(self, tmp_path):
        with pytest.raises(ConversionRefused):
            _convert(tmp_path / "a.zarr", raise_on_failure=True, bogus=1)
        assert _convert(tmp_path / "b.zarr", bogus=1) is False

    def test_a_second_failure_does_not_carry_the_first_reason(self, tmp_path):
        with pytest.raises(ConversionRefused, match="'first'"):
            _convert(tmp_path / "a.zarr", raise_on_failure=True, first=1)
        with pytest.raises(ConversionRefused) as raised:
            _convert(tmp_path / "a.zarr", raise_on_failure=True, second=1)
        assert "'first'" not in str(raised.value)
        assert "'second'" in str(raised.value)

    def test_no_recorder_is_left_open(self, tmp_path):
        _convert(tmp_path / "a.zarr", bogus=1)
        assert convert_module._FAILURES.get() is None
        with pytest.raises(ConversionRefused):
            _convert(tmp_path / "a.zarr", raise_on_failure=True, bogus=1)
        assert convert_module._FAILURES.get() is None

    def test_a_success_with_the_flag_returns_true(self, tmp_path):
        out = tmp_path / "out.zarr"
        assert _convert(out, raise_on_failure=True) is True
        assert out.is_dir()

    def test_several_reasons_join_with_a_newline(self):
        record = convert_module._FailureRecord(
            reasons=["one", "two"], kinds=["refused", "refused"]
        )
        with pytest.raises(ConversionRefused) as raised:
            convert_module._raise_recorded(record)
        assert str(raised.value) == "one\ntwo"

    def test_a_mixed_record_is_a_failure(self):
        record = convert_module._FailureRecord(
            reasons=["one", "two"], kinds=["refused", "error"]
        )
        with pytest.raises(ConversionFailed) as raised:
            convert_module._raise_recorded(record)
        assert not isinstance(raised.value, ConversionRefused)
