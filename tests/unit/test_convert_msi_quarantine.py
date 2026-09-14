"""A failed conversion clears its own destination (issue #293).

``convert_msi`` refuses to write to a path that already exists, and a
conversion that fails part-way through writing left an incomplete
``.zarr`` sitting at that path. So the run that failed was the run that
made the same command fail again -- the second time with "Destination
already exists", a message about the destination rather than about
whatever actually went wrong. Clearing it meant deleting a directory by
hand, and the store looks like a plausible artifact while being
unopenable (``spatialdata.read_zarr()`` raises on it).

The move existed, in ``thyra/__main__.py``, and only the CLI called it.
Every library caller -- Ousia converts through ``convert_msi``, never
through the CLI -- got nothing. Four CLI tests covered the exit status
and all four monkeypatch ``convert_msi`` away, so none of them reached
the gap.

Not to be confused with the Ctrl-C work (#245), which is exception
*routing*: three handlers that turn an interrupt into ``False``. They
all return through this same ``finally``, so they are covered by it
rather than duplicated.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from thyra import convert as convert_module
from thyra.convert import convert_msi, quarantine_partial_output
from thyra.errors import ConversionRefused


def _half_written_store(path: Path) -> None:
    """What a conversion that died mid-write leaves behind."""
    path.mkdir(parents=True)
    (path / "zarr.json").write_text('{"zarr_format": 3, "node_type": "group"}')


@pytest.fixture
def imzml(create_minimal_imzml):
    imzml_path, _, _, _ = create_minimal_imzml
    return imzml_path


def _fail_after_writing(out: Path, how):
    """A ``_perform_conversion_with_cleanup`` that writes, then fails."""

    def replacement(converter, reader):
        _half_written_store(out)
        return how()

    return replacement


class TestTheDestinationIsCleared:
    def test_a_failed_conversion_leaves_nothing_at_the_destination(
        self, imzml, tmp_path, monkeypatch
    ):
        out = tmp_path / "out.zarr"
        monkeypatch.setattr(
            convert_module,
            "_perform_conversion_with_cleanup",
            _fail_after_writing(out, lambda: False),
        )

        assert convert_msi(str(imzml), str(out), pixel_size_um=25.0) is False
        assert not out.exists()

    def test_the_partial_store_is_kept_for_diagnosis(
        self, imzml, tmp_path, monkeypatch
    ):
        out = tmp_path / "out.zarr"
        monkeypatch.setattr(
            convert_module,
            "_perform_conversion_with_cleanup",
            _fail_after_writing(out, lambda: False),
        )

        convert_msi(str(imzml), str(out), pixel_size_um=25.0)

        failed = tmp_path / "out.zarr.failed"
        assert failed.is_dir()
        assert (failed / "zarr.json").exists()

    def test_a_retry_is_accepted(self, imzml, tmp_path, monkeypatch):
        """The point of the whole thing: the same command runs again."""
        out = tmp_path / "out.zarr"
        monkeypatch.setattr(
            convert_module,
            "_perform_conversion_with_cleanup",
            _fail_after_writing(out, lambda: False),
        )
        convert_msi(str(imzml), str(out), pixel_size_um=25.0)

        monkeypatch.setattr(
            convert_module, "_perform_conversion_with_cleanup", lambda c, r: True
        )
        assert convert_msi(str(imzml), str(out), pixel_size_um=25.0) is True

    @pytest.mark.parametrize(
        "how",
        [
            pytest.param(lambda: False, id="returns_false"),
            pytest.param(
                lambda: (_ for _ in ()).throw(ConversionRefused("nope")),
                id="refuses",
            ),
            pytest.param(
                lambda: (_ for _ in ()).throw(RuntimeError("boom")),
                id="raises",
            ),
            pytest.param(
                lambda: (_ for _ in ()).throw(KeyboardInterrupt()),
                id="interrupted",
            ),
        ],
    )
    def test_every_failure_route_clears_it(self, imzml, tmp_path, monkeypatch, how):
        out = tmp_path / "out.zarr"
        monkeypatch.setattr(
            convert_module,
            "_perform_conversion_with_cleanup",
            _fail_after_writing(out, how),
        )

        assert convert_msi(str(imzml), str(out), pixel_size_um=25.0) is False
        assert not out.exists()
        assert (tmp_path / "out.zarr.failed").is_dir()


class TestItDoesNotTouchWhatItDoesNotOwn:
    def test_a_pre_existing_destination_is_refused_and_left_alone(
        self, imzml, tmp_path
    ):
        """The refusal fires *before* the destination is this run's, so
        a store written by an earlier successful run must not be renamed
        out from under the user."""
        out = tmp_path / "out.zarr"
        _half_written_store(out)

        assert convert_msi(str(imzml), str(out), pixel_size_um=25.0) is False
        assert out.is_dir()
        assert (out / "zarr.json").exists()
        assert not (tmp_path / "out.zarr.failed").exists()

    def test_a_successful_conversion_keeps_its_output(
        self, imzml, tmp_path, monkeypatch
    ):
        out = tmp_path / "out.zarr"

        def succeed(converter, reader):
            _half_written_store(out)
            return True

        monkeypatch.setattr(convert_module, "_perform_conversion_with_cleanup", succeed)

        assert convert_msi(str(imzml), str(out), pixel_size_um=25.0) is True
        assert out.is_dir()
        assert not (tmp_path / "out.zarr.failed").exists()

    def test_a_missing_input_renames_nothing(self, tmp_path):
        out = tmp_path / "out.zarr"
        assert (
            convert_msi(str(tmp_path / "absent.imzML"), str(out), pixel_size_um=25.0)
            is False
        )
        assert list(tmp_path.iterdir()) == []


class TestTheNumbering:
    def test_a_second_failure_does_not_overwrite_the_first(self, tmp_path):
        out = tmp_path / "out.zarr"

        _half_written_store(out)
        quarantine_partial_output(out)
        _half_written_store(out)
        quarantine_partial_output(out)

        assert (tmp_path / "out.zarr.failed").is_dir()
        assert (tmp_path / "out.zarr.failed2").is_dir()
        assert not out.exists()

    def test_nothing_to_move_is_a_no_op(self, tmp_path):
        quarantine_partial_output(tmp_path / "absent.zarr")
        assert list(tmp_path.iterdir()) == []
