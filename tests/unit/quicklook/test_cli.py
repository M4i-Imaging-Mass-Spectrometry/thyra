"""The ``thyra ion-image`` and ``thyra report`` commands, end to end."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import zarr
from click.testing import CliRunner
from PIL import Image

from tests.fixtures.hand_built_store import MZ, make_table, write_store
from thyra.convert import convert_msi
from thyra.quicklook.cli import ion_image_command, report_command

_FIXTURES = Path(__file__).resolve().parents[2] / "data" / "fixtures"


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


def _files(folder: Path):
    return sorted(str(path.relative_to(folder)) for path in folder.rglob("*"))


def _assert_refused(result, exit_code: int, *phrases: str) -> None:
    """One ``Error:`` line carrying ``phrases``, the exit status, no traceback."""
    assert result.exit_code == exit_code, result.output
    assert "Traceback" not in result.output
    errors = [line for line in result.output.splitlines() if line.startswith("Error:")]
    assert len(errors) == 1, result.output
    for phrase in phrases:
        assert phrase in errors[0]
    if exit_code == 1:
        assert result.output.strip() == errors[0]
        assert isinstance(result.exception, SystemExit)


class TestIonImage:
    def test_it_writes_a_png_beside_the_store_named_after_it(self, runner, store_path):
        before = _files(store_path)

        result = runner.invoke(ion_image_command, [str(store_path), "505"])

        assert result.exit_code == 0, result.output
        written = store_path.parent / "sample_mz505.0000.png"
        image = Image.open(written)
        assert (image.format, image.mode) == ("PNG", "RGB")
        assert image.width > 480 and image.height > 480  # a 3 x 3 raster, enlarged
        assert _files(store_path) == before  # nothing is written into the store
        assert str(written) in result.output

    def test_it_says_what_it_added_up(self, runner, store_path):
        result = runner.invoke(
            ion_image_command, [str(store_path), "505", "--ppm", "2500"]
        )

        assert "m/z 505.0000 +/- 2500 ppm is 3 bin(s), 504.0000 to 506.0000" in (
            result.output
        )
        assert "Intensity 18 to 2120 over 8 pixels" in result.output

    def test_the_output_option_names_the_file_and_makes_its_folder(
        self, runner, store_path, tmp_path
    ):
        target = tmp_path / "figures" / "lipid.png"

        result = runner.invoke(
            ion_image_command, [str(store_path), "505", "-o", str(target)]
        )

        assert result.exit_code == 0, result.output
        assert Image.open(target).format == "PNG"
        assert not (store_path.parent / "sample_mz505.0000.png").exists()

    def test_a_blank_image_is_written_and_explained(self, runner, store_path):
        result = runner.invoke(ion_image_command, [str(store_path), "503"])

        assert result.exit_code == 0, result.output
        assert "No pixel has signal in this window" in result.output
        assert "nearby is at m/z 504.0000" in result.output
        assert (store_path.parent / "sample_mz503.0000.png").exists()

    def test_several_slice_tables_give_one_image_each(self, runner, tmp_path):
        tables = {
            "ds_z0": make_table(),
            "ds_z1": make_table(),
            "ds_z0_mobility": make_table(
                extra_var={"mobility": np.linspace(0.8, 1.2, MZ.size)}
            ),
        }
        path = write_store(tmp_path / "stack.zarr", tables)

        result = runner.invoke(ion_image_command, [str(path), "505"])

        assert result.exit_code == 0, result.output
        assert sorted(p.name for p in tmp_path.glob("*.png")) == [
            "stack_mz505.0000_ds_z0.png",
            "stack_mz505.0000_ds_z1.png",
        ]

    def test_the_main_table_is_drawn_whatever_the_tables_are_called(
        self, runner, tmp_path
    ):
        """One image, from the table with neither sibling column."""
        tables = {
            "a": make_table(extra_var={"precursor_mz": np.full(MZ.size, 600.0)}),
            "b_msms": make_table(),
            "c": make_table(extra_var={"mobility": np.linspace(0.8, 1.2, MZ.size)}),
        }
        path = write_store(tmp_path / "s.zarr", tables)

        result = runner.invoke(ion_image_command, [str(path), "505"])

        assert result.exit_code == 0, result.output
        assert [p.name for p in tmp_path.glob("*.png")] == ["s_mz505.0000.png"]
        assert result.output.startswith("b_msms: ")


class TestIonImageRefusals:
    def test_a_folder_that_is_not_a_store(self, runner, tmp_path):
        result = runner.invoke(ion_image_command, [str(tmp_path), "505"])

        _assert_refused(result, 1, "is not a Thyra store")

    def test_a_zarr_group_that_is_not_a_thyra_store(self, runner, tmp_path):
        zarr.open_group(str(tmp_path / "other.zarr"), mode="w")

        result = runner.invoke(ion_image_command, [str(tmp_path / "other.zarr"), "505"])

        _assert_refused(result, 1, "is not a Thyra store", "no 'tables' group")

    def test_a_path_that_does_not_exist_is_a_usage_error(self, runner, tmp_path):
        result = runner.invoke(ion_image_command, [str(tmp_path / "no.zarr"), "505"])

        _assert_refused(result, 2, "does not exist")

    def test_an_mz_outside_the_axis(self, runner, store_path):
        result = runner.invoke(ion_image_command, [str(store_path), "5000"])

        _assert_refused(result, 1, "m/z 5000.0000 is outside", "500.0000 to 509.0000")
        assert not list(store_path.parent.glob("*.png"))

    def test_a_window_no_bin_falls_in(self, runner, store_path):
        result = runner.invoke(ion_image_command, [str(store_path), "505.4"])

        _assert_refused(
            result,
            1,
            "No m/z bin lies within 10 ppm of 505.4000",
            "nearest bin is at 505.0000",
            "--ppm 792 would reach it",
        )

    def test_a_store_with_only_sibling_tables(self, runner, tmp_path):
        table = make_table(extra_var={"mobility": np.linspace(0.8, 1.2, MZ.size)})
        path = write_store(tmp_path / "s.zarr", {"ds_z0_mobility": table})

        result = runner.invoke(ion_image_command, [str(path), "505"])

        _assert_refused(result, 1, "holds no table of spectra")

    def test_an_output_inside_the_store_is_a_usage_error(self, runner, store_path):
        before = _files(store_path)

        result = runner.invoke(
            ion_image_command,
            [str(store_path), "505", "-o", str(store_path / "ion.png")],
        )

        _assert_refused(result, 2, "inside the store")
        assert _files(store_path) == before

    def test_an_output_that_is_a_folder_is_a_usage_error(
        self, runner, store_path, tmp_path
    ):
        result = runner.invoke(
            ion_image_command, [str(store_path), "505", "-o", str(tmp_path)]
        )

        _assert_refused(result, 2, "it is a folder")

    def test_an_output_inside_the_store_is_caught_however_it_is_spelt(
        self, runner, store_path
    ):
        inside = store_path.parent / "elsewhere" / ".." / "sample.zarr" / "ion.png"

        result = runner.invoke(
            ion_image_command, [f"{store_path}/", "505", "-o", str(inside)]
        )

        _assert_refused(result, 2, "inside the store")

    def test_a_file_that_cannot_be_written_is_one_line(
        self, runner, store_path, monkeypatch
    ):
        def refuse(self, data):
            raise PermissionError(13, "Permission denied")

        monkeypatch.setattr(Path, "write_bytes", refuse)

        result = runner.invoke(ion_image_command, [str(store_path), "505"])

        _assert_refused(result, 1, "Could not write", "Permission denied")

    def test_each_per_table_name_is_checked_before_anything_is_drawn(
        self, runner, tmp_path
    ):
        tables = {"ds_z0": make_table(), "ds_z1": make_table()}
        path = write_store(tmp_path / "stack.zarr", tables)
        (tmp_path / "out_ds_z1.png").mkdir()

        result = runner.invoke(
            ion_image_command, [str(path), "505", "-o", str(tmp_path / "out.png")]
        )

        _assert_refused(result, 2, "out_ds_z1.png", "it is a folder")
        assert not (tmp_path / "out_ds_z0.png").exists()

    @pytest.mark.parametrize("arguments", [["0"], ["nan"], ["505", "--ppm", "0"]])
    def test_a_number_that_cannot_be_a_mass_or_a_tolerance(
        self, runner, store_path, arguments
    ):
        result = runner.invoke(ion_image_command, [str(store_path), *arguments])

        _assert_refused(result, 2, "must be a number above zero")


class TestReport:
    def test_it_writes_one_html_file_beside_the_store(self, runner, store_path):
        before = _files(store_path)

        result = runner.invoke(report_command, [str(store_path)])

        assert result.exit_code == 0, result.output
        written = store_path.parent / "sample_report.html"
        page = written.read_text(encoding="utf-8")
        assert page.startswith("<!DOCTYPE html>") and "<h1>sample</h1>" in page
        assert "http://" not in page and "https://" not in page
        assert str(store_path.parent) not in page
        assert _files(store_path) == before
        assert f"Wrote {written}" in result.output

    def test_the_output_option_names_the_file(self, runner, store_path, tmp_path):
        target = tmp_path / "out" / "overview.html"

        result = runner.invoke(report_command, [str(store_path), "-o", str(target)])

        assert result.exit_code == 0, result.output
        assert target.read_text(encoding="utf-8").startswith("<!DOCTYPE html>")

    def test_a_file_that_cannot_be_written_is_one_line(
        self, runner, store_path, monkeypatch
    ):
        def refuse(self, data):
            raise PermissionError(13, "Permission denied")

        monkeypatch.setattr(Path, "write_bytes", refuse)

        result = runner.invoke(report_command, [str(store_path)])

        _assert_refused(result, 1, "Could not write", "Permission denied")

    def test_a_folder_that_is_not_a_store(self, runner, tmp_path):
        result = runner.invoke(report_command, [str(tmp_path)])

        _assert_refused(result, 1, "is not a Thyra store")
        assert not list(tmp_path.parent.glob("*_report.html"))

    def test_a_store_with_only_sibling_tables(self, runner, tmp_path):
        table = make_table(extra_var={"precursor_mz": np.full(MZ.size, 600.0)})
        path = write_store(tmp_path / "s.zarr", {"ds_z0_msms": table})

        result = runner.invoke(report_command, [str(path)])

        _assert_refused(result, 1, "holds no table of spectra")

    def test_a_path_that_does_not_exist_is_a_usage_error(self, runner, tmp_path):
        result = runner.invoke(report_command, [str(tmp_path / "no.zarr")])

        _assert_refused(result, 2, "does not exist")

    def test_an_output_inside_the_store_is_a_usage_error(self, runner, store_path):
        result = runner.invoke(
            report_command, [str(store_path), "-o", str(store_path / "r.html")]
        )

        _assert_refused(result, 2, "inside the store")


@pytest.fixture(scope="module")
def converted(tmp_path_factory) -> Path:
    """A store converted from a committed fixture: four pixels, one sibling table."""
    output = tmp_path_factory.mktemp("quicklook") / "mobility_continuous.zarr"
    assert convert_msi(
        str(_FIXTURES / "mobility_continuous.imzML"), str(output)
    ), "conversion fixture failed"
    return output


class TestOnAConvertedStore:
    def test_ion_image_draws_the_table_of_spectra_only(self, runner, converted):
        root = zarr.open_group(str(converted), mode="r")
        assert sorted(root["tables"].keys()) == [
            "msi_dataset_z0",
            "msi_dataset_z0_mobility",
        ]
        mz = np.asarray(root["tables/msi_dataset_z0/var/mz"][:])
        mean = np.asarray(root["tables/msi_dataset_z0/uns/average_spectrum"][:])
        target = tuple(converted.parent.glob("*.png"))
        assert not target

        result = runner.invoke(
            ion_image_command, [str(converted), f"{mz[int(np.argmax(mean))]:.5f}"]
        )

        assert result.exit_code == 0, result.output
        assert result.output.startswith("msi_dataset_z0: ")
        assert "over 4 pixels" in result.output
        (written,) = converted.parent.glob("*.png")
        assert Image.open(written).mode == "RGB"

    def test_report_shows_the_store_and_keeps_its_source_path_out(
        self, runner, converted
    ):
        result = runner.invoke(report_command, [str(converted)])

        assert result.exit_code == 0, result.output
        page = (converted.parent / "mobility_continuous_report.html").read_text(
            encoding="utf-8"
        )
        assert "<td>mobility_continuous.imzML</td>" in page
        assert str(_FIXTURES) not in page and _FIXTURES.name not in page
        assert _FIXTURES.as_posix() not in page
        assert "<th>msi_dataset_z0_mobility</th><td>ion mobility</td>" in page
        assert ">TIC image<" in page and ">Mean spectrum<" in page
        assert "http://" not in page and "https://" not in page


class TestDispatch:
    @pytest.mark.parametrize(
        "name, phrase",
        [("ion-image", "Write the ion image"), ("report", "one-file HTML overview")],
    )
    def test_thyra_reaches_the_subcommand(self, monkeypatch, capsys, name, phrase):
        from thyra.__main__ import cli

        monkeypatch.setattr("sys.argv", ["thyra", name, "--help"], raising=False)
        with pytest.raises(SystemExit) as excinfo:
            cli()

        assert excinfo.value.code == 0
        out = capsys.readouterr().out
        assert phrase in out and f"Usage: thyra {name}" in out

    def test_the_conversion_help_names_both(self, monkeypatch, capsys):
        from thyra.__main__ import cli

        monkeypatch.setattr("sys.argv", ["thyra", "--help"], raising=False)
        with pytest.raises(SystemExit):
            cli()

        out = capsys.readouterr().out
        assert "thyra ion-image STORE MZ" in out and "thyra report STORE" in out

    def test_the_default_tolerance_is_stated_once(self):
        """The command repeats it so that building the command imports nothing."""
        from thyra.quicklook import cli as commands
        from thyra.quicklook.ion_image import DEFAULT_PPM

        assert commands._DEFAULT_PPM == DEFAULT_PPM

    def test_dispatching_costs_no_imaging_import(self):
        """``thyra.__main__`` imports the commands' module on every run."""
        code = (
            "import sys, json; import thyra.quicklook.cli; "
            "print(json.dumps(sorted(m for m in sys.modules "
            "if m.split('.')[0] in ('PIL', 'zarr', 'anndata', 'numpy'))))"
        )
        proc = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, check=True
        )

        assert json.loads(proc.stdout) == []
