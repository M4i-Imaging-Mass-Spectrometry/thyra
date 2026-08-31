"""The `thyra validate` / `thyra export-metaspace` subcommands."""

import json

import pytest
from click.testing import CliRunner

from thyra.metadata.schema import build_msi_metadata
from thyra.metadata.schema.cli import export_metaspace_command, validate_command


def _make_runner() -> CliRunner:
    """A runner whose results expose stderr on every supported click.

    click < 8.2 mixes stderr into output unless asked not to (and
    ``result.stderr`` raises); 8.2 removed the ``mix_stderr`` kwarg and
    always separates.  The project supports both (see the click pin in
    pyproject.toml).
    """
    try:
        return CliRunner(mix_stderr=False)  # type: ignore[call-arg]
    except TypeError:
        return CliRunner()


@pytest.fixture
def runner():
    return _make_runner()


def _write_doc(tmp_path, name="meta.json", mutate=None):
    doc = build_msi_metadata(
        None, pixel_size_um=(20.0, 20.0), source_format="imzml"
    ).to_uns_dict()
    if mutate:
        mutate(doc)
    path = tmp_path / name
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path


class TestValidateCommand:
    def test_valid_document_exits_zero(self, runner, tmp_path):
        result = runner.invoke(validate_command, [str(_write_doc(tmp_path))])
        assert result.exit_code == 0
        assert "OK" in result.output

    def test_invalid_document_exits_one_and_names_the_field(self, runner, tmp_path):
        def corrupt(doc):
            doc["ms_analysis"]["pixel_size_um"]["x"] = -5.0

        path = _write_doc(tmp_path, mutate=corrupt)
        result = runner.invoke(validate_command, [str(path)])
        assert result.exit_code == 1
        assert "ms_analysis.pixel_size_um.x" in result.output

    def test_missing_input_exits_two(self, runner, tmp_path):
        result = runner.invoke(validate_command, [str(tmp_path / "nope.json")])
        assert result.exit_code == 2

    def test_merge_overlays_user_fields(self, runner, tmp_path):
        # An overlay can also break the document; that must be caught.
        path = _write_doc(tmp_path)
        overlay = tmp_path / "user.json"
        overlay.write_text(
            json.dumps({"ms_analysis": {"polarity": "sideways"}}),
            encoding="utf-8",
        )
        result = runner.invoke(validate_command, [str(path), "--merge", str(overlay)])
        assert result.exit_code == 1
        assert "polarity" in result.output

    def test_json_report_is_machine_readable(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        result = runner.invoke(validate_command, [str(path), "--json"])
        assert result.exit_code == 0
        report = json.loads(result.stdout)
        assert report[path.name]["valid"] is True
        assert report[path.name]["issues"] == []


class TestExportMetaspaceCommand:
    def test_writes_the_submission_json(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        output = tmp_path / "out.json"
        result = runner.invoke(export_metaspace_command, [str(path), "-o", str(output)])
        assert result.exit_code == 0
        document = json.loads(output.read_text(encoding="utf-8"))
        assert document["Data_Type"] == "Imaging MS"
        assert document["MS_Analysis"]["Pixel_Size"] == {"Xaxis": 20, "Yaxis": 20}

    def test_default_output_lands_next_to_the_input(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        result = runner.invoke(export_metaspace_command, [str(path)])
        assert result.exit_code == 0
        assert (tmp_path / "meta.metaspace.json").exists()

    def test_stdout_output(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        result = runner.invoke(export_metaspace_command, [str(path), "-o", "-"])
        assert result.exit_code == 0
        assert json.loads(result.stdout)["Data_Type"] == "Imaging MS"

    def test_missing_required_fields_are_warned_on_stderr(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        result = runner.invoke(export_metaspace_command, [str(path), "-o", "-"])
        assert result.exit_code == 0
        assert "Organism" in (result.stderr or "")

    def test_merge_completes_the_submission(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        overlay = tmp_path / "user.json"
        overlay.write_text(
            json.dumps(
                {
                    "sample": {
                        "organism": "Mus musculus",
                        "organism_part": "liver",
                        "condition": "wildtype",
                    }
                }
            ),
            encoding="utf-8",
        )
        result = runner.invoke(
            export_metaspace_command,
            [str(path), "--merge", str(overlay), "-o", "-"],
        )
        assert result.exit_code == 0
        document = json.loads(result.stdout)
        assert document["Sample_Information"]["Organism"] == "Mus musculus"

    def test_non_conforming_document_is_refused(self, runner, tmp_path):
        def corrupt(doc):
            doc["ms_analysis"]["pixel_size_um"]["x"] = -5.0

        path = _write_doc(tmp_path, mutate=corrupt)
        result = runner.invoke(export_metaspace_command, [str(path), "-o", "-"])
        assert result.exit_code == 1


class TestDispatcher:
    def test_subcommands_are_dispatched(self, monkeypatch, capsys):
        from thyra.__main__ import cli

        monkeypatch.setattr("sys.argv", ["thyra", "validate", "--help"], raising=False)
        with pytest.raises(SystemExit) as excinfo:
            cli()
        assert excinfo.value.code == 0
        assert "Validate MSI metadata" in capsys.readouterr().out

    def test_conversion_interface_still_owns_bare_help(self, monkeypatch, capsys):
        from thyra.__main__ import cli

        monkeypatch.setattr("sys.argv", ["thyra", "--help"], raising=False)
        with pytest.raises(SystemExit) as excinfo:
            cli()
        assert excinfo.value.code == 0
        out = capsys.readouterr().out
        assert "INPUT" in out and "OUTPUT" in out
        assert "export-metaspace" in out
