"""Every file a Thyra command writes follows the store path's rule (#525).

``--log-file``, ``thyra metadata -o``, ``thyra export-metaspace -o`` and
``thyra-check-ontology --output`` opened their file without a check. A
path in a folder that did not exist yet, or a path that was a folder,
ended in a traceback -- after the metadata read or the whole ontology
check, and for ``--log-file`` before the conversion started. Now a
missing folder is created when the file is written, and a path that is a
folder, or that lies under a regular file, is refused before any work
with one line and exit 2.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Dict

import pytest
from click.testing import CliRunner

from thyra.metadata.schema import build_msi_metadata
from thyra.utils.output_files import output_file_problem


def _runner() -> CliRunner:
    try:
        return CliRunner(mix_stderr=False)  # type: ignore[call-arg]
    except TypeError:
        return CliRunner()


def _text(result) -> str:
    try:
        stderr = result.stderr
    except ValueError:
        stderr = ""
    return result.output + stderr


@pytest.fixture
def imzml(create_minimal_imzml) -> Path:
    path, _, _, _ = create_minimal_imzml
    return path


@pytest.fixture
def document(tmp_path) -> Path:
    doc = build_msi_metadata(
        None, pixel_size_um=(20.0, 20.0), source_format="imzml"
    ).to_uns_dict()
    path = tmp_path / "doc.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path


def _log_file(target: Path, imzml: Path, tmp_path: Path):
    from thyra.__main__ import main

    args = [str(imzml), str(tmp_path / "s.zarr"), "--pixel-size", "25"]
    return _runner().invoke(main, args + ["--log-file", str(target)])


def _metadata(target: Path, imzml: Path, tmp_path: Path):
    from thyra.metadata.schema.cli import metadata_command

    return _runner().invoke(metadata_command, [str(imzml), "-o", str(target)])


def _metaspace(target: Path, document: Path, tmp_path: Path):
    from thyra.metadata.schema.cli import export_metaspace_command

    return _runner().invoke(
        export_metaspace_command, [str(document), "-o", str(target)]
    )


_COMMANDS: Dict[str, Callable] = {
    "log-file": _log_file,
    "metadata": _metadata,
    "export-metaspace": _metaspace,
}


def _source(name: str, imzml: Path, document: Path) -> Path:
    return document if name == "export-metaspace" else imzml


@pytest.mark.parametrize("name", list(_COMMANDS))
class TestTheCommandLineOutputs:
    def test_a_missing_folder_is_created(self, name, imzml, document, tmp_path):
        target = tmp_path / "new" / "deeper" / "out.file"

        result = _COMMANDS[name](target, _source(name, imzml, document), tmp_path)

        assert "Traceback" not in _text(result)
        assert target.is_file(), _text(result)

    def test_an_existing_folder_is_refused(self, name, imzml, document, tmp_path):
        target = tmp_path / "a_folder"
        target.mkdir()

        result = _COMMANDS[name](target, _source(name, imzml, document), tmp_path)

        assert result.exit_code == 2
        assert str(target) in _text(result)
        assert "Traceback" not in _text(result)
        assert not (tmp_path / "s.zarr").exists()

    def test_a_path_under_a_file_is_refused(self, name, imzml, document, tmp_path):
        blocker = tmp_path / "plain_file"
        blocker.write_text("x")
        target = blocker / "out.file"

        result = _COMMANDS[name](target, _source(name, imzml, document), tmp_path)

        assert result.exit_code == 2
        assert "is a file" in _text(result)
        assert "Traceback" not in _text(result)


class TestCheckedBeforeAnyWork:
    def test_metadata_on_an_unreadable_source_creates_no_folder(self, tmp_path):
        from thyra.metadata.schema.cli import metadata_command

        empty = tmp_path / "empty"
        empty.mkdir()
        target = tmp_path / "newX" / "m.json"

        result = _runner().invoke(metadata_command, [str(empty), "-o", str(target)])

        assert result.exit_code == 1
        assert "Could not read" in _text(result)
        assert not (tmp_path / "newX").exists()

    def test_metadata_to_stdout_is_unchanged(self, imzml):
        from thyra.metadata.schema.cli import metadata_command

        result = _runner().invoke(metadata_command, [str(imzml), "-o", "-"])

        assert result.output.lstrip().startswith("{")

    def test_the_export_default_beside_the_store_is_checked(self, document):
        from thyra.metadata.schema.cli import export_metaspace_command

        default = document.with_name(document.stem + ".metaspace.json")
        default.mkdir()

        result = _runner().invoke(export_metaspace_command, [str(document)])

        assert result.exit_code == 2
        assert str(default) in _text(result)

    @pytest.mark.parametrize("where", ["inside", "at"])
    def test_a_log_file_at_or_inside_the_store_is_refused(self, imzml, tmp_path, where):
        from thyra.__main__ import main

        store = tmp_path / "s.zarr"
        log = store / "run.log" if where == "inside" else store
        result = _runner().invoke(
            main,
            [str(imzml), str(store), "--pixel-size", "25"] + ["--log-file", str(log)],
        )

        assert result.exit_code == 2
        assert not store.exists()


class TestCheckOntology:
    @pytest.fixture
    def sample(self, tmp_path) -> Path:
        path = tmp_path / "sample.imzML"
        path.write_text(
            '<?xml version="1.0"?><mzML xmlns="http://psi.hupo.org/ms/mzml">'
            "<fileDescription><fileContent>"
            '<cvParam cvRef="IMS" accession="IMS:1000042" '
            'name="max count of pixels x" value="3"/>'
            "</fileContent></fileDescription></mzML>"
        )
        return path

    def test_a_missing_folder_is_created(self, sample, tmp_path):
        from thyra.tools.check_ontology import main

        target = tmp_path / "no_such" / "r.json"

        assert main([str(sample), "--output", str(target)]) == 0
        assert json.loads(target.read_text(encoding="utf-8"))

    @pytest.mark.parametrize("shape", ["folder", "under_a_file"])
    @pytest.mark.parametrize("input_kind", ["file", "folder"])
    def test_a_refused_path_never_runs_the_check(
        self, sample, tmp_path, monkeypatch, capsys, shape, input_kind
    ):
        from thyra.metadata import validator
        from thyra.tools.check_ontology import main

        def never(*args, **kwargs):
            raise AssertionError("the check ran on a refused output path")

        monkeypatch.setattr(validator.ImzMLOntologyValidator, "validate_file", never)
        monkeypatch.setattr(
            validator.ImzMLOntologyValidator, "validate_directory", never
        )
        if shape == "folder":
            target = tmp_path / "a_folder"
            target.mkdir()
        else:
            (tmp_path / "plain_file").write_text("x")
            target = tmp_path / "plain_file" / "r.json"
        source = sample if input_kind == "file" else tmp_path

        assert main([str(source), "--output", str(target)]) == 2
        captured = capsys.readouterr()
        assert captured.out == ""
        assert str(target) in captured.err


class TestTheRule:
    def test_a_writable_path_passes(self, tmp_path):
        assert output_file_problem(tmp_path / "a" / "b.json") is None

    def test_a_folder_is_refused(self, tmp_path):
        assert "is a folder" in output_file_problem(tmp_path)  # type: ignore[operator]

    def test_a_path_under_a_file_is_refused(self, tmp_path):
        (tmp_path / "f").write_text("x")
        problem = output_file_problem(tmp_path / "f" / "g" / "h.json")
        assert problem is not None and "is a file" in problem
