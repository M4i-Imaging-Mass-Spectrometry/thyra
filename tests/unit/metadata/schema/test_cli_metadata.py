"""The `thyra validate` / `thyra export-metaspace` subcommands."""

import json
import shutil
import sys
import tempfile
from pathlib import Path

import pytest
from click.testing import CliRunner

from thyra.metadata.schema import build_msi_metadata
from thyra.metadata.schema.cli import export_metaspace_command, validate_command
from thyra.utils.windows_paths import WINDOWS_MAX_PATH, to_extended_length_path


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


class TestJsonFilesNotInPlainUtf8:
    """A JSON input saved by Windows PowerShell is read; another encoding is refused.

    PowerShell 5.1 writes ``thyra metadata src > meta.json`` as UTF-16, and
    UTF-8 with a byte-order mark when its Out-File default is utf8. Both
    ended in a UnicodeDecodeError traceback before.
    """

    def _reencode(self, path, encoding):
        path.write_bytes(path.read_text(encoding="utf-8").encode(encoding))
        return path

    @pytest.mark.parametrize("encoding", ["utf-16", "utf-8-sig", "utf-16-be"])
    def test_validate_reads_it(self, runner, tmp_path, encoding):
        path = self._reencode(_write_doc(tmp_path), encoding)
        result = runner.invoke(validate_command, [str(path)])
        assert result.exit_code == 0, result.output
        assert "OK" in result.output

    def test_export_metaspace_reads_it(self, runner, tmp_path):
        path = self._reencode(_write_doc(tmp_path), "utf-16")
        result = runner.invoke(export_metaspace_command, [str(path), "-o", "-"])
        assert result.exit_code == 0
        assert json.loads(result.stdout)["Data_Type"] == "Imaging MS"

    def test_a_merged_overlay_is_read_too(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        overlay = tmp_path / "user.json"
        overlay.write_text(
            json.dumps({"ms_analysis": {"polarity": "sideways"}}), encoding="utf-16"
        )
        result = runner.invoke(validate_command, [str(path), "--merge", str(overlay)])
        assert result.exit_code == 1
        assert "polarity" in result.output

    def test_an_ansi_overlay_is_refused_in_one_line(self, runner, tmp_path):
        """An editor saving in the Windows code page; one degree sign is enough."""
        path = _write_doc(tmp_path)
        overlay = tmp_path / "user.json"
        overlay.write_text(
            json.dumps({"sample": {"condition": "4 °C"}}, ensure_ascii=False),
            encoding="cp1252",
        )
        result = runner.invoke(validate_command, [str(path), "--merge", str(overlay)])
        assert result.exit_code == 1
        assert not isinstance(result.exception, UnicodeDecodeError)
        message = result.output + (result.stderr or "")
        assert "Could not read" in message
        assert "not UTF-8, UTF-16 or UTF-32" in message
        assert "Traceback" not in message


def _deep_json(path, levels=100_000):
    """A JSON object nested past the parser's recursion limit."""
    path.write_text('{"a": ' + "[" * levels + "]" * levels + "}", encoding="utf-8")
    return path


class TestJsonNestedTooDeeply:
    """Issue #494 (D3): a RecursionError traceback, not a one-line refusal."""

    def _assert_refused(self, result):
        assert result.exit_code == 1
        assert not isinstance(result.exception, RecursionError)
        message = result.output + (result.stderr or "")
        assert "nested too deeply" in message
        assert "Traceback" not in message

    def test_validate_refuses_it(self, runner, tmp_path):
        deep = _deep_json(tmp_path / "deep.json")
        self._assert_refused(runner.invoke(validate_command, [str(deep)]))

    def test_validate_refuses_it_as_an_overlay(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        deep = _deep_json(tmp_path / "deep.json")
        result = runner.invoke(validate_command, [str(path), "--merge", str(deep)])
        self._assert_refused(result)

    def test_export_metaspace_refuses_it(self, runner, tmp_path):
        deep = _deep_json(tmp_path / "deep.json")
        out = tmp_path / "out.json"
        result = runner.invoke(export_metaspace_command, [str(deep), "-o", str(out)])
        self._assert_refused(result)
        assert not out.exists()

    def test_a_depth_the_parser_reads_is_refused_too(self, runner, tmp_path):
        """Python 3.14 on Linux parses 100,000 levels; the bound is Thyra's own."""
        path = _write_doc(tmp_path)
        deep = _deep_json(tmp_path / "deep.json", levels=200)
        result = runner.invoke(validate_command, [str(path), "--merge", str(deep)])
        self._assert_refused(result)
        assert "more than 100 levels" in result.output + (result.stderr or "")

    def test_a_document_within_the_bound_is_read(self, runner, tmp_path):
        """98 array levels inside two object levels: 100 in all."""
        path = _write_doc(tmp_path)
        overlay = tmp_path / "deep.json"
        overlay.write_text(
            '{"sample": {"a": ' + "[" * 98 + "]" * 98 + "}}", encoding="utf-8"
        )
        result = runner.invoke(validate_command, [str(path), "--merge", str(overlay)])
        message = result.output + (result.stderr or "")
        assert "nested too deeply" not in message
        assert "sample.a" in message  # the schema has no such field

    def test_export_metaspace_refuses_it_as_an_overlay(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        deep = _deep_json(tmp_path / "deep.json")
        out = tmp_path / "out.json"
        result = runner.invoke(
            export_metaspace_command, [str(path), "--merge", str(deep), "-o", str(out)]
        )
        self._assert_refused(result)
        assert not out.exists()


class TestNonFiniteNumbers:
    """Issue #494 (D2): an infinite value passed and was written as ``Infinity``.

    ``json.loads`` reads ``1e999`` and the tokens ``Infinity`` and ``NaN``
    as floats, and ``json.dumps`` writes them back as tokens that JSON
    does not have. A JSON input holding one is refused in one line that
    names the field.
    """

    def _overlay(self, tmp_path, token):
        overlay = tmp_path / "user.json"
        overlay.write_text(
            '{"ms_analysis": {"pixel_size_um": {"x": %s, "y": 20.0}}}' % token,
            encoding="utf-8",
        )
        return overlay

    @staticmethod
    def _assert_refused(result, location):
        assert result.exit_code == 1
        message = result.output + (result.stderr or "")
        assert "Could not read" in message
        assert f"{location} is not a finite number" in message
        assert "Traceback" not in message

    @pytest.mark.parametrize("token", ["1e999", "Infinity", "-Infinity", "NaN"])
    def test_validate_refuses_an_overlay_holding_one(self, runner, tmp_path, token):
        path = _write_doc(tmp_path)
        overlay = self._overlay(tmp_path, token)
        result = runner.invoke(validate_command, [str(path), "--merge", str(overlay)])
        self._assert_refused(result, "ms_analysis.pixel_size_um.x")

    def test_validate_refuses_a_document_holding_one(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        path.write_text(
            path.read_text(encoding="utf-8").replace('"x": 20.0', '"x": Infinity'),
            encoding="utf-8",
        )
        assert "Infinity" in path.read_text(encoding="utf-8")
        result = runner.invoke(validate_command, [str(path)])
        self._assert_refused(result, "ms_analysis.pixel_size_um.x")

    def test_a_value_in_a_list_is_located_too(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        overlay = tmp_path / "user.json"
        overlay.write_text(
            '{"processing": [{"parameters": {"a": 1.0}}, {"parameters": {"b": NaN}}]}',
            encoding="utf-8",
        )
        result = runner.invoke(validate_command, [str(path), "--merge", str(overlay)])
        self._assert_refused(result, "processing.1.parameters.b")

    def test_export_metaspace_writes_nothing(self, runner, tmp_path):
        path = _write_doc(tmp_path)
        overlay = self._overlay(tmp_path, "1e999")
        out = tmp_path / "out.json"
        result = runner.invoke(
            export_metaspace_command,
            [str(path), "--merge", str(overlay), "-o", str(out)],
        )
        self._assert_refused(result, "ms_analysis.pixel_size_um.x")
        assert not out.exists()

    def test_the_json_writer_refuses_a_document_too_deep_to_write(self, monkeypatch):
        """The encoder's own recursion limit is refused in one line too."""
        import click

        from thyra.metadata.schema import cli

        def too_deep(*args, **kwargs):
            raise RecursionError("maximum recursion depth exceeded")

        monkeypatch.setattr(cli.json, "dumps", too_deep)
        with pytest.raises(click.ClickException, match="nested too deeply"):
            cli._render_json({"x": 1})

    def test_the_json_writer_refuses_them_in_one_line(self):
        """The backstop for a value from a store or a raw source."""
        import click

        from thyra.metadata.schema.cli import _render_json

        assert json.loads(_render_json({"x": 1.5})) == {"x": 1.5}
        for value in (float("inf"), float("-inf"), float("nan")):
            with pytest.raises(click.ClickException, match="x.0 is not a finite"):
                _render_json({"x": [value]})


class TestABinCountThatIsNotANumber:
    """Issue #494 (D1): ``int()`` raised and the verdict was lost."""

    def test_validate_prints_its_verdict(self, runner, tmp_path):
        import anndata as ad
        import numpy as np
        import zarr

        from thyra.metadata.schema import mass_axis_sha256

        mz = np.array([100.0, 200.0, 300.0])
        root = zarr.open_group(str(tmp_path / "store.zarr"), mode="a")
        table = root.create_group("tables").create_group("t")
        table.create_group("var").create_array("mz", data=mz)
        block = _json_doc()
        block["provenance"]["mass_axis"] = {
            "sha256": mass_axis_sha256(mz),
            "n_bins": "many",
        }
        ad.io.write_elem(table.create_group("uns"), "msi_metadata", block)

        result = runner.invoke(
            validate_command, [str(tmp_path / "store.zarr"), "--json"]
        )

        assert result.exit_code == 1
        assert not isinstance(result.exception, ValueError)
        report = json.loads(result.stdout)
        assert report["t"]["valid"] is False
        locations = {issue["location"] for issue in report["t"]["issues"]}
        assert "provenance.mass_axis" in locations


def _json_doc():
    """A conforming document as plain JSON types, ready for ``write_elem``."""
    return json.loads(
        json.dumps(
            build_msi_metadata(
                None, pixel_size_um=(20.0, 20.0), source_format="imzml"
            ).to_uns_dict()
        )
    )


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


class TestTheStorePathIsPreparedInsideTheCommand:
    """Issue #257: click's existence check ran before the path was prepared.

    The converter writes a store to a deep Windows path through an
    extended-length path, so ``os.path.exists`` on the plain spelling is
    False and ``click.Path(exists=True)`` refused a store the CLI itself
    had just written -- while ``validate_store`` on the prepared path
    validated it fine.
    """

    def test_a_missing_path_is_still_a_usage_error(self, runner, tmp_path):
        result = runner.invoke(validate_command, [str(tmp_path / "nope.zarr")])
        assert result.exit_code == 2

    def test_the_prepared_path_is_what_gets_read(self, tmp_path, monkeypatch):
        """Proven by preparation that redirects: the command must follow it."""
        from thyra.metadata.schema import cli as cli_module

        real = _write_doc(tmp_path, name="real.json")
        typed = tmp_path / "typed.json"

        monkeypatch.setattr(
            cli_module, "prepare_zarr_read_path", lambda path: real, raising=True
        )
        assert cli_module._resolve_store_path(typed) == real

    def test_export_prepares_it_too(self, runner, tmp_path, monkeypatch):
        from thyra.metadata.schema import cli as cli_module

        real = _write_doc(tmp_path, name="real.json")
        monkeypatch.setattr(
            cli_module, "prepare_zarr_read_path", lambda path: real, raising=True
        )

        result = runner.invoke(
            export_metaspace_command, [str(tmp_path / "typed.json"), "-o", "-"]
        )
        assert result.exit_code == 0, result.output


@pytest.mark.skipif(sys.platform != "win32", reason="the 260 character limit")
class TestAPathPastTheWindowsLimit:
    """The case the issue reproduced: a path click cannot even stat."""

    @pytest.fixture
    def deep_document(self):
        base = Path(tempfile.mkdtemp(prefix="thyra_deep_cli_"))
        try:
            deep = base
            while len(str(deep / "metadata.json")) <= WINDOWS_MAX_PATH:
                deep = deep / "nested_directory_segment"
                to_extended_length_path(deep).mkdir()
            document = deep / "metadata.json"
            to_extended_length_path(document).write_text(
                json.dumps(
                    build_msi_metadata(
                        None, pixel_size_um=(20.0, 20.0), source_format="imzml"
                    ).to_uns_dict()
                ),
                encoding="utf-8",
            )
            if document.exists():
                # A machine with LongPathsEnabled=1 -- the GitHub Windows
                # runner is one -- resolves the plain path fine, so the
                # case this reproduces cannot occur on it.
                pytest.skip("Windows long-path support is on; the limit does not apply")
            yield document
        finally:
            shutil.rmtree(to_extended_length_path(base), ignore_errors=True)

    def test_validate_reaches_it(self, runner, deep_document):
        result = runner.invoke(validate_command, [str(deep_document)])
        assert result.exit_code == 0, result.output
        assert "OK" in result.output
