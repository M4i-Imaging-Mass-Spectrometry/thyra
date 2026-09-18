"""The ``thyra metadata`` subcommand: a document from a raw source.

Run against the committed synthetic TDF acquisition, which has a raster
and therefore a pixel size, so this is the ordinary case: the document
conforms and the command exits 0.  The other case -- an acquisition with
no raster, whose document cannot conform -- is in
``tests/unit/test_metadata_without_a_raster.py``, beside the tables that
produce it.

No spectra are decoded and no vendor SDK is loaded, so these run
everywhere.
"""

import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from thyra.metadata.schema.cli import metadata_command

FIXTURE = (
    Path(__file__).resolve().parents[4]
    / "tests"
    / "data"
    / "fixtures"
    / "synthetic_tims.d"
)


def _make_runner() -> CliRunner:
    """A runner whose results expose stderr on every supported click.

    ``result.output`` mixes the two on click 8.2 and is stdout alone on
    8.1, so every assertion below reads ``stdout`` or ``stderr`` by name.
    """
    try:
        return CliRunner(mix_stderr=False)  # type: ignore[call-arg]
    except TypeError:
        return CliRunner()


@pytest.fixture
def runner() -> CliRunner:
    return _make_runner()


@pytest.fixture
def source(tmp_path) -> Path:
    """A writable copy of the committed synthetic acquisition."""
    target = tmp_path / "copy.d"
    shutil.copytree(FIXTURE, target)
    return target


class TestTheDocumentReachesStdout:
    def test_the_default_is_stdout(self, runner, source):
        result = runner.invoke(metadata_command, [str(source)])
        assert result.exit_code == 0, result.output
        document = json.loads(result.stdout)
        assert document["schema_version"]

    def test_a_dash_is_stdout_too(self, runner, source):
        result = runner.invoke(metadata_command, [str(source), "-o", "-"])
        assert result.exit_code == 0, result.output
        assert json.loads(result.stdout)["ms_analysis"]["pixel_size_um"]

    def test_nothing_is_written_beside_the_input(self, runner, source):
        before = sorted(p.name for p in source.parent.iterdir())
        runner.invoke(metadata_command, [str(source)])
        assert sorted(p.name for p in source.parent.iterdir()) == before


class TestTheDocumentDescribesTheAcquisition:
    @pytest.fixture
    def document(self, runner, source) -> dict:
        result = runner.invoke(metadata_command, [str(source)])
        assert result.exit_code == 0, result.output
        return json.loads(result.stdout)

    def test_an_imaging_source_has_a_pixel_size(self, document):
        pitch = document["ms_analysis"]["pixel_size_um"]
        assert pitch["x"] > 0 and pitch["y"] > 0

    def test_the_pitch_is_recorded_as_detected(self, document):
        assert document["provenance"]["pixel_size_source"] == "automatic"

    def test_the_source_format_is_recorded(self, document):
        assert document["provenance"]["source_format"] == "bruker"

    def test_nothing_has_been_done_to_the_data(self, document):
        # There is no conversion behind this document, so there is no
        # processing history -- unlike the same block read out of a store.
        assert "processing" not in document


class TestOutputToAFile:
    def test_the_file_holds_the_document(self, runner, source, tmp_path):
        out = tmp_path / "meta.json"
        result = runner.invoke(metadata_command, [str(source), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert json.loads(out.read_text(encoding="utf-8"))["ms_analysis"]

    def test_stdout_stays_clean_for_the_shell(self, runner, source, tmp_path):
        out = tmp_path / "meta.json"
        result = runner.invoke(metadata_command, [str(source), "-o", str(out)])
        assert result.stdout == ""


class TestMerge:
    def test_user_fields_are_overlaid(self, runner, source, tmp_path):
        overlay = tmp_path / "user.json"
        overlay.write_text(
            json.dumps({"sample": {"organism": "Mus musculus"}}), encoding="utf-8"
        )
        result = runner.invoke(metadata_command, [str(source), "--merge", str(overlay)])
        assert result.exit_code == 0, result.output
        assert json.loads(result.stdout)["sample"]["organism"] == "Mus musculus"

    def test_a_merge_that_breaks_the_schema_is_reported(self, runner, source, tmp_path):
        overlay = tmp_path / "user.json"
        overlay.write_text(
            json.dumps({"ms_analysis": {"polarity": "sideways"}}), encoding="utf-8"
        )
        result = runner.invoke(metadata_command, [str(source), "--merge", str(overlay)])
        assert result.exit_code == 1
        assert "ms_analysis.polarity" in result.stderr
        # The document is still written: what it says is still what the
        # file and the overlay say.
        assert json.loads(result.stdout)["ms_analysis"]["polarity"] == "sideways"


class TestUsage:
    def test_a_missing_input_is_a_usage_error(self, runner, tmp_path):
        result = runner.invoke(metadata_command, [str(tmp_path / "nope.d")])
        assert result.exit_code == 2

    def test_an_unreadable_input_names_the_path(self, runner, tmp_path):
        empty = tmp_path / "empty.d"
        empty.mkdir()
        result = runner.invoke(metadata_command, [str(empty)])
        assert result.exit_code == 1
        assert "empty.d" in result.stderr


class TestTheSubcommandIsDispatched:
    def test_the_first_argument_picks_it_off(self, monkeypatch, capsys):
        from thyra.__main__ import cli

        monkeypatch.setattr("sys.argv", ["thyra", "metadata", "--help"], raising=False)
        with pytest.raises(SystemExit) as excinfo:
            cli()
        assert excinfo.value.code == 0
        assert "without converting it" in capsys.readouterr().out

    def test_the_conversion_help_names_it(self, monkeypatch, capsys):
        from thyra.__main__ import cli

        monkeypatch.setattr("sys.argv", ["thyra", "--help"], raising=False)
        with pytest.raises(SystemExit):
            cli()
        assert "thyra metadata" in capsys.readouterr().out
