import json
import sys
from unittest.mock import patch

import pytest

from thyra.tools.check_ontology import main


class TestCheckOntology:
    """Tests for the ontology checking tool."""

    @patch("thyra.tools.check_ontology.ONTOLOGY", autospec=True)
    @patch("thyra.metadata.validator.ImzMLOntologyValidator", autospec=True)
    def test_main_file_check(
        self, mock_validator_class, mock_ontology, capsys, tmp_path, monkeypatch
    ):
        """Verify the CLI output for a single file check with no unknown terms."""
        # Setup mocks
        mock_validator_instance = mock_validator_class.return_value
        # The summary the validator computed is what gets printed; it used
        # to be discarded for the literal "Test summary for file" (#299).
        mock_validator_instance.validate_file.return_value = {
            "unknown_list": [],
            "summary": (
                "Ontology Validation Summary\n"
                "===========================\n"
                "Total CV terms found: 3"
            ),
        }
        mock_ontology.report_unknown_terms.return_value = "Global unknown terms: None"

        # Create a dummy imzML file
        dummy_imzml = tmp_path / "test.imzML"
        dummy_imzml.write_text(
            '<mzML><cvList><cv id="MS" fullName="MS" version="1.0" '
            'uri="http://example.com/ms.obo"/></cvList><run><spectrumList '
            'count="0"/></run></mzML>'
        )

        # Set up command line arguments
        monkeypatch.setattr(sys, "argv", ["check_ontology", str(dummy_imzml)])

        # Run main
        main()

        # Assertions
        mock_validator_class.assert_called_once()  # Ensure the validator class was instantiated
        mock_validator_instance.validate_file.assert_called_once_with(dummy_imzml)
        captured = capsys.readouterr()
        # The summary comes from the validator's return value.
        expected_output = (
            "Ontology Validation Summary\n"
            "===========================\n"
            "Total CV terms found: 3\n\n"
            "No unknown terms encountered.\n\n"
            "Global unknown terms: None\n"
        )
        assert expected_output == captured.out

    @patch("thyra.tools.check_ontology.ONTOLOGY", autospec=True)
    @patch("thyra.metadata.validator.ImzMLOntologyValidator", autospec=True)
    def test_main_file_check_with_unknown_terms(
        self, mock_validator_class, mock_ontology, capsys, tmp_path, monkeypatch
    ):
        """A file with unknown terms lists them rather than crashing.

        Regression test: main() used to read `unknown_terms` -- which
        validate_file reports as a COUNT -- where it needed the
        `unknown_list` entries, so the first real file with any unknown
        term died on len() of an int.  The earlier mocks fed a list
        under the count's key, so the suite agreed with the bug.
        """
        mock_validator_instance = mock_validator_class.return_value
        mock_validator_instance.validate_file.return_value = {
            "unknown_terms": 1,
            "unknown_list": [
                {"accession": "MS:9999999", "name": "mystery", "value": ""}
            ],
        }
        mock_ontology.report_unknown_terms.return_value = "Global unknown terms: Some"

        dummy_imzml = tmp_path / "test.imzML"
        dummy_imzml.write_text("<mzML/>")

        monkeypatch.setattr(sys, "argv", ["check_ontology", str(dummy_imzml)])
        main()
        captured = capsys.readouterr()
        assert "Found 1 unknown terms:" in captured.out
        assert "  - MS:9999999: mystery" in captured.out

    @patch("thyra.tools.check_ontology.ONTOLOGY", autospec=True)
    @patch("thyra.metadata.validator.ImzMLOntologyValidator", autospec=True)
    def test_main_directory_check(
        self, mock_validator_class, mock_ontology, capsys, tmp_path, monkeypatch
    ):
        """Verify the CLI output for a directory check."""
        # Setup mocks
        mock_validator_instance = mock_validator_class.return_value
        mock_validator_instance.validate_directory.return_value = {
            "files_checked": 2,
            "all_unknown_terms": ["term1", "term2"],
            "summary": "Directory summary",
        }
        # This mock is for the global report at the end
        mock_ontology.report_unknown_terms.return_value = "Global unknown terms: Some"

        # Create a dummy directory and files
        dummy_dir = tmp_path / "test_dir"
        dummy_dir.mkdir()
        (dummy_dir / "file1.imzML").write_text("<mzML/>")
        (dummy_dir / "file2.imzML").write_text("<mzML/>")

        # Set up command line arguments
        monkeypatch.setattr(sys, "argv", ["check_ontology", str(dummy_dir)])

        # Run main
        main()

        # Assertions
        mock_validator_class.assert_called_once()
        mock_validator_instance.validate_directory.assert_called_once_with(dummy_dir)
        captured = capsys.readouterr()

        # CORRECTED: Define the expected output precisely.
        expected_output = (
            "Checked 2 files\n"
            "Found 2 unique unknown terms\n\n"
            "Most common unknown terms:\n"
            "  - term1\n"
            "  - term2\n\n"
            "Global unknown terms: Some\n"
        )

        assert expected_output == captured.out

    @patch("thyra.tools.check_ontology.ONTOLOGY", autospec=True)
    @patch("thyra.metadata.validator.ImzMLOntologyValidator", autospec=True)
    def test_main_output_to_json(
        self, mock_validator_class, mock_ontology, capsys, tmp_path, monkeypatch
    ):
        """Verify that results are correctly saved to a JSON file."""
        # Setup mocks
        mock_validator_instance = mock_validator_class.return_value
        validation_result = {
            "summary": "Test summary for JSON",
            "unknown_terms": 1,
            "unknown_list": [
                {"accession": "MS:9999999", "name": "json_term", "value": ""}
            ],
        }
        mock_validator_instance.validate_file.return_value = validation_result
        mock_ontology.report_unknown_terms.return_value = "Global unknown terms: JSON"

        output_json = tmp_path / "output.json"
        dummy_imzml = tmp_path / "test_json.imzML"
        dummy_imzml.write_text("""<mzML/>""")

        # Set up command line arguments
        monkeypatch.setattr(
            sys,
            "argv",
            ["check_ontology", str(dummy_imzml), "--output", str(output_json)],
        )

        # Run main
        main()

        # Assertions
        captured = capsys.readouterr()

        # CORRECTED: Removed the extra newline for a more likely output format
        expected_output = (
            f"Results saved to {output_json}\nGlobal unknown terms: JSON\n"
        )

        assert expected_output == captured.out

        assert output_json.exists()
        with open(output_json, "r") as f:
            content = json.load(f)
        assert content == validation_result

    @patch("thyra.tools.check_ontology.ONTOLOGY", autospec=True)
    @patch("thyra.metadata.validator.ImzMLOntologyValidator", autospec=True)
    def test_main_verbose_output(
        self, mock_validator_class, mock_ontology, capsys, tmp_path, monkeypatch
    ):
        """Verify the verbose output includes the validator's summary."""
        # Setup mocks
        mock_validator_instance = mock_validator_class.return_value
        mock_validator_instance.validate_file.return_value = {
            "summary": "Verbose summary",
            "unknown_list": [],
        }
        mock_ontology.report_unknown_terms.return_value = (
            "Global unknown terms: Verbose"
        )

        # Create a dummy imzML file
        dummy_imzml = tmp_path / "test_verbose.imzML"
        dummy_imzml.write_text("""<mzML/>""")

        # Set up command line arguments
        monkeypatch.setattr(
            sys, "argv", ["check_ontology", str(dummy_imzml), "--verbose"]
        )

        # Run main
        main()

        # Assertions
        captured = capsys.readouterr()

        # The summary is followed by a blank line on both paths now. It
        # used to be printed without one under --verbose and with one
        # otherwise, because the two were separate branches (#299).
        expected_output = (
            "Verbose summary\n\n"
            "No unknown terms encountered.\n\n"
            "Global unknown terms: Verbose\n"
        )

        assert expected_output == captured.out


class TestTheSummaryReachesTheUser:
    """Issue #299: the tool computed a summary and printed a placeholder.

    ``_print_file_results`` printed a banner and the literal string
    "Test summary for file" for every run without ``--verbose``, so the
    term counts, the known/unknown split and the most common terms were
    all discarded after being calculated. These run the real validator
    against the real ontology, so they also cover the term-name lookup
    that the summary is the only consumer of.
    """

    THREE_TERMS = (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<mzML xmlns="http://psi.hupo.org/ms/mzml">\n'
        "  <fileDescription><fileContent>\n"
        '    <cvParam cvRef="IMS" accession="IMS:1000080" '
        'name="universally unique identifier" value="{0}"/>\n'
        '    <cvParam cvRef="IMS" accession="IMS:1000042" '
        'name="max count of pixels x" value="3"/>\n'
        '    <cvParam cvRef="IMS" accession="IMS:9999999" '
        'name="not a real term" value="x"/>\n'
        "  </fileContent></fileDescription>\n"
        "</mzML>\n"
    )

    def _run(self, tmp_path, capsys, monkeypatch, *extra_args):
        """``monkeypatch``, not a try/finally.

        A test that assigns ``sys.argv`` and restores it by hand leaks the
        value for the rest of the process whenever it fails before the
        finally -- which is what #292 removed from the rest of this file.
        """
        imzml = tmp_path / "three_terms.imzML"
        imzml.write_text(self.THREE_TERMS)
        monkeypatch.setattr(sys, "argv", ["check_ontology", str(imzml), *extra_args])
        main()
        return capsys.readouterr().out

    def test_the_placeholder_is_gone(self, tmp_path, capsys, monkeypatch):
        assert "Test summary for file" not in self._run(tmp_path, capsys, monkeypatch)

    def test_the_term_count_is_the_one_the_validator_computed(
        self, tmp_path, capsys, monkeypatch
    ):
        assert "Total CV terms found: 3" in self._run(tmp_path, capsys, monkeypatch)

    def test_the_known_unknown_split_is_reported(self, tmp_path, capsys, monkeypatch):
        out = self._run(tmp_path, capsys, monkeypatch)
        assert "Known terms: 2" in out
        assert "Unknown terms: 1" in out

    def test_it_does_not_need_verbose(self, tmp_path, capsys, monkeypatch):
        """``--verbose`` turns on INFO logging; it is not a request for
        the tool's only output."""
        plain = self._run(tmp_path, capsys, monkeypatch)
        assert "Total CV terms found: 3" in plain

    def test_the_most_common_terms_carry_their_names_not_their_xsd_types(
        self, tmp_path, capsys, monkeypatch
    ):
        """``ONTOLOGY.get_term`` returns ``(name, data_type)``.

        The summary read ``[1]``, so it printed "xsd:string" for
        IMS:1000080 and ``None`` for a term with no declared type, where
        the term's name belongs.
        """
        out = self._run(tmp_path, capsys, monkeypatch)
        assert "IMS:1000080 (universally unique identifier)" in out
        assert "xsd:string" not in out

    def test_the_unknown_term_is_named_inside_the_summary(
        self, tmp_path, capsys, monkeypatch
    ):
        """Within the summary block, not merely somewhere in the output.

        The unknown-terms report printed the accession anyway, so
        asserting it appears at all passes without the summary ever being
        reached.
        """
        out = self._run(tmp_path, capsys, monkeypatch)
        summary = out.split("Found 1 unknown terms:")[0]
        assert "Unknown Terms:" in summary
        assert "IMS:9999999" in summary


def _imzml(*terms):
    """A minimal imzML document holding ``(accession, times)`` cvParams."""
    params = "".join(
        f'<cvParam cvRef="IMS" accession="{accession}" name="term" value="x"/>' * times
        for accession, times in terms
    )
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<mzML xmlns="http://psi.hupo.org/ms/mzml"><fileDescription>'
        f"<fileContent>{params}</fileContent></fileDescription></mzML>\n"
    ).encode()


GOOD = _imzml(("IMS:1000042", 1), ("IMS:9999999", 1))


class TestBadInput:
    """Issue #529: a missing path or one unreadable file."""

    def test_a_path_that_does_not_exist_exits_2(self, tmp_path, capsys):
        """It used to be searched as an empty folder: "Checked 0 files", 0."""
        report = tmp_path / "report.json"
        code = main([str(tmp_path / "smaple.imzML"), "--output", str(report)])

        captured = capsys.readouterr()
        assert code == 2
        assert captured.err.strip().endswith("does not exist.")
        assert len(captured.err.strip().splitlines()) == 1
        assert captured.out == ""
        assert not report.exists()

    @pytest.mark.parametrize(
        "bad",
        [
            pytest.param(b"", id="empty"),
            pytest.param(GOOD[: len(GOOD) // 2], id="cut-short"),
            pytest.param(GOOD.replace(b'value="x"', b'value="5 \xb5m"'), id="latin-1"),
        ],
    )
    def test_one_unreadable_file_alone_is_one_line_and_exit_1(
        self, tmp_path, capsys, bad
    ):
        path = tmp_path / "bad.imzML"
        path.write_bytes(bad)

        code = main([str(path)])

        captured = capsys.readouterr()
        assert code == 1
        assert "Traceback" not in captured.err
        assert captured.err.startswith(f"Error: cannot read {path}:")
        assert len(captured.err.strip().splitlines()) == 1

    @pytest.mark.parametrize(
        "bad",
        [
            pytest.param(b"", id="empty"),
            pytest.param(GOOD[: len(GOOD) // 2], id="cut-short"),
            pytest.param(GOOD.replace(b'value="x"', b'value="5 \xb5m"'), id="latin-1"),
        ],
    )
    def test_in_a_folder_the_bad_file_is_reported_and_the_rest_checked(
        self, tmp_path, capsys, bad
    ):
        (tmp_path / "a_good.imzML").write_bytes(GOOD)
        (tmp_path / "b_good.imzML").write_bytes(GOOD)
        (tmp_path / "z_bad.imzML").write_bytes(bad)
        report = tmp_path / "report.json"

        code = main([str(tmp_path), "--output", str(report)])

        out = capsys.readouterr().out
        assert code == 1
        assert "Checked 2 files" in out
        assert "Could not read 1 files:" in out
        assert "z_bad.imzML" in out
        saved = json.loads(report.read_text(encoding="utf-8"))
        assert len(saved["per_file_results"]) == 2
        assert list(saved["unreadable_files"]) == [str(tmp_path / "z_bad.imzML")]

    def test_the_suffix_is_matched_in_any_letter_case(self, tmp_path, capsys):
        """``sample.imzml`` was checked on Windows and missed on Linux."""
        (tmp_path / "lower.imzml").write_bytes(GOOD)
        (tmp_path / "upper.IMZML").write_bytes(GOOD)

        assert main([str(tmp_path)]) == 0
        assert "Checked 2 files" in capsys.readouterr().out

    def test_a_folder_named_like_an_imzml_is_not_a_file(self, tmp_path, capsys):
        (tmp_path / "a_good.imzML").write_bytes(GOOD)
        (tmp_path / "z_bad.imzML").mkdir()

        code = main([str(tmp_path)])

        assert code == 0
        assert "Checked 1 files" in capsys.readouterr().out


class TestUnknownTermsAreAListInAFixedOrder:
    """Issue #529: ``all_unknown_terms`` was saved as the text of a Python
    set, and printed in hash order."""

    def _folder(self, tmp_path):
        tmp_path.mkdir(exist_ok=True)
        (tmp_path / "a.imzML").write_bytes(_imzml(("IMS:9999997", 1)))
        (tmp_path / "b.imzML").write_bytes(_imzml(("IMS:9999998", 3)))
        (tmp_path / "c.imzML").write_bytes(_imzml(("IMS:9999999", 1)))
        return tmp_path

    def test_the_json_holds_a_list_most_often_seen_first(self, tmp_path):
        report = tmp_path / "report.json"
        main([str(self._folder(tmp_path / "in")), "--output", str(report)])

        saved = json.loads(report.read_text(encoding="utf-8"))
        assert saved["all_unknown_terms"] == [
            "IMS:9999998",
            "IMS:9999997",
            "IMS:9999999",
        ]
        assert saved["unknown_term_counts"] == {
            "IMS:9999998": 3,
            "IMS:9999997": 1,
            "IMS:9999999": 1,
        }

    def test_the_printed_list_has_the_same_order(self, tmp_path, capsys):
        main([str(self._folder(tmp_path))])

        listed = capsys.readouterr().out.split("Most common unknown terms:")[1]
        lines = [line for line in listed.splitlines() if line.startswith("  - ")]
        assert lines[:3] == [
            "  - IMS:9999998 (seen 3 times)",
            "  - IMS:9999997 (seen 1 times)",
            "  - IMS:9999999 (seen 1 times)",
        ]
