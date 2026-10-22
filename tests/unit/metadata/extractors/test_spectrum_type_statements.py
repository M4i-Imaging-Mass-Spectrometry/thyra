"""The spectrum type is the one every converted spectrum states (D81).

On imzML a spectrum states it through its own terms and the parameter groups
it references; a group no spectrum references counts for nothing. On mzPeak
the per-spectrum ``spectrum_representation`` column of the spectra kept
decides, and an archive read from its peaks member stays centroid. Spectra
stating both types record none. The file-level declaration counts only when
no spectrum states a type, and only if it names one.
"""

import logging
from pathlib import Path

import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

import thyra
from tests.fixtures.mzpeak_builder import Spectrum, build_mzpeak, grid_spectra
from thyra.metadata.constants import SpectrumType, spectrum_type_from_statements
from thyra.readers.imzml.imzml_reader import ImzMLReader
from thyra.readers.mzpeak import MzPeakReader

CENTROID = SpectrumType.CENTROID
PROFILE = SpectrumType.PROFILE

PROFILE_TERM = (
    '<cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>'
)
CENTROID_TERM = (
    '<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>'
)
#: A group no spectrum references.
UNUSED_CENTROID_GROUP = (
    f'<referenceableParamGroup id="peaklist">{CENTROID_TERM}'
    "</referenceableParamGroup>"
)
MS2_TERM = '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>'


def _write_imzml(
    path: Path,
    file_content=(PROFILE_TERM,),
    group_term=PROFILE_TERM,
    spectrum_terms=None,
    unused_group=None,
    positions=None,
):
    """A 3 x 2 processed imzML of profile spectra, rewritten as asked.

    ``file_content`` replaces the representation term in ``fileContent``.
    ``group_term`` replaces it in ``spectrum1``, the group every spectrum
    references (``""`` removes it). ``spectrum_terms`` adds one string per
    spectrum on the spectrum itself. ``unused_group`` puts a group no
    spectrum references ``"first"`` or ``"last"`` in the group list.
    """
    positions = positions or [(x, y, 1) for x in (1, 2, 3) for y in (1, 2)]
    with ImzMLWriter(str(path), mode="processed", spec_type="profile") as writer:
        for position in positions:
            writer.addSpectrum([100.0, 150.0, 200.0], [1.0, 2.0, 3.0], position)
    text = path.read_text()
    head, tail = text.split("</fileContent>", 1)
    head = head.replace(PROFILE_TERM, "".join(file_content), 1)
    tail = tail.replace(PROFILE_TERM, group_term, 1)
    if unused_group == "first":
        tail = tail.replace(
            '<referenceableParamGroup id="spectrum1">',
            UNUSED_CENTROID_GROUP + '<referenceableParamGroup id="spectrum1">',
            1,
        )
    elif unused_group == "last":
        tail = tail.replace(
            "</referenceableParamGroupList>",
            UNUSED_CENTROID_GROUP + "</referenceableParamGroupList>",
            1,
        )
    text = head + "</fileContent>" + tail
    if spectrum_terms is not None:
        parts = text.split("<scanList")
        assert len(parts) == len(spectrum_terms) + 1
        text = parts[0] + "".join(
            (term or "") + "<scanList" + part
            for term, part in zip(spectrum_terms, parts[1:])
        )
    path.write_text(text)
    return path


def _conversion(path: Path):
    """The type the conversion route records, and whether it was guessed."""
    with ImzMLReader(path) as reader:
        essential = reader.get_essential_metadata()
    return essential.spectrum_type, essential.spectrum_type_guessed


def _document(path: Path):
    """The type the metadata-only route records: the head and first spectrum."""
    stored = thyra.read_metadata_document(path)["ms_analysis"].get(
        "spectrum_representation"
    )
    # The document stores the vocabulary word, not the CV name.
    return {"centroid": CENTROID, "profile": PROFILE}.get(stored, stored)


def _warnings(records):
    return [r.getMessage() for r in records if r.levelno >= logging.WARNING]


class TestTheRule:
    @pytest.mark.parametrize(
        "centroid, profile, declared, expected, stated",
        [
            (0, 9, [], PROFILE, True),
            (4, 0, [PROFILE], CENTROID, True),
            (4, 0, [CENTROID, PROFILE], CENTROID, True),
            (8, 1, [], None, True),
            (0, 0, [CENTROID], CENTROID, True),
            (0, 0, [CENTROID, PROFILE], None, True),
            (0, 0, [], None, False),
        ],
    )
    def test_one_rule(self, centroid, profile, declared, expected, stated):
        result = spectrum_type_from_statements(centroid, profile, declared, "f")

        assert result.spectrum_type == expected
        assert result.stated is stated

    def test_mixed_spectra_warn_with_the_counts(self):
        result = spectrum_type_from_statements(8, 1, [], "f")

        assert result.warning is not None
        assert "8 centroid, 1 profile" in result.warning

    def test_spectra_against_the_declaration_warn_with_both(self):
        result = spectrum_type_from_statements(0, 6, [CENTROID], "f")

        assert result.warning is not None
        assert CENTROID in result.warning and PROFILE in result.warning


class TestImzML:
    @pytest.mark.parametrize("where", ["first", "last"])
    def test_an_unused_group_counts_for_nothing(self, tmp_path, where):
        """F176: a centroid group no spectrum references, fileContent silent."""
        path = _write_imzml(
            tmp_path / "unused.imzML", file_content=(), unused_group=where
        )

        assert _conversion(path) == (PROFILE, False)
        assert _document(path) == PROFILE

    @pytest.mark.parametrize(
        "declared", [(CENTROID_TERM, PROFILE_TERM), (PROFILE_TERM, CENTROID_TERM)]
    )
    def test_a_declaration_of_both_yields_to_the_spectra(self, tmp_path, declared):
        """F178: the order fileContent lists the two terms in says nothing."""
        path = _write_imzml(tmp_path / "both.imzML", file_content=declared)

        assert _conversion(path) == (PROFILE, False)
        assert _document(path) == PROFILE

    def test_spectra_win_over_a_declaration_of_the_other(self, tmp_path, thyra_logs):
        path = _write_imzml(tmp_path / "other.imzML", file_content=(CENTROID_TERM,))

        with thyra_logs("thyra.metadata", logging.WARNING) as records:
            assert _conversion(path) == (PROFILE, False)

        assert any(
            "declares 'centroid spectrum'" in message for message in _warnings(records)
        )

    def test_spectra_of_both_types_record_none(self, tmp_path, thyra_logs):
        path = _write_imzml(
            tmp_path / "mixed.imzML",
            group_term="",
            spectrum_terms=[PROFILE_TERM] * 5 + [CENTROID_TERM],
        )

        with thyra_logs("thyra.metadata", logging.WARNING) as records:
            assert _conversion(path) == (None, False)

        assert any("1 centroid, 5 profile" in m for m in _warnings(records))

    @pytest.mark.parametrize("declared, expected", [(PROFILE_TERM, PROFILE)])
    def test_silent_spectra_take_a_declaration_of_one(
        self, tmp_path, declared, expected
    ):
        path = _write_imzml(
            tmp_path / "silent.imzML", file_content=(declared,), group_term=""
        )

        assert _conversion(path) == (expected, False)
        assert _document(path) == expected

    def test_silent_spectra_and_a_declaration_of_both_record_none(self, tmp_path):
        """Not the processed-mode guess: the file did state something."""
        path = _write_imzml(
            tmp_path / "silent_both.imzML",
            file_content=(PROFILE_TERM, CENTROID_TERM),
            group_term="",
        )

        assert _conversion(path) == (None, False)
        assert _document(path) is None

    def test_an_msn_spectrum_left_out_takes_no_part(self, tmp_path):
        """Its own centroid term is not a converted spectrum's statement."""
        positions = [(1, 1, 1), (1, 1, 1), (2, 1, 1), (3, 1, 1)]
        path = _write_imzml(
            tmp_path / "msn.imzML",
            positions=positions,
            spectrum_terms=[None, MS2_TERM + CENTROID_TERM, None, None],
        )

        assert _conversion(path) == (PROFILE, False)

    def test_an_override_is_compared_with_the_spectra(self, tmp_path, thyra_logs):
        """fileContent says centroid, the spectra profile: centroid contradicts."""
        path = _write_imzml(tmp_path / "override.imzML", file_content=(CENTROID_TERM,))

        with thyra_logs("thyra.metadata", logging.WARNING) as records:
            with ImzMLReader(path, spectrum_type="centroid") as reader:
                assert reader.get_essential_metadata().spectrum_type == CENTROID

        assert any("CONTRADICTS" in m for m in _warnings(records))


def _archive_type(path: Path):
    with MzPeakReader(path) as reader:
        return reader.get_essential_metadata().spectrum_type


def _contents(*accessions):
    names = {"MS:1000127": CENTROID, "MS:1000128": PROFILE}
    return [
        {"name": names[a], "accession": a, "value": None, "unit": None}
        for a in accessions
    ]


class TestMzPeak:
    def test_one_centroid_label_among_profile_records_none(self, tmp_path, thyra_logs):
        """F179: the column decides by every converted spectrum, not any one."""
        archive = build_mzpeak(
            tmp_path / "mixed.mzpeak",
            grid_spectra(3, 3),
            file_contents=[],
            spectrum_representation=["MS:1000128"] * 8 + ["MS:1000127"],
        )

        with thyra_logs("thyra.metadata", logging.WARNING) as records:
            assert _archive_type(archive) is None

        assert any("1 centroid, 8 profile" in m for m in _warnings(records))

    def test_an_ms2_spectrum_left_out_takes_no_part(self, tmp_path):
        spectra = grid_spectra(3, 3) + [Spectrum(1, 1, [200.0, 300.0], [5.0, 6.0])]
        archive = build_mzpeak(
            tmp_path / "ms2.mzpeak",
            spectra,
            file_contents=[],
            spectrum_representation=["MS:1000128"] * 9 + ["MS:1000127"],
            ms_levels=[1] * 9 + [2],
        )

        assert _archive_type(archive) == PROFILE

    @pytest.mark.parametrize(
        "contents",
        [("MS:1000127", "MS:1000128"), ("MS:1000128", "MS:1000127")],
    )
    def test_contents_naming_both_yield_to_the_spectra(self, tmp_path, contents):
        """F178: the order of file_description.contents says nothing."""
        archive = build_mzpeak(
            tmp_path / "both.mzpeak",
            grid_spectra(2, 2),
            file_contents=_contents(*contents),
            spectrum_representation="MS:1000128",
        )

        assert _archive_type(archive) == PROFILE

    def test_contents_naming_both_over_silent_spectra_record_none(self, tmp_path):
        archive = build_mzpeak(
            tmp_path / "silent.mzpeak",
            grid_spectra(2, 2),
            file_contents=_contents("MS:1000127", "MS:1000128"),
            spectrum_representation="",
        )

        assert _archive_type(archive) is None

    @pytest.mark.parametrize("column", ["", "profile spectrum", "MS:1000128"])
    @pytest.mark.parametrize("declared", [(), ("MS:1000128",)])
    def test_an_archive_read_from_its_peaks_member_is_centroid(
        self, tmp_path, column, declared
    ):
        """F76: the peaks member holds centroid spectra by definition."""
        archive = build_mzpeak(
            tmp_path / "peaks.mzpeak",
            grid_spectra(3, 2),
            signal="centroid",
            empty_peer=True,
            spectrum_representation=column,
            file_contents=_contents(*declared),
        )

        assert _archive_type(archive) == CENTROID
