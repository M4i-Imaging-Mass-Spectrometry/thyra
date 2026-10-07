"""Polarity is the value every statement in a source agrees on (D35).

An imzML can state it in ``fileContent``, in a referenceable parameter group
and on each spectrum; an mzPeak archive in ``file_description.contents`` and
in each spectrum's ``scan_polarity``. Every place is read, on the conversion
route and on the metadata-only route, and a source that states both values
records neither.
"""

import logging
from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

import thyra
from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from thyra.metadata.constants import agreed_polarity
from thyra.readers.imzml.imzml_reader import ImzMLReader

POSITIVE_TERM = (
    '<cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>'
)
NEGATIVE_TERM = (
    '<cvParam cvRef="MS" accession="MS:1000129" name="negative scan" value=""/>'
)


def _write_imzml(
    path: Path,
    polarity=None,
    file_content_term=None,
    spectrum_terms=None,
) -> Path:
    """A 3 x 2 processed imzML.

    ``polarity`` goes where pyimzml's writer puts it, a parameter group
    every spectrum references. ``file_content_term`` is added to
    ``fileContent``. ``spectrum_terms`` gives one term per spectrum, in
    order, written on the spectrum itself (``None`` writes none).
    """
    kwargs = {} if polarity is None else {"polarity": polarity}
    with ImzMLWriter(str(path), mode="processed", **kwargs) as writer:
        for x in (1, 2, 3):
            for y in (1, 2):
                writer.addSpectrum([100.0, 150.0, 200.0], [1.0, 2.0, 3.0], (x, y, 1))
    text = path.read_text()
    if file_content_term is not None:
        text = text.replace("</fileContent>", file_content_term + "</fileContent>", 1)
    if spectrum_terms is not None:
        parts = text.split("<scanList")
        assert len(parts) == len(spectrum_terms) + 1
        text = parts[0] + "".join(
            ("" if term is None else term) + "<scanList" + part
            for term, part in zip(spectrum_terms, parts[1:])
        )
    path.write_text(text)
    return path


def _conversion_polarity(path: Path):
    """What the conversion route records: every spectrum read."""
    with ImzMLReader(path) as reader:
        return reader.get_comprehensive_metadata().acquisition_params["polarity"]


def _document_polarity(path: Path):
    """What the metadata-only route records: the head and the first spectrum."""
    return thyra.read_metadata_document(path)["ms_analysis"].get("polarity")


class TestAgreedPolarity:
    @pytest.mark.parametrize(
        "stated, expected",
        [
            ([], None),
            (["positive"], "positive"),
            (["negative", "negative"], "negative"),
            (["positive", "negative"], None),
        ],
    )
    def test_one_value_or_none(self, stated, expected):
        assert agreed_polarity(stated) == expected


class TestImzML:
    @pytest.mark.parametrize("polarity", ["positive", "negative"])
    def test_a_parameter_group_states_it(self, tmp_path, polarity):
        """Where pyimzml's own writer puts it."""
        path = _write_imzml(tmp_path / "group.imzML", polarity=polarity)

        assert _conversion_polarity(path) == polarity
        assert _document_polarity(path) == polarity

    def test_every_spectrum_states_it(self, tmp_path):
        path = _write_imzml(
            tmp_path / "spectra.imzML", spectrum_terms=[POSITIVE_TERM] * 6
        )

        assert _conversion_polarity(path) == "positive"
        assert _document_polarity(path) == "positive"

    def test_file_content_states_it(self, tmp_path):
        path = _write_imzml(tmp_path / "fc.imzML", file_content_term=NEGATIVE_TERM)

        assert _conversion_polarity(path) == "negative"
        assert _document_polarity(path) == "negative"

    def test_file_content_and_spectra_disagree(self, tmp_path, thyra_logs):
        """The file-level term does not win over the spectra."""
        path = _write_imzml(
            tmp_path / "contradicts.imzML",
            polarity="negative",
            file_content_term=POSITIVE_TERM,
        )

        with thyra_logs("thyra.metadata", logging.WARNING) as records:
            assert _conversion_polarity(path) is None
        assert "both positive and negative" in records.text
        assert _document_polarity(path) is None

    def test_a_later_spectrum_that_disagrees_is_seen_on_conversion(self, tmp_path):
        """The conversion reads every spectrum, not just the first."""
        terms = [POSITIVE_TERM] * 5 + [NEGATIVE_TERM]
        path = _write_imzml(tmp_path / "alternates.imzML", spectrum_terms=terms)

        assert _conversion_polarity(path) is None

    def test_nothing_stated_records_nothing(self, tmp_path):
        path = _write_imzml(tmp_path / "silent.imzML")

        assert _conversion_polarity(path) is None
        assert _document_polarity(path) is None

    def test_the_store_carries_it(self, tmp_path):
        path = _write_imzml(
            tmp_path / "store.imzML", spectrum_terms=[NEGATIVE_TERM] * 6
        )
        out = tmp_path / "store.zarr"

        assert thyra.convert_msi(path, out, pixel_size_um=10.0)

        from thyra.metadata.schema.store_io import read_msi_metadata_blocks

        block = next(iter(read_msi_metadata_blocks(out).values()))["ms_analysis"]
        assert block["polarity"] == "negative"
        assert block["polarity_term"]["accession"] == "MS:1000129"


class TestMzPeak:
    def _polarity(self, path: Path):
        return thyra.read_metadata_document(path)["ms_analysis"].get("polarity")

    def test_the_file_level_term_of_an_imzml_archive(self, tmp_path):
        contents = [
            {"name": "profile spectrum", "accession": "MS:1000128"},
            {"name": "positive scan", "accession": "MS:1000130"},
        ]
        path = build_mzpeak(
            tmp_path / "fc.mzpeak", grid_spectra(2, 2), file_contents=contents
        )

        assert self._polarity(path) == "positive"

    def test_scan_polarity_alone_as_a_bruker_archive_states_it(self, tmp_path):
        path = build_mzpeak(
            tmp_path / "scans.mzpeak", grid_spectra(2, 2), scan_polarity=[-1] * 4
        )

        assert self._polarity(path) == "negative"

    def test_spectra_that_disagree_record_nothing(self, tmp_path):
        path = build_mzpeak(
            tmp_path / "mixed.mzpeak", grid_spectra(2, 2), scan_polarity=[1, 1, -1, 1]
        )

        assert self._polarity(path) is None

    def test_null_and_unknown_values_state_nothing(self, tmp_path):
        path = build_mzpeak(
            tmp_path / "unknown.mzpeak",
            grid_spectra(2, 2),
            scan_polarity=[None, 0, 1, None],
        )

        assert self._polarity(path) == "positive"

    def test_the_file_term_and_the_spectra_disagree(self, tmp_path):
        contents = [{"name": "positive scan", "accession": "MS:1000130"}]
        path = build_mzpeak(
            tmp_path / "contradicts.mzpeak",
            grid_spectra(2, 2),
            file_contents=contents,
            scan_polarity=[-1] * 4,
        )

        assert self._polarity(path) is None

    def test_an_archive_that_states_nothing(self, tmp_path):
        path = build_mzpeak(tmp_path / "silent.mzpeak", grid_spectra(2, 2))

        assert self._polarity(path) is None


def test_arrays_unchanged_by_the_extra_parse_fields(tmp_path):
    """Collecting the terms does not change what the spectra read as."""
    path = _write_imzml(tmp_path / "arrays.imzML", polarity="positive")
    with ImzMLReader(path) as reader:
        spectra = list(reader.iter_spectra())
    assert len(spectra) == 6
    np.testing.assert_array_equal(spectra[0][1], [100.0, 150.0, 200.0])
