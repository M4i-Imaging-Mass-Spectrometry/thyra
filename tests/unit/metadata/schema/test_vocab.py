"""Vocabulary normalisation: vendor spellings land on PSI-MS terms."""

import pytest

from thyra.metadata.ontology.cache import ONTOLOGY
from thyra.metadata.schema.vocab import (
    ANALYZER_ACCESSIONS,
    IONISATION_SOURCE_ACCESSIONS,
    POLARITY_ACCESSIONS,
    normalize_analyzer,
    normalize_ionisation_source,
    normalize_polarity,
    term_from_accession,
)


class TestTermFromAccession:
    def test_name_comes_from_the_local_table(self):
        term = term_from_accession("MS:1000075")
        assert term.name == ONTOLOGY.terms["MS:1000075"][0]

    def test_unknown_accession_raises(self):
        with pytest.raises(KeyError):
            term_from_accession("MS:9999999")

    @pytest.mark.parametrize(
        "accession",
        sorted(
            set(POLARITY_ACCESSIONS.values())
            | set(IONISATION_SOURCE_ACCESSIONS.values())
            | set(ANALYZER_ACCESSIONS.values())
        ),
    )
    def test_every_vocabulary_accession_resolves_locally(self, accession):
        # The vocabularies may only reference accessions the shipped
        # ontology tables actually contain, or written blocks would
        # fail their own validation.
        assert term_from_accession(accession).name


class TestNormalizePolarity:
    @pytest.mark.parametrize(
        "raw", ["positive", "Positive", "POS", "+", "positive scan"]
    )
    def test_positive_spellings(self, raw):
        result = normalize_polarity(raw)
        assert result is not None
        canonical, term = result
        assert canonical == "positive"
        assert term.accession == "MS:1000130"

    @pytest.mark.parametrize("raw", ["negative", "NEG", "-", "negative ion mode"])
    def test_negative_spellings(self, raw):
        result = normalize_polarity(raw)
        assert result is not None
        assert result[0] == "negative"
        assert result[1].accession == "MS:1000129"

    @pytest.mark.parametrize("raw", [None, "", "  ", "alternating", 42, "both"])
    def test_unrecognised_input_normalises_to_none(self, raw):
        assert normalize_polarity(raw) is None


class TestNormalizeIonisationSource:
    @pytest.mark.parametrize(
        ("raw", "label", "accession"),
        [
            ("MALDI", "MALDI", "MS:1000075"),
            ("matrix-assisted laser desorption ionization", "MALDI", "MS:1000075"),
            ("DESI", "DESI", "MS:1002011"),
            ("esi", "ESI", "MS:1000073"),
            ("TOF-SIMS", "SIMS", "MS:1000402"),
            ("secondary ion mass spectrometry", "SIMS", "MS:1000402"),
        ],
    )
    def test_known_spellings(self, raw, label, accession):
        result = normalize_ionisation_source(raw)
        assert result is not None
        assert result[0] == label
        assert result[1].accession == accession

    def test_unrecognised_input_normalises_to_none(self):
        assert normalize_ionisation_source("laser magic") is None


class TestNormalizeAnalyzer:
    @pytest.mark.parametrize(
        ("raw", "label", "accession"),
        [
            ("TOF", "TOF", "MS:1000084"),
            ("time-of-flight", "TOF", "MS:1000084"),
            ("qTOF", "TOF", "MS:1000084"),
            ("Orbitrap", "Orbitrap", "MS:1000484"),
            ("FTICR", "FTICR", "MS:1000079"),
            ("ion trap", "Ion trap", "MS:1000264"),
        ],
    )
    def test_known_spellings(self, raw, label, accession):
        result = normalize_analyzer(raw)
        assert result is not None
        assert result[0] == label
        assert result[1].accession == accession

    def test_unrecognised_input_normalises_to_none(self):
        assert normalize_analyzer("quadrupole array of mystery") is None
