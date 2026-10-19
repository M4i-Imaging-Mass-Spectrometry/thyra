"""The imzML extractor describes only the spectra that are converted (D70).

An MSn spectrum left out beside MS1 on its pixel adds nothing to the
spectrum count, the peak count or the mass range, in either file mode.
Every pixel keeps a spectrum, so the dimensions do not change.
"""

import pytest

from tests.fixtures.imzml_levels import write_imzml_with_levels
from thyra.metadata.extractors.imzml_extractor import ImzMLMetadataExtractor
from thyra.readers.imzml.imzml_reader import ImzMLReader


def _essential(path):
    """The extractor run on the parser the reader built (with its hook)."""
    reader = ImzMLReader(path)
    try:
        reader._ensure_parser_initialized()
        assert reader.parser is not None
        return ImzMLMetadataExtractor(reader.parser, path).get_essential()
    finally:
        reader.close()


@pytest.mark.parametrize("how", ["own", "group"])
def test_processed_range_and_counts_leave_the_msn_spectrum_out(tmp_path, how):
    spectra = [
        (1, 1, 1, [100.0, 110.0], [1.0, 2.0]),
        (1, 1, 1, [50.0, 900.0, 901.0], [5.0, 5.0, 5.0]),
        (2, 1, 1, [120.0], [3.0]),
    ]
    path = write_imzml_with_levels(tmp_path / "p.imzML", spectra, [1, 2, 1], how=how)

    essential = _essential(path)

    assert essential.n_spectra == 2
    assert essential.total_peaks == 3
    assert essential.mass_range == (100.0, 120.0)
    assert essential.dimensions == (2, 1, 1)


def test_continuous_counts_leave_the_msn_spectrum_out(tmp_path):
    mzs = [150.0, 250.0, 350.0]
    spectra = [
        (1, 1, 1, mzs, [1.0] * 3),
        (2, 1, 1, mzs, [2.0] * 3),
        (2, 1, 1, mzs, [9.0] * 3),
    ]
    path = write_imzml_with_levels(
        tmp_path / "c.imzML", spectra, [1, 1, 2], mode="continuous"
    )

    essential = _essential(path)

    assert essential.n_spectra == 2
    assert essential.total_peaks == 6
    assert essential.mass_range == (150.0, 350.0)


def test_a_file_stating_no_levels_counts_every_spectrum(tmp_path):
    spectra = [
        (1, 1, 1, [100.0], [1.0]),
        (1, 1, 1, [900.0], [1.0]),
    ]
    path = write_imzml_with_levels(tmp_path / "n.imzML", spectra, [None, None])

    essential = _essential(path)

    assert essential.n_spectra == 2
    assert essential.mass_range == (100.0, 900.0)
