"""Tests for the synthetic example data generator."""

import numpy as np
import pytest
from pyimzml.ImzMLParser import ImzMLParser

from thyra.tools.make_example_data import generate_example_imzml


@pytest.fixture
def example_dataset(tmp_path):
    """Generate a small example dataset once for the whole module."""
    path = generate_example_imzml(
        tmp_path / "example.imzML", n_x=12, n_y=9, n_mz_bins=400, pixel_size_um=30.0
    )
    return path


@pytest.mark.unit
def test_writes_imzml_and_ibd(example_dataset):
    """Both halves of the imzML pair are written."""
    assert example_dataset.exists()
    assert example_dataset.with_suffix(".ibd").exists()


@pytest.mark.unit
def test_grid_and_axis_dimensions(example_dataset):
    """Spectra cover the full grid and share one m/z axis."""
    parser = ImzMLParser(str(example_dataset))
    assert len(parser.coordinates) == 12 * 9

    mzs, intensities = parser.getspectrum(0)
    assert len(mzs) == 400
    assert len(intensities) == 400
    assert mzs[0] == pytest.approx(250.0)
    assert mzs[-1] == pytest.approx(1200.0)


@pytest.mark.unit
def test_pixel_size_is_written_to_metadata(example_dataset):
    """Pixel size cvParams are present so Thyra can auto-detect them."""
    text = example_dataset.read_text(encoding="ISO-8859-1")
    assert 'accession="IMS:1000046"' in text
    assert 'accession="IMS:1000047"' in text

    parser = ImzMLParser(str(example_dataset))
    assert float(parser.imzmldict["pixel size x"]) == pytest.approx(30.0)
    assert float(parser.imzmldict["pixel size y"]) == pytest.approx(30.0)


@pytest.mark.unit
def test_coordinates_are_one_based(example_dataset):
    """imzML coordinates start at 1, not 0."""
    parser = ImzMLParser(str(example_dataset))
    xs = [c[0] for c in parser.coordinates]
    ys = [c[1] for c in parser.coordinates]
    assert min(xs) == 1 and max(xs) == 12
    assert min(ys) == 1 and max(ys) == 9


@pytest.mark.unit
def test_regions_are_spatially_distinct(example_dataset):
    """The inner structure carries peaks the outer region does not.

    m/z 888.6 belongs to the inner region only, so its intensity must vary
    across pixels rather than being uniform -- this is what makes the ion
    images in the tutorial show contrast.
    """
    parser = ImzMLParser(str(example_dataset))
    mzs, _ = parser.getspectrum(0)
    idx = int(np.abs(mzs - 888.6).argmin())

    values = np.array(
        [parser.getspectrum(i)[1][idx] for i in range(len(parser.coordinates))]
    )
    assert values.max() > 10 * max(values.min(), 1.0)


@pytest.mark.unit
def test_output_is_deterministic(tmp_path):
    """The same seed reproduces identical intensities."""
    kwargs = dict(n_x=6, n_y=5, n_mz_bins=200)
    a = generate_example_imzml(tmp_path / "a" / "d.imzML", seed=7, **kwargs)
    b = generate_example_imzml(tmp_path / "b" / "d.imzML", seed=7, **kwargs)

    spec_a = ImzMLParser(str(a)).getspectrum(3)[1]
    spec_b = ImzMLParser(str(b)).getspectrum(3)[1]
    np.testing.assert_array_equal(spec_a, spec_b)


@pytest.mark.unit
def test_extension_is_normalised(tmp_path):
    """A path without the .imzML suffix still produces a valid pair."""
    path = generate_example_imzml(tmp_path / "noext.imzML", n_x=4, n_y=4, n_mz_bins=100)
    assert path.suffix == ".imzML"
    assert path.with_suffix(".ibd").exists()
