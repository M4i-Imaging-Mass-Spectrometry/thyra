"""Tests for the synthetic example data generator."""

import numpy as np
import pytest

from tests.fixtures.imzml_parser import production_parser
from thyra.tools.make_example_data import generate_example_imzml


@pytest.fixture
def example_dataset(tmp_path):
    """Generate a small example dataset once for the whole module."""
    path = generate_example_imzml(
        tmp_path / "example.imzML", n_x=12, n_y=9, n_mz_bins=400, pixel_size_um=30.0
    )
    return path


def test_writes_imzml_and_ibd(example_dataset):
    """Both halves of the imzML pair are written."""
    assert example_dataset.exists()
    assert example_dataset.with_suffix(".ibd").exists()


def test_grid_and_axis_dimensions(example_dataset):
    """Spectra cover the full grid and share one m/z axis."""
    with production_parser(example_dataset) as parser:
        assert len(parser.coordinates) == 12 * 9

        mzs, intensities = parser.getspectrum(0)
        assert len(mzs) == 400
        assert len(intensities) == 400
        assert mzs[0] == pytest.approx(250.0)
        assert mzs[-1] == pytest.approx(1200.0)


def test_pixel_size_is_written_to_metadata(example_dataset):
    """Pixel size cvParams are present so Thyra can auto-detect them."""
    text = example_dataset.read_text(encoding="ISO-8859-1")
    assert 'accession="IMS:1000046"' in text
    assert 'accession="IMS:1000047"' in text

    with production_parser(example_dataset) as parser:
        assert float(parser.imzmldict["pixel size x"]) == pytest.approx(30.0)
        assert float(parser.imzmldict["pixel size y"]) == pytest.approx(30.0)


def test_coordinates_are_one_based(example_dataset):
    """imzML coordinates start at 1, not 0."""
    with production_parser(example_dataset) as parser:
        xs = [c[0] for c in parser.coordinates]
        ys = [c[1] for c in parser.coordinates]
        assert min(xs) == 1 and max(xs) == 12
        assert min(ys) == 1 and max(ys) == 9


def test_regions_are_spatially_distinct(example_dataset):
    """The inner structure carries peaks the outer region does not.

    m/z 888.6 belongs to the inner region only, so its intensity must vary
    across pixels rather than being uniform -- this is what makes the ion
    images in the tutorial show contrast.
    """
    with production_parser(example_dataset) as parser:
        mzs, _ = parser.getspectrum(0)
        idx = int(np.abs(mzs - 888.6).argmin())

        values = np.array(
            [parser.getspectrum(i)[1][idx] for i in range(len(parser.coordinates))]
        )
    assert values.max() > 10 * max(values.min(), 1.0)


def test_same_seed_reproduces_identical_intensities(tmp_path):
    """The same seed reproduces identical intensities.

    Renamed from ``test_output_is_deterministic``: the docstring was
    already the true statement while the name promised more than the tool
    delivers. Three places claimed the *files* were byte-identical for a
    fixed seed and none of them was right (issue #314).
    """
    kwargs = dict(n_x=6, n_y=5, n_mz_bins=200)
    a = generate_example_imzml(tmp_path / "a" / "d.imzML", seed=7, **kwargs)
    b = generate_example_imzml(tmp_path / "b" / "d.imzML", seed=7, **kwargs)

    with production_parser(a) as parser_a, production_parser(b) as parser_b:
        spec_a = parser_a.getspectrum(3)[1]
        spec_b = parser_b.getspectrum(3)[1]
    np.testing.assert_array_equal(spec_a, spec_b)


def test_the_ibd_past_its_uuid_header_is_byte_identical(tmp_path):
    """What determinism here actually means, stated positively.

    The 16-byte UUID header is the only non-reproducible range in the
    binary; everything after it is the same bytes for the same seed.
    """
    kwargs = dict(n_x=6, n_y=5, n_mz_bins=200)
    a = generate_example_imzml(tmp_path / "a" / "d.imzML", seed=7, **kwargs)
    b = generate_example_imzml(tmp_path / "b" / "d.imzML", seed=7, **kwargs)

    ibd_a = a.with_suffix(".ibd").read_bytes()
    ibd_b = b.with_suffix(".ibd").read_bytes()
    assert ibd_a[16:] == ibd_b[16:]
    assert ibd_a[:16] != ibd_b[:16]  # and the UUID is genuinely fresh


def test_the_output_path_does_not_reach_the_file(tmp_path):
    """pyimzML derives ``<run id=...>`` from the path it was handed.

    So the generated imzML carried the caller's absolute output directory
    -- into a file the tutorial tells people to share, and one that
    therefore differed between two machines running the same command.
    """
    path = generate_example_imzml(
        tmp_path / "nested" / "sample.imzML", n_x=4, n_y=4, n_mz_bins=100
    )
    text = path.read_text(encoding="ISO-8859-1")

    assert str(tmp_path) not in text
    assert "nested" not in text


def test_the_run_id_is_the_stem_alone(tmp_path):
    """No separator of either flavour, so two platforms agree."""
    import xml.etree.ElementTree as ET  # nosec B405 - our own output

    path = generate_example_imzml(
        tmp_path / "nested" / "sample.imzML", n_x=4, n_y=4, n_mz_bins=100
    )
    root = ET.parse(path).getroot()  # nosec B314 - our own output
    run = root.find(".//{http://psi.hupo.org/ms/mzml}run")
    assert run is not None
    run_id = run.get("id")

    assert run_id == "sample"
    assert "/" not in run_id and "\\" not in run_id


def test_extension_is_normalised(tmp_path):
    """A path without the .imzML suffix still produces a valid pair."""
    path = generate_example_imzml(tmp_path / "noext.imzML", n_x=4, n_y=4, n_mz_bins=100)
    assert path.suffix == ".imzML"
    assert path.with_suffix(".ibd").exists()
