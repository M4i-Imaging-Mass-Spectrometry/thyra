"""End-to-end conversion of mzPeak archives, and parity against imzML.

The parity test is the acceptance gate for this reader: the same spectra,
written once as imzML and once as mzPeak, must convert to the same store.
Anything the reader gets wrong about coordinates, ordering, the mass axis or
the payload shows up as a difference here rather than as a plausible-looking
store nobody checks.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from thyra.convert import convert_msi

spatialdata = pytest.importorskip("spatialdata", reason="SpatialData not installed")

DATASET_ID = "parity"
TABLE = f"{DATASET_ID}_z0"
PIXEL_SIZE_UM = 25.0


def _write_imzml(path: Path, spectra) -> Path:
    """Write the same spectra as a processed-mode imzML."""
    from pyimzml.ImzMLWriter import ImzMLWriter

    with ImzMLWriter(str(path), mode="processed") as writer:
        for spectrum in spectra:
            writer.addSpectrum(
                spectrum.mzs,
                spectrum.intensities,
                (spectrum.x, spectrum.y, 1),
            )
    return path


def _matching_mzpeak(path: Path, spectra, **kwargs) -> Path:
    """Build the mzPeak counterpart of an imzML written from ``spectra``.

    The representation is declared centroid because that is what the imzML
    side reports: pyimzml's processed mode is detected as ``centroid
    spectrum``. Representation feeds the resampling decision tree, so leaving
    the two files disagreeing about it would compare a linear axis against a
    nonlinear one and call the reader wrong for it.
    """
    kwargs.setdefault("spectrum_representation", "centroid spectrum")
    return build_mzpeak(path, spectra, **kwargs)


def _assert_same_representation(imzml: Path, mzpeak: Path) -> None:
    """Fail early, and legibly, if the two inputs disagree on representation.

    Spectrum representation selects the mass-axis law, so a mismatch shows up
    downstream as two axes of different shape -- a wall of unequal floats that
    says nothing about the cause. Checking it here names the problem.
    """
    from thyra.readers.imzml import ImzMLReader
    from thyra.readers.mzpeak import MzPeakReader

    with ImzMLReader(imzml) as reader:
        expected = reader.get_essential_metadata().spectrum_type
    with MzPeakReader(mzpeak) as reader:
        actual = reader.get_essential_metadata().spectrum_type

    assert actual == expected, (
        f"fixtures disagree on spectrum representation: imzML says "
        f"{expected!r}, mzPeak says {actual!r}; the resampled axis law is "
        f"chosen from this, so the stores cannot match"
    )


def _convert(source: Path, output: Path, resampling_config=None) -> None:
    """Convert one input, asserting the conversion reported success."""
    assert (
        convert_msi(
            str(source),
            str(output),
            format_type="spatialdata",
            dataset_id=DATASET_ID,
            pixel_size_um=PIXEL_SIZE_UM,
            resampling_config=resampling_config,
        )
        is True
    )


def _table(store: Path):
    """Load the single table out of a converted store."""
    return spatialdata.SpatialData.read(str(store)).tables[TABLE]


def _dense(table):
    """Densify X for comparison; the fixtures are deliberately tiny."""
    matrix = table.X
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


@pytest.mark.integration
class TestImzMLParity:
    """The same data through two readers must land in the same store."""

    @pytest.mark.parametrize(
        "resampling_config",
        [
            pytest.param(None, id="native-axis"),
            pytest.param(
                {"method": "nearest_neighbor", "target_bins": 64},
                id="resampled",
            ),
        ],
    )
    def test_stores_match(self, tmp_path, resampling_config):
        """X, obs coordinates and the m/z axis all agree.

        Run twice: once on the native union axis, where any disagreement
        about which m/z values exist is visible directly, and once through
        resampling, where the axis is rebuilt and the test instead pins the
        binned intensities.
        """
        spectra = grid_spectra(3, 2, n_points=6)

        imzml = _write_imzml(tmp_path / "source.imzML", spectra)
        mzpeak = _matching_mzpeak(tmp_path / "source.mzpeak", spectra)
        _assert_same_representation(imzml, mzpeak)

        imzml_store = tmp_path / "from_imzml.zarr"
        mzpeak_store = tmp_path / "from_mzpeak.zarr"
        _convert(imzml, imzml_store, resampling_config)
        _convert(mzpeak, mzpeak_store, resampling_config)

        reference = _table(imzml_store)
        candidate = _table(mzpeak_store)

        assert candidate.n_obs == reference.n_obs
        assert candidate.n_vars == reference.n_vars
        np.testing.assert_allclose(
            candidate.var["mz"].to_numpy(),
            reference.var["mz"].to_numpy(),
            rtol=0,
            atol=0,
        )
        for column in ("spatial_x", "spatial_y"):
            np.testing.assert_array_equal(
                candidate.obs[column].to_numpy(),
                reference.obs[column].to_numpy(),
            )
        np.testing.assert_allclose(
            _dense(candidate), _dense(reference), rtol=1e-12, atol=0
        )

    def test_sparse_acquisition_matches(self, tmp_path):
        """Missing pixels land in the same places from both formats.

        Neither reader may densify: an unacquired position must stay
        unacquired rather than becoming a row of zeros.
        """
        spectra = grid_spectra(3, 3, n_points=5, skip=[(2, 2), (3, 1)])

        imzml = _write_imzml(tmp_path / "sparse.imzML", spectra)
        mzpeak = _matching_mzpeak(tmp_path / "sparse.mzpeak", spectra)

        imzml_store = tmp_path / "sparse_imzml.zarr"
        mzpeak_store = tmp_path / "sparse_mzpeak.zarr"
        _convert(imzml, imzml_store)
        _convert(mzpeak, mzpeak_store)

        reference = _table(imzml_store)
        candidate = _table(mzpeak_store)

        assert candidate.n_obs == reference.n_obs == 7
        np.testing.assert_array_equal(
            candidate.obs[["spatial_x", "spatial_y"]].to_numpy(),
            reference.obs[["spatial_x", "spatial_y"]].to_numpy(),
        )
        np.testing.assert_allclose(_dense(candidate), _dense(reference))

    def test_row_group_split_does_not_change_the_store(self, tmp_path):
        """Writer row-group choices are invisible in the output.

        Spectra straddling row-group boundaries are the one place the reader
        does non-trivial bookkeeping, so this pins it end to end rather than
        only at the iteration layer.
        """
        spectra = grid_spectra(3, 2, n_points=7)

        whole = _matching_mzpeak(tmp_path / "whole.mzpeak", spectra)
        split = _matching_mzpeak(tmp_path / "split.mzpeak", spectra, row_group_size=2)

        whole_store = tmp_path / "whole.zarr"
        split_store = tmp_path / "split.zarr"
        _convert(whole, whole_store)
        _convert(split, split_store)

        np.testing.assert_allclose(
            _dense(_table(split_store)), _dense(_table(whole_store))
        )


@pytest.mark.integration
class TestConversion:
    """Conversion behaviour specific to mzPeak inputs."""

    def test_pixel_size_is_detected_from_the_archive(self, tmp_path):
        """An archive declaring IMS:1000046/47 needs no --pixel-size."""
        archive = build_mzpeak(
            tmp_path / "sized.mzpeak",
            grid_spectra(2, 2, n_points=4),
            pixel_size=(12.0, 12.0),
        )
        output = tmp_path / "sized.zarr"

        assert (
            convert_msi(
                str(archive),
                str(output),
                format_type="spatialdata",
                dataset_id=DATASET_ID,
            )
            is True
        )
        assert output.exists()

    def test_declared_grid_extent_is_writable_provenance(self, tmp_path):
        """An archive declaring IMS:1000042/43 still saves.

        The extent is recorded in ``uns``, and zarr writes a dict key as a
        directory name. Keying it by accession made the store unwritable on
        Windows, where a colon is not legal in a path -- and no fixture that
        stopped at the metadata object could see it, because the failure is
        in the save.
        """
        archive = build_mzpeak(
            tmp_path / "extent.mzpeak",
            grid_spectra(2, 2, n_points=4),
            grid=(32, 32),
            spectrum_representation="centroid spectrum",
        )
        output = tmp_path / "extent.zarr"
        _convert(archive, output)

        table = _table(output)
        extent = table.uns["acquisition_params"]["declared_grid_extent"]
        assert extent["max_count_of_pixels_x"] == 32
        assert extent["max_count_of_pixels_y"] == 32

    def test_null_padded_archive_converts(self, tmp_path):
        """Null-pair padding does not reach the stored matrix.

        A NaN surviving into X would poison every downstream statistic, so
        this asserts the stored values are finite as well as correct.
        """
        spectra = grid_spectra(2, 2, n_points=6)
        archive = build_mzpeak(tmp_path / "padded.mzpeak", spectra, null_pair_after=3)
        output = tmp_path / "padded.zarr"
        _convert(archive, output)

        table = _table(output)
        values = _dense(table)
        assert np.isfinite(values).all()
        assert np.isfinite(table.var["mz"].to_numpy()).all()
        assert table.n_vars == 6


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("THYRA_MZPEAK_REFERENCE_ARCHIVE"),
    reason=(
        "Set THYRA_MZPEAK_REFERENCE_ARCHIVE to an archive written by the "
        "HUPO-PSI reference converter to run the interoperability check"
    ),
)
def test_reference_archive_converts(tmp_path):
    """Convert an archive written by the reference implementation.

    The constructed fixtures pin the schema this reader was written against,
    but only a file produced by someone else's writer can show the two agree.
    The reference sample archives are not vendored: that repository publishes
    no licence, so redistributing its data in a package that ships to PyPI
    would not be ours to do. The check therefore runs opt-in, against a path
    the operator supplies.
    """
    archive = Path(os.environ["THYRA_MZPEAK_REFERENCE_ARCHIVE"])
    if not archive.exists():
        pytest.skip(f"reference archive not found: {archive}")

    output = tmp_path / "reference.zarr"
    _convert(archive, output)

    table = _table(output)
    assert table.n_obs > 0
    assert table.n_vars > 0
    values = _dense(table)
    assert np.isfinite(values).all()
    assert np.isfinite(table.var["mz"].to_numpy()).all()
