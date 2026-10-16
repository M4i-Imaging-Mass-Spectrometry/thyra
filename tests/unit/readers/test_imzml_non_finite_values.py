"""Non-finite values in an imzML, end to end (D55).

A NaN, infinite or negative intensity is not a measurement and the
converter drops it before resampling; an infinite m/z is no position on
any axis. The intensity threshold and the source mass range used to take
them in anyway: the threshold summed them into a peak the store never
holds, and one infinite m/z made the default axis an infinite bin count.
"""

import logging
from pathlib import Path

import anndata
import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.convert import convert_msi
from thyra.readers.imzml.imzml_reader import ImzMLReader

#: What the command line passes when no resampling flag is given; the API's
#: ``resampling_config=None`` means ``--no-resample``.
CLI_DEFAULT = {"method": "auto", "axis_type": "auto"}


def _write(path: Path, first_mzs, first_intensities) -> Path:
    """Two processed float64 spectra: the given one, and a plain control."""
    with ImzMLWriter(
        str(path), mode="processed", mz_dtype=np.float64, intensity_dtype=np.float64
    ) as writer:
        writer.addSpectrum(np.asarray(first_mzs), np.asarray(first_intensities), (1, 1))
        writer.addSpectrum(np.array([400.0, 500.0]), np.array([100.0, 100.0]), (2, 1))
    return path


def _stored_at(out: Path, mz: float) -> float:
    """Pixel 1's value in the column at ``mz`` of a raw-axis store."""
    table = anndata.read_zarr(str(out / "tables" / "ds_z0"))
    columns = np.flatnonzero(table.var["mz"].to_numpy() == mz)
    if columns.size == 0:
        return 0.0
    return float(table.X.toarray()[0, columns[0]])


class TestTheThresholdSumsMeasurementsOnly:
    @pytest.mark.parametrize(
        "copies, expected",
        [
            # A NaN or a negative copy no longer removes the finite share.
            ([np.nan, 100.0], 100.0),
            ([-30.0, 60.0], 60.0),
            # An infinite copy no longer keeps a share below the threshold.
            ([np.inf, 10.0], 0.0),
        ],
    )
    def test_the_store_holds_what_it_holds_without_a_threshold(
        self, tmp_path, copies, expected
    ):
        src = _write(tmp_path / "copies.imzML", [400.0, 500.0, 500.0], [100.0, *copies])
        out = tmp_path / "copies.zarr"

        assert convert_msi(
            str(src),
            str(out),
            dataset_id="ds",
            pixel_size_um=10,
            resampling_config=None,
            reader_options={"intensity_threshold": 50.0},
        )

        assert _stored_at(out, 500.0) == expected

    def test_an_infinite_mz_does_not_lift_the_largest_finite_peak(self, tmp_path):
        src = _write(
            tmp_path / "inf_mz.imzML", [400.0, 500.0, np.inf], [100.0, 10.0, 100.0]
        )
        out = tmp_path / "inf_mz.zarr"

        assert convert_msi(
            str(src),
            str(out),
            dataset_id="ds",
            pixel_size_um=10,
            resampling_config={"min_mz": 300.0, "max_mz": 600.0},
            reader_options={"intensity_threshold": 50.0},
        )

        table = anndata.read_zarr(str(out / "tables" / "ds_z0"))
        near = np.abs(table.var["mz"].to_numpy() - 500.0) < 0.5
        assert table.X.toarray()[0, near].sum() == 0.0


class TestAnInfiniteMzLeavesTheRangeFinite:
    @pytest.mark.parametrize("bad_mz", [np.inf, -np.inf])
    def test_the_source_range_is_the_finite_mz(self, tmp_path, bad_mz):
        src = _write(
            tmp_path / "range.imzML", [400.0, 500.0, bad_mz], [100.0, 10.0, 100.0]
        )

        with ImzMLReader(src) as reader:
            assert reader.get_essential_metadata().mass_range == (400.0, 500.0)

    def test_a_nan_beside_the_largest_mz_keeps_it(self, tmp_path):
        """The NaN used to drop its whole spectrum out of the running extrema."""
        src = _write(tmp_path / "nan.imzML", [450.0, 900.0, np.nan], [1.0, 1.0, 1.0])

        with ImzMLReader(src) as reader:
            assert reader.get_essential_metadata().mass_range == (400.0, 900.0)

    def test_the_default_conversion_succeeds(self, tmp_path):
        """It ended in ``OverflowError: cannot convert float infinity``."""
        src = _write(
            tmp_path / "default.imzML", [400.0, 500.0, np.inf], [100.0, 10.0, 100.0]
        )

        assert convert_msi(
            str(src),
            str(tmp_path / "default.zarr"),
            dataset_id="ds",
            pixel_size_um=10,
            resampling_config=CLI_DEFAULT,
        )

    def test_the_raw_axis_still_refuses_it(self, tmp_path, thyra_logs):
        """``--no-resample``: the axis is the source's m/z, inf included."""
        src = _write(
            tmp_path / "raw.imzML", [400.0, 500.0, np.inf], [100.0, 10.0, 100.0]
        )

        with thyra_logs("thyra", logging.ERROR) as records:
            assert not convert_msi(
                str(src),
                str(tmp_path / "raw.zarr"),
                dataset_id="ds",
                pixel_size_um=10,
                resampling_config=None,
            )
        assert "non-finite" in records.text
