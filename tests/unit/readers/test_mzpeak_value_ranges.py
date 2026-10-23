"""mzPeak positions and spectrum numbers past what a cast keeps (issue #527).

A float position just below a whole number was cut down onto its neighbour
(F9); a ``uint64`` position or spectrum number of 2**63 or more wrapped
negative (F17, F126); a region number of 2**31 or more overflowed the int32
``region_number`` column (F143).
"""

from __future__ import annotations

import zipfile

import anndata as ad
import numpy as np
import pyarrow as pa
import pytest

from tests.fixtures.mzpeak_builder import (
    POSITION_X_COLUMN,
    POSITION_Y_COLUMN,
    SCANS_MEMBER,
    Spectrum,
    build_mzpeak,
    grid_spectra,
    rewrite_table,
)
from thyra.convert import convert_msi
from thyra.errors import ConversionRefused
from thyra.readers.mzpeak import MzPeakReader


def _coordinates(archive):
    with MzPeakReader(archive) as reader:
        return sorted(coords for coords, _, _ in reader.iter_spectra())


def _refusal(archive) -> str:
    with pytest.raises(ConversionRefused) as refused:
        with MzPeakReader(archive) as reader:
            reader.get_essential_metadata()
            list(reader.iter_spectra())
    return str(refused.value)


def _floats(edit_x):
    """Store both positions as float64, x through ``edit_x``."""

    def edit(table):
        xs = [edit_x(float(x)) for x in table.column(POSITION_X_COLUMN).to_pylist()]
        ys = [float(y) for y in table.column(POSITION_Y_COLUMN).to_pylist()]
        table = table.set_column(
            table.column_names.index(POSITION_X_COLUMN),
            POSITION_X_COLUMN,
            pa.array(xs, pa.float64()),
        )
        return table.set_column(
            table.column_names.index(POSITION_Y_COLUMN),
            POSITION_Y_COLUMN,
            pa.array(ys, pa.float64()),
        )

    return edit


class TestFloatPositions:
    def test_a_float_a_hair_below_a_whole_number_keeps_its_pixel(self, tmp_path):
        """F9: x = 3 written as 0.3 / 0.1 is 2.9999999999999996."""
        plain = build_mzpeak(tmp_path / "plain.mzpeak", grid_spectra(3, 2))
        noisy = build_mzpeak(tmp_path / "noisy.mzpeak", grid_spectra(3, 2))
        rewrite_table(
            noisy, SCANS_MEMBER, _floats(lambda x: 0.3 / 0.1 if x == 3 else x)
        )

        assert _coordinates(noisy) == _coordinates(plain)
        assert len(set(_coordinates(noisy))) == 6

    def test_a_float_between_pixels_is_refused(self, tmp_path):
        archive = build_mzpeak(tmp_path / "a.mzpeak", grid_spectra(3, 2))
        rewrite_table(archive, SCANS_MEMBER, _floats(lambda x: x + 0.5))

        message = _refusal(archive)
        assert "x positions" in message
        assert "not whole numbers within int64" in message


class TestPositionsPastInt64:
    @pytest.mark.parametrize("shift", [2**63 + 4, 2**63 - 2])
    def test_a_uint64_position_of_2_63_or_more_is_refused(self, tmp_path, shift):
        """F17: all at or above 2**63, and straddling it."""
        archive = build_mzpeak(tmp_path / "a.mzpeak", grid_spectra(3, 3))

        def edit(table):
            xs = [shift + x for x in table.column(POSITION_X_COLUMN).to_pylist()]
            return table.set_column(
                table.column_names.index(POSITION_X_COLUMN),
                pa.field(POSITION_X_COLUMN, pa.uint64()),
                pa.array(xs, pa.uint64()),
            )

        rewrite_table(archive, SCANS_MEMBER, edit)

        assert "x positions" in _refusal(archive)


def _shift_spectrum_numbers(k: int):
    """Shift every spectrum number of every member by ``k``, as uint64."""

    def shift(column):
        values = column.to_numpy(zero_copy_only=False).astype(np.uint64)
        return pa.array(values + np.uint64(k), pa.uint64())

    def edit(table):
        columns = []
        for name in table.column_names:
            column = table.column(name).combine_chunks()
            if name in ("index", "source_index"):
                column = shift(column)
            elif pa.types.is_struct(column.type):
                column = pa.StructArray.from_arrays(
                    [
                        (
                            shift(column.field(f.name))
                            if f.name == "spectrum_index"
                            else column.field(f.name)
                        )
                        for f in column.type
                    ],
                    fields=list(column.type),
                )
            columns.append(column)
        return pa.Table.from_arrays(columns, names=table.column_names)

    return edit


class TestSpectrumNumbersPastInt64:
    @pytest.mark.parametrize("layout", ["point", "chunk"])
    def test_a_spectrum_number_of_2_63_or_more_is_refused(self, tmp_path, layout):
        """F126: they were stored as negative acquisition_order values."""
        archive = build_mzpeak(tmp_path / "a.mzpeak", grid_spectra(3, 2), layout=layout)

        with zipfile.ZipFile(archive) as z:
            members = [m for m in z.namelist() if m.endswith(".parquet")]
        for member in members:
            rewrite_table(archive, member, _shift_spectrum_numbers(2**63 - 3))

        message = _refusal(archive)
        assert "not whole numbers within int64" in message


class TestRegionNumberPastInt32:
    @pytest.mark.parametrize("number", [2**31, 2**40])
    def test_the_region_number_is_stored_as_int64(self, tmp_path, number):
        """F143: the int32 column ended in an OverflowError."""

        def region(n, frames, xs):
            return {
                "region_number": n,
                "name": None,
                "frames": frames,
                "raster_step_um": None,
                "x_index": xs,
                "y_index": [11, 12],
            }

        spectra = [
            Spectrum(x, y, [100.0, 101.0, 102.0], [1.0 + x, 2.0 + y, 3.0])
            for y in (1, 2)
            for x in range(1, 7)
        ]
        archive = build_mzpeak(
            tmp_path / "r.mzpeak",
            spectra,
            file_metadata={
                "imaging": {
                    "coordinate_base": 1,
                    "position_offset": {"x": 100, "y": 10},
                },
                "bruker_maldi": {
                    "regions": [region(number, 8, [101, 104]), region(0, 4, [105, 106])]
                },
            },
        )
        out = tmp_path / "r.zarr"
        assert convert_msi(
            archive,
            out,
            resampling_config={"method": "nearest_neighbor", "target_bins": 50},
        )

        obs = ad.read_zarr(out / "tables" / "msi_dataset_z0").obs
        assert str(obs["region_number"].dtype) == "int64"
        assert set(obs["region_number"]) == {number, 0}
