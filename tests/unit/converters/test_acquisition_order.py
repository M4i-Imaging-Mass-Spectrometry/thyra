"""``obs["acquisition_order"]``: when each row was acquired (design decision D26).

A table's rows follow the raster grid, so a drift plotted against them is
plotted against raster order, which a serpentine scan or a slide of several
regions does not follow. A reader that knows the order hands it over with
every spectrum; the converter keeps, per row, the smallest order among the
spectra the row sums, and writes nothing when the reader does not know.
"""

from pathlib import Path
from typing import List, Tuple

import numpy as np
import pytest

from thyra.converters.spatialdata import SpatialDataConverter
from thyra.converters.spatialdata.streaming_converter import _TableUnit
from thyra.core.base_extractor import MetadataExtractor
from thyra.core.base_reader import BaseMSIReader
from thyra.metadata.schema import MSI_OBS_ACQUISITION_ORDER_COLUMN
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata

COLUMN = MSI_OBS_ACQUISITION_ORDER_COLUMN
AXIS = np.linspace(100.0, 1000.0, 100)

Log = List[Tuple[Tuple[int, int, int], int]]

#: A 3x2 raster whose second row was acquired right to left, numbered the
#: way a Bruker file numbers its frames: from 11, with a gap at 14 where a
#: frame of another region would be.
SERPENTINE: Log = [
    ((0, 0, 0), 11),
    ((1, 0, 0), 12),
    ((2, 0, 0), 13),
    ((2, 1, 0), 15),
    ((1, 1, 0), 16),
    ((0, 1, 0), 17),
]


def _intensities(x: int, y: int) -> np.ndarray:
    """Two peaks on the axis, one moving with x and one with y."""
    intensities = np.zeros_like(AXIS)
    intensities[x * 10 + 20] = 100.0
    intensities[y * 10 + 50] = 200.0
    return intensities


class _Extractor(MetadataExtractor):
    def __init__(self, dims: Tuple[int, int, int]):
        super().__init__(None)
        self._dims = dims

    def _extract_essential_impl(self) -> EssentialMetadata:
        n = self._dims[0] * self._dims[1] * self._dims[2]
        return EssentialMetadata(
            dimensions=self._dims,
            coordinate_bounds=(
                0.0,
                float(self._dims[0] - 1),
                0.0,
                float(self._dims[1] - 1),
            ),
            mass_range=(100.0, 1000.0),
            pixel_size=None,
            n_spectra=n,
            total_peaks=2 * n,
            source_path="/mock/ordered",
        )

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        return ComprehensiveMetadata(
            essential=self._extract_essential_impl(),
            format_specific={},
            acquisition_params={},
            instrument_info={},
            raw_metadata={},
        )


class OrderedReader(BaseMSIReader):
    """Yields ``log`` as written: one spectrum per ``(coords, order)``."""

    def __init__(self, dims: Tuple[int, int, int], log: Log):
        super().__init__(Path("/mock/ordered"))
        self._dims = dims
        self._log = list(log)

    def _create_metadata_extractor(self) -> MetadataExtractor:
        return _Extractor(self._dims)

    def get_common_mass_axis(self) -> np.ndarray:
        return AXIS

    @property
    def has_acquisition_order(self) -> bool:
        return True

    def iter_spectra_with_acquisition_order(self):
        for coords, order in self._log:
            yield coords, order, AXIS, _intensities(coords[0], coords[1])

    def iter_spectra(self):
        for coords, _order in self._log:
            yield coords, AXIS, _intensities(coords[0], coords[1])

    def close(self) -> None:
        pass


class UnorderedReader(OrderedReader):
    """The same spectra from a reader that cannot say when they were acquired."""

    @property
    def has_acquisition_order(self) -> bool:
        return False

    def iter_spectra_with_acquisition_order(self):
        raise NotImplementedError("no order here")


def _tables_obs(reader: BaseMSIReader, tmp_path: Path, **kwargs):
    """Run both passes and return ``obs`` per table key, as finalize builds it."""
    converter = SpatialDataConverter(reader, tmp_path / "out.zarr", **kwargs)
    try:
        converter._initialize_conversion()
        state = converter._create_data_structures()
        converter._process_spectra(state)
        return {
            unit.key: converter._table_obs(unit) for unit in state.units if unit.n_rows
        }
    finally:
        converter.siblings.release_scratch()


class TestTheColumn:
    def test_each_row_carries_the_order_its_position_was_acquired_in(self, tmp_path):
        (obs,) = _tables_obs(OrderedReader((3, 2, 1), SERPENTINE), tmp_path).values()

        # The rows stay in grid order; only the column knows the scan went
        # back along the second row.
        assert obs.index.tolist() == ["0", "1", "2", "3", "4", "5"]
        assert obs[COLUMN].tolist() == [11, 12, 13, 17, 16, 15]
        assert obs[COLUMN].dtype == np.int64

    def test_it_is_the_last_column_before_the_instance_key(self, tmp_path):
        (obs,) = _tables_obs(OrderedReader((3, 2, 1), SERPENTINE), tmp_path).values()

        assert list(obs.columns) == [
            "y",
            "x",
            "region",
            "spatial_x",
            "spatial_y",
            "region_number",
            COLUMN,
        ]

    def test_sorting_by_it_replays_the_acquisition(self, tmp_path):
        (obs,) = _tables_obs(OrderedReader((3, 2, 1), SERPENTINE), tmp_path).values()

        replay = obs.sort_values(COLUMN)
        acquired = [(x, y) for (x, y, _z), _order in SERPENTINE]
        assert list(zip(replay["x"], replay["y"])) == acquired

    def test_a_reader_without_an_order_writes_no_column(self, tmp_path):
        reader = UnorderedReader((3, 2, 1), SERPENTINE)
        (obs,) = _tables_obs(reader, tmp_path).values()

        assert COLUMN not in obs.columns
        assert obs.index.tolist() == ["0", "1", "2", "3", "4", "5"]


class TestTheRowsItDescribes:
    def test_a_position_measured_twice_takes_the_earlier_order(self, tmp_path):
        # The repeat of (0, 0) is yielded last but was acquired first: the
        # row sums both spectra and was first measured at order 2.
        log: Log = [((0, 0, 0), 7), ((1, 0, 0), 8), ((0, 0, 0), 2)]
        (obs,) = _tables_obs(OrderedReader((2, 1, 1), log), tmp_path).values()

        assert obs[COLUMN].tolist() == [2, 8]

    def test_positions_without_a_row_do_not_shift_it(self, tmp_path):
        # (1, 0) and (0, 1) were never acquired and (5, 0) is off the grid,
        # so three positions are dropped and one spectrum is skipped. The
        # column is indexed by the same kept grid as x and y.
        log: Log = [((2, 1, 0), 1), ((0, 0, 0), 2), ((5, 0, 0), 3), ((1, 1, 0), 4)]
        (obs,) = _tables_obs(OrderedReader((3, 2, 1), log), tmp_path).values()

        assert obs.index.tolist() == ["0", "4", "5"]
        assert list(zip(obs["x"], obs["y"], obs[COLUMN])) == [
            (0, 0, 2),
            (1, 1, 4),
            (2, 1, 1),
        ]

    def test_plane_tables_and_a_volume_carry_the_same_orders(self, tmp_path):
        log: Log = [((0, 0, 1), 1), ((1, 0, 1), 2), ((1, 0, 0), 3), ((0, 0, 0), 4)]

        planes = _tables_obs(OrderedReader((2, 1, 2), log), tmp_path / "planes")
        assert planes["msi_dataset_z0"][COLUMN].tolist() == [4, 3]
        assert planes["msi_dataset_z1"][COLUMN].tolist() == [1, 2]

        volume = _tables_obs(
            OrderedReader((2, 1, 2), log), tmp_path / "volume", handle_3d=True
        )["msi_dataset"]
        assert volume["z"].tolist() == [0, 0, 1, 1]
        assert volume[COLUMN].tolist() == [4, 3, 1, 2]

    def test_a_spectrum_without_its_order_is_refused_not_left_blank(self):
        unit = _TableUnit("t", "t_pixels", 0, 4, 10, (2, 2), records_order=True)

        with pytest.raises(ValueError, match="without its acquisition order"):
            unit.count(0, np.array([1]), np.array([1.0]))


def test_the_column_reaches_the_store_as_int64(tmp_path):
    spatialdata = pytest.importorskip("spatialdata")
    from thyra.utils.windows_paths import (
        prepare_zarr_output_path,
        prepare_zarr_read_path,
    )

    out = prepare_zarr_output_path(tmp_path / "ordered.zarr", "ordered")
    converter = SpatialDataConverter(
        OrderedReader((3, 2, 1), SERPENTINE),
        out,
        dataset_id="ordered",
        pixel_size_um=10.0,
    )
    assert converter.convert(), "conversion reported failure"

    obs = spatialdata.read_zarr(prepare_zarr_read_path(out)).tables["ordered_z0"].obs
    assert obs[COLUMN].dtype == np.int64
    assert obs[COLUMN].tolist() == [11, 12, 13, 17, 16, 15]
