"""Opening a store and choosing its tables by their columns."""

from __future__ import annotations

import numpy as np
import pytest
import zarr

from tests.fixtures.hand_built_store import (
    L_SHAPE,
    MZ,
    dense_matrix,
    make_table,
    write_store,
)
from thyra.quicklook import store as store_module
from thyra.quicklook.store import (
    MOBILITY,
    MSI,
    MSMS,
    OTHER,
    QuicklookError,
    open_store,
    table_kind,
)


def _siblings():
    """Three tables whose names say the opposite of what they hold."""
    n = MZ.size
    return {
        "a_spectra": make_table(extra_var={"mobility": np.linspace(0.8, 1.2, n)}),
        "the_mobility_one": make_table(),
        "b_spectra": make_table(extra_var={"precursor_mz": np.full(n, 600.0)}),
    }


class TestTableKind:
    @pytest.mark.parametrize(
        "columns, kind",
        [
            (["mz"], MSI),
            (["mz", "tof"], MSI),
            (["mz", "mobility", "mz_index"], MOBILITY),
            (["mz", "precursor_mz", "precursor_index"], MSMS),
            (["_index"], OTHER),
        ],
    )
    def test_the_var_columns_decide(self, columns, kind):
        assert table_kind(columns) == kind


class TestTheMainTableIsChosenByItsColumns:
    def test_sibling_tables_are_not_tables_of_spectra(self, tmp_path):
        store = open_store(write_store(tmp_path / "s.zarr", _siblings()))

        assert [table.name for table in store.msi_tables()] == ["the_mobility_one"]

    def test_every_table_is_still_listed_with_its_kind_and_shape(self, tmp_path):
        store = open_store(write_store(tmp_path / "s.zarr", _siblings()))

        listed = {table.name: (table.kind, table.shape) for table in store.tables}
        shape = (len(L_SHAPE), MZ.size)
        assert listed == {
            "a_spectra": (MOBILITY, shape),
            "b_spectra": (MSMS, shape),
            "the_mobility_one": (MSI, shape),
        }

    def test_each_slice_table_is_a_table_of_spectra(self, tmp_path):
        tables = {"ds_z0": make_table(), "ds_z1": make_table()}
        store = open_store(write_store(tmp_path / "s.zarr", tables))

        assert [table.name for table in store.msi_tables()] == ["ds_z0", "ds_z1"]

    def test_slice_ten_comes_after_slice_two(self, tmp_path):
        tables = {f"ds_z{index}": make_table() for index in (10, 2, 0, 1)}
        store = open_store(write_store(tmp_path / "s.zarr", tables))

        assert [table.name for table in store.msi_tables()] == [
            "ds_z0",
            "ds_z1",
            "ds_z2",
            "ds_z10",
        ]

    def test_a_store_with_only_sibling_tables_is_refused(self, tmp_path):
        tables = _siblings()
        del tables["the_mobility_one"]
        store = open_store(write_store(tmp_path / "s.zarr", tables))

        with pytest.raises(QuicklookError, match="holds no table of spectra"):
            store.msi_tables()


class TestWhatIsNotAStore:
    def test_a_plain_folder_is_refused(self, tmp_path):
        with pytest.raises(QuicklookError, match="is not a Thyra store"):
            open_store(tmp_path)

    def test_a_file_is_refused(self, tmp_path):
        path = tmp_path / "notes.txt"
        path.write_text("not a store", encoding="utf-8")

        with pytest.raises(QuicklookError, match="is not a Thyra store"):
            open_store(path)

    def test_a_zarr_group_without_tables_is_refused(self, tmp_path):
        zarr.open_group(str(tmp_path / "other.zarr"), mode="w")

        with pytest.raises(QuicklookError, match="has no 'tables' group"):
            open_store(tmp_path / "other.zarr")


class TestMsiTable:
    def test_a_mass_axis_that_goes_backwards_is_refused(self, tmp_path):
        """A binary search over it would return the wrong columns silently."""
        backwards = MZ.copy()
        backwards[[2, 3]] = backwards[[3, 2]]
        path = write_store(tmp_path / "s.zarr", {"t": make_table(mz=backwards)})

        with pytest.raises(QuicklookError, match="not strictly increasing"):
            open_store(path).open_table("t")

    def test_positions_carry_z_only_for_a_volume(self, tmp_path):
        volume = [(x, y, z) for z in (0, 1) for y in (0, 1) for x in (0, 1)]
        tables = {"flat": make_table(), "volume": make_table(positions=volume)}
        store = open_store(write_store(tmp_path / "s.zarr", tables))

        assert store.open_table("flat").positions().z is None
        np.testing.assert_array_equal(
            store.open_table("volume").positions().z, [0, 0, 0, 0, 1, 1, 1, 1]
        )

    def test_a_window_reads_its_own_columns_and_no_others(self, tmp_path, monkeypatch):
        """The whole point on a large store: the read follows the window."""
        path = write_store(tmp_path / "s.zarr", {"t": make_table()})
        table = open_store(path).open_table("t")
        asked = []

        class Spy:
            def __init__(self, name, array):
                self.name, self.array = name, array

            def __getitem__(self, key):
                asked.append((self.name, key.start, key.stop))
                return self.array[key]

        real = store_module._array
        monkeypatch.setattr(
            store_module,
            "_array",
            lambda group, name, where: Spy(name, real(group, name, where)),
        )

        values = table.window_sum(4, 7)

        np.testing.assert_array_equal(values, dense_matrix(len(L_SHAPE))[:, 4:7].sum(1))
        # Columns 0, 1 and 2 hold 8 values each and column 3 none: the block
        # starts at stored value 24 and holds 3 x 8.
        assert sorted(asked) == [
            ("data", 24, 48),
            ("indices", 24, 48),
            ("indptr", 4, 8),
        ]

    def test_the_metadata_block_is_read_from_the_table(self, tmp_path):
        block = {"provenance": {"thyra_version": "4.7.0"}}
        tables = {
            "with": make_table(uns={"msi_metadata": block}),
            "without": make_table(),
        }
        store = open_store(write_store(tmp_path / "s.zarr", tables))

        assert store.open_table("with").metadata_block() == block
        assert store.open_table("without").metadata_block() == {}

    def test_a_mean_spectrum_of_the_wrong_length_is_not_used(self, tmp_path):
        table = make_table(uns={"average_spectrum": np.ones(3)})
        path = write_store(tmp_path / "s.zarr", {"t": table})

        assert open_store(path).open_table("t").average_spectrum() is None


class TestImages:
    def test_the_tic_image_is_found_by_its_table_and_read_from_its_metadata(
        self, tmp_path
    ):
        tic = np.arange(12.0).reshape(1, 3, 4)
        path = write_store(
            tmp_path / "s.zarr", {"ds_z0": make_table()}, images={"ds_z0_tic": [tic]}
        )
        store = open_store(path)

        np.testing.assert_array_equal(store.tic_image("ds_z0"), tic[0])
        assert store.tic_image("another_table") is None

    def test_a_preview_reads_the_coarsest_level_that_is_large_enough(self, tmp_path):
        fine = np.zeros((3, 400, 800), dtype=np.uint8)
        middle = np.ones((3, 100, 200), dtype=np.uint8)
        coarse = np.full((3, 25, 50), 2, dtype=np.uint8)
        path = write_store(
            tmp_path / "s.zarr",
            {"ds_z0": make_table()},
            images={"ds_optical": [fine, middle, coarse]},
        )

        preview = open_store(path).image_preview("ds_optical", 100)

        assert preview.shape == (3, 50, 100)  # the middle level, every 2nd pixel
        assert int(preview.max()) == 1

    def test_a_single_level_is_thinned_while_it_is_read(self, tmp_path):
        only = np.zeros((3, 300, 900), dtype=np.uint8)
        path = write_store(
            tmp_path / "s.zarr", {"ds_z0": make_table()}, images={"photo": [only]}
        )

        assert open_store(path).image_preview("photo", 100).shape == (3, 34, 100)


def test_the_store_name_drops_the_zarr_ending(tmp_path):
    path = write_store(tmp_path / "Brain 1.ZARR", {"t": make_table()})

    assert open_store(path).name == "Brain 1"
