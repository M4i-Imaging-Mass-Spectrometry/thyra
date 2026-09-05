"""The mobility-resolved sibling table, end to end from an imzML with a mobility array.

``mobility_continuous.imzML`` shares one (m/z, 1/K0) feature list across
four pixels; two m/z values repeat because mobility splits them. The summed
MSI table must merge those, the sibling table must keep them apart, and the
two must agree on rows, region and provenance. ``mobility_processed.imzML``
has mobility per pixel and gets the summed table only.
"""

import json
from pathlib import Path

import numpy as np
import pytest

spatialdata = pytest.importorskip("spatialdata")

from thyra.convert import convert_msi  # noqa: E402
from thyra.metadata.schema import (  # noqa: E402
    check_store_var_conventions,
    read_msi_metadata_blocks,
)
from thyra.resampling.types import ResamplingConfig, ResamplingMethod  # noqa: E402
from thyra.utils.windows_paths import prepare_zarr_read_path  # noqa: E402

FIXTURES = Path(__file__).resolve().parents[2] / "data" / "fixtures"
CONTINUOUS = FIXTURES / "mobility_continuous.imzML"
PROCESSED = FIXTURES / "mobility_processed.imzML"

SORTED_MZ = [300.0, 300.0, 450.5, 600.25, 600.25]
SORTED_K0 = [0.95, 1.10, 1.02, 1.20, 1.35]
# Per pixel intensities in (mz, mobility) order; the written file lists the
# 600.25 pair mobility-descending, so 20 and 2 swap places here.
SORTED_ROWS = {
    (0, 0): [10.0, 1.0, 5.0, 20.0, 2.0],
    (1, 0): [11.0, 2.0, 6.0, 21.0, 3.0],
    (0, 1): [12.0, 3.0, 7.0, 22.0, 4.0],
    (1, 1): [13.0, 4.0, 8.0, 23.0, 0.0],
}
SUMMED_ROWS = {
    (0, 0): [11.0, 5.0, 22.0],
    (1, 0): [13.0, 6.0, 24.0],
    (0, 1): [15.0, 7.0, 26.0],
    (1, 1): [17.0, 8.0, 23.0],
}


def _convert(fixture: Path, tmp_path: Path, **kwargs):
    out = tmp_path / f"{fixture.stem}.zarr"
    ok = convert_msi(
        str(fixture), str(out), dataset_id="mob", pixel_size_um=10.0, **kwargs
    )
    assert ok, "conversion reported failure"
    # The mobility table's deepest uns key is the longest path Thyra writes;
    # under pytest's temp directory that passes the Windows limit, which the
    # converter handles on the write side and this helper on the read side.
    readable = prepare_zarr_read_path(out)
    return readable, spatialdata.read_zarr(readable)


def _row(table, x: int, y: int) -> np.ndarray:
    obs = table.obs
    mask = (obs["x"].to_numpy().astype(int) == x) & (
        obs["y"].to_numpy().astype(int) == y
    )
    rows = np.flatnonzero(mask)
    assert rows.size == 1, f"pixel ({x}, {y}) has {rows.size} rows"
    X = table.X[rows[0]]
    return (
        np.asarray(X.toarray()).ravel()
        if hasattr(X, "toarray")
        else np.asarray(X).ravel()
    )


@pytest.mark.parametrize("streaming", [False, True], ids=["in-memory", "streaming"])
class TestContinuousExport:
    def test_both_tables_are_written_and_agree_on_rows(self, tmp_path, streaming):
        _, sdata = _convert(CONTINUOUS, tmp_path, streaming=streaming)
        assert set(sdata.tables) == {"mob_z0", "mob_z0_mobility"}
        summed = sdata.tables["mob_z0"]
        mobility = sdata.tables["mob_z0_mobility"]
        assert summed.n_obs == mobility.n_obs == 4
        assert list(summed.obs.index) == list(mobility.obs.index)
        assert set(mobility.obs["region"].astype(str)) == {"mob_z0_pixels"}
        assert "mob_z0_pixels" in sdata.shapes

    def test_mobility_var_is_the_sorted_feature_pairs(self, tmp_path, streaming):
        _, sdata = _convert(CONTINUOUS, tmp_path, streaming=streaming)
        var = sdata.tables["mob_z0_mobility"].var
        np.testing.assert_array_equal(var["mz"].to_numpy(), SORTED_MZ)
        np.testing.assert_array_equal(var["mobility"].to_numpy(), SORTED_K0)
        np.testing.assert_array_equal(var["mz_index"].to_numpy(), [0, 0, 1, 2, 2])
        np.testing.assert_array_equal(var["mobility_index"].to_numpy(), [0, 2, 1, 3, 4])
        assert list(var.index) == [
            "mz0_im0",
            "mz0_im2",
            "mz1_im1",
            "mz2_im3",
            "mz2_im4",
        ]
        # The summed table's axis is strictly increasing and unique.
        np.testing.assert_array_equal(
            sdata.tables["mob_z0"].var["mz"].to_numpy(), [300.0, 450.5, 600.25]
        )
        assert "mobility" not in sdata.tables["mob_z0"].var.columns

    def test_intensities_land_on_their_features(self, tmp_path, streaming):
        _, sdata = _convert(CONTINUOUS, tmp_path, streaming=streaming)
        mobility = sdata.tables["mob_z0_mobility"]
        for (x, y), expected in SORTED_ROWS.items():
            np.testing.assert_array_equal(_row(mobility, x, y), expected)
        # The zero on the last pixel is absent, not stored.
        assert mobility.X.nnz == 19

    def test_summed_table_merges_features_split_by_mobility(self, tmp_path, streaming):
        _, sdata = _convert(CONTINUOUS, tmp_path, streaming=streaming)
        summed = sdata.tables["mob_z0"]
        for (x, y), expected in SUMMED_ROWS.items():
            np.testing.assert_array_equal(_row(summed, x, y), expected)

    def test_uns_links_the_two_tables(self, tmp_path, streaming):
        out, sdata = _convert(CONTINUOUS, tmp_path, streaming=streaming)
        axis = sdata.tables["mob_z0"].uns["mobility_axis"]
        assert axis["present"] is True or axis["present"] == 1
        assert axis["type_accession"] == "MS:1002815"
        assert axis["unit_accession"] == "MS:1002814"
        assert axis["resolved_table"] == "mob_z0_mobility"
        np.testing.assert_array_equal(
            np.asarray(axis["values"]), [0.95, 1.10, 1.02, 1.35, 1.20]
        )
        assert axis["n_scans"] == 5
        np.testing.assert_array_equal(np.asarray(axis["acq_range"]), [0.95, 1.35])

        # The heatmap sits on the summed table only: the sibling holds the
        # very data it summarises.
        heatmap = sdata.tables["mob_z0"].uns["mobility_heatmap"]
        counts = np.asarray(heatmap["counts"])
        assert counts.shape == (3, 256) and counts.dtype == np.float32
        np.testing.assert_allclose(
            counts.sum(axis=1),
            np.asarray(sdata.tables["mob_z0"].uns["average_spectrum"]),
            rtol=1e-6,
        )

        mob_uns = sdata.tables["mob_z0_mobility"].uns
        assert "mobility_heatmap" not in mob_uns
        feature_axis = mob_uns["feature_axis"]
        assert json.loads(feature_axis["dims"]) == ["mz", "mobility"]
        assert feature_axis["summed_table"] == "mob_z0"
        assert mob_uns["mobility_axis"]["resolved_table"] == "mob_z0_mobility"

        blocks = read_msi_metadata_blocks(out)
        assert set(blocks) == {"mob_z0", "mob_z0_mobility"}
        for block in blocks.values():
            mobility = block["ms_analysis"]["ion_mobility"]
            assert mobility["present"] is True
            assert mobility["separation_term"]["accession"] == "MS:1002815"
            assert mobility["resolved_table"] == "mob_z0_mobility"
            assert "grid" not in mobility

    def test_store_validates_under_both_var_contracts(self, tmp_path, streaming):
        out, _ = _convert(CONTINUOUS, tmp_path, streaming=streaming)
        issues = check_store_var_conventions(out)
        assert set(issues) == {"mob_z0", "mob_z0_mobility"}
        assert all(not table_issues for table_issues in issues.values()), issues


class TestResampledSummedTable:
    def test_features_are_summed_into_their_bins(self, tmp_path):
        config = ResamplingConfig(
            method=ResamplingMethod.NEAREST_NEIGHBOR,
            target_bins=400,
            min_mz=250.0,
            max_mz=650.0,
        )
        _, sdata = _convert(CONTINUOUS, tmp_path, resampling_config=config)
        summed = sdata.tables["mob_z0"]
        axis = summed.var["mz"].to_numpy()
        assert axis.size == 400
        row = _row(summed, 0, 0)
        assert row[np.argmin(np.abs(axis - 300.0))] == pytest.approx(11.0)
        assert row[np.argmin(np.abs(axis - 600.25))] == pytest.approx(22.0)
        assert row.sum() == pytest.approx(38.0)

        mobility = sdata.tables["mob_z0_mobility"]
        np.testing.assert_array_equal(mobility.var["mz"].to_numpy(), SORTED_MZ)
        mz_index = mobility.var["mz_index"].to_numpy()
        assert np.all(np.abs(axis[mz_index] - SORTED_MZ) <= np.diff(axis).max())
        np.testing.assert_array_equal(_row(mobility, 0, 0), SORTED_ROWS[(0, 0)])


class TestOptOut:
    def test_no_mobility_table_keeps_the_axis_description(self, tmp_path):
        _, sdata = _convert(CONTINUOUS, tmp_path, write_mobility_table=False)
        assert set(sdata.tables) == {"mob_z0"}
        axis = sdata.tables["mob_z0"].uns["mobility_axis"]
        assert axis["type_accession"] == "MS:1002815"
        assert "resolved_table" not in axis


@pytest.mark.parametrize("streaming", [False, True], ids=["in-memory", "streaming"])
class TestProcessedExport:
    def test_only_the_summed_table_with_repeats_merged(self, tmp_path, streaming):
        out, sdata = _convert(PROCESSED, tmp_path, streaming=streaming)
        assert set(sdata.tables) == {"mob_z0"}
        summed = sdata.tables["mob_z0"]
        axis = summed.var["mz"].to_numpy()
        assert axis.size == 12 and np.all(np.diff(axis) > 0)
        coords = [(0, 0), (1, 0), (0, 1), (1, 1)]
        for i, (x, y) in enumerate(coords):
            row = _row(summed, x, y)
            assert row[np.flatnonzero(axis == 100.0 + i)[0]] == pytest.approx(
                30.0 * (i + 1)
            )
            assert row[np.flatnonzero(axis == 400.0 + i)[0]] == pytest.approx(
                40.0 * (i + 1)
            )
            assert row.sum() == pytest.approx(100.0 * (i + 1))
        axis_block = summed.uns["mobility_axis"]
        assert axis_block["type_accession"] == "MS:1002815"
        assert "values" not in axis_block
        assert "resolved_table" not in axis_block
        # Per-pixel mobility with no shared axis gives no range to bin over.
        assert "mobility_heatmap" not in summed.uns
        block = read_msi_metadata_blocks(out)["mob_z0"]
        assert block["ms_analysis"]["ion_mobility"]["present"] is True
        assert check_store_var_conventions(out) == {"mob_z0": []}
