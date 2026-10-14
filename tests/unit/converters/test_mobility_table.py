"""The mobility-resolved sibling table, end to end from an imzML with a mobility array.

``mobility_continuous.imzML`` shares one (m/z, 1/K0) feature list across
four pixels; two m/z values repeat because mobility splits them. The summed
MSI table must merge those, the sibling table must keep them apart, and the
two must agree on rows, region and provenance. ``mobility_processed.imzML``
has mobility per pixel: by default it gets the summed table and the
mass-mobility heatmap over the values its pixels hold, and with
``mobility_grid`` the same sibling table, binned onto the heatmap's channels.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pytest

spatialdata = pytest.importorskip("spatialdata")

from thyra.convert import convert_msi  # noqa: E402
from thyra.metadata.schema import (  # noqa: E402
    check_store_var_conventions,
    read_msi_metadata_blocks,
)
from thyra.resampling.mobility_grid import MOBILITY_CHANNELS  # noqa: E402
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

    def test_ccs_is_written_for_an_assumed_charge_of_one(self, tmp_path, streaming):
        out, sdata = _convert(CONTINUOUS, tmp_path, streaming=streaming)
        mobility = sdata.tables["mob_z0_mobility"]
        mz, k0 = np.array(SORTED_MZ), np.array(SORTED_K0)
        reduced = mz * 28.0134 / (mz + 28.0134)
        expected = 18509.8632163405 * k0 / np.sqrt(reduced * 305.0)
        np.testing.assert_allclose(mobility.var["ccs"].to_numpy(), expected, rtol=1e-6)
        # Two isomers at one m/z: the later one in mobility is the larger.
        ccs = mobility.var["ccs"].to_numpy()
        assert ccs[1] > ccs[0] and ccs[4] > ccs[3]
        assert int(mobility.uns["ccs"]["charge_assumed"]) == 1
        assert "ccs" not in sdata.tables["mob_z0"].var.columns
        assert check_store_var_conventions(out)["mob_z0_mobility"] == []

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
        block = read_msi_metadata_blocks(out)["mob_z0"]
        assert block["ms_analysis"]["ion_mobility"]["present"] is True
        assert check_store_var_conventions(out) == {"mob_z0": []}

    def test_the_heatmap_bins_over_the_values_the_pixels_hold(
        self, tmp_path, streaming
    ):
        """No shared axis and no declared range: the pixels' own values decide."""
        _, sdata = _convert(PROCESSED, tmp_path, streaming=streaming)
        heat = sdata.tables["mob_z0"].uns["mobility_heatmap"]
        edges = np.asarray(heat["mobility_edges"])
        assert edges.size == MOBILITY_CHANNELS + 1
        assert edges[0] == pytest.approx(0.8) and edges[-1] == pytest.approx(1.1)
        counts = np.asarray(heat["counts"])
        # The mean over four pixels holding 100, 200, 300 and 400.
        assert float(counts.sum()) == pytest.approx(250.0)
        assert float(heat["current_ratio"]) == pytest.approx(1.0)
        # Each 1/K0 lands in its own channel: 0.8 in the first, 1.1 (the
        # upper edge) in the last.
        per_channel = counts.sum(axis=0)
        expected = {0: 25.0, 85: 75.0, 170: 100.0, MOBILITY_CHANNELS - 1: 50.0}
        assert set(np.flatnonzero(per_channel).tolist()) == set(expected)
        for channel, mean in expected.items():
            assert per_channel[channel] == pytest.approx(mean)

    def test_no_table_without_the_grid_and_the_log_says_how(
        self, tmp_path, streaming, thyra_logs
    ):
        with thyra_logs("thyra.converters.spatialdata.sibling_tables") as records:
            _, sdata = _convert(PROCESSED, tmp_path, streaming=streaming)
        assert "mob_z0_mobility" not in sdata.tables
        assert any("--mobility-grid" in r.getMessage() for r in records)


#: Each fixture point as (m/z base index, 1/K0, intensity at pixel 0); the
#: intensity scales with the pixel, as the fixture builder writes it.
PROCESSED_POINTS = [(0, 0.8, 10.0), (0, 1.1, 20.0), (1, 0.9, 30.0), (2, 1.0, 40.0)]
PROCESSED_MZ_BASES = (100.0, 250.5, 400.0)


@pytest.mark.parametrize("streaming", [False, True], ids=["in-memory", "streaming"])
class TestProcessedExportOnAGrid:
    """``mobility_grid`` on a processed imzML: the contract the TDF grid keeps."""

    def test_the_grid_table_is_written_on_the_heatmap_channels(
        self, tmp_path, streaming
    ):
        out, sdata = _convert(
            PROCESSED, tmp_path, streaming=streaming, mobility_grid=True
        )
        assert set(sdata.tables) == {"mob_z0", "mob_z0_mobility"}
        summed = sdata.tables["mob_z0"]
        block = sdata.tables["mob_z0_mobility"].uns["mobility_grid"]
        assert str(block["law"]) == "linear"
        assert int(block["n_channels"]) == MOBILITY_CHANNELS
        assert float(block["lower"]) == pytest.approx(0.8)
        assert float(block["upper"]) == pytest.approx(1.1)
        np.testing.assert_array_equal(
            np.asarray(block["edges"]),
            np.asarray(summed.uns["mobility_heatmap"]["mobility_edges"]),
        )
        assert summed.uns["mobility_axis"]["resolved_table"] == "mob_z0_mobility"
        assert check_store_var_conventions(out) == {
            "mob_z0": [],
            "mob_z0_mobility": [],
        }

    def test_the_grid_table_carries_ccs_like_every_1_k0_table(
        self, tmp_path, streaming
    ):
        """The fixture states 1/K0 in MS:1002814, so D47 applies (D46)."""
        _, sdata = _convert(
            PROCESSED, tmp_path, streaming=streaming, mobility_grid=True
        )
        grid = sdata.tables["mob_z0_mobility"]
        ccs = grid.var["ccs"].to_numpy()
        assert ccs.shape == (grid.n_vars,) and np.all(np.isfinite(ccs))
        assert int(grid.uns["ccs"]["charge_assumed"]) == 1
        # Same m/z: the higher-mobility channel has the larger cross section.
        mz = grid.var["mz"].to_numpy()
        at_100 = np.flatnonzero(mz == 100.0)
        order = np.argsort(grid.var["mobility"].to_numpy()[at_100])
        assert np.all(np.diff(ccs[at_100][order]) > 0)

    def test_each_point_lands_in_its_cell(self, tmp_path, streaming):
        _, sdata = _convert(
            PROCESSED, tmp_path, streaming=streaming, mobility_grid=True
        )
        grid = sdata.tables["mob_z0_mobility"]
        mz = grid.var["mz"].to_numpy()
        channel = grid.var["mobility_index"].to_numpy().astype(int)
        edges = np.asarray(grid.uns["mobility_grid"]["edges"])
        for i, (x, y) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)]):
            row = _row(grid, x, y)
            assert row.sum() == pytest.approx(100.0 * (i + 1))
            for base, k0, intensity in PROCESSED_POINTS:
                want = min(
                    int(np.searchsorted(edges, k0, side="right")) - 1,
                    MOBILITY_CHANNELS - 1,
                )
                cell = np.flatnonzero(
                    (mz == PROCESSED_MZ_BASES[base] + i) & (channel == want)
                )
                assert cell.size == 1
                assert row[cell[0]] == pytest.approx(intensity * (i + 1))

    def test_the_marginal_reproduces_the_summed_table(self, tmp_path, streaming):
        _, sdata = _convert(
            PROCESSED, tmp_path, streaming=streaming, mobility_grid=True
        )
        block = sdata.tables["mob_z0_mobility"].uns["mobility_marginal"]
        assert str(block["summed_table"]) == "mob_z0"
        assert float(block["current_ratio"]) == pytest.approx(1.0)
        assert float(block["current_ratio_pixel_min"]) == pytest.approx(1.0)
        assert float(block["current_ratio_pixel_max"]) == pytest.approx(1.0)

    def test_the_marginal_holds_under_the_cli_resampling(self, tmp_path, streaming):
        """The CLI resamples by default; the file's edge samples stay in.

        The laid axis puts its first and last bin centres inside the
        file's m/z range, so pixel 0's m/z 100.0 points and pixel 3's
        m/z 403.0 point sit within half a bin past the outer centres.
        The summed table keeps them, so the grid and the heatmap must.
        """
        cli_default = {
            "method": "auto",
            "axis_type": "auto",
            "reference_mz": 1000.0,
        }
        _, sdata = _convert(
            PROCESSED,
            tmp_path,
            streaming=streaming,
            mobility_grid=True,
            resampling_config=cli_default,
        )
        summed = sdata.tables["mob_z0"]
        # The premise: the laid axis's end centres sit inside the data.
        axis = summed.var["mz"].to_numpy()
        assert axis[0] > 100.0 and axis[-1] < 403.0
        assert float(summed.X.sum()) == pytest.approx(1000.0)
        block = sdata.tables["mob_z0_mobility"].uns["mobility_marginal"]
        assert float(block["current_ratio"]) == pytest.approx(1.0)
        assert float(block["current_ratio_pixel_min"]) == pytest.approx(1.0)
        assert float(block["current_ratio_pixel_max"]) == pytest.approx(1.0)
        grid = sdata.tables["mob_z0_mobility"]
        for i, (x, y) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)]):
            assert _row(grid, x, y).sum() == pytest.approx(100.0 * (i + 1))
        heatmap = summed.uns["mobility_heatmap"]
        assert float(heatmap["current_ratio"]) == pytest.approx(1.0)

    def test_a_stated_range_overrides_the_values(self, tmp_path, streaming):
        _, sdata = _convert(
            PROCESSED,
            tmp_path,
            streaming=streaming,
            mobility_grid=True,
            mobility_bins=3,
            mobility_min=0.7,
            mobility_max=1.3,
        )
        block = sdata.tables["mob_z0_mobility"].uns["mobility_grid"]
        assert int(block["n_channels"]) == 3
        np.testing.assert_allclose(np.asarray(block["edges"]), [0.7, 0.9, 1.1, 1.3])

    def test_the_grid_needs_no_heatmap(self, tmp_path, streaming):
        """The range is the reader's, not something the heatmap's pass leaves."""
        _, sdata = _convert(
            PROCESSED,
            tmp_path,
            streaming=streaming,
            mobility_grid=True,
            mobility_heatmap=False,
        )
        assert "mobility_heatmap" not in sdata.tables["mob_z0"].uns
        block = sdata.tables["mob_z0_mobility"].uns["mobility_grid"]
        assert float(block["lower"]) == pytest.approx(0.8)
        assert float(block["upper"]) == pytest.approx(1.1)


class TestAnUnreadableRange:
    def test_the_conversion_keeps_its_summed_table(self, tmp_path, monkeypatch):
        """No range costs the heatmap and the grid, not the conversion."""
        from thyra.errors import ConversionRefused
        from thyra.readers.imzml.imzml_reader import ImzMLReader

        def _refuse(self):
            raise ConversionRefused("ion mobility array truncated in the .ibd")

        monkeypatch.setattr(ImzMLReader, "_scan_mobility_span", _refuse)
        _, sdata = _convert(PROCESSED, tmp_path, mobility_grid=True)
        assert set(sdata.tables) == {"mob_z0"}
        summed = sdata.tables["mob_z0"]
        assert "mobility_heatmap" not in summed.uns
        assert "resolved_table" not in summed.uns["mobility_axis"]
        assert summed.X.sum() == pytest.approx(1000.0)


class TestTheGridOnASharedAxis:
    def test_the_grid_is_ignored_and_said_to_be(self, tmp_path, thyra_logs):
        """A continuous export is read off its features; the flag changes nothing."""
        with thyra_logs(
            "thyra.converters.spatialdata.sibling_tables", logging.WARNING
        ) as records:
            _, sdata = _convert(CONTINUOUS, tmp_path, mobility_grid=True)
        assert "mobility_grid" not in sdata.tables["mob_z0_mobility"].uns
        messages = [r.getMessage() for r in records]
        assert sum("--mobility-grid is ignored" in m for m in messages) == 1
