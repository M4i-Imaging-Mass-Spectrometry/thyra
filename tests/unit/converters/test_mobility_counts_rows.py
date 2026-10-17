"""The mobility siblings count the summed table's rows (issue #495).

The heatmap is a mean over the rows ``average_spectrum`` divides by (B1).
A multi-plane conversion does not warn that the other planes' spectra
were skipped (B3), and a 2D conversion plans tables for the planes the
source states, not the z span (L1). The shared-axis mobility table takes
the memory guard the grid route takes (J1).

Built on the committed continuous imzML with a mobility array: four
spectra, one shared m/z list of five values, then five intensities and
five mobilities per spectrum.
"""

from __future__ import annotations

import logging
import re
import shutil
from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

pytest.importorskip("spatialdata")

import spatialdata  # noqa: E402

from thyra.convert import convert_msi  # noqa: E402

_FIXTURE = Path(__file__).resolve().parents[2] / "data" / "fixtures"
_CONTINUOUS = _FIXTURE / "mobility_continuous"

#: Byte offset of spectrum ``i``'s intensities in the ``.ibd``: the 16-byte
#: UUID, the shared m/z list (5 float64), then 80 bytes per spectrum.
_INTENSITIES = 16 + 40


def _copy(tmp_path: Path, name: str, zero_last: bool = False) -> Path:
    """The continuous fixture, its last spectrum all zero when asked."""
    stem = tmp_path / name
    shutil.copy(_CONTINUOUS.with_suffix(".imzML"), stem.with_suffix(".imzML"))
    blob = bytearray(_CONTINUOUS.with_suffix(".ibd").read_bytes())
    if zero_last:
        start = _INTENSITIES + 3 * 80
        blob[start : start + 40] = bytes(40)
    stem.with_suffix(".ibd").write_bytes(bytes(blob))
    return stem.with_suffix(".imzML")


def _two_planes(tmp_path: Path, zs=(1, 1, 2, 2)) -> Path:
    """The continuous fixture with each spectrum given a ``position z``."""
    stem = tmp_path / "planes"
    planes = iter(zs)
    text = _CONTINUOUS.with_suffix(".imzML").read_text(encoding="utf-8")
    text = re.sub(
        r'<cvParam[^>]*name="position y"[^>]*/>',
        lambda m: m.group(0) + '<cvParam cvRef="IMS" accession="IMS:1000052" '
        f'name="position z" value="{next(planes)}" />',
        text,
    )
    stem.with_suffix(".imzML").write_text(text, encoding="utf-8")
    shutil.copy(_CONTINUOUS.with_suffix(".ibd"), stem.with_suffix(".ibd"))
    return stem.with_suffix(".imzML")


def _convert(src: Path, out: Path, **kwargs):
    assert convert_msi(str(src), str(out), dataset_id="ds", pixel_size_um=10, **kwargs)
    return spatialdata.read_zarr(out).tables


def _heatmap_total(table) -> float:
    return float(np.sum(table.uns["mobility_heatmap"]["counts"]))


class TestTheHeatmapCountsRows:
    """B1: an all-zero spectrum has no row, so it is no pixel either."""

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"resampling_config": None},
            {},
            {"resampling_config": None, "handle_3d": True},
        ],
        ids=["no-resample", "resampled", "volume"],
    )
    def test_the_ratio_is_one(self, tmp_path, kwargs):
        src = _copy(tmp_path, "zero", zero_last=True)
        tables = _convert(src, tmp_path / "out.zarr", **kwargs)
        table = tables["ds" if kwargs.get("handle_3d") else "ds_z0"]

        assert table.n_obs == 3
        assert table.uns["mobility_heatmap"]["current_ratio"] == pytest.approx(1.0)
        assert _heatmap_total(table) == pytest.approx(
            float(np.sum(table.uns["average_spectrum"]))
        )

    def test_a_file_with_every_row_is_unchanged(self, tmp_path):
        table = _convert(
            _copy(tmp_path, "control"), tmp_path / "out.zarr", resampling_config=None
        )["ds_z0"]

        assert table.n_obs == 4
        assert table.uns["mobility_heatmap"]["current_ratio"] == pytest.approx(1.0)

    def test_the_zero_spectrum_is_not_called_skipped(self, tmp_path, thyra_logs):
        src = _copy(tmp_path, "zero", zero_last=True)
        with thyra_logs("thyra", logging.WARNING) as logs:
            _convert(src, tmp_path / "out.zarr", resampling_config=None)

        messages = [r.getMessage() for r in logs]
        assert not [m for m in messages if "had no row" in m]
        assert not [m for m in messages if "mass-mobility heatmap holds" in m]


class TestEachPlaneCountsItsOwnSpectra:
    """B3: a multi-plane 2D conversion loses nothing and says nothing."""

    def test_no_plane_warns_about_the_others(self, tmp_path, thyra_logs):
        with thyra_logs("thyra", logging.WARNING) as logs:
            tables = _convert(
                _two_planes(tmp_path), tmp_path / "out.zarr", resampling_config=None
            )

        assert not [r for r in logs if "had no row" in r.getMessage()]
        for name in ("ds_z0", "ds_z0_mobility", "ds_z1", "ds_z1_mobility"):
            assert tables[name].n_obs == 2, name

    def test_a_plane_still_names_its_own_loss(self):
        from thyra.converters.spatialdata.mobility_table import loses_signal

        assert loses_signal((0, 0, 1), 1, [0.0, 3.0])
        assert not loses_signal((0, 0, 2), 1, [0.0, 3.0])
        assert not loses_signal((0, 0, 1), 1, [0.0, 0.0])
        assert loses_signal((0, 0, 7), None, [1.0])


class TestOnlyStatedPlanesGetATable:
    """L1: z 1 and z 120 are two planes, not 120."""

    def test_two_planes_far_apart(self, tmp_path, thyra_logs):
        src = tmp_path / "gap.imzML"
        with ImzMLWriter(str(src), mode="processed") as writer:
            for z in (1, 120):
                for x in (1, 2, 3):
                    writer.addSpectrum(
                        np.array([150.0, 400.0, 800.0]),
                        np.array([5.0, 1.0, 2.0]),
                        (x, 1, z),
                    )

        with thyra_logs("thyra", logging.WARNING) as logs:
            tables = _convert(src, tmp_path / "out.zarr")

        assert sorted(t for t in tables if t.startswith("ds")) == ["ds_z0", "ds_z119"]
        assert not [r for r in logs if "no position carries" in r.getMessage()]
        assert tables["ds_z0"].n_obs == tables["ds_z119"].n_obs == 3


class TestTheSharedAxisTableIsGuarded:
    """J1: the feature list is checked against free memory before it is built."""

    def _refusing(self, monkeypatch, free_gb: float):
        import thyra.converters.spatialdata.mobility_table as mobility_table

        monkeypatch.setattr(mobility_table, "available_memory_gb", lambda: free_gb)

    def test_a_list_too_wide_is_refused_with_its_lever(
        self, tmp_path, monkeypatch, thyra_logs
    ):
        self._refusing(monkeypatch, 1e-9)
        with thyra_logs("thyra", logging.WARNING) as logs:
            tables = _convert(
                _copy(tmp_path, "src"), tmp_path / "out.zarr", resampling_config=None
            )

        assert "ds_z0" in tables
        assert "ds_z0_mobility" not in tables
        refusals = [
            r.getMessage()
            for r in logs
            if r.getMessage().startswith("No mobility-resolved table")
        ]
        assert len(refusals) == 1
        assert "--no-mobility-table" in refusals[0]
        assert "(m/z, mobility) features" in refusals[0]

    def test_the_guard_is_asked_before_anything_is_allocated(self, monkeypatch):
        import thyra.converters.spatialdata.mobility_table as mobility_table

        asked = []

        def guard(n_features, *args, **kwargs):
            asked.append(n_features)
            return "too wide"

        def no_allocation(*args, **kwargs):  # pragma: no cover - must not run
            raise AssertionError("the assembly was built")

        class Reader:
            def get_shared_mobility_features(self):
                # One pair listed twice: the var holds two features, not three.
                return np.array([100.0, 100.0, 200.0]), np.array([1.0, 1.0, 1.1])

        monkeypatch.setattr(mobility_table, "var_ceiling_refusal", guard)
        monkeypatch.setattr(mobility_table, "CscAssembly", no_allocation)

        built = mobility_table._build_from_shared_axis(
            Reader(), lambda coords: 0, np.array([100.0, 200.0]), 1, Path(".")
        )

        assert built is None
        assert asked == [2]
