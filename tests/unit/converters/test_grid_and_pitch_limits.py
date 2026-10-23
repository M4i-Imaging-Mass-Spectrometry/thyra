"""A pitch past what a float holds, and a grid past memory (issue #527).

F39: a finite ``--pixel-size`` of 1e308 um times the third pixel is
infinity, which reached ``obs`` and the pixel boxes at exit 0. F40: one
spectrum at (1000000, 1000000) beside one at (1, 1) asked the per-position
arrays for 10**12 positions and ended in a ``MemoryError`` traceback, and a
position of 2**31 or more failed an int32 cast in the imzML reader first.
Both are refused in one line now, before anything is allocated.
"""

import logging

import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.convert import convert_msi
from thyra.readers.imzml.imzml_reader import ImzMLReader
from thyra.resampling import memory_guard


def _imzml(path, coords):
    with ImzMLWriter(str(path), mode="processed") as writer:
        for c in coords:
            writer.addSpectrum(
                np.array([100.0, 300.0]), np.array([2.0, 3.0], dtype=np.float32), c
            )
    return path


def _refusal(thyra_logs, *args, **kwargs):
    with thyra_logs("thyra", logging.ERROR) as records:
        assert convert_msi(*args, **kwargs) is False
    return records.text


class TestAPitchPastFloat:
    def test_a_pixel_size_whose_far_edge_overflows_is_refused(
        self, tmp_path, thyra_logs
    ):
        src = _imzml(tmp_path / "row.imzML", [(x, 1, 1) for x in (1, 2, 3)])
        out = tmp_path / "row.zarr"
        text = _refusal(thyra_logs, str(src), str(out), pixel_size_um=1e308)
        assert "past the largest number a float holds" in text
        assert not out.exists()

    def test_a_large_pixel_size_that_stays_finite_converts(self, tmp_path):
        src = _imzml(tmp_path / "row.imzML", [(x, 1, 1) for x in (1, 2, 3)])
        assert convert_msi(str(src), str(tmp_path / "row.zarr"), pixel_size_um=1e300)

    def test_a_z_spacing_whose_far_edge_overflows_is_refused(
        self, tmp_path, thyra_logs
    ):
        src = _imzml(tmp_path / "vol.imzML", [(1, 1, z) for z in (1, 2, 3)])
        text = _refusal(
            thyra_logs,
            str(src),
            str(tmp_path / "vol.zarr"),
            handle_3d=True,
            pixel_size_um=10.0,
            z_spacing_um=1e308,
        )
        assert "z spacing of 1e+308 um across 3 pixels along z" in text


class TestAGridPastMemory:
    def test_one_far_spectrum_is_refused_naming_the_grid(
        self, tmp_path, thyra_logs, monkeypatch
    ):
        monkeypatch.setattr(memory_guard, "available_memory_gb", lambda: 16.0)
        src = _imzml(tmp_path / "far.imzML", [(1, 1, 1), (1_000_000, 1_000_000, 1)])
        out = tmp_path / "far.zarr"
        text = _refusal(thyra_logs, str(src), str(out), pixel_size_um=25.0)
        assert "1,000,000 x 1,000,000" in text
        assert "for 2 spectra" in text
        assert not out.exists()

    def test_a_position_past_int32_reaches_the_grid_check(
        self, tmp_path, thyra_logs, monkeypatch
    ):
        monkeypatch.setattr(memory_guard, "available_memory_gb", lambda: 16.0)
        src = _imzml(tmp_path / "wide.imzML", [(1, 1, 1), (2**31 + 5, 1, 1)])
        text = _refusal(
            thyra_logs, str(src), str(tmp_path / "w.zarr"), pixel_size_um=25
        )
        assert "the pixel grid would be 2,147,483,653 x 1" in text

    def test_the_reader_keeps_a_position_past_int32(self, tmp_path):
        src = _imzml(tmp_path / "wide.imzML", [(1, 1, 1), (2**31 + 5, 1, 1)])
        with ImzMLReader(src) as reader:
            coords = sorted(c for c, _, _ in reader.iter_spectra())
        assert coords == [(0, 0, 0), (2**31 + 4, 0, 0)]

    def test_a_grid_near_the_limit_is_warned_about(self, thyra_logs):
        with thyra_logs("thyra", logging.WARNING) as records:
            refusal = memory_guard.pixel_grid_refusal(
                (1000, 1000), 4, records_order=True, available_gb=0.06
            )
        assert refusal is None
        assert "1,000 x 1,000 = 1,000,000 positions" in records.text
