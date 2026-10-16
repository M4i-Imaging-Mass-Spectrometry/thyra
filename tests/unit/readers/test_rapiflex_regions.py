"""Region map of a multi-region Rapiflex acquisition (issue #448).

The poslog names each spot ``R{region}X{x}Y{y}`` in absolute raster units;
the ``.dat`` offset table places the spectra from the header origin. The
synthetic acquisition here has a nonzero origin and a first spot one column
past it, so a map keyed on the smallest poslog value would be one column
off, as it would on real data.
"""

import logging
import struct

import numpy as np
import pytest

from thyra.readers.bruker.rapiflex import RapiflexReader

N_POINTS = 50
ORIGIN = (100, 40)
WIDTH, HEIGHT = 5, 2


def _write(folder, cells, *, areas=("left", "right"), extra_poslog=()):
    """Write a Rapiflex acquisition with a spectrum in each of ``cells``.

    Args:
        folder: Directory to create.
        cells: ``{(x, y): region}``, 0-based cells relative to the origin.
        areas: ``.mis`` Area names, in order.
        extra_poslog: Further ``(x, y, region)`` poslog spots with no
            spectrum.
    """
    folder.mkdir()
    with open(folder / "sample_poslog.txt", "w") as f:
        f.write("#Timestamp Pos X Y Z\n")
        f.write("2023-01-01 12:00:00.000 __ 0.0 0.0 0.0\n")
        spots = [(x, y, r) for (x, y), r in cells.items()] + list(extra_poslog)
        for x, y, region in spots:
            pos_id = f"R{region:02d}X{ORIGIN[0] + x}Y{ORIGIN[1] + y}"
            f.write(f"2023-01-01 12:00:00.000 {pos_id} {x * 20.0} {y * 20.0} 0.0\n")
    with open(folder / "sample_info.txt", "w") as f:
        f.write(f"Number of Spots: {len(cells)}\n")
        f.write("Mass Start: 100.0\nMass End: 500.0\n")
        f.write(f"DataPoints: {N_POINTS}\nRaster: 20,20\n")
    with open(folder / "sample.mis", "w") as f:
        f.write('<?xml version="1.0"?>\n<ImagingSequence>\n')
        f.write("  <Raster>20,20</Raster>\n")
        for i, name in enumerate(areas):
            f.write(
                f'  <Area Name="{name}" Type="0">'
                f"<Point>{i * 100},0</Point><Point>{i * 100 + 50},50</Point>"
                "</Area>\n"
            )
        f.write("</ImagingSequence>\n")

    n_raster = WIDTH * HEIGHT
    data_start = 48 + n_raster * 4
    offsets = np.zeros(n_raster, dtype=np.uint32)
    order = sorted(cells, key=lambda c: c[1] * WIDTH + c[0])
    for k, (x, y) in enumerate(order):
        offsets[y * WIDTH + x] = data_start + k * N_POINTS * 4
    with open(folder / "sample.dat", "wb") as f:
        f.write(
            struct.pack(
                "<12I",
                48,
                256,
                ORIGIN[0],
                ORIGIN[1],
                WIDTH,
                HEIGHT,
                N_POINTS,
                0,
                0,
                0,
                0,
                0,
            )
        )
        f.write(offsets.tobytes())
        for _ in order:
            f.write(np.ones(N_POINTS, dtype=np.float32).tobytes())
    return folder


def _two_regions():
    """Region 0 in columns 1-2, region 1 in column 4; column 0 is empty."""
    cells = {(x, y): 0 for x in (1, 2) for y in range(HEIGHT)}
    cells.update({(4, y): 1 for y in range(HEIGHT)})
    return cells


def test_two_regions_give_a_map_keyed_from_the_header_origin(tmp_path):
    cells = _two_regions()
    reader = RapiflexReader(_write(tmp_path / "rf", cells))
    try:
        assert reader.get_region_map() == cells
        coords = {c[:2] for c, _, _ in reader.iter_spectra()}
        assert coords == set(cells)
    finally:
        reader.close()


def test_region_info_has_counts_bounds_and_area_names(tmp_path):
    reader = RapiflexReader(_write(tmp_path / "rf", _two_regions()))
    try:
        info = sorted(reader.get_region_info(), key=lambda e: e["region_number"])
    finally:
        reader.close()
    assert info == [
        {"region_number": 0, "n_spectra": 4, "bounds": (1, 0, 2, 1), "name": "left"},
        {"region_number": 1, "n_spectra": 2, "bounds": (4, 0, 4, 1), "name": "right"},
    ]


def test_one_region_is_a_map_with_its_poslog_number(tmp_path):
    """A single-region store keeps the poslog's number, not 1 (D58)."""
    cells = {(x, y): 2 for x in range(WIDTH) for y in range(HEIGHT)}
    reader = RapiflexReader(_write(tmp_path / "rf", cells))
    try:
        assert reader.get_region_map() == cells
        (entry,) = reader.get_region_info()
        assert entry["region_number"] == 2
        assert entry["n_spectra"] == len(cells)
        assert entry["bounds"] == (0, 0, WIDTH - 1, HEIGHT - 1)
    finally:
        reader.close()


def test_a_spectrum_without_a_poslog_spot_gives_no_map(tmp_path, thyra_logs):
    """A partial map would write region -1 for the missing pixels."""
    cells = _two_regions()
    folder = _write(tmp_path / "rf", cells)
    lines = (folder / "sample_poslog.txt").read_text().splitlines()
    (folder / "sample_poslog.txt").write_text("\n".join(lines[:-1]) + "\n")
    reader = RapiflexReader(folder)
    try:
        with thyra_logs("thyra.readers", logging.WARNING) as records:
            assert reader.get_region_map() is None
        assert reader.get_region_info() is None
    finally:
        reader.close()
    assert any("no poslog spot" in r.getMessage() for r in records)


def test_a_poslog_spot_without_a_spectrum_is_left_out(tmp_path):
    cells = _two_regions()
    reader = RapiflexReader(_write(tmp_path / "rf", cells, extra_poslog=[(0, 0, 0)]))
    try:
        assert reader.get_region_map() == cells
    finally:
        reader.close()


@pytest.mark.parametrize("areas", [(), ("only",)])
def test_a_region_without_an_area_has_no_name(tmp_path, areas):
    reader = RapiflexReader(_write(tmp_path / "rf", _two_regions(), areas=areas))
    try:
        info = {e["region_number"]: e for e in reader.get_region_info()}
    finally:
        reader.close()
    assert "name" not in info[1]


def test_a_region_without_spectra_is_not_a_second_region(tmp_path):
    """Poslog spots of a region with no spectra leave a one-region store."""
    cells = {(x, y): 0 for x in (1, 2) for y in range(HEIGHT)}
    reader = RapiflexReader(
        _write(tmp_path / "rf", cells, extra_poslog=[(4, 0, 1), (4, 1, 1)])
    )
    try:
        assert reader.get_region_map() == cells
        assert [e["region_number"] for e in reader.get_region_info()] == [0]
    finally:
        reader.close()


def test_a_conflict_on_an_empty_cell_is_ignored(tmp_path):
    """Two regions on a cell that holds no spectrum place no pixel."""
    cells = _two_regions()
    reader = RapiflexReader(
        _write(tmp_path / "rf", cells, extra_poslog=[(0, 0, 0), (0, 0, 1)])
    )
    try:
        assert reader.get_region_map() == cells
    finally:
        reader.close()


def test_a_spectrum_under_two_regions_gives_no_map(tmp_path):
    cells = _two_regions()
    reader = RapiflexReader(_write(tmp_path / "rf", cells, extra_poslog=[(1, 0, 1)]))
    try:
        assert reader.get_region_map() is None
    finally:
        reader.close()


def test_a_closed_reader_refuses(tmp_path):
    reader = RapiflexReader(_write(tmp_path / "rf", _two_regions()))
    reader.close()
    with pytest.raises(RuntimeError, match="closed"):
        reader.get_region_map()
