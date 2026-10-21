"""Every pixel polygon is its TIC cell, regions without an Area included (#520).

End to end on the committed synthetic timsTOF acquisition, split into one
region per row, with a ``.mis`` that has no teaching points: the fallback
D75 replaced. Before #520 the TIC was stretched over the box of all Areas
while each region's polygons used its own Area, so an Area drawn a step too
wide put cell and polygon 0.63 steps apart, and a region with no Area got
no polygon at all.
"""

from __future__ import annotations

import shutil
import sqlite3
from pathlib import Path
from typing import List

import numpy as np
import pytest

FIXTURE = Path(__file__).resolve().parents[2] / "data" / "fixtures" / "synthetic_tims.d"


def _acquisition(tmp_path: Path, areas: List[str]) -> Path:
    """The fixture with row y=10 as region 0 and row y=11 as region 1."""
    from thyra.readers.bruker.timstof.timstof_reader import BrukerReader
    from thyra.utils.bruker_exceptions import SDKError

    source = tmp_path / "run.d"
    shutil.copytree(FIXTURE, source)
    with sqlite3.connect(source / "analysis.tdf") as conn:
        conn.execute(
            "UPDATE MaldiFrameInfo SET RegionNumber = YIndexPos - 10, "
            "SpotName = 'R0' || (YIndexPos - 10) || 'X' || XIndexPos || "
            "'Y' || YIndexPos"
        )
    conn.close()
    for side in ("-wal", "-shm"):
        (source / f"analysis.tdf{side}").unlink(missing_ok=True)
    body = "".join(
        f'<Area Name="A{i}" Type="0"><Point>{a}</Point><Point>{b}</Point></Area>'
        for i, (a, b) in enumerate(area.split(";") for area in areas)
    )
    (tmp_path / "run.mis").write_text(f"<ImagingSequence>{body}</ImagingSequence>")
    try:
        BrukerReader(source).close()
    except (SDKError, OSError) as exc:  # the vendor library is not loadable here
        pytest.skip(f"Bruker library not loadable on this platform: {exc}")
    return source


@pytest.mark.parametrize(
    "areas",
    [
        ["20,20;80,40", "20,40;80,60"],
        ["20,20;100,40", "20,40;80,60"],
        ["20,20;80,40"],
    ],
    ids=["exact Areas", "Area 0 a step too wide", "no Area for region 1"],
)
def test_each_row_has_a_polygon_on_its_tic_cell(
    tmp_path: Path, areas: List[str]
) -> None:
    import spatialdata as sd
    from spatialdata.transformations import get_transformation

    from thyra.convert import convert_msi

    source = _acquisition(tmp_path, areas)
    output = tmp_path / "o.zarr"
    assert convert_msi(str(source), str(output), include_optical=False)

    sdata = sd.read_zarr(str(output))
    obs = sdata.tables["msi_dataset_z0"].obs
    shapes = sdata.shapes["msi_dataset_z0_pixels"]
    assert len(shapes) == len(obs) == 6
    matrix = get_transformation(sdata.images["msi_dataset_z0_tic"], "global")
    affine = matrix.to_affine_matrix(("x", "y"), ("x", "y"))
    cells = np.c_[obs["x"] + 0.5, obs["y"] + 0.5, np.ones(len(obs))] @ affine.T
    centres = shapes.geometry.loc[obs.index].centroid
    np.testing.assert_allclose(centres.x, cells[:, 0], atol=1e-9)
    np.testing.assert_allclose(centres.y, cells[:, 1], atol=1e-9)
