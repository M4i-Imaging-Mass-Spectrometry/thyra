"""Every converted store records a fingerprint of its mass axis (issue #444, D42).

``msi_metadata.provenance.mass_axis`` holds the SHA-256 of ``var["mz"]``
(as little-endian float64 bytes) and its bin count. It is written whether
or not the conversion resampled, and every table of a conversion carries the
same value: the sibling tables' ``mz_index`` points into the main axis, so
they say which axis they are placed on rather than hash their own ``mz``.
"""

import hashlib
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pytest

spatialdata = pytest.importorskip("spatialdata")

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader  # noqa: E402
from thyra.convert import convert_msi  # noqa: E402
from thyra.converters.spatialdata import SpatialDataConverter  # noqa: E402
from thyra.metadata.schema import (  # noqa: E402
    check_store_var_conventions,
    read_msi_metadata_blocks,
)
from thyra.resampling.types import ResamplingConfig, ResamplingMethod  # noqa: E402
from thyra.utils.windows_paths import prepare_zarr_read_path  # noqa: E402

CONTINUOUS = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "fixtures"
    / "mobility_continuous.imzML"
)


def _var_mz(store: Path, table: str) -> np.ndarray:
    import zarr

    root = zarr.open_group(str(prepare_zarr_read_path(store)), mode="r")
    return np.asarray(root["tables"][table]["var"]["mz"])


def _expected(store: Path, table: str) -> Dict[str, Any]:
    mz = _var_mz(store, table)
    digest = hashlib.sha256(mz.astype("<f8").tobytes()).hexdigest()
    return {"sha256": digest, "n_bins": mz.size}


def _recorded(block: Dict[str, Any]) -> Dict[str, Any]:
    recorded = block["provenance"]["mass_axis"]
    return {"sha256": str(recorded["sha256"]), "n_bins": int(recorded["n_bins"])}


def _convert(tmp_path: Path, **kwargs: Any) -> Path:
    out = tmp_path / "out.zarr"
    assert convert_msi(
        str(CONTINUOUS), str(out), dataset_id="mob", pixel_size_um=10.0, **kwargs
    )
    return out


_RESAMPLED = ResamplingConfig(
    method=ResamplingMethod.NEAREST_NEIGHBOR,
    target_bins=400,
    min_mz=250.0,
    max_mz=650.0,
)


@pytest.mark.parametrize(
    "resampling_config", [_RESAMPLED, None], ids=["resampled", "no-resample"]
)
def test_every_table_carries_the_fingerprint_of_the_main_axis(
    tmp_path, resampling_config
):
    out = _convert(tmp_path, resampling_config=resampling_config)
    blocks = read_msi_metadata_blocks(out)
    assert set(blocks) == {"mob_z0", "mob_z0_mobility"}

    main = _expected(out, "mob_z0")
    if resampling_config is not None:
        assert main["n_bins"] == 400
    for name, block in blocks.items():
        assert _recorded(block) == main, name

    # The sibling's own var["mz"] repeats values split by mobility, so it
    # does not hash to the main axis; its block still names the main axis.
    assert _expected(out, "mob_z0_mobility") != main


def test_a_fresh_store_passes_the_fingerprint_check(tmp_path):
    out = _convert(tmp_path)
    # Without a recorded fingerprint the check has nothing to compare, so
    # a store without one would pass for the wrong reason.
    assert all(_recorded(block) for block in read_msi_metadata_blocks(out).values())
    issues = check_store_var_conventions(out)
    assert all(not table_issues for table_issues in issues.values()), issues


def test_a_volume_and_its_slices_share_one_fingerprint(tmp_path):
    config = MockMSIConfig(
        n_x=4, n_y=4, n_z=3, n_mz_bins=64, peaks_per_spectrum=(4, 8), sparsity=0.0
    )
    out = tmp_path / "volume.zarr"
    converter = SpatialDataConverter(
        MockMSIReader(config),
        out,
        handle_3d=False,
        dataset_id="vol",
        pixel_size_um=10.0,
    )
    assert converter.convert() is True

    blocks = read_msi_metadata_blocks(out)
    assert len(blocks) == 3
    values = {name: _recorded(block) for name, block in blocks.items()}
    assert len({tuple(v.items()) for v in values.values()}) == 1
    first = sorted(blocks)[0]
    assert values[first] == _expected(out, first)

    volume = tmp_path / "volume3d.zarr"
    converter = SpatialDataConverter(
        MockMSIReader(config),
        volume,
        handle_3d=True,
        dataset_id="vol",
        pixel_size_um=10.0,
        z_spacing_um=5.0,
    )
    assert converter.convert() is True
    volume_blocks = read_msi_metadata_blocks(volume)
    assert set(volume_blocks) == {"vol"}
    assert _recorded(volume_blocks["vol"]) == _expected(volume, "vol")
