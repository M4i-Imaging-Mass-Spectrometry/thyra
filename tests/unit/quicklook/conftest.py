"""Fixtures shared by the tests of ``thyra ion-image`` and ``thyra report``."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.fixtures.hand_built_store import make_table, write_store


@pytest.fixture
def store_path(tmp_path: Path) -> Path:
    """A store with one table of spectra over the L-shaped region."""
    return write_store(tmp_path / "sample.zarr", {"msi_dataset_z0": make_table()})
