"""A timsTOF reader asked for spectra it cannot read says why.

A reader opened with ``metadata_only=True`` never loads the Bruker
library, and a closed reader has let go of the library and the database.
Either used to fail on ``None`` with an AttributeError that named neither
cause.
"""

import shutil
from pathlib import Path

import pytest

from tests.fixtures.fake_tims_sdk import FakeSDK
from thyra.core.msms import IsolationWindow
from thyra.readers.bruker.timstof.timstof_reader import BrukerReader
from thyra.utils.bruker_exceptions import DataError, SDKError

FIXTURE = Path(__file__).resolve().parents[2] / "data" / "fixtures" / "synthetic_tims.d"


@pytest.fixture
def tdf(tmp_path) -> Path:
    target = tmp_path / "copy.d"
    shutil.copytree(FIXTURE, target)
    return target


@pytest.fixture
def closed_reader(tdf):
    reader = BrukerReader(tdf, metadata_only=True)
    reader.sdk = FakeSDK()
    reader.handle = 1
    reader.close()
    return reader


def test_metadata_only_reader_names_the_mode(tdf):
    reader = BrukerReader(tdf, metadata_only=True)
    try:
        with pytest.raises(SDKError, match="opened in metadata-only mode"):
            reader._read_frame_spectrum(1)
    finally:
        reader.close()


def test_closed_reader_says_it_is_closed(closed_reader):
    with pytest.raises(DataError, match="the reader has been closed"):
        list(closed_reader.iter_spectra())


def test_closed_reader_mass_axis_says_it_is_closed(closed_reader):
    with pytest.raises(DataError, match="the reader has been closed"):
        closed_reader.get_common_mass_axis()


def test_mobility_guard_tells_closed_from_metadata_only(closed_reader):
    with pytest.raises(SDKError, match="Reading mobility spectra .* has been closed"):
        list(closed_reader.iter_mobility_spectra())


def test_a_window_without_a_scan_range_is_named(tdf):
    reader = BrukerReader(tdf, metadata_only=True)
    try:
        windows = (
            IsolationWindow(target=500.0, scan_begin=0, scan_end=10),
            IsolationWindow(target=600.0),
        )
        with pytest.raises(ValueError, match="window 1 carries no mobility scan"):
            reader._precursor_scan_map(windows)
    finally:
        reader.close()
