"""The TDF reader end to end through the real Bruker library.

``tests/data/fixtures/synthetic_tims.d`` is a hand-written TIMS acquisition
(see ``build_tdf_fixture.py``); ``synthetic_tims_expected.json`` is what was
written into it. The library is bundled for Windows and Linux; anywhere it
cannot be loaded these tests skip rather than fail.

A second group runs only against a real acquisition named by
``THYRA_BRUKER_TDF_DATASET`` and checks the reader against the database's own
per-frame ``SummedIntensities``.
"""

from __future__ import annotations

import json
import os
import sqlite3
from pathlib import Path
from typing import Dict

import numpy as np
import pytest

from thyra.readers.bruker.timstof.timstof_reader import BrukerReader
from thyra.utils.bruker_exceptions import SDKError

pytestmark = pytest.mark.integration

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "fixtures" / "synthetic_tims.d"
EXPECTED = FIXTURE.with_name("synthetic_tims_expected.json")


def _open(mode: str, path: Path = FIXTURE) -> BrukerReader:
    try:
        return BrukerReader(path, tdf_spectrum=mode)
    except (SDKError, OSError) as exc:  # the vendor library is not loadable here
        pytest.skip(f"Bruker library not loadable on this platform: {exc}")


@pytest.fixture(scope="module")
def expected() -> Dict:
    return json.loads(EXPECTED.read_text(encoding="utf-8"))


def _spectra_by_frame(reader: BrukerReader, expected: Dict) -> Dict[int, tuple]:
    """Map frame id -> (mz, intensity) using the fixture's grid."""
    offsets = reader.get_essential_metadata().coordinate_offsets or (0, 0, 0)
    by_coord = {
        (f["x"] - offsets[0], f["y"] - offsets[1]): f["frame"]
        for f in expected["frames"]
    }
    out = {}
    for (x, y, _z), mzs, intensities in reader.iter_spectra():
        out[by_coord[(x, y)]] = (mzs, intensities)
    return out


class TestSyntheticFixture:
    def test_scan_sum_reads_every_scan_of_the_right_frame(self, expected):
        with _open("scan_sum") as reader:
            spectra = _spectra_by_frame(reader, expected)

        assert sorted(spectra) == [f["frame"] for f in expected["frames"]]
        for frame in expected["frames"]:
            mzs, intensities = spectra[frame["frame"]]
            # AccumulationTime is 100 ms in the fixture, so the SDK's
            # intensity scale is exactly 1 and the TIC is the written sum.
            assert intensities.sum() == pytest.approx(frame["tic"])
            assert mzs.size == frame["unique_indices"]
            assert np.all(np.diff(mzs) > 0)

    def test_first_frame_is_present(self, expected):
        # The old reader asked the SDK for frame_id - 1, which does not exist
        # for the first frame, and skipped that pixel with a warning.
        with _open("scan_sum") as reader:
            spectra = _spectra_by_frame(reader, expected)
        assert 1 in spectra

    def test_planted_ion_lands_on_its_calibrated_mz(self, expected):
        frame = expected["frames"][0]
        with _open("scan_sum") as reader:
            spectra = _spectra_by_frame(reader, expected)
            planted_mz = reader.sdk._convert_indices_to_mz(
                reader.handle, frame["frame"], np.array([float(frame["planted_index"])])
            )[0]
        mzs, intensities = spectra[frame["frame"]]
        hit = int(np.argmin(np.abs(mzs - planted_mz)))
        assert mzs[hit] == pytest.approx(planted_mz, abs=1e-9)
        assert intensities[hit] == pytest.approx(frame["planted_intensity"])

    def test_vendor_centroid_yields_one_spectrum_per_pixel(self, expected):
        with _open("scan_sum") as reader:
            summed = _spectra_by_frame(reader, expected)
        with _open("vendor_centroid") as reader:
            centroided = _spectra_by_frame(reader, expected)

        assert sorted(centroided) == sorted(summed)
        for frame_id, (mzs, intensities) in centroided.items():
            assert mzs.size > 0
            assert np.all(np.diff(mzs) > 0)
            # The vendor picker merges bins and drops single counts: never
            # more ion current than the lossless sum, never none at all.
            assert 0 < intensities.sum() <= summed[frame_id][1].sum() * (1 + 1e-9)
            assert mzs.size <= summed[frame_id][0].size

    def test_mobility_axis_decreases_with_scan_number(self, expected):
        with _open("scan_sum") as reader:
            scans = np.arange(expected["n_scans"], dtype=np.float64)
            k0 = reader.sdk.scannum_to_oneoverk0(reader.handle, 1, scans)
        assert k0.shape == scans.shape
        assert np.all(np.diff(k0) < 0)
        assert 0.5 < k0.min() < k0.max() < 2.5

    def test_no_per_pixel_peak_counts_for_tdf(self):
        with _open("scan_sum") as reader:
            assert reader.get_peak_counts_per_pixel() is None

    def test_extractor_reports_the_mobility_dimension(self, expected):
        with _open("scan_sum") as reader:
            info = reader.get_comprehensive_metadata().format_specific["ion_mobility"]
        assert info["present"] is True
        assert info["separation_accession"] == "MS:1002815"
        assert info["num_scans_min"] == info["num_scans_max"] == expected["n_scans"]
        assert info["one_over_k0_range"] == [0.9, 1.99]

    def test_conversion_records_the_mobility_block_and_the_summation(self, tmp_path):
        spatialdata = pytest.importorskip("spatialdata")
        from thyra.convert import convert_msi

        _open("scan_sum").close()  # skip early where the library is missing
        out = tmp_path / "synthetic_tims.zarr"
        ok = convert_msi(
            str(FIXTURE),
            str(out),
            dataset_id="tims",
            pixel_size_um=20.0,
            reader_options={"tdf_spectrum": "scan_sum"},
        )
        assert ok

        from thyra.metadata.schema import read_msi_metadata_blocks

        sdata = spatialdata.read_zarr(out)
        table = next(iter(sdata.tables.values()))
        assert table.n_obs == 6
        block = next(iter(read_msi_metadata_blocks(out).values()))
        mobility = block["ms_analysis"]["ion_mobility"]
        assert mobility["present"] is True
        assert mobility["num_scans"] == 240
        assert mobility["separation_term"]["accession"] == "MS:1002815"
        conversion = block["processing"][0]
        assert conversion["name"] == "conversion"
        assert conversion["parameters"]["tdf_spectrum"] == "scan_sum"


REAL_DATASET = os.environ.get("THYRA_BRUKER_TDF_DATASET")


@pytest.mark.skipif(
    not REAL_DATASET,
    reason="Set THYRA_BRUKER_TDF_DATASET to a TIMS .d directory to run",
)
class TestRealAcquisition:
    def test_every_frame_is_read_and_the_lossless_tic_matches_the_database(self):
        path = Path(REAL_DATASET)  # type: ignore[arg-type]
        con = sqlite3.connect(
            f"file:{(path / 'analysis.tdf').as_posix()}?mode=ro&immutable=1", uri=True
        )
        frames = dict(
            con.execute(
                "SELECT m.Frame, f.SummedIntensities FROM MaldiFrameInfo m "
                "JOIN Frames f ON f.Id = m.Frame ORDER BY m.Frame LIMIT 5"
            ).fetchall()
        )
        n_frames = con.execute("SELECT COUNT(*) FROM MaldiFrameInfo").fetchone()[0]
        con.close()

        with _open("scan_sum", path) as reader:
            n_pixels = sum(1 for _ in reader.iter_spectra())
            per_frame = {
                frame_id: reader.sdk.read_spectrum(
                    reader.handle, frame_id, num_scans=reader._frame_num_scans(frame_id)
                )
                for frame_id in frames
            }
        assert n_pixels == n_frames
        for frame_id, summed in frames.items():
            # SummedIntensities is the unrounded scaled sum; the SDK rounds
            # each pair after scaling, so allow a few counts per pair.
            tic = per_frame[frame_id][1].sum()
            assert tic == pytest.approx(summed, rel=5e-3)

        with _open("vendor_centroid", path) as reader:
            for frame_id in frames:
                mzs, intensities = reader.sdk.read_spectrum(
                    reader.handle, frame_id, num_scans=reader._frame_num_scans(frame_id)
                )
                assert mzs.size > 0
                assert (
                    0.5 * frames[frame_id]
                    < intensities.sum()
                    <= frames[frame_id] * 1.01
                )
