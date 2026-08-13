"""Tests for the PHI SmartSoft-TOF .raw reader.

Builds synthetic .raw files matching the documented layout so the reader can
be exercised without a multi-hundred-megabyte vendor file.
"""

import struct

import numpy as np
import pytest

from thyra.core.registry import detect_format
from thyra.readers.phi import PhiReader, build_mass_axis, parse_phi_header, scan_blocks
from thyra.readers.phi.event_stream import decode_events

HEADER_SIZE = 2048

# slope 1.0 / offset 0.0 makes m/z == (t_us)**2, and a 1000 ns channel makes
# one channel exactly one microsecond, so masses land on 1, 4, 9, ... 100.
_HEADER_TEMPLATE = """SOFH\r
Platform: PC\r
SoftwareVersion: SmartSoft-TOF V3.2.0.35\r
Technique: -TofSIMS\r
SecIonPolarity: -\r
DataType: unsigned long\r
SpecBinSize: {bin_size_ns}\r
SpecBinIncr: 0.00100\r
StartFlightTime: {start_us}\r
StopFlightTime: {stop_us}\r
Mass/Time: {slope}\r
MassOffset: {offset}\r
ImagePixels: {pixels}\r
ScanWidthX: 512.00000\r
NoFrames: {frames}\r
Calibrated: no\r
Calibration: 2  (26.003016, C+N, 26.003099)  (41.998143, C+N+O, 41.998001)\r
HeaderSize: {header_size}\r
[Acq Base]\r
Start Mass (amu)={start_mz}\r
End Mass (amu)={stop_mz}\r
Frames={frames}\r
MSMS Active={msms}\r
Raster Pattern=Scatter\r
Raster Resolution={pixels}\r
[PHI LMIG]\r
Gun Particle=Bi3 +\r
Raster Size (um)={raster_um}\r
Raster Size Calibration=1.550, 1.550\r
[Mosaic Area]\r
Number Of Tiles X={tiles_x}\r
Number Of Tiles Y={tiles_y}\r
[Polarity]\r
Polarity=Negative (-) Ions\r
[Main]\r
Project Name=SmartSoft-TOF\r
EOFH\r
"""


def make_header(**overrides) -> bytes:
    """Render an ASCII header, zero-padded to HeaderSize."""
    fields = {
        "bin_size_ns": "1000.00000",
        "start_us": "0.00000",
        "stop_us": "10.00000",
        "slope": "1.0000000",
        "offset": "0.0000000",
        "pixels": 4,
        "frames": 2,
        "header_size": HEADER_SIZE,
        "start_mz": "1.0",
        "stop_mz": "100.0",
        "msms": "Inactive",
        "raster_um": "8.0",
        "tiles_x": 1,
        "tiles_y": 1,
    }
    fields.update(overrides)
    text = _HEADER_TEMPLATE.format(**fields).encode("latin-1")
    if len(text) > HEADER_SIZE:  # pragma: no cover - guards the fixture itself
        raise AssertionError("synthetic header exceeds HeaderSize")
    return text.ljust(HEADER_SIZE, b"\x00")


def event(x: int, y: int, tof_ps: int) -> int:
    """Encode one ion event as its uint64 record."""
    return (1 << 31) | (tof_ps & 0x7FFFFFFF) | (x << 32) | (y << 43)


def block(block_id: int, payload: bytes = b"", record_size: int = None) -> bytes:
    """Render one block header plus payload."""
    size = len(payload) if record_size is None else record_size
    return (
        struct.pack("<H", block_id) + b"\x00" * 10 + struct.pack("<I", size) + payload
    )


def events_block(words, block_id: int = 1) -> bytes:
    """Render a block of ion events."""
    return block(block_id, np.array(words, dtype="<u8").tobytes())


def write_raw(path, blocks: bytes, **header_overrides):
    """Write a complete synthetic .raw file."""
    path.write_bytes(make_header(**header_overrides) + blocks + block(0))
    return path


@pytest.fixture
def simple_raw(tmp_path):
    """Two pixels with a handful of events, split across two frames."""
    blocks = (
        events_block([event(0, 0, 3_000_000), event(0, 0, 3_000_000)])
        + events_block([event(1, 2, 5_000_000)])
        + block(2)
        + events_block([event(0, 0, 3_000_000), event(1, 2, 10_000_000)])
        + block(2)
    )
    return write_raw(tmp_path / "simple.raw", blocks)


class TestHeaderParsing:
    def test_parses_core_fields(self, simple_raw):
        header = parse_phi_header(simple_raw)
        assert header.header_size == HEADER_SIZE
        assert header.image_pixels == 4
        assert header.mass_slope == 1.0
        assert header.mass_offset == 0.0
        assert header.bin_size_ns == 1000.0
        assert header.polarity == "negative"
        assert header.n_frames == 2
        assert header.msms_active is False
        assert header.data_block_id == 1

    def test_parses_sections_and_calibrants(self, simple_raw):
        header = parse_phi_header(simple_raw)
        assert header.sections["PHI LMIG"]["Gun Particle"] == "Bi3 +"
        assert header.sections["Acq Base"]["Raster Pattern"] == "Scatter"
        assert [c.species for c in header.calibrants] == ["C+N", "C+N+O"]
        assert header.calibrants[0].theoretical_mz == pytest.approx(26.003099)

    def test_pixel_size_ignores_raster_calibration(self, simple_raw):
        """Raster Size Calibration trims the scan generator, it is not a scale.

        Applying it would inflate every coordinate by 55%.
        """
        header = parse_phi_header(simple_raw)
        assert header.raster_size_calibration == pytest.approx(1.55)
        # 8 um raster over 4 pixels
        assert header.pixel_size_um == pytest.approx(2.0)

    def test_msms_selects_block_id_eight(self, tmp_path):
        path = write_raw(tmp_path / "msms.raw", b"", msms="Active")
        assert parse_phi_header(path).data_block_id == 8

    def test_rejects_non_phi_file(self, tmp_path):
        path = tmp_path / "other.raw"
        path.write_bytes(b"NOTPHI" + b"\x00" * 100)
        with pytest.raises(ValueError, match="SOFH"):
            parse_phi_header(path)

    def test_rejects_unterminated_header(self, tmp_path):
        path = tmp_path / "trunc.raw"
        path.write_bytes(b"SOFH\r\nImagePixels: 4\r\n")
        with pytest.raises(ValueError, match="EOFH"):
            parse_phi_header(path)

    def test_rejects_empty_flight_time_range(self, tmp_path):
        path = write_raw(tmp_path / "bad.raw", b"", start_us="9.0", stop_us="1.0")
        with pytest.raises(ValueError, match="flight-time range"):
            parse_phi_header(path)


class TestBlockScan:
    def test_indexes_blocks_and_terminates_cleanly(self, simple_raw):
        header = parse_phi_header(simple_raw)
        index = scan_blocks(simple_raw, header)
        assert index.clean is True
        assert index.end_offset == simple_raw.stat().st_size
        assert len(index.data_spans) == 3
        assert index.n_records == 5
        assert index.n_frames == 2
        assert index.tile_grid == (1, 1)

    def test_unknown_blocks_are_skipped_by_record_size(self, tmp_path):
        blocks = (
            block(99, b"\xde\xad\xbe\xef" * 4)
            + events_block([event(1, 1, 4_000_000)])
            + block(2)
        )
        path = write_raw(tmp_path / "unknown.raw", blocks)
        index = scan_blocks(path, parse_phi_header(path))
        assert index.clean is True
        assert index.n_records == 1
        assert index.block_counts[99] == 1

    def test_raises_when_no_event_blocks(self, tmp_path):
        path = write_raw(tmp_path / "empty.raw", block(2))
        with pytest.raises(ValueError, match="No event blocks"):
            scan_blocks(path, parse_phi_header(path))

    def test_appended_calibration_is_recovered(self, tmp_path):
        appended = (
            b"AppendedBlockType: AsciiInfoBlock01\r\n"
            b"BlockAppendedDate: 07/27/2026 17:29:40\r\n"
            b"Mass/Time: 2.0\r\n"
            b"MassOffset: -0.5\r\n"
        )
        blocks = events_block([event(0, 0, 3_000_000)]) + block(14, appended)
        path = write_raw(tmp_path / "appended.raw", blocks)
        index = scan_blocks(path, parse_phi_header(path))
        assert index.appended_calibration() == (2.0, -0.5)

    def test_no_appended_block_returns_none(self, simple_raw):
        header = parse_phi_header(simple_raw)
        assert scan_blocks(simple_raw, header).appended_calibration() is None


class TestEventDecoding:
    def test_round_trips_coordinates_and_flight_time(self):
        words = np.array([event(511, 300, 21_113_342)], dtype=np.uint64)
        x, y, tof = decode_events(words)
        assert (x[0], y[0], tof[0]) == (511, 300, 21_113_342)

    def test_supports_eleven_bit_coordinates(self):
        words = np.array([event(2047, 2047, 1)], dtype=np.uint64)
        x, y, _ = decode_events(words)
        assert (x[0], y[0]) == (2047, 2047)

    def test_drops_records_failing_the_validity_test(self):
        words = np.array(
            [
                event(1, 1, 1000),
                0xADADADADADADADAD,  # uninitialised-heap fill
                event(2, 2, 2000) & ~(1 << 31),  # tag bit clear
            ],
            dtype=np.uint64,
        )
        x, _, _ = decode_events(words)
        assert x.size == 1


class TestMassAxis:
    def test_axis_is_clipped_to_declared_mass_range(self):
        axis = build_mass_axis(1.0, 0.0, 0.0, 10.0, 1000.0, 1.0, 100.0)
        assert len(axis) == 10
        assert axis.channel_offset == 1
        np.testing.assert_allclose(axis.mz, [1, 4, 9, 16, 25, 36, 49, 64, 81, 100])

    def test_channels_map_back_to_the_axis(self):
        axis = build_mass_axis(1.0, 0.0, 0.0, 10.0, 1000.0, 1.0, 100.0)
        idx, ok = axis.to_channel(np.array([3_000_000, 10_000_000], dtype=np.int64))
        assert ok.all()
        np.testing.assert_array_equal(idx, [2, 9])
        np.testing.assert_allclose(axis.mz[idx], [9.0, 100.0])

    def test_out_of_range_flight_times_are_masked(self):
        axis = build_mass_axis(1.0, 0.0, 0.0, 10.0, 1000.0, 1.0, 100.0)
        _, ok = axis.to_channel(np.array([0, 500_000_000], dtype=np.int64))
        assert not ok.any()

    def test_tof_axis_matches_the_mz_axis(self):
        axis = build_mass_axis(1.0, 0.0, 0.0, 10.0, 1000.0, 1.0, 100.0)
        # slope 1, offset 0 -> m/z is exactly t**2
        np.testing.assert_allclose(axis.tof_us, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        np.testing.assert_allclose(axis.tof_us**2, axis.mz)

    def test_calibration_round_trips(self):
        axis = build_mass_axis(0.382674, -1.598302, 0.0, 120.0, 0.128, 0.5, 1850.0)
        np.testing.assert_allclose(axis.mz_to_tof_us(axis.mz), axis.tof_us, rtol=1e-9)

    def test_recalibration_from_stored_times(self):
        """A corrected axis follows from the stored times, without the file.

        The time axis does not depend on the calibration, so two axes over
        the same channels share it exactly and differ only in m/z.
        """
        stale = build_mass_axis(0.3826529, -1.5976032, 0.0, 120.0, 0.128)
        good = build_mass_axis(0.382674, -1.598302, 0.0, 120.0, 0.128)

        # The axes begin at different channels -- the offset moves where the
        # root turns positive -- but both run to the same stop time on the
        # same grid, so they align from the tail.
        n = min(stale.tof_us.size, good.tof_us.size)
        np.testing.assert_allclose(stale.tof_us[-n:], good.tof_us[-n:], rtol=1e-12)

        # Applying the corrected coefficients to the stored times recovers
        # the corrected masses, so nothing is lost by storing the stale ones.
        recomputed = (0.382674 * stale.tof_us[-n:] - 1.598302) ** 2
        np.testing.assert_allclose(recomputed, good.mz[-n:], rtol=1e-9)

    def test_rejects_non_positive_bin_size(self):
        with pytest.raises(ValueError, match="SpecBinSize"):
            build_mass_axis(1.0, 0.0, 0.0, 10.0, 0.0)

    def test_rejects_empty_mass_window(self):
        with pytest.raises(ValueError, match="No time channel"):
            build_mass_axis(1.0, 0.0, 0.0, 10.0, 1000.0, 500.0, 900.0)


class TestPhiReader:
    def test_dimensions_and_axis(self, simple_raw):
        with PhiReader(simple_raw) as reader:
            assert reader.dimensions == (4, 4, 1)
            np.testing.assert_allclose(
                reader.get_common_mass_axis(),
                [1, 4, 9, 16, 25, 36, 49, 64, 81, 100],
            )

    def test_aggregates_events_per_pixel(self, simple_raw):
        with PhiReader(simple_raw) as reader:
            spectra = {
                (x, y): (mz, inten) for (x, y, _), mz, inten in reader.iter_spectra()
            }
        assert set(spectra) == {(0, 0), (1, 2)}
        # (0,0) saw three events all at 3 us -> m/z 9, summed
        np.testing.assert_allclose(spectra[(0, 0)][0], [9.0])
        np.testing.assert_allclose(spectra[(0, 0)][1], [3.0])
        # (1,2) saw one at 5 us (m/z 25) and one at 10 us (m/z 100)
        np.testing.assert_allclose(spectra[(1, 2)][0], [25.0, 100.0])
        np.testing.assert_allclose(spectra[(1, 2)][1], [1.0, 1.0])

    def test_spectra_are_sorted_by_mz(self, simple_raw):
        with PhiReader(simple_raw) as reader:
            for _, mzs, _ in reader.iter_spectra():
                assert np.all(np.diff(mzs) > 0)

    def test_peak_counts_match_iterated_spectra(self, simple_raw):
        with PhiReader(simple_raw) as reader:
            counts = reader.get_peak_counts_per_pixel()
            n_x = reader.dimensions[0]
            expected = np.zeros_like(counts)
            for (x, y, _), mzs, _ in reader.iter_spectra():
                expected[y * n_x + x] = mzs.size
        np.testing.assert_array_equal(counts, expected)
        # (0,0) holds one occupied channel, (1,2) holds two
        assert counts.sum() == 3

    def test_essential_metadata(self, simple_raw):
        with PhiReader(simple_raw) as reader:
            meta = reader.get_essential_metadata()
        assert meta.dimensions == (4, 4, 1)
        assert meta.n_spectra == 2
        assert meta.total_peaks == 3
        assert meta.pixel_size == pytest.approx((2.0, 2.0))
        assert meta.spectrum_type == "profile spectrum"
        assert meta.coordinate_bounds == (0.0, 3.0, 0.0, 3.0)

    def test_comprehensive_metadata_records_calibration_source(self, simple_raw):
        with PhiReader(simple_raw) as reader:
            meta = reader.get_comprehensive_metadata()
        cal = meta.raw_metadata["calibration"]
        assert cal["source"] == "header"
        assert cal["mass_slope_used"] == 1.0
        assert cal["header_calibrated_flag"] == "no"
        assert meta.instrument_info["vendor"] == "Physical Electronics (PHI)"
        assert meta.format_specific["block_chain_complete"] is True

    def test_iteration_is_repeatable(self, simple_raw):
        with PhiReader(simple_raw) as reader:
            first = [(c, m.tolist()) for c, m, _ in reader.iter_spectra()]
            reader.reset()
            second = [(c, m.tolist()) for c, m, _ in reader.iter_spectra()]
        assert first == second

    def test_pixel_size_override(self, simple_raw):
        with PhiReader(simple_raw, pixel_size_um=1.55) as reader:
            assert reader.pixel_size_um == pytest.approx(1.55)
            meta = reader.get_essential_metadata()
        assert meta.pixel_size == pytest.approx((1.55, 1.55))

    def test_intensity_threshold_filters_low_counts(self, simple_raw):
        with PhiReader(simple_raw, intensity_threshold=2.0) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]
        # only pixel (0,0), with three summed counts, clears the threshold
        assert coords == [(0, 0, 0)]

    def test_appended_calibration_is_preferred(self, tmp_path):
        appended = b"AppendedBlockType: AsciiInfoBlock01\r\nMass/Time: 2.0\r\nMassOffset: 0.0\r\n"
        blocks = events_block([event(0, 0, 1_000_000)]) + block(14, appended)
        path = write_raw(tmp_path / "recal.raw", blocks)
        with PhiReader(path) as reader:
            assert reader.calibration_source == "appended"
            assert reader.mass_axis.slope == 2.0
            # 1 us at slope 2 -> sqrt(m) = 2 -> m/z 4
            ((_, mzs, _),) = list(reader.iter_spectra())
            np.testing.assert_allclose(mzs, [4.0])

    def test_appended_calibration_can_be_disabled(self, tmp_path):
        appended = b"AppendedBlockType: AsciiInfoBlock01\r\nMass/Time: 2.0\r\nMassOffset: 0.0\r\n"
        blocks = events_block([event(0, 0, 1_000_000)]) + block(14, appended)
        path = write_raw(tmp_path / "recal.raw", blocks)
        with PhiReader(path, use_appended_calibration=False) as reader:
            assert reader.calibration_source == "header"
            assert reader.mass_axis.slope == 1.0

    def test_mosaic_tiles_offset_coordinates(self, tmp_path):
        blocks = (
            events_block([event(1, 1, 3_000_000)])
            + block(6, struct.pack("<II", 1, 0))
            + events_block([event(1, 1, 3_000_000)])
        )
        path = write_raw(tmp_path / "mosaic.raw", blocks, tiles_x=2, tiles_y=1)
        with PhiReader(path) as reader:
            assert reader.dimensions == (8, 4, 1)
            coords = sorted(c for c, _, _ in reader.iter_spectra())
        # second tile is offset by one raster width in x
        assert coords == [(1, 1, 0), (5, 1, 0)]

    def test_declared_tiles_do_not_inflate_the_grid(self, simple_raw):
        """Tile counts are declared even when mosaic mode is off."""
        path = write_raw(
            simple_raw.parent / "declared.raw",
            events_block([event(0, 0, 3_000_000)]),
            tiles_x=8,
            tiles_y=2,
        )
        with PhiReader(path) as reader:
            assert reader.header.declared_tiles == (8, 2)
            assert reader.dimensions == (4, 4, 1)

    def test_exposes_flight_time_as_a_mass_axis_annotation(self, simple_raw):
        with PhiReader(simple_raw) as reader:
            ann = reader.get_mass_axis_annotations()
            axis = reader.get_common_mass_axis()
        assert set(ann) == {"tof_us"}
        assert len(ann["tof_us"]) == len(axis)
        np.testing.assert_allclose(ann["tof_us"], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    def test_rejects_directory(self, tmp_path):
        directory = tmp_path / "waters.raw"
        directory.mkdir()
        with pytest.raises(ValueError, match="must be a file"):
            PhiReader(directory)


class TestFlightTimeSurvivesConversion:
    """The measured flight time must reach the store on both write routes.

    The streaming route writes ``var`` to Zarr by hand rather than through a
    dataframe, so the two paths have to be checked separately.
    """

    @pytest.fixture
    def raw_with_events(self, tmp_path):
        words = [
            event(x, y, tof * 1_000_000)
            for x in range(4)
            for y in range(4)
            for tof in (3, 5, 9)
        ]
        return write_raw(tmp_path / "conv.raw", events_block(words) + block(2))

    @pytest.mark.parametrize("streaming", [False, True])
    def test_tof_axis_is_written(self, raw_with_events, tmp_path, streaming):
        sd = pytest.importorskip("spatialdata")
        from thyra.convert import convert_msi

        out = tmp_path / f"out_{int(streaming)}.zarr"
        assert convert_msi(raw_with_events, out, dataset_id="d", streaming=streaming)

        table = next(iter(sd.read_zarr(out).tables.values()))
        assert "tof_us" in table.var.columns

        # The stored times reproduce the stored m/z under the same calibration.
        cal = table.uns["raw_metadata"]["calibration"]
        recomputed = (
            cal["mass_slope_used"] * table.var["tof_us"].to_numpy()
            + cal["mass_offset_used"]
        ) ** 2
        np.testing.assert_allclose(recomputed, table.var["mz"].to_numpy(), rtol=1e-9)


class TestRegistryDetection:
    def test_detects_phi_raw_file(self, simple_raw):
        assert detect_format(simple_raw) == "phi"

    def test_detects_waters_raw_directory(self, tmp_path):
        directory = tmp_path / "waters.raw"
        directory.mkdir()
        (directory / "_FUNC001.DAT").write_bytes(b"\x00")
        assert detect_format(directory) == "waters"

    def test_waters_directory_without_func_files_still_errors(self, tmp_path):
        directory = tmp_path / "empty.raw"
        directory.mkdir()
        with pytest.raises(ValueError, match="_FUNC"):
            detect_format(directory)

    def test_unrecognised_raw_file_is_reported(self, tmp_path):
        path = tmp_path / "mystery.raw"
        path.write_bytes(b"XXXX" + b"\x00" * 64)
        with pytest.raises(ValueError, match="Unrecognised .raw file"):
            detect_format(path)
