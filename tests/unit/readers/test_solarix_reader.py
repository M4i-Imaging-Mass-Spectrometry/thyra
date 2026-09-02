"""Tests for the Bruker solariX peaks.sqlite reader.

Builds synthetic .d directories matching the verified solariX layout --
a peaks.sqlite with hand-packed blobs plus ImagingInfo.xml and a sibling
.mis -- so the reader can be exercised without committing any vendor data.
"""

import logging
import sqlite3

import numpy as np
import pytest

from thyra.core.registry import detect_format, get_reader_class
from thyra.preview import preview_msi
from thyra.readers.bruker import BrukerFolderStructure, BrukerFormat, SolarixReader
from thyra.resampling.types import AxisType, ResamplingMethod

# Properties as found in a real acquisition (ftmsControl 2.2, schema
# Imaging 1.2); tests override individual keys to probe the refusals.
DEFAULT_PROPERTIES = {
    "SchemaType": "Imaging",
    "SchemaVersionMajor": "1",
    "SchemaVersionMinor": "2",
    "AcquisitionSoftware": "ftmsControl",
    "AcquisitionSoftwareVendor": "Bruker",
    "AcquisitionSoftwareVersion": "2.2",
    "InstrumentVendor": "Bruker",
    "InstrumentFamily": "513",
    "InstrumentSourceType": "7",
    "InstrumentSerialNumber": "1272500.00145",
    "AcquisitionDateTime": "2019-09-20T14:10:43.435+02:00",
    "OperatorName": "Admin",
    "MzAcqRangeLower": "100.53",
    "MzAcqRangeUpper": "1000.0",
}

# (Id, Polarity, ScanMode, AcquisitionMode, MsLevel) -- a positive full scan.
DEFAULT_ACQUISITION_KEYS = [(0, 0, 0, 33, 1)]

# Absolute stage-raster coordinates, deliberately far from 0 (a real region
# started at X 1239) so the offset-to-0-based behaviour is actually exercised.
# The four pixels cover a 2x2 bounding box minus one corner: sparse raster.
DEFAULT_SPECTRA = [
    {
        "id": 1,
        "x": 1239,
        "y": 166,
        "mzs": [100.967, 250.5, 999.731],
        "intensities": [1000.0, 2500.5, 300.25],
    },
    {"id": 2, "x": 1240, "y": 166, "mzs": [150.25, 250.5], "intensities": [10.5, 20.0]},
    {"id": 3, "x": 1240, "y": 167, "mzs": [500.125], "intensities": [42.0]},
]


def make_peaks_sqlite(
    d_dir,
    spectra=None,
    properties=None,
    acquisition_keys=None,
    blob_overrides=None,
):
    """Write a synthetic peaks.sqlite with hand-packed peak blobs."""
    spectra = DEFAULT_SPECTRA if spectra is None else spectra
    props = dict(DEFAULT_PROPERTIES)
    props.update(properties or {})
    keys = DEFAULT_ACQUISITION_KEYS if acquisition_keys is None else acquisition_keys
    blob_overrides = blob_overrides or {}

    conn = sqlite3.connect(d_dir / "peaks.sqlite")
    conn.execute("CREATE TABLE Properties (Key TEXT, Value TEXT)")
    conn.executemany("INSERT INTO Properties VALUES (?, ?)", props.items())
    conn.execute(
        "CREATE TABLE AcquisitionKeys (Id INTEGER, Polarity INTEGER, "
        "ScanMode INTEGER, AcquisitionMode INTEGER, MsLevel INTEGER)"
    )
    conn.executemany("INSERT INTO AcquisitionKeys VALUES (?, ?, ?, ?, ?)", keys)
    conn.execute(
        "CREATE TABLE Spectra (Id INTEGER, Chip INTEGER, SpotName TEXT, "
        "RegionNumber INTEGER, XIndexPos INTEGER, YIndexPos INTEGER, "
        "AcquisitionKey INTEGER, ParentMass REAL, DateTime TEXT, "
        "CalibrationDateTime TEXT, MotorPositionX REAL, MotorPositionY REAL, "
        "NumSummations INTEGER, LaserPower REAL, LaserRepRate REAL, "
        "NumPeaks INTEGER, PeakMzValues BLOB, PeakIntensityValues BLOB, "
        "PeakFwhmValues BLOB, PeakSnrValues BLOB, PeakFlags BLOB, "
        "BH BLOB, BC BLOB)"
    )
    for spec in spectra:
        mzs = np.asarray(spec["mzs"], dtype="<f8")
        intensities = np.asarray(spec["intensities"], dtype="<f4")
        n = spec.get("num_peaks", mzs.size)
        mz_blob = blob_overrides.get(("mz", spec["id"]), mzs.tobytes())
        int_blob = blob_overrides.get(("intensity", spec["id"]), intensities.tobytes())
        conn.execute(
            "INSERT INTO Spectra (Id, SpotName, RegionNumber, XIndexPos, "
            "YIndexPos, AcquisitionKey, NumSummations, LaserPower, "
            "LaserRepRate, NumPeaks, PeakMzValues, PeakIntensityValues, "
            "PeakFwhmValues, PeakSnrValues) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                spec["id"],
                f"R00X{spec['x']}Y{spec['y']}",
                spec.get("region", 0),
                spec["x"],
                spec["y"],
                spec.get("acquisition_key", 0),
                spec.get("num_summations", 300),
                spec.get("laser_power", 14.0),
                spec.get("laser_rep_rate", 2000.0),
                n,
                mz_blob,
                int_blob,
                np.full(mzs.size, 2e-4, dtype="<f4").tobytes(),
                np.full(mzs.size, 10.0, dtype="<f4").tobytes(),
            ),
        )
    conn.commit()
    conn.close()


def make_imaging_info(d_dir, n_scans):
    """Write a minimal ImagingInfo.xml with one <scan> entry per spectrum."""
    scans = "".join(
        f"<scan><count>{i + 1}</count><spotName>R00X0Y0</spotName>"
        f"<minutes>0.1</minutes><tic>4.4E9</tic><maxpeak>9.8E8</maxpeak></scan>"
        for i in range(n_scans)
    )
    (d_dir / "ImagingInfo.xml").write_text(
        f"<ImagingInfo>{scans}</ImagingInfo>", encoding="utf-8"
    )


def make_solarix_d(
    parent,
    name="sample",
    spectra=None,
    properties=None,
    acquisition_keys=None,
    blob_overrides=None,
    with_mis=True,
    raster=(50, 50),
    n_info_scans=None,
):
    """Build a synthetic solariX imaging .d with its sibling .mis."""
    spectra = DEFAULT_SPECTRA if spectra is None else spectra
    d_dir = parent / f"{name}.d"
    d_dir.mkdir()
    make_peaks_sqlite(
        d_dir,
        spectra=spectra,
        properties=properties,
        acquisition_keys=acquisition_keys,
        blob_overrides=blob_overrides,
    )
    make_imaging_info(d_dir, len(spectra) if n_info_scans is None else n_info_scans)
    (d_dir / "ser").write_bytes(b"\x00" * 64)
    if with_mis:
        (parent / f"{name}.mis").write_text(
            "<ImagingSequence>"
            f"<Raster>{raster[0]},{raster[1]}</Raster>"
            "</ImagingSequence>",
            encoding="utf-8",
        )
    return d_dir


@pytest.fixture
def solarix_d(tmp_path):
    """A default synthetic solariX .d with a matching .mis."""
    return make_solarix_d(tmp_path)


class TestSolarixDetection:
    """Format detection for the solariX layout and its traps."""

    def test_detect_format_solarix(self, solarix_d):
        assert detect_format(solarix_d) == "solarix"

    def test_folder_structure_detects_solarix(self, solarix_d):
        assert BrukerFolderStructure.detect_format(solarix_d) is BrukerFormat.SOLARIX

    def test_parent_folder_containing_solarix_d(self, tmp_path, solarix_d):
        assert detect_format(tmp_path) == "solarix"

    def test_metadata_files_found(self, solarix_d):
        info = BrukerFolderStructure(solarix_d).analyze()

        assert info.format is BrukerFormat.SOLARIX
        assert info.metadata_files["peaks"] == solarix_d / "peaks.sqlite"
        assert info.metadata_files["imaging_info"] == solarix_d / "ImagingInfo.xml"
        assert info.metadata_files["ser"] == solarix_d / "ser"

    def test_prescan_d_still_falls_through(self, tmp_path):
        """flexImaging pre-scan dirs (fid + analysis.baf) are not imaging runs."""
        prescan = tmp_path / "sample_0_C17_000001.d"
        prescan.mkdir()
        (prescan / "fid").write_bytes(b"\x00" * 16)
        (prescan / "analysis.baf").write_bytes(b"\x00" * 16)
        (prescan / "calib.bin").write_bytes(b"\x00" * 16)

        with pytest.raises(ValueError, match="missing analysis files"):
            detect_format(prescan)

    def test_solarix_family_without_peaks_names_the_fallback(self, tmp_path):
        """ser + ImagingInfo.xml but no peaks.sqlite: still solariX, but the
        error must point at the imzML export route, not claim non-detection."""
        d_dir = tmp_path / "old_acquisition.d"
        d_dir.mkdir()
        (d_dir / "ser").write_bytes(b"\x00" * 16)
        (d_dir / "ImagingInfo.xml").write_text("<ImagingInfo/>", encoding="utf-8")

        with pytest.raises(ValueError, match="without peaks.sqlite") as excinfo:
            detect_format(d_dir)
        assert "imzML" in str(excinfo.value)

    def test_timstof_detection_unaffected(self, tmp_path):
        d_dir = tmp_path / "tims.d"
        d_dir.mkdir()
        (d_dir / "analysis.tdf").touch()

        assert detect_format(d_dir) == "bruker"

    def test_registry_round_trip(self):
        assert get_reader_class("solarix") is SolarixReader


class TestSolarixReader:
    """Spectrum iteration and blob decoding."""

    def test_iter_spectra_round_trip(self, solarix_d):
        with SolarixReader(solarix_d) as reader:
            spectra = list(reader.iter_spectra())

        assert len(spectra) == 3
        coords, mzs, intensities = spectra[0]
        assert coords == (0, 0, 0)
        np.testing.assert_array_equal(mzs, np.array([100.967, 250.5, 999.731]))
        np.testing.assert_array_equal(
            intensities, np.float32([1000.0, 2500.5, 300.25]).astype(np.float64)
        )
        assert mzs.dtype == np.float64
        assert intensities.dtype == np.float64

    def test_absolute_indices_offset_to_zero_based(self, solarix_d):
        """XIndexPos starts at 1239, not 0 -- forget the offset and every
        image becomes a sliver in a huge canvas."""
        with SolarixReader(solarix_d) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]

            assert coords == [(0, 0, 0), (1, 0, 0), (1, 1, 0)]
            assert reader.coordinate_offsets == (1239, 166, 0)
            assert reader.dimensions == (2, 2, 1)

    def test_sparse_raster_pixel_absent(self, solarix_d):
        """The (1239, 167) cell was never acquired: no yield, zero count."""
        with SolarixReader(solarix_d) as reader:
            coords = {c[:2] for c, _, _ in reader.iter_spectra()}
            counts = reader.get_peak_counts_per_pixel()

        assert (0, 1) not in coords
        # pixel_idx = y * n_x + x with n_x = 2
        np.testing.assert_array_equal(counts, np.array([3, 2, 0, 1]))

    def test_truncated_mz_blob_is_refused(self, tmp_path):
        """A blob shorter than NumPeaks * itemsize must not decode silently."""
        d_dir = make_solarix_d(
            tmp_path,
            blob_overrides={("mz", 2): np.asarray([150.25], dtype="<f8").tobytes()},
        )

        with SolarixReader(d_dir) as reader:
            with pytest.raises(ValueError, match="Corrupt PeakMzValues.*Id=2"):
                list(reader.iter_spectra())

    def test_truncated_intensity_blob_is_refused(self, tmp_path):
        d_dir = make_solarix_d(
            tmp_path,
            blob_overrides={("intensity", 3): b"\x00"},
        )

        with SolarixReader(d_dir) as reader:
            with pytest.raises(ValueError, match="Corrupt PeakIntensityValues"):
                list(reader.iter_spectra())

    def test_empty_spectrum_is_skipped(self, tmp_path):
        spectra = DEFAULT_SPECTRA + [
            {"id": 4, "x": 1241, "y": 166, "mzs": [], "intensities": []}
        ]
        d_dir = make_solarix_d(tmp_path, spectra=spectra)

        with SolarixReader(d_dir) as reader:
            assert len(list(reader.iter_spectra())) == 3

    def test_unsorted_mz_values_are_sorted(self, tmp_path):
        spectra = [
            {
                "id": 1,
                "x": 10,
                "y": 20,
                "mzs": [500.0, 100.0, 300.0],
                "intensities": [5.0, 1.0, 3.0],
            }
        ]
        d_dir = make_solarix_d(tmp_path, spectra=spectra)

        with SolarixReader(d_dir) as reader:
            _, mzs, intensities = next(reader.iter_spectra())

        np.testing.assert_array_equal(mzs, [100.0, 300.0, 500.0])
        np.testing.assert_array_equal(intensities, [1.0, 3.0, 5.0])

    def test_intensity_threshold_filters(self, solarix_d):
        with SolarixReader(solarix_d, intensity_threshold=100.0) as reader:
            spectra = {c[:2]: (m, i) for c, m, i in reader.iter_spectra()}

        # Pixel (1, 0) had intensities 10.5 and 20.0 -- both filtered out.
        assert (1, 0) not in spectra
        mzs, intensities = spectra[(0, 0)]
        np.testing.assert_array_equal(mzs, [100.967, 250.5, 999.731])

    def test_common_mass_axis_is_sorted_union(self, solarix_d):
        with SolarixReader(solarix_d) as reader:
            axis = reader.get_common_mass_axis()

        np.testing.assert_array_equal(
            axis,
            np.unique(np.float64([100.967, 250.5, 999.731, 150.25, 250.5, 500.125])),
        )

    def test_has_shared_mass_axis_is_false(self, solarix_d):
        with SolarixReader(solarix_d) as reader:
            assert reader.has_shared_mass_axis is False

    def test_region_map_and_info(self, tmp_path):
        spectra = [
            {
                "id": 1,
                "x": 5,
                "y": 5,
                "mzs": [100.0],
                "intensities": [1.0],
                "region": 0,
            },
            {
                "id": 2,
                "x": 6,
                "y": 5,
                "mzs": [100.0],
                "intensities": [1.0],
                "region": 1,
            },
        ]
        d_dir = make_solarix_d(tmp_path, spectra=spectra)

        with SolarixReader(d_dir) as reader:
            assert reader.get_region_map() == {(0, 0): 0, (1, 0): 1}
            assert reader.get_region_info() == [
                {"region_number": 0, "n_spectra": 1},
                {"region_number": 1, "n_spectra": 1},
            ]

    def test_closed_reader_refuses_iteration(self, solarix_d):
        reader = SolarixReader(solarix_d)
        reader.close()

        with pytest.raises(RuntimeError, match="closed"):
            list(reader.iter_spectra())

    def test_scan_count_mismatch_warns(self, tmp_path, caplog):
        d_dir = make_solarix_d(tmp_path, n_info_scans=99)

        with caplog.at_level(logging.WARNING):
            SolarixReader(d_dir).close()

        assert any("truncated or partially copied" in m for m in caplog.messages)


class TestSolarixSchemaRefusal:
    """Unverified schemas must fail loudly with the found values."""

    def test_wrong_schema_type(self, tmp_path):
        d_dir = make_solarix_d(tmp_path, properties={"SchemaType": "LC-MS"})

        with pytest.raises(ValueError, match="SchemaType='LC-MS'"):
            SolarixReader(d_dir)

    def test_wrong_major_version(self, tmp_path):
        d_dir = make_solarix_d(tmp_path, properties={"SchemaVersionMajor": "2"})

        with pytest.raises(ValueError, match="SchemaVersionMajor='2'"):
            SolarixReader(d_dir)

    def test_no_peaks_sqlite_names_the_fallback(self, tmp_path):
        d_dir = tmp_path / "no_peaks.d"
        d_dir.mkdir()
        (d_dir / "ser").write_bytes(b"\x00" * 16)
        (d_dir / "ImagingInfo.xml").write_text("<ImagingInfo/>", encoding="utf-8")

        with pytest.raises(ValueError, match="imzML"):
            SolarixReader(d_dir)

    def test_missing_properties_table(self, tmp_path):
        d_dir = tmp_path / "empty.d"
        d_dir.mkdir()
        sqlite3.connect(d_dir / "peaks.sqlite").close()
        (d_dir / "ImagingInfo.xml").write_text("<ImagingInfo/>", encoding="utf-8")

        with pytest.raises(ValueError, match="Properties"):
            SolarixReader(d_dir)


class TestSolarixPixelSize:
    """The .mis sibling is the only pixel-size source; never guess."""

    def test_raster_from_matching_mis(self, solarix_d):
        with SolarixReader(solarix_d) as reader:
            assert reader.pixel_size_um == (50.0, 50.0)
            assert ".mis" in reader.pixel_size_source

    def test_missing_mis_reports_unknown(self, tmp_path):
        d_dir = make_solarix_d(tmp_path, with_mis=False)

        with SolarixReader(d_dir) as reader:
            assert reader.pixel_size_um is None
            assert "unknown" in reader.pixel_size_source
            assert reader.get_essential_metadata().pixel_size is None

    def test_single_non_matching_mis_is_used(self, tmp_path):
        d_dir = make_solarix_d(tmp_path, with_mis=False)
        (tmp_path / "renamed.mis").write_text(
            "<ImagingSequence><Raster>25,25</Raster></ImagingSequence>",
            encoding="utf-8",
        )

        with SolarixReader(d_dir) as reader:
            assert reader.pixel_size_um == (25.0, 25.0)

    def test_multiple_non_matching_mis_refused(self, tmp_path):
        """Several other regions' .mis files: guessing one would silently
        assign a wrong raster, so the pixel size stays unknown."""
        d_dir = make_solarix_d(tmp_path, with_mis=False)
        for name in ("region_2.mis", "region_3.mis"):
            (tmp_path / name).write_text(
                "<ImagingSequence><Raster>75,75</Raster></ImagingSequence>",
                encoding="utf-8",
            )

        with SolarixReader(d_dir) as reader:
            assert reader.pixel_size_um is None

    def test_user_override_wins(self, solarix_d):
        with SolarixReader(solarix_d, pixel_size_um=30.0) as reader:
            assert reader.pixel_size_um == (30.0, 30.0)
            assert reader.pixel_size_source == "user override"


class TestSolarixMetadata:
    """Essential and comprehensive metadata content."""

    def test_essential_metadata(self, solarix_d):
        with SolarixReader(solarix_d) as reader:
            essential = reader.get_essential_metadata()

        assert essential.dimensions == (2, 2, 1)
        assert essential.coordinate_bounds == (0.0, 1.0, 0.0, 1.0)
        assert essential.mass_range == (100.53, 1000.0)
        assert essential.pixel_size == (50.0, 50.0)
        assert essential.n_spectra == 3
        assert essential.total_peaks == 6
        assert essential.coordinate_offsets == (1239, 166, 0)
        assert essential.spectrum_type == "centroid spectrum"

    def test_metadata_only_matches_full_mode(self, solarix_d):
        with SolarixReader(solarix_d) as full:
            full_essential = full.get_essential_metadata()
        with SolarixReader(solarix_d, metadata_only=True) as cheap:
            cheap_essential = cheap.get_essential_metadata()

        assert cheap_essential == full_essential

    def test_instrument_info_exact_strings(self, solarix_d):
        """FTICRDetector string-matches "FT-ICR"; these are load-bearing."""
        with SolarixReader(solarix_d) as reader:
            info = reader.get_comprehensive_metadata().instrument_info

        assert info["instrument_type"] == "FT-ICR"
        assert info["manufacturer"] == "Bruker"
        assert info["instrument_model"] == "solariX"
        assert info["instrument_name"] == "solariX"
        assert info["instrument_family"] == 513
        assert info["serial_number"] == "1272500.00145"

    def test_unknown_instrument_family_keeps_raw_number_only(self, tmp_path):
        d_dir = make_solarix_d(tmp_path, properties={"InstrumentFamily": "999"})

        with SolarixReader(d_dir) as reader:
            info = reader.get_comprehensive_metadata().instrument_info

        assert info["instrument_type"] == "FT-ICR"
        assert info["instrument_family"] == 999
        assert "instrument_model" not in info

    def test_polarity_mapping(self, tmp_path):
        positive = make_solarix_d(tmp_path, name="pos")
        negative = make_solarix_d(
            tmp_path, name="neg", acquisition_keys=[(0, 1, 0, 33, 1)]
        )

        with SolarixReader(positive) as reader:
            params = reader.get_comprehensive_metadata().acquisition_params
            assert params["polarity"] == "positive"
            assert params["acquisition_key_0_polarity_raw"] == 0
        with SolarixReader(negative) as reader:
            params = reader.get_comprehensive_metadata().acquisition_params
            assert params["polarity"] == "negative"

    def test_mixed_polarity_keys_stay_unnamed(self, tmp_path):
        d_dir = make_solarix_d(
            tmp_path, acquisition_keys=[(0, 0, 0, 33, 1), (1, 1, 0, 33, 1)]
        )

        with SolarixReader(d_dir) as reader:
            params = reader.get_comprehensive_metadata().acquisition_params

        assert params["polarity"] is None
        assert params["acquisition_key_0_polarity_raw"] == 0
        assert params["acquisition_key_1_polarity_raw"] == 1

    def test_acquisition_params_carry_laser_and_raw_enums(self, solarix_d):
        with SolarixReader(solarix_d) as reader:
            params = reader.get_comprehensive_metadata().acquisition_params

        assert params["laser_power"] == 14.0
        assert params["laser_rep_rate"] == 2000.0
        assert params["num_summations"] == 300
        assert params["acquisition_key_0_scan_mode_raw"] == 0
        assert params["acquisition_key_0_acquisition_mode_raw"] == 33
        assert params["acquisition_key_0_ms_level_raw"] == 1

    def test_format_specific_block(self, solarix_d):
        with SolarixReader(solarix_d) as reader:
            fmt = reader.get_comprehensive_metadata().format_specific

        assert fmt["format"] == "Bruker solariX peaks.sqlite"
        assert fmt["schema_version"] == "1.2"
        assert fmt["acquisition_software"] == "ftmsControl"
        assert fmt["raster_x_range"] == [1239, 1240]
        assert fmt["raster_y_range"] == [166, 167]
        assert fmt["regions"] == [{"region_number": 0, "n_spectra": 3}]


class TestSolarixDecisionTree:
    """A native solariX .d must land on nearest_neighbor + fticr, no flags.

    Same surface the Ousia wizard reads, mirroring the solarix_fticr imzML
    fixture test -- the native route and the export route must agree.
    """

    def test_preview_reports_fticr_and_nearest_neighbor(self, solarix_d):
        preview = preview_msi(solarix_d)

        assert preview.readable, preview.error
        assert preview.instrument_type is AxisType.FTICR
        assert preview.resampling_method is ResamplingMethod.NEAREST_NEIGHBOR
        assert preview.n_pixels == 3
        assert preview.grid_dims == (2, 2)
        assert preview.pixel_size_um == 50.0
        assert preview.mz_range == (100.53, 1000.0)
