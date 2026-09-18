"""A Bruker acquisition that imaged nothing, from the database to the CLI.

A ``.d`` written by a run with no raster -- an electrospray acquisition
on a timsTOF-family instrument -- has ``GlobalMetadata``, ``Frames`` and
the MS/MS tables like any other, and neither ``MaldiFrameInfo`` nor
``MaldiFrameLaserInfo``.  Thyra asked for the laser table first, so the
whole metadata block was unreachable for such a file; the mass range was
never even read.

Everything here is built from the schema of a real timsOmni TSF (schema
3.4) and filled with that acquisition's shapes: 120 fragment frames of
one precursor, whose collision energy steps part way through the run.
The file itself is not in the repository; these tables are.

The two outcomes this pins are opposite on purpose. The metadata is
obtainable, and the conversion is refused -- there are no pixels to
place, and no pixel size would make any.
"""

import sqlite3
from pathlib import Path

import pytest

from thyra.errors import ConversionRefused
from thyra.metadata.document import read_metadata_document
from thyra.readers.bruker.timstof.timstof_reader import BrukerReader

# The columns a TSF 3.4 ``Frames`` table has, in order. Only Polarity,
# MsMsType and NumPeaks are read here; the rest are present because a
# fixture that drops columns tests a schema nobody ships.
_FRAMES_DDL = """
CREATE TABLE Frames (
    Id INTEGER PRIMARY KEY, Time REAL NOT NULL, Polarity CHAR(1) NOT NULL,
    ScanMode INTEGER NOT NULL, MsMsType INTEGER NOT NULL, TimsId INTEGER,
    MaxIntensity INTEGER, SummedIntensities INTEGER, NumPeaks INTEGER,
    MzCalibration INTEGER NOT NULL, T1 REAL, T2 REAL, PropertyGroup INTEGER)
"""

# One row per MS/MS frame, which is what a non-PASEF acquisition writes.
_FRAME_MSMS_DDL = """
CREATE TABLE FrameMsMsInfo (
    Frame INTEGER PRIMARY KEY, Parent INTEGER, TriggerMass REAL,
    IsolationWidth REAL, PrecursorCharge INTEGER, CollisionEnergy REAL)
"""

_GLOBAL_METADATA = {
    "AcquisitionDateTime": "2025-12-10T15:21:16.151+02:00",
    "AcquisitionSoftware": "timsTOF",
    "AcquisitionSoftwareVersion": "7.2.0",
    "InstrumentName": "timsOmni",
    "InstrumentSerialNumber": "0000000.00000",
    "InstrumentSourceType": "11",
    "InstrumentVendor": "Bruker",
    "MethodName": "fragment_one_precursor.m",
    "MzAcqRangeLower": "100.000000",
    "MzAcqRangeUpper": "8000.000000",
    "SchemaType": "TSF",
    "SchemaVersionMajor": "3",
    "SchemaVersionMinor": "4",
}

#: Frames 1..34 fragment the precursor at 50 eV and 35..120 at 60 eV --
#: one precursor, two energies, which is two isolation windows.
_N_FRAMES = 120
_ENERGY_STEPS_AT = 35
_PRECURSOR_MZ = 3888.0
_ISOLATION_WIDTH = 10.0


def _write_tsf(directory: Path, *, polarity: str = "+") -> Path:
    """Write a non-imaging TSF acquisition into ``directory``."""
    directory.mkdir(parents=True, exist_ok=True)
    database = directory / "analysis.tsf"
    with sqlite3.connect(database) as conn:
        conn.execute("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
        conn.execute(_FRAMES_DDL)
        conn.execute(_FRAME_MSMS_DDL)
        conn.executemany(
            "INSERT INTO GlobalMetadata (Key, Value) VALUES (?, ?)",
            sorted(_GLOBAL_METADATA.items()),
        )
        conn.executemany(
            "INSERT INTO Frames (Id, Time, Polarity, ScanMode, MsMsType, "
            "TimsId, MaxIntensity, SummedIntensities, NumPeaks, "
            "MzCalibration) VALUES (?, ?, ?, 13, 2, ?, 100, 1000, 500, 1)",
            [
                (frame, float(frame), polarity, frame * 64)
                for frame in range(1, _N_FRAMES + 1)
            ],
        )
        conn.executemany(
            "INSERT INTO FrameMsMsInfo (Frame, Parent, TriggerMass, "
            "IsolationWidth, PrecursorCharge, CollisionEnergy) "
            "VALUES (?, NULL, ?, ?, NULL, ?)",
            [
                (
                    frame,
                    _PRECURSOR_MZ,
                    _ISOLATION_WIDTH,
                    50.0 if frame < _ENERGY_STEPS_AT else 60.0,
                )
                for frame in range(1, _N_FRAMES + 1)
            ],
        )
    # A TSF acquisition always has its binary beside the database. Nothing
    # here reads it -- no spectra are decoded -- but a reader that probed
    # for it would be right to.
    (directory / "analysis.tsf_bin").write_bytes(b"\x00" * 64)
    return directory


@pytest.fixture
def no_raster(tmp_path) -> Path:
    return _write_tsf(tmp_path / "electrospray.d")


@pytest.fixture
def document(no_raster) -> dict:
    return read_metadata_document(no_raster)


class TestTheMetadataIsObtainable:
    """The block the released extractor could not reach at all."""

    def test_essential_metadata_carries_the_mass_range(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            essential = reader.get_essential_metadata()
        assert essential.mass_range == (100.0, 8000.0)

    def test_there_is_no_pixel_size(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            essential = reader.get_essential_metadata()
        assert essential.pixel_size is None
        assert essential.has_pixel_size is False

    def test_the_grid_is_empty_rather_than_invented(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            essential = reader.get_essential_metadata()
        assert essential.dimensions == (0, 0, 1)
        assert essential.coordinate_bounds == (0.0, 0.0, 0.0, 0.0)

    def test_the_frames_are_counted_off_the_frames_table(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            essential = reader.get_essential_metadata()
        assert essential.n_spectra == _N_FRAMES

    def test_the_file_is_not_called_maldi(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            comprehensive = reader.get_comprehensive_metadata()
        assert comprehensive.format_specific["is_maldi"] is False

    def test_the_global_metadata_survives_the_missing_laser_table(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            comprehensive = reader.get_comprehensive_metadata()
        assert (
            comprehensive.raw_metadata["global_metadata"]["InstrumentName"]
            == "timsOmni"
        )
        assert "frame_info" not in comprehensive.raw_metadata


class TestWhatTheDatabaseSaysAboutTheAcquisition:
    def test_polarity_comes_from_the_frames_table(self, document):
        assert document["ms_analysis"]["polarity"] == "positive"
        assert document["ms_analysis"]["polarity_term"] == {
            "accession": "MS:1000130",
            "name": "positive scan",
        }

    def test_an_alternating_acquisition_states_no_polarity(self, tmp_path):
        directory = _write_tsf(tmp_path / "alternating.d")
        with sqlite3.connect(directory / "analysis.tsf") as conn:
            conn.execute("UPDATE Frames SET Polarity = '-' WHERE Id > 60")
        assert "polarity" not in read_metadata_document(directory)["ms_analysis"]

    def test_the_source_type_is_recorded_but_not_interpreted(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            comprehensive = reader.get_comprehensive_metadata()
        assert comprehensive.format_specific["instrument_source_type"] == 11
        # Two codes have been seen on real files, 1 (MALDI imaging) and 11
        # (electrospray), which is not an enum. Nothing is claimed from it.
        assert (
            "ionisation_source" not in read_metadata_document(no_raster)["ms_analysis"]
        )

    def test_the_instrument_model_is_read(self, document):
        assert document["ms_analysis"]["instrument_model"] == "timsOmni"

    def test_the_acquisition_keeps_its_utc_offset(self, document):
        assert (
            document["acquisition"]["acquisition_datetime"]
            == "2025-12-10T15:21:16.151+02:00"
        )

    def test_the_method_is_named_by_its_file(self, document):
        assert document["acquisition"]["method_file"] == "fragment_one_precursor.m"

    def test_a_tsf_records_no_mobility(self, document):
        assert document["ms_analysis"]["ion_mobility"]["present"] is False


class TestFragmentationOnATsf:
    """``FrameMsMsInfo`` is a TSF table, and was never asked about."""

    def test_the_spectra_are_fragment_spectra(self, document):
        fragmentation = document["ms_analysis"]["fragmentation"]
        assert fragmentation["present"] is True
        assert fragmentation["ms_level"] == 2

    def test_one_precursor_at_two_energies_is_two_windows(self, document):
        windows = document["ms_analysis"]["fragmentation"]["windows"]
        assert [window["target"] for window in windows] == [
            _PRECURSOR_MZ,
            _PRECURSOR_MZ,
        ]
        assert sorted(window["collision_energy"] for window in windows) == [50.0, 60.0]

    def test_the_isolation_width_is_split_into_offsets(self, document):
        window = document["ms_analysis"]["fragmentation"]["windows"][0]
        assert window["lower_offset"] == _ISOLATION_WIDTH / 2
        assert window["upper_offset"] == _ISOLATION_WIDTH / 2

    def test_a_schedule_that_changes_mid_run_is_not_called_constant(self, document):
        # The energy steps at frame 35, so neither window is at every
        # frame. Recorded as it is rather than as the common case.
        assert (
            document["ms_analysis"]["fragmentation"]["constant_across_pixels"] is False
        )

    def test_the_schedule_names_the_format_that_reported_it(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            assert reader.get_fragmentation().source == "bruker_tsf"

    def test_a_survey_tsf_says_so(self, tmp_path):
        directory = _write_tsf(tmp_path / "survey.d")
        with sqlite3.connect(directory / "analysis.tsf") as conn:
            conn.execute("UPDATE Frames SET MsMsType = 0")
            conn.execute("DELETE FROM FrameMsMsInfo")
        fragmentation = read_metadata_document(directory)["ms_analysis"][
            "fragmentation"
        ]
        assert fragmentation["present"] is False
        assert fragmentation["ms_level"] == 1
        assert fragmentation["windows"] == []

    def test_precursors_cannot_be_separated_without_a_mobility_ramp(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            assert reader.has_precursor_spectra is False


class TestTheDocumentShape:
    def test_the_pixel_size_is_absent_rather_than_invented(self, document):
        assert "pixel_size_um" not in document["ms_analysis"]

    def test_nothing_has_been_done_to_the_data(self, document):
        assert "processing" not in document

    def test_the_windows_are_a_list_not_a_json_string(self, document):
        # ``to_uns_dict`` packs them for zarr; a document is not a store.
        assert isinstance(document["ms_analysis"]["fragmentation"]["windows"], list)

    def test_the_document_does_not_validate_against_this_schema_version(self, document):
        from thyra.metadata.schema import validate_document

        _, issues = validate_document(document)
        errors = [issue for issue in issues if issue.severity == "error"]
        assert [issue.location for issue in errors] == ["ms_analysis.pixel_size_um"]


class TestTheConversionIsRefused:
    def test_building_a_converting_reader_refuses(self, no_raster):
        with pytest.raises(ConversionRefused) as refusal:
            BrukerReader(no_raster)
        assert "MaldiFrameInfo" in str(refusal.value)

    def test_the_refusal_says_what_to_run_instead(self, no_raster):
        with pytest.raises(ConversionRefused) as refusal:
            BrukerReader(no_raster)
        assert "thyra metadata" in str(refusal.value)

    def test_a_metadata_only_reader_is_not_refused(self, no_raster):
        with BrukerReader(no_raster, metadata_only=True) as reader:
            assert reader.file_type == "tsf"


class TestAnUnreadableDatabaseIsSomeoneElsesError:
    """The check runs before the database is opened properly, and stays quiet.

    It reads the table list through a connection of its own, so it is the
    first thing to touch the file -- ahead of ``_initialize_database``,
    whose message for a locked file names the file and the application
    likely holding it. A refusal from here would put a worse message in
    front of a better one, so a database it cannot read is not refused.
    """

    def test_an_unreadable_database_falls_through(self, tmp_path, monkeypatch):
        import sqlite3 as sqlite

        from thyra.readers.bruker.timstof import timstof_reader

        directory = _write_tsf(tmp_path / "locked.d")

        def refuse_to_open(*args, **kwargs):
            raise sqlite.OperationalError("database is locked")

        monkeypatch.setattr(timstof_reader, "open_read_only", refuse_to_open)
        with pytest.raises(Exception) as failure:
            BrukerReader(directory)
        # Whatever it is, it is not this check's refusal.
        assert "MaldiFrameInfo" not in str(failure.value)
