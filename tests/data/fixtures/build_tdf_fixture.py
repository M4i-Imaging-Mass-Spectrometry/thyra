"""Build ``synthetic_tims.d``: a hand-written Bruker TDF acquisition.

Everything in it is invented. The point of the fixture is that Bruker's own
``timsdata`` library opens it, walks every TIMS scan of every frame, converts
indices to m/z and scan numbers to 1/K0, and runs its frame-level centroid
extraction on it -- so the TDF reader can be tested end to end through the
real SDK without a byte of anyone's acquisition in the repository.

Layout, worked out against the reference readers (AlphaTims, timsrust) and
verified against the SDK on real files:

``analysis.tdf``
    SQLite with the TDF 3.7 schema. ``Frames`` carries one row per pixel
    with ``TimsId`` = byte offset of the frame block in ``analysis.tdf_bin``,
    ``NumScans`` = ramp length, ``NumPeaks`` = number of (index, scan) pairs.
    ``MaldiFrameInfo`` maps frames to raster positions. ``MzCalibration``
    (model 1) and ``TimsCalibration`` (model 2) hold representative
    coefficients so the SDK's conversions have something to evaluate; they
    describe an instrument model, not a measurement.

``analysis.tdf_bin``
    A 64-byte zero header, then one block per frame: ``uint32 block size``
    (including these 8 bytes), ``uint32 scan count``, zstd payload. The
    payload decompresses to a byte-shuffled ``uint32`` array (all first
    bytes, then all second bytes, ...). Word 0 is the scan count; words
    ``1..scan_count-1`` are ``2 * n_peaks`` for scans ``0..scan_count-2``
    (the last scan's length is implied); then, scan by scan, ``(delta, intensity)``
    pairs where the first delta of a scan is ``index + 1`` and later deltas
    are index differences.

The SDK returns intensities scaled to a 100 ms accumulation
(``100 / Frames.AccumulationTime``; measured 2.00 on a 50 ms run and 9.92 on
a 10.08 ms run, with ``SummedIntensities`` equal to the unrounded scaled
sum). ``AccumulationTime`` is 100.0 here so stored and returned intensities
agree exactly and the expected values below need no factor.

Regenerate with::

    PYTHONPATH=. python tests/data/fixtures/build_tdf_fixture.py

Content-identical on every run (seeded); the SQLite bytes can differ between
SQLite versions, the ``.tdf_bin`` bytes between zstd versions. The committed
files are the fixture; this script is its provenance.
"""

from __future__ import annotations

import json
import shutil
import sqlite3
import struct
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

FIXTURE_DIR = Path(__file__).resolve().parent
OUT_DIR = FIXTURE_DIR / "synthetic_tims.d"
EXPECTED_JSON = FIXTURE_DIR / "synthetic_tims_expected.json"

N_SCANS = 240
GRID: List[Tuple[int, int]] = [(x, y) for y in (10, 11) for x in (100, 101, 102)]
# Planted ions: (digitizer index, apex scan). Each is spread over seven scans
# with a Gaussian mobility profile and a half-height neighbour bin, the shape
# of a real TIMS peak in miniature.
IONS: List[Tuple[int, int]] = [
    (60000, N_SCANS // 3),
    (150000, N_SCANS // 2),
    (240000, (2 * N_SCANS) // 3),
]
N_NOISE_EVENTS = 40
ACCUMULATION_TIME_MS = 100.0  # makes the SDK's intensity scale exactly 1

# TDF 3.7 schema, as written by timsTOF acquisition software.
DDL_SQL = """
CREATE TABLE CalibrationInfo (
    KeyPolarity CHAR(1) CHECK (KeyPolarity IN ('+', '-')), KeyName TEXT, Value TEXT,
    PRIMARY KEY (KeyPolarity, KeyName));
CREATE TABLE CollisionEnergySweepingInfo (
    Frame INTEGER NOT NULL, CollisionId INTEGER NOT NULL, CollisionEnergy REAL NOT NULL,
    CollisionEnergyPercent REAL NOT NULL, PRIMARY KEY(Frame, CollisionId),
    FOREIGN KEY (Frame) REFERENCES Frames (Id)) WITHOUT ROWID;
CREATE TABLE DiaFrameMsMsInfo (
    Frame INTEGER PRIMARY KEY, WindowGroup INTEGER NOT NULL,
    FOREIGN KEY (Frame) REFERENCES Frames (Id),
    FOREIGN KEY (WindowGroup) REFERENCES DiaFrameMsMsWindowGroups (Id));
CREATE TABLE DiaFrameMsMsWindowGroups (Id INTEGER PRIMARY KEY);
CREATE TABLE DiaFrameMsMsWindows (
    WindowGroup INTEGER NOT NULL, ScanNumBegin INTEGER NOT NULL, ScanNumEnd INTEGER NOT NULL,
    IsolationMz REAL NOT NULL, IsolationWidth REAL NOT NULL, CollisionEnergy REAL NOT NULL,
    PRIMARY KEY(WindowGroup, ScanNumBegin),
    FOREIGN KEY (WindowGroup) REFERENCES DiaFrameMsMsWindowGroups (Id)) WITHOUT ROWID;
CREATE TABLE ErrorLog (Frame INTEGER NOT NULL, Scan INTEGER, Message TEXT NOT NULL);
CREATE TABLE FrameMsMsInfo (
    Frame INTEGER PRIMARY KEY, Parent INTEGER, TriggerMass REAL NOT NULL,
    IsolationWidth REAL NOT NULL, PrecursorCharge INTEGER, CollisionEnergy REAL NOT NULL,
    FOREIGN KEY (Frame) REFERENCES Frames (Id));
CREATE TABLE FrameProperties (
    Frame INTEGER NOT NULL, Property INTEGER NOT NULL, Value NOT NULL,
    PRIMARY KEY (Frame, Property), FOREIGN KEY (Frame) REFERENCES Frames (Id)
    FOREIGN KEY (Property) REFERENCES PropertyDefinitions (Id)) WITHOUT ROWID;
CREATE TABLE Frames (
    Id INTEGER PRIMARY KEY, Time REAL NOT NULL,
    Polarity CHAR(1) CHECK (Polarity IN ('+', '-')) NOT NULL, ScanMode INTEGER NOT NULL,
    MsMsType INTEGER NOT NULL, TimsId INTEGER, MaxIntensity INTEGER NOT NULL,
    SummedIntensities INTEGER NOT NULL, NumScans INTEGER NOT NULL, NumPeaks INTEGER NOT NULL,
    MzCalibration INTEGER NOT NULL, T1 REAL NOT NULL, T2 REAL NOT NULL,
    TimsCalibration INTEGER NOT NULL, PropertyGroup INTEGER, AccumulationTime REAL NOT NULL,
    RampTime REAL NOT NULL, Pressure REAL,
    FOREIGN KEY (MzCalibration) REFERENCES MzCalibration (Id),
    FOREIGN KEY (TimsCalibration) REFERENCES TimsCalibration (Id),
    FOREIGN KEY (PropertyGroup) REFERENCES PropertyGroups (Id));
CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);
CREATE TABLE GroupProperties (
    PropertyGroup INTEGER NOT NULL, Property INTEGER NOT NULL, Value NOT NULL,
    PRIMARY KEY (PropertyGroup, Property),
    FOREIGN KEY (PropertyGroup) REFERENCES PropertyGroups (Id),
    FOREIGN KEY (Property) REFERENCES PropertyDefinitions (Id)) WITHOUT ROWID;
CREATE TABLE Maldi2FrameLaserInfo (
    Id INTEGER PRIMARY KEY, TriggerDelayMicroseconds REAL NOT NULL,
    PreBurstLengthMilliseconds INTEGER NOT NULL, PowerCompensationLevel INTEGER NOT NULL,
    FOREIGN KEY(Id) REFERENCES MaldiFrameLaserInfo(Id));
CREATE TABLE MaldiFrameInfo (
    Frame INTEGER PRIMARY KEY NOT NULL, Chip INTEGER NOT NULL, SpotName TEXT,
    RegionNumber INTEGER, XIndexPos INTEGER, YIndexPos INTEGER, LaserPower REAL,
    NumLaserShots INTEGER NOT NULL, LaserRepRate REAL, MotorPositionX REAL,
    MotorPositionY REAL, MotorPositionZ REAL, LaserInfo INTEGER NOT NULL,
    FOREIGN KEY(Frame) REFERENCES Frames(Id),
    FOREIGN KEY(LaserInfo) REFERENCES MaldiFrameLaserInfo(Id));
CREATE TABLE MaldiFrameLaserInfo (
    Id INTEGER PRIMARY KEY, LaserApplicationName TEXT, LaserParameterName TEXT,
    LaserBoost REAL NOT NULL, LaserFocus REAL NOT NULL, BeamScan INTEGER NOT NULL,
    BeamScanSizeX REAL, BeamScanSizeY REAL, WalkOnSpotMode INTEGER NOT NULL,
    WalkOnSpotShots INTEGER, SpotSize REAL);
CREATE TABLE MzCalibration (
    Id INTEGER PRIMARY KEY, ModelType INTEGER NOT NULL, DigitizerTimebase REAL NOT NULL,
    DigitizerDelay REAL NOT NULL, T1 REAL NOT NULL, T2 REAL NOT NULL, dC1 REAL NOT NULL,
    dC2 REAL NOT NULL, C0, C1, C2, C3, C4);
CREATE TABLE PrmFrameMeasurementMode (
    Frame INTEGER PRIMARY KEY, MeasurementModeId TEXT,
    FOREIGN KEY (Frame) REFERENCES Frames(Id));
CREATE TABLE PrmFrameMsMsInfo (
    Frame INTEGER NOT NULL, ScanNumBegin INTEGER NOT NULL, ScanNumEnd INTEGER NOT NULL,
    IsolationMz REAL NOT NULL, IsolationWidth REAL NOT NULL, CollisionEnergy REAL NOT NULL,
    Target INTEGER NOT NULL, PRIMARY KEY (Frame, ScanNumBegin),
    FOREIGN KEY (Frame) REFERENCES Frames(Id),
    FOREIGN KEY (Target) REFERENCES PrmTargets(Id)) WITHOUT ROWID;
CREATE TABLE PrmTargets (
    Id INTEGER PRIMARY KEY,
    ExternalId TEXT CHECK(ExternalId IS NULL OR LENGTH(ExternalId) > 0) UNIQUE,
    Time REAL NOT NULL, OneOverK0 REAL NOT NULL, MonoisotopicMz REAL NOT NULL,
    Charge INTEGER NOT NULL, Description TEXT NOT NULL);
CREATE TABLE PropertyDefinitions (
    Id INTEGER PRIMARY KEY, PermanentName TEXT NOT NULL, Type INTEGER NOT NULL,
    DisplayGroupName TEXT NOT NULL, DisplayName TEXT NOT NULL, DisplayValueText TEXT NOT NULL,
    DisplayFormat TEXT NOT NULL, DisplayDimension TEXT NOT NULL, Description TEXT NOT NULL);
CREATE TABLE PropertyGroups (Id INTEGER PRIMARY KEY) WITHOUT ROWID;
CREATE TABLE Segments (
    Id INTEGER PRIMARY KEY, FirstFrame INTEGER NOT NULL, LastFrame INTEGER NOT NULL,
    IsCalibrationSegment BOOLEAN NOT NULL,
    FOREIGN KEY (FirstFrame) REFERENCES Frames (Id),
    FOREIGN KEY (LastFrame) REFERENCES Frames (Id));
CREATE TABLE TimsCalibration (
    Id INTEGER PRIMARY KEY, ModelType INTEGER NOT NULL, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9);
"""
DDL = [statement.strip() for statement in DDL_SQL.split(";") if statement.strip()]

# Representative TOF (model 1) and TIMS (model 2) calibration coefficients.
# They make index 10,000..300,000 span roughly m/z 60..940 and scan 0..N_SCANS
# span roughly 1/K0 2.0..0.9, decreasing with scan number as on a real ramp.
MZ_CALIBRATION = (
    1,
    1,
    0.2,
    18319.0,
    25.989297270499147,
    24.13372951260117,
    22.0,
    0.0,
    318.2681471005531,
    154305.24250607283,
    -0.0002600306374565382,
    0.0,
    0.0,
)
TIMS_CALIBRATION = (
    1,
    2,
    1,
    N_SCANS - 1,
    263.286764268143,
    91.66124001671024,
    42.857142857142854,
    1,
    0.06474908082027102,
    129.22284861490712,
    9.56073642101304,
    6917.454037200332,
)

GLOBAL_METADATA = {
    "SchemaType": "TDF",
    "SchemaVersionMajor": "3",
    "SchemaVersionMinor": "7",
    "TimsCompressionType": "2",
    "AcquisitionSoftware": "synthetic",
    "AcquisitionSoftwareVendor": "Bruker",
    "AcquisitionSoftwareVersion": "0",
    "AcquisitionDateTime": "2026-01-01T00:00:00.000+00:00",
    "InstrumentName": "timsTOF fleX",
    "InstrumentFamily": "9",
    "InstrumentRevision": "2",
    "InstrumentSerialNumber": "0",
    "InstrumentSourceType": "1",
    "InstrumentVendor": "Bruker",
    "MaldiApplicationType": "Imaging",
    "Geometry": "Imaging_Run",
    "MzAcqRangeLower": "50.000000",
    "MzAcqRangeUpper": "1000.000000",
    "OneOverK0AcqRangeLower": "0.900000",
    "OneOverK0AcqRangeUpper": "1.990000",
    "DigitizerNumSamples": "312509",
    "DigitizerType": "synthetic",
    "MaxNumPeaksPerScan": "5",
    "PeakListIndexScaleFactor": "1",
    "PeakWidthEstimateType": "1",
    "PeakWidthEstimateValue": "0.000025",
    "DenoisingEnabled": "0",
    "ClosedProperly": "1",
    "SampleName": "synthetic_tims",
    "MethodName": "synthetic.m",
    "OperatorName": "thyra",
    "Description": "hand-written test fixture; no acquisition data",
    "AnalysisId": "00000000-0000-0000-0000-000000000000",
    "ImagingAreaMinXIndexPos": "100",
    "ImagingAreaMaxXIndexPos": "102",
    "ImagingAreaMinYIndexPos": "10",
    "ImagingAreaMaxYIndexPos": "11",
}

Scan = Tuple[np.ndarray, np.ndarray]  # (sorted uint32 indices, uint32 intensities)


def make_frame(k: int, rng: np.random.Generator) -> List[Scan]:
    """Per-scan point lists for pixel ``k``: three planted ions plus noise."""
    per_scan: Dict[int, Dict[int, int]] = {}
    for tof, apex in IONS:
        tof_k = tof + 7 * k  # a slightly different index per pixel
        for ds in range(-3, 4):
            s = apex + ds
            if 0 <= s < N_SCANS - 1:  # the last scan stays empty on purpose
                inten = int(200 * np.exp(-0.5 * (ds / 1.5) ** 2)) + 5
                d = per_scan.setdefault(s, {})
                d[tof_k] = d.get(tof_k, 0) + inten
                d[tof_k + 1] = d.get(tof_k + 1, 0) + inten // 2
    for _ in range(N_NOISE_EVENTS):
        s = int(rng.integers(0, N_SCANS - 1))
        t = int(rng.integers(5000, 300000))
        d = per_scan.setdefault(s, {})
        d[t] = d.get(t, 0) + int(rng.integers(1, 30))

    empty = (np.array([], dtype=np.uint32), np.array([], dtype=np.uint32))
    scans: List[Scan] = [empty] * N_SCANS
    for s, d in per_scan.items():
        tof = np.array(sorted(d), dtype=np.uint32)
        scans[s] = (tof, np.array([d[t] for t in tof], dtype=np.uint32))
    return scans


def encode_frame(scans: List[Scan]) -> bytes:
    """One ``analysis.tdf_bin`` frame block (compression type 2)."""
    import zstandard  # build-time dependency only

    n = len(scans)
    words: List[int] = [n]
    for s in range(n - 1):
        words.append(2 * len(scans[s][0]))
    for tof, it in scans:
        if len(tof):
            deltas = np.diff(np.concatenate([[-1], tof.astype(np.int64)])).astype(
                np.uint32
            )
            words.extend(np.column_stack([deltas, it]).ravel().tolist())
    arr = np.asarray(words, dtype=np.uint32)
    shuffled = np.ascontiguousarray(arr.view(np.uint8).reshape(-1, 4).T).tobytes()
    payload = zstandard.ZstdCompressor(level=3).compress(shuffled)
    return struct.pack("<II", 8 + len(payload), n) + payload


def build() -> None:
    rng = np.random.default_rng(7)
    frames = [make_frame(k, rng) for k in range(len(GRID))]

    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    OUT_DIR.mkdir()

    offsets: List[int] = []
    pos = 64
    with open(OUT_DIR / "analysis.tdf_bin", "wb") as f:
        f.write(bytes(64))
        for scans in frames:
            block = encode_frame(scans)
            offsets.append(pos)
            f.write(block)
            pos += len(block)

    con = sqlite3.connect(OUT_DIR / "analysis.tdf")
    con.execute("PRAGMA page_size = 512")
    for statement in DDL:
        con.execute(statement)
    con.executemany(
        "INSERT INTO GlobalMetadata VALUES (?, ?)", list(GLOBAL_METADATA.items())
    )
    con.execute(
        "INSERT INTO MzCalibration VALUES ("
        + ",".join("?" * len(MZ_CALIBRATION))
        + ")",
        MZ_CALIBRATION,
    )
    con.execute(
        "INSERT INTO TimsCalibration VALUES ("
        + ",".join("?" * len(TIMS_CALIBRATION))
        + ")",
        TIMS_CALIBRATION,
    )
    con.execute("INSERT INTO PropertyGroups VALUES (1)")
    con.execute("INSERT INTO Segments VALUES (1, 1, ?, 0)", (len(GRID),))
    con.execute(
        "INSERT INTO MaldiFrameLaserInfo VALUES "
        "(1, 'Imaging 20um', 'Single', 0.0, 88.5, 1, 16.0, 16.0, 0, 200, 20.0)"
    )

    expected = {
        "n_scans": N_SCANS,
        "accumulation_time_ms": ACCUMULATION_TIME_MS,
        "frames": [],
    }
    for frame_id, (scans, offset, (x, y)) in enumerate(
        zip(frames, offsets, GRID), start=1
    ):
        pairs = [
            (int(t), s, int(i))
            for s, (tof, it) in enumerate(scans)
            for t, i in zip(tof.tolist(), it.tolist())
        ]
        num_peaks = len(pairs)
        summed = sum(i for _, _, i in pairs)
        max_intensity = max(i for _, _, i in pairs)
        con.execute(
            "INSERT INTO Frames VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (
                frame_id,
                float(frame_id),
                "+",
                20,
                0,
                offset,
                max_intensity,
                summed,
                N_SCANS,
                num_peaks,
                1,
                25.99,
                24.11,
                1,
                1,
                ACCUMULATION_TIME_MS,
                200.0,
                2.7,
            ),
        )
        con.execute(
            "INSERT INTO MaldiFrameInfo VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (
                frame_id,
                0,
                f"R00X{x}Y{y}",
                0,
                x,
                y,
                70.0,
                50,
                1000.0,
                1000.0 * x,
                -1000.0 * y,
                0.0,
                1,
            ),
        )
        by_index: Dict[int, int] = {}
        for t, _, i in pairs:
            by_index[t] = by_index.get(t, 0) + i
        expected["frames"].append(
            {
                "frame": frame_id,
                "x": x,
                "y": y,
                "num_pairs": num_peaks,
                "tic": summed,
                "unique_indices": len(by_index),
                "planted_index": IONS[1][0] + 7 * (frame_id - 1),
                "planted_intensity": by_index[IONS[1][0] + 7 * (frame_id - 1)],
                "pairs": pairs,
            }
        )
    con.commit()
    con.execute("VACUUM")
    con.close()

    EXPECTED_JSON.write_text(
        json.dumps(expected, indent=1) + "\n", encoding="utf-8", newline="\n"
    )
    sizes = {p.name: p.stat().st_size for p in OUT_DIR.iterdir()}
    print(f"Wrote {OUT_DIR} {sizes} and {EXPECTED_JSON.name}")


if __name__ == "__main__":
    build()
