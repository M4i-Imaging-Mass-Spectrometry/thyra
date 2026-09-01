"""Construct minimal valid ``.mzpeak`` archives for tests.

There is no mzPeak writer in Thyra and no Rust toolchain in CI, so the test
corpus is built here, in pure Python, from ``zipfile`` and ``pyarrow``.

Everything this module emits is written out literally -- the index JSON is a
dict spelled out below, the column names are string constants, the member
names are hardcoded -- and nothing is derived from
:mod:`thyra.readers.mzpeak`. That separation is the point. The hand-authored
imzML corpus in ``tests/data/fixtures`` exists because pyimzml's parser and
writer "agree on each other's mistakes"; a builder that asked the reader how
to spell things would reintroduce exactly that loop, and every test on top of
it would pass while real archives failed.

The shapes here were measured against the reference implementation
(HUPO-PSI/mzPeak @ 502c3a4) and its shipped sample archives, not inferred
from the draft prose.
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

# Member names as the reference converter writes them. Real readers resolve
# these through the index rather than by name; the fixtures use the reference
# spelling so a reader that (wrongly) hardcodes them still sees valid input
# and its bug shows up somewhere more informative than a missing file.
DATA_MEMBER = "spectra_data.parquet"
METADATA_MEMBER = "spectra_metadata.parquet"
SCANS_MEMBER = "spectra_metadata_scans.parquet"
INDEX_MEMBER = "mzpeak_index.json"

POSITION_X_COLUMN = "opt_IMS_1000050_position_x"
POSITION_Y_COLUMN = "opt_IMS_1000051_position_y"

#: Micrometre unit accession, as the reference archive declares it.
UNIT_MICROMETRE = "UO:0000017"


class Spectrum:
    """One spectrum destined for a fixture archive.

    Attributes:
        x: Position along x, in the file's own (typically 1-based) frame.
        y: Position along y.
        mzs: m/z values, ascending.
        intensities: Intensity values, same length as ``mzs``.
    """

    def __init__(
        self,
        x: int,
        y: int,
        mzs: Sequence[float],
        intensities: Sequence[float],
    ):
        """Store one spectrum's coordinates and arrays."""
        self.x = int(x)
        self.y = int(y)
        self.mzs = np.asarray(mzs, dtype=np.float64)
        self.intensities = np.asarray(intensities, dtype=np.float64)
        if self.mzs.size != self.intensities.size:
            raise ValueError("mzs and intensities must be the same length")


def _point_table(
    spectra: Sequence[Spectrum],
    null_pair_after: Optional[int] = None,
) -> pa.Table:
    """Build the point-layout signal table.

    Args:
        spectra: Spectra in ``spectrum_index`` order.
        null_pair_after: When given, insert a null pair (two rows whose m/z
            and intensity are both null) after this many real points of every
            spectrum, imitating the padding the reference converter writes
            where it has removed a run of zeros.

    Returns:
        A table with the single struct column ``point``.
    """
    indices: List[int] = []
    mzs: List[Optional[float]] = []
    intensities: List[Optional[float]] = []

    for index, spectrum in enumerate(spectra):
        for position, (mz, intensity) in enumerate(
            zip(spectrum.mzs, spectrum.intensities)
        ):
            indices.append(index)
            mzs.append(float(mz))
            intensities.append(float(intensity))
            if null_pair_after is not None and position == null_pair_after - 1:
                # A null *pair* -- the reference writer never emits a lone
                # null, and its reader raises on an unpaired one.
                for _ in range(2):
                    indices.append(index)
                    mzs.append(None)
                    intensities.append(None)

    point = pa.StructArray.from_arrays(
        [
            pa.array(indices, type=pa.uint64()),
            pa.array(mzs, type=pa.float64()),
            pa.array(intensities, type=pa.float32()),
        ],
        fields=[
            pa.field("spectrum_index", pa.uint64()),
            pa.field("mz", pa.float64()),
            pa.field("intensity", pa.float32()),
        ],
    )
    return pa.table({"point": point})


def _chunk_table(spectra: Sequence[Spectrum]) -> pa.Table:
    """Build a chunked-layout signal table.

    Only the schema matters: this exists so a reader can be shown refusing the
    layout by name, and no test reads the values back.
    """
    chunk = pa.StructArray.from_arrays(
        [
            pa.array(range(len(spectra)), type=pa.uint64()),
            pa.array([float(s.mzs[0]) for s in spectra], type=pa.float64()),
            pa.array([float(s.mzs[-1]) for s in spectra], type=pa.float64()),
            pa.array(
                [[float(v) for v in s.mzs] for s in spectra],
                type=pa.large_list(pa.float64()),
            ),
            pa.array(["basic" for _ in spectra], type=pa.string()),
            pa.array(
                [[float(v) for v in s.intensities] for s in spectra],
                type=pa.large_list(pa.float32()),
            ),
        ],
        fields=[
            pa.field("spectrum_index", pa.uint64()),
            pa.field("mz_chunk_start", pa.float64()),
            pa.field("mz_chunk_end", pa.float64()),
            pa.field("mz_chunk_values", pa.large_list(pa.float64())),
            pa.field("chunk_encoding", pa.string()),
            pa.field("intensity", pa.large_list(pa.float32())),
        ],
    )
    return pa.table({"chunk": chunk})


def _metadata_table(
    spectra: Sequence[Spectrum],
    null_pair_after: Optional[int],
    spectrum_representation: str,
) -> pa.Table:
    """Build the one-row-per-spectrum metadata table.

    ``number_of_data_points`` counts stored rows, padding included -- which is
    what the reference archive does, and the reason the reader has to correct
    it before reporting a peak count.
    """
    padding = 0 if null_pair_after is None else 2
    return pa.table(
        {
            "index": pa.array(range(len(spectra)), type=pa.uint64()),
            "id": pa.array(
                [f"Scan={i + 1}" for i in range(len(spectra))],
                type=pa.large_string(),
            ),
            "ms_level": pa.array([1] * len(spectra), type=pa.uint8()),
            "time": pa.array(
                [float(i) for i in range(len(spectra))], type=pa.float64()
            ),
            "spectrum_representation": pa.array(
                [spectrum_representation] * len(spectra), type=pa.string()
            ),
            "lowest_observed_mz": pa.array(
                [float(s.mzs.min()) for s in spectra], type=pa.float64()
            ),
            "highest_observed_mz": pa.array(
                [float(s.mzs.max()) for s in spectra], type=pa.float64()
            ),
            "number_of_data_points": pa.array(
                [int(s.mzs.size) + padding for s in spectra], type=pa.uint64()
            ),
            "total_ion_current": pa.array(
                [float(s.intensities.sum()) for s in spectra], type=pa.float32()
            ),
        }
    )


def _scans_table(spectra: Sequence[Spectrum], include_positions: bool) -> pa.Table:
    """Build the scans table, with or without the imaging position columns.

    ``source_index`` is the join key back to the spectrum; ``scan_index``
    numbers scans within a spectrum and is deliberately not the same column,
    because a reader that joins on row order passes until it meets a file
    where they differ.
    """
    columns: Dict[str, pa.Array] = {
        "source_index": pa.array(range(len(spectra)), type=pa.uint64()),
        "scan_index": pa.array([0] * len(spectra), type=pa.uint64()),
        "scan_start_time": pa.array(
            [float(i) for i in range(len(spectra))], type=pa.float32()
        ),
    }
    if include_positions:
        columns[POSITION_X_COLUMN] = pa.array([s.x for s in spectra], type=pa.uint32())
        columns[POSITION_Y_COLUMN] = pa.array([s.y for s in spectra], type=pa.uint32())
    return pa.table(columns)


def _scan_settings(
    pixel_size: Optional[Tuple[float, float]],
    grid: Optional[Tuple[int, int]],
    unit: Optional[str],
) -> List[dict]:
    """Build ``scan_settings_list``, imitating the reference spelling.

    The CV names are copied verbatim from the reference archive, parentheses
    and all: IMS:1000046 is written ``"pixel size (x)"`` and IMS:1000047
    ``"pixel size y"``. A reader matching on name rather than accession finds
    one axis and misses the other, and these fixtures are what makes that
    visible.
    """
    parameters: List[dict] = []
    if grid is not None:
        parameters += [
            {
                "name": "max count of pixels x",
                "accession": "IMS:1000042",
                "value": grid[0],
                "unit": None,
            },
            {
                "name": "max count of pixels y",
                "accession": "IMS:1000043",
                "value": grid[1],
                "unit": None,
            },
        ]
    if pixel_size is not None:
        parameters += [
            {
                "name": "pixel size (x)",
                "accession": "IMS:1000046",
                "value": pixel_size[0],
                "unit": unit,
            },
            {
                "name": "pixel size y",
                "accession": "IMS:1000047",
                "value": pixel_size[1],
                "unit": unit,
            },
        ]
    if not parameters:
        return []
    return [
        {
            "id": "scansettings1",
            "source_file_refs": [],
            "targets": [],
            "parameters": parameters,
        }
    ]


def _index_document(
    include_positions: bool,
    data_kind: str,
    column_mapping_key: Optional[str],
    index_metadata: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Spell out ``mzpeak_index.json``.

    Args:
        include_positions: Whether the scans entry binds the position terms.
        data_kind: How to spell the signal member's kind. The reference
            parser accepts ``"data_arrays"`` and ``"data arrays"`` alike, so
            both must resolve.
        column_mapping_key: Which key carries the CV bindings --
            ``"column_mapping"``, its serde alias ``"metadata_mapping"``, or
            ``None`` to omit the bindings entirely, which the reference
            struct permits via ``serde(default)``.
        index_metadata: Contents of the index's ``metadata`` object. The
            reference imaging archive leaves this empty and puts everything
            in the Parquet footer; other archives populate both.
    """
    scans_bindings = [
        {
            "name": "scan start time",
            "path": "scan_start_time",
            "accession": "MS:1000016",
            "unit": "UO:0000031",
        }
    ]
    if include_positions:
        scans_bindings += [
            {
                "name": "position x",
                "path": POSITION_X_COLUMN,
                "accession": "IMS:1000050",
                "unit": None,
            },
            {
                "name": "position y",
                "path": POSITION_Y_COLUMN,
                "accession": "IMS:1000051",
                "unit": None,
            },
        ]

    def entry(name: str, kind: str, bindings: List[dict]) -> Dict[str, Any]:
        item: Dict[str, Any] = {
            "name": name,
            "entity_type": "spectrum",
            "data_kind": kind,
        }
        if column_mapping_key is not None:
            item[column_mapping_key] = bindings
        item["parameters"] = []
        return item

    return {
        "files": [
            entry(DATA_MEMBER, data_kind, []),
            entry(METADATA_MEMBER, "metadata", []),
            entry(SCANS_MEMBER, "scans", scans_bindings),
        ],
        "metadata": dict(index_metadata) if index_metadata else {},
    }


def build_mzpeak(
    path: Path,
    spectra: Sequence[Spectrum],
    *,
    pixel_size: Optional[Tuple[float, float]] = (25.0, 25.0),
    pixel_size_unit: Optional[str] = UNIT_MICROMETRE,
    grid: Optional[Tuple[int, int]] = None,
    include_positions: bool = True,
    row_group_size: Optional[int] = None,
    layout: str = "point",
    data_kind: str = "data_arrays",
    column_mapping_key: Optional[str] = "column_mapping",
    index_metadata: Optional[Dict[str, Any]] = None,
    null_pair_after: Optional[int] = None,
    spectrum_representation: str = "profile spectrum",
    footer_metadata: bool = True,
) -> Path:
    """Write one ``.mzpeak`` archive and return its path.

    Args:
        path: Destination file. Parent directories must exist.
        spectra: Spectra, in the order they should be indexed.
        pixel_size: ``(x, y)`` pixel size, or ``None`` to omit the terms and
            exercise the "pixel size not found" path.
        pixel_size_unit: Unit accession for the pixel size terms.
        grid: Declared grid extent (IMS:1000042/43), or ``None`` to omit.
        include_positions: When ``False``, the scans member carries no
            position columns, which is what a non-imaging archive looks like.
        row_group_size: Rows per Parquet row group. Set this below a
            spectrum's point count to force spectra across row-group
            boundaries.
        layout: ``"point"`` or ``"chunk"``.
        data_kind: Spelling for the signal member's ``data_kind``.
        column_mapping_key: Key carrying the CV bindings, or ``None``.
        index_metadata: Contents of the index ``metadata`` object.
        null_pair_after: Insert a null pair after this many points of each
            spectrum.
        spectrum_representation: Value for the per-spectrum column.
        footer_metadata: Whether to write the file-level JSON blobs into the
            metadata member's Parquet key-value footer. ``False`` leaves the
            index ``metadata`` object as the only source, which is how a
            reader that reads just one of the two gets caught.

    Returns:
        ``path``, for convenience.
    """
    if layout == "point":
        data_table = _point_table(spectra, null_pair_after)
    elif layout == "chunk":
        data_table = _chunk_table(spectra)
    else:  # pragma: no cover - programming error in a test
        raise ValueError(f"unknown layout {layout!r}")

    scan_settings = _scan_settings(pixel_size, grid, pixel_size_unit)
    file_metadata = {
        "file_description": {
            "contents": [
                {
                    "name": spectrum_representation,
                    "accession": (
                        "MS:1000127"
                        if spectrum_representation == "centroid spectrum"
                        else "MS:1000128"
                    ),
                    "value": None,
                    "unit": None,
                }
            ],
            "source_files": [],
        },
        "scan_settings_list": scan_settings,
        "instrument_configuration_list": [],
        "software_list": [],
        "run": {"id": "fixture", "start_time": None},
    }

    metadata_table = _metadata_table(spectra, null_pair_after, spectrum_representation)
    if footer_metadata:
        metadata_table = metadata_table.replace_schema_metadata(
            {key: json.dumps(value) for key, value in file_metadata.items()}
        )

    path = Path(path)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        for member, table, group_size in (
            (DATA_MEMBER, data_table, row_group_size),
            (METADATA_MEMBER, metadata_table, None),
            (SCANS_MEMBER, _scans_table(spectra, include_positions), None),
        ):
            sink = pa.BufferOutputStream()
            kwargs: Dict[str, Any] = {"write_statistics": True}
            if group_size is not None:
                kwargs["row_group_size"] = group_size
            pq.write_table(table, sink, **kwargs)
            archive.writestr(member, sink.getvalue().to_pybytes())

        archive.writestr(
            INDEX_MEMBER,
            json.dumps(
                _index_document(
                    include_positions,
                    data_kind,
                    column_mapping_key,
                    index_metadata,
                ),
                indent=2,
            ),
        )
    return path


def grid_spectra(
    n_x: int,
    n_y: int,
    n_points: int = 8,
    *,
    base: int = 1,
    skip: Sequence[Tuple[int, int]] = (),
) -> List[Spectrum]:
    """Generate a rectangular acquisition with deterministic arrays.

    Args:
        n_x: Pixels along x.
        n_y: Pixels along y.
        n_points: Points per spectrum.
        base: Coordinate origin as written into the file. The reference
            archives are 1-based; ``0`` exercises the normalisation.
        skip: Positions to leave unacquired, so the grid is sparse.

    Returns:
        Spectra in raster order.
    """
    skipped = {(int(x), int(y)) for x, y in skip}
    spectra: List[Spectrum] = []
    for y in range(n_y):
        for x in range(n_x):
            if (x + base, y + base) in skipped:
                continue
            offset = y * n_x + x
            mzs = 100.0 + np.arange(n_points, dtype=np.float64) * 0.5
            intensities = np.arange(1, n_points + 1, dtype=np.float64) + offset * 10.0
            spectra.append(Spectrum(x + base, y + base, mzs, intensities))
    return spectra
