# thyra/readers/mzpeak/mzpeak_reader.py
"""Experimental reader for mzPeak archives.

mzPeak is a ZIP container of Parquet members plus a JSON index. Every member
is stored uncompressed, so each Parquet file is a contiguous byte range that
pyarrow can read through a plain file object.

Container facts this module relies on, all measured against the reference
implementation (HUPO-PSI/mzPeak @ 502c3a4) rather than inferred from the
draft prose:

* ``mzpeak_index.json`` is ``{"files": [...], "metadata": {...}}``. Each file
  entry carries ``name``, ``entity_type``, ``data_kind`` and -- optionally --
  ``column_mapping`` (alias ``metadata_mapping``) and ``parameters``. Members
  are resolved through this index, never by hardcoded filename.
* The reference parser lowercases and trims both vocabularies and accepts the
  space spelling as an alias for the underscore one, falling through to an
  ``Other`` variant rather than erroring. This reader matches that tolerance:
  an unrecognised role is skipped, not fatal.
* Signal data is one struct column, sorted by ``spectrum_index``. In the
  point layout it is
  ``point: struct<spectrum_index: uint64, mz: double, intensity: float>``,
  one point per row. In the chunked layout it is ``chunk``, one run of
  points per row, with the m/z values under the encoding the row names;
  :mod:`.chunk_decoding` turns those rows back into points. The converter
  that builds the public corpus (okohlbacher/mzPeakConverter) writes the
  chunked layout unless told otherwise.
* Signal lives in two members, split by representation. Profile spectra go
  in the ``data_arrays`` member and centroid spectra in the ``peaks``
  member, always (specification, ``docs/schemas/spectra.md``). The reference
  converter writes both members for every input and leaves the unused one
  with zero rows. One member is read per archive; see
  :meth:`MzPeakArchive.signal_kind`.
* There is no spectrum -> row-group map in the container. The KV key
  ``spectrum_array_index`` sounds like one but describes which column holds
  which CV array type. Per-spectrum sizes come instead from the metadata
  member: ``number_of_data_points`` for the ``data_arrays`` member and
  ``number_of_peaks`` for the ``peaks`` member. Each is null on a spectrum
  the member does not hold.
* File-level metadata lives in the Parquet key-value footer of the metadata
  member *and*, inconsistently, in the index's ``metadata`` object. The
  footer was populated on every reference archive; the index object was empty
  on the only imaging one. The footer wins, the index fills gaps.
* Positions are ``opt_IMS_1000050_position_x`` / ``opt_IMS_1000051_position_y``
  and exist only for imaging acquisitions. Thyra is MSI-only, so an archive
  without them is refused.

The format is a v0.9 draft. Everything here is written to fail loudly with a
named file and a named cause rather than to guess.
"""

from __future__ import annotations

import json
import logging
import shutil
import tempfile
import zipfile
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Generator,
    List,
    NamedTuple,
    Optional,
    Set,
    Tuple,
)

import numpy as np
from numpy.typing import NDArray

from ...core.base_reader import BaseMSIReader, OpticalImageLabel
from ...core.mass_axis import MassAxisAccumulator
from ...core.registry import register_reader
from ...errors import ConversionRefused
from .chunk_decoding import (
    ENCODING_GRID,
    LOSSY_ENCODINGS,
    decode_chunks,
    describe_term,
    validate_chunk_columns,
    validate_encodings,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ...alignment import AreaAlignmentResult
    from ...core.base_extractor import MetadataExtractor

logger = logging.getLogger(__name__)

#: Name of the index member. Fixed by the specification; it is the one
#: filename in the container that may be looked up directly, because it is
#: what tells you the name of everything else.
INDEX_MEMBER = "mzpeak_index.json"

#: ZIP local-file-header magic. Used by format detection so a mislabelled
#: ``.mzpeak`` fails detection instead of failing later inside pyarrow.
ZIP_MAGIC = b"PK\x03\x04"

# Imaging CV terms. Resolved by accession rather than by column or parameter
# name: the reference archives spell IMS:1000046 "pixel size (x)" with
# parentheses and IMS:1000047 "pixel size y" without, so name matching finds
# one axis and misses the other.
IMS_POSITION_X = "IMS:1000050"
IMS_POSITION_Y = "IMS:1000051"
IMS_POSITION_Z = "IMS:1000052"

#: Conventional column names for the position terms. Used only as a fallback
#: when the index carries no column mapping, which the schema permits.
DEFAULT_POSITION_X_COLUMN = "opt_IMS_1000050_position_x"
DEFAULT_POSITION_Y_COLUMN = "opt_IMS_1000051_position_y"
DEFAULT_POSITION_Z_COLUMN = "opt_IMS_1000052_position_z"

#: Base of the positions when the archive does not declare one in
#: ``imaging.coordinate_base``: the imzML convention, which every archive
#: seen so far follows.
SPEC_BASE = 1

#: ``data_kind`` of the member holding profile signal.
PROFILE_KIND = "data_arrays"

#: ``data_kind`` of the member holding centroid signal.
CENTROID_KIND = "peaks"

#: Points a batch of chunk rows is sized to hold. A row group of the chunked
#: layout counts rows, and a row is a run of points, so one row group can
#: hold the whole archive (60 million points in one, on a real file).
CHUNK_BATCH_POINTS = 2_000_000

#: Bounds on the rows of one batch, whatever the points per row.
CHUNK_BATCH_ROWS = (16, 65536)

#: How mzpeak-convert marks an image affine fitted to the FlexImaging
#: teaching points, and the direction that affine maps.
TEACH_POINTS_QUALITY = "teach_points"
IMAGE_TO_PIXEL_MAP = "image_px -> ms_px"

#: Image formats the optical image loader reads.
OPTICAL_IMAGE_SUFFIXES = (".tif", ".tiff", ".jpg", ".jpeg", ".png", ".bmp")

#: The suffix an embedded image is copied out under, by the media type the
#: archive declares for it. The loader picks its decoder by suffix, and a
#: slide scanner's TIFF is a TIFF under a name of its own: the glioma example
#: embeds an Aperio ``.svs`` declared ``image/tiff``.
MEDIA_TYPE_SUFFIXES = {
    "image/tiff": ".tif",
    "image/jpeg": ".jpg",
    "image/png": ".png",
    "image/bmp": ".bmp",
}

#: Metadata column giving each spectrum's MS level.
MS_LEVEL_COLUMN = "ms_level"

#: Scans column holding an ion mobility value, which this reader does not
#: convert.
ION_MOBILITY_COLUMN = "ion_mobility_value"

#: Scans column holding the ion injection time, in milliseconds. For a
#: timsTOF frame mzpeak-convert writes the frame's accumulation time there.
INJECTION_TIME_COLUMN = "ion_injection_time"

#: Scans column holding the user parameters of each scan.
SCAN_PARAMETERS_COLUMN = "parameters"

#: The scan parameter in which mzpeak-convert names the acquisition region
#: of a Bruker frame. It has no accession.
REGION_PARAMETER = "acquisition region"

#: Source file format of a timsTOF TDF run.
BRUKER_TDF_FORMAT = "MS:1002817"

#: Bruker's library returns a TDF intensity as the stored count times this,
#: over the frame's accumulation time in milliseconds, rounded half up.
TDF_INTENSITY_NUMERATOR = 100.0

#: Two m/z of one TDF frame this close, relative, are one TOF bin. Bins of a
#: timsTOF are some 1e-6 apart; one bin's copies differ in the last bit.
SAME_BIN_RELATIVE = 1e-12

#: Metadata column recording how many points each spectrum has in a signal
#: member. The specification requires the matching column for whichever
#: member a writer fills.
COUNT_COLUMNS = {
    PROFILE_KIND: "number_of_data_points",
    CENTROID_KIND: "number_of_peaks",
}


def _normalise_token(value: Any) -> str:
    """Fold an index vocabulary token the way the reference parser does.

    ``DataKind::from_str`` and ``EntityType::from_str`` both lowercase and
    trim before matching, and both accept a space where the serialised form
    uses an underscore (``"data arrays"`` for ``"data_arrays"``,
    ``"mass spectrum"`` for ``"spectrum"``). Folding the same way here means
    an archive written by any conforming writer resolves identically.

    Args:
        value: Raw token from the index, or anything else.

    Returns:
        The folded token, or ``""`` when the value is not a string.
    """
    if not isinstance(value, str):
        return ""
    return value.strip().lower().replace(" ", "_")


def _lazy_pyarrow() -> Tuple[Any, Any]:
    """Import pyarrow on first use.

    pyarrow reaches every Thyra install already: it is a hard dependency of
    spatialdata, which is a hard dependency of Thyra. It is imported lazily
    anyway so that ``import thyra`` -- which imports every reader package to
    trigger registration -- does not pay for a 50 MB extension module that
    only mzPeak archives need.

    Returns:
        The ``pyarrow`` and ``pyarrow.parquet`` modules.

    Raises:
        ImportError: If pyarrow is missing or too old.
    """
    try:
        import pyarrow  # noqa: WPS433 - deliberate lazy import
        import pyarrow.parquet as pq  # noqa: WPS433
    except ImportError as exc:  # pragma: no cover - install-shape dependent
        raise ImportError(
            "Reading mzPeak archives requires pyarrow, which is normally "
            "installed as a dependency of spatialdata. Install it with "
            "`pip install pyarrow>=20`."
        ) from exc

    major = int(str(pyarrow.__version__).split(".")[0])
    if major < 20:
        raise ImportError(
            f"Reading mzPeak archives requires pyarrow >= 20, found "
            f"{pyarrow.__version__}. Earlier versions mis-handle the "
            "large_list and struct columns these archives use."
        )
    return pyarrow, pq


def _positions(column: Any) -> Tuple[NDArray[np.int64], NDArray[np.bool_]]:
    """Read one position column, keeping apart the scans that have none.

    A null is a scan that belongs to no pixel. Cast straight to an integer
    it becomes the smallest one there is, and the minimum of the column
    would then put the origin of the grid at that scan.

    Args:
        column: The position column of the scans member, as pyarrow reads it.

    Returns:
        The positions as integers, 0 where there is none, and a mask that is
        ``True`` where the scan has a position.
    """
    if column.null_count == len(column):
        # Nothing to fill, and a column of the null type cannot be filled.
        return np.zeros(len(column), dtype=np.int64), np.zeros(len(column), dtype=bool)
    present = np.asarray(column.is_valid().to_numpy(), dtype=bool)
    values = np.asarray(column.fill_null(0).to_numpy())
    if values.dtype.kind == "f":
        finite = np.isfinite(values)
        present = present & finite
        values = np.where(finite, values, 0)
    return values.astype(np.int64), present


def _whole(value: Any) -> Optional[int]:
    """``value`` as an int when it is a whole number and not a bool."""
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return int(value)


def _box_bounds(box: List[int], offsets: Tuple[int, int]) -> Tuple[int, int, int, int]:
    """A listed box, ``[x_lo, x_hi, y_lo, y_hi]`` in raster indices, as bounds.

    Bounds are ``(x_lo, y_lo, x_hi, y_hi)`` in 0-based pixel coordinates,
    the order the timsTOF reader gives them in.
    """
    return (
        int(box[0] - offsets[0]),
        int(box[2] - offsets[1]),
        int(box[1] - offsets[0]),
        int(box[3] - offsets[1]),
    )


def _parse_region(
    region: Any,
) -> Optional[Tuple[int, Optional[str], int, List[int]]]:
    """One ``bruker_maldi.regions`` entry, or ``None`` when it is incomplete."""
    if not isinstance(region, dict):
        return None
    number = _whole(region.get("region_number"))
    frames = _whole(region.get("frames"))
    box: List[int] = []
    for key in ("x_index", "y_index"):
        span = region.get(key)
        if not isinstance(span, list) or len(span) != 2:
            return None
        ends = [_whole(end) for end in span]
        if ends[0] is None or ends[1] is None:
            return None
        box.extend([ends[0], ends[1]])
    if number is None or frames is None:
        return None
    name = region.get("name")
    return number, name if isinstance(name, str) and name else None, frames, box


def _on_tdf_scale(
    mzs: NDArray[np.float64], counts: NDArray[np.float64], accumulation_ms: float
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """One TDF frame's raw points as the ``.d`` reader yields that frame.

    Bruker's library scales every point of a scan, ``floor(count * 100 /
    accumulation_ms + 0.5)``, and the ``.d`` reader sums the points of all
    scans that fall on one TOF bin. An archive keeps every scan's point, and
    points of one bin decode to one m/z, so they are summed here. "One m/z"
    is to within :data:`SAME_BIN_RELATIVE`: a chunk's first or last value is
    the bound the converter states, which can differ from the decoded value
    of the same bin in the last bit, while neighbouring bins are some 1e-6
    apart.

    Args:
        mzs: m/z of the frame's points.
        counts: The raw count of each point.
        accumulation_ms: The frame's accumulation time.

    Returns:
        Distinct m/z, ascending, and the summed scaled intensity of each.
    """
    scaled = np.floor(counts * TDF_INTENSITY_NUMERATOR / accumulation_ms + 0.5)
    order = np.argsort(mzs, kind="stable")
    ordered = mzs[order]
    starts = np.concatenate(
        ([True], np.diff(ordered) > SAME_BIN_RELATIVE * np.abs(ordered[1:]))
    )
    if starts.all():
        return ordered, scaled[order]
    group = np.cumsum(starts) - 1
    summed = np.bincount(group, weights=scaled[order])
    return ordered[starts], np.asarray(summed, dtype=np.float64)


class MzPeakSpatialIndex(NamedTuple):
    """Everything about the archive that is one value per spectrum.

    Small even for a very large archive -- one row per spectrum, not per
    point -- so it is read once and kept, which is what lets the reader cut a
    row group into spectra without consulting the file again.

    Attributes:
        spectrum_indices: ``spectrum_index`` of each positioned spectrum,
            ascending.
        coordinates: ``(n, 2)`` array of 0-based ``(x, y)`` pixel coordinates.
        raw_positions: ``(n, 2)`` array of the positions as written, before
            normalisation. Kept because coordinate bounds are reported in the
            file's own frame.
        point_counts: Number of points each spectrum has in the signal
            member that is read; 0 for a spectrum that member does not hold.
        offsets: Where index 0 sits in the source's own frame: the base
            subtracted, plus any ``imaging.position_offset`` the writer
            took out.
        complete: Whether every spectrum the archive lists is kept. When
            one is left out (no position, or an MSn spectrum beside MS1),
            the axis, the mass range and the peak count skip it too.
        z_offset: The position z the archive states for its one plane, 0
            when it states none (D30).
        scan_rows: Row of the scans member each spectrum's position was
            read from, so that other values of the same scan can be read.
    """

    spectrum_indices: NDArray[np.int64]
    coordinates: NDArray[np.int64]
    raw_positions: NDArray[np.int64]
    point_counts: NDArray[np.int64]
    offsets: Tuple[int, int]
    complete: bool = True
    z_offset: int = 0
    scan_rows: Optional[NDArray[np.int64]] = None


class MzPeakArchive:
    """Resolved view over one ``.mzpeak`` container.

    Owns the :class:`zipfile.ZipFile` handle and the parsed index, and hands
    out Parquet members by role. Shared by the reader and the metadata
    extractor so the archive is opened and validated exactly once.
    """

    def __init__(self, path: Path):
        """Open an archive and resolve its members.

        Args:
            path: Path to the ``.mzpeak`` file.

        Raises:
            ValueError: If the file is not a ZIP, carries no index, or the
                index is unreadable.
        """
        self.path = Path(path)
        self._zip: Optional[zipfile.ZipFile] = None
        self._parquet_cache: Dict[str, Any] = {}
        self._file_metadata: Optional[Dict[str, Any]] = None
        self._spatial_index: Optional[MzPeakSpatialIndex] = None
        #: The registered image once looked for; ``()`` when there is none.
        self._registered: Optional[Tuple[Any, ...]] = None
        self._null_count: Optional[int] = None
        self._null_count_cached = False
        self._signal_kind: Optional[str] = None
        self._chunk_encodings: Optional[List[str]] = None

        if not zipfile.is_zipfile(self.path):
            raise ConversionRefused(
                f"Not an mzPeak archive (not a ZIP container): {self.path}"
            )

        self._zip = zipfile.ZipFile(self.path)
        self._members = set(self._zip.namelist())
        if INDEX_MEMBER not in self._members:
            raise ConversionRefused(
                f"Not an mzPeak archive (no {INDEX_MEMBER} member): {self.path}"
            )

        try:
            index = json.loads(self._zip.read(INDEX_MEMBER))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise ConversionRefused(
                f"Malformed {INDEX_MEMBER} in {self.path}: {exc}"
            ) from exc
        if not isinstance(index, dict):
            raise ConversionRefused(
                f"Malformed {INDEX_MEMBER} in {self.path}: expected a JSON "
                f"object, found {type(index).__name__}"
            )

        self.index: Dict[str, Any] = index
        self._roles = self._resolve_roles(index)

    def _resolve_roles(self, index: Dict[str, Any]) -> Dict[Tuple[str, str], dict]:
        """Map ``(entity_type, data_kind)`` to the index entry for that member.

        Entries whose member is missing from the ZIP are dropped with a
        warning rather than raising: the index is allowed to describe more
        than a given writer emitted, and only the members this reader
        actually needs are worth failing over.

        Args:
            index: The parsed index document.

        Returns:
            Mapping from folded role pair to the raw index entry.
        """
        roles: Dict[Tuple[str, str], dict] = {}
        entries = index.get("files")
        if not isinstance(entries, list):
            raise ConversionRefused(
                f"Malformed {INDEX_MEMBER} in {self.path}: 'files' must be a "
                f"list, found {type(entries).__name__}"
            )

        for entry in entries:
            if not isinstance(entry, dict):
                continue
            name = entry.get("name")
            if not isinstance(name, str) or name not in self._members:
                logger.warning(
                    "mzPeak index of %s lists member %r which is not in the "
                    "archive; ignoring it",
                    self.path.name,
                    name,
                )
                continue
            key = (
                _normalise_token(entry.get("entity_type")),
                _normalise_token(entry.get("data_kind")),
            )
            roles.setdefault(key, entry)
        return roles

    @property
    def members(self) -> frozenset:
        """Names of every member in the ZIP."""
        return frozenset(self._members)

    def open_member(self, name: str) -> Any:
        """Open one member for reading, as a binary file object."""
        assert self._zip is not None
        return self._zip.open(name)

    def entry(self, entity_type: str, data_kind: str) -> Optional[dict]:
        """Return the index entry for a role, or ``None`` when absent."""
        return self._roles.get(
            (_normalise_token(entity_type), _normalise_token(data_kind))
        )

    def require_entry(self, entity_type: str, data_kind: str) -> dict:
        """Return the index entry for a role, raising when it is absent.

        Raises:
            ValueError: If no member fills the role.
        """
        found = self.entry(entity_type, data_kind)
        if found is None:
            available = sorted(f"{e}/{k}" for e, k in self._roles)
            raise ConversionRefused(
                f"mzPeak archive {self.path} has no "
                f"'{entity_type}/{data_kind}' member. Present roles: "
                f"{', '.join(available) or 'none'}"
            )
        return found

    def parquet(self, entity_type: str, data_kind: str) -> Any:
        """Open the Parquet member filling a role.

        The handle is cached: members are STORED, so pyarrow reads the
        footer once and then seeks within the archive without re-inflating.
        """
        entry = self.require_entry(entity_type, data_kind)
        name = entry["name"]
        if name not in self._parquet_cache:
            _, pq = _lazy_pyarrow()
            assert self._zip is not None
            self._parquet_cache[name] = pq.ParquetFile(self._zip.open(name))
        return self._parquet_cache[name]

    @staticmethod
    def column_for_accession(entry: dict, accession: str) -> Optional[str]:
        """Find the column bound to a CV accession in an index entry.

        ``column_mapping`` carries the binding, but the reference struct
        declares it ``serde(default)`` with the alias ``metadata_mapping``,
        so both spellings occur and both may be absent entirely.

        Args:
            entry: An index file entry.
            accession: The CV accession to look for, e.g. ``"IMS:1000050"``.

        Returns:
            The column path, or ``None`` when the entry binds no such term.
        """
        mapping = entry.get("column_mapping")
        if not isinstance(mapping, list):
            mapping = entry.get("metadata_mapping")
        if not isinstance(mapping, list):
            return None
        for binding in mapping:
            if not isinstance(binding, dict):
                continue
            if binding.get("accession") == accession:
                path = binding.get("path")
                if isinstance(path, str) and path:
                    return path
        return None

    def _signal_rows(self, data_kind: str) -> Optional[int]:
        """Rows a signal member holds, or ``None`` when it is absent."""
        if self.entry("spectrum", data_kind) is None:
            return None
        return int(self.parquet("spectrum", data_kind).metadata.num_rows)

    def signal_kind(self) -> str:
        """Decide which signal member this archive is read from.

        The profile member is read when it holds rows. The centroid member
        is read when the profile member is empty or absent, which is how
        the reference converter writes every centroid input.

        An archive may fill both, for instance profile spectra with their
        picked peaks. The profile member is read then and the centroid
        member is left alone: adding the two would count the same signal
        twice. A spectrum held only in the centroid member is not converted
        in that case. The log gives the row count of each member.

        Returns:
            ``"data_arrays"`` or ``"peaks"``.

        Raises:
            ValueError: If the archive has neither member.
        """
        if self._signal_kind is not None:
            return self._signal_kind

        profile_rows = self._signal_rows(PROFILE_KIND)
        centroid_rows = self._signal_rows(CENTROID_KIND)
        if profile_rows is None and centroid_rows is None:
            # Raises, naming the roles the archive does have.
            self.require_entry("spectrum", PROFILE_KIND)

        if profile_rows:
            kind = PROFILE_KIND
            if centroid_rows:
                logger.warning(
                    "%s holds signal in both members: %d rows of profile "
                    "data and %d rows of centroid data. Reading the profile "
                    "member only; the centroid member is not converted.",
                    self.path.name,
                    profile_rows,
                    centroid_rows,
                )
        elif centroid_rows:
            kind = CENTROID_KIND
            logger.info(
                "%s holds centroid data only; reading its peaks member",
                self.path.name,
            )
        else:
            # Nothing to read either way. Name a member that exists, so
            # the refusal comes from the empty signal and not from here.
            kind = PROFILE_KIND if profile_rows is not None else CENTROID_KIND

        self._signal_kind = kind
        return kind

    def signal(self) -> Any:
        """Open the signal member that :meth:`signal_kind` selected."""
        return self.parquet("spectrum", self.signal_kind())

    def holds_both_representations(self) -> bool:
        """Whether the profile and the centroid member both hold rows."""
        return bool(self._signal_rows(PROFILE_KIND)) and bool(
            self._signal_rows(CENTROID_KIND)
        )

    def _read_point_counts(
        self,
    ) -> Tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]:
        """Read each spectrum's index, its size in the member read, and its level.

        A null count means the spectrum has no points in that member, so it
        becomes 0. Filled before the conversion to numpy: a null in an
        integer column would otherwise arrive as NaN and fail the cast.

        Returns:
            ``(spectrum_index, counts, ms_levels)``, in the metadata
            member's row order. A level the archive does not give is 0.

        Raises:
            ValueError: If the metadata member lacks the count column.
        """
        kind = self.signal_kind()
        column = COUNT_COLUMNS[kind]
        member = self.parquet("spectrum", "metadata")
        names = member.schema_arrow.names
        if column not in names:
            raise ConversionRefused(
                f"{self.path} is read from its '{kind}' member, but its "
                f"spectrum metadata has no '{column}' column to size the "
                f"spectra from."
            )
        has_levels = MS_LEVEL_COLUMN in names
        metadata = member.read(
            columns=["index", column] + ([MS_LEVEL_COLUMN] if has_levels else [])
        )
        spectrum_index = np.asarray(metadata.column("index").to_numpy(), dtype=np.int64)
        counts = np.asarray(
            metadata.column(column).fill_null(0).to_numpy(), dtype=np.int64
        )
        levels = (
            np.asarray(
                metadata.column(MS_LEVEL_COLUMN).fill_null(0).to_numpy(),
                dtype=np.int64,
            )
            if has_levels
            else np.zeros(spectrum_index.size, dtype=np.int64)
        )
        return spectrum_index, counts, levels

    def _first_level_only(self, levels: NDArray[np.int64]) -> NDArray[np.bool_]:
        """Which placed spectra to keep: all but MSn, when MS1 is placed too.

        Spectra that share a pixel are summed into it. An MS2 spectrum
        summed with the MS1 spectrum of its pixel makes a spectrum the
        instrument never measured, so when MS1 spectra are placed, spectra
        of level 2 and up are left out, and the log says how many. Level 0
        means the source stated none (mzpeak-convert writes 0 then); such
        spectra are kept. Without a placed MS1 spectrum all are kept.
        """
        higher = levels >= 2
        if not (levels == 1).any() or not higher.any():
            return np.ones(levels.size, dtype=bool)
        first = ~higher
        other = levels[higher]
        logger.warning(
            "%s places %d MS1 spectra and %d MSn spectra (%s). The MSn "
            "spectra are left out: summed into the same pixels, they would "
            "give spectra the instrument never measured.",
            self.path.name,
            int((levels == 1).sum()),
            int(other.size),
            ", ".join(
                f"level {level}: {count}"
                for level, count in zip(*np.unique(other, return_counts=True))
            ),
        )
        return first

    def layout(self) -> str:
        """Return the physical layout of the signal member that is read.

        Returns:
            ``"point"`` or ``"chunk"``.

        Raises:
            ValueError: If the member carries neither top-level column.
        """
        fields = list(self.signal().schema_arrow.names)
        if "point" in fields:
            return "point"
        if "chunk" in fields:
            return "chunk"
        raise ConversionRefused(
            f"{self.path} has an unrecognised mzPeak data layout: expected a "
            f"top-level 'point' or 'chunk' column, found {fields}."
        )

    def chunk_columns(self) -> List[str]:
        """Children of the ``chunk`` struct of the member that is read."""
        chunk = self.signal().schema_arrow.field("chunk").type
        return [chunk.field(i).name for i in range(chunk.num_fields)]

    def _distinct(self, column: str) -> List[str]:
        """Distinct values of one string child of the ``chunk`` struct."""
        import pyarrow.compute as pc  # noqa: WPS433 - deliberate lazy import

        table = self.signal().read(columns=[column])
        values = table.column("chunk").combine_chunks()
        if hasattr(values, "num_chunks"):
            values = values.chunk(0)
        for name in column.split(".")[1:]:
            values = values.field(name)
        return sorted(str(value or "") for value in pc.unique(values).to_pylist())

    def chunk_encodings(self) -> List[str]:
        """Chunk encodings the member that is read uses, as CV accessions.

        Read from the ``chunk_encoding`` column alone, which is one short
        string per row, so the whole member is covered before a conversion
        starts. Empty for the point layout.

        Raises:
            ConversionRefused: If the member uses an encoding, a grid model
                or an intensity transform that is not decoded.
        """
        if self._chunk_encodings is not None:
            return self._chunk_encodings
        if self.layout() != "chunk":
            self._chunk_encodings = []
            return self._chunk_encodings

        names = self.chunk_columns()
        validate_chunk_columns(names, str(self.path))
        encodings = self._distinct("chunk.chunk_encoding")
        grid_types: List[str] = []
        if ENCODING_GRID in encodings and "mz_grid" in names:
            # A row under another encoding has no grid, and no grid type.
            found = self._distinct("chunk.mz_grid.grid_type")
            grid_types = [grid_type for grid_type in found if grid_type]
        validate_encodings(encodings, grid_types, names, str(self.path))
        self._chunk_encodings = encodings
        return encodings

    def position_columns(self) -> Optional[Tuple[str, str]]:
        """Resolve the scan columns holding the imaging positions.

        Preference is the CV binding in the index, because the column name is
        a convention while the accession is the contract. Falls back to the
        conventional names when the entry carries no mapping, which the
        reference struct permits.

        Returns:
            ``(x_column, y_column)``, or ``None`` for a non-imaging archive.
        """
        entry = self.entry("spectrum", "scans")
        if entry is None:
            return None
        names = set(self.parquet("spectrum", "scans").schema_arrow.names)
        x_column = (
            self.column_for_accession(entry, IMS_POSITION_X)
            or DEFAULT_POSITION_X_COLUMN
        )
        y_column = (
            self.column_for_accession(entry, IMS_POSITION_Y)
            or DEFAULT_POSITION_Y_COLUMN
        )
        if x_column not in names or y_column not in names:
            return None
        return (x_column, y_column)

    def position_z_column(self) -> Optional[str]:
        """Resolve the scan column holding position z, the way x and y are.

        mzpeak-convert writes one when its source states z (an imzML with
        ``IMS:1000052`` on its spectra) and none otherwise.

        Returns:
            The column name, or ``None`` when the archive has none.
        """
        entry = self.entry("spectrum", "scans")
        if entry is None:
            return None
        names = set(self.parquet("spectrum", "scans").schema_arrow.names)
        z_column = (
            self.column_for_accession(entry, IMS_POSITION_Z)
            or DEFAULT_POSITION_Z_COLUMN
        )
        return z_column if z_column in names else None

    def spatial_index(self) -> MzPeakSpatialIndex:
        """Read positions and per-spectrum point counts.

        Scans join to spectra on ``source_index``, not on row order: the
        reference schema numbers scans within a spectrum separately
        (``scan_index``) and nothing promises the two members are ordered
        alike.

        A spectrum whose scan has a null position belongs to no pixel and is
        left out, as a spectrum with no scan row is.

        Raises:
            ConversionRefused: If the archive is not an imaging acquisition,
                carries no positioned spectra, or gives a scan one position
                without the other.
        """
        if self._spatial_index is not None:
            return self._spatial_index

        columns = self.position_columns()
        if columns is None:
            raise ConversionRefused(
                f"{self.path} is not an imaging mzPeak archive: its scans "
                f"member declares no {IMS_POSITION_X}/{IMS_POSITION_Y} "
                f"position columns. Thyra converts imaging acquisitions only."
            )
        x_column, y_column = columns
        z_column = self.position_z_column()

        scans_member = self.parquet("spectrum", "scans")
        scans = scans_member.read(
            columns=[
                name for name in ("source_index", x_column, y_column, z_column) if name
            ]
        )
        self._report_ion_mobility(scans_member)
        source = np.asarray(scans.column("source_index").to_numpy(), dtype=np.int64)
        xs, has_x = _positions(scans.column(x_column))
        ys, has_y = _positions(scans.column(y_column))

        # The imaging profile sets both positions of a scan or neither. A
        # scan that belongs to no pixel, a calibration scan for instance, has
        # both null. One without the other places the scan nowhere.
        half = has_x != has_y
        if half.any():
            raise ConversionRefused(
                f"{self.path} gives {int(half.sum())} of {half.size} scans one "
                f"position but not the other ({IMS_POSITION_X} without "
                f"{IMS_POSITION_Y}, or the reverse). A scan on a pixel has "
                f"both, and a scan on no pixel has neither."
            )
        placed = has_x & has_y
        # Which row of the scans member each entry came from, carried through
        # the reordering below so that position z can be read for the scans
        # that end up on a pixel.
        scan_row = np.arange(source.size, dtype=np.int64)

        # One spectrum may own several scans; imaging acquisitions write one.
        # Keeping the first per spectrum means a multi-scan file resolves to a
        # single pixel rather than silently overwriting itself. A scan with a
        # position goes before one without, so a spectrum is only unplaced
        # when none of its scans is on a pixel.
        first_placed = np.argsort(~placed, kind="stable")
        source, xs, ys, placed, scan_row = (
            source[first_placed],
            xs[first_placed],
            ys[first_placed],
            placed[first_placed],
            scan_row[first_placed],
        )
        _, first = np.unique(source, return_index=True)
        source, xs, ys, placed, scan_row = (
            source[first],
            xs[first],
            ys[first],
            placed[first],
            scan_row[first],
        )

        spectrum_index, counts, levels = self._read_point_counts()
        order = np.argsort(spectrum_index, kind="stable")
        spectrum_index, counts, levels = (
            spectrum_index[order],
            counts[order],
            levels[order],
        )
        n_listed = int(spectrum_index.size)

        # A spectrum with no scan row has no pixel, so it is dropped rather
        # than placed at the origin.
        lookup = {int(s): i for i, s in enumerate(source)}
        keep = np.array([int(s) in lookup for s in spectrum_index], dtype=bool)
        if not keep.all():
            logger.warning(
                "%d of %d spectra in %s have no scan row and therefore no "
                "position; they are skipped",
                int((~keep).sum()),
                keep.size,
                self.path.name,
            )
        spectrum_index, counts, levels = (
            spectrum_index[keep],
            counts[keep],
            levels[keep],
        )
        rows = np.array([lookup[int(s)] for s in spectrum_index], dtype=np.int64)

        # A spectrum whose scan has no position belongs to no pixel. It is
        # left out as an unacquired pixel is, and must not reach the minimum
        # below, where it would set the origin of the grid.
        on_pixel = placed[rows]
        if not on_pixel.all():
            logger.warning(
                "%d of %d spectra in %s have no position and belong to no "
                "pixel; they are left out",
                int((~on_pixel).sum()),
                on_pixel.size,
                self.path.name,
            )
        spectrum_index, counts, rows, levels = (
            spectrum_index[on_pixel],
            counts[on_pixel],
            rows[on_pixel],
            levels[on_pixel],
        )
        if spectrum_index.size == 0:
            raise ConversionRefused(f"{self.path} contains no positioned spectra.")

        # After the position filter, so that an MS1 scan on no pixel, a
        # survey scan for one, cannot cost the placed spectra their place.
        first = self._first_level_only(levels)
        spectrum_index, counts, rows = (
            spectrum_index[first],
            counts[first],
            rows[first],
        )
        complete = spectrum_index.size == n_listed

        raw = np.stack([xs[rows], ys[rows]], axis=1)
        z_offset = self._z_offset(scans, z_column, scan_row[rows])
        self._report_shared_pixels(raw)
        bases = self._position_bases(raw)
        coordinates = raw - np.array(bases, dtype=np.int64)
        shift = self._position_offset()

        self._spatial_index = MzPeakSpatialIndex(
            spectrum_indices=spectrum_index,
            coordinates=coordinates,
            raw_positions=raw,
            point_counts=counts,
            offsets=(bases[0] + shift[0], bases[1] + shift[1]),
            complete=complete,
            z_offset=z_offset,
            scan_rows=scan_row[rows],
        )
        return self._spatial_index

    def _placed_scan_column(self, name: str) -> Optional[Any]:
        """One scans column, for the scan of each kept spectrum in order.

        ``None`` when the scans member has no such column.
        """
        scans = self.parquet("spectrum", "scans")
        if name not in scans.schema_arrow.names:
            return None
        rows = self.spatial_index().scan_rows
        if rows is None:
            return None
        import pyarrow as pa  # noqa: WPS433 - deliberate lazy import

        column = scans.read(columns=[name]).column(0).combine_chunks()
        return column.take(pa.array(rows))

    def scan_regions(self) -> Optional[NDArray[np.int64]]:
        """The acquisition region each kept spectrum's scan names, or ``None``.

        mzpeak-convert 0.17 states a Bruker frame's ``RegionNumber`` as the
        scan parameter :data:`REGION_PARAMETER`. ``None`` unless every kept
        spectrum states one as a whole number: a region read for some
        pixels and guessed for others would be neither.
        """
        column = self._placed_scan_column(SCAN_PARAMETERS_COLUMN)
        if column is None or len(column) == 0:
            return None
        import pyarrow as pa  # noqa: WPS433 - deliberate lazy import
        import pyarrow.compute as pc  # noqa: WPS433 - deliberate lazy import

        if not pa.types.is_list(column.type) and not pa.types.is_large_list(
            column.type
        ):
            return None
        flat = column.flatten()
        if not pa.types.is_struct(flat.type):
            return None
        fields = {flat.type.field(i).name for i in range(flat.type.num_fields)}
        if not {"name", "value"} <= fields:
            return None
        value = flat.field("value")
        if (
            not pa.types.is_struct(value.type)
            or value.type.get_field_index("integer") < 0
        ):
            return None
        named = pc.fill_null(pc.equal(flat.field("name"), REGION_PARAMETER), False)
        rows = np.asarray(
            pc.list_parent_indices(column).filter(named).to_numpy(), dtype=np.int64
        )
        numbers = value.field("integer").filter(named)
        numbers = np.asarray(numbers.fill_null(-1).to_numpy(), dtype=np.int64)
        # The first such parameter of a scan is its region.
        rows, first = np.unique(rows, return_index=True)
        regions = np.full(len(column), -1, dtype=np.int64)
        regions[rows] = numbers[first]
        if (regions < 0).any():
            return None
        return regions

    def from_tdf(self) -> bool:
        """Whether the archive names a timsTOF TDF run as its source."""
        description = self.file_level_metadata().get("file_description")
        files = (
            description.get("source_files") if isinstance(description, dict) else None
        )
        for source in files if isinstance(files, list) else ():
            parameters = source.get("parameters") if isinstance(source, dict) else None
            for parameter in parameters if isinstance(parameters, list) else ():
                if (
                    isinstance(parameter, dict)
                    and parameter.get("accession") == BRUKER_TDF_FORMAT
                ):
                    return True
        return False

    def accumulation_times(self) -> Optional[NDArray[np.float64]]:
        """Each kept spectrum's ion injection time in ms, or ``None``.

        ``None`` unless every kept spectrum states a positive one.
        """
        column = self._placed_scan_column(INJECTION_TIME_COLUMN)
        if column is None:
            return None
        times = np.asarray(
            column.cast("double").to_numpy(zero_copy_only=False), dtype=np.float64
        )
        if times.size == 0 or not (np.isfinite(times) & (times > 0)).all():
            return None
        return times

    def tdf_accumulation_times(self) -> Optional[NDArray[np.float64]]:
        """The times that put a TDF run's raw counts on the ``.d``'s scale.

        mzpeak-convert stores a TDF frame's counts as the run holds them,
        and gives the frame's accumulation time as its scan's ion injection
        time. ``None`` unless the source is a TDF run and every kept
        spectrum gives one (D33).
        """
        return self.accumulation_times() if self.from_tdf() else None

    def _z_offset(
        self, scans: Any, z_column: Optional[str], placed_rows: NDArray[np.int64]
    ) -> int:
        """The position z of the archive's plane, 0 when it states none.

        Recorded as the z offset, the way D14 records the smallest z of an
        imzML. An archive without z records 0, as every format without z
        does (D30). The spectra on a pixel all sit on index 0 of z either
        way: this reader yields one plane.

        Args:
            scans: The scans member as read, with the z column when there is
                one.
            z_column: The column holding position z, or ``None``.
            placed_rows: The rows of ``scans`` whose spectra are on a pixel.

        Raises:
            ConversionRefused: If some of those spectra state z and the
                others do not, or if they state more than one z. Spectra on
                several planes would otherwise be summed into one.
        """
        if z_column is None:
            return 0
        zs, has_z = _positions(scans.column(z_column))
        zs, has_z = zs[placed_rows], has_z[placed_rows]
        if not has_z.any():
            return 0
        if not has_z.all():
            raise ConversionRefused(
                f"{self.path} gives {int(has_z.sum())} of {has_z.size} spectra "
                f"on a pixel a position z ({IMS_POSITION_Z}) and the others "
                f"none. A plane is stated for every spectrum or for none."
            )
        planes = np.unique(zs)
        if planes.size > 1:
            raise ConversionRefused(
                f"{self.path} places its spectra on {planes.size} planes "
                f"(position z from {int(planes[0])} to {int(planes[-1])}). "
                f"Thyra reads one plane of an mzPeak archive."
            )
        return int(planes[0])

    def _report_shared_pixels(self, raw: NDArray[np.int64]) -> None:
        """Warn when several spectra sit on one pixel; they are summed there."""
        pixels = np.unique(raw, axis=0).shape[0]
        if pixels < raw.shape[0]:
            logger.warning(
                "%s places %d spectra on %d pixels, so some pixels hold more "
                "than one spectrum. Those are summed into one row of the "
                "store.",
                self.path.name,
                raw.shape[0],
                pixels,
            )

    def _report_ion_mobility(self, scans: Any) -> None:
        """Say so when the archive carries ion mobility, which is not read.

        Either as a value per scan or as an array beside m/z in the signal
        member.
        """
        struct = self.signal().schema_arrow.field(self.layout()).type
        arrays = [struct.field(i).name for i in range(struct.num_fields)]
        mobility = [name for name in arrays if "mobility" in name.lower()]
        if mobility:
            logger.warning(
                "%s holds ion mobility arrays (%s). Thyra reads m/z and "
                "intensity only, so the store holds no ion mobility.",
                self.path.name,
                ", ".join(mobility),
            )
        if ION_MOBILITY_COLUMN not in scans.schema_arrow.names:
            return
        column = scans.read(columns=[ION_MOBILITY_COLUMN]).column(0)
        with_value = len(column) - column.null_count
        if with_value:
            logger.warning(
                "%d of %d scans in %s give an ion mobility value. Thyra "
                "reads m/z and intensity only, so the store holds no ion "
                "mobility.",
                with_value,
                len(column),
                self.path.name,
            )

    def imaging_metadata(self) -> Dict[str, Any]:
        """The archive's ``imaging`` block, or an empty one."""
        imaging = self.file_level_metadata().get("imaging")
        return imaging if isinstance(imaging, dict) else {}

    def registered_entry(self) -> Optional[Dict[str, Any]]:
        """The ``imaging.images`` entry fitted to the teaching points, or ``None``.

        mzpeak-convert 0.17 fits a FlexImaging run's sequence image to its
        spots through the ``.mis`` teaching points, and marks that affine
        ``registration_quality: teach_points``. The first image so marked
        is the one; ``None`` when none is, or its affine is not six finite
        numbers mapping :data:`IMAGE_TO_PIXEL_MAP` with an inverse.
        """
        registered = self.registered_image()
        if registered is None:
            return None
        for entry in self.imaging_metadata().get("images") or ():
            if isinstance(entry, dict) and entry.get("archive_path") == registered[0]:
                return entry
        return None

    def registered_image(self) -> Optional[Tuple[str, NDArray[np.float64]]]:
        """The member registered to the pixels, and its 3x3 affine.

        The affine maps image pixels to pixel positions, ``position_x`` and
        ``position_y``, pixel centre to pixel centre. See
        :meth:`registered_entry` for which image it is.
        """
        if self._registered is not None:
            return self._registered or None
        self._registered = ()
        for entry in self.imaging_metadata().get("images") or ():
            affine = entry.get("affine") if isinstance(entry, dict) else None
            if not isinstance(affine, dict):
                continue
            if affine.get("registration_quality") != TEACH_POINTS_QUALITY:
                continue
            path = entry.get("archive_path")
            matrix = affine.get("matrix")
            full = None
            if (
                affine.get("maps") == IMAGE_TO_PIXEL_MAP
                and isinstance(path, str)
                and path in self.members
                and isinstance(matrix, list)
                and len(matrix) == 6
                and all(
                    isinstance(v, (int, float)) and not isinstance(v, bool)
                    for v in matrix
                )
            ):
                full = np.array(
                    [matrix[0:3], matrix[3:6], [0.0, 0.0, 1.0]], dtype=np.float64
                )
                if not np.isfinite(full).all() or abs(np.linalg.det(full)) < 1e-300:
                    full = None
            if full is None:
                logger.warning(
                    "%s marks the affine of image %r as fitted to the "
                    "teaching points, but does not give it as six numbers "
                    "mapping %s that can be inverted; the image is carried "
                    "over unaligned.",
                    self.path.name,
                    path,
                    IMAGE_TO_PIXEL_MAP,
                )
                return None
            self._registered = (path, full)
            return self._registered
        return None

    def _position_bases(self, raw: NDArray[np.int64]) -> Tuple[int, int]:
        """The ``(x, y)`` positions that become index 0, as D14 rules for imzML.

        The base is the one the archive declares in
        ``imaging.coordinate_base``, 1 when it declares none, or the
        smallest position when that is lower. So a cropped image keeps its
        place, and lands on the grid its imzML lands on. Rebasing on the
        smallest position instead moved it to the corner and cut the grid
        to its bounding box.
        """
        declared = self.imaging_metadata().get("coordinate_base")
        if not isinstance(declared, int) or isinstance(declared, bool):
            declared = SPEC_BASE
        smallest = (int(raw[:, 0].min()), int(raw[:, 1].min()))
        bases = (min(smallest[0], declared), min(smallest[1], declared))
        if bases != (declared, declared):
            logger.info(
                "%s declares its positions start at %d, but the smallest is "
                "x=%d, y=%d. Rebasing on that, so no row or column is lost.",
                self.path.name,
                declared,
                smallest[0],
                smallest[1],
            )
        return bases

    def _position_offset(self) -> Tuple[int, int]:
        """How far the archive shifted its positions, ``(0, 0)`` if not at all.

        A writer may shift positions to start at the base and keep the
        shift in ``imaging.position_offset``: mzpeak-convert 0.16.0 does so
        for Bruker raster indices. Adding it back puts the image in the
        frame of the instrument's raster, where a conversion of the same
        run from its ``.d`` puts it.
        """
        offset = self.imaging_metadata().get("position_offset")
        if not isinstance(offset, dict):
            return (0, 0)
        values = (offset.get("x", 0), offset.get("y", 0))
        if not all(
            isinstance(value, int) and not isinstance(value, bool) for value in values
        ):
            logger.warning(
                "%s gives a position_offset that is not two whole numbers "
                "(%r); it is not applied.",
                self.path.name,
                offset,
            )
            return (0, 0)
        return (int(values[0]), int(values[1]))

    def null_count(self) -> Optional[int]:
        """How many points carry a null m/z, or ``None`` when unknowable.

        Read from the Parquet column statistics of the signal member that
        is read, so it costs a footer lookup rather than a pass over the
        data.

        mzPeak compresses profile spectra by dropping interior runs of zero
        intensity and marking each gap with a *null pair*: two adjacent rows
        whose m/z **and** intensity are both null. On the reference imaging
        archive that is 13,286 of 36,856 rows, in 512 pairs per spectrum.

        The reference reader regenerates the missing m/z from the
        per-spectrum ``mz_delta_model`` polynomial plus a locally estimated
        median spacing. Those regenerated values are extrapolations rather
        than the instrument's own numbers, and every one of them carries
        zero intensity.

        Thyra writes a sparse matrix, in which a zero-intensity point
        occupies no storage and contributes nothing to any downstream
        analysis. Reconstructing the pairs would therefore only add
        approximate channels to the common mass axis that can never hold a
        value, so the reader drops them and reports the count instead.

        In the chunked layout the count is taken from the intensity list.
        Padding is null there in every encoding, while the m/z of a chunk
        may sit in a byte buffer that has no null to count.
        """
        if self._null_count_cached:
            return self._null_count
        self._null_count_cached = True
        chunked = self.layout() == "chunk"

        metadata = self.signal().metadata
        total = 0
        seen = False
        for group in range(metadata.num_row_groups):
            row_group = metadata.row_group(group)
            for column in range(row_group.num_columns):
                chunk = row_group.column(column)
                path = chunk.path_in_schema
                if chunked and not path.startswith("chunk.intensity."):
                    continue
                if not chunked and not path.endswith(".mz"):
                    continue
                statistics = chunk.statistics
                if statistics is None or statistics.null_count is None:
                    self._null_count = None
                    return None
                seen = True
                total += int(statistics.null_count)
        self._null_count = total if seen else None
        return self._null_count

    def file_level_metadata(self) -> Dict[str, Any]:
        """Merge the two places file-level metadata is written.

        The Parquet key-value footer of the metadata member was populated on
        every reference archive. The index's ``metadata`` object was
        populated on two of three -- and empty on the only imaging one, which
        is exactly the file whose pixel size matters. Reading either alone
        loses information, so the footer is taken as authoritative and the
        index fills keys the footer lacks (notably ``version``, which appears
        only in the index).

        Returns:
            Mapping of metadata key to its decoded JSON value.
        """
        if self._file_metadata is not None:
            return self._file_metadata

        merged: Dict[str, Any] = {}

        index_metadata = self.index.get("metadata")
        if isinstance(index_metadata, dict):
            merged.update(index_metadata)

        footer = self.parquet("spectrum", "metadata").metadata.metadata or {}
        for raw_key, raw_value in footer.items():
            key = raw_key.decode("utf-8", "replace")
            if key == "ARROW:schema":
                continue
            try:
                merged[key] = json.loads(raw_value)
            except (json.JSONDecodeError, UnicodeDecodeError):
                # Scalars such as spectrum_count are written as bare text.
                merged[key] = raw_value.decode("utf-8", "replace")

        self._file_metadata = merged
        return merged

    def close(self) -> None:
        """Release the ZIP handle and any cached Parquet readers."""
        self._parquet_cache.clear()
        if self._zip is not None:
            self._zip.close()
            self._zip = None


@register_reader("mzpeak")
class MzPeakReader(BaseMSIReader):
    """Experimental reader for mzPeak (HUPO-PSI) imaging archives.

    Experimental because the container is a v0.9 draft: column names, the
    index vocabulary and the placement of file-level metadata have all moved
    between prototype revisions, and are expected to move again before v1.0.
    The reader is written against the reference implementation at 502c3a4 and
    validates what it depends on, so a drifted archive produces a named error
    rather than silently wrong pixels.

    Both layouts are read, and both come out as the same stream of
    spectra. The chunked layout is decoded in :mod:`.chunk_decoding`, which
    names the encodings it reads and refuses the rest.

    Data is read as processed imzML is: one m/z per point, per-spectrum
    axes. :attr:`has_shared_mass_axis` is always ``False`` and the
    resampling decision tree treats these files exactly as it treats
    processed imzML. That holds for the grid encoding too; see the
    property.

    Profile and centroid archives are both read, each from its own member;
    see :meth:`MzPeakArchive.signal_kind`.
    """

    def __init__(
        self,
        data_path: Path,
        intensity_threshold: Optional[float] = None,
        **kwargs: object,
    ):
        """Initialise the reader.

        Args:
            data_path: Path to the ``.mzpeak`` archive.
            intensity_threshold: Minimum intensity to keep; see
                :class:`~thyra.core.base_reader.BaseMSIReader`.
            **kwargs: Accepted and ignored, for signature parity with the
                other readers.
        """
        if kwargs.get("region") is not None:
            # Accepted and ignored before, which converted every region of
            # an archive that now reports several.
            raise ConversionRefused(
                f"Selecting a region (--region {kwargs['region']}) is not "
                f"supported for mzPeak archives. Convert the whole archive; "
                f"obs['region_number'] tells the regions apart."
            )
        super().__init__(data_path, intensity_threshold=intensity_threshold, **kwargs)
        self._archive: Optional[MzPeakArchive] = None
        self._coordinates: Optional[NDArray[np.int64]] = None
        self._point_counts: Optional[NDArray[np.int64]] = None
        self._spectrum_indices: Optional[NDArray[np.int64]] = None
        self._offsets: Optional[Tuple[int, int]] = None
        #: Accumulation time of each positioned spectrum, set when a TDF
        #: run's raw counts are put on the scale of Bruker's library.
        self._accumulation: Optional[NDArray[np.float64]] = None
        self._common_axis: Optional[NDArray[np.float64]] = None
        #: Folder the embedded images are copied to, made on first use.
        self._image_dir: Optional[Path] = None
        self._image_paths: Optional[List[Path]] = None
        #: What the store calls each copied image, by the copy's path.
        self._image_labels: Dict[Path, OpticalImageLabel] = {}
        self._unregistered_warned = False
        self._regions_read = False
        self._region_cache: Optional[
            Tuple[Dict[Tuple[int, int], int], List[Dict[str, Any]]]
        ] = None
        self._announced = False
        #: Null-pair padding points dropped by the iteration in
        #: progress; ``iter_spectra`` clears it as it starts.
        self._dropped_points = 0

    # ------------------------------------------------------------------
    # Archive setup
    # ------------------------------------------------------------------

    @property
    def archive(self) -> MzPeakArchive:
        """The open archive, opened and validated on first access."""
        if self._archive is None:
            self._archive = MzPeakArchive(self.data_path)
            self._validate_layout()
            self._load_spatial_index()
            if not self._announced:
                version = self._archive.file_level_metadata().get("version", "unknown")
                logger.warning(
                    "mzPeak support is EXPERIMENTAL: %s declares container "
                    "version %s, a draft format. Verify converted output "
                    "before relying on it.",
                    self.data_path.name,
                    version,
                )
                self._announced = True
        return self._archive

    def _validate_layout(self) -> None:
        """Refuse what this reader cannot honestly read, before reading.

        A chunked member is checked for the encodings it uses here, so an
        encoding that is not decoded stops the conversion at the start and
        not after the first pass.

        Raises:
            ConversionRefused: If the member has neither layout, or uses a
                chunk encoding, grid model or intensity transform that is
                not decoded.
        """
        assert self._archive is not None
        self._archive.chunk_encodings()

    def _load_spatial_index(self) -> None:
        """Cache the archive's per-spectrum positions and point counts.

        The archive owns the resolution so the metadata extractor sees the
        same pixels the iteration does; this only keeps the arrays where the
        hot loop can reach them without a dict lookup.
        """
        assert self._archive is not None
        index = self._archive.spatial_index()
        self._spectrum_indices = index.spectrum_indices
        self._coordinates = index.coordinates
        self._point_counts = index.point_counts
        self._offsets = index.offsets
        self._accumulation = self._tdf_accumulation()

    def _tdf_accumulation(self) -> Optional[NDArray[np.float64]]:
        """The accumulation times that put a TDF run on the ``.d``'s scale.

        mzpeak-convert stores a TDF frame's counts as the run holds them.
        Bruker's library, which the ``.d`` reader uses, returns each count
        times 100 over the frame's accumulation time in ms, rounded half up
        (D33). The archive gives that time per scan, so the reader applies
        the same rule. ``None`` for any other source; also ``None``, with a
        warning, for a TDF archive that does not give every time.
        """
        assert self._archive is not None
        times = self._archive.tdf_accumulation_times()
        if times is None:
            if not self._archive.from_tdf():
                return None
            logger.warning(
                "%s was made from a timsTOF TDF run but does not give every "
                "frame's accumulation time (%s). Its intensities are stored "
                "as the raw counts, not on the scale of a conversion from "
                "the .d.",
                self.data_path.name,
                INJECTION_TIME_COLUMN,
            )
            return None
        logger.info(
            "%s holds the raw counts of a timsTOF TDF run; each is put on the "
            "scale of a conversion from the .d (x %g / accumulation time in "
            "ms, rounded half up), and points of one frame that share an m/z "
            "are summed.",
            self.data_path.name,
            TDF_INTENSITY_NUMERATOR,
        )
        return times

    # ------------------------------------------------------------------
    # BaseMSIReader contract
    # ------------------------------------------------------------------

    def _create_metadata_extractor(self) -> "MetadataExtractor":
        """Create the mzPeak metadata extractor."""
        from ...metadata.extractors.mzpeak_extractor import MzPeakMetadataExtractor

        return MzPeakMetadataExtractor(self.archive, self.data_path)

    @property
    def has_shared_mass_axis(self) -> bool:
        """Always ``False``.

        A shared axis means every spectrum holds the same m/z values, so
        the converter reads the first spectrum's axis and applies it to
        every pixel. No mzPeak layout states that.

        The grid encoding comes closest. It gives each chunk the model that
        turns a grid index into m/z, but a spectrum still lists the indices
        it holds, and two spectra on one grid may hold different ones.
        Knowing that they do not takes a pass over every index, which is
        the pass a shared axis exists to save.
        """
        return False

    def get_common_mass_axis(self) -> NDArray[np.float64]:
        """Build the union of every m/z value in the archive.

        Accumulated row group by row group rather than over the whole column
        at once: the m/z column is the bulk of the archive (~80% on real
        vendor data, where per-frame calibration makes almost every value
        unique), so materialising it entirely would cost more memory than the
        conversion that follows.
        """
        if self._common_axis is not None:
            return self._common_axis

        data = self.archive.signal()
        lossy = [e for e in self.archive.chunk_encodings() if e in LOSSY_ENCODINGS]
        if lossy:
            # Measured on a centroid image: 331,701 distinct m/z in the
            # point layout, 48,116,750 under MS-Numpress.
            logger.warning(
                "%s stores m/z under %s, which gives one m/z back a little "
                "differently from pixel to pixel. The axis of every distinct "
                "m/z is therefore far longer than the source's. Convert with "
                "resampling to avoid that.",
                self.data_path.name,
                ", ".join(describe_term(e) for e in lossy),
            )
        # ``np.union1d`` per row group re-copied the entire axis every
        # group, which is O(unique) memory but quadratic work over the
        # archive. The shared accumulator's buffer capacity tracks the axis
        # length instead, so the number of merges is logarithmic in the
        # input (#294).
        #
        # Honoured if the caller passes one, but no default here: only
        # imzML sets one, so no archive that converts today starts being
        # refused.
        #
        # ``total_spectra`` is left unset rather than given the row-group
        # count: the refusal message counts *spectra*, and a row group holds
        # many. Saying "after 3 of 12" about row groups would be a wrong
        # denominator rather than a missing one.
        accumulator = MassAxisAccumulator(max_length=self.max_mass_axis_length)
        for mzs in self._iter_mz_blocks(data):
            # Null-pair padding carries no intensity; excluded so the
            # axis holds only channels that can actually take a value.
            accumulator.add(mzs[~np.isnan(mzs)])

        try:
            axis = accumulator.finish()
        except ConversionRefused as e:
            # Keep this reader's own wording, which names the archive --
            # but only for the empty-source refusal. An unconditional
            # rewrite would relabel a max_mass_axis_length refusal, raised
            # on the final fold, as "the archive has no usable signal
            # data": the opposite of what happened, on an archive holding
            # too much. Same guard solariX uses.
            if "No spectra found" not in str(e) and "Failed to extract" not in str(e):
                raise
            raise ConversionRefused(
                f"{self.data_path} yielded no m/z values; the archive has no "
                f"usable signal data."
            ) from e
        self._common_axis = axis.astype(np.float64, copy=False)
        return self._common_axis

    def _chunk_rows_per_batch(self, data: Any) -> int:
        """Rows of the chunked member to decode at a time.

        Sized from the archive's own average, so that a batch holds about
        :data:`CHUNK_BATCH_POINTS` points whether a row is a handful of
        centroids or a long stretch of a profile spectrum.
        """
        low, high = CHUNK_BATCH_ROWS
        rows = int(data.metadata.num_rows)
        points = int(self._require(self._point_counts).sum())
        if rows == 0 or points == 0:
            return high
        return int(min(max(CHUNK_BATCH_POINTS * rows // points, low), high))

    def _iter_chunk_blocks(self, data: Any) -> Generator[
        Tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Decode the chunked member, a batch of rows at a time."""
        batches = data.iter_batches(
            batch_size=self._chunk_rows_per_batch(data), columns=["chunk"]
        )
        for batch in batches:
            if batch.num_rows == 0:
                continue
            decoded = decode_chunks(batch.column(0), str(self.data_path))
            yield decoded.spectrum_index, decoded.mz, decoded.intensity

    def _iter_mz_blocks(self, data: Any) -> Generator[NDArray[Any], None, None]:
        """Yield the m/z of the member that is read, in stored order.

        Only the m/z of spectra that are converted: one left out, on no
        pixel or of another MS level, would add channels no pixel fills.
        """
        if not self.archive.spatial_index().complete:
            kept = self._require(self._spectrum_indices)
            for indices, mzs, _ in self._iter_blocks(data):
                yield mzs[np.isin(indices, kept, assume_unique=False)]
            return
        if self.archive.layout() == "chunk":
            for _, mzs, _ in self._iter_chunk_blocks(data):
                yield mzs
            return
        for group in range(data.metadata.num_row_groups):
            table = data.read_row_group(group, columns=["point"])
            yield self._point_field(table, "mz")

    def _iter_blocks(
        self, data: Any
    ) -> Generator[Tuple[NDArray[np.int64], NDArray[Any], NDArray[Any]], None, None]:
        """Yield the member that is read as runs of points, in stored order.

        Both layouts come out alike: the ``spectrum_index``, m/z and
        intensity of each point, with ``NaN`` m/z on padding.
        """
        if self.archive.layout() == "chunk":
            yield from self._iter_chunk_blocks(data)
            return
        for group in range(data.metadata.num_row_groups):
            table = data.read_row_group(group, columns=["point"])
            if table.num_rows == 0:
                continue
            yield (
                self._point_field(table, "spectrum_index").astype(np.int64, copy=False),
                self._point_field(table, "mz"),
                self._point_field(table, "intensity"),
            )

    @staticmethod
    def _point_field(table: Any, field: str) -> NDArray[Any]:
        """Pull one child out of the ``point`` struct column as numpy."""
        column = table.column("point").combine_chunks()
        # ChunkedArray.combine_chunks() yields a ChunkedArray of one chunk on
        # some pyarrow versions and a StructArray on others; normalise.
        if hasattr(column, "num_chunks"):
            column = column.chunk(0)
        return column.field(field).to_numpy(zero_copy_only=False)

    def iter_spectra(self) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Yield every positioned spectrum in ``spectrum_index`` order.

        Reads one block at a time and cuts it into spectra in memory, so
        the cost is one pass over the archive regardless of spectrum count.
        A block is a row group of the point layout, or a batch of decoded
        rows of the chunked layout. Blocks split spectra at their
        boundaries, so a spectrum's points are carried across iterations
        and emitted only once the next ``spectrum_index`` appears -- which
        is why the emit happens on transition rather than per block.

        Yields:
            ``((x, y, z), mzs, intensities)`` with 0-based coordinates,
            ``z`` always 0, and both arrays float64.
        """
        for coords, _, mzs, intensities in self._iter_spectra_with_order():
            yield coords, mzs, intensities

    @property
    def has_acquisition_order(self) -> bool:
        """Always: every spectrum carries its ``spectrum_index``."""
        return True

    def iter_spectra_with_acquisition_order(self) -> Generator[
        Tuple[
            Tuple[int, int, int],
            int,
            NDArray[np.float64],
            NDArray[np.float64],
        ],
        None,
        None,
    ]:
        """:meth:`iter_spectra` with each spectrum's ``spectrum_index``.

        The index is the archive's own: the spectrum's 0-based position in
        its spectrum list, as mzML numbers it, and the key the signal rows
        are sorted by. Unpositioned spectra keep theirs, so the numbers can
        have gaps.
        """
        yield from self._iter_spectra_with_order()

    def _iter_spectra_with_order(
        self,
    ) -> Generator[
        Tuple[
            Tuple[int, int, int],
            int,
            NDArray[np.float64],
            NDArray[np.float64],
        ],
        None,
        None,
    ]:
        """The block walk behind both spectrum iterators, with indices."""
        # Per iteration, not per reader. Every conversion reads the source
        # twice (issue #226), and a counter carried across the passes made
        # the second pass report the sum of both -- 12 dropped points, then
        # 24, then 36 on a third iteration of the same file. Resetting here
        # rather than in ``reset()`` also covers a caller that iterates
        # again without one.
        self._dropped_points = 0

        data = self.archive.signal()
        positioned = {
            int(s): i for i, s in enumerate(self._require(self._spectrum_indices))
        }

        pending_index: Optional[int] = None
        pending_mz: List[NDArray[Any]] = []
        pending_intensity: List[NDArray[Any]] = []

        for indices, mzs, intensities in self._iter_blocks(data):
            if indices.size == 0:
                continue
            # Cut at every change of spectrum_index. The column is sorted, so
            # a change is a boundary and nothing needs grouping.
            cuts = np.flatnonzero(np.diff(indices)) + 1
            for start, stop in zip(
                np.concatenate(([0], cuts)), np.concatenate((cuts, [indices.size]))
            ):
                index = int(indices[start])
                if pending_index is not None and index != pending_index:
                    emitted = self._emit(
                        pending_index, pending_mz, pending_intensity, positioned
                    )
                    if emitted is not None:
                        yield emitted[0], pending_index, emitted[1], emitted[2]
                    pending_mz, pending_intensity = [], []
                pending_index = index
                pending_mz.append(mzs[start:stop])
                pending_intensity.append(intensities[start:stop])

        if pending_index is not None:
            emitted = self._emit(
                pending_index, pending_mz, pending_intensity, positioned
            )
            if emitted is not None:
                yield emitted[0], pending_index, emitted[1], emitted[2]

        if self._dropped_points:
            logger.info(
                "Dropped %d null-pair padding points from %s; they carry no "
                "intensity and exist only to mark removed zero runs",
                self._dropped_points,
                self.data_path.name,
            )

    def _emit(
        self,
        spectrum_index: int,
        mz_parts: List[NDArray[Any]],
        intensity_parts: List[NDArray[Any]],
        positioned: Dict[int, int],
    ) -> Optional[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]]
    ]:
        """Assemble one spectrum's carried parts into a yieldable tuple.

        Returns ``None`` for a spectrum with no position, with no points left
        after intensity filtering, or with an empty payload -- all three are
        ordinary rather than exceptional, and the converter's sparse grid
        handles the resulting gaps.
        """
        row = positioned.get(spectrum_index)
        if row is None:
            return None

        mzs = np.concatenate(mz_parts).astype(np.float64, copy=False)
        intensities = np.concatenate(intensity_parts).astype(np.float64, copy=False)

        # Drop null-pair padding before anything else sees it. Both columns
        # are null on those rows, so a surviving NaN would propagate into the
        # mass axis and into every intensity statistic downstream. See
        # MzPeakArchive.null_count() for why filling them is the wrong repair.
        valid = ~np.isnan(mzs)
        if not valid.all():
            self._dropped_points += int((~valid).sum())
            mzs = mzs[valid]
            intensities = intensities[valid]

        if self._accumulation is not None:
            mzs, intensities = _on_tdf_scale(
                mzs, intensities, float(self._accumulation[row])
            )
        mzs, intensities = self._apply_intensity_filter(mzs, intensities)
        if mzs.size == 0:
            return None

        coordinates = self._require(self._coordinates)
        return (
            (int(coordinates[row, 0]), int(coordinates[row, 1]), 0),
            mzs,
            intensities,
        )

    @staticmethod
    def _require(value: Optional[NDArray[Any]]) -> NDArray[Any]:
        """Assert that the spatial index has been loaded."""
        if value is None:
            raise RuntimeError("mzPeak spatial index not loaded; access .archive first")
        return value

    @property
    def n_spectra(self) -> int:
        """Number of positioned spectra in the archive."""
        _ = self.archive
        return int(self._require(self._spectrum_indices).size)

    @property
    def mass_range(self) -> Tuple[float, float]:
        """Observed m/z range, from the per-spectrum metadata columns."""
        return self.get_essential_metadata().mass_range

    def get_total_peak_count(self) -> int:
        """Points carrying a value across all positioned spectra.

        Excludes null-pair padding; see :meth:`MzPeakArchive.null_count`.
        """
        return self.get_essential_metadata().total_peaks

    def get_region_map(self) -> Optional[Dict[Tuple[int, int], int]]:
        """Each pixel's acquisition region, or ``None`` for one region.

        mzPeak has no region column. A Bruker archive from mzpeak-convert
        lists its regions under ``bruker_maldi.regions``; see
        :meth:`_regions`.
        """
        regions = self._regions()
        return regions[0] if regions is not None else None

    def get_region_info(self) -> Optional[list]:
        """Region number, spectrum count and box of each region.

        The same keys, box frame and order as the timsTOF reader gives for
        the ``.d`` the archive was made from. ``None`` for one region.
        """
        regions = self._regions()
        return regions[1] if regions is not None else None

    def _regions(
        self,
    ) -> Optional[Tuple[Dict[Tuple[int, int], int], List[Dict[str, Any]]]]:
        """Read the regions of a Bruker archive.

        From mzpeak-convert 0.17 every scan names its region, and those
        names are read (:meth:`_scan_regions`). An older archive only lists
        each region with the raster indices it spans and its frame count.
        Then a pixel is given the region whose box holds it, and that is
        only trusted when every pixel falls in exactly one box and each
        region holds as many pixels as it lists frames. Otherwise the
        regions are not read and a warning says why.
        """
        if self._regions_read:
            return self._region_cache
        self._regions_read = True
        listed = self._listed_regions()
        named = self.archive.scan_regions()
        if named is not None:
            self._region_cache = self._scan_regions(named, listed)
            return self._region_cache
        if len(listed) <= 1:
            return None

        coordinates = self._require(self._coordinates)
        offsets = self._require_offsets()
        assigned = np.full(coordinates.shape[0], -1, dtype=np.int64)
        overlaps = 0
        info: List[Dict[str, Any]] = []
        for number, name, frames, box in listed:
            x_lo, y_lo, x_hi, y_hi = _box_bounds(box, offsets)
            inside = (
                (coordinates[:, 0] >= x_lo)
                & (coordinates[:, 0] <= x_hi)
                & (coordinates[:, 1] >= y_lo)
                & (coordinates[:, 1] <= y_hi)
            )
            overlaps += int((inside & (assigned >= 0)).sum())
            assigned[inside] = number
            entry: Dict[str, Any] = {
                "region_number": number,
                "n_spectra": frames,
                "bounds": (x_lo, y_lo, x_hi, y_hi),
            }
            if name:
                entry["name"] = name
            info.append(entry)

        counts = {n: int((assigned == n).sum()) for n, _, _, _ in listed}
        mismatched = [n for n, _, frames, _ in listed if counts[n] != frames]
        unplaced = int((assigned < 0).sum())
        if overlaps or unplaced or mismatched:
            logger.warning(
                "%s lists %d acquisition regions, but their boxes do not "
                "place each pixel once (%d pixels in two boxes, %d in none, "
                "pixel counts that differ from the listed frames in regions "
                "%s). All pixels are stored as one region.",
                self.data_path.name,
                len(listed),
                overlaps,
                unplaced,
                mismatched or "none",
            )
            return None

        region_map = {
            (int(x), int(y)): int(n)
            for (x, y), n in zip(coordinates.tolist(), assigned.tolist())
        }
        info.sort(key=lambda entry: (-entry["n_spectra"], entry["region_number"]))
        self._region_cache = (region_map, info)
        return self._region_cache

    def _scan_regions(
        self,
        named: NDArray[np.int64],
        listed: List[Tuple[int, Optional[str], int, List[int]]],
    ) -> Optional[Tuple[Dict[Tuple[int, int], int], List[Dict[str, Any]]]]:
        """Regions from the region each scan names; ``None`` for one region.

        The box and the name of a region come from its listing when it has
        one, as for the ``.d``. A region no listing describes is boxed by
        its own pixels.
        """
        numbers = sorted({int(n) for n in named.tolist()})
        if len(numbers) <= 1:
            return None
        coordinates = self._require(self._coordinates)
        offsets = self._require_offsets()
        boxes = {number: (name, box) for number, name, _, box in listed}
        info: List[Dict[str, Any]] = []
        for number in numbers:
            mine = coordinates[named == number]
            entry: Dict[str, Any] = {
                "region_number": number,
                "n_spectra": int(mine.shape[0]),
            }
            name, box = boxes.get(number, (None, None))
            if box is not None:
                entry["bounds"] = _box_bounds(box, offsets)
            else:
                low, high = mine.min(axis=0), mine.max(axis=0)
                entry["bounds"] = (int(low[0]), int(low[1]), int(high[0]), int(high[1]))
            if name:
                entry["name"] = name
            info.append(entry)
        region_map = {
            (int(x), int(y)): int(n)
            for (x, y), n in zip(coordinates.tolist(), named.tolist())
        }
        info.sort(key=lambda entry: (-entry["n_spectra"], entry["region_number"]))
        return region_map, info

    def _single_region(self) -> Tuple[int, Optional[str]]:
        """Number and name of an archive's one region, for its alignment.

        The region map is ``None`` for one region, as from the ``.d``, but
        the archive still names it: in every scan, or in its one listing.
        """
        listed = self._listed_regions()
        named = self.archive.scan_regions()
        number = int(named[0]) if named is not None and named.size else None
        if number is None and len(listed) == 1:
            number = listed[0][0]
        if number is None:
            return 0, None
        names = {n: name for n, name, _, _ in listed}
        return number, names.get(number)

    def _listed_regions(self) -> List[Tuple[int, Optional[str], int, List[int]]]:
        """``(number, name, frames, [x_lo, x_hi, y_lo, y_hi])`` per listed region.

        The box is in raw raster indices, as mzpeak-convert writes it. A
        region that does not give all of these is skipped with a warning,
        and then no region is read at all, since the rest cannot be checked
        against the pixels.
        """
        block = self.archive.file_level_metadata().get("bruker_maldi")
        regions = block.get("regions") if isinstance(block, dict) else None
        if not isinstance(regions, list):
            return []
        listed = []
        for region in regions:
            parsed = _parse_region(region)
            if parsed is None:
                logger.warning(
                    "%s lists a region without a number, a frame count and "
                    "its x and y index ranges (%r); regions are not read.",
                    self.data_path.name,
                    region,
                )
                return []
            listed.append(parsed)
        return listed

    def _require_offsets(self) -> Tuple[int, int]:
        """The spatial index's offsets, once it is loaded."""
        _ = self.archive
        if self._offsets is None:
            raise RuntimeError("mzPeak spatial index not loaded; access .archive first")
        return self._offsets

    def get_optical_image_paths(self) -> List[Path]:
        """The optical images the archive embeds, extracted to files.

        mzpeak-convert stores an image as an ``image`` member, verbatim.
        The optical image loader reads files, so each one is copied to a
        folder of this reader's own, which :meth:`close` removes. The copy
        takes the suffix of the media type the archive declares, so a slide
        scanner's TIFF is read as the TIFF it is. A member in a format the
        loader does not read is skipped with a warning.

        An image whose affine is marked ``teach_points`` is the alignment
        image, and the pixels are placed on it (:meth:`get_image_alignment`).
        Any other affine is not a registration: ``assumed_full_extent``
        stretches the image over the whole acquisition. Such an image is
        carried over unaligned, and its affine kept, with the rest of the
        ``imaging`` block, in the store's raw metadata.
        """
        if self._image_paths is None:
            self._image_paths = self._extract_images()
        return list(self._image_paths)

    def get_optical_image_label(self, path: Path) -> Optional[OpticalImageLabel]:
        """Name a copied image after its member, not by the vendor rule.

        mzpeak-convert numbers its members ``images/image_0000.<ext>``, and
        the rule for vendor folders would call the first one the high
        resolution scan and the second the derived image, whatever they
        show.
        """
        return self._image_labels.get(Path(path))

    def _image_members(self) -> List[str]:
        """Names of the members the index lists as images, in index order."""
        entries = self.archive.index.get("files")
        names = []
        for entry in entries if isinstance(entries, list) else []:
            if not isinstance(entry, dict):
                continue
            if _normalise_token(entry.get("entity_type")) != "image":
                continue
            name = entry.get("name")
            if isinstance(name, str) and name in self.archive.members:
                names.append(name)
        return names

    def _image_media_types(self) -> Dict[str, str]:
        """The media type the archive declares for each image member.

        From ``imaging.images``, where each entry names its member in
        ``archive_path``.
        """
        entries = self.archive.imaging_metadata().get("images")
        media_types: Dict[str, str] = {}
        for entry in entries if isinstance(entries, list) else []:
            if not isinstance(entry, dict):
                continue
            path, media_type = entry.get("archive_path"), entry.get("media_type")
            if isinstance(path, str) and isinstance(media_type, str):
                media_types[path] = media_type.strip().lower()
        return media_types

    @staticmethod
    def _copy_suffix(name: str, media_type: Optional[str]) -> Optional[str]:
        """The suffix the loader reads a member under, or ``None`` for none.

        The declared media type decides, so a TIFF named ``.svs`` is read as
        a TIFF. Without one, the member's own suffix does.
        """
        if media_type in MEDIA_TYPE_SUFFIXES:
            return MEDIA_TYPE_SUFFIXES[media_type]
        suffix = Path(name).suffix.lower()
        return suffix if suffix in OPTICAL_IMAGE_SUFFIXES else None

    def _extract_images(self) -> List[Path]:
        """Copy every readable image member out of the archive."""
        members = self._image_members()
        if not members:
            return []
        media_types = self._image_media_types()
        paths: List[Path] = []
        taken: Set[str] = set()
        for name in members:
            suffix = self._copy_suffix(name, media_types.get(name))
            if suffix is None:
                logger.warning(
                    "%s embeds the image %s (%s), in a format the optical "
                    "image loader does not read; it is not carried over.",
                    self.data_path.name,
                    name,
                    media_types.get(name, "no media type declared"),
                )
                continue
            if self._image_dir is None:
                self._image_dir = Path(tempfile.mkdtemp(prefix="thyra-mzpeak-images-"))
            # The base name only: a member name is a path inside the ZIP
            # and must not choose where on this disk the copy goes. Unique
            # whatever the case, because the store names the element after
            # it, lowercased.
            stem = candidate = Path(name).stem
            number = 1
            while candidate.lower() in taken:
                candidate = f"{stem}_{number}"
                number += 1
            taken.add(candidate.lower())
            target = self._image_dir / f"{candidate}{suffix}"
            with self.archive.open_member(name) as source, target.open("wb") as sink:
                shutil.copyfileobj(source, sink)
            paths.append(target)
            self._image_labels[target] = OpticalImageLabel(candidate, name)
        if paths:
            registered = self._registration()
            logger.info(
                "%s embeds %d optical image(s); %s",
                self.data_path.name,
                len(paths),
                (
                    f"the pixels are placed on {registered[0]} by its "
                    f"teach_points affine."
                    if registered is not None
                    else "they are carried over without an alignment to the pixels."
                ),
            )
        return paths

    def _registration(self) -> Optional[Tuple[str, NDArray[np.float64]]]:
        """The registered image the store carries, and its 3x3 affine.

        The archive says which image is registered
        (:meth:`MzPeakArchive.registered_image`). The pixels are placed on
        it only when it is an image member in a format the loader reads, so
        that the store holds the image its pixels are placed on.
        """
        registered = self.archive.registered_image()
        if registered is None:
            return None
        path, matrix = registered
        media_type = self._image_media_types().get(path)
        if (
            path not in self._image_members()
            or self._copy_suffix(path, media_type) is None
        ):
            if not self._unregistered_warned:
                self._unregistered_warned = True
                logger.warning(
                    "%s registers the image %s to its pixels, but does not "
                    "embed it as an image Thyra reads; the pixels are not "
                    "placed on it.",
                    self.data_path.name,
                    path,
                )
            return None
        return path, matrix

    def get_primary_optical_image_path(self) -> Optional[Path]:
        """The copy of the image the archive registers to the pixels."""
        registration = self._registration()
        if registration is None:
            return None
        for path in self.get_optical_image_paths():
            label = self._image_labels.get(path)
            if label is not None and label.source_file == registration[0]:
                return path
        return None

    def get_image_alignment(self) -> Optional["AreaAlignmentResult"]:
        """Place the pixels on the registered image, as the ``.d`` route does.

        The archive's affine takes an image pixel to a pixel position. The
        TIC image's cell ``(i, j)`` holds the pixel at position
        ``(base_x + i, base_y + j)``, centred on ``(i + 0.5, j + 0.5)``, so
        the matrix from a cell to the image is the inverse of the archive's
        after that shift (D33). It is the lattice of a FlexImaging
        alignment, and the regions are the ones :meth:`get_region_map`
        reads.
        """
        registration = self._registration()
        if registration is None:
            return None
        from ...alignment.teaching_points import (  # noqa: WPS433 - lazy import
            AreaAlignmentResult,
            LatticeFit,
            RegionMapping,
        )

        _, image_to_position = registration
        index = self.archive.spatial_index()
        coordinates = self._require(self._coordinates)
        base = index.raw_positions[0] - coordinates[0]
        cell_to_position = np.array(
            [
                [1.0, 0.0, float(base[0]) - 0.5],
                [0.0, 1.0, float(base[1]) - 0.5],
                [0.0, 0.0, 1.0],
            ]
        )
        cell_to_image = np.linalg.inv(image_to_position) @ cell_to_position

        first_x, first_y = self._require_offsets()
        regions = self.get_region_map()
        if regions is None:
            only, only_name = self._single_region()
            names = {only: only_name or ""}
            numbers = np.full(coordinates.shape[0], only, dtype=np.int64)
        else:
            names = {
                int(entry["region_number"]): str(entry.get("name") or "")
                for entry in self.get_region_info() or []
            }
            numbers = np.array(
                [regions[(int(x), int(y))] for x, y in coordinates.tolist()],
                dtype=np.int64,
            )
        mappings = []
        for number in sorted(set(numbers.tolist())):
            mine = coordinates[numbers == number]
            low, high = mine.min(axis=0), mine.max(axis=0)
            corners = (
                np.array(
                    [
                        [low[0], low[1], 1.0],
                        [high[0] + 1, low[1], 1.0],
                        [high[0] + 1, high[1] + 1, 1.0],
                        [low[0], high[1] + 1, 1.0],
                    ]
                )
                @ cell_to_image.T
            )
            mappings.append(
                RegionMapping(
                    region_id=int(number),
                    name=names.get(int(number)) or str(number),
                    raster_min_x=int(low[0]) + first_x,
                    raster_max_x=int(high[0]) + first_x,
                    raster_min_y=int(low[1]) + first_y,
                    raster_max_y=int(high[1]) + first_y,
                    image_min_x=int(np.floor(corners[:, 0].min())),
                    image_max_x=int(np.ceil(corners[:, 0].max())),
                    image_min_y=int(np.floor(corners[:, 1].min())),
                    image_max_y=int(np.ceil(corners[:, 1].max())),
                )
            )
        pixel_size = self.get_essential_metadata().pixel_size or (0.0, 0.0)
        lattice = LatticeFit(
            cell_to_image=cell_to_image,
            reference_node=(first_x, first_y),
            axis_signs=(
                1 if cell_to_image[0, 0] >= 0 else -1,
                1 if cell_to_image[1, 1] >= 0 else -1,
            ),
            raster_step_um=(float(pixel_size[0]), float(pixel_size[1])),
            spots_outside=0,
            nodes_unmeasured=0,
            reference_from="the teach_points affine of the mzPeak archive",
        )
        return AreaAlignmentResult(
            region_mappings=mappings,
            first_raster_x=first_x,
            first_raster_y=first_y,
            pos_to_region={
                (int(x) + first_x, int(y) + first_y): int(n)
                for (x, y), n in zip(coordinates.tolist(), numbers.tolist())
            },
            lattice=lattice,
        )

    def close(self) -> None:
        """Close the archive and remove the extracted images."""
        if self._archive is not None:
            self._archive.close()
            self._archive = None
        if self._image_dir is not None:
            shutil.rmtree(self._image_dir, ignore_errors=True)
            self._image_dir = None
            self._image_paths = None
            self._image_labels = {}
