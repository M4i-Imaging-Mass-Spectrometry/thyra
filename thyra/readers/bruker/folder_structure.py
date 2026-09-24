# thyra/readers/bruker/folder_structure.py
"""Bruker MSI folder structure abstraction.

This module provides a lightweight, pure Python abstraction for analyzing
Bruker MSI folder structures. It handles format detection and file discovery
without requiring any SDK dependencies.

Supported formats:
- timsTOF: .d folders containing analysis.tdf or analysis.tsf
- Rapiflex: Folders with .dat, _poslog.txt, and _info.txt files
- solariX: .d folders containing peaks.sqlite and ImagingInfo.xml
"""

import logging
import os
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path, PureWindowsPath
from typing import List, Optional, Tuple

from ...errors import ConversionRefused
from .mis_parser import find_mis_file_for_d_folder, parse_mis_file

logger = logging.getLogger(__name__)


def _identity(path: Path) -> str:
    """A key that is the same for two spellings of one file.

    The folders optical discovery searches overlap -- a ``.d``'s parent is
    usually the analyzed root -- so the same image arrives twice under
    different-looking paths and has to be recognised.

    ``abspath`` and ``normcase``, not ``Path.resolve()``: this is string
    work that cannot fail, where ``resolve()`` opens the file to read its
    real name back, which on this platform is exactly the call that trips
    over a path past the 259-character ceiling (see
    ``docs/troubleshooting.md``). Discovery must not be able to abort a
    conversion over a bystander image. The cost is that two paths reaching
    one file through different links stay distinct; the store then carries
    the image twice, which is a great deal better than not converting.
    """
    return os.path.normcase(os.path.abspath(path))


class BrukerFormat(Enum):
    """Bruker MSI data formats."""

    TIMSTOF = "timstof"
    RAPIFLEX = "rapiflex"
    SOLARIX = "solarix"
    UNKNOWN = "unknown"


@dataclass
class BrukerFolderInfo:
    """Information about a Bruker MSI folder structure.

    Attributes:
        path: Root path of the Bruker data
        format: Detected Bruker format
        data_path: Path to the main data (e.g., .d folder or data folder)
        optical_images: Optical image paths, the one the .mis names first
        primary_optical_image: The image the .mis ``<ImageFile>`` names, if
            it was named and found. Always the head of ``optical_images``.
        teaching_points_file: Path to teaching points file (e.g., .mis)
        metadata_files: Dictionary of metadata file paths
    """

    path: Path
    format: BrukerFormat
    data_path: Path
    optical_images: List[Path] = field(default_factory=list)
    primary_optical_image: Optional[Path] = None
    teaching_points_file: Optional[Path] = None
    metadata_files: dict = field(default_factory=dict)


class BrukerFolderStructure:
    """Lightweight analyzer for Bruker MSI folder structures.

    This class provides format detection and file discovery for Bruker MSI
    data without requiring any SDK. It's designed to be used for:
    - Automatic format detection
    - Finding optical images
    - Locating metadata and alignment files

    No SDK is required - all operations are pure Python file system checks.

    Example:
        >>> folder = BrukerFolderStructure(Path("/path/to/data"))
        >>> info = folder.analyze()
        >>> print(f"Format: {info.format.value}")
        >>> print(f"Optical images: {info.optical_images}")
    """

    # File patterns for Rapiflex format
    RAPIFLEX_PATTERNS = {
        "data": "*.dat",
        "poslog": "*_poslog.txt",
        "info": "*_info.txt",
        "mis": "*.mis",
    }

    # File patterns for timsTOF format
    TIMSTOF_PATTERNS = {
        "tdf": "analysis.tdf",
        "tsf": "analysis.tsf",
        "tdf_bin": "analysis.tdf_bin",
        "tsf_bin": "analysis.tsf_bin",
    }

    # File patterns for solariX (FT-ICR / MRMS) format. Detection keys on
    # the processed peak store plus the per-scan index; ``ser`` (the raw
    # transient block) marks the solariX family but is never read.
    SOLARIX_PATTERNS = {
        "peaks": "peaks.sqlite",
        "imaging_info": "ImagingInfo.xml",
        "ser": "ser",
    }

    # Optical image suffixes Thyra reads, lowercase. Order is not
    # significant: matching is a membership test, and the images a folder
    # yields are sorted by path.
    #
    # These are what FlexImaging actually writes, not what Thyra would like
    # to read. TIFF is the common export, but an acquisition's .mis may name
    # a .jpg -- a real Rapiflex slide does -- and .png and .bmp turn up as
    # well. Filtering on `Path.suffix.lower()` rather than globbing one
    # pattern per spelling is deliberate: a glob list has to carry `*.tif`,
    # `*.TIF`, `*.Tif`, ... to be case-insensitive on Linux, and on Windows,
    # where globs already ignore case, every one of those patterns returns
    # the same file again.
    OPTICAL_IMAGE_SUFFIXES = (".tif", ".tiff", ".jpg", ".jpeg", ".png", ".bmp")

    def __init__(self, path: Path):
        """Initialize folder structure analyzer.

        Args:
            path: Path to analyze (can be .d folder or parent folder)
        """
        self.path = Path(path)
        self._info: Optional[BrukerFolderInfo] = None

    def analyze(self) -> BrukerFolderInfo:
        """Analyze the folder structure and return information.

        Returns:
            BrukerFolderInfo with detected format and file paths

        Raises:
            ValueError: If path doesn't exist
        """
        if not self.path.exists():
            raise ConversionRefused(f"Path does not exist: {self.path}")

        # Cache the result
        if self._info is None:
            self._info = self._analyze_structure()

        return self._info

    def _analyze_structure(self) -> BrukerFolderInfo:
        """Perform the actual folder analysis."""
        # First, detect the format
        fmt, data_path = self._detect_format()

        # Find teaching points file. Before the optical images, because the
        # .mis is what names the one that matters.
        teaching_points_file = self._find_teaching_points_file(data_path)

        # Find optical images
        primary_optical = self._find_mis_optical_image(data_path, teaching_points_file)
        optical_images = self._find_optical_images(data_path, primary_optical)

        # Find other metadata files
        metadata_files = self._find_metadata_files(data_path, fmt)

        return BrukerFolderInfo(
            path=self.path,
            format=fmt,
            data_path=data_path,
            optical_images=optical_images,
            primary_optical_image=primary_optical,
            teaching_points_file=teaching_points_file,
            metadata_files=metadata_files,
        )

    def _detect_format(self) -> Tuple["BrukerFormat", Path]:
        """Detect the Bruker format and return (format, data_path)."""
        # Check if this is a .d folder (timsTOF or solariX)
        if self.path.suffix.lower() == ".d":
            if self._is_timstof_folder(self.path):
                return BrukerFormat.TIMSTOF, self.path
            if self._is_solarix_folder(self.path):
                return BrukerFormat.SOLARIX, self.path

        # Check if this folder contains Rapiflex data
        if self._is_rapiflex_folder(self.path):
            return BrukerFormat.RAPIFLEX, self.path

        # Sorted: this decides WHICH dataset gets converted, and what
        # format it is called. ``glob`` order is the filesystem's
        # business (issue #303), and ``__main__``'s own multi-dataset
        # prompt already sorts the same shape.
        d_folders = sorted(self.path.glob("*.d"))
        for d_folder in d_folders:
            if self._is_timstof_folder(d_folder):
                return BrukerFormat.TIMSTOF, d_folder
            if self._is_solarix_folder(d_folder):
                return BrukerFormat.SOLARIX, d_folder

        # Check subfolders for Rapiflex, in the same sorted order.
        for subdir in sorted(self.path.iterdir()):
            if subdir.is_dir() and self._is_rapiflex_folder(subdir):
                return BrukerFormat.RAPIFLEX, subdir

        return BrukerFormat.UNKNOWN, self.path

    def _is_timstof_folder(self, path: Path) -> bool:
        """Check if path is a timsTOF .d folder."""
        if not path.is_dir():
            return False

        has_tdf = (path / self.TIMSTOF_PATTERNS["tdf"]).exists()
        has_tsf = (path / self.TIMSTOF_PATTERNS["tsf"]).exists()

        return has_tdf or has_tsf

    def _is_solarix_folder(self, path: Path) -> bool:
        """Check if path is a solariX imaging .d folder.

        Requires BOTH the processed peak store and the per-scan index.
        FlexImaging pre-scan directories (``fid`` + ``analysis.baf``) carry
        neither and must fall through to UNKNOWN.
        """
        if not path.is_dir():
            return False

        has_peaks = (path / self.SOLARIX_PATTERNS["peaks"]).exists()
        has_imaging_info = (path / self.SOLARIX_PATTERNS["imaging_info"]).exists()

        return has_peaks and has_imaging_info

    @classmethod
    def is_solarix_without_peaks(cls, path: Path) -> bool:
        """Check if path is a solariX-family .d that lacks processed peaks.

        A ``.d`` holding raw transients (``ser``) and ``ImagingInfo.xml``
        but no ``peaks.sqlite`` is a solariX acquisition Thyra cannot read
        natively; callers use this to name the imzML-export fallback instead
        of reporting a generic detection failure.
        """
        path = Path(path)
        if not path.is_dir():
            return False

        has_ser = (path / cls.SOLARIX_PATTERNS["ser"]).exists()
        has_imaging_info = (path / cls.SOLARIX_PATTERNS["imaging_info"]).exists()
        has_peaks = (path / cls.SOLARIX_PATTERNS["peaks"]).exists()

        return has_ser and has_imaging_info and not has_peaks

    def _is_rapiflex_folder(self, path: Path) -> bool:
        """Check if path is a Rapiflex data folder."""
        if not path.is_dir():
            return False

        has_dat = bool(list(path.glob(self.RAPIFLEX_PATTERNS["data"])))
        has_poslog = bool(list(path.glob(self.RAPIFLEX_PATTERNS["poslog"])))
        has_info = bool(list(path.glob(self.RAPIFLEX_PATTERNS["info"])))

        return has_dat and has_poslog and has_info

    def _search_paths(self, data_path: Path) -> List[Path]:
        """Folders an optical image may sit in, nearest the data first."""
        search_paths = [data_path]
        if data_path != self.path:
            search_paths.append(self.path)
        if data_path.parent != data_path:
            search_paths.append(data_path.parent)
        return [p for p in search_paths if p.exists()]

    def _find_mis_optical_image(
        self, data_path: Path, mis_file: Optional[Path]
    ) -> Optional[Path]:
        """Resolve the optical image the .mis names, if it names one.

        FlexImaging records the alignment image twice. ``<ImageFile>`` is a
        bare filename, relative to the acquisition; ``<OriginalImage>`` is an
        absolute path on the machine that acquired the data, drive letter and
        all, and is provenance only -- it will not exist anywhere else, so it
        is never opened here. Only ``<ImageFile>`` is resolved, and only
        against the folders the acquisition itself occupies; its directory
        part, if it somehow has one, is dropped for the same reason.

        This is what makes the right image win when a folder holds several --
        a real Rapiflex acquisition ships a slide overview and two more
        ``.jpg`` files alongside the one the .mis names.

        Args:
            data_path: Path to the data folder
            mis_file: The .mis found for this data folder, or None

        Returns:
            Path to the named image, or None when there is no .mis, it names
            no image, or the named file is not next to the data.
        """
        if mis_file is None:
            return None

        # parse_mis_file refuses a document that declares XML entities. Let
        # that travel: every Bruker reader parses the same file from its own
        # __init__ and would raise first, so this is only reachable for a
        # caller that went to the folder structure directly -- and the answer
        # it gives has to be the same one.
        image_file = parse_mis_file(mis_file).get("ImageFile", "")
        if not image_file:
            return None

        # PureWindowsPath, not Path: on Linux a backslash is an ordinary
        # filename character, so Path("img\\scan.jpg").name is the whole
        # string and the lookup silently misses. .mis files are written by
        # Windows software and spell paths the Windows way wherever Thyra
        # runs.
        name = PureWindowsPath(str(image_file).strip()).name
        if not name:
            return None

        for search_path in self._search_paths(data_path):
            candidate = search_path / name
            if candidate.is_file():
                logger.debug(f"Optical image named by {mis_file.name}: {candidate}")
                return candidate

        logger.warning(
            f"{mis_file.name} names optical image '{image_file}', which is not "
            f"in the acquisition folder; no image gets the alignment's "
            f"coordinate space, and any others found are carried as they are"
        )
        return None

    def _find_optical_images(
        self, data_path: Path, primary: Optional[Path] = None
    ) -> List[Path]:
        """Find optical images in the folder structure.

        Searches the data folder, the analyzed root and the data folder's
        parent, keeping any file whose suffix is in
        :data:`OPTICAL_IMAGE_SUFFIXES`.

        Args:
            data_path: Path to the data folder
            primary: The image the .mis names, from
                :meth:`_find_mis_optical_image`. It is returned first, and
                is included even when the search would not have reached it.

        Returns:
            Paths to optical images: ``primary`` first if there is one, then
            everything else in sorted order.
        """
        optical_images: List[Path] = []
        seen = set()
        if primary is not None:
            optical_images.append(primary)
            seen.add(_identity(primary))

        found: List[Path] = []
        for search_path in self._search_paths(data_path):
            for candidate in search_path.glob("*"):
                if candidate.suffix.lower() not in self.OPTICAL_IMAGE_SUFFIXES:
                    continue
                identity = _identity(candidate)
                if identity in seen or not candidate.is_file():
                    continue
                seen.add(identity)
                found.append(candidate)
                logger.debug(f"Found optical image: {candidate}")

        return optical_images + sorted(found)

    def _find_teaching_points_file(self, data_path: Path) -> Optional[Path]:
        """Find the teaching points / alignment file.

        For Rapiflex, this is the .mis file.
        For timsTOF, teaching points may be in other locations (TBD).

        When multiple .mis files exist (e.g., multi-dataset slides where
        each .d has its own .mis), the file whose stem matches the .d
        folder stem is preferred. This ensures each dataset uses its own
        Area definitions for correct optical alignment.

        Args:
            data_path: Path to the data folder

        Returns:
            Path to teaching points file, or None if not found
        """
        # One locator, shared with the pixel-pitch lookup in
        # BrukerMetadataExtractor. They used to be separate and searched
        # different directories, so a .mis inside the .d could give this
        # one the areas while the other took the pitch from a different
        # file entirely (issue #303). ``self.path`` is the extra this
        # caller contributes; the first and last entries are the default.
        return find_mis_file_for_d_folder(
            data_path, search_paths=[data_path, self.path, data_path.parent]
        )

    def _find_metadata_files(self, data_path: Path, fmt: BrukerFormat) -> dict:
        """Find metadata files based on format.

        Args:
            data_path: Path to the data folder
            fmt: Detected Bruker format

        Returns:
            Dictionary of metadata file paths
        """
        metadata = {}

        if fmt == BrukerFormat.RAPIFLEX:
            # Rapiflex metadata files
            for name, pattern in self.RAPIFLEX_PATTERNS.items():
                files = list(data_path.glob(pattern))
                if files:
                    metadata[name] = files[0] if len(files) == 1 else files

        elif fmt == BrukerFormat.TIMSTOF:
            # timsTOF metadata files
            for name, pattern in self.TIMSTOF_PATTERNS.items():
                file_path = data_path / pattern
                if file_path.exists():
                    metadata[name] = file_path

        elif fmt == BrukerFormat.SOLARIX:
            # solariX metadata files inside the .d
            for name, pattern in self.SOLARIX_PATTERNS.items():
                file_path = data_path / pattern
                if file_path.exists():
                    metadata[name] = file_path
            # Acquisition method (XML) inside the *.m method directory
            method_files = sorted(data_path.glob("*.m/apexAcquisition.method"))
            if method_files:
                metadata["method"] = method_files[0]

        return metadata

    @classmethod
    def detect_format(cls, path: Path) -> BrukerFormat:
        """Quick format detection without full analysis.

        Args:
            path: Path to check

        Returns:
            Detected BrukerFormat
        """
        analyzer = cls(path)
        fmt, _ = analyzer._detect_format()
        return fmt

    @classmethod
    def is_bruker_data(cls, path: Path) -> bool:
        """Check if path contains Bruker MSI data.

        Args:
            path: Path to check

        Returns:
            True if Bruker data is detected
        """
        fmt = cls.detect_format(path)
        return fmt != BrukerFormat.UNKNOWN
