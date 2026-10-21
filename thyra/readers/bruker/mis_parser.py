"""Standalone parser for Bruker FlexImaging .mis XML files.

The .mis file contains teaching point calibration, acquisition area definitions,
and optical image references. It is used by both Rapiflex and timsTOF workflows
for aligning MSI data with optical images.

Parsing goes through ``defusedxml``, which is a hard dependency: see the
import below for why there is no stdlib fallback.
"""

import logging
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Union

# defusedxml, unconditionally, and never xml.etree here. The stdlib parser
# expands entity declarations -- measured on a .mis carrying
# ``<!ENTITY r "5,5">``: xml.etree hands back the expansion and reports a
# 5x5 raster, defusedxml raises EntitiesForbidden -- so a fallback to it is
# not a degraded parse, it is the hole the defusedxml swap was made to
# close. defusedxml is declared in ``[project] dependencies`` with no
# optional-dependencies table anywhere in pyproject.toml, so an install
# without it is broken rather than a supported configuration, and a broken
# install is entitled to the ImportError and its traceback. This is the
# same principle the spatialdata import follows (issue #310); a
# try/except ImportError that warns and carries on is the mirror image of
# the machinery that commit deleted.
import defusedxml.ElementTree as ET
from defusedxml.common import DefusedXmlException

from ...errors import ConversionRefused

if TYPE_CHECKING:
    from xml.etree.ElementTree import Element  # nosec B405 - type hint only

logger = logging.getLogger(__name__)

#: The largest coordinate a ``.mis`` number may hold, in image pixels or
#: stage micrometres. A metre of stage, or a photo a million pixels wide, is
#: far beyond any slide; a value past this is a damaged file, and one past
#: what a float holds ended in an OverflowError inside the alignment.
MAX_COORDINATE = 1e9


def _number(text: str) -> Optional[float]:
    """``text`` as a finite number within :data:`MAX_COORDINATE`, or ``None``."""
    try:
        value = float(text.strip())
    except (ValueError, OverflowError):
        return None
    if not math.isfinite(value) or abs(value) > MAX_COORDINATE:
        return None
    return value


def _whole(value: float) -> Union[int, float]:
    """A whole number as an ``int``, so ``20`` stays ``20`` in the store."""
    return int(value) if value.is_integer() else value


def _pair(text: Optional[str]) -> Optional[List[Union[int, float]]]:
    """``"x,y"`` as two checked numbers, or ``None`` when it is not that."""
    parts = (text or "").split(",")
    if len(parts) != 2:
        return None
    values = [_number(part) for part in parts]
    if any(value is None for value in values):
        return None
    return [_whole(value) for value in values if value is not None]


def area_is_usable(area: Dict[str, Any]) -> bool:
    """Whether a parsed Area has an outline, or is a placeholder.

    An Area that cannot be read keeps its place in the list, so the Areas
    after it still pair with their own regions (region ``N`` is Area ``N``),
    but it has no ``p1``/``p2``: its region has no Area to be placed by.
    """
    return "p1" in area and "p2" in area


def find_mis_file_for_d_folder(
    data_path: Path, search_paths: Optional[Sequence[Path]] = None
) -> Optional[Path]:
    """Locate the .mis file that corresponds to a Bruker .d folder.

    A FlexImaging acquisition typically writes the .mis next to the .d
    folder using the same stem, and Rapiflex writes it inside the data
    folder, so both are searched -- the data folder first.

    This is the one locator. There were four, and they disagreed in two
    ways that changed what a store records (issue #303):

    - **Order.** Three of them took ``list(glob("*.mis"))[0]`` when no
      stem matched, so which file was picked depended on directory
      listing order. The pick decides the pixel pitch
      (``_resolve_pixel_size_um`` prefers the .mis ``<Raster>`` over
      ``BeamScanSize``) and the acquisition areas that ``--region``
      resolves against, so the same dataset could convert two ways on
      two machines.
    - **Where they looked.** This function searched only
      ``data_path.parent`` while ``_find_teaching_points_file`` searched
      the data folder first. A .mis sitting inside the .d was therefore
      visible to the areas lookup and invisible to the pitch lookup, and
      one dataset could take its pitch from one file and its regions
      from another. That needed no unlucky listing order at all.

    Several non-matching candidates are refused rather than guessed at,
    which is what the solariX reader already did. An arbitrary pick
    writes a *wrong* pitch into the store silently; returning ``None``
    falls back to ``BeamScanSize``, which is at least the instrument's
    own answer. A lone non-matching candidate is still accepted, since
    there is nothing to guess between.

    Never raises: the timsTOF reader wraps its caller in
    ``except (ValueError, OSError)``, and :class:`ConversionRefused` is a
    ``ValueError``, so a refusal here would be swallowed and read as
    "no .mis".

    Args:
        data_path: Path to the Bruker .d directory
        search_paths: Directories to search, in order. Defaults to the
            data folder then its parent.

    Returns:
        Path to the matching .mis file, or None if none found.
    """
    if search_paths is None:
        search_paths = [data_path, data_path.parent]
    present = [p for p in search_paths if p.exists()]

    # The stem match is looked for across *every* search path before any
    # fallback runs. Doing it path-by-path let strays in an earlier
    # directory veto the exact match in a later one: two unrelated .mis
    # inside a .d would refuse, and the `sample.mis` sitting next to
    # `sample.d` -- the pick this function is named for -- was never
    # reached. That is a worse answer than the one it replaced.
    for search_path in present:
        matching = search_path / f"{data_path.stem}.mis"
        if matching.exists():
            return matching

    for search_path in present:
        candidates = sorted(search_path.glob("*.mis"))
        if not candidates:
            continue
        if len(candidates) == 1:
            logger.info(
                "No .mis matching stem '%s'; using the only candidate %s in %s",
                data_path.stem,
                candidates[0].name,
                search_path,
            )
            return candidates[0]
        logger.warning(
            "No .mis matching stem '%s' and %d non-matching candidates "
            "in %s -- refusing to guess which one describes this "
            "acquisition",
            data_path.stem,
            len(candidates),
            search_path,
        )
        return None
    return None


def parse_mis_file(path: Path) -> Dict[str, Any]:
    """Parse a Bruker FlexImaging .mis XML file.

    Extracts teaching points, area definitions, raster info, and image
    references from the XML structure.

    A document defusedxml refuses and a document that is simply not
    well-formed are answered differently, on purpose. A malformed or
    truncated .mis costs the optical alignment, which every caller already
    treats as optional -- most acquisitions have no .mis at all -- so it
    stays a warning and an empty result. A refused document is a different
    claim: nothing about that file was read *because Thyra would not read
    it*, and saying so with an empty dict makes a security decision look
    exactly like a .mis that happened to hold nothing. The visible
    consequence was downstream and misleading: no areas, no teaching
    points, no raster, and a later ``--region <name>`` failing as "no such
    region" while the file that caused it went unnamed.

    Args:
        path: Path to the .mis file

    Returns:
        Dictionary with keys: teaching_points, areas, raster, ImageFile,
        OriginalImage, BaseGeometry (all optional depending on file
        content). Empty when the document is not well-formed XML.

    Raises:
        ConversionRefused: If the document declares XML entities or reaches
            for an external reference. Where that lands differs by
            consumer, and only four of the six let it travel untouched.

            The solariX, Rapiflex and timsTOF readers each parse from
            ``__init__``, so it leaves the constructor and reaches
            ``convert_msi``, which prints it once and stops. The timsTOF
            reader's ``except (ValueError, OSError)`` sits above the parse
            and covers the folder-layout lookup only, so it does not
            swallow this.

            ``BrukerMetadataExtractor._resolve_pixel_size_um`` is not a
            constructor and is not uncaught: ``_extract_essential_impl``
            catches it with ``except ConversionRefused: raise``, which is
            what keeps the broad handler under that clause from re-logging
            the refusal at ERROR as "Unexpected error extracting essential
            metadata". It is re-raised unchanged.

            :func:`thyra.preview.preview_msi` is a fifth surface, and the
            commit that wrote this refusal did not enumerate it -- its
            "checked for each consumer" list named four, all of them
            library callers. Preview builds its reader inside ``except
            Exception``, so nothing propagates: the refusal is turned into
            ``MsiPreview.error`` ("Reader construction failed: ...") and
            the preview comes back with ``readable=False``. No reader's
            construction can skip the parse: solariX and timsTOF both run
            it outside their ``metadata_only`` guard, and Rapiflex has no
            such parameter to sit outside of -- ``metadata_only`` is
            absorbed by its ``**kwargs`` and never read -- so it parses
            whenever the .mis is present. A preview therefore reaches this
            function even though it decodes no spectra and, for timsTOF,
            loads no vendor library.

            ``BrukerFolderStructure._find_mis_optical_image`` is a sixth
            surface, added when optical discovery started resolving the
            ``<ImageFile>`` the .mis names instead of globbing for whatever
            images were in the folder. It lets the refusal travel untouched,
            and in practice it is not the one that raises: ``analyze()`` is
            lazy, so every Bruker reader has already parsed the same file
            from its own ``__init__`` and stopped the conversion before
            discovery is ever asked for the optical images.
    """
    metadata: Dict[str, Any] = {}
    root = _parse_root(path)
    if root is None:
        return metadata
    _extract_basic_elements(root, metadata)
    _extract_teaching_points(root, metadata, path.name)
    _extract_reference_point(root, metadata)
    _extract_raster_info(root, metadata, path.name)
    _extract_areas(root, metadata, path.name)
    return metadata


def read_mis_image_file(path: Path) -> Optional[str]:
    """The ``<ImageFile>`` a ``.mis`` names, read without the rest of it.

    Optical discovery needs only this, and a whole parse would log every
    warning about the file's Areas and points a second time: the reader
    that adopted the file has already parsed and reported it.

    Raises:
        ConversionRefused: As :func:`parse_mis_file` does.
    """
    root = _parse_root(path)
    if root is None:
        return None
    metadata: Dict[str, Any] = {}
    _extract_basic_elements(root, metadata)
    image_file = metadata.get("ImageFile")
    return str(image_file) if image_file is not None else None


def _parse_root(path: Path) -> Optional["Element"]:
    """The document's root, or ``None`` (with a warning) when it is not XML.

    Raises:
        ConversionRefused: If the document declares XML entities or reaches
            for an external reference; see :func:`parse_mis_file`.
    """
    try:
        tree = ET.parse(path)  # nosec B314
        root: Optional["Element"] = tree.getroot()
        if root is None:
            # A completed parse always has a root; read like a parse error.
            logger.warning(f"Failed to parse .mis file: {path} has no root element")
        return root

    # DefusedXmlException is a ValueError, ET.ParseError a SyntaxError, so
    # the two clauses are disjoint and their order is presentation only.
    except DefusedXmlException as e:
        raise ConversionRefused(
            f"Refused to read the FlexImaging sequence file {path}: {e}. "
            "The document declares XML entities, or points at an external "
            "resource, and Thyra does not expand either -- an entity can "
            "pull a file off this machine into the acquisition metadata, "
            "or expand until the parser runs out of memory. Nothing was "
            "read from the file, so the acquisition areas, the teaching "
            "points and the raster step it carries are all unavailable and "
            "the conversion cannot use it. Re-export the imaging sequence "
            "from FlexImaging, or -- after reading what the declaration "
            "actually does -- remove the DOCTYPE from the file."
        ) from e

    except ET.ParseError as e:
        logger.warning(f"Failed to parse .mis file: {e}")
        return None


def _extract_basic_elements(root: "Element", metadata: Dict[str, Any]) -> None:
    """Extract basic text elements from .mis XML."""
    for elem_name in ["Method", "ImageFile", "OriginalImage", "BaseGeometry"]:
        elem = root.find(f".//{elem_name}")
        if elem is not None and elem.text:
            metadata[elem_name] = elem.text


def _extract_teaching_points(
    root: "Element", metadata: Dict[str, Any], source: str = "the .mis"
) -> None:
    """Extract teaching point calibration data from .mis XML.

    A point is ``"image x,y;stage x,y"``, each number finite and within
    :data:`MAX_COORDINATE`; a decimal is read as one. A point that is not
    that is left out with a warning naming it, and the others are kept: the
    fit then has one point fewer, which it checks for itself.
    """
    teaching_points: List[Dict[str, List[Union[int, float]]]] = []
    for tp in root.findall(".//TeachPoint"):
        text = (tp.text or "").strip()
        if not text:
            continue
        halves = text.split(";")
        image = _pair(halves[0]) if len(halves) == 2 else None
        stage = _pair(halves[1]) if len(halves) == 2 else None
        if image is None or stage is None:
            logger.warning(
                f"Unreadable TeachPoint {text!r} in {source}: left out of the "
                "optical alignment (expected two finite numbers, a ';', and "
                "two more)"
            )
            continue
        teaching_points.append({"image": image, "stage": stage})
    if teaching_points:
        metadata["teaching_points"] = teaching_points


def _extract_reference_point(root: "Element", metadata: Dict[str, Any]) -> None:
    """Extract the ``<ReferencePoint>``: the image pixel a raster node sits on.

    FlexImaging lays its spots on one lattice through this point (D31), so
    it fixes where the raster falls between whole steps. A value that is not
    a finite number (``inf``, ``nan``) is unreadable too: the lattice fit
    then uses the first teaching point, as for a ``.mis`` without one.
    """
    elem = root.find(".//ReferencePoint")
    if elem is None or not elem.text:
        return
    parts = elem.text.split(",")
    if len(parts) != 2:
        return
    try:
        point = [float(parts[0]), float(parts[1])]
    except ValueError:
        point = []
    if len(point) == 2 and all(math.isfinite(v) for v in point):
        metadata["reference_point"] = point
    else:
        logger.warning(f"Unreadable ReferencePoint in .mis: {elem.text!r}")


def _extract_raster_info(
    root: "Element", metadata: Dict[str, Any], source: str = "the .mis"
) -> None:
    """Extract the raster step, two finite positive numbers, from .mis XML.

    A step that is not that -- zero, negative, one number or three, text --
    is left out with a warning naming it; the pitch then comes from
    elsewhere, and the metadata extractor says from where.
    """
    raster_elem = root.find(".//Raster")
    if raster_elem is None:
        return
    step = _pair(raster_elem.text)
    if step is not None and all(v > 0 for v in step):
        metadata["raster"] = step
        return
    logger.warning(
        f"Unusable Raster {(raster_elem.text or '').strip()!r} in {source}: "
        "a raster step is two numbers above zero, so this one is not used"
    )


def _area_points(
    point_elems: List["Element"],
) -> Optional[List[List[Union[int, float]]]]:
    """An Area's points, or ``None`` when one of them cannot be read."""
    points = [_pair(p.text) for p in point_elems]
    if any(point is None for point in points):
        return None
    return [point for point in points if point is not None]


def _outline_problem(points: List[List[Union[int, float]]]) -> Optional[str]:
    """Why an Area's points outline no region, or ``None`` when they do."""
    if len(points) < 2:
        return f"{len(points)} point(s), and an Area needs at least two"
    xs = [x for x, _ in points]
    ys = [y for _, y in points]
    if max(xs) == min(xs) or max(ys) == min(ys):
        return "no width or no height"
    if len(points) >= 3:
        # A polygon whose points lie on one line encloses nothing, even
        # when its bounding box does.
        twice_area = sum(
            x0 * y1 - x1 * y0
            for (x0, y0), (x1, y1) in zip(points, points[1:] + points[:1])
        )
        if twice_area == 0:
            return "its points lie on one line"
    return None


def _extract_areas(
    root: "Element", metadata: Dict[str, Any], source: str = "the .mis"
) -> None:
    """Extract Area definitions from .mis XML.

    Areas define the image pixel coordinates for each acquisition region.
    Areas may be rectangular (Type=0, 2 points) or polygon (Type=3, N
    points). The bounding box of all points is stored as p1 and p2. The
    outline itself is kept as ``points`` with its ``type``, and the area's
    own ``<Raster>`` step as ``raster``: the lattice fit (D31) needs the
    outline to tell which raster node each region starts on.

    Every ``<Area>`` gets an entry, in file order, because region ``N`` is
    Area ``N``. One that outlines nothing -- a point that is not two finite
    numbers, fewer than two points, no width or height -- keeps its name
    and its place but no outline, with a warning naming it (see
    :func:`area_is_usable`). Dropping it moved every later Area onto the
    region before it.

    Args:
        root: XML root element
        metadata: Dictionary to update with area info
        source: The file's name, for the warnings
    """
    areas: List[Dict[str, Any]] = []
    for area_elem in root.findall(".//Area"):
        area_name = area_elem.get("Name", "")
        points = _area_points(area_elem.findall("Point"))
        problem = (
            "a point that is not two finite numbers"
            if points is None
            else _outline_problem(points)
        )
        if problem is not None or points is None:
            logger.warning(
                f"Area {len(areas)} ('{area_name}') in {source} cannot be used "
                f"({problem}): region {len(areas)} has no Area to be placed by"
            )
            areas.append({"name": area_name, "unusable": problem})
            continue
        xs = [x for x, _ in points]
        ys = [y for _, y in points]
        area: Dict[str, Any] = {
            "name": area_name,
            "p1": [min(xs), min(ys)],
            "p2": [max(xs), max(ys)],
            "points": points,
        }
        area_type = area_elem.get("Type")
        if area_type is not None and area_type.strip().isdigit():
            area["type"] = int(area_type)
        raster_elem = area_elem.find("Raster")
        if raster_elem is not None and raster_elem.text:
            step = raster_elem.text.split(",")
            try:
                values = [float(v) for v in step]
            except ValueError:
                # The outline is still usable without its own step.
                values = []
            if len(values) == 2 and all(math.isfinite(v) and v > 0 for v in values):
                area["raster"] = values
        areas.append(area)

    if areas:
        metadata["areas"] = areas
        logger.debug(f"Parsed {len(areas)} area definitions from .mis")
