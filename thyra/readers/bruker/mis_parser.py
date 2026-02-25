"""Standalone parser for Bruker FlexImaging .mis XML files.

The .mis file contains teaching point calibration, acquisition area definitions,
and optical image references. It is used by both Rapiflex and timsTOF workflows
for aligning MSI data with optical images.
"""

import logging
import xml.etree.ElementTree as ET  # nosec B405 - parsing trusted local instrument files
from pathlib import Path
from typing import Any, Dict, List

logger = logging.getLogger(__name__)


def parse_mis_file(path: Path) -> Dict[str, Any]:
    """Parse a Bruker FlexImaging .mis XML file.

    Extracts teaching points, area definitions, raster info, and image
    references from the XML structure.

    Args:
        path: Path to the .mis file

    Returns:
        Dictionary with keys: teaching_points, areas, raster, ImageFile,
        OriginalImage, BaseGeometry (all optional depending on file content)
    """
    metadata: Dict[str, Any] = {}

    try:
        tree = ET.parse(path)  # nosec B314
        root = tree.getroot()

        _extract_basic_elements(root, metadata)
        _extract_teaching_points(root, metadata)
        _extract_raster_info(root, metadata)
        _extract_areas(root, metadata)

    except ET.ParseError as e:
        logger.warning(f"Failed to parse .mis file: {e}")

    return metadata


def _extract_basic_elements(root: ET.Element, metadata: Dict[str, Any]) -> None:
    """Extract basic text elements from .mis XML."""
    for elem_name in ["Method", "ImageFile", "OriginalImage", "BaseGeometry"]:
        elem = root.find(f".//{elem_name}")
        if elem is not None and elem.text:
            metadata[elem_name] = elem.text


def _extract_teaching_points(root: ET.Element, metadata: Dict[str, Any]) -> None:
    """Extract teaching point calibration data from .mis XML."""
    teaching_points: List[Dict[str, List[int]]] = []
    for tp in root.findall(".//TeachPoint"):
        if tp.text and ";" in tp.text:
            img_coords, stage_coords = tp.text.split(";")
            img_x, img_y = map(int, img_coords.split(","))
            stage_x, stage_y = map(int, stage_coords.split(","))
            teaching_points.append(
                {"image": [img_x, img_y], "stage": [stage_x, stage_y]}
            )
    if teaching_points:
        metadata["teaching_points"] = teaching_points


def _extract_raster_info(root: ET.Element, metadata: Dict[str, Any]) -> None:
    """Extract raster dimensions from .mis XML."""
    raster_elem = root.find(".//Raster")
    if raster_elem is not None and raster_elem.text:
        parts = raster_elem.text.split(",")
        if len(parts) == 2:
            metadata["raster"] = [int(parts[0]), int(parts[1])]


def _extract_areas(root: ET.Element, metadata: Dict[str, Any]) -> None:
    """Extract Area definitions from .mis XML.

    Each Area defines an acquisition region. Areas may be defined by:
    - Two points (bounding box corners) for rectangular regions
    - Multiple points (polygon) for irregular tissue shapes

    In both cases, the bounding box (min/max of all points) is stored
    as p1 and p2, since alignment only needs the enclosing rectangle.

    Args:
        root: XML root element
        metadata: Dictionary to update with area info
    """
    areas: List[Dict[str, Any]] = []
    for area_elem in root.findall(".//Area"):
        area_name = area_elem.get("Name", "")
        points = area_elem.findall("Point")
        if len(points) >= 2:
            try:
                all_x: List[int] = []
                all_y: List[int] = []
                for p in points:
                    p_text = p.text or ""
                    p_parts = p_text.split(",")
                    all_x.append(int(p_parts[0]))
                    all_y.append(int(p_parts[1]))

                # Bounding box from all points
                areas.append(
                    {
                        "name": area_name,
                        "p1": [min(all_x), min(all_y)],
                        "p2": [max(all_x), max(all_y)],
                    }
                )
            except (ValueError, IndexError) as e:
                logger.warning(f"Failed to parse Area '{area_name}': {e}")
                continue

    if areas:
        metadata["areas"] = areas
        logger.debug(f"Parsed {len(areas)} area definitions from .mis")
