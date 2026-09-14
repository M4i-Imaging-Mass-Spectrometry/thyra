# thyra/readers/bruker/base_bruker_reader.py
"""Base class for Bruker MSI readers.

This module provides a common base class for all Bruker MSI readers,
implementing shared functionality for optical image handling and
folder structure analysis.
"""

import logging
from abc import abstractmethod
from pathlib import Path
from typing import List, Optional

from ...core.base_reader import BaseMSIReader
from .folder_structure import BrukerFolderInfo, BrukerFolderStructure

logger = logging.getLogger(__name__)


class BrukerBaseMSIReader(BaseMSIReader):
    """Abstract base class for Bruker MSI data readers.

    This class extends BaseMSIReader with Bruker-specific functionality:
    - Automatic folder structure analysis
    - Optical image discovery
    - Teaching points file location

    Subclasses must implement the abstract methods for reading spectral data.

    Attributes:
        folder_info: BrukerFolderInfo with analyzed folder structure
    """

    def __init__(self, data_path: Path, **kwargs):
        """Initialize Bruker reader with folder analysis.

        Args:
            data_path: Path to Bruker data (folder or .d directory)
            **kwargs: Additional reader-specific parameters
        """
        super().__init__(data_path, **kwargs)

        # Analyze folder structure
        self._folder_structure = BrukerFolderStructure(self.data_path)
        self._folder_info: Optional[BrukerFolderInfo] = None

    @property
    def folder_info(self) -> BrukerFolderInfo:
        """Get analyzed folder information (lazy loaded)."""
        if self._folder_info is None:
            self._folder_info = self._folder_structure.analyze()
        return self._folder_info

    def get_optical_image_paths(self) -> List[Path]:
        """Get paths to optical/microscopy images associated with this data.

        Uses the BrukerFolderStructure to find images in the folder
        hierarchy. Searches the data folder and parent folders for every
        format FlexImaging exports (see
        ``BrukerFolderStructure.OPTICAL_IMAGE_SUFFIXES``), not TIFF alone.

        Returns:
            Paths to optical images, the one the .mis names first
        """
        return self.folder_info.optical_images

    def get_primary_optical_image_path(self) -> Optional[Path]:
        """Get the optical image the .mis names as the alignment image.

        This is the file ``<ImageFile>`` points at, resolved against the
        acquisition folder. It is the image the teaching points are stated
        in, so it is the one whose pixel coordinates the alignment puts the
        MSI shapes into; any other image in the folder is a bystander.

        Returns:
            Path to the named image, or None when no .mis names one that is
            actually there.
        """
        return self.folder_info.primary_optical_image

    def get_teaching_points_file(self) -> Optional[Path]:
        """Get path to the teaching points / alignment file.

        For Rapiflex data, this is typically the .mis file containing
        teaching point calibration data.

        Returns:
            Path to teaching points file, or None if not found
        """
        return self.folder_info.teaching_points_file

    @abstractmethod
    def _create_metadata_extractor(self):
        """Create format-specific metadata extractor.

        Must be implemented by subclasses.
        """
        pass

    @abstractmethod
    def get_common_mass_axis(self):
        """Return the common mass axis for all spectra.

        Must be implemented by subclasses.
        """
        pass

    @abstractmethod
    def iter_spectra(self):
        """Iterate through spectra.

        Must be implemented by subclasses.
        """
        pass

    @abstractmethod
    def close(self):
        """Close all open file handles.

        Must be implemented by subclasses.
        """
        pass
