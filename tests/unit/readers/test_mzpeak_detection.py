"""Format detection for ``.mzpeak`` archives.

Detection is deliberately container-level: extension, ZIP signature, index
member. It does not open the Parquet members, so a non-imaging or
chunked-layout archive still detects as mzPeak and is refused later by the
reader, with an error that can say which of the two it was.
"""

from __future__ import annotations

import zipfile

import pytest

from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from thyra.core.registry import MSIRegistry


@pytest.fixture
def registry():
    """A registry instance for detection calls."""
    return MSIRegistry()


class TestAccepts:
    """Archives that should be recognised."""

    def test_detects_a_valid_archive(self, registry, tmp_path):
        """Extension, ZIP magic and index member all present."""
        archive = build_mzpeak(tmp_path / "valid.mzpeak", grid_spectra(2, 2))
        assert registry.detect_format(archive) == "mzpeak"

    def test_detection_is_case_insensitive(self, registry, tmp_path):
        """``.MZPEAK`` is the same extension."""
        archive = build_mzpeak(tmp_path / "upper.MZPEAK", grid_spectra(2, 2))
        assert registry.detect_format(archive) == "mzpeak"

    def test_non_imaging_archive_still_detects(self, registry, tmp_path):
        """Detection identifies the container; the reader judges the content.

        Refusing here would report "unsupported format" for a file that is a
        perfectly valid mzPeak archive, which is a worse diagnosis than the
        reader's "not an imaging archive".
        """
        archive = build_mzpeak(
            tmp_path / "nonimaging.mzpeak",
            grid_spectra(2, 1),
            include_positions=False,
        )
        assert registry.detect_format(archive) == "mzpeak"


class TestRejects:
    """Things wearing the extension that are not mzPeak archives."""

    def test_rejects_a_file_that_is_not_a_zip(self, registry, tmp_path):
        """The extension alone is not evidence."""
        impostor = tmp_path / "impostor.mzpeak"
        impostor.write_bytes(b"this is not a zip archive at all")

        with pytest.raises(ValueError, match="missing ZIP signature"):
            registry.detect_format(impostor)

    def test_rejects_a_zip_without_the_index_member(self, registry, tmp_path):
        """A ZIP of Parquet files is not an archive without its index.

        The index is what maps members to roles; without it there is nothing
        to resolve against and hardcoding filenames would be a guess.
        """
        bare = tmp_path / "bare.mzpeak"
        with zipfile.ZipFile(bare, "w", compression=zipfile.ZIP_STORED) as handle:
            handle.writestr("spectra_data.parquet", b"not really parquet")

        with pytest.raises(ValueError, match="no mzpeak_index.json member"):
            registry.detect_format(bare)

    def test_rejects_a_directory(self, registry, tmp_path):
        """mzPeak is a single file; the loose Parquet layout is not it."""
        directory = tmp_path / "unpacked.mzpeak"
        directory.mkdir()

        with pytest.raises(ValueError, match="requires a .mzpeak archive file"):
            registry.detect_format(directory)

    def test_rejects_a_missing_path(self, registry, tmp_path):
        """A path that does not exist fails before any format guessing."""
        with pytest.raises(ValueError, match="does not exist"):
            registry.detect_format(tmp_path / "absent.mzpeak")


class TestRegistration:
    """The reader is wired into the registry under its format name."""

    def test_reader_class_is_registered(self):
        """``convert_msi`` resolves the reader through this lookup.

        Deliberately the module-level accessor rather than a fresh registry:
        the decorator registers against the singleton, so a new instance
        would prove nothing about what conversion actually sees.
        """
        from thyra.core.registry import get_reader_class
        from thyra.readers.mzpeak import MzPeakReader

        assert get_reader_class("mzpeak") is MzPeakReader

    def test_extractor_is_registered_for_the_format(self):
        """The dynamic extractor lookup knows the format too."""
        from thyra.metadata.extractors import get_extractor_for_format
        from thyra.metadata.extractors.mzpeak_extractor import MzPeakMetadataExtractor

        assert get_extractor_for_format("mzpeak") is MzPeakMetadataExtractor

    def test_unsupported_format_message_lists_mzpeak(self, registry, tmp_path):
        """A user with an unreadable file is told mzPeak is an option."""
        unknown = tmp_path / "mystery.xyz"
        unknown.write_bytes(b"\x00\x01\x02")

        with pytest.raises(ValueError, match=r"\.mzpeak"):
            registry.detect_format(unknown)
