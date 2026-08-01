"""
Tests for the imzML reader.
"""

from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.readers.imzml import imzml_reader as imzml_reader_module
from thyra.readers.imzml.imzml_reader import ImzMLReader


def write_crlf_unindented_imzml(directory: Path, n_spectra: int) -> Path:
    """Write an imzML laid out the way IONTOF SurfaceLab writes it.

    pyimzml's writer indents its output and terminates lines with LF. IONTOF
    emits neither, so the only text between a spectrum's close tag and the
    next one's open tag is a bare carriage-return/newline pair. That exact
    layout is what trips libxml2 (see ``_initialize_parser``), so a fixture
    only guards the bug if it reproduces the formatting.

    Rewriting the whitespace is semantically a no-op -- imzML carries no mixed
    content -- and the .ibd is left untouched, so every offset stays valid.
    """
    imzml_path = directory / f"crlf_{n_spectra}.imzML"
    mzs = np.linspace(100.0, 500.0, 5)
    intensities = np.arange(1.0, 6.0)

    with ImzMLWriter(str(imzml_path), mode="processed") as writer:
        for i in range(n_spectra):
            writer.addSpectrum(mzs, intensities, (i % 8 + 1, i // 8 + 1, 1))

    normalised = imzml_path.read_bytes().replace(b"\r\n", b"\n")
    lines = [line.lstrip(b" \t") for line in normalised.split(b"\n")]
    imzml_path.write_bytes(b"\r\n".join(lines))
    return imzml_path


class TestImzMLReader:
    """Test the ImzML reader functionality."""

    def test_initialization(self, create_minimal_imzml):
        """Test initializing the reader with a valid file."""
        imzml_path, _, _, _ = create_minimal_imzml

        # Test initialization with path as string
        reader1 = ImzMLReader(str(imzml_path))
        assert hasattr(reader1, "parser")

        # Test initialization with path as Path
        reader2 = ImzMLReader(imzml_path)
        assert hasattr(reader2, "parser")

        # Clean up
        reader1.close()
        reader2.close()

    def test_missing_ibd(self, temp_dir):
        """Test error when .ibd file is missing."""
        # Create imzML file without ibd
        imzml_path = temp_dir / "missing.imzML"
        with open(imzml_path, "w") as f:
            f.write("dummy content")

        reader = ImzMLReader(imzml_path)
        with pytest.raises(ValueError):
            # Error should be raised when parser is accessed due to lazy initialization
            reader.get_essential_metadata()

    def test_get_metadata(self, create_minimal_imzml):
        """Test getting metadata from imzML file."""
        imzml_path, _, _, _ = create_minimal_imzml

        reader = ImzMLReader(imzml_path)

        # Test essential metadata
        essential = reader.get_essential_metadata()
        assert essential.source_path == str(imzml_path)
        assert essential.dimensions == (2, 2, 1)
        assert essential.n_spectra == 4

        # Test comprehensive metadata. ``file_mode`` is asserted against the
        # mode the fixture was written in, not against the set of legal
        # values: the extractor used to answer "processed" unconditionally
        # (it read an attribute pyimzml does not define), and a membership
        # test is satisfied by that bug. The discriminating case -- a
        # continuous file -- lives in
        # tests/unit/metadata/extractors/test_imzml_extractor.py.
        comprehensive = reader.get_comprehensive_metadata()
        assert comprehensive.format_specific.get("file_mode") == "processed"
        assert str(imzml_path) in comprehensive.essential.source_path

        reader.close()

    def test_get_dimensions(self, create_minimal_imzml):
        """Test getting dimensions from imzML file."""
        imzml_path, _, _, _ = create_minimal_imzml

        reader = ImzMLReader(imzml_path)
        essential = reader.get_essential_metadata()
        dimensions = essential.dimensions

        # Our test imzML has a 2x2 grid
        assert len(dimensions) == 3  # (x, y, z)
        assert dimensions[0] == 2  # 2 pixels in x
        assert dimensions[1] == 2  # 2 pixels in y
        assert dimensions[2] == 1  # 1 plane in z

        reader.close()

    def test_get_common_mass_axis(self, create_minimal_imzml):
        """Test getting common mass axis from imzML file."""
        imzml_path, _, mzs, _ = create_minimal_imzml

        reader = ImzMLReader(imzml_path)
        mass_axis = reader.get_common_mass_axis()

        # Check that we got a valid mass axis
        assert len(mass_axis) > 0

        # The values should match our input mzs for a 'processed' imzML
        np.testing.assert_allclose(mass_axis, mzs)

        reader.close()

    def test_iter_spectra(self, create_minimal_imzml):
        """Test iterating through spectra."""
        imzml_path, _, mzs, intensities = create_minimal_imzml

        reader = ImzMLReader(imzml_path)

        # Count spectra and check data
        count = 0
        for (
            coords,
            spectrum_mzs,
            spectrum_intensities,
        ) in reader.iter_spectra():
            # Check coordinates format
            assert len(coords) == 3
            x, y, z = coords
            assert x >= 0 and y >= 0 and z >= 0

            # Check that mz values and intensities are provided
            assert len(spectrum_mzs) > 0
            assert len(spectrum_intensities) > 0

            # Verify arrays are valid
            assert isinstance(spectrum_mzs, np.ndarray)
            assert isinstance(spectrum_intensities, np.ndarray)

            # Don't try to index into common_axis which might cause errors
            # if the implementation changed

            count += 1

        # We should have 4 spectra (2x2 grid)
        assert count == 4

        reader.close()

    def test_iter_and_reconstruct(self, create_minimal_imzml):
        """Test iterating through spectra and reconstructing full data."""
        imzml_path, _, mzs, _ = create_minimal_imzml

        reader = ImzMLReader(imzml_path)

        # Get common mass axis (test that it works)
        reader.get_common_mass_axis()

        # Manually collect data similar to what the former 'read' method would do
        coordinates = []
        intensities = []

        # Iterate through all spectra
        for (
            coords,
            spectrum_mzs,
            spectrum_intensities,
        ) in reader.iter_spectra():
            coordinates.append(coords)

            # In a real application, you might need to map these to the common axis
            # For the test, we'll just collect the data
            intensities.append(spectrum_intensities)

        # We should have 4 spectra (2x2 grid)
        assert len(coordinates) == 4
        assert len(intensities) == 4

        # Get dimensions
        essential = reader.get_essential_metadata()
        dimensions = essential.dimensions
        assert dimensions[0] == 2  # width
        assert dimensions[1] == 2  # height

        reader.close()

    def test_close(self, create_minimal_imzml):
        """Test closing the reader."""
        imzml_path, _, _, _ = create_minimal_imzml

        reader = ImzMLReader(imzml_path)
        # Close should work without errors
        reader.close()

    def test_intensity_threshold_filtering(self, create_minimal_imzml):
        """Test that intensity_threshold filters low values during iteration."""
        imzml_path, _, _, intensities = create_minimal_imzml

        # The test data has intensities [100, 200, 300, 400, 500, ...]
        # Set threshold to filter out values below 250
        threshold = 250.0

        # First, read without threshold to get baseline
        reader_no_thresh = ImzMLReader(imzml_path)
        total_values_no_thresh = 0
        for _, _, spectrum_intensities in reader_no_thresh.iter_spectra():
            total_values_no_thresh += len(spectrum_intensities)
        reader_no_thresh.close()

        # Now read with threshold
        reader_with_thresh = ImzMLReader(imzml_path, intensity_threshold=threshold)
        total_values_with_thresh = 0
        for _, _, spectrum_intensities in reader_with_thresh.iter_spectra():
            # All returned intensities should be >= threshold
            assert np.all(spectrum_intensities >= threshold), (
                f"Found intensities below threshold: "
                f"{spectrum_intensities[spectrum_intensities < threshold]}"
            )
            total_values_with_thresh += len(spectrum_intensities)
        reader_with_thresh.close()

        # With threshold, we should have fewer values
        assert total_values_with_thresh < total_values_no_thresh, (
            f"Expected fewer values with threshold. "
            f"No threshold: {total_values_no_thresh}, "
            f"With threshold: {total_values_with_thresh}"
        )

    def test_intensity_threshold_none_returns_all(self, create_minimal_imzml):
        """Test that intensity_threshold=None returns all values."""
        imzml_path, _, _, _ = create_minimal_imzml

        # Read with explicit None threshold
        reader = ImzMLReader(imzml_path, intensity_threshold=None)

        # Should return all 4 spectra with all their values
        count = 0
        for _, _, spectrum_intensities in reader.iter_spectra():
            assert len(spectrum_intensities) > 0
            count += 1

        assert count == 4, f"Expected 4 spectra, got {count}"
        reader.close()


class TestCrlfUnindentedImzML:
    """Unindented CRLF imzML (IONTOF SurfaceLab) must parse.

    pyimzml removes each <spectrum> from the tree while iterparse is still
    streaming. Under lxml that leaves libxml2's text-node accelerator pointing
    at a freed node, and on this one layout the parser eventually dies with
    ``XMLSyntaxError: xmlSAX2Characters`` -- an out-of-memory error dressed up
    as a syntax error, on a file that is perfectly well-formed.
    """

    def test_parser_is_not_constructed_with_lxml(
        self, create_minimal_imzml, monkeypatch
    ):
        """Pin the parser choice: lxml is unsafe for pyimzml's pruning loop.

        This asserts on the constructor argument rather than on behaviour
        because reproducing the corruption needs a ~16k-spectrum file, and the
        exact threshold shifts with heap layout. The behavioural test below
        covers it where it reproduces; this one holds the line everywhere else.
        """
        imzml_path, _, _, _ = create_minimal_imzml
        recorded = {}
        real_parser = imzml_reader_module.ImzMLParser

        def recording_parser(*args, **kwargs):
            recorded.update(kwargs)
            return real_parser(*args, **kwargs)

        monkeypatch.setattr(imzml_reader_module, "ImzMLParser", recording_parser)

        reader = ImzMLReader(imzml_path)
        reader.get_essential_metadata()
        reader.close()

        assert recorded.get("parse_lib") != "lxml"

    def test_reads_crlf_unindented_file(self, temp_dir):
        """The IONTOF layout round-trips through the reader."""
        imzml_path = write_crlf_unindented_imzml(temp_dir, 16)

        # Guard the fixture itself: if the writer or the reformatting changes
        # so the file is no longer unindented CRLF, this stops testing the
        # thing it claims to test.
        raw = imzml_path.read_bytes()
        assert b"\r\n<spectrum" in raw, "fixture is not unindented CRLF"

        reader = ImzMLReader(imzml_path)
        coordinates = [coords for coords, _, _ in reader.iter_spectra()]
        reader.close()

        assert len(coordinates) == 16
        # 16 spectra written across an 8-wide grid, 1-based in the file and
        # 0-based once the reader has converted them.
        assert coordinates[0] == (0, 0, 0)
        assert coordinates[-1] == (7, 1, 0)

    def test_reads_crlf_file_large_enough_to_corrupt_lxml(self, temp_dir):
        """The size that actually reproduces the libxml2 failure.

        Below roughly 16k spectra the stale-offset writes stay inside the
        buffer libxml2 already had, so the parse survives and proves nothing.
        The exact threshold moves with heap layout, so this may not reproduce
        on every platform -- but it can only ever under-report, never fail on
        correct code, since ElementTree does not go near libxml2. Deliberately
        NOT marked integration: CI runs -m "not integration", and this is the
        only test that exercises the actual defect. Costs ~2s and a ~33 MB
        temporary file.
        """
        imzml_path = write_crlf_unindented_imzml(temp_dir, 16000)

        reader = ImzMLReader(imzml_path)
        essential = reader.get_essential_metadata()
        reader.close()

        assert essential.n_spectra == 16000


def _write_two_pixel_imzml(directory: Path, name: str, spec_type: str) -> Path:
    """A minimal processed imzML declaring the given representation."""
    path = directory / f"{name}.imzML"
    mzs = np.linspace(100.0, 1000.0, 50)
    intensities = np.zeros_like(mzs)
    intensities[10] = 100.0
    with ImzMLWriter(str(path), mode="processed", spec_type=spec_type) as writer:
        writer.addSpectrum(mzs, intensities, (1, 1, 1))
        writer.addSpectrum(mzs, intensities, (2, 1, 1))
    return path


class TestSpectrumTypeReaderOption:
    """``reader_options={"spectrum_type": ...}`` has to reach the extractor.

    The extractor is where the decision is made, but it is constructed by the
    reader, so the value has to be threaded. This is the join that a unit test
    on the extractor alone would not catch.
    """

    def test_absent_means_detect(self, tmp_path):
        path = _write_two_pixel_imzml(tmp_path, "detect", "profile")
        reader = ImzMLReader(path)
        try:
            assert reader.spectrum_type is None
            assert reader.get_essential_metadata().spectrum_type == "profile spectrum"
        finally:
            reader.close()

    @pytest.mark.parametrize(
        "declared,override,expected",
        [
            ("profile", "centroid", "centroid spectrum"),
            ("centroid", "profile", "profile spectrum"),
        ],
    )
    def test_override_reaches_the_stored_metadata(
        self, tmp_path, declared, override, expected
    ):
        path = _write_two_pixel_imzml(tmp_path, f"ovr_{declared}", declared)
        reader = ImzMLReader(path, spectrum_type=override)
        try:
            assert reader.get_essential_metadata().spectrum_type == expected
        finally:
            reader.close()

    def test_reader_normalises_before_storing(self, tmp_path):
        path = _write_two_pixel_imzml(tmp_path, "norm", "profile")
        reader = ImzMLReader(path, spectrum_type="CENTROID")
        try:
            assert reader.spectrum_type == "centroid spectrum"
        finally:
            reader.close()

    def test_a_bad_value_fails_at_construction(self, tmp_path):
        """Fail while the caller is still looking at its own arguments."""
        path = _write_two_pixel_imzml(tmp_path, "bad", "profile")
        with pytest.raises(ValueError, match="Unknown spectrum_type"):
            ImzMLReader(path, spectrum_type="profil")
