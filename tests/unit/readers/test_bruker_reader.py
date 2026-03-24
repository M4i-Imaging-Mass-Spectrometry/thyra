"""
Tests for the Bruker reader.
Note: Full testing requires actual Bruker files and the timsdata DLL.
These tests focus on structure and interface validation.
"""

from thyra.readers.bruker.timstof.timstof_reader import BrukerReader


class TestBrukerReaderStructure:
    """Test the structure and interface of the Bruker reader."""

    def test_class_registration(self):
        """Test that the BrukerReader class is properly registered."""
        from thyra.core.registry import get_reader_class

        # Test that we can get the bruker reader class
        reader_class = get_reader_class("bruker")
        assert reader_class == BrukerReader

    def test_interface_implementation(self):
        """Test that the BrukerReader implements the BaseMSIReader interface."""
        from thyra.core.base_reader import BaseMSIReader

        assert issubclass(BrukerReader, BaseMSIReader)

        # Check that it implements all required methods
        required_methods = [
            "get_essential_metadata",
            "get_comprehensive_metadata",
            "get_common_mass_axis",
            "iter_spectra",
            "close",
            "_create_metadata_extractor",
        ]

        for method in required_methods:
            assert hasattr(BrukerReader, method)
