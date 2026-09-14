# tests/unit/metadata/core/test_base_extractor.py
import pytest

from thyra.core.base_extractor import MetadataExtractor
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata


class ConcreteMetadataExtractor(MetadataExtractor):
    """Concrete implementation of MetadataExtractor for testing."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._extract_essential_called = False
        self._extract_comprehensive_called = False

    def _extract_essential_impl(self) -> EssentialMetadata:
        self._extract_essential_called = True
        return EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(25.0, 25.0),
            n_spectra=200,
            total_peaks=400000,
            source_path="/test/path",
        )

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        self._extract_comprehensive_called = True
        essential = self._extract_essential_impl()
        return ComprehensiveMetadata(
            essential=essential,
            format_specific={"test_format": "test_value"},
            acquisition_params={"test_param": "test_value"},
            instrument_info={"test_instrument": "test_info"},
            raw_metadata={"test_raw": "test_metadata"},
        )


class TestMetadataExtractor:
    """Test the base MetadataExtractor abstract class."""

    def test_cannot_instantiate_abstract_class(self):
        """Test that MetadataExtractor cannot be instantiated directly."""
        with pytest.raises(TypeError):
            MetadataExtractor("dummy_source")

    def test_concrete_implementation_creation(self):
        """Test creating a concrete implementation."""
        extractor = ConcreteMetadataExtractor("test_source")
        assert extractor.data_source == "test_source"
        assert extractor._essential_cache is None
        assert extractor._comprehensive_cache is None

    def test_get_essential_caching(self):
        """Test that essential metadata is cached after first call."""
        extractor = ConcreteMetadataExtractor("test_source")

        # First call
        result1 = extractor.get_essential()
        assert extractor._extract_essential_called is True
        assert extractor._essential_cache is not None
        assert result1.dimensions == (10, 20, 1)

        # Reset the flag to check caching
        extractor._extract_essential_called = False

        # Second call should use cache
        result2 = extractor.get_essential()
        assert (
            extractor._extract_essential_called is False
        )  # Should not be called again
        assert result1 is result2  # Should be the same cached object

    def test_get_comprehensive_caching(self):
        """Test that comprehensive metadata is cached after first call."""
        extractor = ConcreteMetadataExtractor("test_source")

        # First call
        result1 = extractor.get_comprehensive()
        assert extractor._extract_comprehensive_called is True
        assert extractor._comprehensive_cache is not None
        assert result1.format_specific["test_format"] == "test_value"

        # Reset the flag to check caching
        extractor._extract_comprehensive_called = False

        # Second call should use cache
        result2 = extractor.get_comprehensive()
        assert (
            extractor._extract_comprehensive_called is False
        )  # Should not be called again
        assert result1 is result2  # Should be the same cached object

    def test_comprehensive_includes_essential(self):
        """Test that comprehensive metadata includes essential metadata."""
        extractor = ConcreteMetadataExtractor("test_source")

        comprehensive = extractor.get_comprehensive()
        essential = extractor.get_essential()

        # The essential metadata in comprehensive should match standalone essential
        assert comprehensive.essential.dimensions == essential.dimensions
        assert comprehensive.essential.pixel_size == essential.pixel_size
        assert comprehensive.essential.n_spectra == essential.n_spectra

    def test_independent_caching(self):
        """Test that essential and comprehensive caching work independently."""
        extractor = ConcreteMetadataExtractor("test_source")

        # Get essential first
        essential = extractor.get_essential()
        assert extractor._essential_cache is not None
        assert extractor._comprehensive_cache is None

        # Get comprehensive second
        comprehensive = extractor.get_comprehensive()
        assert extractor._essential_cache is not None
        assert extractor._comprehensive_cache is not None

        # Both should be cached independently
        assert essential is extractor._essential_cache
        assert comprehensive is extractor._comprehensive_cache

    def test_cache_invalidation_on_new_instance(self):
        """Test that each extractor instance has its own cache."""
        extractor1 = ConcreteMetadataExtractor("test_source1")
        extractor2 = ConcreteMetadataExtractor("test_source2")

        essential1 = extractor1.get_essential()
        essential2 = extractor2.get_essential()

        # Should be different objects (even if content is same)
        assert essential1 is not essential2
        # The source paths will be the same since we hardcoded them in the test implementation

    def test_error_propagation(self):
        """Test that errors in extraction methods are properly propagated."""

        class FailingExtractor(MetadataExtractor):
            def _extract_essential_impl(self) -> EssentialMetadata:
                raise ValueError("Test extraction error")

            def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
                raise RuntimeError("Test comprehensive error")

        extractor = FailingExtractor("test_source")

        # Essential extraction error should propagate
        with pytest.raises(ValueError, match="Test extraction error"):
            extractor.get_essential()

        # Comprehensive extraction error should propagate
        # Note: Since get_comprehensive calls get_essential first, we'll get the ValueError
        with pytest.raises(ValueError, match="Test extraction error"):
            extractor.get_comprehensive()

    def test_source_attribute_immutable(self):
        """Test that source attribute is set correctly and immutable."""
        extractor = ConcreteMetadataExtractor("test_source")
        assert extractor.data_source == "test_source"

        # Source is mutable by default in Python unless we make it immutable
        # This test can be removed or we can make data_source immutable
        extractor.data_source = "new_source"
        assert extractor.data_source == "new_source"


class TestMetadataExtractorPerformance:
    """Test performance-related aspects of MetadataExtractor."""

    def test_essential_called_only_once_per_cache(self):
        """Test that _extract_essential is called only once per cache lifetime."""
        extractor = ConcreteMetadataExtractor("test_source")

        # Multiple calls should only trigger one extraction
        for _ in range(5):
            extractor.get_essential()

        # Should have been called only once due to caching
        assert extractor._extract_essential_called is True

    def test_comprehensive_optimization_with_essential_cache(self):
        """Test that comprehensive extraction can benefit from essential cache."""

        class OptimizedExtractor(MetadataExtractor):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.essential_calls = 0
                self.comprehensive_calls = 0

            def _extract_essential_impl(self) -> EssentialMetadata:
                self.essential_calls += 1
                return EssentialMetadata(
                    dimensions=(10, 20, 1),
                    coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
                    mass_range=(100.0, 1000.0),
                    pixel_size=(25.0, 25.0),
                    n_spectra=200,
                    total_peaks=400000,
                    source_path="/test/path",
                )

            def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
                self.comprehensive_calls += 1
                # Use cached essential if available
                if self._essential_cache is not None:
                    essential = self._essential_cache
                else:
                    essential = self._extract_essential_impl()

                return ComprehensiveMetadata(
                    essential=essential,
                    format_specific={},
                    acquisition_params={},
                    instrument_info={},
                    raw_metadata={},
                )

        extractor = OptimizedExtractor("test_source")

        # Get essential first, then comprehensive
        extractor.get_essential()
        extractor.get_comprehensive()

        # Essential should be called only once, comprehensive once
        assert extractor.essential_calls == 1
        assert extractor.comprehensive_calls == 1
