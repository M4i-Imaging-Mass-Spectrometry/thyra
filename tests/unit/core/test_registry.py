"""Unit tests for MSIRegistry format detection and class registration."""

import threading
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from thyra.core.base_converter import BaseMSIConverter
from thyra.core.base_reader import BaseMSIReader
from thyra.core.registry import (
    MSIRegistry,
    _registry,
    detect_format,
    get_converter_class,
    get_reader_class,
    register_converter,
    register_reader,
)

# ---------------------------------------------------------------------------
# Helpers: minimal concrete subclasses used across tests
# ---------------------------------------------------------------------------


class _StubReader(BaseMSIReader):
    def get_metadata(self):
        pass

    def get_dimensions(self):
        pass

    def get_common_mass_axis(self):
        pass

    def iter_spectra(self, batch_size=None):
        pass

    def close(self):
        pass


class _StubConverter(BaseMSIConverter):
    def _create_data_structures(self):
        pass

    def _save_output(self, data_structures):
        pass


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _isolated_registry():
    """Snapshot the global registry before each test and restore it after."""
    with _registry._lock:
        saved_readers = _registry._readers.copy()
        saved_converters = _registry._converters.copy()

    yield

    with _registry._lock:
        _registry._readers.clear()
        _registry._readers.update(saved_readers)
        _registry._converters.clear()
        _registry._converters.update(saved_converters)


# ---------------------------------------------------------------------------
# detect_format
# ---------------------------------------------------------------------------


def test_detect_format_raises_for_nonexistent_path(tmp_path):
    """detect_format raises ValueError when path does not exist."""
    missing = tmp_path / "does_not_exist.imzml"
    with pytest.raises(ValueError, match="Input path does not exist"):
        detect_format(missing)


def test_detect_format_imzml_mock():
    """detect_format returns 'imzml' for .imzml path (using mock)."""
    mock_path = MagicMock(spec=Path)
    mock_path.exists.return_value = True
    mock_path.suffix = ".imzml"
    mock_path.is_dir.return_value = False
    # Paired .ibd file must exist
    ibd_mock = MagicMock(spec=Path)
    ibd_mock.exists.return_value = True
    mock_path.with_suffix.return_value = ibd_mock

    with patch.object(
        _registry.__class__, "detect_format", wraps=_registry.detect_format
    ):
        registry = MSIRegistry()
        result = registry.detect_format(mock_path)

    assert result == "imzml"


def test_detect_format_raises_for_unsupported_extension(tmp_path):
    """detect_format raises ValueError for an unknown file extension."""
    unknown = tmp_path / "data.xyz"
    unknown.touch()
    with pytest.raises(ValueError, match="Unsupported format"):
        detect_format(unknown)


# ---------------------------------------------------------------------------
# register_reader / get_reader_class round-trip
# ---------------------------------------------------------------------------


def test_register_and_get_reader_class():
    """register_reader followed by get_reader_class returns the registered class."""
    register_reader("stub_fmt")(_StubReader)
    assert get_reader_class("stub_fmt") is _StubReader


def test_get_reader_class_raises_for_unknown_format():
    """get_reader_class raises ValueError for an unregistered format name."""
    with pytest.raises(ValueError, match="No reader for format"):
        get_reader_class("totally_unknown_format_xyz")


# ---------------------------------------------------------------------------
# register_converter / get_converter_class round-trip
# ---------------------------------------------------------------------------


def test_register_and_get_converter_class():
    """register_converter followed by get_converter_class returns the registered class."""
    register_converter("stub_fmt")(_StubConverter)
    assert get_converter_class("stub_fmt") is _StubConverter


def test_get_converter_class_raises_for_unknown_format():
    """get_converter_class raises ValueError for an unregistered format name."""
    with pytest.raises(ValueError, match="No converter for format"):
        get_converter_class("totally_unknown_format_xyz")


# ---------------------------------------------------------------------------
# Thread safety
# ---------------------------------------------------------------------------


def test_concurrent_registration_is_thread_safe():
    """Registering readers from multiple threads produces consistent results."""
    n_threads = 20
    errors = []

    def register_and_verify(fmt_name, reader_cls):
        try:
            register_reader(fmt_name)(reader_cls)
            retrieved = get_reader_class(fmt_name)
            assert retrieved is reader_cls
        except Exception as exc:
            errors.append(exc)

    # Create unique reader subclasses per thread to avoid cross-contamination
    classes = []
    for i in range(n_threads):
        cls = type(f"_ThreadReader{i}", (_StubReader,), {})
        classes.append(cls)

    threads = [
        threading.Thread(
            target=register_and_verify,
            args=(f"thread_fmt_{i}", classes[i]),
        )
        for i in range(n_threads)
    ]

    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert errors == [], f"Thread safety errors: {errors}"

    # Verify all registrations persisted
    for i in range(n_threads):
        assert get_reader_class(f"thread_fmt_{i}") is classes[i]
