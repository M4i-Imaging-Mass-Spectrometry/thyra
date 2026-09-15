"""
Common test fixtures for thyra tests.
"""

import logging
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, List, Union

import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

# Root directory of the tests
TEST_DIR = Path(__file__).parent.resolve()
# Test data directory
DATA_DIR = TEST_DIR / "data"

# The two lanes, as directories. Membership is decided from these and from
# nothing else -- see pytest_collection_modifyitems below.
UNIT_DIR = TEST_DIR / "unit"
INTEGRATION_DIR = TEST_DIR / "integration"


def pytest_configure(config):
    """Give the suite the same Windows rename retry a conversion gets.

    Zarr renames each metadata document and each rewritten chunk over its
    destination (``zarr/storage/_local.py::_atomic_write``), and on Windows
    that rename has no atomic form: ``MoveFileEx`` must delete the
    destination first, and for well under a millisecond after a name was
    last used NTFS still holds it, so the delete comes back
    ``[WinError 5] Access is denied``. Measured on this project's Windows 11
    machine with pure stdlib and nothing else imported: a replace onto an
    EXISTING file fails 720 times in 20,000 (3.6%); onto an absent one,
    never (0 in 3,000); and onto a FRESH name, never either, which is what
    rules the on-access scanner out. See ``thyra.utils.zarr_atomic_write``
    for the full measurement.

    ``thyra.utils.zarr_atomic_write`` already fixes this for anything a
    converter writes, and :class:`StreamedOpticalImage` now installs it
    too. What is left uncovered is the suite's own raw
    ``SpatialData(...).write(...)`` calls -- the reference stores a test
    builds to compare Thyra's output against. Those are plain spatialdata,
    so no Thyra entry point runs first, and they flake:
    ``tests/unit/converters/test_optical_image_streaming.py`` failed 4 runs
    out of 5 before this hook and 0 out of 8 after it.

    This installs the retry, it does not skip anything: every assertion
    still runs on Windows. CI never needed it -- 0 ``test (windows-latest)``
    failures in 60 runs, the Actions Windows images ship with real-time
    monitoring off -- this is for local Windows development.
    """
    from thyra.utils.zarr_atomic_write import install_windows_atomic_write_retry

    install_windows_atomic_write_retry()


def pytest_collection_modifyitems(config, items):
    """Stamp each test's lane marker from the directory it lives in.

    ``-m integration`` used to mean "the files somebody remembered to
    decorate", which is not the same set as ``tests/integration/``: 14 of
    that directory's 65 tests carried no marker, so the lane CI presents as
    the integration lane silently skipped every end-to-end CLI test. The
    inverse held too -- 14 tests under ``tests/unit/`` claimed the
    integration lane by decorator.

    Stamping by location makes the two agree by construction, and makes a
    misfiled test the only way to get it wrong -- which
    ``tests/unit/test_lane_markers.py`` then catches. The hook only adds and
    never removes, and it does not look at what an item already carries: a
    hand-written ``@pytest.mark.integration`` under ``tests/unit/`` would end
    up with both markers rather than with the one it asked for. There are zero
    such decorators today and ``test_lane_markers.py`` is what keeps it that
    way, so the case is a hazard rather than a behaviour. It runs before ``-m``
    deselection, so the stamp is honoured by the same run that applies it.

    Args:
        config: The pytest config (unused; part of the hook signature).
        items: The collected items, modified in place.
    """
    for item in items:
        raw = getattr(item, "path", None)
        if raw is None:
            continue
        path = Path(raw).resolve()
        if path.is_relative_to(INTEGRATION_DIR):
            item.add_marker(pytest.mark.integration)
        elif path.is_relative_to(UNIT_DIR):
            item.add_marker(pytest.mark.unit)


@contextmanager
def restored_process_globals() -> Iterator[None]:
    """Undo what invoking the CLI does to the process.

    ``setup_logging`` (thyra/utils/logging_config.py:45-53) sets
    ``propagate = False`` on the ``thyra`` logger, clears its handlers and
    sets its level, and ``thyra.__main__.main`` calls it on every
    invocation. All three are process-global, so without this one CLI test
    decides what every test after it can capture -- and ``sys.argv``, which
    ``tests/integration/test_cli.py`` used to assign and never restore, is
    the same class of leak. Restoring here means no test has to remember
    to, and which pytest version is installed stops mattering: 9.x hides
    the symptom by walking ``loggerDict`` for non-propagating loggers,
    8.x does not.

    Snapshotting the root logger as well is insurance rather than a fix
    for anything measured: ``logging.basicConfig`` is the usual way to
    leak into it, and under pytest it is a no-op, because pytest has
    already attached its own handlers to root and CPython only applies
    ``basicConfig``'s ``level=`` when root has none.

    A context manager and not only a fixture, so that
    ``tests/unit/utils/test_logging_config.py`` can assert the round trip
    inside one test rather than across two whose order would decide the
    answer.
    """
    argv_object = sys.argv
    argv_contents = list(sys.argv)
    snapshots = [
        (logger, list(logger.handlers), logger.propagate, logger.level)
        for logger in (logging.getLogger("thyra"), logging.getLogger())
    ]
    try:
        yield
    finally:
        # Both spellings: a module that did ``from sys import argv`` holds
        # the object, everything else reads the attribute.
        sys.argv = argv_object
        sys.argv[:] = argv_contents
        for logger, handlers, propagate, level in snapshots:
            added = [h for h in logger.handlers if h not in handlers]
            logger.handlers = handlers
            logger.propagate = propagate
            logger.setLevel(level)
            for handler in added:
                # Close only what owns an OS handle. A RotatingFileHandler
                # opened by --log-file keeps a tmp_path file open, and
                # Windows will not delete a file that is open. Everything
                # else added during the test may be pytest's own: while a
                # logger is non-propagating, catching_logs attaches the
                # plugin's session-shared caplog_handler and report_handler
                # to it directly, and closing those closes them for the
                # rest of the session.
                if isinstance(handler, logging.FileHandler):
                    handler.close()


@pytest.fixture(autouse=True)
def _restore_process_globals() -> Iterator[None]:
    """Wrap every test in :func:`restored_process_globals`.

    Autouse and declared here rather than per-file, because the tests that
    need it are not the tests that cause it: the CLI test poisons the
    process and some unrelated reader test is what fails.
    """
    with restored_process_globals():
        yield


class _Records(List[logging.LogRecord]):
    """The records ``thyra_logs`` collected, shaped the way ``caplog`` is.

    ``messages`` and ``text`` exist so that moving an assertion off
    ``caplog`` is a rename and not a rewrite, and ``text`` is formatted
    the way caplog formats its own so an assertion carried across still
    matches the same substring.
    """

    @property
    def messages(self) -> List[str]:
        return [record.getMessage() for record in self]

    @property
    def text(self) -> str:
        return "\n".join(
            f"{record.levelname:<8} {record.name}:{record.filename}:"
            f"{record.lineno} {record.getMessage()}"
            for record in self
        )


@pytest.fixture
def thyra_logs():
    """A context manager collecting records from a named Thyra logger.

    Deliberately not ``caplog``. ``setup_logging`` sets
    ``propagate = False`` on the ``thyra`` logger and that is
    process-global: once any test in the session has invoked the CLI,
    caplog's root handler never sees another Thyra record. A caplog
    assertion on Thyra's own logging therefore passes alone and fails in
    the full suite, which is the worst way for a test to be wrong.
    Attaching a handler to the named logger sidesteps propagation, so the
    result does not depend on which tests ran first.

    ``_restore_process_globals`` above now puts the logger back after
    every test, which closes the same hole from the other side. Both are
    wanted: the fixture keeps one test from poisoning the next, this keeps
    a test that invokes the CLI *inside* its own capture block working --
    ``setup_logging`` clears ``thyra``'s handlers, so name a child logger
    in that case, which it does not touch.

    Usage::

        with thyra_logs("thyra.cli", logging.WARNING) as records:
            ...
        assert [r.getMessage() for r in records] == [...]
    """

    @contextmanager
    def capture(
        logger_name: str = "thyra", level: Union[int, str] = logging.INFO
    ) -> Iterator[_Records]:
        # ``caplog.at_level`` takes "WARNING" as happily as logging.WARNING
        # and call sites migrated off it did. Without this the level
        # comparison below raises TypeError on a str.
        if isinstance(level, str):
            level = logging.getLevelNamesMapping()[level]
        collected = _Records()

        class _Collector(logging.Handler):
            def emit(self, record: logging.LogRecord) -> None:
                # What logging.Formatter.format sets, and what caplog's
                # handler therefore leaves behind, because it formats on
                # emit. Assertions reading record.message work without
                # this only while pytest's own root handler happens to
                # have formatted the same record first -- which it does
                # not do for a logger it never sees.
                record.message = record.getMessage()
                collected.append(record)

        logger = logging.getLogger(logger_name)
        handler = _Collector(level=level)
        previous = logger.level
        logger.addHandler(handler)
        if previous > level or previous == logging.NOTSET:
            logger.setLevel(level)
        try:
            yield collected
        finally:
            logger.removeHandler(handler)
            logger.setLevel(previous)

    return capture


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def mock_reader():
    """Create a mock MSI reader for testing converters."""
    from thyra.core.base_reader import BaseMSIReader

    class MockMSIReader(BaseMSIReader):
        def __init__(self, data_path=None, **kwargs):
            super().__init__(data_path or Path("/mock/path"), **kwargs)
            self.closed = False

        def _create_metadata_extractor(self):
            # Create a mock metadata extractor
            from thyra.core.base_extractor import MetadataExtractor
            from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata

            class MockExtractor(MetadataExtractor):
                def _extract_essential_impl(self):
                    return EssentialMetadata(
                        dimensions=(3, 3, 1),
                        coordinate_bounds=(0.0, 2.0, 0.0, 2.0),
                        mass_range=(100.0, 1000.0),
                        pixel_size=None,
                        n_spectra=9,
                        total_peaks=900,  # 9 spectra * 100 peaks each
                        source_path="/mock/path",
                    )

                def _extract_comprehensive_impl(self):
                    return ComprehensiveMetadata(
                        essential=self._extract_essential_impl(),
                        format_specific={"format": "mock"},
                        acquisition_params={},
                        instrument_info={"instrument": "test_instrument"},
                        raw_metadata={"source": "mock"},
                    )

            return MockExtractor(None)

        def get_common_mass_axis(self):
            return np.linspace(100, 1000, 100)  # 100 mass values

        def iter_spectra(self):
            mass_axis = self.get_common_mass_axis()
            for x in range(3):
                for y in range(3):
                    # Create simple synthetic spectrum
                    intensities = np.zeros_like(mass_axis)
                    # Add a few peaks
                    intensities[x * 10 + 20] = 100.0  # Peak varying by x position
                    intensities[y * 10 + 50] = 200.0  # Peak varying by y position
                    yield ((x, y, 0), mass_axis, intensities)

        def close(self):
            self.closed = True

    return MockMSIReader()


@pytest.fixture
def create_minimal_imzml(temp_dir):
    """
    Create a minimal imzML file for testing.
    Returns a tuple of (imzml_path, ibd_path, mzs, intensities)
    """
    imzml_path = temp_dir / "minimal.imzML"
    ibd_path = temp_dir / "minimal.ibd"

    # Create small sample data
    coordinates = [(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)]  # 2x2 grid
    mzs = np.linspace(100, 1000, 50)  # 50 m/z values

    # Create different intensities for each pixel
    all_intensities = []
    for i, (x, y, z) in enumerate(coordinates):
        intensities = np.zeros_like(mzs)
        # Create a few peaks with position-dependent intensity
        intensities[10] = 100.0 * x  # Peak intensity depends on x
        intensities[30] = 150.0 * y  # Peak intensity depends on y
        all_intensities.append(intensities)

    # Write imzML file
    with ImzMLWriter(str(imzml_path), mode="processed") as writer:
        for i, (x, y, z) in enumerate(coordinates):
            writer.addSpectrum(mzs, all_intensities[i], (x, y, z))

    return imzml_path, ibd_path, mzs, all_intensities


def _write_square_imzml(path, origin, side):
    """Write a ``side`` x ``side`` imzML whose lowest coordinate is ``origin``.

    Every pixel carries one peak whose intensity encodes its position, so a
    row that moved can be told from a row that went missing.

    Args:
        path: Where to write the ``.imzML`` (the ``.ibd`` goes beside it).
        origin: The smallest x and y written to the file. ``1`` is what the
            specification says; ``0`` is what some exporters write.
        side: Pixels per side.

    Returns:
        ``(path, mzs, expected)`` where ``expected`` maps the 0-based
        ``(x, y)`` the store should hold to that pixel's intensity.
    """
    mzs = np.linspace(100.0, 200.0, 5)
    expected = {}

    with ImzMLWriter(str(path), mode="processed") as writer:
        for row in range(side):
            for col in range(side):
                intensity = np.zeros_like(mzs)
                # Distinct per pixel, and never zero: a pixel whose only
                # peak is zero gets no row at all, which would hide the
                # very loss this fixture exists to detect.
                intensity[2] = 100.0 + 10.0 * row + col
                writer.addSpectrum(mzs, intensity, (origin + col, origin + row, 1))
                expected[(col, row)] = intensity[2]

    return path, mzs, expected


@pytest.fixture
def zero_based_imzml(temp_dir):
    """A 3x3 imzML numbered from 0 -- the shape of issue #244.

    The imzML specification numbers pixels from 1, and the reader used to
    subtract a constant 1. On a file written 0-based that produced
    ``x = -1`` for the first column, which the converter's grid guard
    dropped: the store came out 4 rows of a 3x3 acquisition, with a
    warning naming a 2x2 grid the file never declared, and exit 0.
    """
    return _write_square_imzml(temp_dir / "zero_based.imzML", origin=0, side=3)


@pytest.fixture
def one_based_imzml(temp_dir):
    """The same 3x3 acquisition written the way the specification says."""
    return _write_square_imzml(temp_dir / "one_based.imzML", origin=1, side=3)


@pytest.fixture
def cropped_one_based_imzml(temp_dir):
    """A 1-based acquisition whose leftmost column is x = 5, not x = 1.

    The reason x and y do not simply rebase on their observed minimum the
    way z does: this file is a legitimately cropped region of the slide,
    and sliding it to the origin would change its grid width, every
    ``obs`` coordinate and its pixel footprint -- on a file that converts
    correctly today.
    """
    return _write_square_imzml(temp_dir / "cropped.imzML", origin=5, side=3)


@pytest.fixture
def mock_waters_data(temp_dir):
    """Create a mock Waters .raw directory structure.

    Creates the minimal structure needed for WatersReader validation:
    a directory with a _FUNC001.DAT file.
    """
    raw_dir = temp_dir / "mock.raw"
    raw_dir.mkdir()
    (raw_dir / "_FUNC001.DAT").write_bytes(b"\x00" * 64)
    return raw_dir
