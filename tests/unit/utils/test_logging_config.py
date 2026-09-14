import logging
import sys
from pathlib import Path

import pytest

from tests.conftest import restored_process_globals
from thyra.utils.logging_config import setup_logging


@pytest.fixture
def logger():
    # Ensure the logger is clean for each test. Nothing puts the handlers
    # back here on purpose: tests/conftest.py's autouse
    # _restore_process_globals does it for every test in the suite.
    logger = logging.getLogger("thyra")
    logger.handlers.clear()
    return logger


def test_setup_logging_console(logger):
    """Test that console logging is set up correctly."""
    setup_logging(log_level=logging.DEBUG)

    assert logger.level == logging.DEBUG
    assert len(logger.handlers) == 1
    assert isinstance(logger.handlers[0], logging.StreamHandler)
    assert logger.handlers[0].level == logging.DEBUG


def test_setup_logging_file(logger, tmp_path):
    """Test that file logging is set up correctly when a log file is provided."""
    log_file = tmp_path / "test.log"
    setup_logging(log_level=logging.INFO, log_file=str(log_file))

    assert len(logger.handlers) == 2

    file_handler = None
    for handler in logger.handlers:
        if isinstance(handler, logging.handlers.RotatingFileHandler):
            file_handler = handler
            break

    assert file_handler is not None
    assert file_handler.level == logging.INFO
    assert Path(file_handler.baseFilename).name == "test.log"

    # Test that a message is written to the log file
    test_message = "This is a test log message."
    logger.info(test_message)

    # Close the handler to ensure the message is flushed to the file
    file_handler.close()

    with open(log_file, "r") as f:
        log_content = f.read()
        assert test_message in log_content


def test_log_level_filtering(logger, tmp_path):
    """Test that log messages are filtered based on the configured log level."""
    log_file = tmp_path / "test_filtering.log"
    setup_logging(log_level=logging.WARNING, log_file=str(log_file))

    debug_message = "This is a debug message."
    warning_message = "This is a warning message."

    logger.debug(debug_message)
    logger.warning(warning_message)

    # Close file handler to ensure logs are written
    for handler in logger.handlers:
        if isinstance(handler, logging.handlers.RotatingFileHandler):
            handler.close()

    with open(log_file, "r") as f:
        log_content = f.read()
        assert debug_message not in log_content
        assert warning_message in log_content


class TestSetupLoggingDoesNotOutliveTheTest:
    """``setup_logging`` reconfigures the process, so something must undo it.

    Every one of these three is process-global and none of them is
    restored by the function that sets them: ``propagate = False`` on the
    ``thyra`` logger, the handler list it clears, and the level. Left
    behind, they decide what every later test can capture -- 23 assertions
    across 12 files went empty behind a single CLI test, and they went
    empty rather than failing loudly, which is the worst way for a test to
    be wrong. ``sys.argv`` is the same class of leak from the same tests.
    """

    def test_the_round_trip_puts_all_four_back(self, monkeypatch):
        """The guard proper: assert inside the block, then outside it.

        Both halves matter. The first says the poisoning is real, so the
        second is not asserting something that was already true; the
        second says the context manager undid it. The suite-wide autouse
        fixture wraps this test too, but it tears down after the body, so
        it cannot make the outer assertions pass on its own -- and neither
        can ``monkeypatch``, whose undo runs at teardown for the same
        reason. Setting argv through it rather than assigning is the
        convention ``tests/unit/test_log_capture_convention.py`` enforces.
        """
        thyra = logging.getLogger("thyra")
        before = (thyra.propagate, list(thyra.handlers), thyra.level)

        with restored_process_globals():
            setup_logging(log_level=logging.DEBUG)
            monkeypatch.setattr(sys, "argv", ["thyra", "in.imzML", "out.zarr"])

            assert thyra.propagate is False
            assert thyra.level == logging.DEBUG
            assert thyra.handlers != before[1]
            assert sys.argv[0] == "thyra"

        assert thyra.propagate is before[0]
        assert thyra.handlers == before[1]
        assert thyra.level == before[2]
        assert sys.argv[0] != "thyra"

    def test_a_log_file_handler_is_closed_so_windows_can_delete_it(self, tmp_path):
        """The one handler that must be closed, and the only one.

        A RotatingFileHandler holds the file open, and Windows refuses to
        delete a file that something has open -- which is why the fixture
        this replaced closed handlers at all. It closed every added
        handler though, and pytest attaches its own session-shared
        capture handlers to a non-propagating logger, so that also closed
        the plugin's.
        """
        log_file = tmp_path / "run.log"
        with restored_process_globals():
            setup_logging(log_level=logging.INFO, log_file=str(log_file))
            opened = [
                h
                for h in logging.getLogger("thyra").handlers
                if isinstance(h, logging.FileHandler)
            ]
            assert len(opened) == 1

        assert opened[0].stream is None
        log_file.unlink()

    def test_a_handler_pytest_owns_is_left_open(self):
        """Restoring must not close what it did not open.

        ``catching_logs`` attaches ``LoggingPlugin.caplog_handler`` and
        ``report_handler`` -- created once per session -- directly to any
        logger that is not propagating, which is exactly the state
        ``setup_logging`` leaves behind. Closing those closed them for the
        rest of the run.
        """
        thyra = logging.getLogger("thyra")
        borrowed = logging.StreamHandler()
        with restored_process_globals():
            setup_logging(log_level=logging.INFO)
            thyra.addHandler(borrowed)

        assert getattr(borrowed, "_closed", False) is False
