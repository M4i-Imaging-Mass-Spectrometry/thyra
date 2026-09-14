"""A failed SDK load must not be cached (issue #298).

``DLLManager.__new__`` published the new instance to ``cls._instance``
*before* calling ``_initialize``. So the first attempt raised its
``SDKError`` -- the long message naming the four places the SDK can go
-- and left a dead manager cached. Every later attempt in that process
returned that dead manager without raising anything, and the failure
resurfaced only at ``manager.dll`` as the bare "No Bruker SDK library
loaded", at a call site with no idea what went wrong.

Measured before the fix, in one process:

    attempt 1: SDKError, 740 chars, "Checked 0 locations."
    attempt 2: raised nothing; is_loaded False; same cached object
               m2.dll -> SDKError: 'No Bruker SDK library loaded'
    attempt 3 (force_reload=True): is_loaded True

So installing the SDK and retrying in the same session could not work,
and ``force_reload=True``, the escape hatch that would have fixed it,
has no caller anywhere in the tree.

``_instance`` is class state that never resets between tests, so every
test here isolates it.
"""

from __future__ import annotations

from unittest import mock

import pytest

from thyra.readers.bruker.timstof.sdk import dll_manager as dll_manager_module
from thyra.readers.bruker.timstof.sdk.dll_manager import DLLManager
from thyra.utils.bruker_exceptions import SDKError


@pytest.fixture(autouse=True)
def _isolate_the_singleton():
    """The cache is class state; do not inherit or leak one."""
    saved = (DLLManager._instance, DLLManager._dll, DLLManager._library_path)
    DLLManager._instance = None
    DLLManager._dll = None
    DLLManager._library_path = None
    try:
        yield
    finally:
        (
            DLLManager._instance,
            DLLManager._dll,
            DLLManager._library_path,
        ) = saved


def _discovery_finds_nothing():
    """Patch discovery so the load fails the way a missing SDK does.

    ``_load_library_by_name`` is the last resort after discovery turns
    up nothing; it must fail too, or the system PATH decides the test.
    """
    return mock.patch.multiple(
        dll_manager_module,
        get_dll_paths=mock.DEFAULT,
        validate_library_path=mock.DEFAULT,
    ), mock.patch.object(
        DLLManager, "_load_library_by_name", side_effect=OSError("no such library")
    )


def _sdk_is_present():
    """Discovery finds one path and loading it works, without a real DLL."""
    fake = mock.MagicMock(name="fake_cdll")
    path = mock.MagicMock(name="libtimsdata")
    return (
        mock.patch.object(dll_manager_module, "get_dll_paths", return_value=[path]),
        mock.patch.object(
            dll_manager_module, "validate_library_path", return_value=True
        ),
        mock.patch.object(DLLManager, "_load_library_at_path", return_value=fake),
        fake,
    )


class _Missing:
    """``with missing():`` -- discovery empty and the fallback failing."""

    def __enter__(self):
        self._multi, self._by_name = _discovery_finds_nothing()
        patched = self._multi.__enter__()
        patched["get_dll_paths"].return_value = []
        patched["validate_library_path"].return_value = False
        self._by_name.__enter__()
        return self

    def __exit__(self, *exc):
        self._by_name.__exit__(*exc)
        return self._multi.__exit__(*exc)


class _Present:
    """``with present() as fake:`` -- one discoverable, loadable library."""

    def __enter__(self):
        *self._patchers, self._fake = _sdk_is_present()
        self._entered = [p.__enter__() for p in self._patchers]
        return self._fake

    def __exit__(self, *exc):
        for p in reversed(self._patchers):
            p.__exit__(*exc)
        return False

    @property
    def load(self):
        return self._entered[-1]


class TestAFailedLoadIsNotCached:
    def test_the_first_attempt_raises(self):
        with _Missing():
            with pytest.raises(SDKError):
                DLLManager()

    def test_nothing_is_cached_after_it(self):
        with _Missing():
            with pytest.raises(SDKError):
                DLLManager()
        assert DLLManager._instance is None

    def test_the_second_attempt_raises_too(self):
        """It used to return the dead manager and raise nothing."""
        with _Missing():
            with pytest.raises(SDKError):
                DLLManager()
            with pytest.raises(SDKError) as excinfo:
                DLLManager()
        # And with the real message, not "No Bruker SDK library loaded".
        assert "BRUKER_SDK_PATH" in str(excinfo.value)

    def test_a_later_attempt_succeeds_once_the_sdk_is_there(self):
        """Installing the SDK and retrying in the same process works."""
        with _Missing():
            with pytest.raises(SDKError):
                DLLManager()

        with _Present() as fake:
            manager = DLLManager()

        assert manager.is_loaded is True
        assert manager.dll is fake


class TestASuccessfulLoadIsStillCached:
    def test_the_same_instance_comes_back(self):
        with _Present():
            first = DLLManager()
            second = DLLManager()

        assert first is second

    def test_the_library_is_loaded_once(self):
        present = _Present()
        with present:
            DLLManager()
            DLLManager()

        assert present.load.call_count == 1


class TestAStrandedSingletonRecovers:
    def test_a_failed_reload_does_not_strand_the_cache(self):
        """``reload()`` clears ``_dll`` before loading and mutates in
        place, so a local-then-publish cannot protect it. The
        ``is_loaded`` clause is what covers this route."""
        with _Present():
            manager = DLLManager()
        assert manager.is_loaded

        with _Missing():
            with pytest.raises(SDKError):
                manager.reload()
            assert manager.is_loaded is False
            # The next construction must not hand the stranded one back.
            with pytest.raises(SDKError):
                DLLManager()
