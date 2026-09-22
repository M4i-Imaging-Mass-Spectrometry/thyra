"""
Tests for the simplified format registry system.
"""

import ast
import subprocess
import sys
import threading
from pathlib import Path
from typing import Dict, Optional

import pytest

import thyra
from thyra.core import registry as registry_module
from thyra.core.base_converter import BaseMSIConverter
from thyra.core.base_reader import BaseMSIReader
from thyra.core.registry import (
    _get_bruker_folder_structure,
    _registry,
    detect_format,
    get_converter_class,
    get_reader_class,
    register_converter,
    register_reader,
)
from thyra.errors import ConversionRefused


class TestRegistry:
    """Test the registry functionality."""

    def setup_method(self):
        """Set up each test by clearing the registries."""
        # Store original registry values to restore later
        with _registry._lock:
            self.original_readers = _registry._readers.copy()
            self.original_converters = _registry._converters.copy()

            # Clear registries for testing
            _registry._readers.clear()
            _registry._converters.clear()

    def teardown_method(self):
        """Restore original registry values after each test."""
        with _registry._lock:
            _registry._readers.clear()
            _registry._readers.update(self.original_readers)

            _registry._converters.clear()
            _registry._converters.update(self.original_converters)

    def test_register_reader(self):
        """Test registering a reader class."""

        # Create a test reader class
        class TestReader(BaseMSIReader):
            def get_metadata(self):
                pass

            def get_dimensions(self):
                pass

            def get_common_mass_axis(self):
                pass

            def iter_spectra(self):
                pass

            def close(self):
                pass

        # Register it
        register_reader("test_format")(TestReader)

        # Check if it was properly registered
        with _registry._lock:
            assert "test_format" in _registry._readers
            assert _registry._readers["test_format"] == TestReader

        # Test getting the reader class
        retrieved_class = get_reader_class("test_format")
        assert retrieved_class == TestReader

    def test_register_converter(self):
        """Test registering a converter class."""

        # Create a test converter class
        class TestConverter(BaseMSIConverter):
            def _create_data_structures(self):
                pass

            def _save_output(self, data_structures):
                pass

        # Register it
        register_converter("test_format")(TestConverter)

        # Check if it was properly registered
        with _registry._lock:
            assert "test_format" in _registry._converters
            assert _registry._converters["test_format"] == TestConverter

        # Test getting the converter class
        retrieved_class = get_converter_class("test_format")
        assert retrieved_class == TestConverter

    def test_detect_format_imzml(self, tmp_path):
        """Test ImzML format detection via extension."""
        # Create test files
        imzml_file = tmp_path / "test.imzml"
        ibd_file = tmp_path / "test.ibd"
        imzml_file.touch()
        ibd_file.touch()

        assert detect_format(imzml_file) == "imzml"

    def test_detect_format_bruker(self, tmp_path):
        """Test Bruker format detection via extension."""
        # Create test directory
        bruker_dir = tmp_path / "test.d"
        bruker_dir.mkdir()
        (bruker_dir / "analysis.tsf").touch()

        assert detect_format(bruker_dir) == "bruker"

    def test_detect_format_bruker_tdf(self, tmp_path):
        """Test Bruker format detection with .tdf file."""
        # Create test directory
        bruker_dir = tmp_path / "test.d"
        bruker_dir.mkdir()
        (bruker_dir / "analysis.tdf").touch()

        assert detect_format(bruker_dir) == "bruker"

    def test_unsupported_extension(self, tmp_path):
        """Test error for unsupported extension."""
        unknown_file = tmp_path / "test.xyz"
        unknown_file.touch()

        with pytest.raises(ValueError, match="Unsupported format"):
            detect_format(unknown_file)

    def test_shimadzu_imdx_recognised_but_in_development(self, tmp_path):
        """Shimadzu .imdx files get a guidance error, not 'unsupported'."""
        imdx_file = tmp_path / "test.imdx"
        imdx_file.touch()

        with pytest.raises(ValueError, match="in development"):
            detect_format(imdx_file)

    def test_shimadzu_kbd_recognised_but_in_development(self, tmp_path):
        """Shimadzu .kbd files get a guidance error, not 'unsupported'."""
        kbd_file = tmp_path / "test.kbd"
        kbd_file.touch()

        with pytest.raises(ValueError, match="export the dataset as imzML"):
            detect_format(kbd_file)

    def test_missing_ibd_file(self, tmp_path):
        """Test error for ImzML without .ibd file."""
        imzml_file = tmp_path / "test.imzml"
        imzml_file.touch()

        with pytest.raises(ValueError, match="requires corresponding .ibd file"):
            detect_format(imzml_file)

    def test_bruker_missing_analysis_files(self, tmp_path):
        """Test error for Bruker .d directory without analysis files."""
        bruker_dir = tmp_path / "test.d"
        bruker_dir.mkdir()

        with pytest.raises(ValueError, match="missing analysis files"):
            detect_format(bruker_dir)

    def test_bruker_not_directory(self, tmp_path):
        """Test error for .d file instead of directory."""
        fake_bruker = tmp_path / "test.d"
        fake_bruker.touch()  # Create as file, not directory

        with pytest.raises(ValueError, match="requires .d directory"):
            detect_format(fake_bruker)

    def test_nonexistent_path(self, tmp_path):
        """Test error for non-existent path."""
        nonexistent = tmp_path / "nonexistent.imzml"

        with pytest.raises(ValueError, match="Input path does not exist"):
            detect_format(nonexistent)

    def test_detect_format_waters(self, tmp_path):
        """Test Waters format detection via .raw directory with _FUNC*.DAT."""
        raw_dir = tmp_path / "test.raw"
        raw_dir.mkdir()
        (raw_dir / "_FUNC001.DAT").write_bytes(b"\x00" * 16)

        assert detect_format(raw_dir) == "waters"

    def test_raw_file_without_phi_magic_is_rejected(self, tmp_path):
        """A .raw file is a PHI candidate, so the error names both vendors.

        Waters .raw is a directory and PHI .raw is a file, so a plain file
        that lacks the PHI SOFH magic belongs to neither.
        """
        fake_raw = tmp_path / "test.raw"
        fake_raw.touch()

        with pytest.raises(ValueError, match="Unrecognised .raw file"):
            detect_format(fake_raw)

    def test_detect_format_phi(self, tmp_path):
        """Test PHI format detection via .raw file with SOFH magic."""
        phi_raw = tmp_path / "test.raw"
        phi_raw.write_bytes(b"SOFH\r\nImagePixels: 4\r\nEOFH\r\n")

        assert detect_format(phi_raw) == "phi"

    def test_waters_missing_func_files(self, tmp_path):
        """Test error for .raw directory without _FUNC*.DAT files."""
        raw_dir = tmp_path / "test.raw"
        raw_dir.mkdir()

        with pytest.raises(ValueError, match="_FUNC.*DAT"):
            detect_format(raw_dir)

    def test_detect_format_waters_generic_directory(self, tmp_path):
        """Test Waters detection from a generic directory (no .raw extension)."""
        some_dir = tmp_path / "my_data"
        some_dir.mkdir()
        (some_dir / "_FUNC001.DAT").write_bytes(b"\x00" * 16)

        assert detect_format(some_dir) == "waters"

    def test_get_nonexistent_reader(self):
        """A reader miss is a refusal, not a bare ValueError.

        Asserted as ``ConversionRefused`` rather than ``ValueError``
        because the type is what ``convert_msi`` dispatches on: the
        refusal handler prints the message once and keeps the traceback
        for DEBUG, the generic handler prints a traceback at ERROR. A
        ``ValueError`` assertion passes under either, so it could not tell
        the two apart -- and the registry is the one place every caller of
        a format name goes through, which is why the convention belongs
        here rather than in each caller's error mapping.
        """
        with pytest.raises(ConversionRefused, match="No reader for format"):
            get_reader_class("nonexistent_format")

    def test_get_nonexistent_converter(self):
        """A converter miss is a refusal too, for the same reason.

        ``thyra/convert.py`` used to route this lookup through a wrapper,
        ``_resolve_converter_class``, that caught the registry's error and
        re-raised ``ConversionRefused`` for any format name containing
        "spatialdata". The wrapper is gone and ``_create_converter`` calls
        ``get_converter_class`` directly.

        Two things hold the convention on that path, not one. This
        assertion holds it at the registry, where every caller of a format
        name arrives. ``TestAnUnregisteredOutputFormat`` in
        ``tests/unit/test_cli_refusals.py``, added by the same commit,
        holds the other end: it drives the lookup through ``convert_msi``
        and asserts the registry's own message reaches the user as a
        single ERROR line naming the formats there are, with the traceback
        kept for DEBUG. Either can fail while the other passes -- this one
        if the registry stops refusing, that one if something between the
        registry and the user starts rewriting the refusal again -- which
        is why both are here.
        """
        with pytest.raises(ConversionRefused, match="No converter for format"):
            get_converter_class("nonexistent_format")


# The two below live outside TestRegistry on purpose: its setup_method
# clears _registry._converters and teardown_method restores it, so the same
# assertions inside the class would only ever test that fixture.
#
# Both were measured to pass against the pre-#310 code as well, and that is
# not a flaw in them: on an install where spatialdata imports, the flag they
# used to be gated behind was True and the behaviour was already correct.
# They state the invariant so a future reviewer can see it asserted
# somewhere; the test that actually separates the two trees is
# tests/unit/test_hard_dependency.py, which fails before the fix.


def test_spatialdata_converter_is_registered_on_import():
    """Importing thyra registers the one output format the docs describe.

    Registration used to sit behind ``if SPATIALDATA_AVAILABLE:``, a flag
    set by a try/except that swallowed the ImportError. When spatialdata
    could not be imported, ``import thyra`` still succeeded and this lookup
    raised "No converter for format 'spatialdata'. Available: []" -- an
    empty registry, naming neither the missing package nor the cause
    (issue #310). spatialdata is a hard dependency, so registration is now
    unconditional and a broken install fails at ``import thyra`` instead.
    """
    assert get_converter_class("spatialdata") is thyra.SpatialDataConverter


class TestTheLazyImportHoldsNoState:
    """Issue #284: the lock's blind spot is deleted, not locked.

    ``_get_bruker_folder_structure`` used to memoise into the module global
    ``_bruker_folder_structure_module`` under an unsynchronised
    check-then-set -- the only mutable global state in ``core.registry``,
    and the only part ``MSIRegistry._lock`` did not cover. The memo is gone:
    ``sys.modules`` already caches the import, and the alternative the issue
    floated -- extending the registry lock over it -- is a deadlock, not a
    fix. These pin both halves of that.
    """

    def test_the_memoising_global_is_not_back(self):
        assert not hasattr(registry_module, "_bruker_folder_structure_module")

    def test_calling_it_adds_no_module_state(self):
        """Any reintroduced memo would show up as a new module attribute."""
        before = set(vars(registry_module))
        _get_bruker_folder_structure()
        assert set(vars(registry_module)) == before

    def test_every_call_yields_the_same_two_objects(self):
        """``sys.modules`` is the cache, so identity holds without a memo."""
        first = _get_bruker_folder_structure()
        second = _get_bruker_folder_structure()
        assert first[0] is second[0]
        assert first[1] is second[1]

    def test_the_lazy_import_does_not_wait_on_the_registry_lock(self):
        """The shape that deadlocks, encoded so it fails instead of hanging.

        Holding ``_registry._lock`` across the import means a caller can own
        the lock while waiting for the ``thyra.readers`` import lock, which a
        thread part-way through importing that package holds while waiting
        for the registry lock inside a ``@register_reader``. Measured on
        CPython 3.13.3: both threads hang for good.

        Here the main thread plays the registering thread and simply holds
        the lock. If the lazy import is ever put under that lock, the worker
        cannot finish and the wait expires -- a failure, not a hung run,
        which is why the thread is a daemon and the wait has a timeout.
        """
        finished = threading.Event()

        def call_it() -> None:
            _get_bruker_folder_structure()
            finished.set()

        with _registry._lock:
            worker = threading.Thread(target=call_it, daemon=True)
            worker.start()
            assert finished.wait(timeout=10), (
                "_get_bruker_folder_structure blocked while another thread "
                "held the registry lock: it is being taken under that lock, "
                "which is the deadlock shape issue #284 measured"
            )


def test_spatialdata_converter_is_always_a_class():
    """``thyra.SpatialDataConverter`` is a class, never None.

    ``thyra/__init__.py`` carried an ``except ImportError`` branch rebinding
    this name to ``None``, and issue #282 item 3 described the name as "a
    class or None depending on installed extras". It was measured to be the
    class even with spatialdata blocked -- the branch was dead, because the
    swallow one layer down meant the import it guarded never raised. Both
    the branch and the flag are gone; the name has one type.
    """
    assert isinstance(thyra.SpatialDataConverter, type)


class TestTheModuleTableMatchesTheDecorators:
    """The one thing lazy registration can get silently wrong (issue #381).

    Registration is still a decorator; what changed is that nothing imports
    every reader up front to run them, so the registry carries a table
    saying which module to import for a given format. A reader added with
    a decorator and not with a table entry would not be found at all, and a
    table entry naming the wrong module would import something that
    registers a different format -- both of which look like "that format is
    not supported" at the far end of a conversion the user has already
    started.

    So the table is checked against the source rather than against a second
    hand-written list. The decorator calls are read out of every module
    under ``thyra/readers/`` and ``thyra/converters/`` with ``ast``, which
    imports nothing: a test that imported them to find the registrations
    would populate the tables as a side effect and could no longer tell
    whether the table or the import had done it.
    """

    @staticmethod
    def _registered_format(node: ast.AST, function: str) -> Optional[str]:
        """The format name in ``function("name")``, or ``None``.

        Covers both spellings in the tree: ``@register_reader("imzml")`` on
        a class, and ``register_converter("spatialdata")(cls)`` called
        after the class body. ``ast.walk`` reaches the inner call of the
        second, so one test against ``ast.Call`` catches both.
        """
        if not isinstance(node, ast.Call):
            return None
        if not isinstance(node.func, ast.Name) or node.func.id != function:
            return None
        if len(node.args) != 1 or not isinstance(node.args[0], ast.Constant):
            return None
        value = node.args[0].value
        return value if isinstance(value, str) else None

    @classmethod
    def _declared(cls, package: str, function: str) -> Dict[str, str]:
        """``format -> dotted module path`` for every decorator in ``package``."""
        root = Path(thyra.__file__).resolve().parent.parent
        declared: Dict[str, str] = {}
        for path in sorted((root / "thyra" / package).rglob("*.py")):
            dotted = ".".join(path.relative_to(root).with_suffix("").parts)
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                format_name = cls._registered_format(node, function)
                if format_name is None:
                    continue
                assert format_name not in declared, (
                    f"{function}({format_name!r}) appears in both "
                    f"{declared[format_name]} and {dotted}; the table can "
                    f"only name one of them"
                )
                declared[format_name] = dotted
        return declared

    def test_every_reader_decorator_has_a_table_entry(self):
        assert self._declared("readers", "register_reader") == dict(
            registry_module._READER_MODULES
        )

    def test_every_converter_decorator_has_a_table_entry(self):
        assert self._declared("converters", "register_converter") == dict(
            registry_module._CONVERTER_MODULES
        )

    @pytest.mark.parametrize("format_name", sorted(registry_module._READER_MODULES))
    def test_each_entry_resolves_to_a_reader(self, format_name):
        """The table is only right if importing what it names registers.

        Equality with the decorator calls above is a source-level check;
        this one drives the lookup, so an entry naming a module that fails
        to import, or that registers under a different name than the
        decorator literal suggests, fails here.
        """
        assert issubclass(get_reader_class(format_name), BaseMSIReader)

    def test_the_converter_entry_resolves(self):
        assert issubclass(get_converter_class("spatialdata"), BaseMSIConverter)

    def test_a_miss_names_every_format_the_table_knows(self):
        """``Available:`` must not depend on what the process has touched.

        The refusal used to list the registered keys, which was every
        format while every reader was imported at startup. Deferred, that
        would be whatever happened to have been looked up already -- so a
        user who mistyped a format before converting anything would be told
        ``Available: []``, which is both wrong and unhelpful.
        """
        with pytest.raises(ConversionRefused) as excinfo:
            get_reader_class("mzml")

        message = str(excinfo.value)
        for format_name in registry_module._READER_MODULES:
            assert format_name in message


def test_a_format_is_imported_only_when_it_is_looked_up():
    """The deferral is real, not a table that something else fills in.

    In-process this cannot be seen: the test session has imported every
    reader long before this runs. The child imports ``thyra``, checks that
    the PHI reader's module is absent, resolves the format, and checks that
    it is now there -- which is the whole of issue #381 in four lines.

    PHI is the subject because nothing else in the package imports it, so
    an incidental import cannot make this pass.
    """
    source = """
import sys

import thyra
from thyra.core.registry import get_reader_class

module = "thyra.readers.phi.phi_reader"
assert module not in sys.modules, "importing thyra imported a reader"
reader = get_reader_class("phi")
assert module in sys.modules, "the lookup did not import the reader"
print("OK", reader.__name__)
"""
    proc = subprocess.run(
        [sys.executable, "-c", source],
        capture_output=True,
        text=True,
        cwd=Path(thyra.__file__).resolve().parent.parent,
    )

    assert proc.returncode == 0, f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
    assert proc.stdout.startswith("OK PhiReader")
