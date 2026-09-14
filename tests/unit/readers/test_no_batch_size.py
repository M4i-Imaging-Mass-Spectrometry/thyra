"""``batch_size`` is gone from the iteration surface (issue #309).

Sixteen ``iter_*`` methods across nine modules declared it and fourteen
never referenced the name at all. The two that did:

- **mzPeak** logged a debug line saying it was ignoring it.
- **imzML** chose between ``_iter_spectra_single`` and
  ``_iter_spectra_batch``, which walked the same flat ``range`` and called
  the same per-spectrum function -- the "batch" form was a nested loop over
  one flat one. Proved on a real fixture before the branch was removed:
  ``single(batch_size=1) == batch(1000) == default``.

Nothing reached it from the outside either. ``_buffer_size`` is assigned
once from ``config.DEFAULT_BUFFER_SIZE``, never reassigned, and no CLI flag
touches it.

The 16 methods are a plain delete, deliberately without a refusal: none of
them takes ``**kwargs``, so a stale caller already gets a loud
``TypeError`` naming the argument. **The constructor is the opposite
case** -- see :class:`TestTheConstructorRefusesIt`.
"""

from __future__ import annotations

import ast
import inspect
import pathlib

import pytest

from thyra.errors import ConversionRefused

_ITERATORS = (
    "iter_spectra",
    "iter_frame_scans",
    "iter_mobility_spectra",
    "iter_precursor_spectra",
)
_THYRA = pathlib.Path(__file__).resolve().parents[3] / "thyra"


def _declaring_modules():
    """Every module defining one of the four iterators, discovered."""
    out = []
    for path in sorted(_THYRA.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8-sig"))
        except SyntaxError:  # pragma: no cover - a BOM'd module, not a reader
            continue
        for node in ast.walk(tree):
            if (
                isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                and node.name in _ITERATORS
            ):
                out.append((path.relative_to(_THYRA.parent).as_posix(), node))
    return out


def test_the_iterators_were_found_at_all():
    """Guard the guard: an empty sweep would pass everything below."""
    assert len(_declaring_modules()) >= 16


@pytest.mark.parametrize(
    "path, node",
    _declaring_modules(),
    ids=[f"{p}:{n.name}" for p, n in _declaring_modules()],
)
def test_no_iterator_declares_batch_size(path, node):
    args = [a.arg for a in node.args.args] + [a.arg for a in node.args.kwonlyargs]
    assert "batch_size" not in args, f"{path}:{node.lineno} still declares it"


class TestAStaleCallerIsTold:
    """No refusal machinery here, and none is needed: the methods take no
    ``**kwargs``, so Python's own error already names the argument."""

    def test_passing_it_raises_type_error(self):
        from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader

        reader = MockMSIReader(MockMSIConfig(n_x=2, n_y=2, n_z=1, n_mz_bins=64))
        with pytest.raises(TypeError, match="batch_size"):
            list(reader.iter_spectra(batch_size=10))


class TestTheConstructorRefusesIt:
    """``ImzMLReader.__init__`` cannot just drop it.

    It takes ``**kwargs``, which forwards to ``BaseMSIReader``, which never
    reads them -- so a plain delete would make ``batch_size=50`` silently
    accepted and no longer documented. D10's Python-API half says a keyword
    that is accepted and does nothing must be answered, and ``BrukerReader``
    already refuses this exact name, so before this the same keyword had two
    answers depending on which reader you held.
    """

    def test_the_keyword_is_refused(self):
        from thyra.readers.imzml.imzml_reader import ImzMLReader

        with pytest.raises(ConversionRefused, match="batch_size"):
            ImzMLReader(pathlib.Path("nonexistent.imzML"), batch_size=50)

    def test_the_message_says_there_is_no_replacement(self):
        from thyra.readers.imzml.imzml_reader import ImzMLReader

        with pytest.raises(ConversionRefused) as excinfo:
            ImzMLReader(pathlib.Path("nonexistent.imzML"), batch_size=50)
        assert "no replacement keyword" in str(excinfo.value)

    def test_a_positional_call_written_against_the_old_signature(self):
        """``cache_coordinates`` moved into the vacated slot.

        ``ImzMLReader(path, 50)`` used to mean ``batch_size=50``; without
        this guard it would now mean ``cache_coordinates=50``, which is
        truthy and silently wrong rather than an error. A keyword refusal
        cannot see a positional call.
        """
        from thyra.readers.imzml.imzml_reader import ImzMLReader

        with pytest.raises(ConversionRefused, match="cache_coordinates"):
            ImzMLReader(pathlib.Path("nonexistent.imzML"), 50)

    def test_a_real_bool_still_works(self):
        """The guard must not refuse the call it was written to protect."""
        from thyra.readers.imzml.imzml_reader import _refuse_shifted_cache_coordinates

        _refuse_shifted_cache_coordinates(True)
        _refuse_shifted_cache_coordinates(False)

    def test_it_is_a_conversion_refused_not_a_bare_value_error(self):
        """Read from the raise, never inferred from a test: ``pytest.raises``
        on ``ValueError`` would pass either way, since ``ConversionRefused``
        subclasses it."""
        from thyra.readers.imzml import imzml_reader

        source = inspect.getsource(imzml_reader._refuse_retired_keywords)
        assert "raise ConversionRefused(" in source


class TestTheImzmlBranchIsGone:
    def test_the_batch_helper_no_longer_exists(self):
        from thyra.readers.imzml.imzml_reader import ImzMLReader

        assert not hasattr(ImzMLReader, "_iter_spectra_batch")

    def test_the_single_loop_remains(self):
        from thyra.readers.imzml.imzml_reader import ImzMLReader

        assert hasattr(ImzMLReader, "_iter_spectra_single")
