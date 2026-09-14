# tests/unit/readers/test_reader_conformance.py
"""One contract, every reader, one set of assertions (issue #274).

Seven readers implement :class:`~thyra.core.base_reader.BaseMSIReader`, and
until this file each was tested individually and differently. Nothing tested
the contract *as a contract*, so nothing noticed when a reader answered the
shared interface in a new way -- which is the mechanism behind the three
different spellings of "unsupported" filed as #275. The interface could drift
because nothing held it still.

**What holds it still.** The case table below is keyed by the registry's own
format names, and :func:`test_every_registered_format_is_covered` compares it
against ``MSIRegistry._readers``. Registering a reader without either a
conformance case or an explicit exemption fails that test, so a new reader is
covered the day it lands rather than the day someone remembers.

**Why the factories build their own data.** Each case writes a synthetic
acquisition into ``tmp_path`` and opens a real reader on it. No mocks: a mock
that satisfies the contract proves the mock conforms, not the reader. The
builders are the ones the per-reader suites already use, imported rather than
duplicated so the two cannot drift apart.

**The one exemption.** Waters has no synthetic path. Its data source *is* the
MassLynx native library -- the reader hands it an opaque file handle and reads
everything back through it -- so there is nothing to write into ``tmp_path``
that would exercise the real code. Its own suite patches ``MassLynxLib``,
which is the right test for that reader and the wrong one for this file.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Callable, Dict

import numpy as np
import pytest

from thyra.core.base_reader import BaseMSIReader
from thyra.core.registry import _registry
from thyra.readers.bruker import SolarixReader
from thyra.readers.bruker.rapiflex import RapiflexReader
from thyra.readers.bruker.timstof.timstof_reader import BrukerReader
from thyra.readers.imzml.imzml_reader import ImzMLReader
from thyra.readers.mzpeak import MzPeakReader
from thyra.readers.phi import PhiReader
from thyra.utils.bruker_exceptions import SDKError

from ...fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from .test_phi_reader import block as phi_block
from .test_phi_reader import event, events_block, write_raw
from .test_rapiflex_reader import build_rapiflex_dataset
from .test_solarix_reader import make_solarix_d

_TDF_FIXTURE = (
    Path(__file__).resolve().parents[2] / "data" / "fixtures" / "synthetic_tims.d"
)


# ----------------------------------------------------------------------
# Factories: tmp_path -> an open reader on synthetic data
# ----------------------------------------------------------------------


def _imzml(tmp_path: Path) -> ImzMLReader:
    from pyimzml.ImzMLWriter import ImzMLWriter

    path = tmp_path / "conformance.imzML"
    mzs = np.linspace(100.0, 400.0, 12)
    with ImzMLWriter(str(path), mode="processed") as writer:
        for y in (1, 2):
            for x in (1, 2):
                writer.addSpectrum(mzs, np.full(mzs.size, float(x + y)), (x, y, 1))
    return ImzMLReader(path)


def _solarix(tmp_path: Path) -> SolarixReader:
    return SolarixReader(make_solarix_d(tmp_path))


def _mzpeak(tmp_path: Path) -> MzPeakReader:
    path = build_mzpeak(tmp_path / "conformance.mzpeak", grid_spectra(2, 2))
    return MzPeakReader(path)


def _phi(tmp_path: Path) -> PhiReader:
    blocks = (
        events_block([event(0, 0, 3_000_000), event(0, 0, 4_000_000)])
        + events_block([event(1, 2, 5_000_000)])
        + phi_block(2)
    )
    return PhiReader(write_raw(tmp_path / "conformance.raw", blocks))


def _rapiflex(tmp_path: Path) -> RapiflexReader:
    folder, *_ = build_rapiflex_dataset(tmp_path)
    return RapiflexReader(folder)


def _bruker_tdf(tmp_path: Path) -> BrukerReader:
    # Copied rather than opened in place: the reader takes SQLite handles on
    # the acquisition, and the fixture is committed to the repository.
    d_dir = tmp_path / _TDF_FIXTURE.name
    shutil.copytree(_TDF_FIXTURE, d_dir)
    try:
        return BrukerReader(d_dir, tdf_spectrum="scan_sum")
    except (SDKError, OSError) as exc:  # the vendor library is not loadable here
        pytest.skip(f"Bruker library not loadable on this platform: {exc}")


CASES: Dict[str, Callable[[Path], BaseMSIReader]] = {
    "imzml": _imzml,
    "solarix": _solarix,
    "mzpeak": _mzpeak,
    "phi": _phi,
    "rapiflex": _rapiflex,
    "bruker": _bruker_tdf,
}

EXEMPT: Dict[str, str] = {
    "waters": (
        "The MassLynx native library is the data source, so there is no "
        "synthetic acquisition to write; see the module docstring."
    ),
}


@pytest.fixture(params=sorted(CASES))
def reader(request, tmp_path: Path):
    """An open reader on synthetic data, one per registered format."""
    instance = CASES[request.param](tmp_path)
    try:
        yield instance
    finally:
        instance.close()


class TestTheCaseTableCoversTheRegistry:
    """The net: a reader cannot be registered without being covered here."""

    def test_every_registered_format_is_covered(self):
        registered = set(_registry._readers)
        covered = set(CASES) | set(EXEMPT)

        assert registered == covered, (
            "every registered reader needs a conformance case or an explicit "
            f"exemption; uncovered: {sorted(registered - covered)}, "
            f"stale: {sorted(covered - registered)}"
        )

    def test_every_case_builds_a_reader_subclass(self, reader):
        assert isinstance(reader, BaseMSIReader)

    @pytest.mark.parametrize("format_name", sorted(EXEMPT))
    def test_exemptions_name_a_registered_format(self, format_name):
        assert format_name in _registry._readers


class TestSpectrumIteration:
    """``iter_spectra`` yields ``((x, y, z), mzs, intensities)``."""

    def test_yields_matching_arrays_at_integer_coordinates(self, reader):
        seen = 0
        for coords, mzs, intensities in reader.iter_spectra():
            seen += 1
            x, y, z = coords
            assert all(isinstance(v, (int, np.integer)) for v in (x, y, z))
            assert min(x, y, z) >= 0, "coordinates are 0-based"
            assert mzs.shape == intensities.shape
            assert mzs.ndim == 1
            if seen >= 4:  # the contract is per-spectrum; four is enough
                break

        assert seen > 0, "a reader on a non-empty acquisition yielded nothing"

    def test_coordinates_fall_inside_the_declared_grid(self, reader):
        n_x, n_y, n_z = reader.get_essential_metadata().dimensions

        for (x, y, z), _, _ in reader.iter_spectra():
            assert x < n_x and y < n_y and z < n_z

    def test_values_are_finite(self, reader):
        for _, mzs, intensities in reader.iter_spectra():
            assert np.all(np.isfinite(mzs))
            assert np.all(np.isfinite(intensities))
            break


class TestTheCommonMassAxis:
    """``get_common_mass_axis`` must always return a usable axis."""

    def test_is_finite_and_strictly_increasing(self, reader):
        axis = reader.get_common_mass_axis()

        assert axis.ndim == 1
        assert axis.size > 0
        assert np.all(np.isfinite(axis))
        assert np.all(np.diff(axis) > 0), "axis is not strictly increasing"

    def test_covers_the_declared_mass_range(self, reader):
        axis = reader.get_common_mass_axis()
        low, high = reader.get_essential_metadata().mass_range

        assert axis[0] >= low - abs(low) * 1e-6
        assert axis[-1] <= high + abs(high) * 1e-6


class TestCapabilityConventions:
    """A predicate that says yes must be backed by a method that works.

    These pin the conventions as they stand. #275 proposes replacing them
    with capability ``Protocol``s; this class is the net under that change,
    and the assertions it has to rewrite are exactly the inconsistency.
    """

    def test_frame_scans_predicate_matches_the_method(self, reader):
        if reader.has_frame_scans:
            next(iter(reader.iter_frame_scans()), None)
        else:
            with pytest.raises(NotImplementedError):
                next(iter(reader.iter_frame_scans()), None)

    def test_ion_mobility_predicate_matches_the_methods(self, reader):
        if reader.has_ion_mobility:
            assert reader.get_mobility_axis() is not None
            next(iter(reader.iter_mobility_spectra()), None)
        else:
            assert reader.get_mobility_axis() is None

    def test_fragmentation_predicate_matches_the_method(self, reader):
        """``get_fragmentation`` is not the demultiplexing predicate.

        #275 lists this as "returns ``None``, and ``iter_precursor_spectra``
        raises", which reads as one capability. It is two. A non-``None``
        schedule means *the MS level is known*: an MS1 acquisition returns
        ``FragmentationSchedule(ms_level=1)`` with no windows, which is a
        real answer and not a promise that precursors can be separated.
        The demultiplexing capability is ``schedule.windows`` being
        non-empty. The synthetic TDF fixture is MS1, so it takes the second
        branch -- and asserting the naive reading here fails against it.
        """
        schedule = reader.get_fragmentation()

        if schedule is not None and schedule.windows:
            next(iter(reader.iter_precursor_spectra()), None)
        else:
            with pytest.raises(NotImplementedError):
                next(iter(reader.iter_precursor_spectra()), None)

    def test_shared_mass_axis_is_a_bool(self, reader):
        assert isinstance(reader.has_shared_mass_axis, bool)
        assert isinstance(reader.has_ion_mobility, bool)
        assert isinstance(reader.has_frame_scans, bool)


class TestLifetime:
    """``reset`` rewinds, ``close`` is idempotent, ``with`` closes."""

    def test_reset_restores_iteration_from_the_first_spectrum(self, reader):
        first = next(iter(reader.iter_spectra()))

        reader.reset()
        again = next(iter(reader.iter_spectra()))

        assert first[0] == again[0]
        np.testing.assert_array_equal(first[1], again[1])

    def test_close_is_idempotent(self, reader):
        reader.close()
        reader.close()  # must not raise

    @pytest.mark.parametrize("format_name", sorted(CASES))
    def test_context_manager_closes(self, format_name, tmp_path):
        # Built directly rather than through the `reader` fixture: this test
        # owns the lifetime it is asserting on.
        instance = CASES[format_name](tmp_path)

        with instance as entered:
            assert entered is instance

        instance.close()  # already closed by __exit__; must still not raise
