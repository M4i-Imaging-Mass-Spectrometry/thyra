"""The z plane a spectrum lands on comes from the file, not from a constant.

imzML does not pin down a base for ``IMS:1000052`` and pyimzml is
inconsistent about it: when the term is absent it synthesises ``z = 1``,
but an explicit ``z = 0`` is passed through verbatim -- contradicting its
own docstring, which promises zero.

``np.maximum(z - 1, 0)`` folded both onto plane 0. On a file written
0-based that merges planes 0 and 1: a two-plane acquisition converts to a
two-row table whose rows are the *union* of both planes, the conversion
returns ``True`` and nothing warns. Rebasing on the smallest z present is
correct for either convention, and answers 1 -- the pre-existing value --
for every file in the corpus, none of which declares ``IMS:1000052``.

The same constant appeared one layer up in the extractor, deciding both
``dimensions[2]`` and which row of the CSR ``indptr`` a pixel's peak count
lands on, so both layers are exercised here.
"""

from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.readers.imzml.imzml_reader import ImzMLReader

_MZS = np.linspace(100.0, 200.0, 5)
_INTENSITIES = np.arange(1.0, 6.0)


def _write_two_planes(directory: Path, name: str, z_values) -> Path:
    """A 2x1x2 acquisition whose z values are exactly ``z_values``.

    Every plane carries the same m/z values with plane-dependent
    intensities, so a merge shows up as summed intensity rather than as a
    shape change alone.
    """
    path = directory / f"{name}.imzML"
    with ImzMLWriter(str(path), mode="processed") as writer:
        for plane, z in enumerate(z_values):
            for x in (1, 2):
                writer.addSpectrum(_MZS, _INTENSITIES * (plane + 1), (x, 1, z))
    return path


@pytest.mark.parametrize(
    "z_values, description",
    [((0, 1), "written 0-based"), ((1, 2), "written 1-based")],
)
class TestZBaseIsReadOffTheFile:
    """Both conventions must describe the same acquisition identically."""

    def test_two_planes_are_two_planes(self, tmp_path, z_values, description):
        """The headline: 0-based files reported ``n_z = 1``."""
        path = _write_two_planes(tmp_path, f"dims_{z_values[0]}", z_values)
        reader = ImzMLReader(path)
        try:
            assert reader.get_essential_metadata().dimensions == (2, 1, 2), description
        finally:
            reader.close()

    def test_coordinates_span_both_planes(self, tmp_path, z_values, description):
        """The cached array is what the converters index rows by."""
        path = _write_two_planes(tmp_path, f"coords_{z_values[0]}", z_values)
        reader = ImzMLReader(path)
        try:
            coords = [c for c, _, _ in reader.iter_spectra()]
        finally:
            reader.close()

        assert sorted(coords) == [(0, 0, 0), (0, 0, 1), (1, 0, 0), (1, 0, 1)]

    def test_no_two_spectra_share_a_position(self, tmp_path, z_values, description):
        """The property the merge violated, stated directly."""
        path = _write_two_planes(tmp_path, f"unique_{z_values[0]}", z_values)
        reader = ImzMLReader(path)
        try:
            coords = [c for c, _, _ in reader.iter_spectra()]
        finally:
            reader.close()

        assert len(set(coords)) == len(coords)

    def test_peak_counts_land_on_their_own_row(self, tmp_path, z_values, description):
        """These counts become the CSR ``indptr``.

        With the base wrong, plane 1's counts overwrote plane 0's and the
        second half of the array stayed zero -- so the matrix was built
        with the wrong number of entries reserved per row.
        """
        path = _write_two_planes(tmp_path, f"counts_{z_values[0]}", z_values)
        reader = ImzMLReader(path)
        try:
            counts = reader.get_peak_counts_per_pixel()
        finally:
            reader.close()

        assert counts is not None
        assert list(counts) == [len(_MZS)] * 4


class TestZBaseWithoutTheCoordinateCache:
    """``cache_coordinates=False`` takes a different branch; it must agree."""

    def test_the_uncached_path_rebases_too(self, tmp_path):
        path = _write_two_planes(tmp_path, "uncached", (0, 1))
        reader = ImzMLReader(path, cache_coordinates=False)
        try:
            coords = [c for c, _, _ in reader.iter_spectra()]
        finally:
            reader.close()

        assert sorted(coords) == [(0, 0, 0), (0, 0, 1), (1, 0, 0), (1, 0, 1)]


class TestSinglePlaneFilesAreUnchanged:
    """The corpus shape: no ``IMS:1000052``, so pyimzml synthesises z = 1.

    This is the regression guard for the fix itself -- every real file
    goes down this path and must come out exactly as before.
    """

    @pytest.mark.parametrize("z", [0, 1])
    def test_one_plane_is_plane_zero(self, tmp_path, z):
        path = _write_two_planes(tmp_path, f"single_{z}", (z,))
        reader = ImzMLReader(path)
        try:
            assert reader.get_essential_metadata().dimensions == (2, 1, 1)
            assert sorted(c for c, _, _ in reader.iter_spectra()) == [
                (0, 0, 0),
                (1, 0, 0),
            ]
        finally:
            reader.close()
