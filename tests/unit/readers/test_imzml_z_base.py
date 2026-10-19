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

import re
from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from thyra import convert_msi
from thyra.errors import ConversionRefused
from thyra.readers.imzml.imzml_reader import ImzMLReader
from thyra.readers.mzpeak import MzPeakReader
from thyra.utils.imzml_coordinate_base import PLANE_RULE

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


_GRID = ((1, 1), (2, 1), (1, 2), (2, 2))
_Z_TERM = re.compile(r"(<cvParam [^>]*IMS:1000052[^>]*/>)")


def _write_grid(directory: Path, name: str, edit_z) -> Path:
    """Four spectra that all state z = 1, with ``edit_z`` applied to the terms.

    ``edit_z`` gets the list of the four position z terms as written and
    returns what each is replaced with.
    """
    path = directory / f"{name}.imzML"
    with ImzMLWriter(str(path), mode="continuous") as writer:
        for x, y in _GRID:
            writer.addSpectrum(_MZS, _INTENSITIES, (x, y, 1))
    parts = _Z_TERM.split(path.read_text())
    at = [i for i, part in enumerate(parts) if "IMS:1000052" in part]
    assert len(at) == 4
    for i, new in zip(at, edit_z([parts[i] for i in at])):
        parts[i] = new
    path.write_text("".join(parts))
    return path


def _without_z(unstated):
    """Drop the z term of the spectra at these positions in the grid."""
    return lambda terms: ["" if k in unstated else t for k, t in enumerate(terms)]


def _essential(path: Path):
    reader = ImzMLReader(path)
    try:
        return reader.get_essential_metadata()
    finally:
        reader.close()


class TestSomeSpectraStateZ:
    """Four spectra, three with z and one without: a plane cannot be stated.

    The spectrum without z used to land on a plane made from pyimzml's
    stand-in of 1; the recorded z offset then depended on which spectrum
    came first.
    """

    @pytest.mark.parametrize("lacking", [0, 1, 3], ids=["first", "middle", "last"])
    def test_the_reader_refuses_it_with_the_plane_rule(self, tmp_path, lacking):
        path = _write_grid(tmp_path, f"lacks_{lacking}", _without_z({lacking}))

        with pytest.raises(ConversionRefused, match="3 of 4 spectra a position z") as e:
            _essential(path)

        assert PLANE_RULE in str(e.value)
        assert "IMS:1000052" in str(e.value)

    def test_the_conversion_returns_false_and_writes_no_store(self, tmp_path):
        path = _write_grid(tmp_path, "lacks_first", _without_z({0}))
        store = tmp_path / "out.zarr"

        assert convert_msi(path, store, pixel_size_um=10.0) is False

        assert not store.exists()

    def test_imzml_and_mzpeak_refuse_with_the_same_sentence(self, tmp_path):
        """One rule, one sentence, whichever route the file takes."""
        path = _write_grid(tmp_path, "lacks_last", _without_z({3}))
        archive = build_mzpeak(
            tmp_path / "somez.mzpeak", grid_spectra(2, 2), position_z=[1, 1, None, 1]
        )
        marker = "a position z (IMS:1000052) and the others none."

        with pytest.raises(ConversionRefused) as imzml_refusal:
            _essential(path)
        with pytest.raises(ConversionRefused) as mzpeak_refusal:
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

        imzml_tail = str(imzml_refusal.value).split(marker)
        mzpeak_tail = str(mzpeak_refusal.value).split(marker)
        assert len(imzml_tail) == len(mzpeak_tail) == 2
        assert imzml_tail[1] == mzpeak_tail[1]


class TestAZTermWithoutAValue:
    """A z term that carries no value states no z, so the offset is 0.

    It used to count as stated, recording 1 for what was only pyimzml's
    stand-in.
    """

    def test_every_term_valueless_records_zero(self, tmp_path):
        path = _write_grid(
            tmp_path,
            "valueless",
            lambda terms: [t.replace(' value="1"', "") for t in terms],
        )
        assert " value=" not in "".join(_Z_TERM.findall(path.read_text()))

        assert _essential(path).coordinate_offsets == (1, 1, 0)

    def test_the_same_terms_with_a_value_record_one(self, tmp_path):
        """Control: the valueless edit is the only difference."""
        path = _write_grid(tmp_path, "valued", lambda terms: terms)

        assert _essential(path).coordinate_offsets == (1, 1, 1)


def _z_through_groups(path: Path, planes) -> None:
    """Move each scan's z into a group the scan references, as a writer may."""
    text = path.read_text()
    pattern = (
        r'(<referenceableParamGroupRef ref="scan1"/>)(\s*<cvParam[^>]*IMS:1000050'
        r"[^>]*/>\s*<cvParam[^>]*IMS:1000051[^>]*/>)\s*<cvParam[^>]*IMS:1000052"
        r'[^>]*value="(\d+)"/>'
    )
    text, n = re.subn(pattern, r'\1<referenceableParamGroupRef ref="z\3"/>\2', text)
    assert n == len(planes) * 4
    groups = "".join(
        f'<referenceableParamGroup id="z{z}"><cvParam cvRef="IMS" '
        f'accession="IMS:1000052" name="position z" value="{z}"/>'
        f"</referenceableParamGroup>"
        for z in sorted(set(planes))
    )
    path.write_text(
        text.replace(
            "</referenceableParamGroupList>", groups + "</referenceableParamGroupList>"
        )
    )


def _write_grouped_planes(directory: Path, planes=(1, 2)) -> Path:
    path = directory / "grouped.imzML"
    with ImzMLWriter(str(path), mode="continuous") as writer:
        for z in planes:
            for x, y in _GRID:
                writer.addSpectrum(_MZS, _INTENSITIES * 10**z, (x, y, z))
    _z_through_groups(path, planes)
    return path


class TestZStatedThroughAGroup:
    """pyimzml does not read z from a referenced group.

    Two planes stated that way used to be summed into one table, and a
    single plane at z = 3 recorded offset 0.
    """

    def test_the_reader_refuses_and_names_the_group(self, tmp_path):
        path = _write_grouped_planes(tmp_path)

        with pytest.raises(ConversionRefused) as e:
            _essential(path)

        assert "through a referenced parameter group" in str(e.value)
        assert "parameter group 'z1'" in str(e.value)

    def test_one_plane_through_a_group_is_refused_too(self, tmp_path):
        path = _write_grouped_planes(tmp_path, planes=(3,))

        with pytest.raises(ConversionRefused, match="parameter group 'z3'"):
            _essential(path)

    def test_the_conversion_returns_false_and_writes_no_store(self, tmp_path):
        path = _write_grouped_planes(tmp_path)
        store = tmp_path / "out.zarr"

        assert convert_msi(path, store, pixel_size_um=10.0) is False

        assert not store.exists()
