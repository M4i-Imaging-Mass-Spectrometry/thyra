"""The CV constants moved, and both spellings still reach the same objects.

:mod:`thyra.metadata.constants` was carved out of
``thyra/resampling/constants.py`` under issue #381. Three vendor metadata
extractors and the imzML reader imported accession codes and spectrum-type
names from the resampling package, which is the wrong direction --
resampling consumes a description of the data, it does not define one --
and it is why importing a metadata extractor loaded all 25 resampling
modules.

A move like that breaks quietly in two ways, and both are asserted here:
an importer left pointing at the old module, and a re-export that is a
copy rather than the same object, so ``is`` comparisons and ``in``
membership on the alias table start disagreeing between call sites.
"""

import pytest

from thyra.metadata import constants as new_home
from thyra.resampling import constants as old_home

_MOVED = [
    "SPECTRUM_TYPE_ALIASES",
    "BinaryDataType",
    "ImzMLAccessions",
    "SpectrumType",
    "normalize_spectrum_type",
]


@pytest.mark.parametrize("name", _MOVED)
def test_the_old_module_still_serves_the_same_object(name):
    assert getattr(old_home, name) is getattr(new_home, name)


@pytest.mark.parametrize("name", ["BinaryDataType", "ImzMLAccessions", "SpectrumType"])
def test_the_resampling_package_still_re_exports_it(name):
    """``from thyra.resampling import SpectrumType`` is the documented spelling.

    These three are the moved names that ``thyra/resampling/__init__.py``
    lists in its ``__all__``; the package is a public surface and the move
    is not a breaking change. ``SPECTRUM_TYPE_ALIASES`` and
    ``normalize_spectrum_type`` were never re-exported at that level, and
    are not now.
    """
    import thyra.resampling as package

    assert getattr(package, name) is getattr(new_home, name)


def test_the_resampling_threshold_did_not_move():
    """``Thresholds`` is the one number in that module that is about resampling.

    It is the peak density above which a spectrum is called profile, read
    by ``data_characteristics`` and by nothing under ``thyra/metadata/``.
    Moving it too would have traded one wrong home for another.
    """
    assert old_home.Thresholds.PROFILE_PEAK_DENSITY == 5000
    assert not hasattr(new_home, "Thresholds")


class TestNormalizeSpectrumType:
    """Behaviour carried over unchanged, pinned at its new address."""

    @pytest.mark.parametrize(
        "given, expected",
        [
            ("profile", new_home.SpectrumType.PROFILE),
            ("PROFILE", new_home.SpectrumType.PROFILE),
            ("  profile spectrum ", new_home.SpectrumType.PROFILE),
            ("centroid", new_home.SpectrumType.CENTROID),
            ("centroided", new_home.SpectrumType.CENTROID),
            ("centroid spectrum", new_home.SpectrumType.CENTROID),
        ],
    )
    def test_it_maps_every_accepted_spelling(self, given, expected):
        assert new_home.normalize_spectrum_type(given) == expected

    def test_none_means_no_override(self):
        assert new_home.normalize_spectrum_type(None) is None

    def test_an_unknown_spelling_is_refused(self):
        """A typo must not be indistinguishable from "no override"."""
        from thyra.errors import ConversionRefused

        with pytest.raises(ConversionRefused, match="Unknown spectrum_type"):
            new_home.normalize_spectrum_type("proflie")
