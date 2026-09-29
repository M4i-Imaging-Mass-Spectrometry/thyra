"""Reader-level tests for the experimental mzPeak reader.

Every archive here is built by :mod:`tests.fixtures.mzpeak_builder`, which
spells the container out literally rather than asking the reader how to write
it -- see that module's docstring for why that separation matters.
"""

from __future__ import annotations

import json
import warnings
import zipfile

import numpy as np
import pytest

from tests.fixtures.mzpeak_builder import Spectrum, build_mzpeak, grid_spectra
from thyra.errors import ConversionRefused
from thyra.readers.mzpeak import MzPeakReader


@pytest.fixture
def simple_archive(tmp_path):
    """A 3x2 acquisition, 1-based positions, eight points per spectrum."""
    return build_mzpeak(tmp_path / "simple.mzpeak", grid_spectra(3, 2))


class TestIteration:
    """Coordinates, payloads and ordering coming out of iter_spectra."""

    def test_yields_every_spectrum_in_index_order(self, simple_archive):
        """All spectra arrive, ascending, with 0-based coordinates."""
        with MzPeakReader(simple_archive) as reader:
            emitted = list(reader.iter_spectra())

        assert len(emitted) == 6
        assert [coords for coords, _, _ in emitted] == [
            (0, 0, 0),
            (1, 0, 0),
            (2, 0, 0),
            (0, 1, 0),
            (1, 1, 0),
            (2, 1, 0),
        ]

    def test_payload_matches_what_was_written(self, tmp_path):
        """m/z and intensity survive the round trip, as float64."""
        spectra = [Spectrum(1, 1, [100.0, 200.5, 300.25], [7.0, 8.0, 9.0])]
        archive = build_mzpeak(tmp_path / "one.mzpeak", spectra)

        with MzPeakReader(archive) as reader:
            ((_, mzs, intensities),) = list(reader.iter_spectra())

        np.testing.assert_array_equal(mzs, [100.0, 200.5, 300.25])
        np.testing.assert_allclose(intensities, [7.0, 8.0, 9.0])
        assert mzs.dtype == np.float64
        assert intensities.dtype == np.float64

    def test_spectrum_straddling_row_groups_is_emitted_once(self, tmp_path):
        """A spectrum split across row groups is reassembled, not duplicated.

        Row groups cut at a fixed row count with no regard for spectrum
        boundaries, so with eight points per spectrum and three rows per group
        every spectrum spans several groups.
        """
        archive = build_mzpeak(
            tmp_path / "straddle.mzpeak",
            grid_spectra(2, 2, n_points=8),
            row_group_size=3,
        )

        with MzPeakReader(archive) as reader:
            emitted = list(reader.iter_spectra())

        assert len(emitted) == 4
        assert [coords for coords, _, _ in emitted] == [
            (0, 0, 0),
            (1, 0, 0),
            (0, 1, 0),
            (1, 1, 0),
        ]
        for _, mzs, intensities in emitted:
            assert mzs.size == 8
            assert intensities.size == 8
            # Ascending within the spectrum proves the pieces were joined in
            # order rather than concatenated arbitrarily.
            assert np.all(np.diff(mzs) > 0)

    def test_row_group_size_does_not_change_the_result(self, tmp_path):
        """Reading is invariant to how the writer chose its row groups."""
        spectra = grid_spectra(3, 2, n_points=7)
        whole = build_mzpeak(tmp_path / "whole.mzpeak", spectra)
        split = build_mzpeak(tmp_path / "split.mzpeak", spectra, row_group_size=2)

        with MzPeakReader(whole) as reader:
            reference = list(reader.iter_spectra())
        with MzPeakReader(split) as reader:
            candidate = list(reader.iter_spectra())

        assert len(reference) == len(candidate)
        for (coords_a, mz_a, int_a), (coords_b, mz_b, int_b) in zip(
            reference, candidate
        ):
            assert coords_a == coords_b
            np.testing.assert_array_equal(mz_a, mz_b)
            np.testing.assert_array_equal(int_a, int_b)

    def test_zero_based_positions_are_normalised(self, tmp_path):
        """A 0-based file lands on the same grid as a 1-based one.

        Normalisation subtracts the observed minimum rather than a constant,
        because the format promises nothing about the origin.
        """
        archive = build_mzpeak(tmp_path / "zero.mzpeak", grid_spectra(2, 2, base=0))

        with MzPeakReader(archive) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]

        assert coords == [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)]

    def test_offset_positions_are_normalised(self, tmp_path):
        """A file whose positions start at 5 still converts to a 0-based grid."""
        spectra = [
            Spectrum(5, 9, [100.0, 101.0], [1.0, 2.0]),
            Spectrum(6, 9, [100.0, 101.0], [3.0, 4.0]),
            Spectrum(5, 10, [100.0, 101.0], [5.0, 6.0]),
        ]
        archive = build_mzpeak(tmp_path / "offset.mzpeak", spectra)

        with MzPeakReader(archive) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]
            essential = reader.get_essential_metadata()

        assert coords == [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
        assert essential.coordinate_offsets == (5, 9, 0)

    def test_missing_pixels_are_left_missing(self, tmp_path):
        """A sparse acquisition yields only acquired pixels.

        Missing pixels are ordinary in imaging mzPeak -- they were 38% of the
        real dataset this reader was designed against -- and the converter's
        sparse grid handles the gaps. Densifying here would invent spectra.
        """
        spectra = grid_spectra(3, 3, skip=[(2, 2), (3, 1)])
        archive = build_mzpeak(tmp_path / "sparse.mzpeak", spectra)

        with MzPeakReader(archive) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]
            essential = reader.get_essential_metadata()

        assert len(coords) == 7
        assert (1, 1, 0) not in coords  # the skipped (2, 2) in file frame
        assert (2, 0, 0) not in coords  # the skipped (3, 1) in file frame
        assert essential.dimensions == (3, 3, 1)
        assert essential.n_spectra == 7

    def test_the_acquisition_order_is_the_spectrum_index(self, tmp_path):
        """The archive's own index, not the raster: a serpentine second row."""
        first, second = grid_spectra(2, 1), grid_spectra(2, 2)[2:]
        spectra = first + list(reversed(second))
        archive = build_mzpeak(tmp_path / "serpentine.mzpeak", spectra)

        with MzPeakReader(archive) as reader:
            assert reader.has_acquisition_order
            ordered = [
                (coords, order)
                for coords, order, _, _ in reader.iter_spectra_with_acquisition_order()
            ]

        assert ordered == [
            ((0, 0, 0), 0),
            ((1, 0, 0), 1),
            ((1, 1, 0), 2),
            ((0, 1, 0), 3),
        ]


class TestMassAxis:
    """The reader's view of the m/z axis."""

    def test_never_claims_a_shared_axis(self, simple_archive):
        """The point layout, the one Thyra reads, has an axis per spectrum."""
        with MzPeakReader(simple_archive) as reader:
            assert reader.has_shared_mass_axis is False

    def test_common_axis_is_the_sorted_union(self, tmp_path):
        """Per-spectrum axes combine into one ascending, deduplicated axis."""
        spectra = [
            Spectrum(1, 1, [100.0, 102.0], [1.0, 2.0]),
            Spectrum(2, 1, [101.0, 102.0], [3.0, 4.0]),
        ]
        archive = build_mzpeak(tmp_path / "union.mzpeak", spectra)

        with MzPeakReader(archive) as reader:
            axis = reader.get_common_mass_axis()

        np.testing.assert_array_equal(axis, [100.0, 101.0, 102.0])


class TestNullPairPadding:
    """Null pairs mark removed zero runs and must not reach the converter."""

    def test_padding_is_dropped_from_the_payload(self, tmp_path):
        """Rows whose m/z and intensity are both null never surface."""
        spectra = grid_spectra(2, 1, n_points=6)
        archive = build_mzpeak(tmp_path / "padded.mzpeak", spectra, null_pair_after=3)

        with MzPeakReader(archive) as reader:
            emitted = list(reader.iter_spectra())

        for _, mzs, intensities in emitted:
            assert mzs.size == 6
            assert intensities.size == 6
            assert not np.isnan(mzs).any()
            assert not np.isnan(intensities).any()

    def test_the_dropped_count_is_per_iteration(self, tmp_path):
        """Every conversion iterates twice, and the count must not accumulate.

        ``_dropped_points`` was set once in ``__init__`` and only ever
        incremented, so the second pass logged the sum of both passes --
        12 dropped points, then 24, then 36 on a third iteration of the
        same file (issue #238). Cleared as ``iter_spectra`` starts, which
        also covers a caller that iterates again without ``reset()``.
        """
        spectra = grid_spectra(3, 2, n_points=6)
        archive = build_mzpeak(tmp_path / "counted.mzpeak", spectra, null_pair_after=3)

        counts = []
        with MzPeakReader(archive) as reader:
            for _ in range(3):
                reader.reset()
                list(reader.iter_spectra())
                counts.append(reader._dropped_points)
            list(reader.iter_spectra())  # no reset() at all
            counts.append(reader._dropped_points)

        assert counts[0] > 0
        assert counts == [counts[0]] * 4

    def test_padding_is_excluded_from_the_mass_axis(self, tmp_path):
        """The axis holds only channels that can take a value."""
        spectra = grid_spectra(2, 1, n_points=6)
        archive = build_mzpeak(
            tmp_path / "padded_axis.mzpeak", spectra, null_pair_after=3
        )

        with MzPeakReader(archive) as reader:
            axis = reader.get_common_mass_axis()

        assert not np.isnan(axis).any()
        assert np.all(np.diff(axis) > 0)

    def test_peak_counts_are_corrected_for_padding(self, tmp_path):
        """Recorded point counts include padding; the reported total does not.

        ``number_of_data_points`` counts stored rows, so a file with a null
        pair per spectrum records two more per spectrum than it can deliver.
        Reporting the recorded figure would have the converter pre-allocate
        for points that never arrive.
        """
        spectra = grid_spectra(2, 1, n_points=6)
        archive = build_mzpeak(
            tmp_path / "padded_counts.mzpeak", spectra, null_pair_after=3
        )

        with MzPeakReader(archive) as reader:
            essential = reader.get_essential_metadata()
            delivered = sum(mzs.size for _, mzs, _ in reader.iter_spectra())

        assert essential.total_peaks == delivered == 12


class TestIndexTolerance:
    """The index vocabulary is parsed as tolerantly as the reference parser."""

    @pytest.mark.parametrize("spelling", ["data_arrays", "data arrays", "DATA ARRAYS"])
    def test_data_kind_spellings_all_resolve(self, tmp_path, spelling):
        """Underscore, space and case variants name the same role.

        ``DataKind::from_str`` lowercases and trims before matching and
        carries ``data arrays`` as an explicit alias, so a reader that
        compares the raw string rejects valid archives.
        """
        archive = build_mzpeak(
            tmp_path / f"kind_{abs(hash(spelling))}.mzpeak",
            grid_spectra(2, 1),
            data_kind=spelling,
        )

        with MzPeakReader(archive) as reader:
            assert len(list(reader.iter_spectra())) == 2

    def test_metadata_mapping_alias_is_accepted(self, tmp_path):
        """``metadata_mapping`` is a serde alias for ``column_mapping``."""
        archive = build_mzpeak(
            tmp_path / "alias.mzpeak",
            grid_spectra(2, 1),
            column_mapping_key="metadata_mapping",
        )

        with MzPeakReader(archive) as reader:
            assert len(list(reader.iter_spectra())) == 2

    def test_absent_column_mapping_falls_back_to_conventional_names(self, tmp_path):
        """Bindings are ``serde(default)``, so an archive may omit them."""
        archive = build_mzpeak(
            tmp_path / "nomapping.mzpeak",
            grid_spectra(2, 1),
            column_mapping_key=None,
        )

        with MzPeakReader(archive) as reader:
            assert [c for c, _, _ in reader.iter_spectra()] == [
                (0, 0, 0),
                (1, 0, 0),
            ]


class TestRefusals:
    """Layouts and acquisitions the reader will not pretend to handle."""

    def test_chunked_layout_is_refused_by_name(self, tmp_path):
        """The chunked encoding is a different layout, not a variant."""
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak", grid_spectra(2, 1), layout="chunk"
        )

        with pytest.raises(NotImplementedError, match="chunked layout"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_non_imaging_archive_is_refused(self, tmp_path):
        """No positions means no pixels, and Thyra is MSI-only.

        A valid mzPeak file can carry no positions at all -- the reference
        converter only emits them when it happens to see imaging input.
        """
        archive = build_mzpeak(
            tmp_path / "nonimaging.mzpeak",
            grid_spectra(2, 1),
            include_positions=False,
        )

        with pytest.raises(ConversionRefused, match="not an imaging mzPeak archive"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_region_map_is_none(self, simple_archive):
        """mzPeak carries no region or ROI identity of any kind."""
        with MzPeakReader(simple_archive) as reader:
            assert reader.get_region_map() is None
            assert reader.get_region_info() is None


def _peak_lists(spectra, keep=3):
    """A shorter list per spectrum, on m/z values the profile never holds."""
    return [
        Spectrum(s.x, s.y, s.mzs[:keep] + 0.125, s.intensities[:keep] * 2.0)
        for s in spectra
    ]


def _drop_signal_members(archive, tmp_path):
    """Rewrite an archive's index so it names no signal member at all."""
    stripped = tmp_path / f"stripped_{archive.name}"
    with (
        zipfile.ZipFile(archive) as source,
        zipfile.ZipFile(stripped, "w", compression=zipfile.ZIP_STORED) as target,
    ):
        for name in source.namelist():
            payload = source.read(name)
            if name == "mzpeak_index.json":
                index = json.loads(payload)
                index["files"] = [
                    entry
                    for entry in index["files"]
                    if entry["data_kind"] not in ("data_arrays", "peaks")
                ]
                payload = json.dumps(index).encode("utf-8")
            target.writestr(name, payload)
    return stripped


class TestSignalMember:
    """Profile data lives in one member and centroid data in another.

    The specification puts profile spectra in the ``data_arrays`` member and
    centroid spectra in the ``peaks`` member, always. The reference converter
    writes both members for every input and leaves the unused one with zero
    rows, so ``empty_peer=True`` is the shape a real archive has.
    """

    def test_profile_archive_is_read_from_the_data_member(self, tmp_path):
        """An empty peaks member beside the data changes nothing."""
        spectra = grid_spectra(3, 2, n_points=6)
        archive = build_mzpeak(tmp_path / "profile.mzpeak", spectra, empty_peer=True)

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with MzPeakReader(archive) as reader:
                kind = reader.archive.signal_kind()
                total = reader.get_total_peak_count()
                emitted = list(reader.iter_spectra())

        assert kind == "data_arrays"
        assert len(emitted) == 6
        assert total == sum(mzs.size for _, mzs, _ in emitted) == 36
        for spectrum, (_, mzs, intensities) in zip(spectra, emitted):
            np.testing.assert_array_equal(mzs, spectrum.mzs)
            np.testing.assert_allclose(intensities, spectrum.intensities)

    @pytest.mark.parametrize("empty_peer", [True, False], ids=["empty", "absent"])
    def test_centroid_archive_is_read_from_the_peaks_member(self, tmp_path, empty_peer):
        """The data member holds no rows, or is not there: read the peaks.

        With the data member empty every ``number_of_data_points`` is null.
        Read as the point count, that column failed a cast with a numpy
        RuntimeWarning and the archive was then refused for holding no m/z
        values. Warnings are errors here so the cast cannot come back.
        """
        spectra = grid_spectra(3, 2, n_points=6)
        archive = build_mzpeak(
            tmp_path / "centroid.mzpeak",
            spectra,
            signal="centroid",
            empty_peer=empty_peer,
        )

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with MzPeakReader(archive) as reader:
                kind = reader.archive.signal_kind()
                axis = reader.get_common_mass_axis()
                total = reader.get_total_peak_count()
                emitted = list(reader.iter_spectra())

        assert kind == "peaks"
        assert [coords for coords, _, _ in emitted] == [
            (0, 0, 0),
            (1, 0, 0),
            (2, 0, 0),
            (0, 1, 0),
            (1, 1, 0),
            (2, 1, 0),
        ]
        assert total == sum(mzs.size for _, mzs, _ in emitted) == 36
        np.testing.assert_array_equal(axis, spectra[0].mzs)
        for spectrum, (_, mzs, intensities) in zip(spectra, emitted):
            np.testing.assert_array_equal(mzs, spectrum.mzs)
            np.testing.assert_allclose(intensities, spectrum.intensities)

    def test_centroid_archive_split_across_row_groups(self, tmp_path):
        """The peaks member is cut into spectra the way the data member is."""
        spectra = grid_spectra(2, 2, n_points=8)
        archive = build_mzpeak(
            tmp_path / "centroid_split.mzpeak",
            spectra,
            signal="centroid",
            empty_peer=True,
            row_group_size=3,
        )

        with MzPeakReader(archive) as reader:
            emitted = list(reader.iter_spectra())

        assert len(emitted) == 4
        for spectrum, (_, mzs, _) in zip(spectra, emitted):
            np.testing.assert_array_equal(mzs, spectrum.mzs)

    def test_both_members_filled_reads_the_data_member_only(self, tmp_path, thyra_logs):
        """Profile and picked peaks of the same spectra are not added up."""
        spectra = grid_spectra(3, 2, n_points=6)
        archive = build_mzpeak(
            tmp_path / "both.mzpeak",
            spectra,
            signal="both",
            centroids=_peak_lists(spectra),
        )

        with thyra_logs("thyra.readers.mzpeak.mzpeak_reader", "WARNING") as logs:
            with MzPeakReader(archive) as reader:
                kind = reader.archive.signal_kind()
                axis = reader.get_common_mass_axis()
                total = reader.get_total_peak_count()
                emitted = list(reader.iter_spectra())

        assert kind == "data_arrays"
        assert total == sum(mzs.size for _, mzs, _ in emitted) == 36
        # The peak lists sit 0.125 above the profile points, so one of their
        # m/z values on the axis would mean the peaks member was read too.
        np.testing.assert_array_equal(axis, spectra[0].mzs)
        for spectrum, (_, mzs, intensities) in zip(spectra, emitted):
            np.testing.assert_array_equal(mzs, spectrum.mzs)
            np.testing.assert_allclose(intensities, spectrum.intensities)

        said = [r.getMessage() for r in logs if "both members" in r.getMessage()]
        assert len(said) == 1
        assert "36 rows of profile data" in said[0]
        assert "18 rows of centroid data" in said[0]
        assert "Reading the profile member only" in said[0]

    def test_single_member_archives_say_nothing_about_both(self, tmp_path, thyra_logs):
        """The notice is for archives that fill both members, and only them."""
        spectra = grid_spectra(2, 1)
        archives = [
            build_mzpeak(tmp_path / "p.mzpeak", spectra, empty_peer=True),
            build_mzpeak(
                tmp_path / "c.mzpeak", spectra, signal="centroid", empty_peer=True
            ),
        ]

        with thyra_logs("thyra.readers.mzpeak.mzpeak_reader", "INFO") as logs:
            for archive in archives:
                with MzPeakReader(archive) as reader:
                    list(reader.iter_spectra())

        assert not [r for r in logs if "both members" in r.getMessage()]

    def test_spectrum_held_only_in_the_peaks_member_is_not_converted(self, tmp_path):
        """With both members filled, a null point count means no points.

        The reference converter splits a run by representation: a profile
        spectrum has a null ``number_of_peaks`` and a centroid one a null
        ``number_of_data_points``. The data member is the one read, so the
        centroid-only spectrum is not emitted, as an unacquired pixel is not.
        """
        spectra = [
            Spectrum(1, 1, [100.0, 101.0, 102.0], [1.0, 2.0, 3.0]),
            Spectrum(2, 1, [], []),
            Spectrum(3, 1, [100.0, 101.0], [4.0, 5.0]),
        ]
        centroids = [None, Spectrum(2, 1, [150.0, 151.0], [9.0, 9.0]), None]
        archive = build_mzpeak(
            tmp_path / "split.mzpeak", spectra, signal="both", centroids=centroids
        )

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with MzPeakReader(archive) as reader:
                counts = reader.archive.spatial_index().point_counts
                essential = reader.get_essential_metadata()
                emitted = list(reader.iter_spectra())

        np.testing.assert_array_equal(counts, [3, 0, 2])
        assert [coords for coords, _, _ in emitted] == [(0, 0, 0), (2, 0, 0)]
        assert essential.n_spectra == 3
        assert essential.dimensions == (3, 1, 1)
        assert essential.total_peaks == 5

    def test_layout_is_read_from_the_member_that_holds_the_data(self, tmp_path):
        """A chunked peaks member is refused though the data member is not."""
        archive = build_mzpeak(
            tmp_path / "chunked_peaks.mzpeak",
            grid_spectra(2, 1),
            signal="centroid",
            empty_peer=True,
            layout="chunk",
        )

        with pytest.raises(NotImplementedError, match="chunked layout"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_padding_is_counted_in_the_member_that_is_read(self, tmp_path):
        """Null pairs in the data member are not charged to a peaks read.

        Both archives hold six real points per spectrum. The profile one
        pads each spectrum with a null pair; the centroid one has no padding
        to correct for, whatever its empty data member's statistics say.
        """
        spectra = grid_spectra(2, 1, n_points=6)
        padded = build_mzpeak(
            tmp_path / "padded.mzpeak", spectra, null_pair_after=3, empty_peer=True
        )
        centroid = build_mzpeak(
            tmp_path / "centroid.mzpeak", spectra, signal="centroid", empty_peer=True
        )

        with MzPeakReader(padded) as reader:
            assert reader.archive.null_count() == 4
            assert reader.get_total_peak_count() == 12
        with MzPeakReader(centroid) as reader:
            assert not reader.archive.null_count()
            assert reader.get_total_peak_count() == 12

    def test_both_members_filled_corrects_for_the_data_member_padding(self, tmp_path):
        """The padding correction follows the data member when both hold rows."""
        spectra = grid_spectra(2, 1, n_points=6)
        archive = build_mzpeak(
            tmp_path / "both_padded.mzpeak",
            spectra,
            signal="both",
            centroids=_peak_lists(spectra),
            null_pair_after=3,
        )

        with MzPeakReader(archive) as reader:
            total = reader.get_total_peak_count()
            delivered = sum(mzs.size for _, mzs, _ in reader.iter_spectra())

        assert total == delivered == 12

    def test_missing_count_column_is_refused_by_name(self, tmp_path):
        """The column that sizes the spectra is named when it is not there."""
        archive = build_mzpeak(
            tmp_path / "uncounted.mzpeak",
            grid_spectra(2, 1),
            signal="centroid",
            empty_peer=True,
            declare_counts=False,
        )

        with pytest.raises(ConversionRefused, match="no 'number_of_peaks' column"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_two_empty_members_are_refused_for_holding_no_signal(self, tmp_path):
        """No rows in either member is an empty archive, said plainly."""
        spectra = [Spectrum(1, 1, [], []), Spectrum(2, 1, [], [])]
        archive = build_mzpeak(
            tmp_path / "empty.mzpeak",
            spectra,
            signal="both",
            centroids=[None, None],
        )

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with pytest.raises(ConversionRefused, match="yielded no m/z values"):
                with MzPeakReader(archive) as reader:
                    reader.get_common_mass_axis()

    def test_archive_without_a_signal_member_is_refused(self, tmp_path):
        """Neither member in the index: the refusal lists what is there."""
        whole = build_mzpeak(tmp_path / "whole.mzpeak", grid_spectra(2, 1))
        archive = _drop_signal_members(whole, tmp_path)

        with pytest.raises(ConversionRefused, match="no 'spectrum/data_arrays'"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()
