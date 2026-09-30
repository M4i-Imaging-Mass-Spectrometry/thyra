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

from tests.fixtures.mzpeak_builder import (
    DELTA,
    GRID,
    NUMPRESS_LINEAR,
    PLAIN,
    POSITION_X_COLUMN,
    POSITION_Y_COLUMN,
    SCANS_MEMBER,
    SQUARE_ROOT_GRID,
    Spectrum,
    build_mzpeak,
    grid_spectra,
)
from thyra.errors import ConversionRefused
from thyra.readers.mzpeak import MzPeakReader, mzpeak_reader

#: Every chunk encoding the reader decodes.
ENCODINGS = [PLAIN, DELTA, NUMPRESS_LINEAR, GRID]


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


def _calibration_scan(n_points=8):
    """A spectrum with no position: a scan that belongs to no pixel."""
    mzs = 100.0 + np.arange(n_points, dtype=np.float64) * 0.5
    return Spectrum(None, None, mzs, np.full(n_points, 999.0))


def _with_scans(archive, tmp_path, rows):
    """Rewrite an archive's scans member from ``(source_index, x, y)`` rows.

    The builder writes one scan per spectrum. A spectrum with several scans
    has to be written here.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    table = pa.table(
        {
            "source_index": pa.array([r[0] for r in rows], type=pa.uint64()),
            "scan_index": pa.array(range(len(rows)), type=pa.uint64()),
            "scan_start_time": pa.array([0.0] * len(rows), type=pa.float32()),
            POSITION_X_COLUMN: pa.array([r[1] for r in rows], type=pa.uint32()),
            POSITION_Y_COLUMN: pa.array([r[2] for r in rows], type=pa.uint32()),
        }
    )
    sink = pa.BufferOutputStream()
    pq.write_table(table, sink)
    rewritten = tmp_path / f"scans_{archive.name}"
    with (
        zipfile.ZipFile(archive) as source,
        zipfile.ZipFile(rewritten, "w", compression=zipfile.ZIP_STORED) as target,
    ):
        for name in source.namelist():
            payload = source.read(name)
            if name == SCANS_MEMBER:
                payload = sink.getvalue().to_pybytes()
            target.writestr(name, payload)
    return rewritten


class TestScansWithoutAPosition:
    """A scan on no pixel is left out, and must not move the grid.

    The imaging profile gives a scan that belongs to no pixel, a calibration
    scan for instance, a null in both position columns. Read as an integer
    a null is the smallest one there is, and the grid's origin went there.
    """

    LAYOUTS = [("point", DELTA)] + [("chunk", encoding) for encoding in ENCODINGS]

    @pytest.mark.parametrize(("layout", "encoding"), LAYOUTS)
    def test_the_scan_is_left_out_and_the_others_are_read(
        self, tmp_path, layout, encoding
    ):
        """The positioned spectra arrive as they do without the scan."""
        placed = grid_spectra(3, 2)
        twin = build_mzpeak(
            tmp_path / "twin.mzpeak", placed, layout=layout, chunk_encoding=encoding
        )
        archive = build_mzpeak(
            tmp_path / "scan.mzpeak",
            placed[:2] + [_calibration_scan()] + placed[2:],
            layout=layout,
            chunk_encoding=encoding,
        )

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with MzPeakReader(archive) as reader:
                emitted = list(reader.iter_spectra())
        with MzPeakReader(twin) as reader:
            expected = list(reader.iter_spectra())

        assert len(emitted) == len(expected) == 6
        for (coords, mzs, intensities), (coords_e, mzs_e, intensities_e) in zip(
            emitted, expected
        ):
            assert coords == coords_e
            np.testing.assert_array_equal(mzs, mzs_e)
            np.testing.assert_array_equal(intensities, intensities_e)

    def test_the_grid_is_the_one_of_the_positioned_spectra(self, tmp_path):
        """Origin, extent and count ignore the scan on no pixel."""
        placed = [
            Spectrum(5, 9, [100.0, 101.0], [1.0, 2.0]),
            Spectrum(6, 9, [100.0, 101.0], [3.0, 4.0]),
            Spectrum(5, 10, [100.0, 101.0], [5.0, 6.0]),
        ]
        archive = build_mzpeak(
            tmp_path / "offset.mzpeak", [_calibration_scan(2)] + placed
        )

        with MzPeakReader(archive) as reader:
            essential = reader.get_essential_metadata()

        assert essential.coordinate_offsets == (5, 9, 0)
        assert essential.dimensions == (2, 2, 1)
        assert essential.n_spectra == 3

    def test_the_acquisition_order_keeps_the_gap(self, tmp_path):
        """The scan keeps its place in the spectrum list; the rest keep theirs."""
        placed = grid_spectra(2, 1)
        archive = build_mzpeak(
            tmp_path / "gap.mzpeak", [placed[0], _calibration_scan(), placed[1]]
        )

        with MzPeakReader(archive) as reader:
            ordered = [
                (coords, order)
                for coords, order, _, _ in reader.iter_spectra_with_acquisition_order()
            ]

        assert ordered == [((0, 0, 0), 0), ((1, 0, 0), 2)]

    def test_the_log_counts_what_was_left_out(self, tmp_path, thyra_logs):
        """One line says how many spectra are on no pixel."""
        placed = grid_spectra(2, 2)
        archive = build_mzpeak(
            tmp_path / "two.mzpeak",
            [_calibration_scan()] + placed + [_calibration_scan()],
        )

        with thyra_logs("thyra.readers.mzpeak.mzpeak_reader", "WARNING") as logs:
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

        assert any(
            "2 of 6 spectra" in message and "have no position" in message
            for message in logs.messages
        )

    @pytest.mark.parametrize(("x", "y"), [(3, None), (None, 3)])
    def test_one_position_without_the_other_is_refused(self, tmp_path, x, y):
        """Half a position places a scan nowhere, and is not guessed at."""
        mzs = [100.0, 101.0]
        archive = build_mzpeak(
            tmp_path / "half.mzpeak",
            grid_spectra(2, 1) + [Spectrum(x, y, mzs, [1.0, 2.0])],
        )

        with pytest.raises(ConversionRefused, match="one position but not the other"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_an_archive_with_no_scan_on_a_pixel_is_refused(self, tmp_path):
        """Position columns that hold nothing give no image."""
        archive = build_mzpeak(
            tmp_path / "none.mzpeak", [_calibration_scan(), _calibration_scan()]
        )

        with pytest.raises(ConversionRefused, match="contains no positioned spectra"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    @pytest.mark.parametrize(
        ("values", "arrow_type", "expected", "present"),
        [
            ([1, None, 3], "uint32", [1, 0, 3], [True, False, True]),
            ([None, None], "null", [0, 0], [False, False]),
            ([1.0, float("nan"), None], "float64", [1, 0, 0], [True, False, False]),
            ([], "uint32", [], []),
        ],
    )
    def test_a_position_column_is_read_with_its_gaps(
        self, values, arrow_type, expected, present
    ):
        """Nulls, a column of the null type and a float NaN all mean no position."""
        import pyarrow as pa

        column = pa.chunked_array([pa.array(values, type=getattr(pa, arrow_type)())])

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            positions, has_position = mzpeak_reader._positions(column)

        assert positions.tolist() == expected
        assert positions.dtype == np.int64
        assert has_position.tolist() == present

    def test_a_spectrum_is_placed_by_its_scan_that_has_a_position(self, tmp_path):
        """A first scan on no pixel does not unplace the spectrum."""
        built = build_mzpeak(tmp_path / "multi.mzpeak", grid_spectra(2, 1))
        archive = _with_scans(
            built, tmp_path, [(0, None, None), (0, 1, 1), (1, 2, 1), (1, None, None)]
        )

        with MzPeakReader(archive) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]

        assert coords == [(0, 0, 0), (1, 0, 0)]


class TestMassAxis:
    """The reader's view of the m/z axis."""

    def test_never_claims_a_shared_axis(self, simple_archive):
        """The point layout has an axis per spectrum."""
        with MzPeakReader(simple_archive) as reader:
            assert reader.has_shared_mass_axis is False

    def test_one_grid_for_every_spectrum_is_not_a_shared_axis(self, tmp_path):
        """A grid says where a point may lie, not which points a pixel has.

        Every row of this archive carries the same grid model, and every
        spectrum the same m/z. The reader still reports no shared axis:
        it would take a pass over every index to know the second part.
        """
        archive = build_mzpeak(
            tmp_path / "grid.mzpeak",
            grid_spectra(2, 2),
            layout="chunk",
            chunk_encoding=GRID,
        )

        with MzPeakReader(archive) as reader:
            assert reader.has_shared_mass_axis is False
            axis = reader.get_common_mass_axis()

        np.testing.assert_array_equal(axis, 100.0 + 0.5 * np.arange(8))

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

    def test_chunk_encoding_that_is_not_decoded_is_refused_by_name(self, tmp_path):
        """The refusal gives the term, and comes before any spectrum."""
        archive = build_mzpeak(
            tmp_path / "slof.mzpeak",
            grid_spectra(2, 1),
            layout="chunk",
            chunk_encoding="MS:1002314",
        )

        with pytest.raises(ConversionRefused) as refusal:
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

        message = str(refusal.value)
        assert "MS:1002314 (MS-Numpress short logged float compression)" in message
        assert "slof.mzpeak" in message

    def test_one_refused_row_among_many_is_enough(self, tmp_path):
        """Every row is checked at the start, not the first batch alone."""
        spectra = grid_spectra(3, 3, n_points=8)
        rows = 9 * 2
        archive = build_mzpeak(
            tmp_path / "late.mzpeak",
            spectra,
            layout="chunk",
            chunk_encoding=[DELTA] * (rows - 1) + ["MS:4000000"],
        )

        with pytest.raises(ConversionRefused, match="chunk encoding MS:4000000,"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_grid_model_that_is_not_decoded_is_refused_by_name(self, tmp_path):
        """The open models are decoded; a vendor's model is named and left."""
        archive = build_mzpeak(
            tmp_path / "vendor_grid.mzpeak",
            grid_spectra(2, 1),
            layout="chunk",
            chunk_encoding=GRID,
            grid_type="MS:9999002",
        )

        with pytest.raises(ConversionRefused, match=r"grid of type MS:9999002 \("):
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
        """A chunked peaks member beside an empty point-layout data member.

        The converter writes this pair when it is asked for the point layout
        and finds centroid m/z on a lattice: the peaks member is chunked all
        the same.
        """
        spectra = grid_spectra(2, 1)
        archive = build_mzpeak(
            tmp_path / "chunked_peaks.mzpeak",
            spectra,
            signal="centroid",
            empty_peer=True,
            layout="chunk",
            chunk_encoding=GRID,
        )

        with MzPeakReader(archive) as reader:
            assert reader.archive.layout() == "chunk"
            assert reader.get_essential_metadata().spectrum_type == "centroid spectrum"
            emitted = list(reader.iter_spectra())

        assert [coords for coords, _, _ in emitted] == [(0, 0, 0), (1, 0, 0)]
        for (_, mzs, intensities), spectrum in zip(emitted, spectra):
            np.testing.assert_array_equal(mzs, spectrum.mzs)
            np.testing.assert_array_equal(intensities, spectrum.intensities)

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

    @pytest.mark.parametrize("layout", ["point", "chunk"])
    def test_two_empty_members_are_refused_for_holding_no_signal(
        self, tmp_path, layout
    ):
        """No rows in either member is an empty archive, said plainly."""
        spectra = [Spectrum(1, 1, [], []), Spectrum(2, 1, [], [])]
        archive = build_mzpeak(
            tmp_path / "empty.mzpeak",
            spectra,
            signal="both",
            centroids=[None, None],
            layout=layout,
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


def _read(archive):
    """Everything a conversion takes from a reader, in one pass each."""
    with MzPeakReader(archive) as reader:
        essential = reader.get_essential_metadata()
        return {
            "axis": reader.get_common_mass_axis(),
            "spectra": list(reader.iter_spectra_with_acquisition_order()),
            "n_spectra": essential.n_spectra,
            "total_peaks": essential.total_peaks,
            "mass_range": essential.mass_range,
            "spectrum_type": essential.spectrum_type,
            "dimensions": essential.dimensions,
        }


def _assert_same_reading(candidate, reference):
    """Two archives gave the reader the same data, to the bit."""
    for key in ("n_spectra", "total_peaks", "mass_range", "spectrum_type"):
        assert candidate[key] == reference[key], key
    assert candidate["dimensions"] == reference["dimensions"]
    np.testing.assert_array_equal(candidate["axis"], reference["axis"])
    assert len(candidate["spectra"]) == len(reference["spectra"])
    for ours, theirs in zip(candidate["spectra"], reference["spectra"]):
        assert ours[0] == theirs[0]
        assert ours[1] == theirs[1]
        np.testing.assert_array_equal(ours[2], theirs[2])
        np.testing.assert_array_equal(ours[3], theirs[3])
        assert ours[2].dtype == np.float64
        assert ours[3].dtype == np.float64


class TestChunkedLayout:
    """A chunked archive reads as its point twin does.

    The twin is the same spectra written in the point layout, so whatever
    the point layout is known to give, the chunked one is held to.
    """

    @pytest.mark.parametrize("encoding", ENCODINGS)
    @pytest.mark.parametrize("chunk_points", [1, 3, 100])
    def test_reads_as_the_point_twin(self, tmp_path, encoding, chunk_points):
        """Pixels, order, m/z, intensities, axis and counts all agree."""
        spectra = grid_spectra(3, 2, n_points=7, skip=[(2, 1)])
        twin = build_mzpeak(tmp_path / "point.mzpeak", spectra)
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak",
            spectra,
            layout="chunk",
            chunk_encoding=encoding,
            chunk_points=chunk_points,
        )

        _assert_same_reading(_read(archive), _read(twin))

    def test_square_root_grid_reads_as_the_point_twin(self, tmp_path):
        """The second open grid model."""
        spectra = [
            Spectrum(x, 1, [(10.0 + k / 4.0) ** 2 for k in range(x, x + 6)], range(6))
            for x in (1, 2, 3)
        ]
        twin = build_mzpeak(tmp_path / "point.mzpeak", spectra)
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak",
            spectra,
            layout="chunk",
            chunk_encoding=GRID,
            grid_type=SQUARE_ROOT_GRID,
        )

        _assert_same_reading(_read(archive), _read(twin))

    def test_mixed_encodings_read_as_the_point_twin(self, tmp_path):
        """A member may change encoding from one row to the next."""
        spectra = grid_spectra(3, 2, n_points=9)
        twin = build_mzpeak(tmp_path / "point.mzpeak", spectra)
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak",
            spectra,
            layout="chunk",
            chunk_encoding=ENCODINGS,
            chunk_points=2,
        )

        with MzPeakReader(archive) as reader:
            assert reader.archive.chunk_encodings() == sorted(ENCODINGS)
        _assert_same_reading(_read(archive), _read(twin))

    @pytest.mark.parametrize("encoding", ENCODINGS)
    @pytest.mark.parametrize("signal", ["centroid", "both"])
    def test_member_choice_holds_for_chunked_members(self, tmp_path, encoding, signal):
        """Peaks when the data member is empty; the data member otherwise."""
        spectra = grid_spectra(2, 2, n_points=6)
        centroids = [Spectrum(s.x, s.y, s.mzs[:2], s.intensities[:2]) for s in spectra]
        options = {
            "signal": signal,
            "empty_peer": signal == "centroid",
            "centroids": centroids if signal == "both" else None,
        }
        twin = build_mzpeak(tmp_path / "point.mzpeak", spectra, **options)
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak",
            spectra,
            layout="chunk",
            chunk_encoding=encoding,
            **options,
        )

        _assert_same_reading(_read(archive), _read(twin))

    @pytest.mark.parametrize("encoding", ENCODINGS)
    @pytest.mark.parametrize("rows", [1, 2, 5])
    def test_batch_size_does_not_change_the_result(
        self, tmp_path, monkeypatch, encoding, rows
    ):
        """Spectra cut by a batch boundary are joined, not doubled."""
        spectra = grid_spectra(3, 2, n_points=8)
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak",
            spectra,
            layout="chunk",
            chunk_encoding=encoding,
            chunk_points=3,
        )
        whole = _read(archive)

        monkeypatch.setattr(mzpeak_reader, "CHUNK_BATCH_ROWS", (rows, rows))

        _assert_same_reading(_read(archive), whole)

    @pytest.mark.parametrize("row_group_size", [1, 2, 7])
    def test_row_group_size_does_not_change_the_result(self, tmp_path, row_group_size):
        """Nor are those cut by a row group of the archive."""
        spectra = grid_spectra(3, 2, n_points=8)
        whole = build_mzpeak(
            tmp_path / "whole.mzpeak", spectra, layout="chunk", chunk_points=3
        )
        split = build_mzpeak(
            tmp_path / "split.mzpeak",
            spectra,
            layout="chunk",
            chunk_points=3,
            row_group_size=row_group_size,
        )

        _assert_same_reading(_read(split), _read(whole))

    def test_batches_are_sized_by_points(self, tmp_path):
        """Long rows make small batches; the bounds hold either way."""
        spectra = grid_spectra(2, 2, n_points=8)
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak", spectra, layout="chunk", chunk_points=4
        )

        with MzPeakReader(archive) as reader:
            data = reader.archive.signal()
            low, high = mzpeak_reader.CHUNK_BATCH_ROWS
            assert reader._chunk_rows_per_batch(data) == high
            reader._point_counts = reader._point_counts * 10**9
            assert reader._chunk_rows_per_batch(data) == low

    @pytest.mark.parametrize("encoding", [PLAIN, DELTA, NUMPRESS_LINEAR])
    @pytest.mark.parametrize("chunk_points", [1, 2, 3, 4, 5, 100])
    def test_padding_reads_as_the_point_twin(self, tmp_path, encoding, chunk_points):
        """Null pairs are dropped and counted as the point layout's are."""
        spectra = grid_spectra(2, 2, n_points=6)
        twin = build_mzpeak(tmp_path / "point.mzpeak", spectra, null_pair_after=3)
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak",
            spectra,
            null_pair_after=3,
            layout="chunk",
            chunk_encoding=encoding,
            chunk_points=chunk_points,
        )

        with MzPeakReader(twin) as reader:
            list(reader.iter_spectra())
            expected = reader._dropped_points
            assert reader.archive.null_count() == expected == 8
        with MzPeakReader(archive) as reader:
            list(reader.iter_spectra())
            assert reader._dropped_points == expected
            assert reader.archive.null_count() == expected
        _assert_same_reading(_read(archive), _read(twin))

    def test_lossy_m_z_stay_inside_the_declared_range(self, tmp_path):
        """The ends of every spectrum are exact, so the range holds them.

        The scale makes every decoded m/z a little off, as MS-Numpress does
        on real data. The declared range comes from the exact values, and
        the resampled axis is built on it.
        """
        spectra = [
            Spectrum(
                x,
                1,
                [100.1234567 + 0.3 * x, 250.7654321, 399.9876543 - 0.3 * x],
                [5, 6, 7],
            )
            for x in (1, 2, 3)
        ]
        scaled = np.concatenate([s.mzs for s in spectra]) * 1e5
        assert not np.any(scaled == np.round(scaled)), "the fixture is not lossy"
        archive = build_mzpeak(
            tmp_path / "lossy.mzpeak",
            spectra,
            layout="chunk",
            chunk_encoding=NUMPRESS_LINEAR,
            numpress_fixed_point=1e5,
        )

        with MzPeakReader(archive) as reader:
            low, high = reader.get_essential_metadata().mass_range
            emitted = list(reader.iter_spectra())

        assert low == min(s.mzs[0] for s in spectra)
        assert high == max(s.mzs[-1] for s in spectra)
        for (_, mzs, _), spectrum in zip(emitted, spectra):
            assert mzs[0] == spectrum.mzs[0]
            assert mzs[-1] == spectrum.mzs[-1]
            assert mzs[1] != spectrum.mzs[1]
            assert mzs[1] == pytest.approx(spectrum.mzs[1], abs=1e-5)
            assert low <= mzs.min() and mzs.max() <= high

    @pytest.mark.parametrize(
        ("encoding", "warned"),
        [(PLAIN, False), (DELTA, False), (NUMPRESS_LINEAR, True), (GRID, True)],
    )
    def test_the_axis_of_every_m_z_warns_on_a_lossy_archive(
        self, tmp_path, thyra_logs, encoding, warned
    ):
        """Without resampling a lossy archive gives a very long axis."""
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak",
            grid_spectra(2, 1),
            layout="chunk",
            chunk_encoding=encoding,
        )

        with MzPeakReader(archive) as reader:
            reader.get_essential_metadata()
            with thyra_logs("thyra.readers.mzpeak.mzpeak_reader", "WARNING") as logs:
                reader.get_common_mass_axis()

        messages = [record.getMessage() for record in logs]
        assert [encoding in message for message in messages] == [True] * warned
        assert all("resampling" in message for message in messages)

    def test_the_store_is_told_which_encodings_it_came_from(self, tmp_path):
        """Two of them are lossy, so the record is worth keeping."""
        spectra = grid_spectra(2, 1)
        twin = build_mzpeak(tmp_path / "point.mzpeak", spectra)
        archive = build_mzpeak(
            tmp_path / "chunked.mzpeak",
            spectra,
            layout="chunk",
            chunk_encoding=NUMPRESS_LINEAR,
        )

        with MzPeakReader(archive) as reader:
            chunked = reader.get_comprehensive_metadata().format_specific
        with MzPeakReader(twin) as reader:
            point = reader.get_comprehensive_metadata().format_specific

        assert chunked["layout"] == "chunk"
        assert chunked["chunk_encodings"] == [NUMPRESS_LINEAR]
        assert point["layout"] == "point"
        assert "chunk_encodings" not in point

    def test_intensity_filter_applies_to_chunked_archives(self, tmp_path):
        """The threshold is the reader's, whatever the layout."""
        spectra = grid_spectra(2, 1, n_points=6)
        twin = build_mzpeak(tmp_path / "point.mzpeak", spectra)
        archive = build_mzpeak(tmp_path / "chunked.mzpeak", spectra, layout="chunk")

        with MzPeakReader(twin, intensity_threshold=4.0) as reader:
            expected = list(reader.iter_spectra())
        with MzPeakReader(archive, intensity_threshold=4.0) as reader:
            emitted = list(reader.iter_spectra())

        assert len(emitted) == len(expected) == 2
        for ours, theirs in zip(emitted, expected):
            np.testing.assert_array_equal(ours[1], theirs[1])
            np.testing.assert_array_equal(ours[2], theirs[2])
