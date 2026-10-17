"""The sibling tables fed from the summed table's own passes (design decision D5).

A stand-in source hands its frames over as records, the way a Bruker TDF
does, and also answers the three iterators the standalone passes use.
Converted both ways on the streaming route, the stores must be identical
in every table -- and the fused way must not have touched the iterators
the sinks used to read from, which is the whole saving.
"""

import logging
from pathlib import Path
from typing import Generator, List, Optional, Tuple

import numpy as np
import pytest
from numpy.typing import NDArray

pytest.importorskip("spatialdata")

import spatialdata  # noqa: E402

from thyra.core.base_extractor import MetadataExtractor  # noqa: E402
from thyra.core.base_reader import BaseMSIReader  # noqa: E402
from thyra.core.mobility import MobilityAxis  # noqa: E402
from thyra.core.msms import FragmentationSchedule, IsolationWindow  # noqa: E402
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata  # noqa: E402

PIXELS = [(0, 0), (1, 0), (0, 1), (1, 1)]
#: A decreasing 1/K0 ramp of 8 scans, as a TDF's is; window 0 owns scans
#: 0..4 (1/K0 1.5 .. 1.27), window 1 scans 4..8 (1.21 .. 1.1).
SCAN_K0 = np.linspace(1.5, 1.1, 8)
WINDOWS = (
    IsolationWindow(313.275, 0.5, 0.5, 34.3, scan_begin=0, scan_end=4),
    IsolationWindow(936.578, 0.5, 0.5, 49.6, scan_begin=4, scan_end=8),
)
SCHEDULE = FragmentationSchedule(ms_level=2, windows=WINDOWS, source="stub")

#: One point cloud per pixel: (m/z, scan, intensity). Two points at one
#: m/z on different scans, one m/z in both windows, one pixel sparse.
CLOUDS = {
    0: [(100.0, 1, 10.0), (100.0, 5, 4.0), (300.0, 2, 20.0), (200.0, 6, 30.0)],
    1: [(300.0, 0, 40.0), (200.0, 7, 50.0), (500.0, 5, 60.0)],
    2: [(100.0, 3, 5.0), (400.0, 4, 7.0)],
    3: [(500.0, 2, 1.0), (100.0, 6, 2.0), (500.0, 7, 3.0)],
}


def _summed(p: int) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    totals: dict = {}
    for mz, _scan, intensity in CLOUDS[p]:
        totals[mz] = totals.get(mz, 0.0) + intensity
    mzs = np.array(sorted(totals), dtype=np.float64)
    return mzs, np.array([totals[m] for m in mzs], dtype=np.float64)


def _cloud(p: int):
    points = sorted(CLOUDS[p], key=lambda t: (t[1], t[0]))
    mzs = np.array([t[0] for t in points], dtype=np.float64)
    mobility = np.array([SCAN_K0[t[1]] for t in points], dtype=np.float64)
    intensities = np.array([t[2] for t in points], dtype=np.float64)
    return mzs, mobility, intensities


def _precursors(p: int) -> List[Tuple[int, NDArray[np.float64], NDArray[np.float64]]]:
    out = []
    for w, window in enumerate(WINDOWS):
        totals: dict = {}
        for mz, scan, intensity in CLOUDS[p]:
            if window.scan_begin <= scan < window.scan_end:
                totals[mz] = totals.get(mz, 0.0) + intensity
        if totals:
            mzs = np.array(sorted(totals), dtype=np.float64)
            out.append((w, mzs, np.array([totals[m] for m in mzs], dtype=np.float64)))
    return out


class _Frame:
    def __init__(self, p: int, coords):
        self.coords = coords
        self._p = p

    def spectrum(self):
        return _summed(self._p)

    def mobility_points(self):
        return _cloud(self._p)

    def precursor_spectra(self):
        return _precursors(self._p)


class _IndexedFrame(_Frame):
    """A record that offers the indexed view, as the Bruker TDF record does.

    Its flat view raises, so a sink reached through it would fail the
    conversion: whatever the store holds came through the indexed view.
    """

    def mobility_points(self):
        raise AssertionError("the indexed view was offered and not taken")

    def mobility_points_indexed(self):
        mzs, mobility, intensities = _cloud(self._p)
        unique, inverse = np.unique(mzs, return_inverse=True)
        return unique, np.asarray(inverse).ravel(), mobility, intensities


class _StubExtractor(MetadataExtractor):
    def __init__(self):
        super().__init__(data_source=None)

    def _extract_essential_impl(self) -> EssentialMetadata:
        return EssentialMetadata(
            dimensions=(2, 2, 1),
            coordinate_bounds=(0.0, 1.0, 0.0, 1.0),
            mass_range=(100.0, 500.0),
            pixel_size=(10.0, 10.0),
            n_spectra=4,
            total_peaks=12,
            source_path="stub_fused",
            spectrum_type="centroid spectrum",
        )

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        return ComprehensiveMetadata(
            essential=self._extract_essential_impl(),
            format_specific={"format": "stub"},
            acquisition_params={},
            instrument_info={},
            raw_metadata={},
        )


class FusedStubReader(BaseMSIReader):
    """A TIMS PASEF source in miniature that also hands its frames over as records."""

    frame_scans = True

    def __init__(self):
        super().__init__(Path("stub_fused"))
        self.mobility_passes = 0
        self.precursor_passes = 0
        self.frame_passes = 0

    def _create_metadata_extractor(self) -> MetadataExtractor:
        return _StubExtractor()

    @property
    def has_shared_mass_axis(self) -> bool:
        return False

    def get_common_mass_axis(self) -> NDArray[np.float64]:
        return np.array([100.0, 200.0, 300.0, 400.0, 500.0])

    def iter_spectra(self) -> Generator:
        for p, (x, y) in enumerate(PIXELS):
            mzs, intensities = _summed(p)
            yield (x, y, 0), mzs, intensities

    # -- mobility ----------------------------------------------------------

    @property
    def has_ion_mobility(self) -> bool:
        return True

    @property
    def has_shared_mobility_axis(self) -> bool:
        return False

    def get_mobility_axis(self) -> Optional[MobilityAxis]:
        return MobilityAxis(
            kind_accession="MS:1002815",
            kind_name="inverse reduced ion mobility",
            unit_accession="MS:1002814",
            unit_name="volt-second per square centimeter",
            values=SCAN_K0.copy(),
            acq_range=(1.1, 1.5),
            source="stub",
        )

    def iter_mobility_spectra(self) -> Generator:
        self.mobility_passes += 1
        for p, (x, y) in enumerate(PIXELS):
            mzs, mobility, intensities = _cloud(p)
            yield (x, y, 0), mzs, mobility, intensities

    # -- MS/MS -------------------------------------------------------------

    def get_fragmentation(self) -> Optional[FragmentationSchedule]:
        return SCHEDULE

    @property
    def has_precursor_spectra(self) -> bool:
        return True

    def iter_precursor_spectra(self) -> Generator:
        self.precursor_passes += 1
        for p, (x, y) in enumerate(PIXELS):
            for window, mzs, intensities in _precursors(p):
                yield (x, y, 0), window, mzs, intensities

    # -- the records --------------------------------------------------------

    @property
    def has_frame_scans(self) -> bool:
        return self.frame_scans

    def iter_frame_scans(self) -> Generator:
        self.frame_passes += 1
        for p, (x, y) in enumerate(PIXELS):
            yield _Frame(p, (x, y, 0))

    def close(self) -> None:
        pass


class UnfusedStubReader(FusedStubReader):
    frame_scans = False


#: When each pixel of ``PIXELS`` was acquired: the second row right to
#: left, numbered as Bruker numbers frames, with a gap where another
#: region's frame would be (design decision D26).
ORDER = {0: 1, 1: 2, 2: 5, 3: 4}


class _OrderedFrame(_Frame):
    @property
    def acquisition_order(self) -> int:
        return ORDER[self._p]


class OrderedStubReader(FusedStubReader):
    """The same source, which knows when each pixel was acquired."""

    @property
    def has_acquisition_order(self) -> bool:
        return True

    def iter_spectra_with_acquisition_order(self) -> Generator:
        for p, (x, y) in enumerate(PIXELS):
            mzs, intensities = _summed(p)
            yield (x, y, 0), ORDER[p], mzs, intensities

    def iter_frame_scans(self) -> Generator:
        self.frame_passes += 1
        for p, (x, y) in enumerate(PIXELS):
            yield _OrderedFrame(p, (x, y, 0))


class UnfusedOrderedStubReader(OrderedStubReader):
    frame_scans = False


#: The four spectra, then pixel 0 measured a second time, last and so out
#: of order: every entry of the repeat lands on a (row, column) the first
#: measurement already holds, in every table.
REMEASURED = list(enumerate(PIXELS)) + [(0, PIXELS[0])]


class RemeasuredStubReader(FusedStubReader):
    """The same source with a fifth frame repeating pixel 0 on its own position."""

    def iter_spectra(self) -> Generator:
        for p, (x, y) in REMEASURED:
            mzs, intensities = _summed(p)
            yield (x, y, 0), mzs, intensities

    def iter_mobility_spectra(self) -> Generator:
        self.mobility_passes += 1
        for p, (x, y) in REMEASURED:
            mzs, mobility, intensities = _cloud(p)
            yield (x, y, 0), mzs, mobility, intensities

    def iter_precursor_spectra(self) -> Generator:
        self.precursor_passes += 1
        for p, (x, y) in REMEASURED:
            for window, mzs, intensities in _precursors(p):
                yield (x, y, 0), window, mzs, intensities

    def iter_frame_scans(self) -> Generator:
        self.frame_passes += 1
        for p, (x, y) in REMEASURED:
            yield _Frame(p, (x, y, 0))


class UnfusedRemeasuredStubReader(RemeasuredStubReader):
    frame_scans = False


class _SharedAxis:
    """Every pixel's points listed as one ``(m/z, mobility)`` feature axis."""

    @property
    def has_shared_mobility_axis(self) -> bool:
        return True

    def get_shared_mobility_features(self):
        pairs = sorted(
            {(mz, SCAN_K0[scan]) for cloud in CLOUDS.values() for mz, scan, _ in cloud}
        )
        return (
            np.array([mz for mz, _ in pairs], dtype=np.float64),
            np.array([k0 for _, k0 in pairs], dtype=np.float64),
        )


class SharedAxisStubReader(_SharedAxis, UnfusedStubReader):
    pass


class SharedAxisRemeasuredStubReader(_SharedAxis, UnfusedRemeasuredStubReader):
    pass


class _OutOfRangeRepeat(SharedAxisStubReader):
    """Pixel 0 measured again, its summed spectrum wholly outside the mass range.

    The summed table never sees the repeat -- the spectrum resamples to
    nothing, so the position is measured once as far as it knows -- but
    the shared feature axis keeps every listed feature whatever the range,
    so the mobility table still reads pixel 0's cloud twice.
    """

    def iter_spectra(self) -> Generator:
        yield from super().iter_spectra()
        yield (0, 0, 0), np.array([900.0]), np.array([5.0])

    def iter_mobility_spectra(self) -> Generator:
        yield from super().iter_mobility_spectra()
        mzs, mobility, intensities = _cloud(0)
        yield (0, 0, 0), mzs, mobility, intensities


class _GridOutOfRangeRepeat(UnfusedStubReader):
    """:class:`_OutOfRangeRepeat` on the grid route (issue #495, finding B2).

    The summed table resamples the repeat to nothing and so never flags
    pixel 0 as measured twice; only the grid builder's own repeated-row
    check (D40) can sum the second cloud into row 0.
    """

    def iter_spectra(self) -> Generator:
        yield from super().iter_spectra()
        yield (0, 0, 0), np.array([900.0]), np.array([5.0])

    def iter_mobility_spectra(self) -> Generator:
        yield from super().iter_mobility_spectra()
        mzs, mobility, intensities = _cloud(0)
        yield (0, 0, 0), mzs, mobility, intensities


class IndexedStubReader(FusedStubReader):
    """The same source whose records hand their points over indexed."""

    def iter_frame_scans(self) -> Generator:
        self.frame_passes += 1
        for p, (x, y) in enumerate(PIXELS):
            yield _IndexedFrame(p, (x, y, 0))


#: The four spectra, then pixel 1's cloud measured again on pixel 0's
#: position, last and so out of order: the summed table sums it into
#: that position's row, which is one pixel, not two (D38).
REPEATED = list(enumerate(PIXELS)) + [(1, PIXELS[0])]


class RepeatedStubReader(FusedStubReader):
    """The same source with a fifth frame on a position already measured."""

    def iter_spectra(self) -> Generator:
        for p, (x, y) in REPEATED:
            mzs, intensities = _summed(p)
            yield (x, y, 0), mzs, intensities

    def iter_mobility_spectra(self) -> Generator:
        self.mobility_passes += 1
        for p, (x, y) in REPEATED:
            mzs, mobility, intensities = _cloud(p)
            yield (x, y, 0), mzs, mobility, intensities

    def iter_precursor_spectra(self) -> Generator:
        self.precursor_passes += 1
        for p, (x, y) in REPEATED:
            for window, mzs, intensities in _precursors(p):
                yield (x, y, 0), window, mzs, intensities

    def iter_frame_scans(self) -> Generator:
        self.frame_passes += 1
        for p, (x, y) in REPEATED:
            yield _Frame(p, (x, y, 0))


class UnfusedRepeatedStubReader(RepeatedStubReader):
    frame_scans = False


RESAMPLED = {
    "method": "nearest_neighbor",
    "axis_type": "constant",
    "target_bins": 43,
    "min_mz": 90.0,
    "max_mz": 510.0,
}


def _convert(reader, out: Path, resampling=RESAMPLED, **kwargs) -> Path:
    from thyra.converters.spatialdata.streaming_converter import (
        StreamingSpatialDataConverter,
    )
    from thyra.utils.windows_paths import prepare_zarr_output_path

    out = prepare_zarr_output_path(out, "stub")
    converter = StreamingSpatialDataConverter(
        reader,
        out,
        dataset_id="stub",
        pixel_size_um=10.0,
        resampling_config=resampling,
        mobility_grid=True,
        **kwargs,
    )
    assert converter.convert(), "conversion reported failure"
    return out


def _read(out: Path):
    from thyra.utils.windows_paths import prepare_zarr_read_path

    return spatialdata.read_zarr(prepare_zarr_read_path(out))


def _dense(table) -> np.ndarray:
    X = table.X
    return np.asarray(X.toarray() if hasattr(X, "toarray") else X)


def _distinct_entries(table) -> int:
    merged = table.X.copy()
    merged.sum_duplicates()
    return int(merged.nnz)


class TestARemeasuredPosition:
    """A position measured twice is one summed row in every table (#241).

    The summed table always summed it. The siblings did on the fused
    route only: the standalone passes stored both measurements' entries
    side by side, a matrix that densifies to the same values but holds
    more entries than cells and is not canonical.
    """

    @pytest.mark.parametrize(
        "plain, remeasured",
        [
            (FusedStubReader, RemeasuredStubReader),
            (UnfusedStubReader, UnfusedRemeasuredStubReader),
            (SharedAxisStubReader, SharedAxisRemeasuredStubReader),
        ],
        ids=["fused", "unfused", "shared-axis"],
    )
    def test_every_table_sums_the_repeat_into_one_entry_per_cell(
        self, tmp_path, plain, remeasured
    ):
        once = _read(_convert(plain(), tmp_path / "once.zarr")).tables
        twice = _read(_convert(remeasured(), tmp_path / "twice.zarr")).tables

        assert set(twice) == {"stub_z0", "stub_z0_mobility", "stub_z0_msms"}
        for key, table in twice.items():
            assert table.n_obs == 4, key
            assert table.X.nnz == _distinct_entries(table), key
            assert table.X.nnz == once[key].X.nnz, key
            assert list(table.var.index) == list(once[key].var.index), key
            # Pixel 0 is row 0, and its row is its spectrum twice over.
            expected = _dense(once[key])
            expected[0] *= 2
            np.testing.assert_allclose(_dense(table), expected, err_msg=key)

    def test_a_repeat_the_summed_table_never_saw_is_summed_too(self, tmp_path):
        tables = _read(_convert(_OutOfRangeRepeat(), tmp_path / "range.zarr")).tables
        once = _read(_convert(SharedAxisStubReader(), tmp_path / "once.zarr")).tables

        table = tables["stub_z0_mobility"]
        assert table.X.nnz == _distinct_entries(table)
        expected = _dense(once["stub_z0_mobility"])
        expected[0] *= 2
        np.testing.assert_allclose(_dense(table), expected)

    def test_the_grid_sums_a_repeat_the_summed_table_never_saw(self, tmp_path):
        tables = _read(
            _convert(_GridOutOfRangeRepeat(), tmp_path / "range.zarr")
        ).tables
        once = _read(_convert(UnfusedStubReader(), tmp_path / "once.zarr")).tables

        table = tables["stub_z0_mobility"]
        assert table.X.nnz == _distinct_entries(table)
        expected = _dense(once["stub_z0_mobility"])
        expected[0] *= 2
        np.testing.assert_allclose(_dense(table), expected)

    def test_both_routes_write_the_same_store(self, tmp_path):
        a = _read(_convert(RemeasuredStubReader(), tmp_path / "fused.zarr")).tables
        b = _read(
            _convert(UnfusedRemeasuredStubReader(), tmp_path / "unfused.zarr")
        ).tables

        for key in ("stub_z0_mobility", "stub_z0_msms"):
            assert a[key].X.nnz == b[key].X.nnz, key
            np.testing.assert_array_equal(_dense(a[key]), _dense(b[key]))
        for key, block in (
            ("stub_z0_mobility", "mobility_marginal"),
            ("stub_z0_msms", "demultiplexed_current"),
        ):
            assert float(a[key].uns[block]["current_ratio"]) == pytest.approx(1.0)
            assert float(b[key].uns[block]["current_ratio"]) == pytest.approx(1.0)


class TestFusedPasses:
    def test_the_fused_store_is_the_unfused_store(self, tmp_path):
        fused, unfused = FusedStubReader(), UnfusedStubReader()
        a = _read(_convert(fused, tmp_path / "fused.zarr"))
        b = _read(_convert(unfused, tmp_path / "unfused.zarr"))

        assert set(a.tables) == set(b.tables)
        assert set(a.tables) == {"stub_z0", "stub_z0_mobility", "stub_z0_msms"}
        for key in a.tables:
            ta, tb = a.tables[key], b.tables[key]
            np.testing.assert_array_equal(_dense(ta), _dense(tb))
            assert list(ta.var.index) == list(tb.var.index)
            assert list(ta.obs.index) == list(tb.obs.index)
            for column in ta.var.columns:
                np.testing.assert_array_equal(
                    ta.var[column].to_numpy(), tb.var[column].to_numpy()
                )
        heat_a = a.tables["stub_z0"].uns["mobility_heatmap"]["counts"]
        heat_b = b.tables["stub_z0"].uns["mobility_heatmap"]["counts"]
        np.testing.assert_array_equal(np.asarray(heat_a), np.asarray(heat_b))

    def test_the_fused_route_reads_the_source_twice_and_never_the_iterators(
        self, tmp_path
    ):
        fused = FusedStubReader()
        _convert(fused, tmp_path / "fused.zarr")
        # Count and scatter: one read each, serving the summed table, the
        # heatmap, the grid and the split alike.
        assert fused.frame_passes == 2
        assert fused.mobility_passes == 0
        assert fused.precursor_passes == 0

        unfused = UnfusedStubReader()
        _convert(unfused, tmp_path / "unfused.zarr")
        assert unfused.frame_passes == 0
        # Heatmap + discovery fused into one, then the grid's scatter.
        assert unfused.mobility_passes == 2
        assert unfused.precursor_passes == 2

    def test_the_indexed_view_gives_the_same_store(self, tmp_path):
        # The record offering (unique m/z, inverse) instead of the points'
        # m/z changes how many binary searches the mapping does, nothing
        # else. Its flat view raises, so the conversion completing is
        # itself the proof that the indexed one was preferred.
        a = _read(_convert(FusedStubReader(), tmp_path / "flat.zarr"))
        b = _read(_convert(IndexedStubReader(), tmp_path / "indexed.zarr"))

        assert set(a.tables) == set(b.tables)
        for key in a.tables:
            ta, tb = a.tables[key], b.tables[key]
            np.testing.assert_array_equal(_dense(ta), _dense(tb))
            assert list(ta.var.index) == list(tb.var.index)
            assert list(ta.obs.index) == list(tb.obs.index)
        np.testing.assert_array_equal(
            np.asarray(a.tables["stub_z0"].uns["mobility_heatmap"]["counts"]),
            np.asarray(b.tables["stub_z0"].uns["mobility_heatmap"]["counts"]),
        )

    def test_a_record_without_the_indexed_view_still_converts(self, tmp_path):
        # The fallback is not incidental: a record that cannot factor its
        # m/z (every source but the TDF one) offers no such method.
        assert not hasattr(_Frame(0, (0, 0, 0)), "mobility_points_indexed")
        sdata = _read(_convert(FusedStubReader(), tmp_path / "flat_only.zarr"))
        assert "mobility_heatmap" in sdata.tables["stub_z0"].uns

    def test_both_routes_give_every_table_the_acquisition_order(self, tmp_path):
        # The records carry the order on the fused route and the ordered
        # iterator on the other; the siblings copy the summed table's obs,
        # so all three tables say the same thing either way.
        a = _read(_convert(OrderedStubReader(), tmp_path / "fused.zarr"))
        b = _read(_convert(UnfusedOrderedStubReader(), tmp_path / "unfused.zarr"))

        assert set(a.tables) == {"stub_z0", "stub_z0_mobility", "stub_z0_msms"}
        for sdata in (a, b):
            for key in sdata.tables:
                column = sdata.tables[key].obs["acquisition_order"]
                assert column.dtype == np.int64
                # Rows in grid order: (0, 0), (1, 0), (0, 1), (1, 1).
                assert column.tolist() == [1, 2, 5, 4]
        for key in a.tables:
            np.testing.assert_array_equal(_dense(a.tables[key]), _dense(b.tables[key]))

    def test_a_source_without_an_order_gives_no_table_the_column(self, tmp_path):
        sdata = _read(_convert(FusedStubReader(), tmp_path / "fused.zarr"))

        for table in sdata.tables.values():
            assert "acquisition_order" not in table.obs.columns

    def test_the_marginals_are_exact_either_way(self, tmp_path):
        for reader, name in ((FusedStubReader(), "f"), (UnfusedStubReader(), "u")):
            sdata = _read(_convert(reader, tmp_path / f"{name}.zarr"))
            grid = sdata.tables["stub_z0_mobility"].uns["mobility_marginal"]
            split = sdata.tables["stub_z0_msms"].uns["demultiplexed_current"]
            assert float(grid["current_ratio"]) == pytest.approx(1.0)
            assert float(split["current_ratio"]) == pytest.approx(1.0)


class TestRepeatedPosition:
    @pytest.mark.parametrize(
        "reader_type", [RepeatedStubReader, UnfusedRepeatedStubReader]
    )
    def test_the_heatmap_divides_by_positions_not_spectra(
        self, tmp_path, thyra_logs, reader_type
    ):
        with thyra_logs("thyra.converters", logging.WARNING) as records:
            out = _convert(reader_type(), tmp_path / "repeated.zarr")
        table = _read(out).tables["stub_z0"]
        messages = [r.getMessage() for r in records]

        assert table.n_obs == 4
        counts = np.asarray(table.uns["mobility_heatmap"]["counts"])
        every_point = sum(_cloud(p)[2].sum() for p, _ in REPEATED)
        assert counts.sum() == pytest.approx(every_point / 4)
        assert float(table.uns["mobility_heatmap"]["current_ratio"]) == (
            pytest.approx(1.0)
        )
        assert not any("mass-mobility heatmap" in m for m in messages)


class _EmptySummedFrame(_Frame):
    """Pixel 3's record: raw points, but a summed spectrum with nothing in it.

    What a TDF frame under the vendor centroid gives when the picker finds
    no peak in raw scans that carry current.
    """

    def spectrum(self):
        return None if self._p == 3 else super().spectrum()


class EmptySummedStubReader(FusedStubReader):
    """Pixel 3 gets no row in the summed table; its raw cloud is still read."""

    def iter_spectra(self) -> Generator:
        for p, (x, y) in enumerate(PIXELS):
            if p != 3:
                mzs, intensities = _summed(p)
                yield (x, y, 0), mzs, intensities

    def iter_frame_scans(self) -> Generator:
        self.frame_passes += 1
        for p, (x, y) in enumerate(PIXELS):
            yield _EmptySummedFrame(p, (x, y, 0))


class UnfusedEmptySummedStubReader(EmptySummedStubReader):
    frame_scans = False


class TestAPositionWithoutARow:
    """The heatmap is a mean over the summed table's rows (issue #495, D58)."""

    @pytest.mark.parametrize(
        "reader_type", [EmptySummedStubReader, UnfusedEmptySummedStubReader]
    )
    def test_its_points_and_pixel_stay_out(self, tmp_path, reader_type):
        table = _read(_convert(reader_type(), tmp_path / "out.zarr")).tables["stub_z0"]

        assert table.n_obs == 3
        counts = np.asarray(table.uns["mobility_heatmap"]["counts"])
        rowed = sum(_cloud(p)[2].sum() for p in range(3))
        assert counts.sum() == pytest.approx(rowed / 3)
        assert float(table.uns["mobility_heatmap"]["current_ratio"]) == (
            pytest.approx(1.0)
        )

    def test_both_routes_write_the_same_heatmap(self, tmp_path):
        a = _read(_convert(EmptySummedStubReader(), tmp_path / "a.zarr")).tables
        b = _read(_convert(UnfusedEmptySummedStubReader(), tmp_path / "b.zarr")).tables

        np.testing.assert_array_equal(
            np.asarray(a["stub_z0"].uns["mobility_heatmap"]["counts"]),
            np.asarray(b["stub_z0"].uns["mobility_heatmap"]["counts"]),
        )


class _ZeroCloudFrame(_EmptySummedFrame):
    """Pixel 3 as before, its raw points all zero: it carries no signal."""

    def mobility_points(self):
        mzs, mobility, intensities = _cloud(self._p)
        if self._p == 3:
            intensities = np.zeros_like(intensities)
        return mzs, mobility, intensities


class ZeroCloudStubReader(EmptySummedStubReader):
    def iter_frame_scans(self) -> Generator:
        self.frame_passes += 1
        for p, (x, y) in enumerate(PIXELS):
            yield _ZeroCloudFrame(p, (x, y, 0))


class TestTheFusedGridCountsWhatIsLost:
    """B3 on the fused route: a rowless frame is skipped only with signal (D58)."""

    @pytest.mark.parametrize(
        "reader_type, skipped",
        [(EmptySummedStubReader, 1), (ZeroCloudStubReader, 0)],
        ids=["with-signal", "all-zero"],
    )
    def test_the_skip_warning(self, tmp_path, thyra_logs, reader_type, skipped):
        with thyra_logs("thyra.converters", logging.WARNING) as records:
            _convert(reader_type(), tmp_path / "out.zarr")
        said = [
            r.getMessage()
            for r in records
            if "mobility spectra had no row" in r.getMessage()
        ]

        assert said == (
            [f"{skipped} mobility spectra had no row in the MSI table and were skipped"]
            if skipped
            else []
        )


#: A laid non-constant axis over the range the clouds span (100 to 500).
#: Its first and last centre sit inside that range, so the clouds' m/z 100
#: and 500 points are in the half-bin margin the summed table keeps.
PHYSICS = {
    "method": "nearest_neighbor",
    "axis_type": "linear_tof",
    "target_bins": 40,
    "min_mz": 100.0,
    "max_mz": 500.0,
}


class TestPointsInTheEdgeMargin:
    """Every sibling keeps the points the summed table keeps (the kept range).

    On a constant axis the first and last centre are the range, so the
    margin is empty. On a physics axis they are half a bin inside it, and a
    sink that tested the strict axis span dropped the file's first and last
    samples: its marginals then fell short of the summed table.
    """

    @pytest.mark.parametrize(
        "reader_type",
        [FusedStubReader, IndexedStubReader, UnfusedStubReader],
        ids=["fused", "indexed", "unfused"],
    )
    def test_the_marginals_and_the_summed_table_hold_every_point(
        self, tmp_path, reader_type
    ):
        sdata = _read(_convert(reader_type(), tmp_path / "s.zarr", resampling=PHYSICS))
        summed = sdata.tables["stub_z0"]
        axis = summed.var["mz"].to_numpy()
        # The premise: both end samples are past the end centres.
        assert axis[0] > 100.0 and axis[-1] < 500.0

        every_point = sum(_cloud(p)[2].sum() for p in range(len(PIXELS)))
        assert float(_dense(summed).sum()) == pytest.approx(every_point)

        heatmap = summed.uns["mobility_heatmap"]
        assert float(heatmap["current_ratio"]) == pytest.approx(1.0)
        assert np.asarray(heatmap["counts"]).sum() * len(PIXELS) == pytest.approx(
            every_point, rel=1e-5
        )

        grid = sdata.tables["stub_z0_mobility"]
        assert float(_dense(grid).sum()) == pytest.approx(every_point)
        marginal = grid.uns["mobility_marginal"]
        assert float(marginal["current_ratio"]) == pytest.approx(1.0)
        assert float(marginal["current_ratio_pixel_min"]) == pytest.approx(1.0)
        assert float(marginal["current_ratio_pixel_max"]) == pytest.approx(1.0)

        split = sdata.tables["stub_z0_msms"].uns["demultiplexed_current"]
        assert float(split["current_ratio"]) == pytest.approx(1.0)
