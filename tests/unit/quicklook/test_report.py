"""The one-file HTML overview of a store."""

from __future__ import annotations

import base64
import io
import re
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from tests.fixtures.hand_built_store import L_SHAPE, MZ, make_table, write_store
from thyra.quicklook.report import (
    _as_photo,
    _tic_raster,
    build_report,
    file_name_only,
    pick_peaks,
    spectrum_svg,
    table_facts,
)
from thyra.quicklook.store import Positions, open_store

#: Where the source was read from, as a store records it: the whole path.
WINDOWS_SOURCE = "\\".join(["D:", "acquisitions", "jdoe", "slide 7", "run 7.imzML"])
POSIX_SOURCE = "/".join(["", "home", "jdoe", "slide 7", "run 7.imzML"])


def _metadata(source_path: str = WINDOWS_SOURCE, **analysis):
    """A metadata block with fields the page prints and fields it must not."""
    return {
        "schema_version": "0.14.0",
        "ms_analysis": {
            "pixel_size_um": {"x": 25.0, "y": 25.0},
            "instrument_model": "timsTOF fleX",
            "manufacturer": "Bruker",
            "serial_number": "SN-0042-SECRET",
            **analysis,
        },
        "sample": {"sample_name": "patient Jane Doe"},
        "provenance": {
            "thyra_version": "4.7.0",
            "source_path": source_path,
            "pixel_size_source": "automatic",
        },
    }


def _store(tmp_path: Path, images=None, attrs=None, **metadata):
    n = MZ.size
    tables = {
        "ds_z0": make_table(uns={"msi_metadata": _metadata(**metadata)}),
        "ds_z0_mobility": make_table(extra_var={"mobility": np.linspace(0.8, 1.2, n)}),
        "ds_z0_msms": make_table(
            mz=np.arange(100.0, 104.0), extra_var={"precursor_mz": np.full(4, 600.0)}
        ),
    }
    return write_store(tmp_path / "Brain 1.zarr", tables, images=images, attrs=attrs)


def _embedded_images(page: str):
    return [
        Image.open(io.BytesIO(base64.b64decode(data)))
        for data in re.findall(r'src="data:image/png;base64,([^"]+)"', page)
    ]


class TestPickPeaks:
    AXIS = np.arange(100.0, 110.0, 0.5)  # twenty bins, half an m/z apart

    def _spectrum(self, **heights):
        spectrum = np.zeros(self.AXIS.size)
        for mz, height in heights.items():
            spectrum[
                int(np.flatnonzero(self.AXIS == float(mz[1:].replace("_", ".")))[0])
            ] = height
        return spectrum

    def test_the_tallest_come_first(self):
        spectrum = self._spectrum(m101=5.0, m104=9.0, m107=7.0)

        peaks = pick_peaks(self.AXIS, spectrum, count=3, separation=1.5)

        assert [self.AXIS[peak] for peak in peaks] == [104.0, 107.0, 101.0]

    def test_a_peak_too_close_to_a_taller_one_is_skipped(self):
        """101 is the second tallest, and 1 m/z from the tallest."""
        spectrum = self._spectrum(m100=9.0, m101=8.0, m105=3.0)

        peaks = pick_peaks(self.AXIS, spectrum, count=2, separation=1.5)

        assert [self.AXIS[peak] for peak in peaks] == [100.0, 105.0]

    def test_the_shoulder_of_a_peak_is_not_a_peak(self):
        spectrum = self._spectrum(m103_5=4.0, m104=9.0, m104_5=6.0, m108=1.0)

        peaks = pick_peaks(self.AXIS, spectrum, count=5, separation=0.1)

        assert [self.AXIS[peak] for peak in peaks] == [104.0, 108.0]

    def test_no_more_than_asked_for(self):
        spectrum = self._spectrum(m100=1.0, m102=2.0, m104=3.0, m106=4.0)

        assert len(pick_peaks(self.AXIS, spectrum, count=2, separation=1.5)) == 2

    def test_a_spectrum_without_signal_has_no_peaks(self):
        assert pick_peaks(self.AXIS, np.zeros(self.AXIS.size)) == []


class TestFileNameOnly:
    @pytest.mark.parametrize("path", [WINDOWS_SOURCE, POSIX_SOURCE])
    def test_only_the_last_part_is_kept_whichever_separator(self, path):
        assert file_name_only(path) == "run 7.imzML"

    def test_a_folder_source_keeps_its_name(self):
        assert (
            file_name_only(WINDOWS_SOURCE.replace("run 7.imzML", "run.d\\")) == "run.d"
        )

    @pytest.mark.parametrize("missing", [None, "", "   ", 7])
    def test_nothing_recorded_is_none(self, missing):
        assert file_name_only(missing) is None


class TestSpectrumSvg:
    def test_a_one_bin_peak_among_many_bins_keeps_its_full_height(self):
        """Thinning by taking every n-th bin would lose it."""
        mz = np.linspace(100.0, 1000.0, 200_001)
        spectrum = np.ones(mz.size)
        spectrum[123_457] = 1000.0

        svg = spectrum_svg(mz, spectrum)

        heights = [float(y) for y in re.findall(r"L[\d.]+,([\d.]+)", svg)]
        assert len(heights) <= 852
        assert min(heights) == 15.0  # the top of the plot
        assert ">1000<" in svg

    def test_it_draws_with_the_colour_of_the_page(self):
        svg = spectrum_svg(MZ, np.arange(1.0, MZ.size + 1))

        assert "currentColor" in svg
        assert "xmlns" not in svg and "#" not in svg

    def test_a_step_with_no_bin_in_it_is_drawn_at_zero(self):
        """Two bins far apart: the gap between them carries no intensity."""
        svg = spectrum_svg(np.array([100.0, 101.0, 900.0]), np.array([5.0, 5.0, 5.0]))

        heights = [float(y) for y in re.findall(r"L[\d.]+,([\d.]+)", svg)]
        assert heights.count(15.0) == 2
        assert heights.count(260.0) == 2  # the floor, plus the closing corner


class TestTheTicRaster:
    """The stored TIC image, laid out again on the raster of the ion images."""

    #: Two slices of two by two; the last pixel of the second was not measured.
    VOLUME = [(x, y, z) for z in (0, 1) for y in (0, 1) for x in (0, 1)][:-1]

    def _positions(self, points):
        points = np.asarray(points, dtype=np.int64)
        z = points[:, 2] if points.shape[1] == 3 else None
        return Positions(points[:, 0], points[:, 1], z)

    def test_a_volume_is_read_slice_row_column(self):
        """Each voxel holds its own address, so a swapped axis cannot pass."""
        z, y, x = np.indices((2, 2, 2))
        tic = (1 + 100 * z + 10 * y + x).astype(float)

        raster = _tic_raster(tic, self._positions(self.VOLUME))

        nan = np.nan
        expected = [[1.0, 2.0, nan, 101.0, 102.0], [11.0, 12.0, nan, 111.0, nan]]
        np.testing.assert_array_equal(raster, expected)

    def test_a_volume_of_one_plane_still_marks_what_was_not_measured(self):
        """Its table carries z, all of one value; its TIC is a single plane."""
        points = [(x, y, 0) for x, y in L_SHAPE]
        tic = np.full((8, 13), 7.0)

        raster = _tic_raster(tic, self._positions(points))

        assert raster.shape == (3, 3)
        assert np.isnan(raster[2, 2]) and np.count_nonzero(np.isnan(raster)) == 1

    def test_an_image_the_positions_do_not_index_is_shown_as_stored(self):
        tic = np.full((2, 2), 7.0)

        too_small = _tic_raster(tic, self._positions(L_SHAPE))
        several_slices = _tic_raster(tic, self._positions(self.VOLUME))

        assert too_small is tic and several_slices is tic

    def test_an_image_with_no_pixels_is_left_out(self):
        assert _tic_raster(np.zeros((0, 0)), self._positions(L_SHAPE)) is None

    def test_a_volume_goes_through_the_whole_page(self, tmp_path):
        z, y, x = np.indices((2, 2, 2))
        tic = (1 + 100 * z + 10 * y + x).astype(float)[np.newaxis]
        table = make_table(positions=self.VOLUME)
        path = write_store(tmp_path / "v.zarr", {"ds": table}, images={"ds_tic": [tic]})

        page = build_report(open_store(path))

        assert "<td>2 x 2 pixels, 2 slices</td>" in page
        tic_image = np.asarray(_embedded_images(page)[0])
        assert tic_image.shape[1] / tic_image.shape[0] == 2.5  # 5 wide, 2 high
        assert tuple(tic_image[-1, -1]) == (190, 190, 190)


class TestOpticalPixels:
    def test_sixteen_bit_values_are_brought_to_eight(self):
        pixels = np.array([[[0, 32768, 65535]]] * 3, dtype=np.uint16)

        np.testing.assert_array_equal(_as_photo(pixels)[0, :, 0], [0, 128, 255])

    def test_signed_values_do_not_wrap_round(self):
        pixels = np.array([[[-100, 0, 100]]], dtype=np.int16)

        low, middle, high = (int(value) for value in _as_photo(pixels)[0])
        assert (low, high) == (0, 255) and middle in (127, 128)

    def test_fewer_than_three_channels_is_the_first_in_grey(self):
        pixels = np.stack([np.full((2, 2), 9), np.zeros((2, 2))]).astype(np.uint8)

        photo = _as_photo(pixels)

        assert photo.shape == (2, 2) and int(photo.max()) == 9

    def test_a_fourth_channel_is_dropped(self):
        assert _as_photo(np.zeros((4, 2, 5), dtype=np.uint8)).shape == (2, 5, 3)


class TestTableFacts:
    def _facts(self, tmp_path, **metadata):
        table = open_store(_store(tmp_path, **metadata)).open_table("ds_z0")
        block = _metadata(**metadata)
        return dict(table_facts(table, table.positions(), block))

    def test_the_facts_come_from_the_table_and_its_metadata(self, tmp_path):
        assert self._facts(tmp_path) == {
            "Instrument": "Bruker timsTOF fleX",
            "Pixel size": "25 µm",
            "Raster": "3 x 3 pixels",
            "Measured pixels": str(len(L_SHAPE)),
            "m/z range": "500.0000 to 509.0000",
            "m/z bins": "10",
            "Converted with": "Thyra 4.7.0",
            "Source file": "run 7.imzML",
        }

    def test_a_pixel_that_is_not_square_states_both_sides(self, tmp_path):
        facts = self._facts(tmp_path, pixel_size_um={"x": 20.0, "y": 50.0})

        assert facts["Pixel size"] == "20 x 50 µm"

    def test_a_block_with_nothing_in_it_says_so(self, tmp_path):
        table = open_store(_store(tmp_path)).open_table("ds_z0")

        facts = dict(table_facts(table, table.positions(), {}))

        assert facts["Instrument"] == "not stated in the source"
        assert facts["Pixel size"] == "not stated"
        assert facts["Source file"] == "not recorded"


class TestTheReport:
    @pytest.fixture
    def page(self, tmp_path):
        tic = np.zeros((1, 8, 13))
        tic[0, 5:8, 10:13] = 50.0
        photo = np.zeros((3, 40, 60), dtype=np.uint8)
        attrs = {
            "optical_images": {
                "elements": {"ds_optical_highres": {"source_file": "slide_0000.tif"}}
            }
        }
        path = _store(
            tmp_path,
            images={"ds_z0_tic": [tic], "ds_optical_highres": [photo]},
            attrs=attrs,
        )
        return build_report(open_store(path))

    def test_it_refers_to_nothing_outside_itself(self, page):
        assert "http://" not in page and "https://" not in page
        assert "<script" not in page.lower()
        assert "<link" not in page.lower() and "@import" not in page
        assert "url(" not in page
        sources = re.findall(r'(?:src|href)="([^"]*)"', page)
        assert sources and all(
            source.startswith("data:image/png;base64,") for source in sources
        )

    def test_it_holds_no_path_and_no_name_from_the_source_metadata(
        self, page, tmp_path
    ):
        assert "run 7.imzML" in page
        for private in ("jdoe", "acquisitions", "slide 7", "D:", "SECRET", "Jane Doe"):
            assert private not in page
        assert str(tmp_path) not in page and tmp_path.name not in page

    def test_a_posix_source_path_is_cut_down_too(self, tmp_path):
        page = build_report(open_store(_store(tmp_path, source_path=POSIX_SOURCE)))

        assert "run 7.imzML" in page
        assert "jdoe" not in page and "/home" not in page

    def test_the_sections_come_in_the_stated_order(self, page):
        marks = [
            "<h1>Brain 1</h1>",
            "Measured pixels",
            ">TIC image<",
            ">Mean spectrum<",
            "<svg",
            ">Ion images of the tallest peaks<",
            ">Optical image<",
            ">Tables<",
        ]
        found = [page.index(mark) for mark in marks]

        assert found == sorted(found)
        assert len(re.findall(r"<h[23]>", page)) == 5

    def test_every_image_is_an_embedded_png(self, page):
        images = _embedded_images(page)

        assert {image.format for image in images} == {"PNG"}
        # The TIC, one image for each of the two peaks, the optical image.
        assert len(images) == 4

    def test_the_ion_images_are_of_the_peaks_of_the_mean_spectrum(self, page):
        """The mean spectrum rises to 509, with a gap at 503: two peaks."""
        captions = re.findall(r"<figcaption>m/z ([\d.]+),", page)

        assert captions == ["509.0000", "502.0000"]

    def test_the_command_it_suggests_can_be_pasted(self, page, tmp_path):
        """A store name with a space in it is quoted."""
        assert "<code>thyra ion-image &quot;Brain 1.zarr&quot; 509.0000</code>" in page

        plain = write_store(tmp_path / "plain.zarr", {"t": make_table()})
        assert "<code>thyra ion-image plain.zarr 509.0000</code>" in build_report(
            open_store(plain)
        )

    def test_the_sibling_tables_are_listed_with_their_shapes(self, page):
        rows = re.findall(r"<tr><th>(ds_[^<]*)</th>((?:<td>[^<]*</td>)+)</tr>", page)
        listed = {name: re.findall(r"<td>([^<]*)</td>", cells) for name, cells in rows}

        assert listed == {
            "ds_z0": ["spectra", "8", "10"],
            "ds_z0_mobility": ["ion mobility", "8", "10"],
            "ds_z0_msms": ["MS/MS", "8", "4"],
        }

    def test_the_optical_image_is_named_by_its_file_and_reduced(self, page):
        assert "ds_optical_highres, from slide_0000.tif" in page

    def test_the_tic_image_marks_what_was_not_measured(self, page):
        """The stored TIC is a full grid; the page draws the measured pixels."""
        tic = np.asarray(_embedded_images(page)[0])

        assert tic.shape[0] == tic.shape[1]  # the 3 x 3 box of the L, enlarged
        block = tic.shape[0] // 3
        assert tuple(tic[-1, -1]) == (190, 190, 190)
        assert tuple(tic[0, 0]) != (190, 190, 190)
        assert block * 3 == tic.shape[0]

    def test_it_reads_in_a_light_and_in_a_dark_browser(self, page):
        assert '<meta name="color-scheme" content="light dark">' in page
        assert "color-scheme: light dark" in page
        style = page[page.index("<style>") : page.index("</style>")]
        assert "color:" not in style.replace("color-scheme:", "")

    def test_a_store_without_an_optical_image_has_no_such_section(self, tmp_path):
        page = build_report(open_store(_store(tmp_path)))

        assert "Optical image" not in page
        assert ">TIC image<" not in page  # nor a TIC section without a TIC image
        assert ">Tables<" in page

    def test_names_from_the_store_are_escaped(self, tmp_path):
        path = write_store(tmp_path / "R&D.zarr", {"a&b": make_table()})

        page = build_report(open_store(path))

        assert "<th>a&amp;b</th>" in page and "a&b" not in page
        assert "<h1>R&amp;D</h1>" in page and "R&D" not in page

    def test_a_matrix_stored_by_row_still_gets_the_rest_of_the_page(self, tmp_path):
        """Only the ion images need the matrix; the reason stands in for them."""
        tic = np.ones((1, 8, 13))
        path = write_store(
            tmp_path / "old.zarr", {"t": make_table(csr=True)}, images={"t_tic": [tic]}
        )

        page = build_report(open_store(path))

        assert ">TIC image<" in page and ">Mean spectrum<" in page
        assert "does not store its matrix column by column" in page
        assert len(_embedded_images(page)) == 1  # the TIC image only

    def test_a_damaged_mean_spectrum_does_not_blank_the_drawing(self, tmp_path):
        spectrum = np.arange(1.0, MZ.size + 1)
        spectrum[4] = np.nan
        table = make_table(uns={"average_spectrum": spectrum})
        path = write_store(tmp_path / "s.zarr", {"t": table})

        page = build_report(open_store(path))

        assert "nan" not in page.lower()
        assert ">10<" in page  # the tallest value still labels the axis

    def test_each_slice_table_gets_its_own_sections(self, tmp_path):
        tables = {"ds_z0": make_table(), "ds_z1": make_table()}
        page = build_report(open_store(write_store(tmp_path / "s.zarr", tables)))

        assert "<h2>ds_z0</h2>" in page and "<h2>ds_z1</h2>" in page
        assert page.count("<h3>Mean spectrum</h3>") == 2
