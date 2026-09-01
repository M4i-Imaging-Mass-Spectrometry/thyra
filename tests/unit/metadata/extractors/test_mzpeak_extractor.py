"""Metadata extraction from mzPeak archives."""

from __future__ import annotations

import pytest

from tests.fixtures.mzpeak_builder import Spectrum, build_mzpeak, grid_spectra
from thyra.readers.mzpeak import MzPeakReader


def _essential(path):
    """Read essential metadata and close the archive."""
    with MzPeakReader(path) as reader:
        return reader.get_essential_metadata()


class TestPixelSize:
    """Pixel size is matched on accession and normalised to micrometres."""

    def test_micrometre_pixel_size_is_read(self, tmp_path):
        """The common case: IMS:1000046/47 in UO:0000017."""
        archive = build_mzpeak(
            tmp_path / "um.mzpeak", grid_spectra(2, 2), pixel_size=(30.0, 40.0)
        )
        assert _essential(archive).pixel_size == (30.0, 40.0)

    @pytest.mark.parametrize(
        ("unit", "value", "expected"),
        [
            ("UO:0000017", 25.0, 25.0),  # micrometre
            ("UO:0000018", 4406.25, 4.40625),  # nanometre
            ("UO:0000016", 0.05, 50.0),  # millimetre
        ],
    )
    def test_units_are_converted(self, tmp_path, unit, value, expected):
        """Everything downstream is micrometres, so units are folded here.

        The nanometre case is the one that bit the imzML reader: a 4.40625 um
        pixel read as 4406.25 um if the unit is dropped.
        """
        archive = build_mzpeak(
            tmp_path / f"unit_{unit.replace(':', '_')}.mzpeak",
            grid_spectra(2, 2),
            pixel_size=(value, value),
            pixel_size_unit=unit,
        )
        pixel_size = _essential(archive).pixel_size
        assert pixel_size is not None
        assert pixel_size[0] == pytest.approx(expected)
        assert pixel_size[1] == pytest.approx(expected)

    def test_missing_pixel_size_returns_none(self, tmp_path):
        """Real Bruker exports carry no scan-settings block at all.

        The converter then falls through to ``--pixel-size``, exactly as it
        does for an imzML that omits the terms.
        """
        archive = build_mzpeak(
            tmp_path / "nopixel.mzpeak", grid_spectra(2, 2), pixel_size=None
        )
        essential = _essential(archive)
        assert essential.pixel_size is None
        assert essential.has_pixel_size is False

    def test_unsupported_unit_is_refused_not_assumed(self, tmp_path):
        """An unknown unit yields no pixel size rather than a silent scale.

        Passing the number through unconverted would be a scale error that
        nothing downstream could detect.
        """
        archive = build_mzpeak(
            tmp_path / "badunit.mzpeak",
            grid_spectra(2, 2),
            pixel_size=None,
            footer_metadata=False,
            index_metadata={
                "scan_settings_list": [
                    {
                        "id": "scansettings1",
                        "parameters": [
                            {
                                "name": "pixel size (x)",
                                "accession": "IMS:1000046",
                                "value": 25.0,
                                "unit": "UO:9999999",
                            }
                        ],
                    }
                ]
            },
        )
        assert _essential(archive).pixel_size is None

    def test_single_declared_axis_is_treated_as_square(self, tmp_path):
        """A file declaring only IMS:1000046 describes a square pixel."""
        archive = build_mzpeak(
            tmp_path / "oneaxis.mzpeak",
            grid_spectra(2, 2),
            pixel_size=None,
            footer_metadata=False,
            index_metadata={
                "scan_settings_list": [
                    {
                        "id": "scansettings1",
                        "parameters": [
                            {
                                "name": "pixel size (x)",
                                "accession": "IMS:1000046",
                                "value": 12.5,
                                "unit": "UO:0000017",
                            }
                        ],
                    }
                ]
            },
        )
        assert _essential(archive).pixel_size == (12.5, 12.5)


class TestMetadataPlacement:
    """File-level metadata lives in two places and either may be empty."""

    def test_read_from_the_parquet_footer(self, tmp_path):
        """The reference imaging archive puts everything here."""
        archive = build_mzpeak(
            tmp_path / "footer.mzpeak",
            grid_spectra(2, 2),
            pixel_size=(15.0, 15.0),
            index_metadata=None,
            footer_metadata=True,
        )
        assert _essential(archive).pixel_size == (15.0, 15.0)

    def test_read_from_the_index_metadata_object(self, tmp_path):
        """Other reference archives populate the index instead.

        A reader that consults only one of the two locations finds nothing on
        half the corpus, which is why both are merged.
        """
        archive = build_mzpeak(
            tmp_path / "indexmeta.mzpeak",
            grid_spectra(2, 2),
            pixel_size=None,
            footer_metadata=False,
            index_metadata={
                "version": "0.9.0",
                "scan_settings_list": [
                    {
                        "id": "scansettings1",
                        "parameters": [
                            {
                                "name": "pixel size (x)",
                                "accession": "IMS:1000046",
                                "value": 18.0,
                                "unit": "UO:0000017",
                            },
                            {
                                "name": "pixel size y",
                                "accession": "IMS:1000047",
                                "value": 18.0,
                                "unit": "UO:0000017",
                            },
                        ],
                    }
                ],
            },
        )
        assert _essential(archive).pixel_size == (18.0, 18.0)

    def test_container_version_is_reported(self, tmp_path):
        """The draft version appears only in the index metadata object."""
        archive = build_mzpeak(
            tmp_path / "version.mzpeak",
            grid_spectra(2, 2),
            index_metadata={"version": "0.9.0"},
        )
        with MzPeakReader(archive) as reader:
            comprehensive = reader.get_comprehensive_metadata()
        assert comprehensive.format_specific["container_version"] == "0.9.0"


class TestEssentialShape:
    """Dimensions, bounds, counts and representation."""

    def test_dimensions_and_bounds(self, tmp_path):
        """Grid comes from the positions actually present."""
        archive = build_mzpeak(tmp_path / "shape.mzpeak", grid_spectra(4, 3))
        essential = _essential(archive)

        assert essential.dimensions == (4, 3, 1)
        assert essential.coordinate_bounds == (1.0, 4.0, 1.0, 3.0)
        assert essential.n_spectra == 12

    def test_declared_extent_does_not_override_observed_positions(self, tmp_path):
        """A declared grid larger than the data does not pad the output.

        IMS:1000042/43 are kept as provenance only: trusting them would size
        the grid from an assertion rather than from the pixels present.
        """
        archive = build_mzpeak(
            tmp_path / "extent.mzpeak", grid_spectra(2, 2), grid=(64, 64)
        )
        with MzPeakReader(archive) as reader:
            essential = reader.get_essential_metadata()
            comprehensive = reader.get_comprehensive_metadata()

        assert essential.dimensions == (2, 2, 1)
        # Plain names, not accessions: these become uns keys and zarr writes
        # a dict key as a directory name, where a colon is illegal on Windows.
        assert comprehensive.acquisition_params["declared_grid_extent"] == {
            "max_count_of_pixels_x": 64,
            "max_count_of_pixels_y": 64,
        }

    def test_mass_range_spans_every_spectrum(self, tmp_path):
        """Range is the min of the lows and the max of the highs."""
        spectra = [
            Spectrum(1, 1, [100.0, 150.0], [1.0, 2.0]),
            Spectrum(2, 1, [120.0, 900.0], [3.0, 4.0]),
        ]
        archive = build_mzpeak(tmp_path / "range.mzpeak", spectra)
        assert _essential(archive).mass_range == (100.0, 900.0)

    @pytest.mark.parametrize(
        "representation", ["profile spectrum", "centroid spectrum"]
    )
    def test_spectrum_type_is_reported(self, tmp_path, representation):
        """Representation drives the resampling decision tree."""
        archive = build_mzpeak(
            tmp_path / f"rep_{representation.split()[0]}.mzpeak",
            grid_spectra(2, 2),
            spectrum_representation=representation,
        )
        assert _essential(archive).spectrum_type == representation

    def test_source_path_is_the_archive(self, tmp_path):
        """Provenance points at the file that was read."""
        archive = build_mzpeak(tmp_path / "source.mzpeak", grid_spectra(2, 2))
        assert _essential(archive).source_path == str(archive)
