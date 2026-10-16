"""The region option and the region listing (issue #448).

Only the timsTOF reader can select a region. Every other reader refuses
``region`` now, through the API as well as the CLI: solariX, rapifleX,
imzML, PHI and Waters dropped it in silence, and on mzPeak the CLI logged
it as ignored just before the reader refused it. :func:`thyra.list_regions`
lists the regions of any input without converting it.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

import thyra
from tests.unit.readers.test_mzpeak_archive_metadata import _region
from tests.unit.readers.test_mzpeak_archive_metadata import (
    _two_regions as _mzpeak_two_regions,
)
from tests.unit.readers.test_rapiflex_regions import _two_regions as _rf_cells
from tests.unit.readers.test_rapiflex_regions import _write as _write_rapiflex
from tests.unit.readers.test_solarix_reader import make_solarix_d
from thyra.__main__ import main
from thyra.convert import convert_msi
from thyra.core.registry import _READER_MODULES, get_reader_class
from thyra.errors import ConversionRefused

_FIXTURES = Path(__file__).resolve().parents[1] / "data" / "fixtures"
_IMZML = _FIXTURES / "timstof_flex_export.imzML"
_TDF = _FIXTURES / "synthetic_tims.d"


class TestTheRegionOptionIsRefused:
    @pytest.mark.parametrize("format_name", sorted(set(_READER_MODULES) - {"bruker"}))
    def test_every_reader_but_timstof_refuses_it(self, format_name, tmp_path):
        """The refusal comes before any file is opened."""
        reader_class = get_reader_class(format_name)
        with pytest.raises(ConversionRefused, match="only supported for Bruker"):
            reader_class(tmp_path / "missing", region=1)

    def test_the_api_returns_false_with_one_error(self, tmp_path, thyra_logs):
        with thyra_logs("thyra", logging.WARNING) as records:
            result = convert_msi(_IMZML, tmp_path / "out.zarr", region="01")

        assert result is False
        assert not (tmp_path / "out.zarr").exists()
        errors = [r.getMessage() for r in records if r.levelno >= logging.ERROR]
        assert len(errors) == 1
        assert "only supported for Bruker timsTOF" in errors[0]

    @pytest.mark.parametrize("build", ["imzml", "mzpeak"])
    def test_the_cli_says_one_thing(self, build, tmp_path):
        """No "ignored" line before the refusal."""
        if build == "imzml":
            source = _IMZML
        else:
            source = _mzpeak_two_regions(
                tmp_path,
                [
                    _region(0, 4, (101, 102), (11, 12)),
                    _region(1, 4, (104, 105), (11, 12)),
                ],
            )
        result = CliRunner().invoke(
            main, [str(source), str(tmp_path / "out.zarr"), "--region", "1"]
        )

        assert result.exit_code != 0
        assert "only supported for Bruker timsTOF" in result.output
        assert "ignored" not in result.output.lower()


class TestListRegions:
    def test_it_is_public(self):
        assert "list_regions" in thyra.__all__
        assert thyra.list_regions is not None

    def test_an_input_without_regions_is_one_whole_run_entry(self):
        assert thyra.list_regions(_IMZML) == [
            {
                "region_number": 1,
                "name": None,
                "n_spectra": thyra.preview_msi(_IMZML).n_pixels,
                "bounds": (
                    0,
                    0,
                    thyra.preview_msi(_IMZML).grid_dims[0] - 1,
                    thyra.preview_msi(_IMZML).grid_dims[1] - 1,
                ),
            }
        ]

    def test_rapiflex_regions_carry_names_and_boxes(self, tmp_path):
        folder = _write_rapiflex(tmp_path / "rf", _rf_cells())
        assert thyra.list_regions(folder) == [
            {
                "region_number": 0,
                "name": "left",
                "n_spectra": 4,
                "bounds": (1, 0, 2, 1),
            },
            {
                "region_number": 1,
                "name": "right",
                "n_spectra": 2,
                "bounds": (4, 0, 4, 1),
            },
        ]

    def test_solarix_regions_carry_boxes(self, tmp_path):
        spectra = [
            {
                "id": i + 1,
                "x": 10 + x,
                "y": 20,
                "mzs": [100.0],
                "intensities": [1.0],
                "region": region,
            }
            for i, (x, region) in enumerate([(0, 0), (1, 0), (5, 1)])
        ]
        d_dir = make_solarix_d(tmp_path, spectra=spectra)
        assert thyra.list_regions(d_dir) == [
            {"region_number": 0, "name": None, "n_spectra": 2, "bounds": (0, 0, 1, 0)},
            {"region_number": 1, "name": None, "n_spectra": 1, "bounds": (5, 0, 5, 0)},
        ]

    def test_mzpeak_regions_sorted_by_number(self, tmp_path):
        archive = _mzpeak_two_regions(
            tmp_path,
            [
                _region(1, 4, (104, 105), (11, 12)),
                _region(0, 4, (101, 102), (11, 12), name="left"),
            ],
        )
        assert thyra.list_regions(archive) == [
            {
                "region_number": 0,
                "name": "left",
                "n_spectra": 4,
                "bounds": (0, 0, 1, 1),
            },
            {"region_number": 1, "name": None, "n_spectra": 4, "bounds": (3, 0, 4, 1)},
        ]

    def test_timstof_one_region(self, tmp_path):
        """The DB RegionNumber, 0 here, not the placeholder 1 (D57)."""
        d_dir = tmp_path / _TDF.name
        shutil.copytree(_TDF, d_dir)
        preview = thyra.preview_msi(d_dir)
        assert preview.readable, preview.error

        regions = thyra.list_regions(d_dir)

        assert [r["region_number"] for r in regions] == [0]
        assert regions[0]["bounds"] == (
            0,
            0,
            preview.grid_dims[0] - 1,
            preview.grid_dims[1] - 1,
        )
        assert preview.regions is None

    def test_the_documented_loop_works_on_a_single_region_run(self, tmp_path):
        """docs/api.md: one store per region, from list_regions (#496 G3)."""
        import spatialdata

        d_dir = tmp_path / _TDF.name
        shutil.copytree(_TDF, d_dir)

        stores = []
        for region in thyra.list_regions(d_dir):
            number = region["region_number"]
            out = tmp_path / f"region_{number}.zarr"
            assert convert_msi(
                str(d_dir), str(out), region=number, include_optical=False
            )
            stores.append((number, out))

        assert [number for number, _ in stores] == [0]
        table = spatialdata.read_zarr(stores[0][1]).tables["msi_dataset_z0"]
        assert set(table.obs["region_number"]) == {0}
        assert "average_spectrum_per_region" not in table.uns

    def test_an_unreadable_input_is_refused(self, tmp_path):
        with pytest.raises(ConversionRefused, match="could not be read"):
            thyra.list_regions(tmp_path / "nothing.imzML")
