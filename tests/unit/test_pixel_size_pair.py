"""The pixel size override takes an x and y pair (issue #570).

``--pixel-size`` and ``convert_msi(pixel_size_um=...)`` took one number and
applied it to both axes. A raster whose pixels are not square, a line scan
of 40 by 150 um for one, could then be corrected by hand only when the
source already stated both values. The override now takes the pair:
``--pixel-size 40x150`` on the command line, ``(40, 150)`` in the API.

The fixture states 50 by 50 um, so every pair below contradicts the source
and has to win. The assertions are on the stored artefacts, through the API
and through the command line.
"""

from __future__ import annotations

import base64
import io
import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest
import spatialdata
from click.testing import CliRunner
from PIL import Image
from spatialdata.transformations import get_transformation

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.__main__ import PixelSizeParam, _validate_basic_params, main
from thyra.convert import convert_msi, pixel_size_problem
from thyra.errors import ConversionRefused
from thyra.metadata.schema import MSI_METADATA_UNS_KEY

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "fixtures"
#: An imzML of 3 by 2 pixels that states 50 by 50 um.
DECLARED = FIXTURES / "pixel_size_unit_declared.imzML"
#: A timsTOF folder, for a second format.
TIMS = FIXTURES / "synthetic_tims.d"

X_UM, Y_UM = 40.0, 150.0
N_X, N_Y = 3, 2
KEY = "msi_dataset_z0"

#: What a store records about a size given by hand, for one value and a pair.
MANUAL_RECORD = {
    "method": "manual",
    "source_format": "imzml",
    "detection_successful": False,
    "note": "Pixel size manually specified via --pixel-size parameter",
}


def _root_attrs(store: Path) -> Dict[str, Any]:
    return json.loads((store / "zarr.json").read_text())["attributes"]


def _by_api(folder: Path, pixel_size_um: Any, source: Path = DECLARED) -> Path:
    store = folder / "api.zarr"
    assert convert_msi(
        source, store, pixel_size_um=pixel_size_um, include_optical=False
    )
    return store


def _by_cli(folder: Path, token: str) -> Path:
    store = folder / "cli.zarr"
    result = CliRunner().invoke(
        main, [str(DECLARED), str(store), "--pixel-size", token, "--no-optical"]
    )
    assert result.exit_code == 0, result.output
    return store


@pytest.fixture(scope="module")
def api_store(tmp_path_factory) -> Path:
    return _by_api(tmp_path_factory.mktemp("pair-api"), (X_UM, Y_UM))


@pytest.fixture(scope="module")
def cli_store(tmp_path_factory) -> Path:
    return _by_cli(tmp_path_factory.mktemp("pair-cli"), "40x150")


@pytest.fixture(scope="module", params=["api", "cli"])
def store(request, api_store, cli_store) -> Path:
    """The same conversion, by each of its two routes."""
    return api_store if request.param == "api" else cli_store


class TestAStatedPairIsInEveryBlock:
    """40 by 150 um, wherever the store states or uses a pixel size."""

    def test_the_root_attributes(self, store):
        attrs = _root_attrs(store)
        assert attrs["pixel_size_x_um"] == X_UM
        assert attrs["pixel_size_y_um"] == Y_UM
        cs = attrs["coordinate_systems"]["global"]
        assert (cs["pixel_size_um_x"], cs["pixel_size_um_y"]) == (X_UM, Y_UM)

    def test_the_pixel_shapes(self, store):
        shapes = spatialdata.SpatialData.read(str(store)).shapes[f"{KEY}_pixels"]
        assert len(shapes) == N_X * N_Y
        for polygon in shapes.geometry:
            xmin, ymin, xmax, ymax = polygon.bounds
            assert xmax - xmin == pytest.approx(X_UM)
            assert ymax - ymin == pytest.approx(Y_UM)

    def test_the_positions_step_by_each_size(self, store):
        obs = spatialdata.SpatialData.read(str(store)).tables[KEY].obs
        np.testing.assert_array_equal(obs["spatial_x"], obs["x"] * X_UM)
        np.testing.assert_array_equal(obs["spatial_y"], obs["y"] * Y_UM)
        assert set(np.diff(np.unique(obs["spatial_x"]))) == {X_UM}
        assert set(np.diff(np.unique(obs["spatial_y"]))) == {Y_UM}

    def test_the_tic_image_lies_exactly_on_the_shapes(self, store):
        """Each TIC cell covers its pixel's shape, edge for edge."""
        sdata = spatialdata.SpatialData.read(str(store))
        matrix = np.asarray(
            get_transformation(
                sdata.images[f"{KEY}_tic"], to_coordinate_system="global"
            ).to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y"))
        )
        np.testing.assert_array_equal(
            matrix,
            [[X_UM, 0.0, -0.5 * X_UM], [0.0, Y_UM, -0.5 * Y_UM], [0.0, 0.0, 1.0]],
        )
        obs = sdata.tables[KEY].obs
        shapes = sdata.shapes[f"{KEY}_pixels"]
        for label in obs.index:
            x, y = obs.loc[label, "x"], obs.loc[label, "y"]
            low = matrix @ [x, y, 1.0]
            high = matrix @ [x + 1, y + 1, 1.0]
            assert (*low[:2], *high[:2]) == pytest.approx(
                shapes.geometry.loc[label].bounds
            )

    def test_the_image_and_the_shapes_cover_the_same_rectangle(self, store):
        sdata = spatialdata.SpatialData.read(str(store))
        for element in (sdata.images[f"{KEY}_tic"], sdata.shapes[f"{KEY}_pixels"]):
            extent = spatialdata.get_extent(element, coordinate_system="global")
            assert extent["x"] == pytest.approx((-0.5 * X_UM, (N_X - 0.5) * X_UM))
            assert extent["y"] == pytest.approx((-0.5 * Y_UM, (N_Y - 0.5) * Y_UM))

    def test_the_metadata_block(self, store):
        block = (
            spatialdata.SpatialData.read(str(store))
            .tables[KEY]
            .uns[MSI_METADATA_UNS_KEY]
        )
        assert block["ms_analysis"]["pixel_size_um"] == {"x": X_UM, "y": Y_UM}
        assert block["provenance"]["pixel_size_source"] == "manual"

    def test_it_is_recorded_as_given_by_hand_like_one_value(self, store, tmp_path):
        """The same record, key for key, as a single value leaves."""
        pair = _root_attrs(store)
        one = _root_attrs(_by_api(tmp_path, 40))

        assert pair["pixel_size_detection_info"] == MANUAL_RECORD
        assert one["pixel_size_detection_info"] == MANUAL_RECORD
        assert pair["conversion_options"] == {
            **one["conversion_options"],
            "pixel_size_y_um": Y_UM,
        }
        assert one["conversion_options"]["pixel_size_y_um"] == X_UM


class TestTheSpellings:
    def test_a_capital_x_reads_the_same(self, tmp_path):
        capital = _root_attrs(_by_cli(tmp_path, "40X150"))
        assert capital["pixel_size_x_um"] == X_UM
        assert capital["pixel_size_y_um"] == Y_UM

    @pytest.mark.parametrize("token", ["40", "40.0", "4e1"])
    def test_one_number_still_means_both_axes(self, tmp_path, token):
        attrs = _root_attrs(_by_cli(tmp_path, token))
        assert attrs["pixel_size_x_um"] == attrs["pixel_size_y_um"] == X_UM

    @pytest.mark.parametrize(
        "token, read",
        [("40", 40.0), ("40x150", (40.0, 150.0)), ("0.5X2e1", (0.5, 20.0))],
    )
    def test_what_the_option_hands_on(self, token, read):
        value = PixelSizeParam().convert(token, None, None)
        assert value == read
        assert type(value) is type(read)

    @pytest.mark.parametrize(
        "pair",
        [
            (40, 150),
            [40, 150],
            (np.float32(40), np.int64(150)),
            [np.float64(40.0), 150],
        ],
        ids=["tuple", "list", "numpy", "mixed"],
    )
    def test_the_api_takes_a_tuple_or_a_list_of_any_real_numbers(self, tmp_path, pair):
        store = _by_api(tmp_path, pair)
        attrs = _root_attrs(store)
        assert attrs["pixel_size_x_um"] == X_UM
        assert attrs["pixel_size_y_um"] == Y_UM
        obs = spatialdata.SpatialData.read(str(store)).tables[KEY].obs
        assert obs["spatial_x"].dtype == np.float64
        assert obs["spatial_y"].dtype == np.float64

    def test_a_second_format_takes_the_pair(self, tmp_path):
        from thyra.readers.bruker.timstof.timstof_reader import BrukerReader
        from thyra.utils.bruker_exceptions import SDKError

        try:
            with BrukerReader(TIMS):
                pass
        except (SDKError, OSError) as exc:  # pragma: no cover - platform
            pytest.skip(f"Bruker library not loadable on this platform: {exc}")

        store = _by_api(tmp_path, (X_UM, Y_UM), source=TIMS)
        attrs = _root_attrs(store)
        assert attrs["pixel_size_x_um"] == X_UM
        assert attrs["pixel_size_y_um"] == Y_UM
        assert attrs["pixel_size_detection_info"]["method"] == "manual"
        sdata = spatialdata.SpatialData.read(str(store))
        key = next(k for k in sdata.tables if not k.endswith(("_mobility", "_msms")))
        obs = sdata.tables[key].obs
        np.testing.assert_array_equal(obs["spatial_x"], obs["x"] * X_UM)
        np.testing.assert_array_equal(obs["spatial_y"], obs["y"] * Y_UM)


#: Each is refused by the command line before the source is opened.
REFUSED_TOKENS = [
    "40x",
    "x150",
    "40x150x10",
    "40,150",
    "40 150",
    "40 x 150",
    "forty",
    "0",
    "-5",
    "nan",
    "inf",
    "0x150",
    "40x0",
    "40x-150",
    "40xnan",
    "infx150",
]


class TestTheCommandLineRefusals:
    @pytest.mark.parametrize("token", REFUSED_TOKENS)
    def test_it_is_refused_in_one_line_before_the_source_is_opened(
        self, tmp_path, monkeypatch, token
    ):
        def never(*_args, **_kwargs):  # pragma: no cover - must not run
            raise AssertionError("the conversion started")

        monkeypatch.setattr("thyra.convert.convert_msi", never)
        monkeypatch.setattr("thyra.convert._create_reader", never)
        store = tmp_path / "out.zarr"
        result = CliRunner().invoke(
            main, [str(DECLARED), str(store), "--pixel-size", token]
        )

        assert result.exit_code == 2, result.output
        errors = [
            line for line in result.output.splitlines() if line.startswith("Error")
        ]
        assert len(errors) == 1, result.output
        assert "'--pixel-size'" in errors[0]
        assert "--pixel-size 40)" in errors[0]
        assert "--pixel-size 40x150)" in errors[0]
        assert "Traceback" not in result.output
        assert list(tmp_path.iterdir()) == []

    @pytest.mark.parametrize(
        "pair", [(40.0, 0.0), (-1.0, 150.0), (40.0, float("nan")), (float("inf"), 1.0)]
    )
    def test_the_value_check_reads_both_halves(self, pair):
        import click

        with pytest.raises(click.BadParameter, match="finite positive"):
            _validate_basic_params(pair, "dataset")

    def test_a_usable_pair_passes_the_value_check(self):
        _validate_basic_params((40.0, 150.0), "dataset")


#: Each is refused by ``convert_msi``, with the part of the reason to find.
REFUSED_VALUES = [
    (True, "one number or an (x, y) pair of two numbers, not bool"),
    ((40, 150, 10), "got 3 values"),
    ([40], "got 1 values"),
    ((), "got 0 values"),
    ((40, 0), "Pixel size y must be a finite positive number"),
    ((-40, 150), "Pixel size x must be a finite positive number"),
    ((40, float("nan")), "Pixel size y must be a finite positive number"),
    ((float("inf"), 150), "Pixel size x must be a finite positive number"),
    ((True, 150), "Pixel size x must be a number, not bool"),
    ((40, "150"), "Pixel size y must be a number, not str"),
    ((40, None), "Pixel size y must be a number, not NoneType"),
    ("40x150", "one number or an (x, y) pair of two numbers, not str"),
    ({"x": 40, "y": 150}, "an (x, y) pair of two numbers, not dict"),
    (np.array([40.0, 150.0]), "an (x, y) pair of two numbers, not ndarray"),
    (np.bool_(True), "an (x, y) pair of two numbers, not bool"),
]


class TestTheApiRefusals:
    @pytest.fixture
    def no_reader(self, monkeypatch):
        def never(*_args, **_kwargs):  # pragma: no cover - must not run
            raise AssertionError("the source was opened")

        monkeypatch.setattr("thyra.convert._create_reader", never)

    @pytest.mark.parametrize("value, said", REFUSED_VALUES)
    def test_it_is_refused_in_one_line_before_the_source_is_opened(
        self, tmp_path, thyra_logs, no_reader, value, said
    ):
        store = tmp_path / "out.zarr"
        with thyra_logs("thyra.convert", logging.ERROR) as records:
            assert convert_msi(DECLARED, store, pixel_size_um=value) is False

        assert len(records) == 1, [r.getMessage() for r in records]
        assert said in records[0].getMessage()
        assert "\n" not in records[0].getMessage()
        assert list(tmp_path.iterdir()) == []

    @pytest.mark.parametrize("value, said", REFUSED_VALUES)
    def test_it_is_raised_as_a_refusal_when_asked(
        self, tmp_path, no_reader, value, said
    ):
        with pytest.raises(ConversionRefused, match=re.escape(said)):
            convert_msi(
                DECLARED,
                tmp_path / "out.zarr",
                pixel_size_um=value,
                raise_on_failure=True,
            )
        assert list(tmp_path.iterdir()) == []

    @pytest.mark.parametrize(
        "value", [40, 40.0, np.float32(40), (40, 150), [40.0, 150.0]]
    )
    def test_a_usable_size_has_no_problem(self, value):
        assert pixel_size_problem(value) is None


class TestAVolume:
    """``handle_3d`` takes the pair; the z fallback keeps the x pitch."""

    @pytest.fixture
    def volume(self, tmp_path, monkeypatch, thyra_logs):
        reader = MockMSIReader(
            MockMSIConfig(
                n_x=N_X,
                n_y=N_Y,
                n_z=2,
                n_mz_bins=64,
                peaks_per_spectrum=(3, 5),
                seed=5,
            )
        )
        monkeypatch.setattr(
            "thyra.convert._create_reader", lambda *_a, **_k: (reader, "mock")
        )
        source = tmp_path / "volume.imzML"
        source.write_bytes(b"")
        store = tmp_path / "volume.zarr"
        with thyra_logs("thyra.core.base_converter", logging.WARNING) as records:
            assert convert_msi(
                source,
                store,
                dataset_id="v",
                pixel_size_um=(X_UM, Y_UM),
                handle_3d=True,
                include_optical=False,
            )
        return store, [r.getMessage() for r in records]

    def test_both_sizes_reach_the_volume(self, volume):
        store, _ = volume
        attrs = _root_attrs(store)
        assert attrs["pixel_size_x_um"] == X_UM
        assert attrs["pixel_size_y_um"] == Y_UM

        sdata = spatialdata.SpatialData.read(str(store))
        obs = sdata.tables["v"].obs
        np.testing.assert_array_equal(obs["spatial_x"], obs["x"] * X_UM)
        np.testing.assert_array_equal(obs["spatial_y"], obs["y"] * Y_UM)
        matrix = np.asarray(
            get_transformation(
                sdata.images["v_tic"], to_coordinate_system="global"
            ).to_affine_matrix(input_axes=("x", "y", "z"), output_axes=("x", "y", "z"))
        )
        assert [matrix[0, 0], matrix[1, 1]] == [X_UM, Y_UM]
        for polygon in sdata.shapes["v_pixels"].geometry:
            xmin, ymin, xmax, ymax = polygon.bounds
            assert (xmax - xmin, ymax - ymin) == pytest.approx((X_UM, Y_UM))

    def test_the_z_spacing_falls_back_to_the_x_size_and_says_so(self, volume):
        store, warnings = volume
        attrs = _root_attrs(store)
        cs = attrs["coordinate_systems"]["global"]
        assert cs["z_spacing_um"] == X_UM
        assert cs["z_spacing_source"] == "assumed_isotropic"
        assert attrs["conversion_options"]["z_spacing_um"] == X_UM

        sdata = spatialdata.SpatialData.read(str(store))
        np.testing.assert_array_equal(
            sdata.tables["v"].obs["spatial_z"], sdata.tables["v"].obs["z"] * X_UM
        )
        assumed = [m for m in warnings if "No z spacing was supplied" in m]
        assert len(assumed) == 1, warnings
        assert "(40 um)" in assumed[0]
        assert "--z-spacing" in assumed[0]


def _embedded_sizes(page: str) -> List[Tuple[int, int]]:
    """The size of each image embedded in a report, in page order."""
    found = re.findall(r'src="data:image/png;base64,([^"]+)"', page)
    return [Image.open(io.BytesIO(base64.b64decode(data))).size for data in found]


class TestTheQuickLooksDrawItToScale:
    """3 by 2 pixels of 40 by 150 um cover 120 by 300 um: taller than wide."""

    def test_the_report(self, api_store, tmp_path):
        from thyra.quicklook.cli import report_command

        page = tmp_path / "report.html"
        result = CliRunner().invoke(report_command, [str(api_store), "-o", str(page)])
        assert result.exit_code == 0, result.output

        sizes = _embedded_sizes(page.read_text(encoding="utf-8"))
        assert sizes, "the report holds no image"
        for width, height in sizes:
            assert abs(height * (N_X * X_UM) / (N_Y * Y_UM) - width) <= 1, sizes

    def test_the_ion_image(self, api_store, tmp_path):
        from thyra.quicklook.cli import ion_image_command

        sizes = {}
        for name, store in (("pair", api_store), ("square", _by_api(tmp_path, 40))):
            target = tmp_path / f"{name}.png"
            result = CliRunner().invoke(
                ion_image_command, [str(store), "184", "-o", str(target)]
            )
            assert result.exit_code == 0, result.output
            sizes[name] = Image.open(target).size

        assert sizes["pair"][1] > sizes["pair"][0], sizes
        assert sizes["square"][0] > sizes["square"][1], sizes
