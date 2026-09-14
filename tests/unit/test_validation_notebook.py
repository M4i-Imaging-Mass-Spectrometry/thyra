"""The Colab validation notebook, checked as source (issues #296, #316).

`README.md` and `docs/index.md` both badge this notebook as the way to try
Thyra, and nothing tested it -- `grep -rn ipynb tests/` was empty. Three
things were wrong with it and none of them would have been caught:

- **Cell 14 picked an image by position.** `sdata.images` is sorted and an
  optical image sorts before the TIC (`..._optical_highres` <
  `..._z0_tic`), so on any store converted with `--include-optical` -- the
  default -- it showed the optical photo captioned as an ion image, or
  crashed with `AttributeError: 'DataTree' object has no attribute 'data'`
  when that photo was large enough to be written multiscale.
- **Cell 4 reached for a POSIX shell.** `!mkdir -p data` under cmd.exe
  makes two directories, one called `-p`, and reports success.
- **Cell 9 printed a boolean.** `convert_msi` reports every failure by
  returning False, so the notebook announced "Conversion succeeded: False"
  and carried on into a `FileNotFoundError` from inside zarr.

These are source-shape assertions rather than an executed notebook: running
it needs Colab, a 30 MB generated dataset and a full conversion. They are
cheap and they each fail on the shape that was there.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

_NOTEBOOK = (
    Path(__file__).resolve().parents[2]
    / "notebooks"
    / "Thyra_Validation_Workflow.ipynb"
)


@pytest.fixture(scope="module")
def cells():
    nb = json.loads(_NOTEBOOK.read_text(encoding="utf-8"))
    return [
        (i, "".join(c["source"]))
        for i, c in enumerate(nb["cells"])
        if c["cell_type"] == "code"
    ]


def test_the_notebook_is_where_the_badges_point():
    assert _NOTEBOOK.is_file()


class TestItRunsOffColab:
    #: Cell 5 is the Zenodo download, behind ``DOWNLOAD_FULL_DATASET = False``
    #: and explicitly a Colab-only path. Everything else must not need a shell.
    COLAB_ONLY = {5}

    def test_no_cell_shells_out(self, cells):
        offenders = [
            index
            for index, source in cells
            if index not in self.COLAB_ONLY
            and any(line.lstrip().startswith("!") for line in source.split("\n"))
        ]
        assert not offenders, f"cells {offenders} use the IPython shell escape"

    def test_the_dataset_cell_uses_the_running_interpreter(self, cells):
        source = dict(cells)[4]
        assert "sys.executable" in source
        assert "thyra.tools.make_example_data" in source

    def test_the_dataset_cell_fails_loudly(self, cells):
        """``check=True``: a generator failure must not be printed and passed."""
        assert "check=True" in dict(cells)[4]


class TestAFailedConversionStopsTheNotebook:
    def test_the_conversion_cell_raises_on_false(self, cells):
        source = dict(cells)[9]
        assert "ok = convert_msi(" in source
        assert "raise SystemExit(" in source

    def test_it_raises_after_reporting(self, cells):
        """The message is still printed; the raise follows it."""
        source = dict(cells)[9]
        assert source.index("Conversion succeeded") < source.index("raise SystemExit")

    def test_the_message_names_the_usual_causes(self, cells):
        source = dict(cells)[9]
        assert ".ibd" in source and "pixel_size_um" in source


class TestTheImageIsChosenByName:
    def test_no_cell_picks_an_image_positionally(self, cells):
        offenders = [
            index for index, source in cells if "list(sdata.images.keys())[0]" in source
        ]
        assert not offenders, f"cells {offenders} pick an image by position"

    def test_the_plot_cell_selects_the_tic(self, cells):
        source = dict(cells)[14]
        assert 'endswith("_tic")' in source

    def test_it_handles_a_multiscale_image(self, cells):
        """An optical image with a short side of 2000 px or more is written
        multiscale, and a DataTree has no ``.data``."""
        source = dict(cells)[14]
        assert "DataTree" in source

    def test_it_says_so_when_there_is_no_tic(self, cells):
        source = dict(cells)[14]
        assert "raise RuntimeError(" in source


class TestTheTableCellsAreLeftAlone:
    """Cells 13 and 16 pick a *table* positionally, which is correct.

    The sibling tables are ``{key}_mobility`` and ``{key}_msms``, and
    ``{key}`` is a strict prefix of both, so sorted order always puts the
    spectral table first. Not the same defect -- this pins that nobody
    "fixes" it by analogy.
    """

    def test_they_still_index_position_zero(self, cells):
        source = dict(cells)[13]
        assert "list(sdata.tables.keys())[0]" in source
