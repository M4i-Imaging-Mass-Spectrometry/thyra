"""The streaming converter names a skipped setup step instead of crashing on None.

Each guard covers state that only exists after an earlier step ran: the
grid dimensions (initialization), a table's row layout (pass 1's
``finish_counting``) and the in-memory ``obs`` frame. No conversion
route reaches them out of order; before the guards, calling them early
failed with a bare ``TypeError`` on ``None`` that named no step.
"""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
    _obs_frame,
    _TableUnit,
)


def _unit() -> _TableUnit:
    return _TableUnit("t", "t_pixels", 0, n_grid=4, n_cols=3, tic_shape=(2, 2))


def test_locate_before_dimensions_names_the_step() -> None:
    converter = StreamingSpatialDataConverter.__new__(StreamingSpatialDataConverter)
    converter._dimensions = None
    with pytest.raises(ValueError, match="Dimensions are not initialized"):
        converter._locate([_unit()], 0, 0, 0)


def test_row_of_before_finish_counting_names_the_step() -> None:
    with pytest.raises(RuntimeError, match="row layout is not decided yet"):
        _unit().row_of(0)


def test_row_of_after_finish_counting_reads_the_layout() -> None:
    unit = _unit()
    unit.count(2, np.array([0, 1]), np.array([1.0, 2.0]))
    unit.finish_counting()
    assert [unit.row_of(g) for g in range(4)] == [-1, -1, 0, -1]


def test_obs_frame_returns_the_dataframe() -> None:
    adata = AnnData(obs=pd.DataFrame({"x": [0, 1]}, index=["0", "1"]))
    assert _obs_frame(adata) is adata.obs


def test_obs_frame_refuses_anything_else() -> None:
    with pytest.raises(TypeError, match="Expected obs as a DataFrame"):
        _obs_frame(SimpleNamespace(obs=object()))
