"""What the four stages of a conversion hand each other.

``BaseMSIConverter.convert()`` is a template method whose stages
communicate through one object:

    state = self._create_data_structures()
    self._process_spectra(state)
    self._finalize_data(state)
    success = self._save_output(state)

That object used to be a ``Dict[str, Any]``. Nothing declared what a stage
required or produced, the twelve keys in use could only be learned by
grepping for bracket literals, and no checker could verify any of it in a
package that ships the ``Typing :: Typed`` classifier (issue #273).

The fields below are that dictionary's keys, typed. Two things the dict
could not express and this does:

* **Which fields are optional.** ``region_*`` exist only for a dataset
  with a region map, and the old dict said so by simply not having the
  key -- so some call sites used ``state["region_row_count"]`` and others
  ``state.get("region_total_intensity")``, and only the author knew which
  was safe where. They are ``None`` here, once, for everyone.
* **Which are filled in later.** ``tables``, ``shapes`` and ``images``
  start empty and are filled by the finalize stage; the counters start at
  zero and are set by the scatter pass. Their defaults say so.

``units`` and ``passes`` are deliberately loose. A table unit is private
to the streaming converter and a fused-pass record is built per reader, so
typing them here would make :mod:`thyra.core` depend on a converter's
internals to describe a container it only passes along.
"""

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, List, Optional

import numpy as np
import pandas as pd
from numpy.typing import NDArray

if TYPE_CHECKING:  # pragma: no cover - typing only
    from anndata import AnnData


@dataclass
class ConversionState:
    """The working state of one conversion, threaded through its stages.

    Attributes:
        units: The tables being built -- one per z plane, or one for the
            whole volume. Opaque here; see the streaming converter.
        var_df: The feature frame every table copies for its ``var``.
        passes: The fused sibling passes when the reader hands its frames
            over as records (design decision D5), else ``None``.
        tables: Element key -> table, filled by the finalize stage.
        shapes: Element key -> pixel shapes, filled alongside the tables.
        images: Element key -> image, the TIC rasters and any optical
            images.
        pixel_count: Grid positions that carried a spectrum, counted by
            the scatter pass.
        input_peaks: Peaks read from the source, before binning. Compared
            against what reached the store to report what was dropped.
        out_of_grid_spectra: Spectra whose coordinates fell outside the
            declared grid, counted rather than silently discarded.
        region_total_intensity: Region -> summed intensity per m/z bin, or
            ``None`` when the dataset declares no regions.
        region_row_count: Region -> rows contributing to that sum, or
            ``None`` for the same reason.
        avg_spectrum_per_region: Region -> mean spectrum, derived from the
            two fields above once the counting pass is done. Keyed by
            ``str``, unlike the two ``int``-keyed fields it is computed
            from: it goes straight into ``uns`` and a Zarr group's keys
            are strings. Typing the dict is what surfaced the difference
            -- under ``Dict[str, Any]`` the three region fields looked
            alike and mypy had nothing to check.
    """

    units: List[Any]
    var_df: pd.DataFrame
    passes: Optional[Any] = None
    tables: Dict[str, "AnnData"] = field(default_factory=dict)
    shapes: Dict[str, Any] = field(default_factory=dict)
    images: Dict[str, Any] = field(default_factory=dict)
    pixel_count: int = 0
    input_peaks: int = 0
    out_of_grid_spectra: int = 0
    region_total_intensity: Optional[Dict[int, NDArray[np.float64]]] = None
    region_row_count: Optional[Dict[int, int]] = None
    avg_spectrum_per_region: Optional[Dict[str, NDArray[np.float64]]] = None

    @property
    def has_regions(self) -> bool:
        """Whether this dataset declared a region map."""
        return self.region_row_count is not None
