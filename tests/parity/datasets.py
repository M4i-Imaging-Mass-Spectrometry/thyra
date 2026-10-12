"""Where the parity reference data lives, and what each dataset folder holds.

Reference data never enters the repository. One environment variable,
:data:`PARITY_ENV`, names a folder; every subfolder holding a
``parity.json`` is a dataset. A dataset folder looks like::

    <name>/
        parity.json      what to convert and how SCiLS Lab imported it
        source/          the source data, copied unchanged
        scils/           the .slx, its import log or .sip, and the
                         reference.npz / reference.json that
                         tests/parity/scils_extract.py writes
        thyra/           Thyra stores written by the harness (regenerable)

``parity.json`` keys:

``source``
    Path of the input Thyra converts, relative to the dataset folder.
``spectrum``
    ``"profile"`` or ``"centroid"``: what SCiLS Lab imported. It decides
    the resampling method the harness asks Thyra for (see
    :func:`tests.parity.thyra_runner.matched_resampling_config`).
``reader_options``
    Optional; passed to ``convert_msi`` unchanged.
``pixel_size_um``
    Optional; passed to ``convert_msi`` for a source that does not state
    one Thyra will read, set to what was entered into SCiLS Lab.
``method``
    Optional; overrides the resampling method the harness would choose.
``axis_alignment``
    Optional; ``"scils_values"`` widens the requested range by half a bin
    so Thyra's bin midpoints land on SCiLS Lab's axis values (see
    :func:`tests.parity.thyra_runner.matched_resampling_config`).
``resampling_overrides``
    Optional dict merged into the ``resampling_config`` last.
``thyra_supported``
    Optional, default true. False for a format Thyra cannot read: the
    SCiLS side is still extracted and reported, the comparison is skipped.
``variants``
    Optional, ``{label: {"reader_options": {...}, "method": ...}}``. Each
    is converted and compared too and reported on its own row, so a
    reader option that changes the result can be shown beside Thyra's
    default. Only the main entry is asserted on by the tests.
``scils_import``, ``notes``
    Free text copied into the report.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

#: The environment variable naming the folder of dataset folders.
PARITY_ENV = "THYRA_PARITY_DATA"


@dataclass(frozen=True)
class ParityDataset:
    """One dataset folder.

    Attributes:
        name: The folder name, which is the dataset's public name.
        root: The dataset folder.
        manifest: The parsed ``parity.json``.
    """

    name: str
    root: Path
    manifest: Dict[str, Any] = field(repr=False)

    @property
    def source(self) -> Path:
        """The input Thyra converts."""
        return self.root / self.manifest["source"]

    @property
    def scils_dir(self) -> Path:
        """The folder holding the SCiLS Lab outputs."""
        return self.root / "scils"

    @property
    def thyra_dir(self) -> Path:
        """The folder the harness writes Thyra stores into."""
        return self.root / "thyra"

    @property
    def spectrum(self) -> str:
        """``"profile"`` or ``"centroid"``."""
        return str(self.manifest.get("spectrum", "profile"))

    @property
    def reader_options(self) -> Dict[str, Any]:
        """Options passed to ``convert_msi`` as ``reader_options``."""
        return dict(self.manifest.get("reader_options") or {})

    @property
    def pixel_size_um(self) -> Optional[float]:
        """The pixel size to pass, when the source's own is not read."""
        value = self.manifest.get("pixel_size_um")
        return None if value is None else float(value)

    @property
    def thyra_supported(self) -> bool:
        """Whether Thyra has a reader for this source."""
        return bool(self.manifest.get("thyra_supported", True))

    @property
    def has_reference(self) -> bool:
        """Whether the SCiLS reference has been extracted."""
        return (self.scils_dir / "reference.npz").is_file()

    def variants(self) -> Dict[str, "ParityDataset"]:
        """The manifest's ``variants``, each as a dataset with its overrides applied."""
        out = {}
        for label, override in (self.manifest.get("variants") or {}).items():
            manifest = {k: v for k, v in self.manifest.items() if k != "variants"}
            manifest.update(override)
            out[label] = ParityDataset(
                name=self.name, root=self.root, manifest=manifest
            )
        return out


def parity_root() -> Optional[Path]:
    """The folder named by :data:`PARITY_ENV`, or None when unset or missing."""
    value = os.environ.get(PARITY_ENV, "").strip()
    if not value:
        return None
    root = Path(value)
    return root if root.is_dir() else None


def discover(root: Optional[Path] = None) -> List[ParityDataset]:
    """Every dataset folder under ``root``, sorted by name.

    Args:
        root: The folder to search; defaults to :func:`parity_root`.

    Returns:
        The datasets, or an empty list when there is no root.
    """
    root = parity_root() if root is None else root
    if root is None:
        return []
    found = []
    for manifest_path in sorted(root.glob("*/parity.json")):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        found.append(
            ParityDataset(
                name=manifest_path.parent.name,
                root=manifest_path.parent,
                manifest=manifest,
            )
        )
    return found
