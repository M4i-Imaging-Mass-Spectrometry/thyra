# Coordinate Systems

Every SpatialData zarr Thyra writes carries a small but crucial promise:
**all elements within the zarr resolve to the same physical frame at the
``"global"`` coordinate system, and that frame is documented in the
zarr metadata.** This page explains the contract so consumers (Ousia,
registration tooling, custom analysis scripts) can render and compute
without guessing what the numbers mean.

---

## The contract in one sentence

In every Thyra-produced zarr,

> the ``"global"`` coordinate system is one self-consistent frame, its
> unit and pixel-size calibration are recorded under
> ``zarr.attrs["coordinate_systems"]["global"]``, and every element
> registers a transform that lands in that frame.

If you only remember one thing, remember that.

---

## The schema

At the zarr top level Thyra writes:

```python
zarr.attrs["coordinate_systems"] = {
    "global": {
        "unit": "micrometer" | "pixel",
        "pixel_size_um_x": float | None,
        "pixel_size_um_y": float | None,
        "reference_element": str | None,
        "convention_version": 1,
        "produced_by": "thyra/<version>",
    }
}
```

| Field | Meaning |
|-------|---------|
| ``unit`` | Unit of one step in ``"global"`` -- either ``"micrometer"`` or ``"pixel"``. |
| ``pixel_size_um_x``, ``pixel_size_um_y`` | Conversion factor: micrometers per one ``"global"`` unit. ``1.0`` when ``unit="micrometer"``; physical pixel size when ``unit="pixel"``. ``None`` if the producer cannot calibrate (e.g. an uncalibrated optical photo). |
| ``reference_element`` | Name of the canonical raster element that defines pixel space, when ``unit="pixel"``. ``None`` otherwise. |
| ``convention_version`` | Schema version; bump when the shape of this attr changes. Currently ``1``. |
| ``produced_by`` | ``"thyra/<version>"`` for Thyra-produced zarrs. |

Consumers that don't know the schema can still open the zarr; they
just won't be able to render distances in micrometers without
making assumptions. Consumers that do know it (Ousia, the EscDat
registration script, this project's tests) read this attr first and
trust it.

---

## What ``"global"`` is in each Thyra mode

Thyra operates in two modes depending on whether the input came with
FlexImaging optical alignment data. The two modes pick different
``"global"`` conventions because they have different *canonical
references*.

### Mode A: standalone (no optical alignment)

When there is no optical photo to align against, the only naturally
meaningful frame is **physical micrometers of the imaged tissue**.

- TIC / ion images: stored intrinsically in raster pixel indices,
  with a transform ``Scale([pixel_size_um, pixel_size_um])`` to
  ``"global"``.
- Pixel-polygon shapes: stored intrinsically in micrometers
  (``spatial_x = x * pixel_size_um``), with ``Identity`` to
  ``"global"``.
- ``zarr.attrs["coordinate_systems"]["global"]`` declares
  ``unit="micrometer"`` and fills ``pixel_size_um_x/y`` with the
  MSI grid pixel size.

Both elements resolve to the same micrometer extent at ``"global"``.

### Mode B: with FlexImaging optical alignment

When the input includes a ``.mis`` ``ImageFile`` and registration
landmarks, Thyra picks a different canonical reference: **the
optical photo itself**. The natural frame is then *that image's
pixel grid*, because the photo can be drawn with no transform.

- TIC / ion images: stored intrinsically in raster pixel indices,
  with an ``Affine`` transform mapping into the optical pixel
  grid.
- Pixel-polygon shapes: stored directly in optical pixels, with
  ``Identity`` to ``"global"``.
- Primary optical image: ``Identity`` to ``"global"``.
- Other optical images (overview, etc.): a ``Scale`` to align with
  the primary image's pixel grid.
- ``zarr.attrs["coordinate_systems"]["global"]`` declares
  ``unit="pixel"`` with ``reference_element`` set to the primary
  optical filename; ``pixel_size_um_x/y`` is typically ``None``
  because FlexImaging photos do not generally carry a um-per-pixel
  calibration.

Note that the two modes pick the right convention for what is
actually known about the data; consumers should look at
``unit`` rather than assuming Thyra always uses one or the other.

---

## Reading the coordinate-system metadata

```python
import json
from pathlib import Path

zarr_path = Path("output.zarr")
with open(zarr_path / "zarr.json") as f:
    attrs = json.load(f)["attributes"]

cs = attrs["coordinate_systems"]["global"]
print(f"global unit:        {cs['unit']}")
print(f"pixel size (um):    {cs['pixel_size_um_x']} x {cs['pixel_size_um_y']}")
print(f"reference element:  {cs['reference_element']}")
print(f"schema version:     {cs['convention_version']}")
```

To resolve a ``"global"`` coordinate to micrometers, multiply by the
``pixel_size_um_x/y`` factors. When ``unit="micrometer"`` those
factors are ``1.0`` and the multiplication is a no-op; when
``unit="pixel"`` they perform the px-to-um conversion.

---

## Why two conventions?

Picking ``unit="micrometer"`` for everything would be simpler, but
also wrong:

- In **Mode B** there is a canonical raster image (the optical
  photo). Putting ``"global"`` in micrometers would force that
  image to carry a transform, since its data is intrinsically in
  pixels. With ``unit="pixel"``, the canonical image carries
  ``Identity`` and downstream tools can blit it without thinking.
- In **Mode A** there is no canonical raster image - the MSI
  itself is the data. There's no privileged pixel grid to use as
  ``"global"``, so ``unit="micrometer"`` is the natural and
  unambiguous choice.

The rule generalises: **``"global"`` is whichever frame is most
natural for that dataset's canonical reference, and the metadata
declares which one was picked.**

This matches how the broader spatial-omics ecosystem works:
``spatialdata-io``'s Xenium reader, for example, uses
``unit="pixel"`` because Xenium has a canonical morphology image;
its MERFISH-style readers use ``unit="micrometer"`` because they
do not.

---

## Verifying the contract

Because the contract is just metadata, it can drift if a producer
has a bug. Every Thyra release has a unit-test guard that:

1. Runs a representative conversion end-to-end.
2. Reads the resulting zarr.
3. Asserts the bbox of the TIC image and the pixel-polygon shapes
   resolve to the same ``"global"`` extent within a half-pixel
   tolerance.
4. Asserts the ``coordinate_systems.global`` attr is present and
   self-consistent.

See ``tests/unit/converters/test_coordinate_systems.py``.

A consumer can run the same check at load time. The recommended
pattern is to compute the bbox at ``"global"`` for every element
and warn if any pair differs by more than ~10x; that catches the
``Scale(1/pixel_size)`` vs identity unit-confusion class of bug
without flagging legitimate small-vs-large element pairs (e.g.
a small ROI within a full-tissue dataset).

---

## Cross-modality registration

When a Thyra zarr is later registered against another modality
(e.g. Xenium via the EscDat registration pipeline), the
registration step rewrites every Thyra element's transform to
``"global"`` so that the merged zarr's ``"global"`` matches the
other modality's convention. The merged zarr's
``coordinate_systems.global`` attr is rewritten to declare the new
unit and reference element.

That is, **the coordinate-system contract is per-zarr, not
per-element-or-modality**. A standalone Thyra zarr can have
``unit="micrometer"`` while a merged Thyra+Xenium zarr produced
from the same data will have ``unit="pixel"``. Both are correct
for their respective zarrs.

---

## What if I open a zarr that doesn't have this attr?

Older Thyra outputs (pre-schema) do not carry the
``coordinate_systems`` attr. They are still readable, but you have
to infer the convention from element layouts and accept some
risk that two elements disagree silently.

If you control the data: **regenerate**. The schema is cheap to
write and the resulting zarrs are self-describing forever after.
If you do not control the data, the best fallback is to write a
producer-specific heuristic (see Ousia's coordinate-system loader
for an example: it recognises Xenium-from-spatialdata-io zarrs by
their ``morphology_focus`` image and assumes Xenium's documented
pixel size).
