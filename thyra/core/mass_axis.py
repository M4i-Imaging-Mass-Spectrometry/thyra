"""Building a raw mass axis: the union of every spectrum's m/z values.

Five readers needed this and five wrote it, four of them badly (issue
#294). The shapes were:

- **imzML** folded each batch into a running axis, so its transient was
  O(unique). This module is that implementation, lifted.
- **Waters** and **solariX** appended every spectrum's array to a list,
  then paid ``np.concatenate`` and ``np.unique`` -- three payload-sized
  allocations and a transient of O(*total peaks*), which is the whole
  file rather than its distinct values. Worst exactly when the axis is
  small: a shared profile grid (the Waters default on an MRT) has few
  unique values and many spectra, so the concatenate is pure waste.
- **mzPeak** folded with ``np.union1d`` per row group, which re-copies
  the whole axis every group -- O(unique) memory but quadratic work.
- **timsTOF** accumulated into a ``set[float]``, which is O(unique) in
  elements but ~10x that in bytes: every value becomes a boxed CPython
  float (24 bytes, 32 with pymalloc's block rounding), plus the set
  table, plus ``sorted()``'s pointer list. It reported its own cost as
  ``len(unique) * 8``, the output array alone.

Nothing downstream could make up for it. The converter refuses a
too-wide axis, but only *after* ``get_common_mass_axis()`` has returned,
so a build whose transient is O(total) has already peaked by then.

## ``finish()`` takes no argument, and that is load-bearing

The accumulator formats a refusal naming how far the build got, and in
its first home that index arrived as a parameter: ``add(mzs, idx)`` and
``finish(idx)``. Its one caller seeded ``idx = -1`` before the loop so
the call was safe on an empty source.

Every new caller would have had to know that. Three of the four readers
adopting it here loop over something that can be empty -- Waters over
functions and scans, solariX over a cursor, mzPeak over row groups -- and
the natural ``acc.finish(scan)`` is an ``UnboundLocalError`` on a source
with no spectra, which is the shape a reader is *most* likely to be
handed by a broken file. Measured on all three before the change.

So the count lives in the accumulator: ``add()`` increments it, and
``finish()`` needs nothing. The hazard is gone structurally rather than
by asking four call sites to remember a sentinel.
"""

import logging
from typing import Any, List, Optional

import numpy as np
from numpy.typing import NDArray

from ..errors import ConversionRefused

logger = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_MAX_MASS_AXIS_LENGTH",
    "MASS_AXIS_BATCH_VALUES",
    "MASS_AXIS_PROBE_BLOCK",
    "MassAxisAccumulator",
    "dedupe_sorted",
    "filter_absent",
    "merge_disjoint",
    "validate_max_mass_axis_length",
]


# 1Mi elements is 8 MiB at float64. Measured at a 512 MiB m/z payload on a
# shared grid (peak extra RSS / wall clock): 2^14 -> 7.14 MiB / 3.04 s,
# 2^17 -> 6.90 / 3.06, 2^20 -> 15.30 / 2.30, 2^22 -> 40.29 / 2.11,
# 2^24 -> 148.38 / 2.06. A 1024x sweep of the knob moves wall clock 1.5x with
# no cliff anywhere, so this is a knee rather than a fragile constant.
MASS_AXIS_BATCH_VALUES = 1 << 20

# Block size for the searchsorted probe in ``filter_absent``, in ELEMENTS.
# Bounds that function's temporaries to O(block) rather than O(batch).
MASS_AXIS_PROBE_BLOCK = 1 << 20

# Default cap on the number of unique m/z values the processed-mode raw mass
# axis may reach, in ELEMENTS. Taken from SCiLS Lab, which applies the same
# limit to the same quantity: "Data sets in SCiLS Lab are limited to a maximum
# of 10 million bins on the common mass axis" (SCiLS Lab 2026b User Guide,
# p.76). The cap previously defaulted to None because no defensible number was
# available; this is one, from the tool Thyra's resampling was designed
# against. Override with ``reader_options={"max_mass_axis_length": N}``, or
# pass None for the old unlimited behaviour.
DEFAULT_MAX_MASS_AXIS_LENGTH = 10_000_000


def dedupe_sorted(a: NDArray[Any]) -> NDArray[Any]:
    """Deduplicate an ALREADY-SORTED array, as ``np.unique`` would.

    Reproduces numpy's ``_unique1d`` mask including its ``equal_nan=True``
    behaviour, which collapses a trailing run of NaNs to its first element.
    Unlike ``np.unique`` it does not make the extra full-size ``flatten()``
    copy, which is one of the three payload-sized allocations the old
    collect-everything-then-unique build paid for.

    Args:
        a: A sorted array.

    Returns:
        The distinct values of ``a``, in order.
    """
    n = a.size
    if n == 0:
        return a[:0].copy()

    mask = np.empty(n, dtype=bool)
    mask[0] = True

    if a.dtype.kind == "f" and bool(np.isnan(a[-1])):
        first_nan = int(np.searchsorted(a, a[-1], side="left"))
        # The ``first_nan > 0`` guard is load-bearing: an all-NaN input makes
        # the slices empty and ``np.not_equal`` raises on the shape mismatch.
        if first_nan > 0:
            np.not_equal(a[1:first_nan], a[: first_nan - 1], out=mask[1:first_nan])
        mask[first_nan] = True
        mask[first_nan + 1 :] = False
    else:
        np.not_equal(a[1:], a[:-1], out=mask[1:])

    return a[mask]


def filter_absent(
    acc: NDArray[Any],
    s: NDArray[Any],
    block: int = MASS_AXIS_PROBE_BLOCK,
) -> NDArray[Any]:
    """Return the elements of ``s`` that are not already in ``acc``.

    Both arrays must be sorted and free of duplicates. ``s`` is walked in
    blocks so the search temporaries stay O(block) rather than O(len(s)).
    Returns ``s`` itself, uncopied, when none of its values are present --
    the all-distinct fast path, worth half a payload at peak.

    Args:
        acc: The running axis, sorted and unique.
        s: A candidate batch, sorted and unique.
        block: Probe block size in elements.

    Returns:
        The subset of ``s`` absent from ``acc``.
    """
    n = s.size
    # ``acc.size == 0`` must be handled here: the ``np.minimum`` below would
    # otherwise index ``acc[-1]`` on an empty array and raise IndexError.
    if n == 0 or acc.size == 0:
        return s

    keep = np.empty(n, dtype=bool)
    n_keep = 0
    is_float = acc.dtype.kind == "f"
    a_len = acc.size

    for lo in range(0, n, block):
        hi = min(lo + block, n)
        sb = s[lo:hi]
        pos = np.searchsorted(acc, sb, side="left")
        out_of_range = pos >= a_len
        np.minimum(pos, a_len - 1, out=pos)
        va = acc[pos]
        eq = va == sb
        if is_float:
            # searchsorted routes a NaN key onto acc's first NaN, but NaN does
            # not compare equal to itself, so without this the same NaN would
            # be re-inserted on every fold.
            np.logical_or(eq, np.isnan(va) & np.isnan(sb), out=eq)
        np.logical_and(eq, ~out_of_range, out=eq)
        np.logical_not(eq, out=eq)
        keep[lo:hi] = eq
        n_keep += int(np.count_nonzero(eq))
        del pos, out_of_range, va, eq, sb

    if n_keep == n:
        return s
    if n_keep == 0:
        return s[:0]
    return s[keep]


def merge_disjoint(box_a: List[Any], box_b: List[Any]) -> NDArray[Any]:
    """Merge two DISJOINT sorted arrays, releasing each once it is copied.

    The inputs arrive boxed in one-element lists so this function can drop the
    caller's reference as soon as each side has been copied. Passing them as
    plain parameters would keep both alive for the whole call and add a full
    copy of the larger side to the peak.

    The result needs no second deduplication because ``filter_absent`` has
    already removed the overlap, which saves another full-size copy.

    Args:
        box_a: One-element list holding the first sorted array.
        box_b: One-element list holding the second sorted array.

    Returns:
        The sorted union of the two inputs.
    """
    a, b = box_a[0], box_b[0]
    out = np.empty(a.size + b.size, dtype=a.dtype)

    n_a = a.size
    out[:n_a] = a
    box_a[0] = None
    del a

    out[n_a:] = b
    box_b[0] = None
    del b

    # Exactly two ascending runs, which is the case quicksort handles best.
    out.sort(kind="quicksort")
    return out


class MassAxisAccumulator:
    """Builds a sorted, unique m/z axis without holding every spectrum.

    Values are copied into a scratch buffer; when that buffer fills, it is
    sorted, deduplicated, and merged into the running axis. The buffer's
    capacity tracks the axis length, so the number of merges is logarithmic in
    the input rather than linear -- a fixed batch size would re-copy the
    accumulator once per batch, which is the quadratic behaviour that makes a
    naive ``np.union1d`` fold far slower than the code it replaces.

    Peak memory therefore depends on the number of *unique* m/z values rather
    than on the size of the file.
    """

    def __init__(
        self, total_spectra: Optional[int] = None, max_length: Optional[int] = None
    ) -> None:
        """Initialize the accumulator.

        Args:
            total_spectra: Spectrum count, used only in the error message,
                and ``None`` when the reader cannot say cheaply -- solariX
                is walking a cursor and would need a second ``COUNT(*)``
                over the table to answer. The message then reports how far
                the build got without claiming a denominator it does not
                have.
            max_length: Cap on unique m/z values, or None for unlimited.
        """
        self._total_spectra = total_spectra
        self._max_length = max_length
        self._acc: Optional[NDArray[Any]] = None
        self._buf: Optional[NDArray[Any]] = None
        self._cap = 0
        self._cap_target = MASS_AXIS_BATCH_VALUES
        self._n = 0
        self.saw_any = False
        #: Spectra handed to :meth:`add`, for the refusal message. Tracked
        #: here so no caller has to keep a loop variable alive for it.
        self.n_seen = 0

    @property
    def n_unique(self) -> int:
        """Distinct m/z values folded in so far.

        Only what has been *folded*: values still sitting in the scratch
        buffer are not counted, so this is a lower bound during a build and
        exact once :meth:`finish` has run. It is for progress reporting,
        which is the only thing that asks mid-build.
        """
        return 0 if self._acc is None else int(self._acc.size)

    def add(self, mzs: NDArray[Any]) -> None:
        """Buffer one spectrum's m/z values, folding first if the buffer is full.

        Args:
            mzs: The spectrum's m/z values. Not modified.
        """
        m = mzs.size
        self.n_seen += 1
        if not m:
            return
        self.saw_any = True

        if self._buf is None:
            self._allocate(m, mzs.dtype)
        elif self._n + m > self._cap:
            self._fold()
            if self._buf is None:
                self._allocate(m, mzs.dtype)
            elif m > self._cap:
                self._buf = None
                self._allocate(m, mzs.dtype, exact=True)

        assert self._buf is not None
        self._buf[self._n : self._n + m] = mzs
        self._n += m

    def finish(self) -> NDArray[Any]:
        """Fold whatever is buffered and return the axis.

        Returns:
            Sorted, unique m/z values.

        Raises:
            ConversionRefused: If no spectrum yielded any m/z values.
        """
        self._fold()
        self._buf = None

        if not self.saw_any or self._acc is None:
            raise ConversionRefused("No spectra found to build common mass axis")
        if self._acc.size == 0:
            raise ConversionRefused("Failed to extract any m/z values")
        return self._acc

    def _allocate(self, m: int, dtype: Any, exact: bool = False) -> None:
        """Allocate the scratch buffer, at least large enough for one spectrum."""
        self._cap = m if exact else max(self._cap_target, m)
        self._buf = np.empty(self._cap, dtype=dtype)

    def _fold(self) -> None:
        """Sort, deduplicate and merge the buffered batch into the axis."""
        if self._n == 0 or self._buf is None:
            return

        view = self._buf[: self._n]
        view.sort(kind="quicksort")
        run = dedupe_sorted(view)
        self._n = 0

        if self._cap > MASS_AXIS_BATCH_VALUES:
            # The scratch has grown to axis size, so release it before the
            # merge allocates. Gating on the batch floor keeps this from firing
            # on every one of the thousands of folds a shared-grid dataset
            # performs, where re-faulting the buffer each time costs more than
            # the memory it frees.
            view = None
            self._buf = None
            self._cap = 0

        if self._acc is None:
            self._acc = run
        else:
            new = filter_absent(self._acc, run)
            run = None  # drop the batch before the merge allocates
            if new.size:
                box_a, box_b = [self._acc], [new]
                self._acc = None
                new = None
                self._acc = merge_disjoint(box_a, box_b)

        self._check_limit()

        # Record the next capacity but do NOT allocate it here. On the final
        # fold the stream is already exhausted, so allocating would commit a
        # whole axis-sized buffer that is never written.
        self._cap_target = max(MASS_AXIS_BATCH_VALUES, self._acc.size)
        if self._cap and self._cap != self._cap_target:
            self._buf = None
            self._cap = 0

    def _check_limit(self) -> None:
        """Stop if the axis has outgrown ``max_length``."""
        if self._max_length is None or self._acc is None:
            return
        if self._acc.size <= self._max_length:
            return
        progress = (
            f"after {self.n_seen:,} of {self._total_spectra:,} spectra"
            if self._total_spectra is not None
            else f"after {self.n_seen:,} spectra"
        )
        raise ConversionRefused(
            f"Common mass axis exceeded {self._max_length:,} unique m/z values "
            f"{progress} "
            f"({self._acc.size:,} so far). The peak lists in this dataset do "
            "not share m/z values, so a raw axis grows to roughly one column "
            "per peak, which is not usable downstream. Convert with resampling "
            "instead (it is the default; --no-resample disables it), or raise "
            "max_mass_axis_length."
        )


def validate_max_mass_axis_length(value: Any) -> Optional[int]:
    """Check the raw-axis cap, which is a count of unique m/z values.

    Refused here rather than at extraction time, so a bad value fails while
    the caller is still looking at their own arguments -- the same reason
    ``spectrum_type`` is normalised in the constructor. Until then ``-1``,
    ``0`` and ``2.5`` were all accepted and then compared against a growing
    axis, so every file was refused for "exceeding" a cap it could not
    possibly satisfy, in a message that named the nonsense value back at the
    person who set it (issue #261).

    ``bool`` is rejected explicitly: it is an ``int`` subclass, so
    ``max_mass_axis_length=True`` would otherwise mean a cap of one.

    Args:
        value: What the caller passed, or the default.

    Returns:
        The value unchanged, once it is ``None`` or a positive int.

    Raises:
        ConversionRefused: On anything else.
    """
    if value is None or (
        isinstance(value, int) and not isinstance(value, bool) and value > 0
    ):
        return value

    raise ConversionRefused(
        f"max_mass_axis_length must be a positive integer or None (no limit), "
        f"got {value!r}. It caps how many unique m/z values a processed-mode "
        f"raw axis may reach before the build gives up; the default is "
        f"{DEFAULT_MAX_MASS_AXIS_LENGTH:,}."
    )
