"""Counting the spectra a read drops, instead of narrating each one.

Every reader's frame loop announced each dropped spectrum with its own
``logger.warning`` and counted none of them. On a 17,423-pixel
acquisition where the source is systematically unreadable -- a truncated
``.tdf_bin``, a region whose ``MaldiFrameInfo`` rows are missing -- that
is 17,423 warning lines, twice, because the converter reads the source
in two passes. The one number a user needs, *how much of my data did not
make it*, was in none of them; the run then either succeeded with a
short store or hit the empty-conversion refusal, neither of which names
the drops as the cause.

So: name the first one, count them all, and say the total once when the
pass is over. That is the shape the sibling tables already use
(``msms_table``'s ``n_dropped``, ``mobility_table``'s merged-feature and
empty-row counts) -- this is the reader-side version of it, not a new
convention.

The first warning is kept in full because it carries the *cause*, which
is what makes the total actionable; the rest of the causes are logged at
DEBUG, where a user chasing a specific frame can still find them.
"""

import logging
from typing import Optional

__all__ = ["DropTally"]


class DropTally:
    """Warn once with the cause, count the rest, summarise at the end.

    Args:
        logger: The module logger to report through, so lines keep the
            name of the code that dropped the spectrum.
        what: What is being dropped, as it should read mid-sentence --
            e.g. ``"frames with no coordinates"``.

    Example:
        >>> tally = DropTally(logger, "frames with no coordinates")
        >>> for frame_id in frame_ids:            # doctest: +SKIP
        ...     coords = lookup(frame_id)
        ...     if coords is None:
        ...         tally.drop(f"frame {frame_id}")
        ...         continue
        ...     yield frame_id, coords
        >>> tally.summarise(len(frame_ids))       # doctest: +SKIP
    """

    def __init__(self, logger: logging.Logger, what: str) -> None:
        """Start an empty tally."""
        self._logger = logger
        self._what = what
        self.n_dropped = 0

    def drop(self, detail: str, cause: Optional[BaseException] = None) -> None:
        """Record one dropped spectrum.

        Args:
            detail: What was dropped, e.g. ``"frame 41"``.
            cause: The exception, when there was one.
        """
        self.n_dropped += 1
        suffix = f": {cause}" if cause is not None else ""
        if self.n_dropped == 1:
            self._logger.warning(
                "Dropping %s%s. Further occurrences are logged at DEBUG; "
                "the total is reported when the pass ends.",
                detail,
                suffix,
            )
        else:
            self._logger.debug("Dropping %s%s", detail, suffix)

    def summarise(self, n_total: Optional[int] = None) -> None:
        """Report the total once, if anything was dropped.

        Args:
            n_total: How many were attempted, when it is known. The
                fraction is what tells a user whether they lost a corner
                or the acquisition.
        """
        if not self.n_dropped:
            return
        if n_total:
            self._logger.warning(
                "Dropped %d of %d %s", self.n_dropped, n_total, self._what
            )
        else:
            self._logger.warning("Dropped %d %s", self.n_dropped, self._what)
