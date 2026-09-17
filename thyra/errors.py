"""The exception a refusal is spelled with.

Thyra refuses a great deal on purpose: a ``.d`` directory with no
analysis files, a Waters raster whose stage never moved, a ``tof`` axis
with no width law, an output path that already holds a store. Each of
those refusals is a sentence someone wrote for the person who hit it,
and each used to reach that person as a thirty-line traceback --
``convert_msi`` caught every exception alike and logged the message
*and* the traceback at ERROR, so a refusal Thyra planned for was
presented exactly like a crash it did not (issue #234).

:class:`ConversionRefused` is how the two are told apart. Raising it
says "this is the whole explanation": the CLI prints the message once at
ERROR and keeps the traceback for ``--log-level DEBUG``. Anything else
reaching the same handler is unexpected, and keeps its traceback at
ERROR, because for those the traceback *is* the explanation.

It subclasses :class:`ValueError` because every one of these sites
raised a ``ValueError`` before, and callers -- Thyra's own tests
included -- catch that. Nothing that used to work stops working.

The distinction is about who the message is for, not about where it is
raised. An internal invariant ("Common mass axis is not initialized")
stays a plain ``ValueError``: nobody can act on it, and its traceback is
the only useful part.
"""


class ConversionRefused(ValueError):
    """Thyra declined to do something, and the message says why.

    The message is addressed to the person who ran the conversion and is
    expected to name what to do about it. No traceback is shown for it
    at the default log level.
    """


#: What a provenance block is allowed to fail with (issue #280).
#:
#: The blocks that carry it are assembled from whatever shape a reader's
#: extractor produced -- a missing attribute, a key the vendor did not
#: write, a value of the wrong type -- so these four are the expected
#: outcome of an unfamiliar source and cost the store one section.
#: Everything else is Thyra breaking its own invariant, and a store missing
#: a section it was asked to write is a worse outcome than a traceback, so
#: the rest propagates.
#:
#: :class:`ConversionRefused` subclasses ``ValueError`` and is therefore
#: *inside* this tuple. That is deliberate: none of the sites that use it
#: calls anything that refuses. A site that does must re-raise the refusal
#: ahead of the catch, because a refusal is addressed to the person who ran
#: the conversion and a log line is not delivery.
#:
#: It lives here, beside that refusal, because both the converter and
#: :mod:`thyra.metadata.uns_assembler` catch it, and neither module is a
#: home for it: the converter imports the assembler, so the assembler
#: cannot import the tuple back without a cycle, and the converter's error
#: policy does not belong inside :mod:`thyra.metadata`.
MALFORMED_METADATA = (AttributeError, KeyError, TypeError, ValueError)


__all__ = ["ConversionRefused", "MALFORMED_METADATA"]
