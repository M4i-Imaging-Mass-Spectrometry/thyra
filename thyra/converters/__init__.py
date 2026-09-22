"""Converter plugins for MSI data output formats.

Importing this package is what pulls in ``spatialdata``, and through it
``dask`` and ``anndata``. It used to happen during ``import thyra``, for
the ``@register_converter`` side effect; the registry now imports
:mod:`thyra.converters.spatialdata.converter` the first time the
``spatialdata`` output format is resolved (issue #381), so this runs when
a conversion is actually being set up.

The Dask configuration below moved here from ``thyra/__init__.py`` for the
same reason: reading it there meant importing ``dask`` to silence a
``dask`` warning, on every import of anything under ``thyra``. It has to
be set before ``dask.dataframe`` is first imported, which is why it comes
before the subpackage import rather than after.
"""

try:
    import dask

    # Configure Dask to use new query planning (silences legacy DataFrame
    # warning)
    dask.config.set({"dataframe.query-planning": True})
except ImportError:
    pass

from . import spatialdata  # noqa: E402,F401
