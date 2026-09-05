"""Make pyimzml survive a cvParam that has a numeric type but no value.

``pyimzml.ontology.ontology.convert_xml_value`` converts a cvParam's value
to the type its ontology entry declares. It catches ``ValueError`` -- an
empty ``value=""`` on a float-typed term becomes ``None`` -- but not
``TypeError``, so a term with **no** ``value`` attribute at all
(``float(None)``) aborts the whole metadata parse, before a single spectrum
is read.

That is exactly what the pyimzML fork behind TIMSCONVERT and TIMSImaging
writes for the ion mobility array it adds: a ``referenceableParamGroup``
with ``MS:1003006 mean inverse reduced ion mobility array`` carrying a unit
and no value, and pyimzml's table types ``MS:1003006`` as ``xsd:float``.
Upstream pyimzml therefore cannot open any such export. Measured on the
hand-authored ``mobility_continuous.imzML``: ``TypeError: float() argument
must be a string or a real number, not 'NoneType'`` from
``ImzMLParser.__init__``.

The patch is as narrow as the defect: a missing value on a typed term
converts to ``None``, the same result an empty value already produces.
Every other conversion is untouched. It is installed once, the first time
Thyra builds a parser, and is idempotent.
"""

import logging

logger = logging.getLogger(__name__)

_INSTALLED = False


def ensure_lenient_cv_param_values() -> None:
    """Patch pyimzml so a value-less typed cvParam parses as ``None``."""
    global _INSTALLED
    if _INSTALLED:
        return

    from pyimzml import metadata as pyimzml_metadata
    from pyimzml.ontology import ontology as pyimzml_ontology

    original = pyimzml_ontology.convert_xml_value

    def convert_xml_value(dtype, value):
        try:
            return original(dtype, value)
        except TypeError:
            if value is None:
                return None
            raise

    convert_xml_value.__doc__ = original.__doc__
    convert_xml_value.__wrapped__ = original  # type: ignore[attr-defined]

    pyimzml_ontology.convert_xml_value = convert_xml_value
    # metadata.py imports the function by name, so its own binding must be
    # replaced too; lookup_and_convert_cv_param resolves it through the
    # ontology module's globals, which the assignment above already covers.
    if getattr(pyimzml_metadata, "convert_xml_value", None) is original:
        pyimzml_metadata.convert_xml_value = convert_xml_value

    _INSTALLED = True
    logger.debug("Installed the lenient cvParam value conversion for pyimzml")
