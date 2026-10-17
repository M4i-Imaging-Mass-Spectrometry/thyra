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


_SPECTRUM_TERMS_INSTALLED = False


def ensure_spectrum_term_fix() -> None:
    """Patch pyimzml to correct a mislabelled spectrum term as it does elsewhere.

    pyimzml corrects a cvParam an exporter wrote under the wrong accession
    ("positive scan" as ``MS:1000128``, the profile spectrum term) where it
    parses terms itself: ``fileContent`` and the parameter groups. The
    per-spectrum values asked for with ``include_spectra_metadata`` are
    looked up by accession alone, so a file stating that pair on each
    spectrum had no polarity at all, while the same pair in a group gave
    "positive" (D55). After the patch, ``spectrum_metadata_fields`` reads as
    if the file were written correctly. Nothing else the parser records
    changes, and a parser built without ``include_spectra_metadata`` is
    untouched.

    The hook is pyimzml's per-spectrum step, which is private.
    ``test_polarity_statements`` converts a file that states the pair on
    each spectrum, so a pyimzml release that renames the step fails there
    rather than losing the polarity without a word. Idempotent.
    """
    global _SPECTRUM_TERMS_INSTALLED
    if _SPECTRUM_TERMS_INSTALLED:
        return

    from pyimzml import ImzMLParser as parser_module

    from ...metadata.constants import MISLABELLED_POLARITY_TERMS

    parser_class = parser_module.ImzMLParser
    # Name-mangled, as pyimzml's own calls of it are.
    step = "_ImzMLParser__process_spectrum"
    original = getattr(parser_class, step)
    prefix = parser_module.XMLNS_PREFIX
    corrections = [
        (f'.//{prefix}cvParam[@accession="{accession}"][@name="{name}"]', fixed)
        for (accession, name), fixed in MISLABELLED_POLARITY_TERMS.items()
    ]

    def process_spectrum(self, elem, include_spectra_metadata):
        original(self, elem, include_spectra_metadata)
        fields = getattr(self, "spectrum_metadata_fields", None)
        if not isinstance(fields, dict):
            return
        for path, fixed in corrections:
            values = fields.get(fixed)
            # The step above appended this spectrum's value last; it is None
            # when the spectrum states no correct term.
            if values and values[-1] is None and elem.find(path) is not None:
                values[-1] = True

    process_spectrum.__wrapped__ = original  # type: ignore[attr-defined]
    setattr(parser_class, step, process_spectrum)
    _SPECTRUM_TERMS_INSTALLED = True
    logger.debug("Installed the spectrum term correction for pyimzml")
