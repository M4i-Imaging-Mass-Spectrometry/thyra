"""Construct an ``ImzMLParser`` the way production does.

Production has exactly one ``ImzMLParser(`` construction site,
``ImzMLReader._initialize_parser``, and it passes ``ibd_file=`` an already-open
handle. Omitting the argument sends pyimzml down ``_infer_bin_filename``
instead, which is a different rule: it globs the parent directory and keeps
anything matching an unanchored ``.+\\.ibd`` whose ``Path.stem`` compares equal
to the imzML's, case-sensitively.

The two rules disagree in three ways that matter:

* **Case.** pyimzml's ``re.IGNORECASE`` is decorative, because the stem compare
  beside it is a plain ``==``. ``SAMPLE.IBD`` next to ``sample.imzML`` raises
  ``IndexError: list index out of range`` in pyimzml and resolves fine through
  ``with_suffix('.ibd').exists()``, which is case-insensitive on Windows.
* **Directories.** The regex is matched against the full path, so a *parent*
  directory named ``proj.ibd`` makes every same-stem sibling eligible --
  including the imzML itself.
* **Partial downloads.** The same unanchored match accepts ``sample.ibdtmp``.

Thyra's rule is the more robust of the two; see ``docs/imzml-parser-notes.md``
before changing either. What matters here is that tests exercising the parser
should exercise the rule production uses, otherwise no test can catch a
resolution bug on the path that ships.
"""

from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Union

from pyimzml.ImzMLParser import ImzMLParser


@contextmanager
def production_parser(imzml_path: Union[str, Path]) -> Iterator[ImzMLParser]:
    """Yield a parser built exactly as ``ImzMLReader._initialize_parser`` builds it.

    Mirrors production on both counts that have bitten before: the ``.ibd`` is
    opened here and handed over as ``ibd_file=``, and ``parse_lib`` is pinned to
    ``ElementTree`` rather than left to pyimzml's default.

    Args:
        imzml_path: Path to the imzML. Its ``.ibd`` sibling must exist.

    Yields:
        The open parser. The handle is closed on the way out, which the
        callers that used to construct a bare parser did not do.
    """
    imzml_path = Path(imzml_path)
    ibd_path = imzml_path.with_suffix(".ibd")
    if not ibd_path.exists():
        raise ValueError(f"Corresponding .ibd file not found for {imzml_path}")

    ibd_file = open(ibd_path, mode="rb")
    try:
        yield ImzMLParser(
            filename=str(imzml_path),
            parse_lib="ElementTree",
            ibd_file=ibd_file,
        )
    finally:
        ibd_file.close()
