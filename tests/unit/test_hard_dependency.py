"""A declared dependency's absence must abort, with the real ImportError.

``pyproject.toml`` declares ``spatialdata`` and ``defusedxml`` in
``[project] dependencies`` with no ``[project.optional-dependencies]``
table anywhere, so an install without either is broken, not a supported
configuration.

spatialdata used to be treated as optional in three independent places,
and the one that fired swallowed the cause: ``import thyra`` succeeded,
every SpatialData name in the base converter was rebound to ``None``, the
converter was never registered, and the first thing the user saw was a
registry miss naming no package and no reason. The real ``ImportError``
sat in a module-level string one module away and reached the user on no
path (issue #310).

defusedxml is the same shape of mistake with a sharper edge: the fallback
was to ``xml.etree``, which is the parser the defusedxml swap was made to
get away from, so a warning nobody reads was all that separated a
hardened parse from an entity-expanding one.

What a broken install is entitled to, in both cases, is the exception and
its traceback, pointing at the import that failed. That is what these
tests hold.

**Where it surfaces moved under issue #381, and the exception did not.**
The fix for #310 was to import both unconditionally during ``import
thyra``, which is exactly the 3.5 s that made the metadata layer unusable
on its own: importing ``thyra.metadata.schema`` initialises ``thyra``
first, so requiring ``spatialdata`` there would have meant an
implementation that only writes the metadata document still had to have
the imaging stack installed to import a pydantic model. So neither package
is checked at the front door any more. ``spatialdata`` is imported by
:mod:`thyra.converters`, which a conversion reaches while resolving its
output format, before it opens the input; ``defusedxml`` by
``thyra.readers.bruker.mis_parser``, which every Bruker path reaches
during format detection. Both raise where they are imported. What #310
actually forbade -- a swallowed ImportError, a name rebound to ``None``, a
registry miss standing in for a missing package -- is still gone, and the
last test here is what says so.
"""

import subprocess
import sys

# The child installs this finder ahead of every other one and refuses one
# root package -- ``fullname.split(".")[0]`` keeps anndata, geopandas,
# shapely and thyra's own ``thyra.converters.spatialdata`` subpackage
# importable, so what is measured is the one dependency.
_BLOCK_TEMPLATE = """
import importlib.abc
import sys


class _Blocker(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] == "{package}":
            raise ImportError("blocked by tests/unit/test_hard_dependency.py")
        return None


sys.meta_path.insert(0, _Blocker())

import thyra

{body}
"""

# Resolving the output format is what a conversion does before it opens the
# input, and it is what imports the converter package. ``print`` is reached
# only if the import it depends on did not raise.
_BLOCK_SPATIALDATA = _BLOCK_TEMPLATE.format(
    package="spatialdata",
    body=(
        "from thyra.core.registry import get_converter_class\n"
        'print("resolving the converter SUCCEEDED", '
        'get_converter_class("spatialdata"))'
    ),
)

# The witness for defusedxml is the module object the parser binds ``ET``
# to: on the fallback the name existed and was ``xml.etree.ElementTree``,
# so printing it is what distinguishes "imported defusedxml" from
# "imported something". Format detection is the path that reaches it -- a
# ``.d`` directory sends the registry to BrukerFolderStructure, whose
# module imports the ``.mis`` parser.
_BLOCK_DEFUSEDXML = _BLOCK_TEMPLATE.format(
    package="defusedxml",
    body=(
        "from thyra.readers.bruker import folder_structure\n"
        'print("importing the folder structure SUCCEEDED", '
        "folder_structure.parse_mis_file)"
    ),
)

# Neither package is needed to describe a dataset, and after #381 neither
# is loaded to do it. The schema is the import an outside implementation of
# the metadata document would make.
_BLOCK_BOTH_FOR_THE_SCHEMA = """
import importlib.abc
import sys


class _Blocker(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] in ("spatialdata", "defusedxml"):
            raise ImportError("blocked by tests/unit/test_hard_dependency.py")
        return None


sys.meta_path.insert(0, _Blocker())

import thyra
from thyra.metadata.schema import MSIMetadata, validate_document

print("OK", MSIMetadata.__name__, validate_document.__name__)
"""


def _run_child(source: str) -> "subprocess.CompletedProcess[str]":
    """Import thyra in a fresh interpreter that cannot import one package.

    A child is required: ``thyra`` and its dependencies are all already in
    this process's ``sys.modules``, so a finder installed here would never
    be consulted for any of them.
    """
    return subprocess.run(
        [sys.executable, "-c", source],
        capture_output=True,
        text=True,
    )


def test_resolving_the_converter_fails_when_spatialdata_is_missing():
    """The conversion must die at the lookup, not warn and carry on.

    The exit status is the load-bearing assertion. stderr alone does not
    discriminate: before the #310 fix the child also wrote "SpatialData
    dependencies not available: ..." to stderr, because the swallowing
    handler logged a warning and logging's lastResort handler printed it --
    so the substring "spatialdata" was present either way. What changed is
    that the process dies.
    """
    proc = _run_child(_BLOCK_SPATIALDATA)

    assert proc.returncode != 0, (
        "the spatialdata converter resolved without spatialdata; the "
        "absence of a hard dependency was swallowed.\n"
        f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
    )
    assert "resolving the converter SUCCEEDED" not in proc.stdout


def test_the_failure_names_the_import_that_failed():
    """The user gets the real exception, not a paraphrase of it.

    "ImportError" is the discriminator here: the message the swallowing
    path produced ("SpatialData dependencies not available", later
    "SpatialData converter unavailable") named neither the exception type
    nor the import site. So is the absence of "No converter for format" --
    a registry miss is what #310 was filed about, and the registry must not
    turn a broken install back into one now that it is the thing doing the
    importing.
    """
    proc = _run_child(_BLOCK_SPATIALDATA)

    assert "ImportError" in proc.stderr, proc.stderr
    assert "spatialdata" in proc.stderr, proc.stderr
    # The traceback points at the import statement, which is the part of a
    # broken install that is actually actionable.
    assert "Traceback" in proc.stderr, proc.stderr
    assert "No converter for format" not in proc.stderr, proc.stderr


def test_bruker_detection_fails_when_defusedxml_is_missing():
    """No silent downgrade to the parser that expands entities.

    ``mis_parser`` carried ``try: import defusedxml ... except ImportError:
    import xml.etree.ElementTree``, which logged a warning and carried on.
    Measured on a ``.mis`` holding ``<!ENTITY r "5,5">``: the stdlib parser
    returns the expansion and reports a 5x5 raster, defusedxml raises
    EntitiesForbidden. So the fallback was not a reduced-function parse, it
    was the unhardened one -- reachable on any install where the wheel
    happened not to be there, announced by a single WARNING line during
    ``import thyra`` that nothing is watching.

    Every Bruker path goes through this parser (the Rapiflex, timsTOF and
    solariX readers, and BrukerMetadataExtractor) and reaches it through
    ``folder_structure``, which the registry imports to tell one Bruker
    format from another -- so the abort is at detection, before a single
    spectrum is read.
    """
    proc = _run_child(_BLOCK_DEFUSEDXML)

    assert proc.returncode != 0, (
        "the .mis parser imported without defusedxml; it fell back to the "
        "entity-expanding stdlib parser.\n"
        f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
    )
    assert "importing the folder structure SUCCEEDED" not in proc.stdout
    assert "ImportError" in proc.stderr, proc.stderr
    assert "defusedxml" in proc.stderr, proc.stderr
    assert "mis_parser" in proc.stderr, proc.stderr
    assert "Traceback" in proc.stderr, proc.stderr


def test_the_metadata_schema_needs_neither():
    """And this is what the move bought (issue #381).

    A program that writes the metadata document and converts nothing --
    the reason the schema is versioned and ontology-mapped at all -- gets
    ``MSIMetadata`` and ``validate_document`` without the imaging stack.
    Checking the two dependencies during ``import thyra``, which is how
    #310 was fixed, would have made this fail: a parent package always
    initialises before its subpackage.
    """
    proc = _run_child(_BLOCK_BOTH_FOR_THE_SCHEMA)

    assert proc.returncode == 0, f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
    assert proc.stdout.startswith("OK MSIMetadata validate_document")
