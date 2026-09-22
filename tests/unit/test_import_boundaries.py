"""What each front door is allowed to import, measured in a fresh process.

Issue #381: ``thyra/__init__.py`` imported the reader and converter
packages for their registration side effects, and a parent package always
initialises first, so ``import thyra.metadata.schema`` -- a pydantic model
and a controlled vocabulary -- loaded every vendor reader, the SpatialData
converter, ``spatialdata``, ``dask`` and ``anndata``: 3.5 s and 3,580
modules before the caller's first line.

Every assertion here needs its own interpreter. By the time a test session
reaches this file it has imported most of the package, so ``sys.modules``
in this process says nothing about what any single import pulls in.

These are boundaries, not budgets. They name modules that must be absent,
not a module count or a time, because a count moves with every
dependency's own imports and would fail for reasons that have nothing to
do with the boundary. What made the old cost invisible was that no test
could see the difference between "the schema imports spatialdata" and
"something else in the process did".
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

import pytest

import thyra

#: The checkout, so the child imports the tree under test rather than
#: whatever ``site-packages`` holds.
_REPO_ROOT = Path(thyra.__file__).resolve().parent.parent

_CENSUS = """
import json
import sys

{statement}

modules = sorted(sys.modules)
json.dump(
    {{
        "thyra": [m for m in modules if m == "thyra" or m.startswith("thyra.")],
        "roots": sorted({{m.split(".")[0] for m in modules}}),
    }},
    sys.stdout,
)
"""


def _census(statement: str) -> Dict[str, List[str]]:
    """Run ``statement`` in a fresh interpreter and report its ``sys.modules``."""
    proc = subprocess.run(
        [sys.executable, "-c", _CENSUS.format(statement=statement)],
        capture_output=True,
        text=True,
        cwd=_REPO_ROOT,
    )
    assert proc.returncode == 0, (
        f"child failed running {statement!r}\n"
        f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
    )
    return json.loads(proc.stdout)


def _under(prefix: str, modules: List[str]) -> List[str]:
    return [m for m in modules if m == prefix or m.startswith(f"{prefix}.")]


#: The smallest complete document: no reader, no store, no pixels. The
#: assertion is there so the census cannot pass on a document that failed
#: to fill the provenance block the version goes into.
_BUILD_A_DOCUMENT = """
import thyra.metadata.schema as schema

document = schema.build_msi_metadata(
    None, pixel_size_um=(10.0, 10.0), source_format="imzml"
)
assert document.provenance.thyra_version
"""


class TestTheMetadataSchemaStandsAlone:
    """The acceptance criterion of issue #381.

    A second implementation of the metadata document -- the schema is a
    candidate to be cited by an external guideline, which implies one --
    reaches ``MSIMetadata`` and ``validate_document`` and nothing else. It
    has no business loading an imaging converter to do it, and for
    instrument families this package does not read, it could not.
    """

    @pytest.fixture(scope="class")
    def census(self) -> Dict[str, List[str]]:
        return _census("import thyra.metadata.schema")

    def test_it_loads_no_converter(self, census):
        assert _under("thyra.converters", census["thyra"]) == []

    def test_it_loads_no_reader(self, census):
        assert _under("thyra.readers", census["thyra"]) == []

    def test_it_loads_no_resampling(self, census):
        """The extractors are the path this used to arrive by.

        ``thyra/metadata/__init__.py`` re-exported the six vendor
        extractors eagerly, and three of them took ``SpectrumType`` or
        ``ImzMLAccessions`` from ``thyra.resampling.constants``, so
        importing the schema imported all 25 resampling modules. Both ends
        are fixed: the re-export is lazy, and those constants moved to
        ``thyra.metadata.constants``, which is where a description of the
        data belongs.
        """
        assert _under("thyra.resampling", census["thyra"]) == []

    @pytest.mark.parametrize("package", ["spatialdata", "anndata", "dask"])
    def test_it_loads_no_output_stack(self, census, package):
        assert package not in census["roots"]

    def test_the_schema_itself_is_there(self, census):
        assert "thyra.metadata.schema.models" in census["thyra"]


class TestBuildingADocumentStaysInsideTheBoundary:
    """The boundary has to hold at call time, not only at import time.

    ``build_msi_metadata`` reads the package version with ``from thyra
    import __version__`` while it fills the provenance block, so building
    a document re-enters the package root. That was a real trap while the
    root imported the readers and the converter: importing the schema
    could be made cheap and the first document would pay the cost anyway.

    It costs nothing now, because the root it re-enters is empty and
    already imported -- but that is a property of ``thyra/__init__.py``
    that nothing else states, and moving one import back there would
    restore the trap silently. This asserts the outcome rather than the
    mechanism, so a future ``_version.py`` or ``importlib.metadata``
    reading would satisfy it unchanged.
    """

    @pytest.fixture(scope="class")
    def census(self) -> Dict[str, List[str]]:
        return _census(_BUILD_A_DOCUMENT)

    def test_no_converter_is_imported(self, census):
        assert _under("thyra.converters", census["thyra"]) == []

    @pytest.mark.parametrize("package", ["spatialdata", "anndata", "dask"])
    def test_no_output_stack_is_imported(self, census, package):
        assert package not in census["roots"]


class TestTheFrontDoorRegistersNothing:
    """``import thyra`` costs what the caller asked for: nothing.

    Registration moved into the registry, which imports a format's module
    the first time that format is looked up. Anything that used to come in
    through the front door for its side effect now arrives when it is
    needed.
    """

    @pytest.fixture(scope="class")
    def census(self) -> Dict[str, List[str]]:
        return _census("import thyra")

    @pytest.mark.parametrize(
        "module",
        ["thyra.convert", "thyra.preview", "thyra.converters", "thyra.readers"],
    )
    def test_it_imports_no_heavy_submodule(self, census, module):
        assert module not in census["thyra"]

    @pytest.mark.parametrize("package", ["spatialdata", "anndata", "dask", "pandas"])
    def test_it_imports_no_heavy_dependency(self, census, package):
        assert package not in census["roots"]


class TestTheCommandLineDoesNotPayBeforeItParses:
    """``thyra --version`` and ``thyra validate`` reach no converter.

    ``thyra/__main__.py`` bound ``convert_msi`` at module scope, so every
    invocation -- ``--version``, ``--help``, ``validate`` on a standalone
    JSON document, ``metadata`` -- imported the converter and ``pandas``
    before click looked at the first argument. The conversion command
    imports it where it calls it.
    """

    @pytest.fixture(scope="class")
    def census(self) -> Dict[str, List[str]]:
        return _census("import thyra.__main__")

    def test_the_converter_is_not_imported(self, census):
        assert _under("thyra.converters", census["thyra"]) == []

    @pytest.mark.parametrize("package", ["spatialdata", "anndata"])
    def test_the_output_stack_is_not_imported(self, census, package):
        assert package not in census["roots"]

    def test_the_conversion_entry_point_is_not_imported(self, census):
        assert "thyra.convert" not in census["thyra"]


class TestTheLazyNamesStillResolve:
    """Deferring an import must not remove a name.

    ``thyra.convert_msi``, ``thyra.preview_msi``, ``thyra.MsiPreview``,
    ``thyra.read_metadata_document`` and ``thyra.SpatialDataConverter``
    were module attributes bound by the eager imports; they are :pep:`562`
    lookups now. ``thyra.readers`` and ``thyra.converters`` were bound as
    a side effect of importing them for registration, and resolve the same
    way.
    """

    @pytest.mark.parametrize(
        "name",
        [
            "MsiPreview",
            "SpatialDataConverter",
            "convert_msi",
            "preview_msi",
            "read_metadata_document",
        ],
    )
    def test_the_public_api_resolves(self, name):
        assert getattr(thyra, name) is not None

    @pytest.mark.parametrize("name", ["converters", "readers"])
    def test_the_subpackages_resolve(self, name):
        assert getattr(thyra, name).__name__ == f"thyra.{name}"

    def test_dir_lists_what_getattr_serves(self):
        """``dir(thyra)`` is the discoverable surface, and it must not shrink.

        A module with ``__getattr__`` and no ``__dir__`` lists only what is
        already in its ``__dict__``, so tab completion on a fresh import
        would have shown ``__version__`` and nothing else.
        """
        listed = dir(thyra)
        for name in thyra.__all__:
            assert name in listed
        assert "readers" in listed and "converters" in listed

    def test_an_unknown_name_is_still_an_attribute_error(self):
        with pytest.raises(AttributeError, match="no attribute 'not_a_name'"):
            thyra.not_a_name
