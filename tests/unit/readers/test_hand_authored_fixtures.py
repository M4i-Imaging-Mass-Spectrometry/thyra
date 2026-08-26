"""Tests over the hand-authored imzML corpus in ``tests/data/fixtures``.

Every other imzML the suite reads was produced by pyimzml's own
``ImzMLWriter``, so parser and writer agree on each other's mistakes and nine
structural features of real vendor files are unreachable. These four pairs were
written out literally to break that loop; see the corpus README for what each
one carries and for the two ways they can be destroyed without a test noticing.

A number of assertions below pin behaviour that is **wrong**. Each carries a
comment naming the audit finding it characterises and what a fix would turn it
into. ``assert imzmldict["pixel size x"] == 4406.25`` is a record of pyimzml
discarding a nanometre unit, not a claim that 4406.25 um is the pixel size --
the extractor compensates, which is what ``TestUnitNanometre`` now asserts.
"""

import shutil
import subprocess
import warnings
from pathlib import Path
from typing import List

import numpy as np
import pytest
from pyimzml.ImzMLParser import ImzMLParser

from tests.fixtures.imzml_parser import production_parser
from thyra.convert import convert_msi
from thyra.readers.imzml.imzml_reader import ImzMLReader

REPO_ROOT = Path(__file__).resolve().parents[3]
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "fixtures"

STEMS = [
    "iontof_sparse",
    "unit_nanometre",
    "two_scansettings",
    "two_precision_terms",
]

# two_precision_terms declares its precision twice, which _validate_parser_state
# refuses -- so it cannot be opened through the production reader. The corpus-wide
# tests below use this list; the refusal itself is asserted in
# TestTwoPrecisionTerms.
READABLE_STEMS = [stem for stem in STEMS if stem != "two_precision_terms"]

# Measured on the committed bytes. iontof_sparse.imzML is the only fixture
# terminated with CRLF, and every one of its line breaks is one.
IONTOF_CRLF_COUNT = 174
IONTOF_BYTE_COUNT = 11593


def fixture_path(stem: str) -> Path:
    """Return the imzML of a corpus pair, failing loudly if it is missing."""
    path = FIXTURE_DIR / f"{stem}.imzML"
    assert path.exists(), f"{path} is missing -- did .gitignore eat the corpus?"
    return path


def open_reader(stem: str) -> ImzMLReader:
    """Open a corpus fixture through the production reader."""
    reader = ImzMLReader(fixture_path(stem))
    reader._ensure_parser_initialized()
    return reader


def raw_bytes(stem: str) -> bytes:
    """Read a fixture's imzML with no text-mode translation whatsoever."""
    return fixture_path(stem).read_bytes()


class TestCorpusParsesThroughTheReader:
    """Every pair must survive the path production actually uses."""

    @pytest.mark.parametrize("stem", READABLE_STEMS)
    def test_reader_opens_the_file(self, stem):
        """The fixture parses and reports coherent essential metadata."""
        reader = open_reader(stem)
        try:
            essential = reader.get_essential_metadata()

            assert essential.n_spectra > 0
            assert essential.n_spectra == len(reader.parser.coordinates)
            assert len(essential.dimensions) == 3
            assert all(d > 0 for d in essential.dimensions)
            assert essential.mass_range[0] < essential.mass_range[1]
        finally:
            reader.close()

    @pytest.mark.parametrize("stem", READABLE_STEMS)
    def test_spectra_round_trip_out_of_the_ibd(self, stem):
        """The hand-packed .ibd decodes to the values it was packed from.

        The offsets, the IMS:1000103 lengths and the byte layout were all
        written by hand, so this is the check that the corpus is internally
        consistent rather than merely parseable.
        """
        reader = open_reader(stem)
        try:
            for index in range(len(reader.parser.coordinates)):
                mzs, intensities = reader.parser.getspectrum(index)

                assert mzs.size == intensities.size
                assert mzs.size > 0
                assert np.all(np.diff(mzs) > 0), "m/z must ascend"
                assert np.all(intensities > 0)
        finally:
            reader.close()

    @pytest.mark.parametrize("stem", STEMS)
    def test_index_attribute_is_zero_based(self, stem):
        """All three real files number spectra from 0; ImzMLWriter numbers from 1.

        Nothing in Thyra reads ``spectrum/@index`` today, so this is a property
        of the corpus rather than a bug guard -- but it is one of the nine
        divergences that made a generated fixture unable to stand in for a real
        file, and it is free to keep true.
        """
        text = raw_bytes(stem).decode("iso-8859-1")
        indices = sorted(
            int(chunk.split('"')[0]) for chunk in text.split('index="')[1:]
        )

        assert indices == list(range(len(indices)))
        assert indices[0] == 0


class TestIontofSparse:
    """The bellini class: everything an IONTOF export does that pyimzml does not."""

    def test_declares_iso_8859_1(self):
        """Real IONTOF exports are Latin-1; ImzMLWriter only ever emits UTF-8."""
        assert raw_bytes("iontof_sparse").startswith(
            b'<?xml version="1.0" encoding="ISO-8859-1"?>'
        )

    def test_every_line_break_is_crlf(self):
        """CRLF is the trait, and it is the one most easily lost.

        Checked here against the worktree file. ``TestCommittedBytes`` makes the
        same assertion against the git blob, which is the copy a normalising
        ``.gitattributes`` would silently rewrite.
        """
        data = raw_bytes("iontof_sparse")

        assert len(data) == IONTOF_BYTE_COUNT
        assert data.count(b"\r\n") == IONTOF_CRLF_COUNT
        assert data.count(b"\n") == data.count(b"\r\n"), "a bare LF crept in"

    def test_is_unindented(self):
        """No leading whitespace on any line.

        Unindented CRLF is precisely the layout that makes lxml raise
        ``XMLSyntaxError: xmlSAX2Characters`` on a perfectly good file, which is
        why the reader pins ``parse_lib="ElementTree"``. An indented fixture
        does not reproduce it.
        """
        lines = raw_bytes("iontof_sparse").split(b"\r\n")

        assert not [line for line in lines if line[:1] in (b" ", b"\t")]

    def test_has_no_mzml_version_attribute(self):
        """mzML 1.1.0 requires @version. IONTOF omits it; pyimzml never looks."""
        text = raw_bytes("iontof_sparse").decode("iso-8859-1")
        opening_tag = text[text.index("<mzML") : text.index(">", text.index("<mzML"))]

        assert "version=" not in opening_tag

    def test_has_no_position_z_term(self):
        """IMS:1000052 appears zero times in all three real files.

        Every ImzMLWriter fixture has it, so the branch where pyimzml
        synthesises z is exercised only by real data -- and by this file.
        """
        assert b"IMS:1000052" not in raw_bytes("iontof_sparse")

        reader = open_reader("iontof_sparse")
        try:
            assert {c[2] for c in reader.parser.coordinates} == {1}
        finally:
            reader.close()

    def test_misspelled_names_are_corrected_with_a_warning(self):
        """Three ontology names are wrong, and pyimzml repairs each with a UserWarning.

        The accessions are correct throughout; only the names are off, the way
        real bellini has them. The repair is silent apart from the warning, and
        ``thyra/convert.py`` suppresses the IMS:1000046 one by module name --
        so this is the only place the name-correction path is observable.
        """
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            reader = open_reader("iontof_sparse")
            reader.close()

        corrections = [
            str(w.message)
            for w in caught
            if "found with incorrect name" in str(w.message)
        ]
        corrected_accessions = {
            accession
            for accession in ("IMS:1000042", "IMS:1000043", "IMS:1000046")
            if any(accession in message for message in corrections)
        }

        assert corrected_accessions == {"IMS:1000042", "IMS:1000043", "IMS:1000046"}

    def test_max_dimension_contradicts_pixel_size(self):
        """The two spatial cvParams disagree and neither carries a unit.

        4 pixels x 4.40625 um is 17.625 um of raster; ``max dimension`` says 15.
        Deriving um/pixel from max dimension gives 3.75, reading ``pixel size``
        gives 4.40625. Real bellini is 12.8% apart in exactly this way, while
        pea and xenium are self-consistent -- so code that reads the wrong
        cvParam is invisible on SCiLS files and silently wrong on IONTOF ones.
        """
        reader = open_reader("iontof_sparse")
        try:
            meta = reader.parser.imzmldict
            derived = meta["max dimension x"] / meta["max count of pixels x"]
            declared = meta["pixel size x"]

            assert derived == pytest.approx(3.75)
            assert declared == pytest.approx(4.40625)
            assert abs(derived - declared) / declared > 0.1
        finally:
            reader.close()

    def test_acquisition_is_sparse(self):
        """Six of sixteen declared grid positions were acquired.

        A sparse acquisition is what makes the write routes' disagreement about
        empty pixel rows observable at all -- three of four drop them, the PCS
        route keeps them as all-zero observations.
        """
        reader = open_reader("iontof_sparse")
        try:
            meta = reader.parser.imzmldict
            declared_positions = (
                meta["max count of pixels x"] * meta["max count of pixels y"]
            )

            assert len(reader.parser.coordinates) == 6
            assert declared_positions == 16
        finally:
            reader.close()


class TestUnitNanometre:
    """Audit #9, closed: the discarded ``unitAccession``, worth a factor of 1000.

    ``_extract_pixel_size_fast`` now pairs the bare ``imzmldict`` number with
    the unit preserved on the ParamGroup path and converts to micrometres, so
    the assertions below state the truth rather than characterise the error.
    """

    def _variant_with_unit(self, tmp_path, unit_accession: str, unit_name: str):
        """Copy the fixture pair into ``tmp_path`` with the unit swapped."""
        text = raw_bytes("unit_nanometre").decode("utf-8")
        text = text.replace(
            'unitAccession="UO:0000018"', f'unitAccession="{unit_accession}"'
        )
        text = text.replace('unitName="nanometer"', f'unitName="{unit_name}"')

        imzml = tmp_path / "unit_variant.imzML"
        imzml.write_text(text, encoding="utf-8")
        shutil.copyfile(
            fixture_path("unit_nanometre").with_suffix(".ibd"),
            imzml.with_suffix(".ibd"),
        )
        return imzml

    def test_the_document_declares_nanometre(self):
        """The unit is in the file. Whether anything reads it is the next test."""
        text = raw_bytes("unit_nanometre").decode("utf-8")

        assert 'accession="IMS:1000046"' in text
        assert 'unitAccession="UO:0000018"' in text
        assert 'unitName="nanometer"' in text

    def test_imzmldict_drops_the_unit(self):
        """pyimzml's dict is still unit-blind; the extractor must not trust it.

        ``convert_cv_param`` takes no unit argument, so ``__readimzmlmeta``
        cannot carry one and ``imzmldict['pixel size x']`` is a bare number in
        whatever unit the vendor declared. This pins the reason
        ``_extract_pixel_size_fast`` reads the unit off the ParamGroup path
        rather than taking this number at face value.
        """
        reader = open_reader("unit_nanometre")
        try:
            assert reader.parser.imzmldict["pixel size x"] == pytest.approx(4406.25)
            assert reader.parser.imzmldict["pixel size y"] == pytest.approx(4406.25)
        finally:
            reader.close()

    def test_nanometre_pixel_size_reaches_essential_metadata_in_um(self):
        """4406.25 nanometre reads as 4.40625 um, not 4406.25 um."""
        reader = open_reader("unit_nanometre")
        try:
            assert reader.get_essential_metadata().pixel_size == pytest.approx(
                (4.40625, 4.40625)
            )
        finally:
            reader.close()

    def test_converted_coordinates_reach_disk_in_um(self, temp_dir):
        """End to end, because the point of #9 was that ``convert_msi`` returned
        ``True`` and wrote coordinates 1000x too large rather than failing. For
        a 2x2 acquisition at 4.40625 um the x positions are ``[0, 4.40625]``.
        """
        import spatialdata

        output = temp_dir / "unit_nanometre.zarr"

        assert convert_msi(fixture_path("unit_nanometre"), output, streaming=False)

        table = next(iter(spatialdata.read_zarr(output).tables.values()))

        assert sorted(set(table.obs["spatial_x"])) == pytest.approx([0.0, 4.40625])

    def test_millimetre_pixel_size_is_scaled_up(self, tmp_path):
        """UO:0000016 goes the other way: x1000 rather than /1000."""
        imzml = self._variant_with_unit(tmp_path, "UO:0000016", "millimeter")
        reader = ImzMLReader(imzml)
        try:
            assert reader.get_essential_metadata().pixel_size == pytest.approx(
                (4406250.0, 4406250.0)
            )
        finally:
            reader.close()

    def test_an_unrecognised_unit_is_refused(self, tmp_path):
        """A unit outside mm/um/nm raises instead of guessing a factor.

        UO:0000015 is centimetre -- a real length unit deliberately not in the
        table, so accepting it unconverted would be a silent x10000 error.
        """
        imzml = self._variant_with_unit(tmp_path, "UO:0000015", "centimeter")
        reader = ImzMLReader(imzml)
        try:
            with pytest.raises(ValueError, match="UO:0000015"):
                reader.get_essential_metadata()
        finally:
            reader.close()

    def test_the_unit_survives_on_the_param_group_path(self):
        """``cv_params`` keeps what ``imzmldict`` threw away.

        This is the path the fix reads: the unit is in the parsed metadata,
        one attribute along a different access path from the bare number.
        """
        reader = open_reader("unit_nanometre")
        try:
            scan_settings = next(iter(reader.parser.metadata.scan_settings.values()))
            units = {
                param[1]: param[6]
                for param in scan_settings.cv_params
                if param[1] in ("IMS:1000046", "IMS:1000047")
            }

            assert units == {
                "IMS:1000046": "UO:0000018",
                "IMS:1000047": "UO:0000018",
            }
        finally:
            reader.close()


class TestTwoScanSettings:
    """Audit #9's neighbour: first-match-anywhere across two blocks."""

    def test_two_blocks_are_declared(self):
        """Both blocks survive into ``metadata.scan_settings``.

        All three real files declare ``count="1"``, so nothing in the corpus
        reaches the multi-block path. This fixture is the discriminator for
        lane I's check (h), ``len(scan_settings) == 1``.
        """
        reader = open_reader("two_scansettings")
        try:
            settings = reader.parser.metadata.scan_settings

            assert len(settings) == 2
            assert set(settings) == {"scanSettings1", "scanSettings2"}
        finally:
            reader.close()

    def test_imzmldict_is_a_chimera_of_both_blocks(self):
        """CHARACTERISATION of audit #9's neighbour -- lane I check (h) warns on this.

        ``__readimzmlmeta`` resolves each accession independently by
        first-match-anywhere below ``scanSettingsList``, so the dict mixes the
        blocks: pixel counts and pixel size from block 1, max dimensions from
        block 2, which is the only block declaring them. The resulting triple
        describes neither block, and block 2's own 25.0 um pixel size is never
        seen by anything.

        Lane I adds a warning here rather than a refusal; a full fix belongs
        with the #9 work, which no lane in handout I owns.
        """
        reader = open_reader("two_scansettings")
        try:
            meta = reader.parser.imzmldict

            assert meta["max count of pixels x"] == 2  # block 1 only
            assert meta["max dimension x"] == 200  # block 2 only
            assert meta["pixel size x"] == pytest.approx(10.0)  # block 1 wins

            block_two = reader.parser.metadata.scan_settings["scanSettings2"]
            block_two_pixel_size = {
                param[2] for param in block_two.cv_params if param[1] == "IMS:1000046"
            }

            assert block_two_pixel_size == {25.0}
        finally:
            reader.close()

    def test_the_reader_silently_picks_one(self):
        """CHARACTERISATION of audit #9's neighbour -- lane I check (h) warns on this.

        Thyra's grid model cannot represent per-region pixel sizes, so picking
        one of two contradictory declarations without saying so is the defect.
        Today the file converts as though block 2 did not exist.
        """
        reader = open_reader("two_scansettings")
        try:
            assert reader.get_essential_metadata().pixel_size == pytest.approx(
                (10.0, 10.0)
            )
        finally:
            reader.close()


class TestTwoPrecisionTerms:
    """Audit #7: precision resolved by dict order rather than by declaration."""

    def test_the_group_declares_both_precisions(self):
        """Two contradictory precision terms as siblings in one param group.

        Schema-valid, and the discriminator for lane I's check (f). That check
        must read ``param_by_name`` rather than the group's own ``cv_params``,
        because a legal file may declare its precision in a *referenced* group
        -- so this fixture puts both terms where either form would see them.
        """
        with production_parser(fixture_path("two_precision_terms")) as parser:
            group = parser.metadata.referenceable_param_groups[parser.mzGroupId]
            declared = {
                name
                for name in ("32-bit float", "64-bit float")
                if name in group.param_by_name
            }

            assert declared == {"32-bit float", "64-bit float"}

    def test_pyimzml_resolves_it_by_dict_insertion_order(self):
        """CHARACTERISATION of audit #7 -- lane I refuses this file.

        ``__process_metadata`` keeps the *last* match while scanning
        ``PRECISION_DICT`` in insertion order (32-bit float, 64-bit float,
        32-bit integer, 64-bit integer), so the winner is "64-bit float"
        whatever the document order says. Nothing warns.

        Lane I's ``_validate_parser_state()`` turns this into a raise; this
        assertion becomes ``pytest.raises`` in that diff.
        """
        with production_parser(fixture_path("two_precision_terms")) as parser:
            assert parser.mzPrecision == "d"

    def test_the_ambiguous_file_is_refused(self):
        """The ambiguous declaration is refused. Inverted from characterisation.

        The .ibd here really is float64, so pyimzml's guess happens to be
        right and the file reads correctly. That is deliberate: the fixture
        isolates the ambiguous *declaration*, so a validator that refuses it is
        refusing the ambiguity and not reacting to corrupt data. Get the
        declaration wrong on a genuinely float32 file and it decodes as float64
        garbage at the correct length, with ``convert_msi`` returning ``True``.
        """
        with pytest.raises(ValueError, match="declares 2 precision terms"):
            open_reader("two_precision_terms")

        # The bytes really are float64, so what is refused is the ambiguous
        # declaration and not unreadable data: read the same file through
        # pyimzml directly and the values are exactly what was packed.
        with production_parser(fixture_path("two_precision_terms")) as parser:
            mzs, _ = parser.getspectrum(0)

            assert mzs.dtype == np.float64
            assert mzs[0] == pytest.approx(100.0)
            assert mzs[-1] == pytest.approx(500.0)


def _git(*args: str) -> str:
    """Run git in the repository root and return stdout, or skip if it cannot."""
    if shutil.which("git") is None:
        pytest.skip("git is not on PATH")

    result = subprocess.run(
        ["git", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        pytest.skip(f"git {' '.join(args)} failed: {result.stderr.decode().strip()}")

    return result.stdout.decode("utf-8", errors="replace")


def _blob(repo_path: str) -> bytes:
    """Return the raw staged bytes of a tracked path.

    ``git cat-file -p`` reads the object exactly as stored -- no smudge filter,
    no eol conversion. Reading the worktree file instead would answer a
    different question, and it is the question that cannot fail.
    """
    if shutil.which("git") is None:
        pytest.skip("git is not on PATH")

    result = subprocess.run(
        ["git", "cat-file", "-p", f":{repo_path}"],
        cwd=REPO_ROOT,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        pytest.skip(f"{repo_path} is not in the git index")

    return result.stdout


def _is_ignored(repo_path: str) -> bool:
    """Report whether .gitignore would exclude a path.

    ``--no-index`` is load-bearing. Without it ``check-ignore`` declines to
    answer for anything already in the index, which is every file in the
    corpus -- so the question that matters, *would ``git add`` still take this
    file*, silently becomes unaskable and every answer comes back "no".
    """
    if shutil.which("git") is None:
        pytest.skip("git is not on PATH")

    result = subprocess.run(
        ["git", "check-ignore", "--no-index", "--", repo_path],
        cwd=REPO_ROOT,
        capture_output=True,
        check=False,
    )
    if result.returncode not in (0, 1):
        pytest.skip(f"git check-ignore failed: {result.stderr.decode().strip()}")

    return result.returncode == 0


class TestCommittedBytes:
    """The corpus is only useful if what is committed is what was authored.

    ``.gitattributes:15`` is ``* text=auto eol=lf``. Without an exemption,
    staging ``iontof_sparse.imzML`` rewrites all 174 CRLFs and the fixture stops
    testing anything -- with the suite green on both sides, because every other
    test in this module reads the worktree file. These read the blob.
    """

    def test_crlf_survives_into_the_index(self):
        """The committed iontof_sparse still has all 174 CRLFs and no bare LF."""
        blob = _blob("tests/data/fixtures/iontof_sparse.imzML")

        assert len(blob) == IONTOF_BYTE_COUNT
        assert blob.count(b"\r\n") == IONTOF_CRLF_COUNT
        assert blob.count(b"\n") == blob.count(b"\r\n")

    @pytest.mark.parametrize("stem", STEMS)
    @pytest.mark.parametrize("suffix", ["imzML", "ibd"])
    def test_blob_is_byte_identical_to_the_worktree(self, stem, suffix):
        """Nothing in the checkout/staging round trip may touch a fixture byte."""
        repo_path = f"tests/data/fixtures/{stem}.{suffix}"
        blob = _blob(repo_path)

        assert blob == (REPO_ROOT / repo_path).read_bytes()

    @pytest.mark.parametrize("stem", STEMS)
    def test_imzml_is_exempt_from_eol_normalisation(self, stem):
        """The exemption itself, so deleting it fails here and not somewhere subtle."""
        attributes = _git(
            "check-attr", "text", "--", f"tests/data/fixtures/{stem}.imzML"
        )

        assert attributes.strip().endswith("text: unset")

    @pytest.mark.parametrize("stem", STEMS)
    def test_ibd_is_marked_binary(self, stem):
        """.ibd files carry packed float64; a text heuristic must never see them."""
        attributes = _git(
            "check-attr", "binary", "--", f"tests/data/fixtures/{stem}.ibd"
        )

        assert attributes.strip().endswith("binary: set")

    def test_the_corpus_is_actually_tracked(self):
        """All eight files are in the index, not merely present in the checkout."""
        tracked = _git("ls-files", "tests/data/fixtures").split()
        expected: List[str] = [
            f"tests/data/fixtures/{stem}.{suffix}"
            for stem in STEMS
            for suffix in ("imzML", "ibd")
        ]

        assert set(expected) <= set(tracked)

    @pytest.mark.parametrize("stem", STEMS)
    @pytest.mark.parametrize("suffix", ["imzML", "ibd"])
    def test_gitignore_still_admits_the_corpus(self, stem, suffix):
        """``*imzML`` and ``*ibd`` are unanchored; the negations are what allow these.

        ``test_the_corpus_is_actually_tracked`` cannot see this. Once a file is
        in the index .gitignore no longer applies to it, so deleting the
        negations leaves that assertion green and breaks only the next person
        to regenerate a fixture and try to add it.
        """
        assert not _is_ignored(f"tests/data/fixtures/{stem}.{suffix}")

    def test_the_unanchored_ignore_rules_are_still_live(self):
        """The control: the same names one directory up are still ignored.

        Without it the test above would pass just as happily if somebody had
        deleted ``*imzML`` and ``*ibd`` from .gitignore outright, which is a
        different repository from the one it means to describe.
        """
        assert _is_ignored("tests/data/scratch.imzML")
        assert _is_ignored("tests/data/scratch.ibd")


class TestProductionResolvesTheIbdExplicitly:
    """The corpus is also the place to pin how the .ibd is found.

    Production passes ``ibd_file=`` explicitly; pyimzml's ``_infer_bin_filename``
    is the other rule, and the two disagree on case, on directories and on an
    unanchored regex. Thyra's ``with_suffix('.ibd').exists()`` is the more
    robust of the two -- see ``docs/imzml-parser-notes.md`` before "fixing" it
    to match.
    """

    def test_the_reader_hands_the_parser_an_open_handle(self):
        """``ImzMLReader`` never lets pyimzml infer the .ibd path."""
        reader = open_reader("iontof_sparse")
        try:
            assert reader.ibd_file is not None
            assert reader.parser.m is reader.ibd_file
            assert reader.ibd_path == fixture_path("iontof_sparse").with_suffix(".ibd")
        finally:
            reader.close()

    def test_an_explicit_handle_and_an_inferred_one_agree_here(self):
        """Both resolution rules find the same file for a well-named pair.

        They diverge only on the awkward cases (``SAMPLE.IBD``, a parent
        directory containing ``.ibd``, a half-finished ``sample.ibdtmp``). This
        pins the easy case so a change to either rule shows up as a diff.
        """
        path = fixture_path("iontof_sparse")
        inferred = ImzMLParser._infer_bin_filename(path)

        assert Path(inferred) == path.with_suffix(".ibd")
