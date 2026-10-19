"""Tests for ``ImzMLReader._validate_parser_state`` and its two consumers.

pyimzml trusts its own offsets absolutely and so, until now, did Thyra: a
declared offset that points past the end of the ``.ibd`` yields an empty array
rather than an error, and the affected pixels leave the store without a word.
These tests pin the refusals, and -- just as importantly -- pin the shapes that
must *not* be refused.

Fault fixtures are built at runtime from ``ImzMLWriter`` plus byte surgery on a
copy. That is deliberately the cheap half: writer output cannot reach the
structural features of real vendor files, which is why the committed
hand-authored corpus exists separately.
"""

import logging
import os
import re
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, List, Optional
from xml.etree import ElementTree  # nosec B405

import numpy as np
import pytest
from pyimzml.compression import ZlibCompression
from pyimzml.ImzMLParser import ImzMLParser

from thyra.convert import _create_converter
from thyra.errors import ConversionRefused
from thyra.metadata.document import read_metadata_document
from thyra.preview import preview_msi
from thyra.readers.imzml import imzml_reader as imzml_reader_module
from thyra.readers.imzml.imzml_reader import ImzMLReader

_MODULE_LOGGER_NAME = "thyra.readers.imzml.imzml_reader"

# The three real files live outside the repository: ``test_data/`` is
# gitignored (``.gitignore:55``), so the "a real file is still accepted" tests
# below can never run in CI and skip cleanly when the directory is absent --
# which it also is in a git worktree. ``THYRA_TEST_DATA`` overrides the
# location.
_REAL_DATA_DIR = Path(
    os.environ.get("THYRA_TEST_DATA", Path(__file__).resolve().parents[3] / "test_data")
)


@contextmanager
def _capture_module_logs(level: int = logging.WARNING) -> Iterator[List[str]]:
    """Collect the reader module's log records at ``level`` and above.

    A handler on the named logger rather than pytest's ``caplog``, which needs
    propagation to the root handler and can be silently switched off by
    whatever ran before this test.
    """
    records: List[str] = []

    class _Capture(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(record.getMessage())

    module_logger = logging.getLogger(_MODULE_LOGGER_NAME)
    handler = _Capture(level=level)
    previous_level = module_logger.level
    module_logger.addHandler(handler)
    module_logger.setLevel(level)
    try:
        yield records
    finally:
        module_logger.removeHandler(handler)
        module_logger.setLevel(previous_level)


def write_imzml(
    directory: Path,
    name: str = "sample",
    n_spectra: int = 6,
    n_peaks: int = 5,
    **writer_kwargs,
) -> Path:
    """Write a small well-formed processed-mode imzML/ibd pair.

    Args:
        directory: Where to write.
        name: File stem.
        n_spectra: How many spectra to write.
        n_peaks: Peaks per spectrum.
        **writer_kwargs: Passed straight to ``ImzMLWriter``.

    Returns:
        Path to the written ``.imzML``.
    """
    # Imported here so the module-level import list stays about what is being
    # tested rather than about how fixtures are made.
    from pyimzml.ImzMLWriter import ImzMLWriter

    path = directory / f"{name}.imzML"
    mzs = np.linspace(100.0, 500.0, n_peaks)
    intensities = np.arange(1.0, n_peaks + 1.0)

    with ImzMLWriter(str(path), mode="processed", **writer_kwargs) as writer:
        for i in range(n_spectra):
            writer.addSpectrum(mzs, intensities, (i % 3 + 1, i // 3 + 1, 1))
    return path


def _nth_binary_array_block(text: str, spectrum: int, array: str) -> "re.Match[str]":
    """Locate one ``<binaryDataArray>`` block in an imzML document.

    Every block is matched first and the wanted array selected afterwards. A
    single pattern anchored on the ``referenceableParamGroupRef`` would match
    from the *previous* array's opening tag, because the lazy ``.*?`` in front
    of the ref happily spans a block boundary.

    Args:
        text: The whole imzML document.
        spectrum: 0-based spectrum index.
        array: ``"mzArray"`` or ``"intensityArray"``.

    Returns:
        The match covering that spectrum's block for that array.
    """
    blocks = [
        match
        for match in re.finditer(
            r"<binaryDataArray\b.*?</binaryDataArray>", text, re.DOTALL
        )
        if f'ref="{array}"' in match.group(0)
    ]
    assert len(blocks) > spectrum, f"only {len(blocks)} {array} blocks found"
    return blocks[spectrum]


def poison_cv_param(
    imzml_path: Path, spectrum: int, array: str, accession: str, value: str
) -> None:
    """Rewrite one cvParam value inside one spectrum's binary array block.

    Args:
        imzml_path: The file to edit in place.
        spectrum: 0-based spectrum index.
        array: ``"mzArray"`` or ``"intensityArray"``.
        accession: The cvParam accession to rewrite.
        value: Its replacement value.
    """
    text = imzml_path.read_text(encoding="utf-8")
    block = _nth_binary_array_block(text, spectrum, array)
    edited = re.sub(
        rf'(accession="{accession}"[^>]*?value=")[^"]*(")',
        rf"\g<1>{value}\g<2>",
        block.group(0),
        count=1,
    )
    assert edited != block.group(0), f"{accession} not found in {array}[{spectrum}]"
    imzml_path.write_text(
        text[: block.start()] + edited + text[block.end() :], encoding="utf-8"
    )


def spectrum_end_byte(imzml_path: Path, spectrum: int) -> int:
    """Return the last byte the given spectrum occupies in the ``.ibd``.

    Args:
        imzml_path: A well-formed imzML whose ``.ibd`` is intact.
        spectrum: 0-based spectrum index.

    Returns:
        The exclusive end byte of that spectrum's intensity array.
    """
    parser = ImzMLParser(str(imzml_path), parse_lib="ElementTree")
    try:
        return int(
            parser.intensityOffsets[spectrum]
            + parser.intensityLengths[spectrum]
            * parser.sizeDict[parser.intensityPrecision]
        )
    finally:
        parser.m.close()


def set_spectrum_list_count(imzml_path: Path, value: int) -> Path:
    """Rewrite the ``<spectrumList>`` count attribute and nothing else.

    Args:
        imzml_path: The file to edit in place.
        value: The count to declare.

    Returns:
        The same path, for chaining onto ``write_imzml``.
    """
    text = imzml_path.read_text(encoding="utf-8")
    edited, n = re.subn(
        r'(<spectrumList\b[^>]*?count=")[^"]*(")',
        rf"\g<1>{value}\g<2>",
        text,
    )
    assert n == 1, f"expected one spectrumList count, rewrote {n}"
    imzml_path.write_text(edited, encoding="utf-8")
    return imzml_path


def empty_the_spectrum_list(imzml_path: Path, declared: Optional[int] = None) -> Path:
    """Drop every ``<spectrum>``, leaving the head and the opening tag intact.

    ``ImzMLWriter`` cannot write this directly -- it sizes the raster with
    ``max()`` over the coordinates and raises on an empty one -- so the two
    halves are made separately here.

    The two counts a real file can carry are genuinely different faults, and
    ``declared`` is which one to build: left alone, the attribute keeps saying
    what the writer wrote, which is a truncated export -- a header sized for an
    acquisition that then produced nothing. ``0`` is a file that is honestly
    empty and says so.

    Args:
        imzml_path: The file to edit in place.
        declared: The count to restate, or None to keep the written one.

    Returns:
        The same path, for chaining onto ``write_imzml``.
    """
    text = imzml_path.read_text(encoding="utf-8")
    edited, n = re.subn(
        r"(<spectrumList\b[^>]*>).*?</spectrumList>",
        r"\g<1></spectrumList>",
        text,
        flags=re.DOTALL,
    )
    assert n == 1, f"expected one <spectrumList>, rewrote {n}"
    imzml_path.write_text(edited, encoding="utf-8")
    if declared is not None:
        return set_spectrum_list_count(imzml_path, declared)
    return imzml_path


def _edit_text(path: Path, old: str, new: str, count: int = 1) -> None:
    """Replace ``old`` by ``new`` in the document, which must contain it."""
    text = path.read_text(encoding="utf-8")
    assert old in text, f"{old!r} not in {path.name}"
    path.write_text(text.replace(old, new, count), encoding="utf-8")


def _refusal(path: Path) -> ConversionRefused:
    """The refusal ``get_essential_metadata`` raises for this file."""
    with pytest.raises(ConversionRefused) as excinfo:
        reader = ImzMLReader(path)
        try:
            reader.get_essential_metadata()
        finally:
            reader.close()
    return excinfo.value


class TestCleanFilesAreAccepted:
    """The validator must be invisible on anything well-formed."""

    def test_writer_output_passes_with_no_warnings(self, temp_dir):
        path = write_imzml(temp_dir)
        with _capture_module_logs() as records:
            reader = ImzMLReader(path)
            reader._ensure_parser_initialized()
            reader.close()
        assert records == []

    def test_32_bit_integer_intensity_is_accepted(self, temp_dir):
        """MS:1000519 is spec-legal and converts correctly today.

        pyimzml's own writer emits ``32-bit integer`` for
        ``intensity_dtype=np.int32``. A precision allow-list of ``("f", "d")``
        -- which is what the audit proposed -- refuses this file.
        """
        path = write_imzml(temp_dir, intensity_dtype=np.int32)
        with _capture_module_logs() as records:
            reader = ImzMLReader(path)
            reader._ensure_parser_initialized()
            assert reader.parser.intensityPrecision == "i"
            reader.close()
        assert records == []

    def test_multiple_spectra_read_back_after_validation(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=6)
        reader = ImzMLReader(path)
        try:
            assert reader.n_spectra == 6
            mzs, intensities = reader.parser.getspectrum(3)
            assert mzs.size == 5 and intensities.size == 5
        finally:
            reader.close()


class TestOffsetsAgainstTheIbd:
    """Audit #4 -- nothing validated offsets or lengths against the .ibd."""

    def test_truncated_ibd_names_the_spectrum_and_both_byte_figures(self, temp_dir):
        """An interrupted copy cut on a spectrum boundary.

        pyimzml returns zero-length arrays for every spectrum past the cut
        without raising, so before this check ``convert_msi`` returned True
        and wrote a store containing only the pixels before it.
        """
        path = write_imzml(temp_dir, n_spectra=6)
        cut = spectrum_end_byte(path, 2)
        ibd = path.with_suffix(".ibd")
        full_size = ibd.stat().st_size
        with open(ibd, "r+b") as handle:
            handle.truncate(cut)
        assert cut < full_size

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as excinfo:
            reader._ensure_parser_initialized()

        message = str(excinfo.value)
        assert "spectrum 3" in message
        assert f"{cut:,}" in message
        assert ibd.name in message

    def test_offset_past_the_end_of_the_ibd_is_refused(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=6)
        ibd_size = path.with_suffix(".ibd").stat().st_size
        poison_cv_param(path, 4, "mzArray", "IMS:1000102", str(ibd_size + 4096))

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused, match=r"spectrum 4 declares a m/z array"):
            reader._ensure_parser_initialized()

    def test_leading_negative_offset_is_refused(self, temp_dir):
        """pyimzml's offset repair is a no-op on a *leading* negative.

        ``__fix_offsets`` seeds ``prev_value`` with NaN and ``nan >= 0`` is
        False, so a negative on spectrum 0 survives the repair untouched and
        every later seek is computed from it.
        """
        path = write_imzml(temp_dir, n_spectra=6)
        poison_cv_param(path, 0, "mzArray", "IMS:1000102", "-16")

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused, match=r"negative m/z offset"):
            reader._ensure_parser_initialized()

    def test_array_length_disagreement_is_refused(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=6, n_peaks=5)
        poison_cv_param(path, 2, "intensityArray", "IMS:1000103", "3")

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as excinfo:
            reader._ensure_parser_initialized()
        message = str(excinfo.value)
        assert "spectrum 2" in message
        assert "5" in message and "3" in message

    def test_non_monotonic_offsets_warn_rather_than_refuse(self, temp_dir):
        """Legal, but it is the precondition pyimzml's repair assumes away."""
        path = write_imzml(temp_dir, n_spectra=6)
        parser = ImzMLParser(str(path), parse_lib="ElementTree")
        first, second = int(parser.mzOffsets[1]), int(parser.mzOffsets[2])
        parser.m.close()
        # Every spectrum here is the same length, so swapping two offsets
        # leaves both reads inside the file.
        poison_cv_param(path, 1, "mzArray", "IMS:1000102", str(second))
        poison_cv_param(path, 2, "mzArray", "IMS:1000102", str(first))

        with _capture_module_logs() as records:
            reader = ImzMLReader(path)
            reader._ensure_parser_initialized()
            reader.close()
        assert any("not monotonically non-decreasing" in r for r in records)

    def test_unaccounted_trailing_bytes_warn_rather_than_refuse(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=6)
        ibd = path.with_suffix(".ibd")
        with open(ibd, "ab") as handle:
            handle.write(b"\x00" * 64)

        with _capture_module_logs() as records:
            reader = ImzMLReader(path)
            reader._ensure_parser_initialized()
            reader.close()
        assert any("unaccounted for" in r for r in records)

    def test_a_later_spectrum_stating_more_values_than_it_stores_is_refused(
        self, temp_dir
    ):
        """The F188 layout: six spectra of three peaks, spectrum 1 says four.

        Its last m/z value was the first bytes of its intensity array, and
        its last intensity the first bytes of the next m/z array: an invented
        m/z of 134,217,760 joined the mass range. The extent check only
        looks at the last array in the file, and spectrum 0 was the only one
        whose lengths were compared.
        """
        path = write_imzml(
            temp_dir,
            n_spectra=6,
            n_peaks=3,
            mz_dtype=np.float64,
            intensity_dtype=np.float32,
        )
        poison_cv_param(path, 1, "mzArray", "IMS:1000103", "4")
        poison_cv_param(path, 1, "intensityArray", "IMS:1000103", "4")

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as excinfo:
            reader._ensure_parser_initialized()

        message = str(excinfo.value)
        assert "spectrum 1" in message
        assert "runs into" in message

    def test_a_continuous_file_sharing_one_mz_block_builds_quietly(self, temp_dir):
        """Identical byte ranges are not an overlap.

        Every spectrum of a continuous file points at the same m/z block, so
        a check that refused any two arrays on the same bytes would refuse
        every continuous file.
        """
        from pyimzml.ImzMLWriter import ImzMLWriter

        path = temp_dir / "continuous.imzML"
        mzs = np.linspace(100.0, 500.0, 5)
        with ImzMLWriter(str(path), mode="continuous") as writer:
            for i in range(6):
                writer.addSpectrum(
                    mzs, np.arange(1.0, 6.0) + i, (i % 3 + 1, i // 3 + 1, 1)
                )

        with _capture_module_logs() as records:
            reader = ImzMLReader(path)
            reader._ensure_parser_initialized()
            reader.close()

        assert records == []


class TestBinaryArrayDeclarations:
    """Audit #6, #7, #8 and #13 -- what the param groups claim."""

    def test_zlib_compression_is_refused_naming_the_accession(self, temp_dir):
        """pyimzml 1.5.5 has no decompression path at all.

        Ungated, it reads ``IMS:1000103 x itemsize`` raw deflate bytes and
        ``np.frombuffer`` turns them into numbers of exactly the declared
        length, so every one of Thyra's emptiness and length gates passes.
        """
        path = write_imzml(
            temp_dir,
            mz_compression=ZlibCompression(),
            intensity_compression=ZlibCompression(),
        )
        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused, match=r"MS:1000574"):
            reader._ensure_parser_initialized()

    def test_two_precision_terms_in_one_group_are_refused(self, temp_dir):
        """Schema-valid, and pyimzml picks between them by dictionary order.

        The ranking is ``32-bit float < 64-bit float < 32-bit integer <
        64-bit integer`` regardless of document order, so a genuinely float32
        array declared alongside ``64-bit float`` decodes as float64 garbage
        at the declared length.
        """
        path = write_imzml(temp_dir)
        text = path.read_text(encoding="utf-8")
        needle = '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"'
        assert needle in text, "writer output changed shape; update the needle"
        extra = (
            '<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" '
            'value=""/>\n      '
        )
        path.write_text(text.replace(needle, extra + needle, 1), encoding="utf-8")

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused, match=r"declares 2 precision terms"):
            reader._ensure_parser_initialized()

    def test_64_bit_integer_is_refused_as_platform_dependent(self, temp_dir):
        """``SIZE_DICT['l']`` is 8 everywhere; ``np.dtype('l')`` is not.

        pyimzml reads N*8 bytes and decodes them as 2N int32 values on
        Windows and N int64 values on Linux. One file, two answers -- which
        is the durable defect here, not corruption.
        """
        path = write_imzml(temp_dir, intensity_dtype=np.int64)
        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused, match=r"64-bit integer"):
            reader._ensure_parser_initialized()

    def test_encoded_length_contradicting_the_precision_is_refused(self, temp_dir):
        """Audit #8 -- a correct accession carrying a wrong name.

        ``ACCESSION_FIX_MAPPING`` rewrites the accession and keeps the raw
        name, and pyimzml derives the precision from the name, so the file
        decodes at the wrong width at exactly the right length. Only
        ``IMS:1000104`` disagrees, and it is checked here on spectrum 0.
        """
        path = write_imzml(temp_dir, n_peaks=5)
        # 5 float64 values are 40 bytes; claim they were encoded as 20.
        poison_cv_param(path, 0, "mzArray", "IMS:1000104", "20")

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as excinfo:
            reader._ensure_parser_initialized()
        message = str(excinfo.value)
        assert "IMS:1000104" in message
        assert "20" in message and "40" in message

    def test_two_scan_settings_blocks_warn_rather_than_refuse(self, temp_dir):
        """Audit #9's neighbour: a per-accession chimera, not an error yet.

        Refusing belongs with the pixel-size unit work; Thyra's grid model
        cannot represent per-region pixel sizes either way.
        """
        path = write_imzml(temp_dir)
        text = path.read_text(encoding="utf-8")
        match = re.search(r"<scanSettings\b.*?</scanSettings>", text, re.DOTALL)
        assert match is not None
        duplicate = match.group(0).replace("scanSettings1", "scanSettings2", 1)
        path.write_text(
            text[: match.end()] + "\n" + duplicate + text[match.end() :],
            encoding="utf-8",
        )

        with _capture_module_logs() as records:
            reader = ImzMLReader(path)
            reader._ensure_parser_initialized()
            reader.close()
        assert any("<scanSettings> blocks" in r for r in records)

    @pytest.mark.parametrize("keeps_encoded_length", [False, True])
    def test_an_array_stating_its_own_precision_against_its_group_is_refused(
        self, temp_dir, keeps_encoded_length
    ):
        """F177 -- pyimzml decodes at the group's width and ignores the array.

        Spectrum 0's intensity array states 32-bit float on itself while its
        group says 64-bit. It used to decode at 64 bits: 514 stored for a row
        that sums to 15. With the encoded length (``IMS:1000104``) left in,
        it agrees with the group and the old encoded-length check passed;
        with it removed there was nothing to compare at all.
        """
        path = write_imzml(temp_dir, mz_dtype=np.float64, intensity_dtype=np.float64)
        text = path.read_text(encoding="utf-8")
        block = _nth_binary_array_block(text, 0, "intensityArray")
        own = block.group(0).replace(
            '<referenceableParamGroupRef ref="intensityArray"/>',
            '<referenceableParamGroupRef ref="intensityArray"/>'
            '<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" '
            'value=""/>',
            1,
        )
        if not keeps_encoded_length:
            own = re.sub(r"<cvParam[^>]*IMS:1000104[^>]*/>", "", own)
            assert "IMS:1000104" not in own
        path.write_text(
            text[: block.start()] + own + text[block.end() :], encoding="utf-8"
        )

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as excinfo:
            reader._ensure_parser_initialized()

        message = str(excinfo.value)
        assert "spectrum 0" in message
        assert "intensity array" in message
        assert "32-bit float" in message
        assert "64-bit float" in message
        assert "'intensityArray'" in message


class TestAnUnreadableDocumentIsOneRefusal:
    """F54, F48 and F192 -- pyimzml's own errors reached the user as tracebacks.

    A cut-short or empty document, a position written ``2.0``, a position in
    a parameter group or a second ``<scan>``, and a second group declaring
    the m/z array each ended in four ERROR lines and a raw ``ParseError``,
    ``ValueError``, ``TypeError`` or ``AttributeError``. Each is now one
    ``ConversionRefused`` naming the file and the cause.

    They use ``temp_dir``, whose removal fails on Windows while a file in it
    is open: a refusal raised inside the parse once left the document's
    handle to the garbage collector.
    """

    def test_a_document_cut_in_half_says_it_is_not_well_formed(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=4)
        text = path.read_bytes()
        path.write_bytes(text[: len(text) // 2])

        refusal = _refusal(path)

        assert "not well-formed XML" in str(refusal)
        assert path.name in str(refusal)

    def test_the_original_error_is_kept_as_the_cause(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=4)
        text = path.read_bytes()
        path.write_bytes(text[: len(text) // 2])

        refusal = _refusal(path)

        assert isinstance(refusal.__cause__, ElementTree.ParseError)

    def test_an_empty_file_says_it_is_not_well_formed(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=4)
        path.write_bytes(b"")

        refusal = _refusal(path)

        assert "not well-formed XML" in str(refusal)

    def test_a_position_written_with_a_decimal_point_names_the_spectrum(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=4)
        text = path.read_text(encoding="utf-8")
        edited, n = re.subn(
            r'(name="position y" value="2)"', r'\g<1>.0"', text, count=1
        )
        assert n == 1
        path.write_text(edited, encoding="utf-8")

        message = str(_refusal(path))

        assert "position y (IMS:1000051)" in message
        assert "'2.0'" in message
        assert "Positions are whole numbers" in message
        assert "spectrum" in message

    def test_a_position_held_in_a_parameter_group_names_spectrum_and_group(
        self, temp_dir
    ):
        """The last spectrum's x and y moved into a group its scan references."""
        path = write_imzml(temp_dir, n_spectra=4)
        text = path.read_text(encoding="utf-8")
        positions = list(
            re.finditer(
                r'<cvParam[^>]*"IMS:1000050"[^>]*/>\s*<cvParam[^>]*"IMS:1000051"'
                r"[^>]*/>",
                text,
            )
        )
        last = positions[-1]
        text = (
            text[: last.start()]
            + '<referenceableParamGroupRef ref="pos"/>'
            + text[last.end() :]
        )
        text = text.replace(
            "</referenceableParamGroupList>",
            f'<referenceableParamGroup id="pos">{last.group(0)}'
            "</referenceableParamGroup></referenceableParamGroupList>",
            1,
        )
        path.write_text(text, encoding="utf-8")

        message = str(_refusal(path))

        assert "spectrum 3" in message
        assert "position x (IMS:1000050)" in message
        assert "parameter group 'pos'" in message

    def test_a_position_in_a_second_scan_behind_an_empty_first_scan(self, temp_dir):
        """The F48 repro: an empty ``<scan/>`` ahead of the one with the position."""
        path = write_imzml(temp_dir, n_spectra=4)
        _edit_text(path, '<scanList count="1">', '<scanList count="2"><scan/>')

        message = str(_refusal(path))

        assert "spectrum 0" in message
        assert "<scan> 2" in message
        assert "reads only the first <scan>" in message

    @pytest.mark.parametrize("used", [False, True], ids=["unused", "used"])
    def test_a_second_mz_group_names_both_groups(self, temp_dir, used):
        """The F192 repro: a second group carrying the m/z array term.

        Declared and used by no spectrum, or used by every other one, the
        parser took one group and raised an ``AttributeError`` on a spectrum
        that referenced the other.
        """
        path = write_imzml(temp_dir, n_spectra=6, mz_dtype=np.float64)
        text = path.read_text(encoding="utf-8")
        group = re.search(
            r'<referenceableParamGroup id="mzArray">.*?</referenceableParamGroup>',
            text,
            re.DOTALL,
        )
        assert group is not None
        text = text.replace(
            group.group(0),
            group.group(0) + group.group(0).replace('id="mzArray"', 'id="mzArrayB"'),
            1,
        )
        if used:
            parts = text.split("<spectrum ")
            parts = [
                (
                    p.replace('ref="mzArray"', 'ref="mzArrayB"', 1)
                    if k % 2 == 0 and k
                    else p
                )
                for k, p in enumerate(parts)
            ]
            text = "<spectrum ".join(parts)
        path.write_text(text, encoding="utf-8")

        message = str(_refusal(path))

        assert "The file declares the m/z array in 2 groups" in message
        assert "'mzArray', 'mzArrayB'" in message


class TestRefusalCleansUp:
    """Refusing must not leave the reader holding the .ibd."""

    def test_the_ibd_handle_is_released_when_validation_refuses(self, temp_dir):
        """Windows will not let the caller delete a file Thyra still holds."""
        path = write_imzml(temp_dir)
        poison_cv_param(path, 0, "mzArray", "IMS:1000102", "-16")

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused):
            reader._ensure_parser_initialized()
        assert reader.ibd_file is None
        path.with_suffix(".ibd").unlink()


class TestARefusedFileIsNotParsedTwice:
    """Initialization is expensive, and every public entry point triggers it.

    ``n_spectra``, ``get_essential_metadata`` and four more all call
    ``_ensure_parser_initialized``, so without a memo a refused file is
    re-parsed once per attempt -- 63 s of XML each time on a 2.1 GB imzML --
    only to fail identically.
    """

    def test_the_second_attempt_re_raises_without_re_parsing(
        self, temp_dir, monkeypatch
    ):
        path = write_imzml(temp_dir, n_spectra=6)
        poison_cv_param(path, 0, "mzArray", "IMS:1000102", "-16")

        constructions: List[str] = []
        real_parser = imzml_reader_module.ImzMLParser

        def _counting_parser(*args, **kwargs):
            constructions.append(kwargs.get("filename", ""))
            return real_parser(*args, **kwargs)

        monkeypatch.setattr(imzml_reader_module, "ImzMLParser", _counting_parser)

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as first:
            reader._ensure_parser_initialized()
        with pytest.raises(ConversionRefused) as second:
            reader._ensure_parser_initialized()

        assert len(constructions) == 1
        assert second.value is first.value

    def test_the_memo_keeps_the_original_type_and_message(self, temp_dir, monkeypatch):
        """Rebuilding the error as ``type(e)(str(e))`` would lose both.

        Not every failure here is a ValueError Thyra raised: pyimzml itself
        raises whatever its XML and its own dependencies raise, and those
        constructors do not all take a single message. The reader names the
        file and the cause in one refusal (D63) and keeps the original as
        its cause; the memo hands back that same refusal.
        """

        class _TwoArgError(Exception):
            def __init__(self, code: int, detail: str) -> None:
                super().__init__(f"{code}: {detail}")
                self.code = code

        def _refusing_parser(*args, **kwargs):
            raise _TwoArgError(7, "the vendor library said no")

        monkeypatch.setattr(imzml_reader_module, "ImzMLParser", _refusing_parser)

        reader = ImzMLReader(write_imzml(temp_dir))
        refusals = []
        for _ in range(2):
            with pytest.raises(ConversionRefused) as excinfo:
                reader._ensure_parser_initialized()
            cause = excinfo.value.__cause__
            assert isinstance(cause, _TwoArgError)
            assert cause.code == 7
            assert str(cause) == "7: the vendor library said no"
            assert "7: the vendor library said no" in str(excinfo.value)
            refusals.append(excinfo.value)
        assert refusals[0] is refusals[1]

    def test_a_clean_file_is_still_parsed_exactly_once(self, temp_dir, monkeypatch):
        path = write_imzml(temp_dir, n_spectra=6)

        constructions: List[str] = []
        real_parser = imzml_reader_module.ImzMLParser

        def _counting_parser(*args, **kwargs):
            constructions.append(kwargs.get("filename", ""))
            return real_parser(*args, **kwargs)

        monkeypatch.setattr(imzml_reader_module, "ImzMLParser", _counting_parser)

        reader = ImzMLReader(path)
        try:
            assert reader.n_spectra == 6
            reader._ensure_parser_initialized()
        finally:
            reader.close()
        assert len(constructions) == 1


class TestConverterCreationDoesNotSwallowRefusals:
    """The first metadata read must fail loudly, not at DEBUG.

    ``_should_use_streaming`` used to make that read and once hid the
    validator's failure behind a size estimate; the size gate is gone,
    and ``_create_converter`` now makes the read itself, outside any
    try, for the same reason: the converter's constructor extracts the
    same metadata inside a try that logs at DEBUG.
    """

    def test_metadata_failure_propagates(self, temp_dir):
        class _RefusingReader:
            def get_essential_metadata(self):
                raise ValueError("imzML spectrum 3 declares a m/z array ending at ...")

        from thyra.core.base_converter import PixelSizeSource

        # Deliberately NOT narrowed to ConversionRefused, unlike every other
        # ValueError assertion in this directory. The stub above raises a plain
        # ValueError from the test body, and the claim is that
        # _create_converter swallows a metadata failure of ANY type -- narrowing
        # this would assert something the reader never raises and test less.
        with pytest.raises(ValueError, match="spectrum 3"):
            _create_converter(
                "spatialdata",
                _RefusingReader(),
                temp_dir / "out.zarr",
                "ds",
                10.0,
                PixelSizeSource.USER_PROVIDED,
                False,
                {},
            )


class TestPreviewReportsARefusedFile:
    """``preview_msi`` drives the Ousia Import Wizard's per-sample card.

    The refusal survives the metadata-only path of issue #360, but not the
    message that used to come with it. A preview no longer parses the
    document, so it no longer holds every spectrum's offsets: it reads the
    last array the file declares out of the tail of the XML and compares
    that one end byte against the ``.ibd``. So a truncated binary is still
    refused, and still before anyone converts it, but the spectrum named is
    the last one rather than the first one affected -- and the byte named is
    what the file claims, not where it was cut.

    The check is one-sided on purpose and the second test below is what
    says so: it can miss a bad file whose offsets are not in document
    order, and it can never refuse a good one.
    """

    def test_refused_file_is_unreadable_with_the_validator_message(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=6)
        cut = spectrum_end_byte(path, 2)
        intact_end = path.with_suffix(".ibd").stat().st_size
        with open(path.with_suffix(".ibd"), "r+b") as handle:
            handle.truncate(cut)

        preview = preview_msi(path)
        assert preview.readable is False
        assert preview.error is not None
        # The last spectrum, whose array is the one read out of the tail.
        assert "spectrum 6" in preview.error
        assert f"{intact_end:,}" in preview.error
        assert f"{cut:,}" in preview.error
        # The card still renders; it just renders the error.
        assert preview.n_pixels == 0
        assert preview.grid_dims == (0, 0)

    def test_a_file_whose_binary_is_whole_is_not_refused(self, temp_dir):
        """Guard the guard: the tail check fires on a shortfall, not on a tail.

        Without this, a check that refused every file would pass the test
        above for the wrong reason -- and it would refuse every acquisition
        in the wizard.
        """
        path = write_imzml(temp_dir, n_spectra=6)

        preview = preview_msi(path)

        assert preview.readable is True, preview.error
        assert preview.error is None

    def test_clean_file_is_still_readable(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=6)
        preview = preview_msi(path)
        assert preview.readable is True
        assert preview.n_pixels == 6

    def test_a_too_long_encoded_length_on_the_last_array_is_not_a_truncation(
        self, temp_dir
    ):
        """F155 -- the last array states 40 bytes, 20 are stored, all values present.

        The preview ended the last array at offset plus its encoded length
        (``IMS:1000104``) and called the file truncated, and ``thyra
        metadata`` exited 1; the conversion read it. The end is now the
        offset plus the value count times the item size, so all three agree.
        """
        path = write_imzml(
            temp_dir,
            n_spectra=4,
            n_peaks=5,
            mz_dtype=np.float64,
            intensity_dtype=np.float32,
        )
        poison_cv_param(path, 3, "intensityArray", "IMS:1000104", "40")

        preview = preview_msi(path)
        document = read_metadata_document(path)

        assert preview.readable is True, preview.error
        assert document["ms_analysis"]["n_spectra"] == 4


class TestIbdUuid:
    """``IMS:1000080`` against the first 16 bytes of the ``.ibd``.

    The pair is the only thing that tells an ``.imzML`` apart from a
    *different* acquisition's ``.ibd`` renamed to sit beside it: every other
    check reads the XML's own offsets and lengths against the binary's size,
    which a wrong-but-similar file satisfies. It warns rather than refuses --
    see ``ImzMLReader._check_ibd_uuid`` and design decision D17.
    """

    @staticmethod
    def _set_declared_uuid(imzml_path: Path, value: str) -> None:
        """Rewrite the ``IMS:1000080`` value in the file description."""
        text = imzml_path.read_text(encoding="utf-8")
        edited, n = re.subn(
            r'(accession="IMS:1000080"[^>]*?value=")[^"]*(")',
            rf"\g<1>{value}\g<2>",
            text,
        )
        assert n == 1, f"expected one IMS:1000080 cvParam, rewrote {n}"
        imzml_path.write_text(edited, encoding="utf-8")

    @staticmethod
    def _read_warnings(path: Path) -> List[str]:
        with _capture_module_logs() as records:
            reader = ImzMLReader(path)
            reader._ensure_parser_initialized()
            reader.close()
        return [r for r in records if "binary-file UUID" in r]

    def test_matching_uuid_is_silent(self, temp_dir):
        """pyimzml's own writer puts the same UUID in both places."""
        path = write_imzml(temp_dir)
        assert self._read_warnings(path) == []

    def test_mismatched_uuid_warns_and_still_reads(self, temp_dir):
        """The bellini case: two different values, a readable file."""
        path = write_imzml(temp_dir)
        self._set_declared_uuid(path, "{FC37F303-A9C0-4CD3-A28E-1D18E523C269}")

        with _capture_module_logs() as records:
            reader = ImzMLReader(path)
            reader._ensure_parser_initialized()
            n = reader.n_spectra
            reader.close()

        assert n == 6, "the warning must not cost the file its spectra"
        said = [r for r in records if "binary-file UUID" in r]
        assert len(said) == 1
        assert "fc37f303a9c04cd3a28e1d18e523c269" in said[0]

    def test_braces_and_case_are_not_a_mismatch(self, temp_dir):
        """IONTOF writes the registry braces, SCiLS does not; neither is data.

        The hyphens are positional too. Only the 32 hex digits are compared,
        so re-spelling the file's own UUID in the other convention must stay
        silent rather than reporting a disagreement with itself.
        """
        path = write_imzml(temp_dir)
        header = path.with_suffix(".ibd").read_bytes()[:16].hex()
        plain = uuid.UUID(hex=header)
        for spelling in (
            str(plain),
            str(plain).upper(),
            "{" + str(plain).upper() + "}",
            "  {" + str(plain) + "}  ",
        ):
            self._set_declared_uuid(path, spelling)
            assert self._read_warnings(path) == [], f"warned on {spelling!r}"

    def test_a_file_declaring_no_uuid_is_silent(self, temp_dir):
        """Nothing to compare is not a disagreement."""
        path = write_imzml(temp_dir)
        text = path.read_text(encoding="utf-8")
        edited, n = re.subn(r"\s*<cvParam[^>]*IMS:1000080[^>]*/>", "", text)
        assert n == 1
        path.write_text(edited, encoding="utf-8")
        assert self._read_warnings(path) == []


class TestADocumentThatDeclaresNoSpectra:
    """A ``<spectrumList>`` holding nothing, refused before the parser exists.

    Every other refusal in this file comes out of ``_validate_parser_state``,
    which runs against a parser that was built. This one cannot: pyimzml's
    constructor ends by taking the z extent off the coordinates it collected,
    ``np.asarray(self.coordinates)[:, 2].max()``, and an empty list makes that
    a 1-D array -- so the constructor dies with ``IndexError: too many indices
    for array`` and leaves nothing to validate. A preview reported exactly that
    sentence as its whole account of the file.

    ``ImzMLMetadataExtractor._extract_essential_impl`` refuses the same file
    with "No coordinates found in ImzML file", so the fact was already known
    one stage later; it was only ever unreachable.
    """

    def test_an_empty_spectrum_list_is_refused_by_name(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=6)
        empty_the_spectrum_list(path, declared=0)

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as excinfo:
            reader.get_essential_metadata()

        message = str(excinfo.value)
        assert "sample.imzML" in message
        assert "declares no spectra" in message
        assert "<spectrumList>" in message
        assert "too many indices" not in message

    def test_a_count_the_list_does_not_back_up_is_named(self, temp_dir):
        """count="6" over nothing is a different fault from count="0"."""
        path = write_imzml(temp_dir, n_spectra=6)
        empty_the_spectrum_list(path)

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as excinfo:
            reader.get_essential_metadata()

        assert "count attribute says 6" in str(excinfo.value)

    def test_a_count_of_zero_over_nothing_is_not_reported_twice(self, temp_dir):
        """The count agrees with the list, so there is one fault to report."""
        path = write_imzml(temp_dir, n_spectra=6)
        empty_the_spectrum_list(path, declared=0)

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as excinfo:
            reader.get_essential_metadata()

        assert "count attribute" not in str(excinfo.value)

    def test_a_document_with_no_spectrum_list_at_all_is_refused_too(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=6)
        text = path.read_text(encoding="utf-8")
        edited, n = re.subn(
            r"<spectrumList\b.*?</spectrumList>", "", text, flags=re.DOTALL
        )
        assert n == 1
        path.write_text(edited, encoding="utf-8")

        reader = ImzMLReader(path)
        with pytest.raises(ConversionRefused) as excinfo:
            reader.get_essential_metadata()

        assert "declares no spectra" in str(excinfo.value)

    def test_the_parser_is_never_built(self, temp_dir, monkeypatch):
        """A guard on the document, not a rewording of what pyimzml raised.

        Pinned because the cheaper fix was to catch the ``IndexError`` and
        rename it, which would go quiet the day pyimzml handles an empty
        coordinate list itself.
        """
        path = write_imzml(temp_dir, n_spectra=6)
        empty_the_spectrum_list(path)

        constructions: List[str] = []

        def _forbidden_parser(*args, **kwargs):
            constructions.append(kwargs.get("filename", ""))
            raise AssertionError("the parser must not be built for this file")

        monkeypatch.setattr(imzml_reader_module, "ImzMLParser", _forbidden_parser)

        with pytest.raises(ConversionRefused):
            ImzMLReader(path)._ensure_parser_initialized()

        assert constructions == []

    def test_the_preview_card_says_what_is_wrong_with_the_file(self, temp_dir):
        """Where this message is actually read: the Ousia Import Wizard."""
        path = write_imzml(temp_dir, n_spectra=6)
        empty_the_spectrum_list(path)

        preview = preview_msi(path)

        assert preview.readable is False
        assert preview.error is not None
        assert "declares no spectra" in preview.error
        assert "too many indices" not in preview.error
        assert preview.n_pixels == 0

    def test_a_file_that_holds_spectra_is_not_refused(self, temp_dir):
        """Guard the guard: a check that refused everything would pass above."""
        path = write_imzml(temp_dir, n_spectra=6)

        reader = ImzMLReader(path)
        try:
            assert reader.n_spectra == 6
        finally:
            reader.close()

    def test_a_count_of_zero_over_a_list_that_holds_spectra_is_not_refused(
        self, temp_dir
    ):
        """The count is read for the message, never for the verdict."""
        path = write_imzml(temp_dir, n_spectra=6)
        set_spectrum_list_count(path, 0)

        reader = ImzMLReader(path)
        try:
            assert reader.n_spectra == 6
        finally:
            reader.close()

    def test_a_document_that_will_not_parse_is_refused_in_one_line(self, temp_dir):
        """This check must never be the one that condemns a broken document.

        It reads the file before pyimzml does, so a file neither can parse
        would otherwise be reported as declaring no spectra -- which is a
        guess, and the wrong one. The reader says it is not well-formed XML,
        in one line, and keeps the parser's error as the cause.
        """
        path = write_imzml(temp_dir, n_spectra=6)
        path.write_text("<mzML><spectrumList count=", encoding="utf-8")

        with pytest.raises(ConversionRefused) as excinfo:
            ImzMLReader(path)._ensure_parser_initialized()

        message = str(excinfo.value)
        assert "declares no spectra" not in message
        assert "not well-formed XML" in message
        assert "\n" not in message


@pytest.mark.skipif(
    not (_REAL_DATA_DIR / "bellini.imzML").exists(),
    reason=f"real MSI corpus not present at {_REAL_DATA_DIR} (test_data/ is gitignored)",
)
class TestRealFilesAreStillAccepted:
    """A validator that refuses real data is worse than no validator.

    ``test_data/`` is gitignored, so this can never run in CI and is skipped
    outright in a worktree. It is the whole point of the lane all the same:
    every check above has to survive contact with the files Thyra actually
    converts. Set ``THYRA_TEST_DATA`` to point it somewhere else.

    ``20240826_xenium_0041899.imzML`` is deliberately absent: its XML alone
    takes ~63 s to parse, which belongs in the integration lane rather than
    here.
    """

    #: What each real file is allowed to say, and nothing else. ``bellini``
    #: is an IONTOF SurfaceLab export whose ``IMS:1000080`` and ``.ibd``
    #: header hold two different UUIDs -- measured, not hypothetical: the XML
    #: says ``{FC37F303-...C269}`` and the binary begins
    #: ``3ad1bacd...f731``, with the first spectrum at byte 16 so the header
    #: slot is genuinely populated. Every other check passes and the file
    #: converts correctly, which is exactly why that check warns instead of
    #: refusing (issue #261, design decision D17). Listing the warning here
    #: rather than dropping the assertion keeps the rest of the guarantee:
    #: any *other* message from these files still fails the test.
    _ALLOWED_WARNINGS = {
        "bellini.imzML": ("imzML declares binary-file UUID",),
        "pea.imzML": (),
    }

    @pytest.mark.parametrize("name", ["bellini.imzML", "pea.imzML"])
    def test_real_file_initialises_with_no_errors_and_no_warnings(self, name):
        path = _REAL_DATA_DIR / name
        if not path.exists():
            pytest.skip(f"{name} not present")

        with _capture_module_logs() as records:
            reader = ImzMLReader(path)
            reader._ensure_parser_initialized()
            n = reader.n_spectra
            reader.close()

        allowed = self._ALLOWED_WARNINGS[name]
        unexpected = [r for r in records if not any(r.startswith(a) for a in allowed)]

        assert n > 0
        assert unexpected == [], f"{name} tripped the validator: {unexpected}"
        for prefix in allowed:
            assert any(
                r.startswith(prefix) for r in records
            ), f"{name} no longer warns {prefix!r}; the check or the file changed"
