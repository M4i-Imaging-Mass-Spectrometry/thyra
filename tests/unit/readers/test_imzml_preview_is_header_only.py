"""``preview_msi`` on an imzML parses no spectrum list and decodes no ``.ibd``.

``preview.py`` promises "No spectra are read" and passes
``metadata_only=True`` to say so. ``ImzMLReader.__init__`` swallowed the
kwarg through ``**kwargs``, so a preview built an ordinary parser: pyimzml
walked every ``<spectrum>`` element for its offsets, and
``_get_mass_range_processed`` then decoded each spectrum's m/z array out of
the binary to find the extrema. Measured on a 2.0 GiB imzML with 918,855
spectra: 64 s, against a documented ``<500 ms`` (issue #360). The same
file's head answers in 0.4 ms.

The assertions here are about which doors are opened, not about how long
it took: a wall-clock threshold on a synthetic file measures the machine.
``ImzMLParser`` is the one door into both the spectrum list and the
``.ibd`` -- ``_initialize_parser`` is its only construction site -- so
counting constructions settles it exactly.

What the head can and cannot answer is the other half, and each answer has
a test below: the raster and the spectrum count come off the declared
geometry, the mass range off the per-spectrum ``MS:1000528``/``MS:1000527``
that most writers record, and a file recording neither is reported as
unknown rather than guessed from its first spectrum.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from thyra.errors import ConversionRefused
from thyra.preview import preview_msi
from thyra.readers.imzml import imzml_reader as imzml_reader_module
from thyra.readers.imzml.imzml_reader import ImzMLReader

from .test_imzml_parser_state_validation import empty_the_spectrum_list, write_imzml

_OBSERVED_MZ_CVPARAM = re.compile(
    r'\s*<cvParam[^>]*accession="MS:100052[78]"[^>]*/>\s*'
)


@pytest.fixture
def parser_builds(monkeypatch):
    """Count every construction of the real pyimzml parser."""
    calls: list[str] = []
    original = imzml_reader_module.ImzMLParser

    def counting(*args, **kwargs):
        calls.append(str(kwargs.get("filename", args[0] if args else "")))
        return original(*args, **kwargs)

    monkeypatch.setattr(imzml_reader_module, "ImzMLParser", counting)
    return calls


def write_zero_based_imzml(directory: Path, n_spectra: int = 6) -> Path:
    """The same 3 x 2 raster as ``write_imzml``, written from origin 0.

    Kept here rather than folded into the shared helper: what it exists to
    produce is a file whose *declaration* is one short of its raster, which
    is a property of pyimzml's writer and belongs next to the test that
    reads it.
    """
    import numpy as np
    from pyimzml.ImzMLWriter import ImzMLWriter

    path = directory / "zero_based.imzML"
    mzs = np.linspace(100.0, 500.0, 5)
    intensities = np.arange(1.0, 6.0)

    with ImzMLWriter(str(path), mode="processed") as writer:
        for i in range(n_spectra):
            writer.addSpectrum(mzs, intensities, (i % 3, i // 3, 1))
    return path


def strip_observed_mz(path: Path) -> Path:
    """Remove every observed-m/z cvParam, as an IONTOF export has none.

    The binary and its offsets are untouched: what changes is only whether
    the document states the range it already contains.
    """
    text = path.read_text(encoding="utf-8")
    stripped, n = _OBSERVED_MZ_CVPARAM.subn("\n", text)
    assert n, "fixture carried no observed-m/z cvParams to strip"
    path.write_text(stripped, encoding="utf-8")
    return path


class TestThePreviewNeverBuildsAParser:
    def test_preview_parses_no_spectra(self, temp_dir, parser_builds):
        """The whole point: a preview reads the head, not the document."""
        path = write_imzml(temp_dir, n_spectra=6)

        preview = preview_msi(path)

        assert preview.readable, preview.error
        assert parser_builds == [], "preview built a pyimzml parser"

    def test_a_real_read_still_builds_one(self, temp_dir, parser_builds):
        """Guard the guard: the counter must be able to see a parse.

        Without this, a fixture that silently stopped intercepting
        ``ImzMLParser`` would make the test above pass for the wrong
        reason.
        """
        path = write_imzml(temp_dir, n_spectra=6)

        reader = ImzMLReader(path)
        try:
            reader.get_essential_metadata()
        finally:
            reader.close()

        assert len(parser_builds) == 1


class TestTheHeadAnswersTheCard:
    def test_the_declared_raster_and_count_are_reported(self, temp_dir):
        path = write_imzml(temp_dir, n_spectra=6)

        preview = preview_msi(path)

        assert preview.grid_dims == (3, 2)
        assert preview.n_pixels == 6

    def test_the_range_is_the_one_a_full_read_finds(self, temp_dir):
        """The recorded extrema are the extrema, not an approximation.

        This is what makes reading them off the XML a substitution rather
        than a trade: the writer computed them from the same arrays the
        ``.ibd`` scan would decode.
        """
        path = write_imzml(temp_dir, n_spectra=6)

        from_head = preview_msi(path).mz_range

        reader = ImzMLReader(path)
        try:
            from_binary = reader.get_essential_metadata().mass_range
        finally:
            reader.close()

        assert from_head == pytest.approx(from_binary)

    def test_a_file_recording_no_range_reports_unknown(self, temp_dir):
        """``None``, not the first spectrum's own extrema.

        IONTOF SurfaceLab writes no observed-m/z terms at all. Filling the
        gap from spectrum one was measured 28 Da narrow on a real file --
        a range narrow enough to be believed and wrong enough to mislead.
        """
        path = strip_observed_mz(write_imzml(temp_dir, n_spectra=6))

        preview = preview_msi(path)

        assert preview.readable, preview.error
        assert preview.mz_range is None
        # Everything the head does know is still reported.
        assert preview.grid_dims == (3, 2)
        assert preview.n_pixels == 6


class TestTheRangeScanIsChunked:
    def test_a_term_split_across_chunks_is_still_found(self, temp_dir, monkeypatch):
        """The scan reads 8 MB at a time, so terms straddle its boundaries.

        A cvParam cut in half by a read boundary matches in neither half,
        which is why each chunk carries the tail of the one before it. The
        overlap is shrunk here so a fixture small enough to commit crosses
        it many times; the file's real chunk size would need a 16 MB one.
        """
        from thyra.readers.imzml import header as header_module

        path = write_imzml(temp_dir, n_spectra=40)
        whole = header_module.scan_observed_mz_range(path)

        monkeypatch.setattr(header_module, "_SCAN_CHUNK_BYTES", 512)
        monkeypatch.setattr(header_module, "_SCAN_OVERLAP_BYTES", 256)
        chunked = header_module.scan_observed_mz_range(path)

        assert chunked == whole
        assert chunked is not None

    def test_the_overlap_does_not_double_a_bound(self, temp_dir, monkeypatch):
        """Matching an element twice is harmless, which is why it is allowed.

        Elements inside the carried-over tail are matched by two chunks in
        a row. A minimum and a maximum do not care, and this is the test
        that says so rather than a comment claiming it.
        """
        from thyra.readers.imzml import header as header_module

        path = write_imzml(temp_dir, n_spectra=40)
        monkeypatch.setattr(header_module, "_SCAN_CHUNK_BYTES", 300)
        monkeypatch.setattr(header_module, "_SCAN_OVERLAP_BYTES", 299)

        assert header_module.scan_observed_mz_range(path) == pytest.approx(
            (100.0, 500.0)
        )


def write_continuous_imzml(directory: Path, n_spectra: int = 40) -> Path:
    """The same raster as ``write_imzml``, in continuous mode.

    Every spectrum shares the one m/z array, which is what ``IMS:1000030``
    declares and what lets the range scan stop at the first spectrum.
    """
    import numpy as np
    from pyimzml.ImzMLWriter import ImzMLWriter

    path = directory / "continuous.imzML"
    mzs = np.linspace(100.0, 500.0, 5)
    intensities = np.arange(1.0, 6.0)

    with ImzMLWriter(str(path), mode="continuous") as writer:
        for i in range(n_spectra):
            writer.addSpectrum(mzs, intensities, (i % 3 + 1, i // 3 + 1, 1))
    return path


_HIGHEST_OBSERVED = re.compile(r'(accession="MS:1000527"[^>]*?value=")([^"]*)(")')


def raise_recorded_highest(path: Path, ordinal: int, value: float) -> Path:
    """Rewrite one spectrum's recorded highest observed m/z.

    ``ordinal`` counts occurrences in document order; ``-1`` is the last
    spectrum.  The binary is untouched: what changes is what the document
    CLAIMS, which is all the scan reads.
    """
    text = path.read_text(encoding="utf-8")
    matches = list(_HIGHEST_OBSERVED.finditer(text))
    assert matches, "fixture carried no highest-observed-m/z cvParams"
    match = matches[ordinal]
    text = f"{text[: match.start(2)]}{value}{text[match.end(2):]}"
    path.write_text(text, encoding="utf-8")
    return path


class TestAContinuousFileStopsAtItsFirstSpectrum:
    """A shared m/z axis is stated once, so the scan need not read it 117,911 times.

    Measured on a 270 MB continuous export over a network share: the full
    pass cost 5.6 s cold and 0.58 s warm for a range the first spectrum
    already stated (issue #371).  The early stop is exact by the
    specification -- continuous mode means one m/z array for every
    spectrum -- and it is checked against the file's own first chunk
    rather than taken on trust.
    """

    def test_the_preview_range_is_the_one_a_full_read_finds(self, temp_dir):
        path = write_continuous_imzml(temp_dir)

        from_head = preview_msi(path).mz_range

        reader = ImzMLReader(path)
        try:
            from_binary = reader.get_essential_metadata().mass_range
        finally:
            reader.close()

        assert from_head == pytest.approx(from_binary)

    def test_the_scan_stops_after_the_first_chunk(self, temp_dir, monkeypatch):
        """Observable only by planting a claim past the first chunk.

        The last spectrum is made to claim a higher m/z than the array it
        shares with every other spectrum can hold.  A full scan sees it; the
        shared-axis scan stops before it, and that is the documented trade:
        the specification says the value cannot differ, so the rest of the
        document is not read to check.  The plain scan is unchanged, which
        is what keeps a processed file honest.
        """
        from thyra.readers.imzml import header as header_module

        path = raise_recorded_highest(write_continuous_imzml(temp_dir), -1, 9999.0)
        monkeypatch.setattr(header_module, "_SCAN_CHUNK_BYTES", 1024)
        monkeypatch.setattr(header_module, "_SCAN_OVERLAP_BYTES", 256)

        assert header_module.scan_observed_mz_range(path, shared_axis=True) == (
            pytest.approx((100.0, 500.0))
        )
        assert header_module.scan_observed_mz_range(path) == pytest.approx(
            (100.0, 9999.0)
        )

    def test_a_disagreement_in_the_first_chunk_scans_in_full(self, temp_dir):
        """The specification is checked, not assumed.

        A writer that declares continuous and still records differing
        extrema on its first spectra is scanned to the end like any other
        file, and the range is the extremes of everything it recorded.
        """
        from thyra.readers.imzml import header as header_module

        path = raise_recorded_highest(write_continuous_imzml(temp_dir), 1, 9999.0)

        assert header_module.scan_observed_mz_range(path, shared_axis=True) == (
            pytest.approx((100.0, 9999.0))
        )

    def test_a_processed_file_is_still_scanned_in_full(self, temp_dir, monkeypatch):
        """The mode gates the early stop, and it is read off the head.

        The same planted claim on a processed file reaches the preview,
        because each of its spectra has its own array and the last one
        really can range higher than the first.
        """
        from thyra.readers.imzml import header as header_module

        path = raise_recorded_highest(write_imzml(temp_dir, n_spectra=40), -1, 9999.0)
        monkeypatch.setattr(header_module, "_SCAN_CHUNK_BYTES", 1024)
        monkeypatch.setattr(header_module, "_SCAN_OVERLAP_BYTES", 256)

        preview = preview_msi(path)

        assert preview.readable, preview.error
        assert preview.mz_range == pytest.approx((100.0, 9999.0))


class TestWhatTheHeadWillNotGuess:
    def test_a_zero_based_file_falls_back_to_the_coordinates(
        self, temp_dir, parser_builds
    ):
        """A 0-based file's declaration is ambiguous, so it is not used.

        ``max count of pixels x`` is a count to the specification and the
        largest coordinate to pyimzml's writer; for a 0-based file those
        differ by one, and nothing in the head says which was meant. The
        coordinate path settles it, so this file pays for a parse and gets
        the right raster rather than saving the parse and being one short.
        """
        path = write_zero_based_imzml(temp_dir)

        preview = preview_msi(path)

        assert preview.grid_dims == (3, 2)
        assert len(parser_builds) == 1, "the ambiguous file skipped the parse"

    def test_a_file_declaring_no_spectra_falls_back_too(self, temp_dir, parser_builds):
        """An empty acquisition keeps being refused, not described.

        The head could describe this file -- it declares a raster and a
        pitch -- and the card would read "a raster this size, zero pixels
        in it", which is a description of something nobody can convert. So
        the head declines and the coordinate path's refusal stands.

        The refusal is now the reader's own, raised off the document before
        the parser is built, which is why the fallback this file takes costs
        no parse at all. What says the head declined is the message: had it
        not, the card would read as a perfectly readable 3 x 2 raster with
        nothing in it.
        """
        path = empty_the_spectrum_list(write_imzml(temp_dir, n_spectra=6), declared=0)

        preview = preview_msi(path)

        assert preview.readable is False
        assert "declares no spectra" in (preview.error or "")
        assert parser_builds == [], "the fallback no longer needs a parse"

    def test_a_missing_ibd_is_still_refused(self, temp_dir):
        """A preview that called this readable would invite a refusal."""
        path = write_imzml(temp_dir, n_spectra=6)
        path.with_suffix(".ibd").unlink()

        preview = preview_msi(path)

        assert preview.readable is False
        assert ".ibd" in (preview.error or "")


class TestAMetadataOnlyReaderCanStillRead:
    def test_spectra_are_available_if_asked_for(self, temp_dir, parser_builds):
        """``metadata_only`` is where the metadata comes from, not a lock.

        Unlike the Bruker reader, which skips an SDK load it cannot undo,
        nothing here is given up at construction: the parser is simply not
        built for metadata. A caller that goes on to read spectra pays for
        it then, which is what ``preview_msi`` never does.
        """
        path = write_imzml(temp_dir, n_spectra=6)
        reader = ImzMLReader(path, metadata_only=True)
        try:
            reader.get_essential_metadata()
            assert parser_builds == []

            spectra = list(reader.iter_spectra())
        finally:
            reader.close()

        assert len(spectra) == 6
        assert len(parser_builds) == 1


def test_refusals_still_carry_their_type(temp_dir):
    """The missing-``.ibd`` refusal is a ``ConversionRefused``, as before.

    ``preview_msi`` turns every exception into ``readable=False``, so the
    type is only visible to a direct caller -- which the converter is.
    """
    path = write_imzml(temp_dir, n_spectra=6)
    path.with_suffix(".ibd").unlink()

    reader = ImzMLReader(path, metadata_only=True)
    with pytest.raises(ConversionRefused):
        reader.get_essential_metadata()
