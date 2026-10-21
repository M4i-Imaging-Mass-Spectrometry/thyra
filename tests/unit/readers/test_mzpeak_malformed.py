"""Malformed mzPeak archives are read or refused by name, never a traceback.

Each archive here is a good one from :mod:`tests.fixtures.mzpeak_builder`
with one thing changed: a column re-typed or dropped, a member cut short,
a field renamed. The reader either reads what is still there, or refuses
the archive with a :class:`ConversionRefused` that names the archive and
the member or column at fault.
"""

from __future__ import annotations

import logging
import zipfile

import numpy as np
import pyarrow as pa
import pytest

from tests.fixtures.mzpeak_builder import (
    DATA_MEMBER,
    INDEX_MEMBER,
    METADATA_MEMBER,
    PEAKS_MEMBER,
    POSITION_X_COLUMN,
    POSITION_Y_COLUMN,
    POSITION_Z_COLUMN,
    SCANS_MEMBER,
    Spectrum,
    build_mzpeak,
    grid_spectra,
    rewrite_member,
    rewrite_table,
)
from thyra.errors import ConversionRefused
from thyra.metadata.constants import SpectrumType
from thyra.readers.mzpeak import MzPeakReader, mzpeak_reader

READER_LOGGER = "thyra.readers.mzpeak.mzpeak_reader"


def _retyped(*columns: str, kind=pa.string()):
    """A table edit that casts some columns to another type."""

    def edit(table):
        for column in columns:
            table = table.set_column(
                table.column_names.index(column),
                column,
                table.column(column).cast(kind),
            )
        return table

    return edit


def _dropped(*columns: str):
    """A table edit that drops some columns."""
    return lambda table: table.drop_columns(list(columns))


def _coordinates(archive):
    """The coordinates of every spectrum the reader yields."""
    with MzPeakReader(archive) as reader:
        return [coords for coords, _, _ in reader.iter_spectra()]


def _refusal(archive) -> str:
    """The message the reader refuses the archive with."""
    with pytest.raises(ConversionRefused) as refused:
        with MzPeakReader(archive) as reader:
            reader.get_essential_metadata()
            list(reader.iter_spectra())
    return str(refused.value)


@pytest.fixture
def archive(tmp_path):
    """A 3x2 profile archive in the point layout."""
    return build_mzpeak(tmp_path / "a.mzpeak", grid_spectra(3, 2))


class TestPositionsStoredAsText:
    """A position column of whole numbers written as text (F42)."""

    def test_whole_numbers_as_text_are_positions(self, archive, thyra_logs):
        """The same pixels as the integer columns, and a warning."""
        expected = _coordinates(archive)
        rewrite_table(
            archive, SCANS_MEMBER, _retyped(POSITION_X_COLUMN, POSITION_Y_COLUMN)
        )

        with thyra_logs(READER_LOGGER, logging.WARNING) as records:
            assert _coordinates(archive) == expected

        assert "as text" in "\n".join(records.messages)

    def test_text_spellings_of_a_whole_number(self, archive):
        """'3.0' and ' 3 ' are 3; a blank is a scan on no pixel."""
        expected = _coordinates(archive)

        def edit(table):
            index = table.column_names.index(POSITION_X_COLUMN)
            values = [
                f" {v}.0 " for v in table.column(POSITION_X_COLUMN).to_pylist()
            ] + [""]
            ys = table.column(POSITION_Y_COLUMN).to_pylist() + [None]
            extra = {name: [None] for name in table.column_names}
            extra["source_index"] = [table.num_rows]
            grown = pa.concat_tables(
                [table, pa.table(extra, schema=table.schema)]
            ).combine_chunks()
            grown = grown.set_column(index, POSITION_X_COLUMN, pa.array(values))
            y_index = grown.column_names.index(POSITION_Y_COLUMN)
            return grown.set_column(
                y_index, POSITION_Y_COLUMN, pa.array(ys, type=pa.uint32())
            )

        rewrite_table(archive, SCANS_MEMBER, edit)

        assert _coordinates(archive) == expected

    def test_a_categorical_text_column_is_read(self, archive):
        """Dictionary-encoded text, as pandas writes a categorical."""
        expected = _coordinates(archive)

        def edit(table):
            index = table.column_names.index(POSITION_X_COLUMN)
            text = table.column(POSITION_X_COLUMN).cast(pa.string())
            return table.set_column(
                index, POSITION_X_COLUMN, text.combine_chunks().dictionary_encode()
            )

        rewrite_table(archive, SCANS_MEMBER, edit)

        assert _coordinates(archive) == expected

    def test_z_as_text_is_read(self, tmp_path):
        """A text z column states its plane as the integer one does."""
        archive = build_mzpeak(
            tmp_path / "z.mzpeak", grid_spectra(3, 2), position_z=[4] * 6
        )
        rewrite_table(archive, SCANS_MEMBER, _retyped(POSITION_Z_COLUMN))

        with MzPeakReader(archive) as reader:
            offsets = reader.get_essential_metadata().coordinate_offsets

        assert offsets is not None
        assert offsets[2] == 4

    def test_text_that_is_no_whole_number_is_refused(self, archive):
        """'1.5' is no pixel. The refusal names the column."""

        def edit(table):
            values = [f"{v}.5" for v in table.column(POSITION_X_COLUMN).to_pylist()]
            index = table.column_names.index(POSITION_X_COLUMN)
            return table.set_column(index, POSITION_X_COLUMN, pa.array(values))

        rewrite_table(archive, SCANS_MEMBER, edit)

        message = _refusal(archive)
        assert POSITION_X_COLUMN in message
        assert SCANS_MEMBER in message
        assert "whole number" in message

    def test_another_type_is_refused_by_name(self, archive):
        """A boolean position is refused, naming the column and its type."""
        rewrite_table(
            archive, SCANS_MEMBER, _retyped(POSITION_Y_COLUMN, kind=pa.bool_())
        )

        message = _refusal(archive)
        assert POSITION_Y_COLUMN in message
        assert "bool" in message


class TestMissingColumns:
    """A column the reader needs is refused by name; one it can do without is not."""

    @pytest.mark.parametrize(
        ("member", "column"),
        [(METADATA_MEMBER, "index"), (SCANS_MEMBER, "source_index")],
    )
    def test_a_join_column_is_required(self, archive, member, column):
        """Without it no spectrum can be placed (F46)."""
        rewrite_table(archive, member, _dropped(column))

        message = _refusal(archive)
        assert f"'{column}'" in message
        assert member in message
        assert str(archive) in message

    def test_no_observed_bounds_reads_them_from_the_points(self, archive):
        """Absent bound columns state no bound, as null ones do (F46)."""
        with MzPeakReader(archive) as reader:
            expected = reader.get_essential_metadata().mass_range
        rewrite_table(
            archive,
            METADATA_MEMBER,
            _dropped("lowest_observed_mz", "highest_observed_mz"),
        )

        with MzPeakReader(archive) as reader:
            assert reader.get_essential_metadata().mass_range == expected

    def test_no_representation_column_states_none(self, tmp_path):
        """No column and no file-level statement: the type is not stated (F44)."""
        archive = build_mzpeak(
            tmp_path / "a.mzpeak",
            grid_spectra(3, 2, n_points=5),
            file_contents=[{"name": "mass spectrum", "accession": "MS:1000294"}],
        )
        rewrite_table(archive, METADATA_MEMBER, _dropped("spectrum_representation"))

        with MzPeakReader(archive) as reader:
            assert reader.get_essential_metadata().spectrum_type is None
            assert len(list(reader.iter_spectra())) == 6


class TestDamagedMembers:
    """A member the ZIP holds intact but that cannot be read (F46)."""

    @pytest.mark.parametrize(
        ("member", "edit"),
        [
            (METADATA_MEMBER, lambda data: data[: len(data) // 2]),
            (DATA_MEMBER, lambda data: b""),
            (SCANS_MEMBER, lambda data: b"not parquet"),
        ],
        ids=["metadata cut in half", "data member empty", "scans not parquet"],
    )
    def test_refused_by_member(self, archive, member, edit):
        """The refusal names the archive and the member."""
        rewrite_member(archive, member, edit)

        message = _refusal(archive)
        assert member in message
        assert "damaged" in message

    def test_an_index_with_a_bad_crc(self, archive):
        """One byte of the stored index changed: its CRC no longer matches."""
        with zipfile.ZipFile(archive) as opened:
            stored = opened.read(INDEX_MEMBER)
        raw = bytearray(archive.read_bytes())
        start = raw.find(stored)
        assert start >= 0
        raw[start + 5] ^= 0x20
        archive.write_bytes(bytes(raw))

        message = _refusal(archive)
        assert INDEX_MEMBER in message


class TestPointLayoutFields:
    """A ``point`` struct whose fields are named otherwise (F45)."""

    @pytest.mark.parametrize(
        ("old", "new"),
        [
            ("intensity", "intensity_array"),
            ("mz", "mz_array"),
            ("spectrum_index", "index"),
        ],
    )
    @pytest.mark.parametrize(
        ("signal", "member"), [("profile", DATA_MEMBER), ("centroid", PEAKS_MEMBER)]
    )
    def test_refused_before_any_spectrum(self, tmp_path, old, new, signal, member):
        """Refused on opening, naming the missing field and the ones found."""
        archive = build_mzpeak(
            tmp_path / "a.mzpeak",
            grid_spectra(3, 2, n_points=4),
            layout="point",
            signal=signal,
        )

        def rename(table):
            kids = [
                pa.field(new if f.name == old else f.name, f.type)
                for f in table.schema.field("point").type
            ]
            return table.cast(pa.schema([pa.field("point", pa.struct(kids))]))

        rewrite_table(archive, member, rename)

        with pytest.raises(ConversionRefused) as refused:
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()
        message = str(refused.value)
        assert f"'{old}'" in message
        assert new in message
        assert member in message

    def test_the_extractor_alone_refuses_it_too(self, tmp_path):
        """The metadata document is made without the reader; it refuses too."""
        from thyra.metadata.extractors.mzpeak_extractor import MzPeakMetadataExtractor

        spectra = grid_spectra(3, 2, n_points=4)
        archive = build_mzpeak(tmp_path / "a.mzpeak", spectra, footer_metadata=False)
        rewrite_table(archive, METADATA_MEMBER, _dropped("lowest_observed_mz"))

        def rename(table):
            kids = [
                pa.field("mz_array" if f.name == "mz" else f.name, f.type)
                for f in table.schema.field("point").type
            ]
            return table.cast(pa.schema([pa.field("point", pa.struct(kids))]))

        rewrite_table(archive, DATA_MEMBER, rename)

        opened = mzpeak_reader.MzPeakArchive(archive)
        try:
            with pytest.raises(ConversionRefused, match="'mz'"):
                MzPeakMetadataExtractor(opened, archive).get_essential()
        finally:
            opened.close()


class TestUnusedPeaksMember:
    """The member that is not read cannot stop a conversion (F47)."""

    def test_an_empty_peaks_member_beside_profile_signal(self, tmp_path, thyra_logs):
        """The profile spectra convert; the log names the peaks member."""
        good = build_mzpeak(
            tmp_path / "good.mzpeak",
            grid_spectra(4, 3),
            signal="profile",
            empty_peer=True,
        )
        expected = _coordinates(good)
        rewrite_member(good, PEAKS_MEMBER, lambda data: b"")

        with thyra_logs(READER_LOGGER, logging.WARNING) as records:
            assert _coordinates(good) == expected
            with MzPeakReader(good) as reader:
                assert (
                    reader.get_essential_metadata().spectrum_type
                    == SpectrumType.PROFILE
                )

        assert PEAKS_MEMBER in "\n".join(records.messages)

    def test_a_damaged_member_that_is_read_is_refused(self, tmp_path):
        """A centroid archive's peaks member is its signal: refused by name."""
        archive = build_mzpeak(
            tmp_path / "c.mzpeak",
            grid_spectra(4, 3),
            signal="centroid",
            empty_peer=True,
        )
        rewrite_member(archive, PEAKS_MEMBER, lambda data: b"")

        assert PEAKS_MEMBER in _refusal(archive)

    def test_a_damaged_profile_member_beside_centroid_signal_is_refused(self, tmp_path):
        """Whether it held the signal cannot be told, so it is not passed over."""
        archive = build_mzpeak(
            tmp_path / "c.mzpeak",
            grid_spectra(4, 3),
            signal="centroid",
            empty_peer=True,
        )
        rewrite_member(archive, DATA_MEMBER, lambda data: b"")

        assert DATA_MEMBER in _refusal(archive)


class TestNoSignal:
    """An archive whose every spectrum is empty (F52)."""

    def test_refused_as_holding_no_signal(self, tmp_path):
        """One refusal, saying so, before any spectrum is read."""
        spectra = [Spectrum(x, y, [], []) for y in (1, 2) for x in (1, 2)]
        archive = build_mzpeak(tmp_path / "empty.mzpeak", spectra)

        with pytest.raises(ConversionRefused, match="holds no signal"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()


class TestChunkBatches:
    """A batch decodes about the target number of points (F142)."""

    @staticmethod
    def _largest_batch(tmp_path, monkeypatch, keep: int) -> int:
        """Rows of the largest batch when every keep-th spectrum has a position."""
        monkeypatch.setattr(mzpeak_reader, "CHUNK_BATCH_POINTS", 10_000)
        rows = []
        decode = mzpeak_reader.decode_chunks

        def counting(column, path):
            rows.append(len(column))
            return decode(column, path)

        monkeypatch.setattr(mzpeak_reader, "decode_chunks", counting)
        mz = 100.0 + np.arange(400) * 0.01
        spectra = [
            Spectrum(
                1 + i % 20 if i % keep == 0 else None,
                1 + i // 20 if i % keep == 0 else None,
                mz,
                np.ones(400),
            )
            for i in range(400)
        ]
        archive = build_mzpeak(
            tmp_path / f"keep{keep}.mzpeak", spectra, layout="chunk", chunk_points=100
        )
        with MzPeakReader(archive) as reader:
            list(reader.iter_spectra())
        return max(rows)

    def test_the_batch_does_not_grow_as_spectra_are_left_out(
        self, tmp_path, monkeypatch
    ):
        """Every spectrum kept, or one in a hundred: the same batch."""
        every = self._largest_batch(tmp_path, monkeypatch, keep=1)
        few = self._largest_batch(tmp_path, monkeypatch, keep=100)

        assert every == 100
        assert few == every


class TestLossyAxisWarning:
    """The warning says what can happen, not what did (F84)."""

    def test_worded_as_a_possibility(self, tmp_path, thyra_logs):
        """The archive does not say how many m/z its source held."""
        archive = build_mzpeak(
            tmp_path / "n.mzpeak",
            grid_spectra(3, 2),
            layout="chunk",
            chunk_encoding="MS:1002312",
        )
        with thyra_logs(READER_LOGGER, logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                reader.get_common_mass_axis()

        logged = "\n".join(records.messages)
        assert "can then grow far longer" in logged
        assert "therefore" not in logged
