"""An imzML leaves an MSn spectrum out only on a pixel that holds MS1 too (D70).

The mzPeak reader has applied this rule since D70; the imzML route reads
``MS:1000511`` on the spectrum or in a group it references, and applies the
same rule, so a file and its mzPeak archive convert to the same store
(issue #522, F26).
"""

import logging
import shutil
from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

from tests.fixtures.imzml_levels import level_term, write_imzml_with_levels
from thyra.errors import ConversionRefused
from thyra.readers.imzml.imzml_reader import ImzMLReader


def _read(path):
    """What the conversion route takes from a file."""
    with ImzMLReader(path) as reader:
        spectra = list(reader.iter_spectra_with_acquisition_order())
        return {
            "spectra": spectra,
            "orders": [order for _, order, _, _ in spectra],
            "n_spectra": reader.n_spectra,
            "peaks": reader.get_total_peak_count(),
            "axis": reader.get_common_mass_axis().tolist(),
            "essential": reader.get_essential_metadata(),
        }


def _give_spectrum_0_its_own_mz_block(path, mzs):
    """Point spectrum 0's m/z array at new float64 values after the .ibd's end.

    A continuous file whose spectrum 0 reads another block than the rest:
    a reader that takes the shared axis from spectrum 0 would show it.
    """
    ibd = path.with_suffix(".ibd")
    offset = ibd.stat().st_size
    with open(ibd, "ab") as handle:
        handle.write(np.asarray(mzs, dtype="<f8").tobytes())
    text = path.read_text()
    start = text.index("<spectrum ")
    marker = 'name="external offset" value="'
    at = text.index(marker, text.index('ref="mzArray"', start)) + len(marker)
    end = text.index('"', at)
    path.write_text(text[:at] + str(offset) + text[end:])


def _msn_messages(records):
    return [r.getMessage() for r in records if "MSn spectra" in r.getMessage()]


@pytest.mark.parametrize("how", ["own", "group"])
def test_only_the_shared_pixel_loses_its_ms2(tmp_path, how, thyra_logs):
    """MS1 and MS2 on one pixel, MS2 alone on the next, MS1 after it."""
    spectra = [
        (1, 1, 1, [100.0, 110.0], [1.0, 2.0]),
        # Outside the kept spectra's range, so the range would show it.
        (1, 1, 1, [900.0, 901.0, 902.0], [1000.0, 1.0, 1.0]),
        (2, 1, 1, [160.0], [10.0]),
        (3, 1, 1, [120.0], [4.0]),
    ]
    path = write_imzml_with_levels(
        tmp_path / "mixed.imzML", spectra, [1, 2, 2, 1], how=how
    )

    with thyra_logs("thyra.readers.imzml", logging.WARNING) as records:
        reading = _read(path)

    # The left-out spectrum keeps its place: the one after it is still 2.
    assert reading["orders"] == [0, 2, 3]
    assert [coords for coords, _, _, _ in reading["spectra"]] == [
        (0, 0, 0),
        (1, 0, 0),
        (2, 0, 0),
    ]
    np.testing.assert_array_equal(reading["spectra"][0][3], [1.0, 2.0])
    assert reading["n_spectra"] == 3
    assert reading["peaks"] == 4
    assert reading["axis"] == [100.0, 110.0, 120.0, 160.0]
    essential = reading["essential"]
    assert essential.n_spectra == 3
    assert essential.total_peaks == 4
    assert essential.mass_range == (100.0, 160.0)
    assert essential.dimensions == (3, 1, 1)
    messages = _msn_messages(records)
    assert len(messages) == 1
    assert "1 MSn spectra (level 2: 1)" in messages[0]
    assert "1 MSn spectra on pixels of their own are kept" in messages[0]


@pytest.mark.parametrize(
    "first_level", [1, 0, None], ids=["level-1", "level-0", "level-none"]
)
def test_an_ms1_spectrum_of_any_stated_kind_keeps_ms2_out(tmp_path, first_level):
    """Level 0 or no level beside level 2 reads as MS1."""
    spectra = []
    for x in (1, 2):
        spectra += [
            (x, 1, 1, [100.0, 101.0], [1.0, 2.0]),
            (x, 1, 1, [150.0], [1000.0]),
        ]
    path = write_imzml_with_levels(
        tmp_path / "pairs.imzML", spectra, [first_level, 2] * 2
    )

    reading = _read(path)

    assert reading["n_spectra"] == 2
    assert reading["axis"] == [100.0, 101.0]
    for _, _, mzs, intensities in reading["spectra"]:
        np.testing.assert_array_equal(mzs, [100.0, 101.0])
        assert float(intensities.sum()) == 3.0


def test_the_writer_s_group_level_0_reads_as_ms1(tmp_path):
    """pyimzml states level 0 in a group; a group of level 2 beside it is MSn."""
    spectra = [
        (1, 1, 1, [100.0], [1.0]),
        (1, 1, 1, [150.0], [1000.0]),
    ]
    path = write_imzml_with_levels(
        tmp_path / "writer.imzML", spectra, [None, 2], how="group"
    )

    reading = _read(path)

    assert reading["orders"] == [0]
    assert reading["axis"] == [100.0]


def test_a_spectrum_s_own_level_comes_before_its_group_s(tmp_path):
    """An own level 1 beside a referenced group of level 2 is MS1."""
    spectra = [
        (1, 1, 1, [100.0], [1.0]),
        (1, 1, 1, [150.0], [1000.0]),
    ]
    path = write_imzml_with_levels(
        tmp_path / "own.imzML", spectra, [None, 2], how="group"
    )
    text = path.read_text()
    group_ref = '<referenceableParamGroupRef ref="level2"/>'
    assert text.count(group_ref) == 1
    path.write_text(text.replace(group_ref, group_ref + level_term(1)))

    assert _read(path)["orders"] == [0, 1]


def test_a_pixel_on_another_plane_is_another_pixel(tmp_path):
    """MS1 on plane 1 does not leave out MS2 at the same x, y on plane 2."""
    spectra = [
        (1, 1, 1, [100.0], [1.0]),
        (1, 1, 2, [150.0], [1000.0]),
    ]
    path = write_imzml_with_levels(tmp_path / "planes.imzML", spectra, [1, 2])

    reading = _read(path)

    assert reading["orders"] == [0, 1]
    assert reading["axis"] == [100.0, 150.0]


def test_continuous_mode_counts_and_sums_the_kept_spectra(tmp_path):
    """The F26 shape of issue #522: TIC 40 and 80, not 40 and 1280."""
    mzs = [150.0, 250.0, 350.0, 450.0]
    spectra = [
        (1, 1, 1, mzs, [10.0] * 4),
        (2, 1, 1, mzs, [20.0] * 4),
        (2, 1, 1, mzs, [0.0, 500.0, 0.0, 700.0]),
    ]
    path = write_imzml_with_levels(
        tmp_path / "f26.imzML", spectra, [1, 1, 2], mode="continuous"
    )

    reading = _read(path)

    assert [float(i.sum()) for _, _, _, i in reading["spectra"]] == [40.0, 80.0]
    assert reading["n_spectra"] == 2
    assert reading["peaks"] == 8
    assert reading["essential"].n_spectra == 2
    assert reading["essential"].total_peaks == 8


def test_a_continuous_axis_comes_from_the_first_kept_spectrum(tmp_path):
    """An MSn spectrum left out at position 0 gives neither axis nor range."""
    spectra = [
        (1, 1, 1, [100.0, 101.0], [2.0, 2.0]),
        (1, 1, 1, [100.0, 101.0], [1.0, 1.0]),
        (2, 1, 1, [100.0, 101.0], [3.0, 3.0]),
    ]
    path = write_imzml_with_levels(
        tmp_path / "first.imzML", spectra, [2, 1, 1], mode="continuous"
    )
    _give_spectrum_0_its_own_mz_block(path, [500.0, 501.0])

    reading = _read(path)

    assert reading["orders"] == [1, 2]
    assert reading["axis"] == [100.0, 101.0]
    assert reading["essential"].n_spectra == 2
    assert reading["essential"].total_peaks == 4


def test_a_left_out_msn_spectrum_on_its_own_block_keeps_the_one_block_read(
    tmp_path, thyra_logs
):
    """The measured block test counts converted spectra only (D73, D82)."""
    spectra = [
        (1, 1, 1, [100.0, 101.0], [2.0, 2.0]),
        (1, 1, 1, [100.0, 101.0], [1.0, 1.0]),
        (2, 1, 1, [100.0, 101.0], [3.0, 3.0]),
    ]
    path = write_imzml_with_levels(
        tmp_path / "first.imzML", spectra, [2, 1, 1], mode="continuous"
    )
    _give_spectrum_0_its_own_mz_block(path, [500.0, 501.0])

    with thyra_logs("thyra.readers.imzml.imzml_reader") as records:
        with ImzMLReader(path) as reader:
            assert reader.has_shared_mass_axis is True
            assert reader._continuous_uniform
            mass_range = reader.get_essential_metadata().mass_range

    assert mass_range == (100.0, 101.0)
    assert not any("declares continuous mode" in r.getMessage() for r in records)


def _mobility_file_with_a_left_out_msn_spectrum(
    tmp_path, drop_its_array=False, shorten_its_array=False, leave_it_in=False
):
    """The processed mobility fixture, its fourth spectrum moved onto the
    second's pixel and stated MS2; optionally without its mobility array, or
    with one shorter than its m/z array. ``leave_it_in`` keeps the spectrum
    where it was and states no level, so it is converted."""
    fixtures = Path(__file__).resolve().parents[2] / "data" / "fixtures"
    path = tmp_path / "mobility.imzML"
    shutil.copyfile(fixtures / "mobility_processed.ibd", path.with_suffix(".ibd"))
    text = (fixtures / "mobility_processed.imzML").read_text()
    head, fourth = text.split('<spectrum id="Scan=4"')
    if not leave_it_in:
        fourth = fourth.replace(
            '<referenceableParamGroupRef ref="spectrum1" />',
            '<referenceableParamGroupRef ref="spectrum1" />' + level_term(2),
            1,
        ).replace('name="position y" value="2"', 'name="position y" value="1"', 1)
    mobility = fourth.index(
        '<binaryDataArray encodedLength="0"><referenceableParamGroupRef '
        'ref="mobilityArray" />'
    )
    if shorten_its_array:
        length = 'name="external array length" value="4"'
        at = fourth.index(length, mobility)
        fourth = fourth[:at] + length.replace('"4"', '"3"') + fourth[at + len(length) :]
    if drop_its_array:
        end = fourth.index("</binaryDataArray>", mobility) + len("</binaryDataArray>")
        fourth = fourth[:mobility] + fourth[end:]
        fourth = fourth.replace(
            '<binaryDataArrayList count="3">', '<binaryDataArrayList count="2">', 1
        )
    path.write_text(head + '<spectrum id="Scan=4"' + fourth)
    return path


def test_a_left_out_msn_spectrum_needs_no_mobility_array(tmp_path):
    """It is never read, so the file is not refused for lacking one there."""
    path = _mobility_file_with_a_left_out_msn_spectrum(tmp_path, drop_its_array=True)

    with ImzMLReader(path) as reader:
        items = list(reader.iter_mobility_spectra())

    assert [coords for coords, _, _, _ in items] == [(0, 0, 0), (1, 0, 0), (0, 1, 0)]


def test_a_left_out_msn_spectrum_may_hold_a_mobility_array_of_another_length(
    tmp_path,
):
    """Its arrays are not checked for being parallel, as it is never read."""
    path = _mobility_file_with_a_left_out_msn_spectrum(tmp_path, shorten_its_array=True)

    with ImzMLReader(path) as reader:
        items = list(reader.iter_mobility_spectra())

    assert len(items) == 3


@pytest.mark.parametrize(
    "change, refusal",
    [
        ({"drop_its_array": True}, "1 of 4 spectra reference no ion mobility array"),
        ({"shorten_its_array": True}, "Spectrum 3 declares 3 ion mobility values"),
    ],
    ids=["missing", "shorter"],
)
def test_a_converted_spectrum_s_mobility_array_is_still_checked(
    tmp_path, change, refusal
):
    """The same file with that spectrum converted is refused, as before."""
    path = _mobility_file_with_a_left_out_msn_spectrum(
        tmp_path, leave_it_in=True, **change
    )

    with pytest.raises(ConversionRefused, match=refusal):
        ImzMLReader(path).n_spectra


def test_mobility_spectra_and_range_leave_the_msn_spectrum_out(tmp_path, monkeypatch):
    """The moved MS2 spectrum gives no point cloud and no range."""
    path = _mobility_file_with_a_left_out_msn_spectrum(tmp_path)

    with ImzMLReader(path) as reader:
        original = reader._read_mobility_array

        def _wide_on_the_msn_spectrum(idx):
            values = original(idx)
            return np.full(values.size, 9.0) if idx == 3 else values

        monkeypatch.setattr(reader, "_read_mobility_array", _wide_on_the_msn_spectrum)
        items = list(reader.iter_mobility_spectra())
        span = reader.get_mobility_value_range()

    assert [coords for coords, _, _, _ in items] == [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    for _, _, mobility, _ in items:
        np.testing.assert_array_equal(mobility, [0.8, 1.1, 0.9, 1.0])
    assert span == (0.8, 1.1)


@pytest.mark.parametrize(
    "levels", [None, [1, 1, 1], [2, 2, 2]], ids=["none", "all-ms1", "all-ms2"]
)
def test_a_file_without_ms1_beside_msn_keeps_every_spectrum(
    tmp_path, levels, thyra_logs
):
    """No levels, only MS1, or only MS2: nothing is left out or logged."""
    spectra = [
        (1, 1, 1, [100.0], [1.0]),
        (1, 1, 1, [150.0], [2.0]),
        (2, 1, 1, [160.0], [3.0]),
    ]
    path = tmp_path / "plain.imzML"
    if levels is None:
        with ImzMLWriter(str(path), mode="processed") as writer:
            for x, y, z, mzs, intensities in spectra:
                writer.addSpectrum(mzs, intensities, (x, y, z))
    else:
        write_imzml_with_levels(path, spectra, levels)

    with thyra_logs("thyra.readers.imzml", logging.WARNING) as records:
        reading = _read(path)

    assert reading["n_spectra"] == 3
    assert reading["orders"] == [0, 1, 2]
    assert reading["axis"] == [100.0, 150.0, 160.0]
    assert reading["essential"].n_spectra == 3
    assert not _msn_messages(records)
