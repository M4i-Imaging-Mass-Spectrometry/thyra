"""The decoders of the mzPeak chunked layout, one encoding at a time.

Two kinds of evidence are kept here, because each covers what the other
cannot.

The rows under ``CONVERTER_*`` were written by someone else's writer:
mzpeak-convert 0.14.0 (okohlbacher/mzPeakConverter @ 0ed311e) converted
synthetic spectra, made for this purpose, once in the chunked layout and
once in the point layout. A row is one chunk as that converter stored it,
beside the m/z the point layout holds for the same points. They show the
decoders agree with a writer that is not ours. The rows under
``CONVERTER_TIMSTOF_*`` are two chunks mzpeak-convert 0.17.2 wrote from the
public timsTOF fleX TDF run of MassIVE MSV000088438, from two frames with
different calibrations; the converter states the m/z of each row's first
and last bin, evaluated by its own reference implementation, as the row's
bounds.

The other tests encode with :mod:`tests.fixtures.mzpeak_builder`, which
follows the specification's text and asks the decoder nothing. They reach
the shapes a short real row does not: nulls, negative residuals, rows of
one point, counts that disagree.
"""

from __future__ import annotations

import struct

import numpy as np
import pyarrow as pa
import pytest

from tests.fixtures.mzpeak_builder import (
    DELTA,
    GRID,
    GRID_PARAMETERS,
    LINEAR_GRID,
    NUMPRESS_LINEAR,
    PLAIN,
    SQUARE_ROOT_GRID,
    Spectrum,
    chunk_table,
    delta_encode,
    grid_indices,
    numpress_linear_encode,
)
from thyra.errors import ConversionRefused
from thyra.readers.mzpeak import chunk_decoding
from thyra.readers.mzpeak.chunk_decoding import (
    decode_chunks,
    decode_delta,
    decode_grid,
    decode_numpress_linear,
    decode_plain,
    running_sums,
    take_bounds,
    timstof_mz,
)

#: The timsTOF m/z grid model of the reference implementation.
TIMSTOF_MZ_GRID = "MS:9999002"

SOURCE = "fixture.mzpeak"

#: One MS-Numpress chunk of 14 points, as the converter stored it.
CONVERTER_NUMPRESS_BYTES = [
    65, 114, 224, 183, 48, 0, 0, 0, 67, 231, 109, 123, 209, 255, 255, 127,
    16, 210, 112, 18, 157, 10, 153, 236, 47, 29, 206, 193, 213, 168, 108, 73,
    52, 138, 221, 122, 198, 144, 78, 40, 38, 237, 50, 212, 141, 125, 41, 195,
    156, 25, 176, 83, 146, 225, 238, 218, 242, 97, 33, 202, 182, 32,
]  # fmt: skip
CONVERTER_NUMPRESS_MZ = [
    104.61332681135059, 108.48724287956014, 114.11072932850564,
    117.11739053647167, 120.80872279913262, 128.5476028824506,
    129.39232936182776, 130.14447617289295, 131.10155636964345,
    132.77324899631446, 134.5395544570974, 134.74717413998872,
    140.19798894100046, 147.7003944338548,
]  # fmt: skip

#: One delta-encoded chunk of 12 points, as the converter stored it.
CONVERTER_DELTA_START = 102.6591
CONVERTER_DELTA_VALUES = [
    2.0752999999999986, 4.961500000000001, 2.2276999999999987,
    2.082400000000007, 0.7592999999999961, 8.659900000000007,
    0.6323999999999899, 4.766200000000012, 7.724199999999996,
    2.1189999999999998, 7.919900000000013,
]  # fmt: skip
CONVERTER_DELTA_MZ = [
    102.6591, 104.7344, 109.6959, 111.9236, 114.006, 114.7653,
    123.4252, 124.0576, 128.8238, 136.548, 138.667, 146.5869,
]  # fmt: skip

#: One grid-encoded chunk of 12 points, as the converter stored it.
CONVERTER_GRID_PARAMETERS = [0.049890254999999994, 5.788144191212053e-12, 0.0005]
CONVERTER_GRID_INDICES = [
    195225786, 36807998, 47908275, 174062699, 535171533, 1110070825,
    735192466, 126724555, 180472353, 161709862, 525652074, 270743083,
]  # fmt: skip
CONVERTER_GRID_MZ = [
    102.0405, 102.4666, 103.0212, 105.0362, 111.2315, 124.082,
    132.5928, 134.0598, 136.149, 138.021, 144.1061, 147.2403,
]  # fmt: skip


#: Two timsTOF m/z grid chunks: parameters, indices, and the bounds the
#: converter evaluated for the first and the last bin.
CONVERTER_TIMSTOF_ROWS = [
    (
        [
            316.8515891632165, 2531.224084396041, 7.555625863514373e-05, 0.0,
            -0.002708343104990962, 0.19999999999999998, 25628.6,
        ],
        [3732, 1347, 847, 2958, 4774, 304, 3275, 4145, 1856, 2892, 1636, 3051],
        (105.98313200667917, 154.62569669420392),
    ),
    (
        [
            316.8515891632165, 2531.2240816152676, 7.555625846913327e-05, 0.0,
            -0.002708343104990962, 0.19999999999999998, 25628.6,
        ],
        [294660, 80, 788, 1906, 324, 1455, 2885, 6, 1311, 1],
        (1107.680159057425, 1154.2099969562803),
    ),
]  # fmt: skip


def _numpress(rows):
    """Encode rows with the builder and lay them end to end."""
    buffers = [numpress_linear_encode(list(row)) for row in rows]
    data = np.frombuffer(b"".join(buffers), dtype=np.uint8)
    lengths = np.array([len(b) for b in buffers], dtype=np.int64)
    counts = np.array([len(row) for row in rows], dtype=np.int64)
    return data, lengths, counts


def _chunk(spectra, **kwargs):
    """The ``chunk`` struct column of a fixture table."""
    return chunk_table(spectra, **kwargs).column("chunk").combine_chunks()


class TestAgainstTheConverter:
    """Chunks written by mzpeak-convert decode to the point layout's m/z."""

    def test_numpress_linear(self):
        """Lossy by design; the converter's own bound here is 1e-6."""
        data = np.array(CONVERTER_NUMPRESS_BYTES, dtype=np.uint8)
        decoded = decode_numpress_linear(
            data, np.array([data.size]), np.array([14]), SOURCE
        )
        np.testing.assert_allclose(decoded, CONVERTER_NUMPRESS_MZ, rtol=0, atol=1e-7)

    def test_delta(self):
        """The sums are taken in stored order, so the values come back whole."""
        decoded = decode_delta(
            np.array([CONVERTER_DELTA_START]),
            np.array(CONVERTER_DELTA_VALUES),
            np.array([11]),
            np.array([12]),
            SOURCE,
        )
        np.testing.assert_array_equal(decoded, CONVERTER_DELTA_MZ)

    def test_linear_grid(self):
        """The converter promises every value within 1e-6 of the source."""
        decoded = decode_grid(
            LINEAR_GRID,
            (np.array(CONVERTER_GRID_PARAMETERS), np.array([3])),
            (np.array(CONVERTER_GRID_INDICES, dtype=np.uint32), np.array([12])),
            np.array([12]),
            SOURCE,
        )
        np.testing.assert_allclose(decoded, CONVERTER_GRID_MZ, rtol=0, atol=1e-7)

    def test_timstof_mz_grid(self):
        """Each row's ends decode to the bounds the converter evaluated.

        The two rows sit side by side with different parameters, so each
        must be decoded under its own.
        """
        parameters = [value for row in CONVERTER_TIMSTOF_ROWS for value in row[0]]
        indices = [value for row in CONVERTER_TIMSTOF_ROWS for value in row[1]]
        counts = np.array([len(row[1]) for row in CONVERTER_TIMSTOF_ROWS])

        decoded = decode_grid(
            TIMSTOF_MZ_GRID,
            (np.array(parameters), np.array([7, 7])),
            (np.array(indices, dtype=np.uint32), counts),
            counts,
            SOURCE,
        )

        ends = [decoded[0], decoded[11], decoded[12], decoded[21]]
        bounds = [b for row in CONVERTER_TIMSTOF_ROWS for b in row[2]]
        np.testing.assert_allclose(ends, bounds, rtol=1e-15, atol=0)
        assert np.all(np.diff(decoded[:12]) > 0)
        assert np.all(np.diff(decoded[12:]) > 0)


class TestNumpressLinear:
    """MS-Numpress linear prediction, encoded by the builder."""

    @pytest.mark.parametrize("count", [0, 1, 2, 3, 4, 17, 300])
    def test_every_length_round_trips(self, count):
        """Buffers of no value, one, two, and an odd or even tail."""
        values = 100.0 + 0.25 * np.arange(count) ** 1.5
        values = np.round(values * 1024) / 1024
        data, lengths, counts = _numpress([values])

        decoded = decode_numpress_linear(data, lengths, counts, SOURCE)

        np.testing.assert_array_equal(decoded, values)

    def test_rows_of_unequal_length_decode_side_by_side(self):
        """Many buffers at once give what each gives alone."""
        rng = np.random.default_rng(188)
        rows = [
            np.sort(np.round(rng.uniform(50, 2000, n) * 4096) / 4096)
            for n in rng.integers(0, 60, 40)
        ]
        data, lengths, counts = _numpress(rows)

        decoded = decode_numpress_linear(data, lengths, counts, SOURCE)

        np.testing.assert_array_equal(decoded, np.concatenate(rows))

    def test_residuals_of_either_sign_and_every_width(self):
        """Steps that grow and shrink, by one unit and by millions."""
        steps = [0, 1, -1, 15, -16, 255, -256, 70000, -70000, 9000000, -9000000]
        values = 1000.0 + np.cumsum(np.cumsum(np.array(steps) / 1024.0))
        data, lengths, counts = _numpress([values])

        decoded = decode_numpress_linear(data, lengths, counts, SOURCE)

        np.testing.assert_array_equal(decoded, values)

    def test_a_stored_zero_is_padding(self):
        """MS-Numpress has no null, so the writer's zero stands for one."""
        data, lengths, counts = _numpress([[100.0, 0.0, 0.0, 101.5, 102.0]])

        decoded = decode_numpress_linear(data, lengths, counts, SOURCE)

        np.testing.assert_array_equal(np.isnan(decoded), [0, 1, 1, 0, 0])
        np.testing.assert_array_equal(decoded[[0, 3, 4]], [100.0, 101.5, 102.0])

    def test_the_scale_is_read_big_endian(self):
        """The first eight bytes are the scale, most significant first."""
        buffer = numpress_linear_encode([3.0, 4.0], fixed_point=1024.0)

        assert struct.unpack(">d", buffer[:8]) == (1024.0,)
        decoded = decode_numpress_linear(
            np.frombuffer(buffer, dtype=np.uint8),
            np.array([len(buffer)]),
            np.array([2]),
            SOURCE,
        )
        np.testing.assert_array_equal(decoded, [3.0, 4.0])

    @pytest.mark.parametrize("claimed", [4, 6])
    def test_a_count_that_disagrees_is_refused(self, claimed):
        """Five values stored, another number of intensities beside them."""
        data, lengths, _ = _numpress([[100.0, 101.0, 102.5, 103.0, 105.0]])

        with pytest.raises(ConversionRefused, match="counts disagree.*MS:1002312"):
            decode_numpress_linear(data, lengths, np.array([claimed]), SOURCE)

    def test_a_buffer_cut_short_is_refused(self):
        """A header that is not all there cannot be read as values."""
        data, _, counts = _numpress([[100.0, 101.0, 102.5]])

        with pytest.raises(ConversionRefused, match="too short"):
            decode_numpress_linear(data[:14], np.array([14]), counts, SOURCE)


class TestDelta:
    """Delta encoding, nulls included."""

    def test_sums_are_taken_in_stored_order(self):
        """Each value is the one before plus its difference, nothing else."""
        values = [100.1, 100.30000000000001, 100.7, 103.9]
        encoded = delta_encode(values)
        expected = np.array(values[:1] + encoded, dtype=np.float64)
        for place in range(1, expected.size):
            expected[place] = expected[place - 1] + expected[place]

        decoded = decode_delta(
            np.array([values[0]]),
            np.array(encoded),
            np.array([3]),
            np.array([4]),
            SOURCE,
        )

        np.testing.assert_array_equal(decoded, expected)

    def test_a_value_after_a_null_is_taken_as_it_is(self):
        """Differences stop at a null pair and resume after it."""
        values = [100.0, 100.5, None, None, 200.0, 200.25]
        encoded = np.array(delta_encode(values), dtype=np.float64)

        decoded = decode_delta(
            np.array([100.0]), encoded, np.array([5]), np.array([6]), SOURCE
        )

        np.testing.assert_array_equal(
            decoded, [100.0, 100.5, np.nan, np.nan, 200.0, 200.25]
        )

    def test_a_row_that_opens_with_a_null_keeps_it_in_the_list(self):
        """The list is then as long as the intensities, start not added."""
        values = [None, 150.0, 150.5]
        encoded = np.array(delta_encode(values), dtype=np.float64)
        assert encoded.size == 3

        decoded = decode_delta(
            np.array([150.0]), encoded, np.array([3]), np.array([3]), SOURCE
        )

        np.testing.assert_array_equal(decoded, [np.nan, 150.0, 150.5])

    def test_a_row_that_opens_with_a_null_pair(self):
        """The reference decoder adds the start here; the count says not to."""
        values = [None, None, 150.0, 150.5]
        encoded = np.array(delta_encode(values), dtype=np.float64)

        decoded = decode_delta(
            np.array([150.0]), encoded, np.array([4]), np.array([4]), SOURCE
        )

        np.testing.assert_array_equal(decoded, [np.nan, np.nan, 150.0, 150.5])

    def test_a_start_followed_by_a_null_pair(self):
        """A lone first point, then a removed run of zeros."""
        values = [100.0, None, None, 150.0]
        encoded = np.array(delta_encode(values), dtype=np.float64)

        decoded = decode_delta(
            np.array([100.0]), encoded, np.array([3]), np.array([4]), SOURCE
        )

        np.testing.assert_array_equal(decoded, [100.0, np.nan, np.nan, 150.0])

    def test_a_row_of_one_point_has_an_empty_list(self):
        """The start is the point."""
        decoded = decode_delta(
            np.array([321.5]), np.array([]), np.array([0]), np.array([1]), SOURCE
        )

        np.testing.assert_array_equal(decoded, [321.5])

    def test_a_count_that_disagrees_is_refused(self):
        """Three differences make four values, not six."""
        with pytest.raises(ConversionRefused, match="counts disagree.*MS:1003089"):
            decode_delta(
                np.array([100.0]),
                np.array([0.5, 0.5, 0.5]),
                np.array([3]),
                np.array([6]),
                SOURCE,
            )

    def test_short_and_long_runs_sum_alike(self, monkeypatch):
        """The two routes through the running sum give the same bits."""
        rng = np.random.default_rng(7)
        lengths = rng.integers(1, 40, 60)
        starts = np.cumsum(lengths) - lengths
        values = rng.uniform(0.001, 3.0, int(lengths.sum()))

        side_by_side = values.copy()
        running_sums(side_by_side, starts, lengths)
        monkeypatch.setattr(chunk_decoding, "SHORT_RUN", 0)
        one_by_one = values.copy()
        running_sums(one_by_one, starts, lengths)

        expected = np.concatenate(
            [np.cumsum(values[s : s + n]) for s, n in zip(starts, lengths)]
        )
        np.testing.assert_array_equal(side_by_side, expected)
        np.testing.assert_array_equal(one_by_one, expected)


class TestPlain:
    """Values stored as they are."""

    def test_the_start_is_put_back_in_front(self):
        """``mz_chunk_values`` leaves the first value to the start."""
        decoded = decode_plain(
            np.array([100.0, 200.0]),
            np.array([100.5, 101.0, 200.5]),
            np.array([2, 1]),
            np.array([3, 2]),
            SOURCE,
        )

        np.testing.assert_array_equal(decoded, [100.0, 100.5, 101.0, 200.0, 200.5])

    def test_a_count_that_disagrees_is_refused(self):
        """Two values and a start make three points, not two."""
        with pytest.raises(ConversionRefused, match="counts disagree.*MS:1000576"):
            decode_plain(
                np.array([100.0]),
                np.array([100.5, 101.0]),
                np.array([2]),
                np.array([2]),
                SOURCE,
            )


class TestGrid:
    """Grid indices and the model that turns them into m/z."""

    def test_linear_grid(self):
        """``(intercept + slope * index) / scale``, on the worked example.

        The specification gives intercept 95, slope 3.75e-7, scale 1 and
        index 2190583200 for 916.4687.
        """
        decoded = decode_grid(
            LINEAR_GRID,
            (np.array([95.0, 3.75e-7, 1.0]), np.array([3])),
            (np.array([2190583200], dtype=np.uint32), np.array([1])),
            np.array([1]),
            SOURCE,
        )

        np.testing.assert_allclose(decoded, [916.4687], rtol=0, atol=1e-9)

    def test_square_root_grid(self):
        """The square of the bracket, over the scale."""
        values = [(10.0 + k / 4.0) ** 2 for k in range(6)]
        indices = grid_indices(values, SQUARE_ROOT_GRID, GRID_PARAMETERS)

        decoded = decode_grid(
            SQUARE_ROOT_GRID,
            (np.array(GRID_PARAMETERS), np.array([3])),
            (np.array(indices, dtype=np.uint32), np.array([6])),
            np.array([6]),
            SOURCE,
        )

        np.testing.assert_array_equal(decoded, values)

    def test_indices_are_steps_and_restart_at_every_row(self):
        """The first index of a row is whole; the rest are differences."""
        parameters = np.tile(np.array(GRID_PARAMETERS), 2)
        first = grid_indices([100.0, 100.5, 102.0], LINEAR_GRID, GRID_PARAMETERS)
        second = grid_indices([50.0, 51.25], LINEAR_GRID, GRID_PARAMETERS)

        decoded = decode_grid(
            LINEAR_GRID,
            (parameters, np.array([3, 3])),
            (np.array(first + second, dtype=np.uint32), np.array([3, 2])),
            np.array([3, 2]),
            SOURCE,
        )

        np.testing.assert_array_equal(decoded, [100.0, 100.5, 102.0, 50.0, 51.25])

    def test_a_scale_left_out_is_one(self):
        """Two parameters are an intercept and a slope."""
        decoded = decode_grid(
            LINEAR_GRID,
            (np.array([10.0, 0.5]), np.array([2])),
            (np.array([4, 2], dtype=np.uint32), np.array([2])),
            np.array([2]),
            SOURCE,
        )

        np.testing.assert_array_equal(decoded, [12.0, 13.0])

    def test_a_count_that_disagrees_is_refused(self):
        """Two indices beside three intensities."""
        with pytest.raises(ConversionRefused, match="counts disagree.*MS:1003826"):
            decode_grid(
                LINEAR_GRID,
                (np.array(GRID_PARAMETERS), np.array([3])),
                (np.array([4, 2], dtype=np.uint32), np.array([2])),
                np.array([3]),
                SOURCE,
            )


class TestTimstofMz:
    """The timsTOF m/z grid: a TOF bin solved through the TDF calibration.

    The flight time of a bin is ``index * timebase + delay``. The model
    states it as ``C0 + beta * u + C2 * u**2 + C3 * u**3`` of
    ``u = sqrt(m/z + C4)``. Each test checks the decoded m/z against that
    statement, not against how it is solved.
    """

    #: C0, beta, C2, C3, C4, timebase, delay, of the shape real runs have.
    PARAMETERS = [316.85, 2531.22, 7.5556e-05, 0.0, -0.0027, 0.2, 25628.6]

    @staticmethod
    def _flight(mz, parameters):
        c0, beta, c2, c3, c4 = parameters[:5]
        u = np.sqrt(np.asarray(mz) + c4)
        return c0 + beta * u + c2 * u**2 + c3 * u**3

    def _check(self, parameters, bins):
        decoded = timstof_mz(np.asarray(bins, dtype=np.float64), np.array(parameters))
        flight = np.asarray(bins, dtype=np.float64) * parameters[5] + parameters[6]
        np.testing.assert_allclose(
            self._flight(decoded, parameters), flight, rtol=1e-14, atol=0
        )
        return decoded

    @pytest.mark.parametrize("c2", [7.5556e-05, 0.0, -2e-05])
    def test_quadratic(self, c2):
        """Closed form, with the quadratic term of either sign or none."""
        parameters = list(self.PARAMETERS)
        parameters[2] = c2
        decoded = self._check(parameters, [0, 1, 3732, 299_999, 400_000])
        assert np.all(np.diff(decoded) > 0)

    def test_cubic(self):
        """A cubic term is solved by Newton's method from the linear solution."""
        parameters = list(self.PARAMETERS)
        parameters[3] = 1e-9
        self._check(parameters, [100, 50_000, 250_000, 400_000])

    def test_rows_take_their_own_parameters(self):
        """A frame's temperature correction changes beta and C2 per row."""
        other = list(self.PARAMETERS)
        other[1] *= 1.00002
        decoded = decode_grid(
            TIMSTOF_MZ_GRID,
            (np.array(self.PARAMETERS + other), np.array([7, 7])),
            (np.array([50_000, 50_000], dtype=np.uint32), np.array([1, 1])),
            np.array([1, 1]),
            SOURCE,
        )

        assert decoded[0] != decoded[1]
        np.testing.assert_array_equal(
            decoded,
            [
                self._check(self.PARAMETERS, [50_000])[0],
                self._check(other, [50_000])[0],
            ],
        )

    def test_too_few_parameters_are_refused(self):
        """The model needs all seven; three are a linear grid's."""
        with pytest.raises(ConversionRefused, match="3 parameters; the model needs 7"):
            decode_grid(
                TIMSTOF_MZ_GRID,
                (np.array(GRID_PARAMETERS), np.array([3])),
                (np.array([4, 2], dtype=np.uint32), np.array([2])),
                np.array([2]),
                SOURCE,
            )


class TestRows:
    """Whole chunk rows, as a signal member holds them."""

    SPECTRA = [
        Spectrum(1, 1, [100.0, 100.5, 101.0, 150.25, 151.0], [1, 2, 3, 4, 5]),
        Spectrum(2, 1, [100.5, 300.0], [6, 7]),
        Spectrum(3, 1, [99.0, 100.0, 101.5, 400.0, 400.5, 401.0, 402.0], range(8, 15)),
    ]

    @pytest.mark.parametrize("encoding", [PLAIN, DELTA, NUMPRESS_LINEAR, GRID])
    @pytest.mark.parametrize("chunk_points", [1, 2, 3, 50])
    def test_every_encoding_gives_the_points_back(self, encoding, chunk_points):
        """However the rows are cut, the points are those that went in."""
        chunk = _chunk(self.SPECTRA, encoding=encoding, chunk_points=chunk_points)

        decoded = decode_chunks(chunk, SOURCE)

        np.testing.assert_array_equal(
            decoded.spectrum_index, np.repeat([0, 1, 2], [5, 2, 7])
        )
        np.testing.assert_array_equal(
            decoded.mz, np.concatenate([s.mzs for s in self.SPECTRA])
        )
        np.testing.assert_array_equal(decoded.intensity, np.arange(1.0, 15.0))

    def test_a_member_may_mix_encodings(self):
        """Each row is decoded by what it names."""
        mixed = [DELTA, NUMPRESS_LINEAR, PLAIN, GRID]
        chunk = _chunk(self.SPECTRA, encoding=mixed, chunk_points=2)
        assert set(chunk.field("chunk_encoding").to_pylist()) == set(mixed)

        decoded = decode_chunks(chunk, SOURCE)

        np.testing.assert_array_equal(
            decoded.mz, np.concatenate([s.mzs for s in self.SPECTRA])
        )

    @pytest.mark.parametrize("encoding", [PLAIN, DELTA, NUMPRESS_LINEAR])
    @pytest.mark.parametrize("chunk_points", [1, 2, 3, 4, 5, 50])
    def test_padding_comes_back_as_nan(self, encoding, chunk_points):
        """Wherever the cut falls, a null pair stays two padding points.

        The cuts include the awkward ones: a row of padding alone, a row
        that opens with one null, and a row that opens with both.
        """
        chunk = _chunk(
            self.SPECTRA,
            null_pair_after=2,
            encoding=encoding,
            chunk_points=chunk_points,
        )

        decoded = decode_chunks(chunk, SOURCE)

        expected = np.concatenate(
            [
                np.concatenate((s.mzs[:2], [np.nan, np.nan], s.mzs[2:]))
                for s in self.SPECTRA
            ]
        )
        np.testing.assert_array_equal(decoded.mz, expected)
        np.testing.assert_array_equal(np.isnan(decoded.intensity), np.isnan(decoded.mz))

    def test_a_row_of_padding_alone_holds_no_m_z(self):
        """Both bounds are zero, and nothing is decoded from the row.

        The reference reader skips such a row whatever its list holds. Here
        the list is missing, which would otherwise be a count that
        disagrees.
        """
        chunk = _chunk(self.SPECTRA[:1], encoding=DELTA, chunk_points=5)
        row = chunk.to_pylist()[0]
        padding = dict(
            row,
            mz_chunk_start=0.0,
            mz_chunk_end=0.0,
            mz_chunk_values=None,
            intensity=[None, None],
        )

        decoded = decode_chunks(pa.array([row, padding], type=chunk.type), SOURCE)

        np.testing.assert_array_equal(
            decoded.mz, [100.0, 100.5, 101.0, 150.25, 151.0, np.nan, np.nan]
        )
        np.testing.assert_array_equal(decoded.spectrum_index, [0] * 7)

    def test_a_sliced_batch_is_read_from_its_own_rows(self):
        """A batch is a window on a larger array, offsets and all."""
        chunk = _chunk(self.SPECTRA, encoding=DELTA, chunk_points=2)

        decoded = decode_chunks(chunk.slice(3, 3), SOURCE)

        np.testing.assert_array_equal(decoded.spectrum_index, [1, 1, 2, 2, 2, 2])
        np.testing.assert_array_equal(
            decoded.mz, [100.5, 300.0, 99.0, 100.0, 101.5, 400.0]
        )


def _numpress_row(start, end, intensities=14):
    """The converter's MS-Numpress chunk as a row, under chosen bounds."""
    fields = [
        pa.field("spectrum_index", pa.uint64()),
        pa.field("mz_chunk_start", pa.float64()),
        pa.field("mz_chunk_end", pa.float64()),
        pa.field("mz_chunk_values", pa.large_list(pa.float64())),
        pa.field("chunk_encoding", pa.string()),
        pa.field("intensity", pa.large_list(pa.float32())),
        pa.field("mz_numpress_linear_bytes", pa.large_list(pa.uint8())),
    ]
    row = {
        "spectrum_index": 0,
        "mz_chunk_start": start,
        "mz_chunk_end": end,
        "mz_chunk_values": None,
        "chunk_encoding": NUMPRESS_LINEAR,
        "intensity": [float(i) for i in range(intensities)],
        "mz_numpress_linear_bytes": CONVERTER_NUMPRESS_BYTES,
    }
    return pa.array([row], type=pa.struct(fields))


class TestBounds:
    """The ends of a row are the bounds the row states.

    A lossy encoding returns the first and last m/z of a spectrum a little
    off. The archive declares the m/z range of the spectrum from the exact
    values, the resampled axis is built on that range, and a peak that
    decodes below its first bin is left out of the store.
    """

    def test_the_ends_of_a_lossy_row_are_exact(self):
        """First and last come back to the bit; the rest to the encoding."""
        expected = np.array(CONVERTER_NUMPRESS_MZ)
        buffer = np.array(CONVERTER_NUMPRESS_BYTES, dtype=np.uint8)
        as_decoded = decode_numpress_linear(
            buffer, np.array([buffer.size]), np.array([14]), SOURCE
        )
        assert as_decoded[0] != expected[0]
        assert as_decoded[-1] != expected[-1]

        decoded = decode_chunks(_numpress_row(expected[0], expected[-1]), SOURCE)

        assert decoded.mz[0] == expected[0]
        assert decoded.mz[-1] == expected[-1]
        np.testing.assert_array_equal(decoded.mz[1:-1], as_decoded[1:-1])
        np.testing.assert_allclose(decoded.mz, expected, rtol=0, atol=1e-7)

    def test_a_bound_that_is_not_the_value_is_left_alone(self):
        """A writer may give the edges of the interval it cut at."""
        buffer = np.array(CONVERTER_NUMPRESS_BYTES, dtype=np.uint8)
        as_decoded = decode_numpress_linear(
            buffer, np.array([buffer.size]), np.array([14]), SOURCE
        )

        decoded = decode_chunks(_numpress_row(100.0, 150.0), SOURCE)

        np.testing.assert_array_equal(decoded.mz, as_decoded)

    def test_padding_at_an_end_is_passed_over(self):
        """The bound belongs to the first point that is not padding."""
        mz = np.array([np.nan, 100.0000001, 101.0, 102.0000001, np.nan, np.nan, 7.0])
        counts = np.array([5, 1, 0, 1])

        take_bounds(
            mz,
            counts,
            np.array([100.0, 0.0, 0.0, 7.0]),
            np.array([102.0, 0.0, 0.0, 7.0]),
        )

        np.testing.assert_array_equal(
            mz, [np.nan, 100.0, 101.0, 102.0, np.nan, np.nan, 7.0]
        )


class TestRefusals:
    """What is not decoded is refused by its term."""

    def test_an_encoding_that_is_not_decoded(self):
        """The accession and the name are both in the message."""
        chunk = _chunk(TestRows.SPECTRA, encoding="MS:1002314")

        with pytest.raises(ConversionRefused) as refusal:
            decode_chunks(chunk, SOURCE)

        message = str(refusal.value)
        assert "MS:1002314 (MS-Numpress short logged float compression)" in message
        assert SOURCE in message

    def test_an_encoding_nobody_knows(self):
        """A term with no name here is still given as written."""
        chunk = _chunk(TestRows.SPECTRA, encoding="MS:4000000")

        with pytest.raises(ConversionRefused, match="chunk encoding MS:4000000,"):
            decode_chunks(chunk, SOURCE)

    def test_a_grid_model_that_is_not_decoded(self):
        """Grid coordinate interpolation, for one; every decoded model is listed."""
        chunk = _chunk(TestRows.SPECTRA, encoding=GRID, grid_type="MS:1003822")

        with pytest.raises(ConversionRefused) as refused:
            decode_chunks(chunk, SOURCE)

        message = str(refused.value)
        assert "grid of type MS:1003822 (grid coordinate interpolation)" in message
        assert "MS:9999002 (timsTOF m/z grid" in message

    def test_an_intensity_array_under_a_transform(self):
        """The column name is the sign; the term is named with it."""
        chunk = _chunk(TestRows.SPECTRA, encoding=DELTA)
        names = [chunk.type.field(i).name for i in range(chunk.type.num_fields)]
        renamed = pa.StructArray.from_arrays(
            [chunk.field(name) for name in names],
            names=[
                "intensity_numpress_slof_bytes" if name == "intensity" else name
                for name in names
            ],
        )

        with pytest.raises(ConversionRefused) as refusal:
            decode_chunks(renamed, SOURCE)

        message = str(refusal.value)
        assert "MS:1002314 (MS-Numpress short logged float compression)" in message
        assert "intensity_numpress_slof_bytes" in message

    def test_chunks_cut_along_another_axis(self):
        """No m/z bounds, no m/z chunks."""
        chunk = _chunk(TestRows.SPECTRA, encoding=DELTA)
        names = [chunk.type.field(i).name for i in range(chunk.type.num_fields)]
        renamed = pa.StructArray.from_arrays(
            [chunk.field(name) for name in names],
            names=[name.replace("mz_chunk", "time_chunk") for name in names],
        )

        with pytest.raises(ConversionRefused, match="lack mz_chunk_start"):
            decode_chunks(renamed, SOURCE)

    def test_an_encoding_without_its_column(self):
        """A row names MS-Numpress and the member has no byte column."""
        chunk = _chunk(TestRows.SPECTRA, encoding=DELTA)
        names = [chunk.type.field(i).name for i in range(chunk.type.num_fields)]
        relabelled = pa.StructArray.from_arrays(
            [
                (
                    pa.array([NUMPRESS_LINEAR] * len(chunk), type=pa.string())
                    if name == "chunk_encoding"
                    else chunk.field(name)
                )
                for name in names
            ],
            names=names,
        )

        with pytest.raises(ConversionRefused, match="no 'mz_numpress_linear_bytes'"):
            decode_chunks(relabelled, SOURCE)
