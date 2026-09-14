# tests/unit/converters/test_phi_default_axis_resolves_peaks.py
"""The PHI default axis has to put enough bins across a peak to keep its shape.

``PhiToFSIMSDetector`` used to inherit the generic ``linear_tof`` default of
17 mDa at m/z 300. On a real nanoTOF acquisition that laid **one bin across
the base peak**: 42% of the peaks measured over a twelve-acquisition corpus
came out under two bins per peak width, and the acquisition's own vendor peak
list -- which carries a 6 mDa window at nominal m/z 27 to separate C15N- from
13CN-, 6.3 mDa apart -- had a single bin inside that window. Per-peak recovery
against the instrument's own ``.bif6`` peak-image export ran 88.6-113.1%.

Nothing downstream could see it. Nearest-neighbour conserves counts exactly,
so the total ion image stayed bit-identical to the vendor's own export on both
the broken and the fixed axis; only *windowed* numbers moved. That is why the
third test here pins the TIC as well: a future change must not be able to buy
peak fidelity back by trading away sum fidelity, or the reverse.

The fix was the axis law, not the bin width. A PHI pixel is a list of ion
arrival times, so re-binning is exact and the only question is how fine a
histogram preserves the peak -- which is what ``AxisType.TOF`` and
``PHI_NANOTOF_LAW`` answer, at three bins per measured peak width everywhere.
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
import pytest

from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)
from thyra.readers.phi import PhiReader
from thyra.resampling.mass_axis.tof_generator import (
    DEFAULT_BINS_PER_FWHM,
    PHI_NANOTOF_LAW,
    tof_fwhm_mda,
)
from thyra.resampling.types import ResamplingConfig

# The reference acquisition's own calibration and channel width.
SLOPE = 0.382653
BIN_SIZE_NS = 0.128
START_MZ = 0.5
STOP_MZ = 200.0
PIXELS = 8
HEADER_SIZE = 2048

#: Peaks planted in the synthetic acquisition, as ``(m/z, width multiple)``.
#: The multiple is of :data:`PHI_NANOTOF_LAW`'s own prediction at that m/z.
#: 1.3x is the median of the real corpus -- the law is fitted to its 5th
#: percentile, so a typical peak is wider than it predicts.
PLANTED = [(26.0, 1.3), (43.0, 1.3), (79.0, 1.3), (147.0, 1.3)]

EVENTS_PER_PEAK_PER_PIXEL = 400

_HEADER = """SOFH\r
Platform: PC\r
SoftwareVersion: SmartSoft-TOF V3.2.0.35\r
Technique: -TofSIMS\r
SecIonPolarity: -\r
DataType: unsigned long\r
SpecBinSize: {bin_size_ns}\r
SpecBinIncr: 0.00100\r
StartFlightTime: {start_us}\r
StopFlightTime: {stop_us}\r
Mass/Time: {slope}\r
MassOffset: 0.0000000\r
ImagePixels: {pixels}\r
ScanWidthX: 128.00000\r
NoFrames: 1\r
Calibrated: yes\r
HeaderSize: {header_size}\r
[Acq Base]\r
Start Mass (amu)={start_mz}\r
End Mass (amu)={stop_mz}\r
Frames=1\r
MSMS Active=Inactive\r
Raster Pattern=Scatter\r
Raster Resolution={pixels}\r
[PHI LMIG]\r
Gun Particle=Bi3 +\r
Raster Size (um)=8.0\r
Raster Size Calibration=1.550, 1.550\r
[Mosaic Area]\r
Number Of Tiles X=1\r
Number Of Tiles Y=1\r
[Polarity]\r
Polarity=Negative (-) Ions\r
[Main]\r
Project Name=SmartSoft-TOF\r
EOFH\r
"""


def _tof_us(mz: float) -> float:
    """Flight time of an m/z under the synthetic calibration (offset 0)."""
    return float(np.sqrt(mz) / SLOPE)


def _sigma_us(mz: float, fwhm_da: float) -> float:
    """Flight-time sigma of a peak of ``fwhm_da`` at ``mz``.

    ``m = (slope t)^2`` gives ``dm/dt = 2 slope sqrt(m)``, so a width in m/z
    divides by that to become a width in time.
    """
    fwhm_us = fwhm_da / (2.0 * SLOPE * np.sqrt(mz))
    return float(fwhm_us / (2.0 * np.sqrt(2.0 * np.log(2.0))))


def planted_fwhm_da(mz: float, multiple: float) -> float:
    """The FWHM a planted peak is drawn with, in Da."""
    return multiple * float(tof_fwhm_mda(mz, *PHI_NANOTOF_LAW)) * 1e-3


def _block(block_id: int, payload: bytes = b"") -> bytes:
    return (
        struct.pack("<H", block_id)
        + b"\x00" * 10
        + struct.pack("<I", len(payload))
        + payload
    )


def _event(x: int, y: int, tof_ps: int) -> int:
    return (1 << 31) | (int(tof_ps) & 0x7FFFFFFF) | (x << 32) | (y << 43)


@pytest.fixture(scope="module")
def synthetic_raw(tmp_path_factory) -> Path:
    """A PHI ``.raw`` whose peaks have a width we chose.

    Every pixel gets the same four peaks, so the store's own mean spectrum
    can be measured against the width that was planted.
    """
    rng = np.random.default_rng(20260914)
    start_us = _tof_us(START_MZ)
    stop_us = _tof_us(STOP_MZ)

    words = []
    for y in range(PIXELS):
        for x in range(PIXELS):
            for mz, multiple in PLANTED:
                times = rng.normal(
                    _tof_us(mz),
                    _sigma_us(mz, planted_fwhm_da(mz, multiple)),
                    EVENTS_PER_PEAK_PER_PIXEL,
                )
                for t in times:
                    words.append(_event(x, y, round(t * 1e6)))

    payload = _block(1, np.array(words, dtype="<u8").tobytes()) + _block(2)
    header = _HEADER.format(
        bin_size_ns=f"{BIN_SIZE_NS:.5f}",
        start_us=f"{start_us:.5f}",
        stop_us=f"{stop_us:.5f}",
        slope=f"{SLOPE:.7f}",
        pixels=PIXELS,
        header_size=HEADER_SIZE,
        start_mz=START_MZ,
        stop_mz=STOP_MZ,
    ).encode("latin-1")
    assert len(header) <= HEADER_SIZE, "synthetic header exceeds HeaderSize"

    path = tmp_path_factory.mktemp("phi") / "synthetic.raw"
    path.write_bytes(header.ljust(HEADER_SIZE, b"\x00") + payload + _block(0))
    return path


def _convert(raw: Path, out: Path, resample: bool):
    """Convert the way the CLI does, with resampling on or off."""
    converter = StreamingSpatialDataConverter(
        reader=PhiReader(raw),
        output_path=out,
        dataset_id="phi",
        include_optical=False,
        resampling_config=ResamplingConfig() if resample else None,
    )
    assert converter.convert() is True
    return out


def _table(store: Path):
    import spatialdata as sd

    return next(iter(sd.read_zarr(store).tables.values()))


@pytest.fixture(scope="module")
def resampled(synthetic_raw, tmp_path_factory):
    """The default route: the detector picks the axis."""
    out = tmp_path_factory.mktemp("r") / "d.zarr"
    return _table(_convert(synthetic_raw, out, resample=True))


@pytest.fixture(scope="module")
def native(synthetic_raw, tmp_path_factory):
    """``--no-resample``: the detector's own channel grid, untouched."""
    out = tmp_path_factory.mktemp("n") / "n.zarr"
    return _table(_convert(synthetic_raw, out, resample=False))


def _mean_spectrum(table) -> np.ndarray:
    return np.asarray(table.X.mean(axis=0)).ravel()


def _measured_sigma_da(
    mz: np.ndarray, spectrum: np.ndarray, centre: float, half_span: float
) -> float:
    """Intensity-weighted RMS width of the peak at ``centre``, in Da.

    The second moment, not half-max crossings. Crossings need the flanks
    resolved finely enough to interpolate between two samples, and the
    whole point of this axis is that it carries only three bins across a
    peak width -- measuring that way scatters by +-15% depending on where
    the bin grid happens to fall. Every count contributes to a moment, so
    it stays accurate at this bin density and the binning's own effect on
    it is exactly predictable.
    """
    lo, hi = np.searchsorted(mz, [centre - half_span, centre + half_span])
    x = mz[lo:hi]
    w = spectrum[lo:hi].astype(np.float64)
    total = w.sum()
    assert total > 0, f"no counts around m/z {centre}"
    mean = float((w * x).sum() / total)
    return float(np.sqrt((w * (x - mean) ** 2).sum() / total))


class TestBinsAcrossAPeak:
    """The axis must resolve the peaks the instrument produces."""

    def test_at_least_three_bins_across_every_planted_peak(self, resampled):
        """The contract: three bins per peak width, at every m/z.

        The old ``linear_tof`` default managed 1.0 at the base peak of the
        reference acquisition, which is what made windowed sums a coin flip.
        """
        mz = resampled.var["mz"].to_numpy()
        for centre, multiple in PLANTED:
            width = planted_fwhm_da(centre, multiple)
            lo, hi = np.searchsorted(mz, [centre - width / 2, centre + width / 2])
            assert hi - lo >= 3, (
                f"m/z {centre}: {hi - lo} bins across a {width * 1e3:.2f} mDa "
                f"FWHM, need at least 3"
            )

    def test_bin_width_is_the_law_over_three(self, resampled):
        """The axis is laid at ``FWHM(m) / 3`` under the declared pair."""
        mz = resampled.var["mz"].to_numpy()
        for centre, _ in PLANTED:
            i = int(np.searchsorted(mz, centre))
            realised = float(mz[i + 1] - mz[i])
            expected = (
                float(tof_fwhm_mda(centre, *PHI_NANOTOF_LAW))
                * 1e-3
                / DEFAULT_BINS_PER_FWHM
            )
            assert realised == pytest.approx(expected, rel=0.02)

    def test_the_stored_peak_keeps_its_planted_width(self, resampled, native):
        """Binning at three per width broadens a peak by 1.5%, and no more.

        Histogramming convolves the peak with one bin, and variances add:
        ``sigma_stored^2 = sigma^2 + w^2 / 12``. At ``w = FWHM / 3`` that is
        a 1.5% broadening, which is the whole price of this axis. The old
        ``linear_tof`` default put ``w`` at roughly one FWHM on this data,
        where the same arithmetic gives 14% -- and 14% is the *average*
        case, with the peak's position against the grid deciding the rest.
        """
        mz_r = resampled.var["mz"].to_numpy()
        mz_n = native.var["mz"].to_numpy()
        stored = _mean_spectrum(resampled)
        source = _mean_spectrum(native)

        for centre, multiple in PLANTED:
            span = 6.0 * planted_fwhm_da(centre, multiple)
            sigma_source = _measured_sigma_da(mz_n, source, centre, span)
            sigma_stored = _measured_sigma_da(mz_r, stored, centre, span)
            i = int(np.searchsorted(mz_r, centre))
            bin_width = float(mz_r[i + 1] - mz_r[i])
            predicted = np.sqrt(sigma_source**2 + bin_width**2 / 12.0)
            assert sigma_stored == pytest.approx(predicted, rel=0.03), (
                f"m/z {centre}: source sigma {sigma_source * 1e3:.3f} mDa, "
                f"stored {sigma_stored * 1e3:.3f}, one bin predicts "
                f"{predicted * 1e3:.3f}"
            )
            assert sigma_stored < 1.05 * sigma_source, (
                f"m/z {centre}: the axis broadened the peak by "
                f"{100 * (sigma_stored / sigma_source - 1):.1f}%"
            )


class TestWindowedRecovery:
    """Sums over a narrow window must survive the resampling.

    This is the check that fails loudly on the old axis: with one bin inside
    a peak window, whether an ion over- or under-collects depends on where
    its peak happens to fall against the bin grid, so it is a coin flip per
    ion -- 88.6% for one and 113.1% for another on the same real file.
    """

    def test_window_sums_match_the_unresampled_store(self, resampled, native):
        """Each planted peak, summed over +-1 FWHM, against the raw channels."""
        mz_r = resampled.var["mz"].to_numpy()
        mz_n = native.var["mz"].to_numpy()
        totals_r = np.asarray(resampled.X.sum(axis=0)).ravel()
        totals_n = np.asarray(native.X.sum(axis=0)).ravel()

        for centre, multiple in PLANTED:
            width = planted_fwhm_da(centre, multiple)
            lo, hi = centre - width, centre + width
            ours = totals_r[slice(*np.searchsorted(mz_r, [lo, hi]))].sum()
            theirs = totals_n[slice(*np.searchsorted(mz_n, [lo, hi]))].sum()
            assert theirs > 0
            assert ours / theirs == pytest.approx(1.0, abs=0.05), (
                f"m/z {centre}: recovered {100 * ours / theirs:.1f}% of the "
                f"counts the unresampled store holds in the same window"
            )


class TestTotalsAreUntouched:
    """Neither route may trade sum fidelity for peak fidelity."""

    def test_every_pixel_total_is_identical_on_both_routes(self, resampled, native):
        """Nearest-neighbour re-histograms counts; it never creates or drops one.

        On the real acquisition this holds to the count against the vendor's
        own total-ion export, on 262,144 of 262,144 pixels, on every axis
        tried -- which is exactly why the windowed numbers above are the only
        thing that can catch a mistuned axis.
        """
        per_pixel_r = np.asarray(resampled.X.sum(axis=1)).ravel()
        per_pixel_n = np.asarray(native.X.sum(axis=1)).ravel()
        assert per_pixel_r.sum() == pytest.approx(per_pixel_n.sum(), rel=0)
        np.testing.assert_allclose(
            np.sort(per_pixel_r), np.sort(per_pixel_n), rtol=0, atol=0
        )

    def test_the_total_is_every_event_that_was_planted(self, resampled):
        """No silent loss between the reader and the store."""
        planted = PIXELS * PIXELS * len(PLANTED) * EVENTS_PER_PEAK_PER_PIXEL
        stored = float(np.asarray(resampled.X.sum()))
        # A few events per peak scatter outside the declared mass range.
        assert 0.99 * planted <= stored <= planted
