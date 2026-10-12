"""What a parity run is held to, as named constants.

These start permissive on purpose: the first job of the harness is to
measure how far Thyra sits from SCiLS Lab, not to decide how close it must
be. Each report prints the measured value beside the threshold it was held
to, so tightening one is a matter of reading the reports and moving the
number. A failure here means "outside the measured envelope", which may as
well be a SCiLS quirk as a Thyra defect; the report is what tells them
apart.
"""

#: Largest axis offset, in ppm, for axes that should be the same grid.
AXIS_MAX_PPM = 50.0

#: Largest axis offset as a fraction of the local bin width. Above 0.5 a
#: Thyra bin no longer covers the SCiLS bin of the same index.
AXIS_MAX_BIN_FRACTION = 0.5

#: Largest difference between the two fitted spacing exponents
#: (:func:`tests.parity.metrics.spacing_exponent`). The laws sit 0.5 apart.
AXIS_MAX_EXPONENT_DIFF = 0.1

#: Mean spectrum, compared bin by bin when the axes match.
MEAN_MAX_REL_L1 = 0.5
MEAN_MIN_PEARSON = 0.9
MEAN_SUM_RATIO_RANGE = (0.8, 1.25)

#: Mean spectrum after both sides are summed into
#: :data:`COARSE_BINS` equal-width bins: insensitive to where each fine
#: bin sits, so it isolates the intensity scale from the axis.
COARSE_BINS = 1000
COARSE_MEAN_MAX_REL_L1 = 0.2
COARSE_MEAN_MIN_PEARSON = 0.98

#: Per-pixel TIC.
TIC_MIN_PEARSON = 0.95
TIC_MEDIAN_RATIO_RANGE = (0.8, 1.25)

#: Below this coefficient of variation of SCiLS's per-spot TIC, the TIC is
#: constant across spots (a TIC-normalised export) and its correlation with
#: Thyra's measures rounding noise only; the ratio is still checked.
TIC_MIN_REFERENCE_CV = 1e-3

#: Fraction of SCiLS spots that must pair with a Thyra pixel.
MIN_MATCHED_FRACTION = 0.95

#: Individual spectra, compared bin by bin when the axes match.
SPOT_MIN_PEARSON = 0.8
SPOT_COARSE_MIN_PEARSON = 0.95
