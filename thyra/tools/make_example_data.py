# thyra/tools/make_example_data.py
"""Generate a small synthetic MSI dataset for tutorials and smoke tests.

The generated imzML/ibd pair is deterministic for a given seed and small
enough to convert in seconds, so it can be used to verify an installation
end to end without downloading vendor data.

The phantom is a sagittal brain-like section with two internal structures
that carry distinct peak sets, so TIC images and single-ion images show
real spatial contrast.
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter

logger = logging.getLogger(__name__)

# Peak sets in m/z, chosen in the phospholipid range covered by the
# example dataset on Zenodo (250-1200 Da). Values are illustrative --
# this is synthetic data and is not intended to model real lipid ratios.
_MATRIX_PEAKS: Dict[float, float] = {
    # Low-mass matrix cluster ions, present everywhere including background.
    296.4: 0.30,
    379.1: 0.22,
    551.2: 0.15,
}
_SHARED_TISSUE_PEAKS: Dict[float, float] = {
    # Abundant in any tissue region.
    760.6: 1.00,
    782.6: 0.85,
    786.6: 0.55,
}
_GREY_MATTER_PEAKS: Dict[float, float] = {
    772.5: 0.90,
    798.5: 0.70,
    826.6: 0.45,
}
_WHITE_MATTER_PEAKS: Dict[float, float] = {
    888.6: 1.00,
    906.6: 0.65,
    834.5: 0.40,
}

_PEAK_SIGMA_DA = 0.5


def _ellipse_mask(
    shape: Tuple[int, int],
    centre: Tuple[float, float],
    radii: Tuple[float, float],
) -> np.ndarray:
    """Boolean mask of an axis-aligned ellipse on a (y, x) grid."""
    n_y, n_x = shape
    # np.mgrid is typed as Any, which propagates through the arithmetic below
    # and makes the returned mask untyped. Broadcasting two typed 1-D ranges
    # produces the same grid while keeping the return type checkable.
    yy = np.arange(n_y, dtype=np.float64).reshape(n_y, 1)
    xx = np.arange(n_x, dtype=np.float64).reshape(1, n_x)
    cy, cx = centre
    ry, rx = radii
    return ((yy - cy) / ry) ** 2 + ((xx - cx) / rx) ** 2 <= 1.0


def _build_region_masks(n_x: int, n_y: int) -> Dict[str, np.ndarray]:
    """Build tissue, grey matter and white matter masks."""
    shape = (n_y, n_x)

    tissue = _ellipse_mask(
        shape, centre=(n_y / 2.0, n_x / 2.0), radii=(n_y * 0.42, n_x * 0.45)
    )

    # White matter: a smaller ellipse offset towards the lower left,
    # loosely standing in for a corpus-callosum-like structure.
    white = _ellipse_mask(
        shape,
        centre=(n_y * 0.60, n_x * 0.42),
        radii=(n_y * 0.16, n_x * 0.26),
    )
    white &= tissue

    # Grey matter is the rest of the tissue.
    grey = tissue & ~white

    return {"tissue": tissue, "grey": grey, "white": white}


def _add_peaks(
    spectrum: np.ndarray,
    mz_axis: np.ndarray,
    peaks: Dict[float, float],
    scale: float,
) -> None:
    """Add Gaussian peaks to a spectrum in place."""
    for centre, amplitude in peaks.items():
        # Only touch the window around the peak -- much faster than
        # evaluating the Gaussian over the whole axis per peak.
        lo, hi = np.searchsorted(
            mz_axis, [centre - 6 * _PEAK_SIGMA_DA, centre + 6 * _PEAK_SIGMA_DA]
        )
        if lo >= hi:
            continue
        window = mz_axis[lo:hi]
        spectrum[lo:hi] += (amplitude * scale) * np.exp(
            -((window - centre) ** 2) / (2.0 * _PEAK_SIGMA_DA**2)
        )


def generate_example_imzml(
    output_path: Path,
    n_x: int = 48,
    n_y: int = 36,
    pixel_size_um: float = 25.0,
    mz_min: float = 250.0,
    mz_max: float = 1200.0,
    n_mz_bins: int = 4000,
    seed: int = 0,
) -> Path:
    """Write a synthetic continuous-mode imzML dataset.

    Args:
        output_path: Path of the .imzML file to write. The matching .ibd
            file is written alongside it by the writer.
        n_x: Grid width in pixels.
        n_y: Grid height in pixels.
        pixel_size_um: Pixel size written into the imzML metadata, so
            that Thyra can auto-detect it during conversion.
        mz_min: Lower bound of the m/z axis.
        mz_max: Upper bound of the m/z axis.
        n_mz_bins: Number of points on the shared m/z axis.
        seed: Seed for the noise generator, for reproducibility.

    Returns:
        The path to the written .imzML file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    mz_axis = np.linspace(mz_min, mz_max, n_mz_bins)
    masks = _build_region_masks(n_x, n_y)

    # A smooth intensity gradient across the section, so that pixels
    # within one region are not all identical.
    yy, xx = np.mgrid[0:n_y, 0:n_x]
    gradient = 0.75 + 0.5 * (xx / max(n_x - 1, 1)) + 0.25 * (yy / max(n_y - 1, 1))

    with ImzMLWriter(
        str(output_path), mode="continuous", spec_type="profile"
    ) as writer:
        for y in range(n_y):
            for x in range(n_x):
                spectrum = np.zeros(n_mz_bins, dtype=np.float64)

                # Matrix signal is present over the whole slide.
                _add_peaks(spectrum, mz_axis, _MATRIX_PEAKS, scale=120.0)

                if masks["tissue"][y, x]:
                    scale = 4000.0 * gradient[y, x]
                    _add_peaks(spectrum, mz_axis, _SHARED_TISSUE_PEAKS, scale)
                    if masks["white"][y, x]:
                        _add_peaks(spectrum, mz_axis, _WHITE_MATTER_PEAKS, scale)
                    else:
                        _add_peaks(spectrum, mz_axis, _GREY_MATTER_PEAKS, scale)

                # Baseline plus shot-like noise.
                spectrum += rng.normal(loc=4.0, scale=1.5, size=n_mz_bins)
                np.clip(spectrum, 0.0, None, out=spectrum)

                # imzML coordinates are 1-based.
                writer.addSpectrum(
                    mz_axis, spectrum.astype(np.float32), (x + 1, y + 1, 1)
                )

    _inject_pixel_size(output_path, pixel_size_um)
    return output_path


def _inject_pixel_size(imzml_path: Path, pixel_size_um: float) -> None:
    """Add pixel size cvParams to the scanSettings element.

    pyimzml's writer does not emit IMS:1000046 / IMS:1000047, which is
    what Thyra reads to auto-detect pixel size. Adding them here keeps
    the generated file representative of real acquisitions.
    """
    text = imzml_path.read_text(encoding="ISO-8859-1")

    anchor = '<cvParam cvRef="IMS" accession="IMS:1000043"'
    idx = text.find(anchor)
    if idx == -1:
        logger.warning(
            "Could not locate scanSettings anchor in %s; pixel size not written. "
            "Pass --pixel-size to thyra when converting.",
            imzml_path,
        )
        return

    line_end = text.find("\n", idx)
    if line_end == -1:
        line_end = len(text)

    # Note the asymmetric names: the imaging MS ontology defines IMS:1000046
    # as "pixel size (x)" but IMS:1000047 as "pixel size y".
    indent = " " * 6
    addition = (
        f'\n{indent}<cvParam cvRef="IMS" accession="IMS:1000046" '
        f'name="pixel size (x)" value="{pixel_size_um}"/>'
        f'\n{indent}<cvParam cvRef="IMS" accession="IMS:1000047" '
        f'name="pixel size y" value="{pixel_size_um}"/>'
    )

    imzml_path.write_text(
        text[:line_end] + addition + text[line_end:], encoding="ISO-8859-1"
    )


def _create_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate a small synthetic MSI dataset (imzML + ibd) for "
            "tutorials and installation checks."
        )
    )
    parser.add_argument(
        "output",
        nargs="?",
        default="example_data/synthetic_brain.imzML",
        help=("Output .imzML path " "(default: example_data/synthetic_brain.imzML)"),
    )
    parser.add_argument(
        "--pixels-x", type=int, default=48, help="Grid width in pixels (default: 48)"
    )
    parser.add_argument(
        "--pixels-y", type=int, default=36, help="Grid height in pixels (default: 36)"
    )
    parser.add_argument(
        "--pixel-size",
        type=float,
        default=25.0,
        help="Pixel size in micrometers written to metadata (default: 25)",
    )
    parser.add_argument(
        "--mz-bins",
        type=int,
        default=4000,
        help="Number of points on the m/z axis (default: 4000)",
    )
    parser.add_argument("--seed", type=int, default=0, help="Random seed (default: 0)")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    return parser


def main() -> int:
    """Entry point for the thyra-example-data command."""
    args = _create_parser().parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(message)s",
    )

    output = Path(args.output)
    if output.suffix.lower() != ".imzml":
        output = output.with_suffix(".imzML")

    path = generate_example_imzml(
        output,
        n_x=args.pixels_x,
        n_y=args.pixels_y,
        pixel_size_um=args.pixel_size,
        n_mz_bins=args.mz_bins,
        seed=args.seed,
    )

    ibd = path.with_suffix(".ibd")
    total_mb = (path.stat().st_size + ibd.stat().st_size) / 1024**2
    n_pixels = args.pixels_x * args.pixels_y

    print(f"Wrote {path}")
    print(f"Wrote {ibd}")
    print(
        f"{n_pixels} pixels ({args.pixels_x} x {args.pixels_y}), "
        f"{args.mz_bins} m/z bins, {args.pixel_size} um, {total_mb:.1f} MB total"
    )
    print()
    print("Convert it with:")
    print(f"  thyra {path} {path.with_suffix('.zarr')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
