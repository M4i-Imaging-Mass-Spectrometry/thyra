"""The two commands that look at a converted store: ``ion-image``, ``report``.

Both take a store that already exists and write one file beside it. They
follow the metadata subcommands in everything a user can see: a missing
path or an unusable ``-o`` is a usage error (exit 2) raised before any
work, and a store the command cannot draw is one line and exit 1, never a
traceback.

Only click is imported at module scope. ``thyra.__main__`` imports this
module to dispatch on the first argument, on every run, so numpy, Zarr and
Pillow are imported inside the commands that use them.
"""

from __future__ import annotations

import os
from math import isfinite
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Tuple

import click

from ..utils.output_files import make_parent_folder, output_file_problem
from ..utils.windows_paths import prepare_zarr_read_path

if TYPE_CHECKING:
    from .store import Store

#: Kept in step with ``thyra.quicklook.ion_image.DEFAULT_PPM``; repeated
#: here so that building the command imports nothing.
_DEFAULT_PPM = 10.0


def _existing_store(path: Path) -> Tuple[Path, Path]:
    """The store's absolute path and the path to read it through.

    A path that does not exist is refused, as the metadata subcommands do.
    The check runs on the prepared path, so a store under a path past the
    Windows limit is found (see ``thyra.metadata.schema.cli``). Preparing
    a path walks the store, so the prepared path is handed on and the walk
    happens once.
    """
    read_path = prepare_zarr_read_path(path)
    if not read_path.exists():
        raise click.BadParameter(
            f"Path {str(path)!r} does not exist.", param_hint="STORE"
        )
    return Path(os.path.abspath(path)), read_path


def _check_output(output: Path, store: Path) -> None:
    """Refuse an output file that cannot be written, or that is in the store.

    A store is one dataset and stays whole: nothing is written into it.
    """
    problem = output_file_problem(output)
    if problem is None and Path(os.path.abspath(output)).is_relative_to(store):
        problem = (
            f"Cannot write {output}: it is inside the store. Give a path "
            "outside the .zarr folder."
        )
    if problem is not None:
        raise click.BadParameter(problem, param_hint="'-o' / '--output'")


def _open(store: Path, read_path: Path) -> "Store":
    """Open the store, turning a refusal into a one-line click error."""
    from .store import QuicklookError, open_store

    try:
        return open_store(store, read_path)
    except QuicklookError as exc:
        raise click.ClickException(str(exc)) from exc


def _write(target: Path, content: bytes) -> None:
    """Write the output file; a failure is one line, not a traceback."""
    try:
        make_parent_folder(target)
        target.write_bytes(content)
    except OSError as exc:
        reason = exc.strerror or type(exc).__name__
        raise click.ClickException(f"Could not write {target}: {reason}.") from exc


def _positive(value: float, hint: str) -> None:
    if not isfinite(value) or value <= 0:
        raise click.BadParameter("must be a number above zero.", param_hint=hint)


def _blank_image_note(peak: Optional[float]) -> str:
    """What to say about an image in which no pixel has signal."""
    note = "No pixel has signal in this window: the image is blank."
    if peak is not None:
        note += f" The tallest peak of the mean spectrum nearby is at m/z {peak:.4f}."
    return note


@click.command("ion-image")
@click.argument("store", type=click.Path(exists=False, path_type=Path))
@click.argument("mz", type=click.FLOAT)
@click.option(
    "--ppm",
    type=click.FLOAT,
    default=_DEFAULT_PPM,
    show_default=True,
    help="Half-width of the m/z window, in parts per million of MZ.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Output file (default: <store>_mz<MZ>.png next to the store).",
)
def ion_image_command(
    store: Path, mz: float, ppm: float, output_path: Optional[Path]
) -> None:
    """Write the ion image of one m/z value as a PNG.

    STORE is a converted .zarr store and MZ the m/z value. The image adds
    up, per pixel, every m/z bin within --ppm of MZ. A position that was
    never measured is drawn grey.

    A store with several tables of spectra, one per slice, gives one image
    per table; each file name then ends in its table's name.
    """
    from .ion_image import make_ion_image, tallest_peak_near
    from .render import format_intensity, ion_image_figure, png_bytes, scale_top
    from .store import QuicklookError, store_name

    _positive(mz, "MZ")
    _positive(ppm, "'--ppm'")
    store, read_path = _existing_store(store)
    if output_path is None:
        output_path = store.with_name(f"{store_name(store)}_mz{mz:.4f}.png")
    _check_output(output_path, store)

    opened = _open(store, read_path)
    try:
        tables = opened.msi_tables()
        targets = {summary.name: output_path for summary in tables}
        if len(tables) > 1:
            for name in targets:
                targets[name] = output_path.with_name(
                    f"{output_path.stem}_{name}{output_path.suffix}"
                )
                _check_output(targets[name], store)
        for name, target in targets.items():
            table = opened.open_table(name)
            ion = make_ion_image(table, mz, ppm)
            figure = ion_image_figure(ion, f"{opened.name}, {name}")
            _write(target, png_bytes(figure))
            window = ion.window
            click.echo(
                f"{name}: m/z {mz:.4f} +/- {ppm:g} ppm is "
                f"{window.n_bins} bin(s), {window.low:.4f} to {window.high:.4f}. "
                f"Intensity {format_intensity(ion.lowest)} to "
                f"{format_intensity(ion.highest)} over {ion.n_pixels} pixels; "
                f"the colours run from 0 to {format_intensity(scale_top(ion.image))}."
            )
            if ion.highest == 0:
                click.echo(_blank_image_note(tallest_peak_near(table, mz)))
            click.echo(f"Wrote {target} ({figure.width} x {figure.height} pixels)")
    except QuicklookError as exc:
        raise click.ClickException(str(exc)) from exc


@click.command("report")
@click.argument("store", type=click.Path(exists=False, path_type=Path))
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Output file (default: <store>_report.html next to the store).",
)
def report_command(store: Path, output_path: Optional[Path]) -> None:
    """Write a one-file HTML overview of a converted store.

    STORE is a converted .zarr store. The page shows what the store says
    about the acquisition, the TIC image, the mean spectrum, ion images of
    its tallest peaks, the optical image when there is one, and the list of
    tables. It opens in any browser, without a network connection.
    """
    from .report import build_report
    from .store import QuicklookError, store_name

    store, read_path = _existing_store(store)
    if output_path is None:
        output_path = store.with_name(f"{store_name(store)}_report.html")
    _check_output(output_path, store)

    opened = _open(store, read_path)
    try:
        page = build_report(opened).encode("utf-8")
    except QuicklookError as exc:
        raise click.ClickException(str(exc)) from exc
    _write(output_path, page)
    click.echo(f"Wrote {output_path} ({len(page) / 1024:.0f} kB)")


# Consumed by thyra.__main__ beside the metadata subcommands.
QUICKLOOK_SUBCOMMANDS: Dict[str, click.Command] = {
    "ion-image": ion_image_command,
    "report": report_command,
}
