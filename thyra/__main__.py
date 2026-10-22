"""The command line, which is argument checking around one call.

Every command here ends in :func:`thyra.convert.convert_msi` or one of its
neighbours; the CLI adds no conversion behaviour of its own, so that the
Python API and the command line cannot describe different products. What
it does add is the part a library cannot: refusing bad input while the
user is still looking at what they typed.

Hence the number of ``_validate_*`` helpers below. A pixel size, an m/z
range, a mobility grid and a ToF law can each be individually well formed
and jointly meaningless, and a conversion is long enough that discovering
this at write time is expensive. The output path is checked for
writability before anything is read, for the same reason.

The import block is ordered deliberately and carries ``# noqa: E402``: the
dependency warning configuration has to run before the libraries that emit
those warnings are imported, so the usual "imports first" rule is
inverted on purpose here.

**What it does not import at module scope is ``thyra.convert``** (issue
#381). Binding ``convert_msi`` here meant every invocation paid for the
converter, ``pandas`` and ``spatialdata`` before click had parsed a single
argument -- including ``thyra --version`` and the two subcommands,
``validate`` and ``metadata``, that never reach a converter at all. The
conversion command imports it where it calls it. Tests that replace the
conversion patch ``thyra.convert.convert_msi``, which is the function
rather than a second name for it.
"""

# Configure dependencies to suppress warnings BEFORE any imports
import logging  # noqa: E402
import os  # noqa: E402
import warnings  # noqa: E402
from math import isfinite  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Dict, List, Literal, Optional, Tuple  # noqa: E402

import click  # noqa: E402

from thyra import __version__  # noqa: E402
from thyra.reader_options import (  # noqa: E402
    READER_OPTION_FORMATS,
    format_name,
    reader_option_reaches,
)
from thyra.resampling.mobility_grid import MOBILITY_CHANNELS  # noqa: E402
from thyra.utils.logging_config import setup_logging  # noqa: E402
from thyra.utils.output_files import (  # noqa: E402
    make_parent_folder,
    output_file_problem,
)

# Not ``__name__``: under ``python -m thyra`` that is "__main__", and
# setup_logging configures the "thyra" logger with propagate=False, so
# every line the CLI itself logged went to Python's last-resort handler
# instead -- absent from --log-file entirely, and "Conversion completed
# successfully" appeared in 0 of 67 successful runs (issue #258). The
# console-script entry point was unaffected, which is why it survived.
logger = logging.getLogger("thyra.cli")

# Configure Dask to use new query planning (silences legacy DataFrame warning)
os.environ["DASK_DATAFRAME__QUERY_PLANNING"] = "True"

# Suppress dependency warnings at the earliest possible moment
warnings.filterwarnings("ignore", category=FutureWarning, module="dask")
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
warnings.filterwarnings(
    "ignore",
    message="The legacy Dask DataFrame implementation is deprecated",
    category=FutureWarning,
)


def _is_usable_number(value: float) -> bool:
    """Whether a numeric option carries a real, positive quantity.

    ``value <= 0`` is False for NaN and for +infinity, so every guard
    written that way admitted both. ``--pixel-size nan`` reached the store
    as the dataset's pixel size, and ``--tof-law nan 1`` was refused only
    by the axis generator, as a traceback, after a 47 s metadata scan
    (issue #231).
    """
    return isfinite(value) and value > 0


def _validate_basic_params(pixel_size: Optional[float], dataset_id: str) -> None:
    """Validate basic conversion parameters."""
    if pixel_size is not None and not _is_usable_number(pixel_size):
        raise click.BadParameter(
            "Pixel size must be a finite positive number", param_hint="pixel_size"
        )
    from thyra.convert import dataset_id_problem

    # Emptiness used to be the whole check, so an id with a space or a
    # slash was refused by SpatialData only at the write, after both
    # passes over the source (issue #250).
    problem = dataset_id_problem(dataset_id)
    if problem is not None:
        raise click.BadParameter(problem, param_hint="dataset_id")


def _validate_positive_int(value: Optional[int], param_name: str, label: str) -> None:
    """Validate that an optional int parameter is positive if provided."""
    if value is not None and value <= 0:
        raise click.BadParameter(f"{label} must be positive", param_hint=param_name)


def _validate_positive_float(
    value: Optional[float], param_name: str, label: str
) -> None:
    """Validate that an optional float parameter is positive if provided."""
    if value is not None and not _is_usable_number(value):
        raise click.BadParameter(
            f"{label} must be a finite positive number", param_hint=param_name
        )


def _validate_mz_range(min_mz: Optional[float], max_mz: Optional[float]) -> None:
    """Validate that min_mz is less than max_mz when both are provided."""
    if min_mz is not None and max_mz is not None and min_mz >= max_mz:
        raise click.BadParameter("Minimum m/z must be less than maximum m/z")


def _validate_resampling_params(
    resample_bins: Optional[int],
    resample_min_mz: Optional[float],
    resample_max_mz: Optional[float],
    resample_width_at_mz: Optional[float],
    resample_reference_mz: float,
) -> None:
    """Validate resampling parameters."""
    _validate_positive_int(resample_bins, "resample_bins", "Number of resampling bins")
    _validate_positive_float(resample_min_mz, "resample_min_mz", "Minimum m/z")
    _validate_positive_float(resample_max_mz, "resample_max_mz", "Maximum m/z")
    _validate_mz_range(resample_min_mz, resample_max_mz)

    if resample_bins is not None and resample_width_at_mz is not None:
        raise click.BadParameter(
            "--resample-bins and --resample-width-at-mz are mutually exclusive"
        )

    _validate_positive_float(resample_width_at_mz, "resample_width_at_mz", "Mass width")

    if not _is_usable_number(resample_reference_mz):
        raise click.BadParameter(
            "Reference m/z must be a finite positive number",
            param_hint="resample_reference_mz",
        )


def _validate_mobility_grid_params(
    mobility_bins: Optional[int],
    mobility_min: Optional[float],
    mobility_max: Optional[float],
) -> None:
    """Validate the mobility grid's own numbers.

    Unchecked, ``--mobility-bins 0`` reached the grid builder and
    ``--mobility-min 2 --mobility-max 1`` an empty range, both on a
    source where the flags do apply (issue #260). The edges are
    positions on the mobility axis rather than quantities, so unlike
    every other float option here they may be zero or negative; only
    their order and finiteness are checked.
    """
    if mobility_bins is not None and mobility_bins <= 0:
        raise click.BadParameter(
            "A mobility grid needs at least one channel", param_hint="mobility_bins"
        )
    for value, hint in ((mobility_min, "mobility_min"), (mobility_max, "mobility_max")):
        if value is not None and not isfinite(value):
            raise click.BadParameter(
                "Mobility grid edges must be finite numbers", param_hint=hint
            )
    if (
        mobility_min is not None
        and mobility_max is not None
        and mobility_min >= mobility_max
    ):
        raise click.BadParameter(
            "--mobility-min must be below --mobility-max: the grid bins the "
            "range between them"
        )


def _validate_tof_law(tof_law: Optional[Tuple[float, float]]) -> None:
    """Validate ``--tof-law A B``: both non-negative, not both zero.

    A pair with no width would lay no axis. The bin width itself is not a
    separate flag: ``--resample-width-at-mz`` at ``--resample-reference-mz``
    fixes the bins per peak width for a ``tof`` axis exactly as it fixes the
    width of any other, and 3 bins per peak width is the default.
    """
    if tof_law is None:
        return
    a, b = tof_law
    if not (isfinite(a) and isfinite(b)):
        raise click.BadParameter(
            "The TOF width law needs two finite numbers", param_hint="tof_law"
        )
    if a < 0 or b < 0 or (a == 0 and b == 0):
        raise click.BadParameter(
            "The TOF width law needs A >= 0 and B >= 0 with at least one of "
            "them positive",
            param_hint="tof_law",
        )


def _validate_input_path(input: Path) -> None:
    """Validate input path and format requirements."""
    if not input.exists():
        raise click.BadParameter(f"Input path does not exist: {input}")

    if input.is_file() and input.suffix.lower() == ".imzml":
        ibd_path = input.with_suffix(".ibd")
        if not ibd_path.exists():
            raise click.BadParameter(
                "ImzML file requires corresponding .ibd file, but not "
                f"found: {ibd_path}"
            )
    elif input.is_dir() and input.suffix.lower() == ".d":
        from thyra.core.registry import detect_format

        # Several Bruker formats share the .d extension (timsTOF via
        # analysis.tsf/.tdf, solariX via peaks.sqlite), so delegate to the
        # registry's detection rather than duplicating the marker files
        # here -- its errors also carry the format-specific guidance
        # (e.g. the imzML-export fallback for a peaks-less solariX .d).
        try:
            detect_format(input)
        except ValueError as e:
            raise click.BadParameter(str(e)) from e


def _validate_output_path(output: Path) -> None:
    """Validate output path.

    Beyond existence, the deepest existing ancestor must be a directory
    Thyra may write into (issues #255, #312). That check is
    :func:`thyra.convert.output_ancestor_problem`, shared with
    ``convert_msi`` so the two routes refuse the same paths (issue #519).
    """
    if output.exists():
        raise click.BadParameter(f"Output path already exists: {output}")

    from thyra.convert import output_ancestor_problem

    problem = output_ancestor_problem(output)
    if problem is not None:
        raise click.BadParameter(problem)


def _validate_log_file(log_file: Optional[Path], output: Path) -> None:
    """Refuse a ``--log-file`` no folder creation can make writable.

    Checked with the other paths, before any work: a folder given as the
    log file ended in a traceback before the conversion started
    (issue #525). The rule is :func:`thyra.utils.output_files.output_file_problem`,
    the one every file a Thyra command writes goes through. A log file
    inside the store's path is refused too: creating its folder would
    create the output, and the store is renamed onto that name at the end.
    """
    if log_file is None:
        return
    problem = output_file_problem(log_file)
    store = output.resolve()
    log = log_file.resolve()
    if problem is None and (log == store or store in log.parents):
        problem = (
            f"Cannot write {log_file}: it is the output store {output} or "
            "lies inside it. Write the log file beside the store instead."
        )
    if problem is not None:
        raise click.BadParameter(problem, param_hint="'--log-file'")


def _display_calibration_info(input: Path, use_recalibrated: bool) -> None:
    """Display calibration information for Bruker datasets.

    Note: This is informational only. Full interactive selection
    will be implemented in the future (see GitHub issue #54).
    """
    from thyra.readers.bruker.timstof.timstof_reader import read_calibration_states

    # Which state is applied is the reader's own rule (issue #518): the
    # newest state that maps exactly this acquisition's frames, which is
    # not always the newest state, and none with --no-recalibrated.
    states = read_calibration_states(input, use_recalibrated)
    if not states:
        return

    # A recalibration adds a state, so the count is a property of the
    # dataset, not of a state. It used to be printed per row as
    # ``state["id"] - 1``, which made a three-state file claim the
    # second state had been recalibrated once and the third twice.
    n_recalibrations = len(states) - 1
    applied = next((state for state in states if state["applied"]), None)

    click.echo("\n" + "=" * 60)
    click.echo("Calibration Information (Display Only)")
    click.echo("=" * 60)
    if n_recalibrations:
        click.echo(
            f"  Recalibrated {n_recalibrations} time"
            f"{'s' if n_recalibrations > 1 else ''} since acquisition"
        )
    for state in states:
        marker = " (will be used)" if state["applied"] else ""
        click.echo(
            f"  State {state['id']}: {state['datetime']} "
            f"[{state['source']}]{marker}"
        )

    click.echo("\n" + _applied_calibration_line(states, applied, use_recalibrated))
    click.echo(
        "\nNote: this lists the states only; one cannot be chosen here. "
        "Thyra uses the newest state that maps this acquisition's frames, "
        "or the original calibration with --no-recalibrated."
    )
    click.echo("=" * 60 + "\n")


def _applied_calibration_line(
    states: List[Dict[str, Any]],
    applied: Optional[Dict[str, Any]],
    use_recalibrated: bool,
) -> str:
    """The one line saying which calibration the conversion applies."""
    if not use_recalibrated:
        return (
            "Using original calibration (--no-recalibrated flag set): "
            "none of these states"
        )
    newest = states[-1]["id"]
    if applied is None:
        return (
            "Using original calibration: no state maps exactly this "
            "acquisition's frames, so Bruker's library applies none of them"
        )
    if applied["id"] == newest:
        return f"Using the newest calibration state (State {newest})"
    return (
        f"Using calibration state {applied['id']}, not the newest: the newer "
        "states do not map exactly this acquisition's frames, so Bruker's "
        "library skips them"
    )


def _select_bruker_dataset(input_path: Path) -> Path:
    """Prompt user to select a dataset when multiple .d folders exist.

    If the input directory contains multiple Bruker .d folders, displays
    them as a numbered list and lets the user pick one interactively.

    Only runs count: ``.d`` folders Thyra reads, timsTOF or solariX. A
    folder holding one run is that run, so everything after this --
    the calibration display included -- sees the ``.d`` (issue #525).
    With no answer to read (nothing on standard input), the run ends
    with the one-line refusal ``convert_msi`` gives and exit 1, where
    click printed a bare ``Aborted!``.

    Args:
        input_path: The user-provided input path

    Returns:
        The selected .d folder path, or the original path if no
        selection is needed
    """
    from thyra.readers.bruker.folder_structure import (
        BrukerFolderStructure,
        several_runs_message,
    )

    d_folders = BrukerFolderStructure.runs_in(input_path)

    if not d_folders:
        return input_path
    if len(d_folders) == 1:
        return d_folders[0]

    click.echo(f"\nFound {len(d_folders)} datasets in {input_path.name}:")
    for i, d_folder in enumerate(d_folders, 1):
        click.echo(f"  [{i}] {d_folder.name}")

    try:
        choice: int = click.prompt(
            "\nSelect dataset to convert",
            type=click.IntRange(1, len(d_folders)),
        )
    except click.Abort:
        logger.error(
            "%s",
            several_runs_message(input_path, d_folders, then=", or pipe its number in"),
        )
        raise SystemExit(1) from None

    selected: Path = d_folders[choice - 1]
    click.echo(f"  -> {selected.name}\n")
    return selected


def _build_resampling_config(
    resample_method: str,
    mass_axis_type: str,
    resample_bins: Optional[int],
    resample_min_mz: Optional[float],
    resample_max_mz: Optional[float],
    resample_width_at_mz: Optional[float],
    resample_reference_mz: float,
    resample_gap_tolerance: Optional[float] = None,
    tof_law: Optional[Tuple[float, float]] = None,
) -> dict:
    """Build resampling configuration dictionary.

    ``tof_law`` is the CLI's one spelling of the two-term width law; the
    config carries it as the ``tof_a`` / ``tof_b`` pair the Python API
    takes. ``bins_per_fwhm`` is API-only: on the command line the width at
    the reference m/z says the same thing.
    """
    return {
        "method": resample_method,
        "axis_type": mass_axis_type,
        "target_bins": resample_bins,
        "min_mz": resample_min_mz,
        "max_mz": resample_max_mz,
        "width_at_mz": resample_width_at_mz,
        "reference_mz": resample_reference_mz,
        "gap_tolerance_da": resample_gap_tolerance,
        "tof_a": None if tof_law is None else tof_law[0],
        "tof_b": None if tof_law is None else tof_law[1],
    }


#: The options that decide the mass axis, which a ``--mass-axis-plan`` has
#: already fixed. Combining one with the plan would either be ignored or
#: contradict it, so it is refused.
_PLAN_FIXES = (
    "--no-resample",
    "--mass-axis-type",
    "--resample-bins",
    "--resample-min-mz",
    "--resample-max-mz",
    "--resample-width-at-mz",
    "--resample-reference-mz",
    "--tof-law",
)


def _refuse_axis_flags_with_plan(given: dict) -> None:
    """Refuse every axis option given beside ``--mass-axis-plan``.

    Args:
        given: Option spelling mapped to whether it was given, for each of
            :data:`_PLAN_FIXES`.

    Raises:
        click.BadParameter: Naming the options, since the plan fixes the axis.
    """
    clash = [flag for flag in _PLAN_FIXES if given.get(flag)]
    if clash:
        raise click.BadParameter(
            f"{', '.join(clash)} cannot be combined with --mass-axis-plan: "
            "the plan already fixes the mass axis.",
            param_hint="--mass-axis-plan",
        )


def _resampling_config_from_plan(
    plan_file: Path,
    resample_method: str,
    resample_gap_tolerance: Optional[float],
) -> dict:
    """The resampling config a ``--mass-axis-plan`` file stands for.

    ``--resample-method`` is accepted unless the plan forces a different
    method, and ``--resample-gap-tolerance`` always: neither changes the
    axis.

    Raises:
        click.BadParameter: If the file is not a plan, or the method clashes
            with the one the plan forces.
    """
    from thyra.errors import ConversionRefused
    from thyra.resampling.plan import MassAxisPlan

    try:
        plan = MassAxisPlan.load(plan_file)
    except ConversionRefused as exc:
        raise click.BadParameter(str(exc), param_hint="--mass-axis-plan") from exc
    config = plan.to_resampling_config()
    if resample_method != "auto":
        if plan.method is not None and plan.method.value != resample_method:
            raise click.BadParameter(
                f"the plan already forces the method {plan.method.value}; "
                "drop --resample-method.",
                param_hint="--resample-method",
            )
        config["method"] = resample_method
    if resample_gap_tolerance is not None:
        config["gap_tolerance_da"] = resample_gap_tolerance
    return config


def _build_reader_options(
    use_recalibrated: bool,
    intensity_threshold: Optional[float],
    spectrum_type: str = "auto",
    tdf_spectrum: Optional[str] = None,
    waters_spectrum: Optional[str] = None,
) -> dict[str, bool | float | str]:
    """Build reader options dictionary from CLI parameters.

    ``spectrum_type="auto"`` is the CLI's way of saying "no override", so it is
    omitted entirely rather than forwarded -- readers take ``None``/absent to
    mean detect, and ``"auto"`` is not a representation. ``tdf_spectrum`` is
    likewise only forwarded when given: it is a Bruker TDF option, and a
    reader for any other format would reject an unexpected keyword.
    ``waters_spectrum`` is the same for Waters .raw, and is spelled as a
    representation rather than as ``use_centroid`` because that is what a
    caller is choosing; the reader keyword it sets is the negation.
    """
    options: dict[str, bool | float | str] = {
        "use_recalibrated_state": use_recalibrated
    }
    if intensity_threshold is not None:
        options["intensity_threshold"] = intensity_threshold
    if spectrum_type != "auto":
        options["spectrum_type"] = spectrum_type
    if tdf_spectrum is not None:
        options["tdf_spectrum"] = tdf_spectrum
    if waters_spectrum is not None:
        options["use_centroid"] = waters_spectrum == "centroid"
    return options


def _reader_formats(key: str) -> Tuple[str, ...]:
    """The formats whose reader uses ``key``, from the shared table."""
    formats = READER_OPTION_FORMATS[key]
    if formats is None:
        raise KeyError(f"{key} reaches every format, so it is never ignored")
    return formats


#: Options that only reach one family of sources, and the detected
#: formats each does reach. docs/cli.md has always said the vendor groups
#: are "ignored on other formats", but nothing said so at the time: every
#: one of these was accepted on an imzML input and produced a store
#: byte-identical to the plain run, with no log line (issue #260).
#: ``--tof-law`` on a non-tof axis and ``--msms-table`` on imzML already
#: set the precedent of saying it out loud; these follow. ``--region`` is
#: not here: the reader refuses it on every other format (issue #448), and
#: an "ignored" line before that refusal said the opposite.
#: The reader options take their formats from ``READER_OPTION_FORMATS``,
#: the table ``convert_msi`` warns from, so both routes name the same
#: options as ignored (issue #519).
IGNORED_ELSEWHERE: dict[str, Tuple[str, Tuple[str, ...]]] = {
    "tdf_spectrum": ("--tdf-spectrum", _reader_formats("tdf_spectrum")),
    "waters_spectrum": ("--waters-spectrum", _reader_formats("use_centroid")),
    "use_recalibrated": (
        "--use-recalibrated/--no-recalibrated",
        _reader_formats("use_recalibrated_state"),
    ),
    "interactive_calibration": ("--interactive-calibration", ("bruker", "solarix")),
    "spectrum_type": ("--spectrum-type", _reader_formats("spectrum_type")),
    "mobility_grid": (
        "--mobility-grid/--no-mobility-grid",
        ("bruker", "imzml", "waters"),
    ),
    "mobility_bins": ("--mobility-bins", ("bruker", "imzml", "waters")),
    "mobility_min": ("--mobility-min", ("bruker", "imzml", "waters")),
    "mobility_max": ("--mobility-max", ("bruker", "imzml", "waters")),
}

#: The three options that size the grid ``--mobility-grid`` builds. They
#: do nothing without it, which docs/cli.md says and nothing else did.
GRID_SIZING_FLAGS = ("mobility_bins", "mobility_min", "mobility_max")


def _given_on_the_command_line(ctx: Optional[click.Context], name: str) -> bool:
    """Whether ``name`` was typed, rather than left at its default.

    Reads click's parameter source rather than comparing against the
    default, so ``--use-recalibrated`` (whose default is on) and
    ``--mobility-bins 256`` (whose default is the value) are recognised
    as given while an untouched option stays quiet.
    """
    if ctx is None:
        return False
    source = ctx.get_parameter_source(name)
    return source is click.core.ParameterSource.COMMANDLINE


def _warn_ignored_flags(
    ctx: Optional[click.Context], input_format: str, mobility_grid: bool
) -> None:
    """Name every option given on a source, or in a run, that ignores it.

    Args:
        ctx: The click context of the running command, or ``None`` when
            the function is called outside one (a direct unit test).
        input_format: What ``detect_format`` said the source is.
        mobility_grid: Whether a mobility grid table was asked for, which
            is what the three grid-sizing options size.
    """
    if ctx is None:
        return
    source_name = format_name(input_format)
    for name, (spelling, formats) in IGNORED_ELSEWHERE.items():
        if input_format in formats or not _given_on_the_command_line(ctx, name):
            continue
        logger.warning(
            "%s is ignored on %s input: it only applies to %s. The store is "
            "written as though it had not been given.",
            spelling,
            source_name,
            " or ".join(format_name(f) for f in formats),
        )

    if mobility_grid:
        return
    for name in GRID_SIZING_FLAGS:
        spelling, formats = IGNORED_ELSEWHERE[name]
        if input_format not in formats or not _given_on_the_command_line(ctx, name):
            continue
        logger.warning(
            "%s is ignored without --mobility-grid: it sizes the mobility "
            "grid, and no grid table is being written.",
            spelling,
        )


#: The options that shape the resampled axis or how spectra are placed on
#: it, by parameter name. ``--no-resample`` builds no such axis, so each one
#: given beside it does nothing.
RESAMPLING_FLAGS = {
    "resample_method": "--resample-method",
    "mass_axis_type": "--mass-axis-type",
    "resample_bins": "--resample-bins",
    "resample_min_mz": "--resample-min-mz",
    "resample_max_mz": "--resample-max-mz",
    "resample_width_at_mz": "--resample-width-at-mz",
    "resample_reference_mz": "--resample-reference-mz",
    "resample_gap_tolerance": "--resample-gap-tolerance",
    "tof_law": "--tof-law",
}


def _warn_resampling_flags_without_resampling(
    ctx: Optional[click.Context], resample: bool
) -> None:
    """Name every resampling option ``--no-resample`` makes dead.

    They were dropped in silence: ``--resample-min-mz 300
    --resample-max-mz 400`` beside it stored the full range and exited 0
    (issue #519). One WARNING per option, as for the vendor options
    (issue #260).
    """
    if resample or ctx is None:
        return
    for name, spelling in RESAMPLING_FLAGS.items():
        if _given_on_the_command_line(ctx, name):
            logger.warning(
                "%s is ignored with --no-resample: no resampled axis is "
                "built. The store is written as though it had not been given.",
                spelling,
            )


def _parse_streaming_option(streaming: str) -> bool | Literal["auto"]:
    """Convert the streaming CLI string to a typed value."""
    if streaming == "true":
        return True
    if streaming == "false":
        return False
    return "auto"


def _quarantine_partial_output(output: Path) -> None:
    """Move a partially written store aside after a failed conversion.

    A thin delegation: the move itself lives in
    :func:`thyra.convert.quarantine_partial_output`, because
    ``convert_msi`` is the public entry point and a library caller --
    Ousia converts through it, never through this CLI -- used to get no
    cleanup at all (issue #293).

    Kept as a call rather than deleted because ``convert_msi`` resolves
    and may extend the output path (``prepare_zarr_output_path``) while
    the CLI still holds the path the user typed. It is idempotent: a
    store already moved aside does not exist here, so this returns
    immediately. The converter writes under a temporary name and sets
    its own partial store aside (issue #525), so on that route there is
    usually nothing here to move.
    """
    from thyra.convert import quarantine_partial_output

    quarantine_partial_output(output)


def _handle_post_conversion(success: bool, output: Path) -> bool:
    """Report the conversion result, quarantining a partial store on failure.

    Returns whether the run succeeded. There are no fallible post-conversion
    steps left, so today this is just ``success`` -- but the caller keys the
    process exit status off the return value rather than off ``convert_msi``'s
    bool, so the next step added here has to report its failure instead of
    dropping it on the floor the way ``--optimize-chunks`` did.
    """
    if success:
        logger.info(f"Conversion completed successfully. Output stored at {output}")
        return True

    logger.error("Conversion failed.")
    _quarantine_partial_output(output)
    return False


class GroupedCommand(click.Command):
    """Command that groups options into sections in --help output.

    Every visible option must appear in exactly one ``GROUPS`` entry.
    Anything missing still shows up -- ``format_help`` sweeps the remainder
    into a trailing ``Options:`` section rather than dropping it -- but that
    section is a defect signal, not a home. It is where
    ``--resample-gap-tolerance`` and ``--spectrum-type`` sat for two
    releases, next to ``--version`` and ``--help``, while docs/cli.md told
    readers the help was organised by category.

    ``tests/unit/test_cli_help_grouping.py`` fails if the trailing section is
    ever non-empty again, so adding an option without classifying it here is
    caught before release rather than by a reader.

    Section order is mirrored by docs/cli.md; keep the two in step.
    """

    GROUPS = {
        "Conversion": [
            "--pixel-size",
            "--region",
            "--no-resample",
            "--resample",
            "--include-optical",
            "--no-optical",
            "--mobility-table",
            "--no-mobility-table",
            "--mobility-heatmap",
            "--no-mobility-heatmap",
            "--mobility-grid",
            "--no-mobility-grid",
            "--msms-table",
            "--no-msms-table",
        ],
        "Logging": ["--log-level", "-v", "--log-file"],
        "Ion mobility grid (advanced)": [
            "--mobility-bins",
            "--mobility-min",
            "--mobility-max",
        ],
        "Resampling (advanced)": [
            "--resample-method",
            "--mass-axis-type",
            "--resample-bins",
            "--resample-min-mz",
            "--resample-max-mz",
            "--resample-width-at-mz",
            "--resample-reference-mz",
            "--resample-gap-tolerance",
            "--tof-law",
            "--mass-axis-plan",
        ],
        "imzML-specific": ["--spectrum-type"],
        "Bruker-specific": [
            "--use-recalibrated",
            "--no-recalibrated",
            "--interactive-calibration",
            "--intensity-threshold",
            "--tdf-spectrum",
        ],
        "Waters-specific": ["--waters-spectrum"],
        "Other": ["--dataset-id", "--sample-id", "--handle-3d", "--z-spacing"],
        "General": ["--version", "--help"],
    }

    def format_help(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        """Custom help that groups options into sections."""
        self.format_usage(ctx, formatter)
        self.format_help_text(ctx, formatter)

        opts, opt_map = self._collect_options(ctx)
        if not opts:
            return

        used: set = set()
        for group_name, group_opts in self.GROUPS.items():
            records = self._records_for_group(group_opts, opt_map, used)
            if records:
                with formatter.section(group_name):
                    formatter.write_dl(records)

        remaining = [rv for param, rv in opts if id(param) not in used]
        if remaining:
            with formatter.section("Options"):
                formatter.write_dl(remaining)

    @staticmethod
    def _collect_options(ctx: click.Context) -> tuple:
        """Collect options and build name-to-record mapping."""
        opts = []
        for param in ctx.command.get_params(ctx):
            rv = param.get_help_record(ctx)
            if rv is not None:
                opts.append((param, rv))

        opt_map: dict = {}
        for param, rv in opts:
            for name in getattr(param, "opts", []) + getattr(
                param, "secondary_opts", []
            ):
                opt_map[name] = (param, rv)

        return opts, opt_map

    @staticmethod
    def _records_for_group(group_opts: list, opt_map: dict, used: set) -> list:
        """Extract help records for a group, tracking used params."""
        records = []
        for opt_name in group_opts:
            if opt_name in opt_map and id(opt_map[opt_name][0]) not in used:
                param, rv = opt_map[opt_name]
                records.append(rv)
                used.add(id(param))
        return records


@click.command(cls=GroupedCommand)
@click.version_option(version=__version__, prog_name="thyra")
@click.argument("input", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
# -- Conversion --
@click.option(
    "--format",
    type=click.Choice(["spatialdata"]),
    default="spatialdata",
    hidden=True,
    help=(
        "The only output format there is. Accepted so existing scripts keep "
        "running; it selects nothing. The Choice stays so an unknown value "
        "still fails before any data is read."
    ),
)
@click.option(
    "--pixel-size",
    type=float,
    default=None,
    help="Pixel size in um (default: auto-detect from metadata)",
)
@click.option(
    "--region",
    type=str,
    default=None,
    help=(
        "Convert one region of a timsTOF .d (default: all regions). "
        "Accepts a .mis Area Name (e.g. '03') or an integer DB "
        "RegionNumber. Refused for every other format."
    ),
)
@click.option(
    "--resample/--no-resample",
    default=True,
    help="Mass axis resampling (default: enabled)",
)
@click.option(
    "--mobility-table/--no-mobility-table",
    default=True,
    help=(
        "When the source shares one set of (m/z, ion mobility) features across "
        "pixels (an imzML export with a mobility array), also write them as a "
        "mobility-resolved sibling table next to the summed MSI table "
        "(default: enabled)"
    ),
)
@click.option(
    "--mobility-heatmap/--no-mobility-heatmap",
    default=True,
    help=(
        "When the source has an ion mobility dimension (Bruker TDF with TIMS "
        "engaged, an imzML export with a mobility array, a Waters run with "
        "drift bins), store the mean "
        "mass-mobility frame on the summed table as uns['mobility_heatmap'] "
        "(default: enabled; fed from the summed table's own passes on a "
        "Bruker TDF, one extra pass over any other source, and on a "
        "processed imzML a first read of its mobility values for the range)"
    ),
)
@click.option(
    "--mobility-grid/--no-mobility-grid",
    default=False,
    help=(
        "When the source carries ion mobility per pixel rather than as a "
        "shared feature list (Bruker TDF, processed imzML, Waters drift "
        "time), bin the point cloud onto a common mobility grid and write "
        "the same mobility-resolved sibling table (default: disabled; a "
        "much larger table, built out of core so any acquisition fits; fed "
        "from the summed table's own passes on a Bruker TDF, extra passes "
        "over other sources). Its marginal over channels reproduces the "
        "summed table exactly under the default --tdf-spectrum scan_sum, "
        "and on Waters under --waters-spectrum profile"
    ),
)
@click.option(
    "--msms-table/--no-msms-table",
    default=None,
    help=(
        "When the source isolates several precursors per pixel in disjoint "
        "mobility slices (Bruker PASEF), also write them split apart as a "
        "demultiplexed sibling table next to the summed MSI table "
        "(default: enabled, since the summed spectrum of such a pixel mixes "
        "unrelated fragment spectra; fed from the summed table's own passes, "
        "nothing on a source that is not MS/MS). The split adds back up to "
        "the summed table "
        "exactly under the default --tdf-spectrum scan_sum"
    ),
)
@click.option(
    "--include-optical/--no-optical",
    default=True,
    help="Include optical images in output (default: True)",
)
# -- Ion mobility grid (advanced) --
@click.option(
    "--mobility-bins",
    type=int,
    default=None,
    help=(
        "Mobility channels the grid divides the range into, for "
        "--mobility-grid. Default: the mass-mobility heatmap's own channels "
        f"over the same edges ({MOBILITY_CHANNELS}, or one per drift bin on a "
        "Waters run with fewer), so a box drawn on the heatmap indexes grid "
        "channels directly; a number gives that up"
    ),
)
@click.option(
    "--mobility-min",
    type=float,
    default=None,
    help=(
        "Lower edge of the mobility grid, in the axis unit (1/K0 for TIMS, "
        "ms for Waters drift time). "
        "Default: the smallest value the source's mobility axis holds, which "
        "is where the heatmap starts"
    ),
)
@click.option(
    "--mobility-max",
    type=float,
    default=None,
    help=(
        "Upper edge of the mobility grid, in the axis unit. Default: the "
        "largest value the source's mobility axis holds"
    ),
)
# -- Logging --
@click.option(
    "-v",
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
    default="INFO",
    help="Logging level (default: INFO)",
)
@click.option(
    "--log-file",
    type=click.Path(path_type=Path),
    default=None,
    help="Write logs to file",
)
# -- Performance --
@click.option(
    "--streaming",
    type=click.Choice(["auto", "true", "false"]),
    default="auto",
    hidden=True,
    help=(
        "Deprecated no-op: every conversion streams since v3.23. Accepted "
        "so existing scripts keep running."
    ),
)
@click.option(
    "--optimize-chunks",
    is_flag=True,
    hidden=True,
    help="Deprecated no-op, accepted so existing scripts keep running.",
)
# -- Resampling (advanced) --
@click.option(
    "--resample-method",
    type=click.Choice(["auto", "nearest_neighbor", "tic_preserving"]),
    default="auto",
    help="Resampling method (default: auto-detect)",
)
@click.option(
    "--mass-axis-type",
    type=click.Choice(
        ["auto", "constant", "linear_tof", "reflector_tof", "tof", "orbitrap", "fticr"]
    ),
    default="auto",
    help=(
        "Mass axis spacing type (default: auto-detect). 'tof' lays bins at a "
        "measured peak width sqrt(A m + B m^2) mDa, of which linear_tof "
        "(B=0) and reflector_tof (A=0) are the limits; the instrument's own "
        "pair applies (SELECT SERIES MRT centroid, timsTOF, PHI nanoTOF) "
        "unless --tof-law "
        "gives one."
    ),
)
@click.option(
    "--tof-law",
    type=float,
    nargs=2,
    default=None,
    metavar="A B",
    help=(
        "Coefficients of the 'tof' width law, A in mDa^2/Da and B "
        "dimensionless (MRT 0.0185 9.1e-6, timsTOF 0.0877 8.74e-4, PHI "
        "nanoTOF 0.454 0.0284). Only "
        "needed for an instrument Thyra has no pair for. Bins per peak "
        "width default to 3; --resample-width-at-mz sets them otherwise."
    ),
)
@click.option(
    "--resample-bins",
    type=int,
    default=None,
    help="Number of bins (mutually exclusive with --resample-width-at-mz)",
)
@click.option(
    "--resample-min-mz",
    type=float,
    default=None,
    help="Minimum m/z (default: auto-detect)",
)
@click.option(
    "--resample-max-mz",
    type=float,
    default=None,
    help="Maximum m/z (default: auto-detect)",
)
@click.option(
    "--resample-width-at-mz",
    type=float,
    default=None,
    help="Mass width in Da at reference m/z for physics-based binning",
)
@click.option(
    "--resample-reference-mz",
    type=float,
    default=1000.0,
    help="Reference m/z for width specification (default: 1000.0)",
)
@click.option(
    "--resample-gap-tolerance",
    type=float,
    default=None,
    help=(
        "For tic_preserving only: discard target bins farther than this many "
        "Da from any measured m/z, instead of interpolating across the gap "
        "(default: no limit)"
    ),
)
@click.option(
    "--mass-axis-plan",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Use the mass axis of a plan file (from MassAxisPlan.save) so several "
        "data sets share one axis. The plan fixes the axis, so the other axis "
        "options and --no-resample cannot be combined with it."
    ),
)
@click.option(
    "--spectrum-type",
    type=click.Choice(["auto", "profile", "centroid"]),
    default="auto",
    help=(
        "Spectrum representation (default: auto-detect from the file). "
        "Overrides what the imzML declares via MS:1000127/MS:1000128; "
        "contradicting a declaration is logged as a warning. SCiLS Lab spells "
        "this --rep_type. imzML only."
    ),
)
# -- Bruker-specific --
@click.option(
    "--use-recalibrated/--no-recalibrated",
    default=True,
    help="Use recalibrated state (default: True)",
)
@click.option(
    "--interactive-calibration",
    is_flag=True,
    help="Display Bruker calibration states",
)
@click.option(
    "--intensity-threshold",
    type=float,
    default=None,
    help="Minimum intensity filter (useful for continuous mode data)",
)
@click.option(
    "--tdf-spectrum",
    type=click.Choice(["vendor_centroid", "scan_sum"]),
    default=None,
    help=(
        "How a TDF (TIMS) frame's mobility scans are collapsed into one "
        "spectrum per pixel: scan_sum (default) sums every scan and keeps all "
        "of the ion current, matching Bruker's own per-frame total; "
        "vendor_centroid is Bruker's frame-level peak picker over the full "
        "ramp, matching TSF line spectra and SCiLS Lab. TDF only."
    ),
)
# -- Waters-specific --
@click.option(
    "--waters-spectrum",
    type=click.Choice(["centroid", "profile"]),
    default=None,
    help=(
        "What MassLynx hands back for one Waters .raw pixel: centroid is "
        "the vendor peak picker, which reports one centroid where the "
        "sampled trace has two maxima a few mDa apart; profile is that "
        "sampled trace, which keeps them apart at the cost of a larger "
        "store. The default is profile on a SELECT SERIES MRT and centroid "
        "on every other Waters instrument. Waters .raw only."
    ),
)
# -- Other --
@click.option(
    "--dataset-id",
    default="msi_dataset",
    help="Dataset identifier (default: msi_dataset)",
)
@click.option(
    "--sample-id",
    type=str,
    default=None,
    help=(
        "Name the sample in the store's metadata (default: none). "
        "Never read from the file."
    ),
)
@click.option(
    "--handle-3d",
    is_flag=True,
    help="Process as 3D data instead of 2D slices (see --z-spacing)",
)
@click.option(
    "--z-spacing",
    type=float,
    default=None,
    help=(
        "Distance between consecutive slices in um, for --handle-3d. "
        "Default: reuse the in-plane pixel size and record it as an "
        "assumption -- correct only if sections happen to be one pixel "
        "width apart. Set this whenever the section thickness is known."
    ),
)
def main(
    input: Path,
    output: Path,
    format: str,
    dataset_id: str,
    sample_id: Optional[str],
    pixel_size: Optional[float],
    handle_3d: bool,
    z_spacing: Optional[float],
    optimize_chunks: bool,
    log_level: str,
    log_file: Optional[Path],
    use_recalibrated: bool,
    interactive_calibration: bool,
    resample: bool,
    resample_method: str,
    resample_bins: Optional[int],
    resample_min_mz: Optional[float],
    resample_max_mz: Optional[float],
    resample_width_at_mz: Optional[float],
    resample_reference_mz: float,
    resample_gap_tolerance: Optional[float],
    mass_axis_type: str,
    tof_law: Optional[Tuple[float, float]],
    mass_axis_plan: Optional[Path],
    spectrum_type: str,
    include_optical: bool,
    mobility_table: bool,
    mobility_heatmap: bool,
    mobility_grid: bool,
    mobility_bins: Optional[int],
    mobility_min: Optional[float],
    mobility_max: Optional[float],
    msms_table: Optional[bool],
    intensity_threshold: Optional[float],
    tdf_spectrum: Optional[str],
    waters_spectrum: Optional[str],
    streaming: str,
    region: Optional[str],
):
    """Convert MSI data to SpatialData format.

    INPUT: Path to input MSI file or directory
    OUTPUT: Path for output file

    Subcommands (each has its own --help): 'thyra metadata INPUT'
    writes a raw source's metadata document without converting it;
    'thyra validate PATH' validates MSI metadata against the schema;
    'thyra export-metaspace PATH' writes the METASPACE submission JSON.
    """
    # Validate all parameters
    _validate_basic_params(pixel_size, dataset_id)
    _validate_resampling_params(
        resample_bins,
        resample_min_mz,
        resample_max_mz,
        resample_width_at_mz,
        resample_reference_mz,
    )
    _validate_positive_float(
        resample_gap_tolerance, "resample_gap_tolerance", "Gap tolerance"
    )
    _validate_tof_law(tof_law)
    plan_config = None
    if mass_axis_plan is not None:
        ctx = click.get_current_context(silent=True)
        _refuse_axis_flags_with_plan(
            {
                "--no-resample": not resample,
                "--mass-axis-type": mass_axis_type != "auto",
                "--resample-bins": resample_bins is not None,
                "--resample-min-mz": resample_min_mz is not None,
                "--resample-max-mz": resample_max_mz is not None,
                "--resample-width-at-mz": resample_width_at_mz is not None,
                "--resample-reference-mz": (
                    resample_reference_mz != 1000.0
                    or _given_on_the_command_line(ctx, "resample_reference_mz")
                ),
                "--tof-law": tof_law is not None,
            }
        )
        plan_config = _resampling_config_from_plan(
            mass_axis_plan, resample_method, resample_gap_tolerance
        )
    _validate_positive_float(
        intensity_threshold, "intensity_threshold", "Intensity threshold"
    )
    _validate_positive_float(z_spacing, "z_spacing", "Z spacing")
    _validate_mobility_grid_params(mobility_bins, mobility_min, mobility_max)
    _validate_input_path(input)
    _validate_output_path(output)
    _validate_log_file(log_file, output)

    # Configure logging. The log file's folder is created here, as the
    # store's is, rather than ending the run in a traceback (issue #525).
    if log_file is not None:
        make_parent_folder(log_file)
    setup_logging(log_level=getattr(logging, log_level), log_file=log_file)

    # Warn before doing any work: on a multi-hour conversion, telling the user
    # at the end that the flag was dead is too late to be useful.
    if optimize_chunks:
        logger.warning(
            "--optimize-chunks is deprecated and does nothing. Chunk sizes for "
            "the converted store are chosen at write time; the post-hoc pass "
            "this flag used to invoke never worked and has been removed. The "
            "flag is still accepted so existing scripts keep running, and will "
            "be dropped in a future release."
        )

    # If input folder has multiple .d datasets, let the user choose
    input = _select_bruker_dataset(input)

    from thyra.convert import convert_msi
    from thyra.core.registry import detect_format

    # Before any work: say which of the given options this source ignores.
    # Detection failures are left to the conversion, which reports them
    # with the guidance the registry writes for each format.
    try:
        detected_format = detect_format(input)
    except Exception:
        detected_format = ""
    if detected_format:
        _warn_ignored_flags(
            click.get_current_context(silent=True), detected_format, mobility_grid
        )
    _warn_resampling_flags_without_resampling(
        click.get_current_context(silent=True), resample
    )

    # Display calibration info if requested (Bruker datasets only)
    if interactive_calibration and input.is_dir() and input.suffix.lower() == ".d":
        _display_calibration_info(input, use_recalibrated)

    # Build resampling config if enabled
    resampling_config = plan_config or (
        _build_resampling_config(
            resample_method,
            mass_axis_type,
            resample_bins,
            resample_min_mz,
            resample_max_mz,
            resample_width_at_mz,
            resample_reference_mz,
            resample_gap_tolerance,
            tof_law,
        )
        if resample
        else None
    )

    # Build reader options for format-specific settings
    reader_options = _build_reader_options(
        use_recalibrated,
        intensity_threshold,
        spectrum_type,
        tdf_spectrum,
        waters_spectrum,
    )
    # Only the options that reach this source: the ones that do not were
    # named above, with the flag's own spelling, and convert_msi would name
    # them again. use_recalibrated_state is always set, so without this an
    # imzML run with no vendor flag got a line about one (issue #519).
    if detected_format:
        reader_options = {
            key: value
            for key, value in reader_options.items()
            if reader_option_reaches(key, detected_format)
        }

    msms_table_option: Dict[str, Any] = (
        {} if msms_table is None else {"msms_table": msms_table}
    )

    # Perform conversion. ``convert_msi`` catches its own interrupt, but
    # this catches one that lands anywhere else between here and the
    # quarantine below -- so however the run is stopped, the exit status
    # and the partial store are the ones the "Exit status" section of
    # docs/cli.md describes (issue #245).
    try:
        success = convert_msi(
            str(input),
            str(output),
            format_type=format,
            dataset_id=dataset_id,
            pixel_size_um=pixel_size,
            handle_3d=handle_3d,
            z_spacing_um=z_spacing,
            resampling_config=resampling_config,
            reader_options=reader_options,
            include_optical=include_optical,
            streaming=_parse_streaming_option(streaming),
            region=region,
            sample_id=sample_id,
            write_mobility_table=mobility_table,
            mobility_heatmap=mobility_heatmap,
            mobility_grid=mobility_grid,
            mobility_bins=mobility_bins,
            mobility_min=mobility_min,
            mobility_max=mobility_max,
            # Passed only when given: the converter's own default applies
            # otherwise, and an explicit request is what the lossless-spectrum
            # check in convert.py reacts to.
            **msms_table_option,
        )
    except KeyboardInterrupt:
        logger.error("Interrupted.")
        success = False

    ok = _handle_post_conversion(success, output)

    # Surface failure to the shell: without this a script or CI wrapper
    # calling `thyra` sees success even though nothing usable was written.
    if not ok:
        raise SystemExit(1)


def cli() -> None:
    """Console-script entry point: dispatch subcommands, else convert.

    ``thyra INPUT OUTPUT`` converts exactly as it always has; the
    metadata subcommands (``thyra metadata``, ``thyra validate``,
    ``thyra export-metaspace``) are picked off by their first argument
    before click sees it.  Dispatch is hand-rolled rather than a
    ``click.Group`` because a group cannot carry the two positional
    arguments the conversion interface is documented with.
    """
    import sys

    from thyra.metadata.schema.cli import METADATA_SUBCOMMANDS

    args = sys.argv[1:]
    if args and args[0] in METADATA_SUBCOMMANDS:
        command = METADATA_SUBCOMMANDS[args[0]]
        command.main(args[1:], prog_name=f"thyra {args[0]}")
        return
    main()


if __name__ == "__main__":
    cli()
