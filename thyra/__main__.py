# thyra/__main__.py

# Configure dependencies to suppress warnings BEFORE any imports
import logging  # noqa: E402
import os  # noqa: E402
import sqlite3  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Literal, Optional  # noqa: E402

import click  # noqa: E402

from thyra.convert import convert_msi  # noqa: E402
from thyra.utils.data_processors import optimize_zarr_chunks  # noqa: E402
from thyra.utils.logging_config import setup_logging  # noqa: E402

logger = logging.getLogger(__name__)

# Configure Dask to use new query planning (silences legacy DataFrame warning)
os.environ["DASK_DATAFRAME__QUERY_PLANNING"] = "True"

# Suppress dependency warnings at the earliest possible moment
warnings.filterwarnings("ignore", category=FutureWarning, module="dask")
warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
warnings.filterwarnings(
    "ignore",
    message="The legacy Dask DataFrame implementation is deprecated",
    category=FutureWarning,
)


def _get_calibration_states(bruker_path: Path) -> list[dict]:
    """Read calibration states from calibration.sqlite.

    Args:
        bruker_path: Path to Bruker .d directory

    Returns:
        List of calibration state dictionaries with id, datetime, and version info
    """
    cal_file = bruker_path / "calibration.sqlite"
    if not cal_file.exists():
        return []

    try:
        conn = sqlite3.connect(str(cal_file))
        cursor = conn.cursor()

        # Query calibration states
        cursor.execute(
            """
            SELECT cs.Id, ci.DateTime
            FROM CalibrationState cs
            LEFT JOIN CalibrationInfo ci ON cs.Id = ci.StateId
            ORDER BY cs.Id
            """
        )

        states = []
        for row in cursor.fetchall():
            state_id, datetime_str = row
            states.append(
                {
                    "id": state_id,
                    "datetime": datetime_str or "Unknown",
                    "version": state_id,
                }
            )

        conn.close()
        return states

    except Exception:
        return []


def _validate_basic_params(pixel_size: Optional[float], dataset_id: str) -> None:
    """Validate basic conversion parameters."""
    if pixel_size is not None and pixel_size <= 0:
        raise click.BadParameter("Pixel size must be positive", param_hint="pixel_size")
    if not dataset_id.strip():
        raise click.BadParameter("Dataset ID cannot be empty", param_hint="dataset_id")


def _validate_positive_int(value: Optional[int], param_name: str, label: str) -> None:
    """Validate that an optional int parameter is positive if provided."""
    if value is not None and value <= 0:
        raise click.BadParameter(f"{label} must be positive", param_hint=param_name)


def _validate_positive_float(
    value: Optional[float], param_name: str, label: str
) -> None:
    """Validate that an optional float parameter is positive if provided."""
    if value is not None and value <= 0:
        raise click.BadParameter(f"{label} must be positive", param_hint=param_name)


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

    if resample_reference_mz <= 0:
        raise click.BadParameter(
            "Reference m/z must be positive", param_hint="resample_reference_mz"
        )


def _validate_input_path(input: Path) -> None:
    """Validate input path and format requirements."""
    if not input.exists():
        raise click.BadParameter(f"Input path does not exist: {input}")

    if input.is_file() and input.suffix.lower() == ".imzml":
        ibd_path = input.with_suffix(".ibd")
        if not ibd_path.exists():
            raise click.BadParameter(
                f"ImzML file requires corresponding .ibd file, but not found: {ibd_path}"
            )
    elif input.is_dir() and input.suffix.lower() == ".d":
        if (
            not (input / "analysis.tsf").exists()
            and not (input / "analysis.tdf").exists()
        ):
            raise click.BadParameter(
                f"Bruker .d directory requires analysis.tsf or analysis.tdf file: {input}"
            )


def _validate_output_path(output: Path) -> None:
    """Validate output path."""
    if output.exists():
        raise click.BadParameter(f"Output path already exists: {output}")


def _display_calibration_info(input: Path, use_recalibrated: bool) -> None:
    """Display calibration information for Bruker datasets.

    Note: This is informational only. Full interactive selection
    will be implemented in the future (see GitHub issue #54).
    """
    states = _get_calibration_states(input)
    if not states:
        return

    click.echo("\n" + "=" * 60)
    click.echo("Calibration Information (Display Only)")
    click.echo("=" * 60)
    for state in states:
        is_active = state["id"] == max(s["id"] for s in states)
        active_marker = " (active/will be used)" if is_active else ""
        recal_info = (
            f" - recalibrated {state['version'] - 1} times"
            if state["version"] > 1
            else ""
        )
        click.echo(
            f"  State {state['id']}: {state['datetime']}{recal_info}{active_marker}"
        )

    if use_recalibrated:
        click.echo(
            f"\nUsing active calibration state (State {max(s['id'] for s in states)})"
        )
    else:
        click.echo("\nUsing original calibration (--no-recalibrated flag set)")

    click.echo("\nNote: Interactive selection not yet available. See GitHub issue #54.")
    click.echo("=" * 60 + "\n")


def _select_bruker_dataset(input_path: Path) -> Path:
    """Prompt user to select a dataset when multiple .d folders exist.

    If the input directory contains multiple Bruker .d folders, displays
    them as a numbered list and lets the user pick one interactively.

    Args:
        input_path: The user-provided input path

    Returns:
        The selected .d folder path, or the original path if no
        selection is needed
    """
    if input_path.suffix.lower() == ".d":
        return input_path

    if not input_path.is_dir():
        return input_path

    d_folders = sorted(
        f for f in input_path.iterdir() if f.is_dir() and f.suffix.lower() == ".d"
    )

    if len(d_folders) <= 1:
        return input_path

    click.echo(f"\nFound {len(d_folders)} datasets in {input_path.name}:")
    for i, d_folder in enumerate(d_folders, 1):
        click.echo(f"  [{i}] {d_folder.name}")

    choice: int = click.prompt(
        "\nSelect dataset to convert",
        type=click.IntRange(1, len(d_folders)),
    )

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
) -> dict:
    """Build resampling configuration dictionary."""
    return {
        "method": resample_method,
        "axis_type": mass_axis_type,
        "target_bins": resample_bins,
        "min_mz": resample_min_mz,
        "max_mz": resample_max_mz,
        "width_at_mz": resample_width_at_mz,
        "reference_mz": resample_reference_mz,
    }


def _build_reader_options(
    use_recalibrated: bool,
    intensity_threshold: Optional[float],
) -> dict[str, bool | float]:
    """Build reader options dictionary from CLI parameters."""
    options: dict[str, bool | float] = {"use_recalibrated_state": use_recalibrated}
    if intensity_threshold is not None:
        options["intensity_threshold"] = intensity_threshold
    return options


def _parse_streaming_option(streaming: str) -> bool | Literal["auto"]:
    """Convert the streaming CLI string to a typed value."""
    if streaming == "true":
        return True
    if streaming == "false":
        return False
    return "auto"


def _handle_post_conversion(
    success: bool,
    optimize_chunks: bool,
    format: str,
    output: Path,
    dataset_id: str,
) -> None:
    """Handle chunk optimization and result logging after conversion."""
    if success and optimize_chunks and format == "spatialdata":
        optimize_zarr_chunks(str(output), f"tables/{dataset_id}/X")

    if success:
        logger.info(f"Conversion completed successfully. Output stored at {output}")
    else:
        logger.error("Conversion failed.")


class GroupedCommand(click.Command):
    """Command that groups options into sections in --help output."""

    GROUPS = {
        "Conversion": [
            "--format",
            "--pixel-size",
            "--region",
            "--no-resample",
            "--resample",
            "--include-optical",
            "--no-optical",
        ],
        "Logging": ["--log-level", "-v", "--log-file"],
        "Resampling (advanced)": [
            "--resample-method",
            "--mass-axis-type",
            "--resample-bins",
            "--resample-min-mz",
            "--resample-max-mz",
            "--resample-width-at-mz",
            "--resample-reference-mz",
        ],
        "Performance": ["--streaming", "--optimize-chunks", "--sparse-format"],
        "Bruker-specific": [
            "--use-recalibrated",
            "--no-recalibrated",
            "--interactive-calibration",
            "--intensity-threshold",
        ],
        "Other": ["--dataset-id", "--handle-3d"],
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
@click.argument("input", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
# -- Conversion --
@click.option(
    "--format",
    type=click.Choice(["spatialdata"]),
    default="spatialdata",
    help="Output format (default: spatialdata)",
)
@click.option(
    "--pixel-size",
    type=float,
    default=None,
    help="Pixel size in um (default: auto-detect from metadata)",
)
@click.option(
    "--region",
    type=int,
    default=None,
    help="Convert specific region number (default: all regions)",
)
@click.option(
    "--resample/--no-resample",
    default=True,
    help="Mass axis resampling (default: enabled)",
)
@click.option(
    "--include-optical/--no-optical",
    default=True,
    help="Include optical images in output (default: True)",
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
    help="Streaming mode for large datasets (default: auto)",
)
@click.option(
    "--optimize-chunks",
    is_flag=True,
    help="Optimize Zarr chunks after conversion",
)
@click.option(
    "--sparse-format",
    type=click.Choice(["csc", "csr"]),
    default="csc",
    help="Sparse matrix format: csc or csr (default: csc)",
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
        ["auto", "constant", "linear_tof", "reflector_tof", "orbitrap", "fticr"]
    ),
    default="auto",
    help="Mass axis spacing type (default: auto-detect)",
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
# -- Other --
@click.option(
    "--dataset-id",
    default="msi_dataset",
    help="Dataset identifier (default: msi_dataset)",
)
@click.option(
    "--handle-3d",
    is_flag=True,
    help="Process as 3D data instead of 2D slices",
)
def main(
    input: Path,
    output: Path,
    format: str,
    dataset_id: str,
    pixel_size: Optional[float],
    handle_3d: bool,
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
    mass_axis_type: str,
    sparse_format: str,
    include_optical: bool,
    intensity_threshold: Optional[float],
    streaming: str,
    region: Optional[int],
):
    """Convert MSI data to SpatialData format.

    INPUT: Path to input MSI file or directory
    OUTPUT: Path for output file
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
        intensity_threshold, "intensity_threshold", "Intensity threshold"
    )
    _validate_input_path(input)
    _validate_output_path(output)

    # Configure logging
    setup_logging(log_level=getattr(logging, log_level), log_file=log_file)

    # If input folder has multiple .d datasets, let the user choose
    input = _select_bruker_dataset(input)

    # Display calibration info if requested (Bruker datasets only)
    if interactive_calibration and input.is_dir() and input.suffix.lower() == ".d":
        _display_calibration_info(input, use_recalibrated)

    # Build resampling config if enabled
    resampling_config = (
        _build_resampling_config(
            resample_method,
            mass_axis_type,
            resample_bins,
            resample_min_mz,
            resample_max_mz,
            resample_width_at_mz,
            resample_reference_mz,
        )
        if resample
        else None
    )

    # Build reader options for format-specific settings
    reader_options = _build_reader_options(use_recalibrated, intensity_threshold)

    # Perform conversion
    success = convert_msi(
        str(input),
        str(output),
        format_type=format,
        dataset_id=dataset_id,
        pixel_size_um=pixel_size,
        handle_3d=handle_3d,
        resampling_config=resampling_config,
        reader_options=reader_options,
        sparse_format=sparse_format,
        include_optical=include_optical,
        streaming=_parse_streaming_option(streaming),
        region=region,
    )

    _handle_post_conversion(success, optimize_chunks, format, output, dataset_id)


if __name__ == "__main__":
    main()
