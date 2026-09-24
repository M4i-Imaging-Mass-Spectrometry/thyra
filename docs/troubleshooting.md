# Troubleshooting

Find the message or problem you see, then follow the fix under it. The
explanations are folded away in the **Advanced** boxes, for when you want
them.

## Installing

### The terminal says `thyra` is not recognized, or "command not found"

1. Close the terminal and open a new one. A new terminal picks up newly
   installed commands.
2. If you installed with uv, run `uv tool update-shell`, then open a new
   terminal again.
3. If you installed with pip, type `python -m thyra` wherever this site says
   `thyra`. It runs the same program:

    ```bash
    python -m thyra --version
    ```

??? advanced "Advanced: Why this happens"
    A terminal finds commands by looking through a list of folders, called
    the PATH. The installer puts `thyra` in a folder that is not always on
    that list yet. `uv tool update-shell` adds uv's folder to it.
    `python -m thyra` avoids the list altogether, because it asks Python to
    run Thyra directly.

### pip says "No matching distribution found for thyra"

Your Python version is not one Thyra supports. Thyra needs Python 3.12 or
3.13; Python 3.14 does not work yet. The simplest fix is to install with uv,
which fetches Python 3.13 for you: see [Install](install.md).

### `thyra --version` shows a version starting with 1

Your Python is 3.11, so pip installed the last old version of Thyra that ran
on it. Install with uv instead (see [Install](install.md)), or use Python
3.12 or 3.13.

### pip says "externally-managed-environment"

Your system's Python does not allow pip to install into it. This is common
on Linux and with Homebrew on macOS. Install with uv instead:

```bash
uv tool install --python 3.13 thyra
```

## Converting

### "Output path already exists"

Thyra never overwrites an earlier result. Give the new result another name,
or delete the old `.zarr` folder first.

### "Pixel size not found in metadata"

Your file does not record its pixel size. Run the conversion again and give
the size in micrometres:

```bash
thyra input.imzML output.zarr --pixel-size 50
```

### "Bruker SDK is not supported on macOS"

Thyra cannot read Bruker timsTOF data on a Mac. Convert it on a Windows or
Linux computer, then copy the `.zarr` result to your Mac. The result opens
anywhere.

On Windows or Linux, an error about the Bruker SDK means the copy of Bruker's
library that comes with Thyra did not load. Reinstall Thyra, and if the error
stays, [ask for help](#ask-for-help).

### An imzML file is refused before the conversion starts

Thyra checks an imzML file before converting it, and stops if it could not
read the data correctly. The message names the problem.

- **The message says the `.ibd` is truncated.** The `.ibd` file is
  incomplete, usually because a copy did not finish. Copy it again and
  compare the file size with the original.
- **Any other message.** Export the data again from your instrument
  software, choosing uncompressed data and 64-bit floating-point m/z values.

??? advanced "Advanced: Everything Thyra checks in an imzML file"
    Thyra compares what the `.imzML` file declares with the `.ibd` file on
    disk, before it reads a single spectrum. It refuses:

    - an `.ibd` shorter than the offsets in the `.imzML` say it should be;
    - `zlib` compression (`MS:1000574`) on either array, which the reader
      cannot decompress;
    - a data type with no precision, with two precisions, or with one that
      disagrees with what the reader resolved;
    - 64-bit integer arrays, whose size differs between Windows and Linux;
    - a first spectrum whose encoded length (`IMS:1000104`) does not match its
      number of values at the declared precision;
    - a negative offset or length, or a spectrum whose m/z and intensity
      arrays have different lengths.

    Some unusual files are fine and only produce a warning: extra bytes at the
    end of the `.ibd`, offsets out of order, and more than one
    `<scanSettings>` block. 32-bit integer arrays are allowed.

### The conversion failed and left a folder ending in `.failed`

When a conversion fails, Thyra renames the unfinished result, for example to
`brain.zarr.failed`. That keeps it out of the way, so you can run the same
command again straight away. The `.failed` folder is only useful for working
out what went wrong: delete it when you no longer need it.

### "WinError 5: Access is denied" (Windows)

Another program has the result open. Close every program that has loaded the
`.zarr` folder, such as napari, a Jupyter notebook or a Python window. Then
run the conversion again, or give the result a new name.

??? advanced "Advanced: Why this happens"
    Windows cannot replace a file in one step, so the storage library deletes
    the old file and renames the new one into its place. Now and then,
    Windows refuses that delete for a moment. The library retries, which
    clears the short refusals on its own. An error that remains after the
    retries means a program really is holding the folder open. Antivirus
    software is not the cause: the same thing happens with real-time
    scanning turned off.

### The computer runs out of memory

Thyra keeps the spectra on disk while it works, but the list of m/z bins
stays in memory. Use fewer bins:

```bash
thyra large.d output.zarr --resample-bins 20000
```

## Opening the result

### Windows: a result looks incomplete in another program { #windows-long-paths }

Windows limits a normal path to 260 characters. A `.zarr` folder holds files
deep inside it, so a result stored in a deep folder can go past that limit.
Thyra writes such results correctly. Other programs may read them with parts
missing, without any error.

The simplest fix is a short location, such as `C:\msi\brain.zarr`. Always
give programs the full path, starting with the drive letter.

??? advanced "Advanced: Other fixes, and why a relative path matters"
    Thyra's own commands (`thyra validate`, `thyra export-metaspace`) and
    `thyra.metadata.schema.read_msi_metadata_blocks` handle long paths by
    themselves. For other tools:

    - **Turn on long paths in Windows.** Set `LongPathsEnabled` to 1 under
      `HKLM\SYSTEM\CurrentControlSet\Control\FileSystem`. This needs
      administrator rights.
    - **Prepare the path in Python.** Pass the store through
      `thyra.utils.windows_paths.prepare_zarr_read_path`, or add the `\\?\`
      prefix yourself:

        ```python
        import spatialdata as sd
        sdata = sd.read_zarr(r"\\?\C:\very\long\path\output.zarr")
        ```

    Without long paths turned on, a relative path such as
    `..\results\brain.zarr` is measured before Windows removes the `..`
    parts. It can go over the limit even when the full path does not. If
    something looks missing, compare the two spellings:

    ```python
    import os
    os.path.exists(p), os.path.exists(os.path.abspath(p))
    ```

    `False, True` means the path is the problem, not the data.

## Getting help

### See more of what Thyra is doing

Add `-v DEBUG` to print every step, and `--log-file` to save it:

```bash
thyra input.imzML output.zarr -v DEBUG --log-file conversion.log
```

### Ask for help

[Open an issue on GitHub](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/issues).
Say what you ran and what you expected, and attach the log file from the
command above.
