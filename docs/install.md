# Install Thyra

You install Thyra once. It takes a few minutes, and you do not need to know
how to program. You will type three commands into a terminal, and this page
shows each one.

If you already use Python, skip to [Already use Python?](#already-use-python).

## Step 1: Open a terminal

A terminal is a window where you type a command and press Enter to run it.

=== "Windows"

    Press the Windows key, type **PowerShell**, and press Enter.

=== "macOS"

    Press Command and Space together, type **Terminal**, and press Enter.

=== "Linux"

    Open the **Terminal** application from your applications menu.

## Step 2: Install uv

uv is a free, small program that installs Thyra for you. It also fetches the
right version of Python, so you do not need to install Python yourself.

Copy the line for your computer, paste it into the terminal, and press Enter:

=== "Windows"

    ```powershell
    powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
    ```

=== "macOS"

    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```

=== "Linux"

    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```

When it finishes, **close the terminal and open a new one**. The new terminal
knows where uv is.

??? advanced "Advanced: Other ways to get uv"
    Some work computers block the command above. If yours does, ask your IT
    department, or install uv another way:

    - Windows: `winget install --id=astral-sh.uv -e`
    - macOS: `brew install uv`

    The [uv installation guide](https://docs.astral.sh/uv/getting-started/installation/)
    lists every option.

## Step 3: Install Thyra

```bash
uv tool install --python 3.13 thyra
```

This downloads Thyra and everything it needs, including Python 3.13 if your
computer does not have it yet. It takes a minute or two.

## Step 4: Check that it works

```bash
thyra --version
```

You should see `thyra, version 4.0.0` or a newer number. Thyra is installed.

**Next:** [convert your first dataset](getting-started.md).

## Already use Python?

If you work in Python 3.12 or 3.13 (in Jupyter, for example), install Thyra
into that environment instead. Then you can use it from Python as well as
from the terminal:

```bash
pip install thyra
```

Python 3.14 is not supported yet. On 3.14, pip reports that it cannot find
Thyra.

??? advanced "Advanced: conda, virtual environments and installing from source"
    Thyra installs with pip into any environment that runs Python 3.12 or
    3.13. With conda:

    ```bash
    conda create -n thyra python=3.13
    conda activate thyra
    pip install thyra
    ```

    To work on Thyra itself, install it from source. The steps are in
    [Contributing](contributing.md).

## What works on which computer

Most data can be converted on any computer. Two formats need Windows or
Linux, because Thyra reads them with the instrument maker's own software,
which Thyra includes only for those systems.

| Data | Windows | Linux | macOS |
|---|---|---|---|
| imzML, Bruker solariX, Bruker rapifleX, PHI, mzPeak | yes | yes | yes |
| Bruker timsTOF | yes | yes | no |
| Waters | yes | yes | no |

On a Mac, convert timsTOF and Waters data on a Windows or Linux computer,
then copy the `.zarr` result back. The result opens on any system.

## Update Thyra

When a new version comes out, update with the same tool you installed with:

=== "uv"

    ```bash
    uv tool upgrade thyra
    ```

=== "pip"

    ```bash
    pip install --upgrade thyra
    ```

## Something went wrong?

[Troubleshooting](troubleshooting.md) covers the problems people meet most
often, starting with a `thyra` command that the terminal does not find.
