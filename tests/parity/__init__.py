"""Parity of Thyra's conversions against SCiLS Lab, measured by outputs only.

Run with ``THYRA_PARITY_DATA`` pointing at a folder of dataset folders (see
:mod:`tests.parity.datasets`)::

    pytest -m parity tests/parity
    python -m tests.parity.report

Without that variable every test here skips. ``scils_extract.py`` makes the
reference files from SCiLS Lab's ``.slx`` and needs the ``scilslab`` client,
which ships with SCiLS Lab and is not a Thyra dependency.
"""
