# Thyra Test Suite

This directory contains the test suite for the `thyra` package, providing both unit tests and integration tests.

## Test Structure

The tests are organized as follows:

- `unit/`: Fast tests for individual components
  - Core functionality (registry, base classes)
  - Readers (imzML and Bruker)
  - Converters (SpatialData)
  - Utility functions

- `integration/`: End-to-end tests for the full conversion workflow
  - imzML format conversion tests
  - Bruker format conversion tests
  - Command-line interface tests

- `data/`: Test data for running the tests
  - Contains minimal test data for both imzML and Bruker formats

- `conftest.py`: Common fixtures for all tests

## Running the Tests

### Running Unit Tests

To run only the unit tests (fast, no external dependencies):

```bash
pytest
```

or explicitly:

```bash
pytest -m "not integration"
```

### Running Integration Tests

To run the integration tests:

```bash
pytest -m integration
```

### Running All Tests

To run both unit and integration tests:

```bash
pytest -m "unit or integration"
```

### Running Specific Test Files

To run tests from a specific file:

```bash
pytest tests/unit/test_registry.py
```

### Running Tests With Coverage Report

To run tests with coverage:

```bash
pytest --cov=thyra
```

For a detailed coverage report:

```bash
pytest --cov=thyra --cov-report=html
```

## Test Dependencies

In addition to the regular package dependencies, the test suite requires:

- pytest
- pytest-cov (optional, for coverage reports)
- mock (for mocking and patching)

These can be installed with:

```bash
pip install pytest pytest-cov mock
```

## Special Test Considerations

- **Bruker Tests**: Some Bruker tests require the Bruker timsdata DLL/shared library to be available. Tests will be skipped if these dependencies are not found.

- **SpatialData Tests**: Tests for SpatialData format conversion require the `spatialdata` package to be installed. Tests will be skipped if this package is not available.

- **Mock Data**: The test suite uses fixtures to create minimal test data rather than requiring large real-world datasets.

## What the local imzML corpus can and cannot exercise

`test_data/` is gitignored, so the three real imzML files live only on
developer machines and CI never sees them. Their properties were measured once,
in full, and are recorded here because re-deriving them costs a streaming scan
of a 2.1 GB XML and an 8-minute conversion — and because **the gaps are the
useful part**. Several code paths look covered and are not.

| | bellini | pea | 20240826_xenium_0041899 |
| --- | --- | --- | --- |
| imzML / `.ibd` | 21 MB / 582 MB | 29 MB / 722 MB | 2.1 GB / 14.0 GB |
| spectra | 16,384 | 12,737 | 918,855 |
| vendor | IONTOF SurfaceLab | SCiLS | SCiLS |
| layout | unindented, CRLF | indented, LF | indented, LF |
| mode | profile | centroid | centroid |
| m/z precision | 64-bit float | 64-bit float | 64-bit float |
| intensity precision | 64-bit float | 32-bit float | 32-bit float |
| convert (default settings) | 5.8 s | 8.7 s | 494 s |

### Nothing here exercises these

Any change to the paths below is untested against real data, whatever the suite
says:

- **Continuous mode.** All three files are `IMS:1000031 processed`. There is no
  continuous file in the corpus, so `get_common_mass_axis`'s spectrum-0 branch
  and `_get_mass_range_continuous` have never met one.
- **Any z dimension.** `IMS:1000052` appears zero times, so pyimzml synthesises
  `z=1` everywhere and the z-handling code only ever sees the one shape real
  data never has.
- **Compression.** All three declare `MS:1000576 no compression`.
- **Integer binary arrays.** No file declares `IMS:1000141/1000142` or
  `MS:1000519/1000522`.
- **More than one `<scanSettings>`.** All three declare `count="1"`.
- **Wrapped 32-bit offsets.** Zero negative offsets and zero inversions across
  **1,896,000** offsets, so pyimzml's `__fix_offsets` repair is dormant on the
  entire corpus. Xenium is the interesting case and it is *correct*: 71.7% of
  its offsets exceed 2^32, the largest is 14,032,038,028, and they are written
  as proper unsigned 64-bit.

Some of these are reachable through the hand-authored corpus in
`data/fixtures/` — see that directory's README for which.

### Two properties worth relying on

- **All three `.ibd` files are packed exactly.** `max(offset + length x itemsize)`
  equals the on-disk size to the byte: 581,940,736 / 721,795,516 /
  14,032,053,520. A file that fails that check is damaged, not unusual.
- **`IMS:1000104` is consistent.** The encoded byte length equals
  `IMS:1000103 x itemsize` for all ~1.9 million arrays, zero violations, which
  is what makes it usable as a cross-check on a declared precision.

### One inconsistency to know about

**bellini's spatial metadata contradicts itself.** `max count of pixels x` x
`pixel size x` is `128 x 4.40625 = 564`, but `max dimension x` says `500`.
Deriving um/pixel from `max dimension / max count` gives `3.90625`; reading
`pixel size` gives `4.40625` — 12.8% apart, with no unit on either cvParam to
arbitrate. pea and xenium are exactly self-consistent, so code that reads the
wrong one is invisible on SCiLS files and silently wrong on IONTOF files.

Unrelated but easy to trip over: bellini's m/z values all lie on a single
global linear-TOF grid (666,526 values from 300 spectra fit
`sqrt(m/z) = s0 + k*g` to within float64 noise), so it is a sparse subset of a
common axis despite being stored as `processed`. pea and xenium are genuinely
arbitrary centroids.

## Adding New Tests

When adding new functionality to the `thyra` package, please follow these guidelines for testing:

1. Add unit tests for new components, functions, or methods
2. Update integration tests if the conversion workflow is affected
3. Add new test fixtures to `conftest.py` if needed
4. Document any special requirements or dependencies for new tests
