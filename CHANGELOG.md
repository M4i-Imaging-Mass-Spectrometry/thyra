# CHANGELOG

<!-- version list -->

## v1.16.0 (2026-03-04)

### Bug Fixes

- Correct multi-brain optical alignment for shared TIFF slides
  ([#78](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/78),
  [`77fed1d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/77fed1de8edc6129a9fb6908cfc3b928810529ec))

### Features

- Interactive dataset selection, grouped CLI help, resample default
  ([#78](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/78),
  [`77fed1d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/77fed1de8edc6129a9fb6908cfc3b928810529ec))

- Multi-brain alignment fix and CLI improvements
  ([#78](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/78),
  [`77fed1d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/77fed1de8edc6129a9fb6908cfc3b928810529ec))


## v1.15.1 (2026-02-26)

### Bug Fixes

- Code quality sweep - mypy, logging, asserts, zarr consolidation
  ([#77](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/77),
  [`fd51e6f`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/fd51e6f8989617a664e9f911502e1d101dfa28f9))


## v1.15.0 (2026-02-25)

### Features

- Multi-region support, optical alignment fixes, and image scaling
  ([#76](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/76),
  [`80c0d0b`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/80c0d0bf89d9cdf7420d6f0266c1448e8ecabece))

- Optical alignment transforms, multi-region support, and streaming fixes
  ([#76](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/76),
  [`80c0d0b`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/80c0d0bf89d9cdf7420d6f0266c1448e8ecabece))

- Optical alignment, multi-region support, and streaming fixes
  ([#76](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/76),
  [`80c0d0b`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/80c0d0bf89d9cdf7420d6f0266c1448e8ecabece))

### Refactoring

- Extract _resolve_imaging_bounds to reduce complexity
  ([#76](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/76),
  [`80c0d0b`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/80c0d0bf89d9cdf7420d6f0266c1448e8ecabece))


## v1.14.1 (2026-02-18)

### Bug Fixes

- Resolve all 287 mypy type errors across codebase
  ([#65](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/65),
  [`899d81f`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/899d81f499fb1218e46fd94a7cf58b22c3768cdf))


## v1.14.0 (2026-02-13)

### Bug Fixes

- **ci**: Fix release workflow to properly detect and publish new versions
  ([`76f9bb0`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/76f9bb08ed26fccd99c5e7fb6bbdf1c53e3546ce))

### Documentation

- Update README with Waters .raw format support
  ([#73](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/73),
  [`d372d3c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/d372d3cac6ad8d573eb9f80ef26c96acaf51900b))

### Features

- Add Waters .raw MSI reader with MassLynx native library support
  ([#73](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/73),
  [`d372d3c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/d372d3cac6ad8d573eb9f80ef26c96acaf51900b))

### Refactoring

- Reduce cyclomatic complexity in main() and _scan_all_ms_spectra
  ([#73](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/73),
  [`d372d3c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/d372d3cac6ad8d573eb9f80ef26c96acaf51900b))

### Testing

- Add comprehensive unit tests for Waters reader
  ([#73](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/73),
  [`d372d3c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/d372d3cac6ad8d573eb9f80ef26c96acaf51900b))


## v1.13.0 (2026-01-24)

### Features

- Add CLI support and tests for intensity threshold
  ([#71](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/71),
  [`ab1dfe4`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/ab1dfe4f9fc9f81c7d78d7c83ddd17c4c534d917))

- Move intensity threshold filtering to reader level
  ([#71](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/71),
  [`ab1dfe4`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/ab1dfe4f9fc9f81c7d78d7c83ddd17c4c534d917))

- Strategy pattern for instrument detection and continuous mode optimization
  ([#72](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/72),
  [`74c1f29`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/74c1f29b8d4dbdfa8a220ca800d5f7ed98fbc181))


## v1.12.1 (2026-01-23)

### Bug Fixes

- Lower PCS threshold from 50 GB to 30 GB for memory efficiency
  ([`74f6700`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/74f6700be667408b07cef2e9f3b78bcb3647a5a8))


## v1.12.0 (2026-01-23)

### Bug Fixes

- Support datasets with >2.1 billion non-zeros in streaming converter
  ([`f111ed8`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f111ed835dd4d34ee0d4c7e44b4a8412e687f99d))


## v1.11.1 (2026-01-23)

### Bug Fixes

- Handle ResamplingConfig dataclass in streaming converter
  ([`abb9d06`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/abb9d061a683be4355bb879a99fbb4dabd0a5d46))


## v1.11.0 (2026-01-23)

### Features

- Add streaming parameter to convert_msi API
  ([#70](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/70),
  [`37c7d2d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/37c7d2df313eea358873312aa51c6fe6fc61b930))


## v1.10.0 (2026-01-23)

### Bug Fixes

- Correct release workflow YAML syntax and job dependencies
  ([`5ccf1ec`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/5ccf1ec66c3f28a4f7bdc068c738ccd08a03df01))

### Code Style

- Apply black formatting to streaming converter
  ([`1432d27`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/1432d277270dc936289d215980b78e7014db1dc4))

### Features

- Add streaming converter for memory-efficient large dataset conversion
  ([`55f7d42`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/55f7d4294d133a10f0281afe2bba49e3ee93a880))

- Add streaming converter for memory-efficient large dataset conversion
  ([`fd7acf5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/fd7acf5a2aedcf31855f7dc7d39f962dc5e2e4a1))

- Implement no-cache CSC streaming for memory-efficient large dataset conversion
  ([`4bac2d1`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/4bac2d1560349bee46dbc0076fddf5a9384a3391))

### Refactoring

- Reduce _get_mass_range_complete complexity from 13 to ~5
  ([`94bcc95`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/94bcc952b90eafd95afd7fd5d67e7c9bacdeba1b))

- Reduce _stream_build_coo complexity from 16 to ~7
  ([`69b5ba9`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/69b5ba92806d7b94c4fb3644039f2e641c520eac))

- Remove dead code (zero_copy parameter and _convert_with_scipy)
  ([`d6affd5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/d6affd520d186fc1e4f21ae97fa20a1ca10d4f68))

- Streamline streaming converter code
  ([`29a834f`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/29a834fa9b1cfaccfbc14ba08168a8724e8faf5c))


## v1.9.0 (2025-12-15)

- Initial Release
