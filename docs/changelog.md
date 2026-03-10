# Changelog

All notable changes to Thyra are documented here. This project uses
[Semantic Versioning](https://semver.org/) and
[Conventional Commits](https://www.conventionalcommits.org/).

<!-- version list -->

## v1.17.2 (2025-03-10)

### Bug Fixes

- Add `region_number` to `.obs` in streaming and 3D converters for consistent
  region tracking across all converter paths

## v1.17.1 (2025-03-10)

### Bug Fixes

- Store region info as JSON string in `uns["regions"]` to preserve dict
  structure through AnnData/Zarr round-trip (previously stored as stringified
  numpy array)

## v1.17.0 (2025-03-10)

### Features

- Include area names from `.mis` file in region info output, so
  `uns["regions"]` now contains human-readable region names alongside numbers

## v1.16.0 (2025-03-09)

### Features

- Interactive dataset selection when a folder contains multiple `.d` datasets
- Grouped `--help` output organising CLI options by category
- Resample enabled by default (`--resample / --no-resample`)

### Bug Fixes

- Correct multi-brain optical alignment for shared TIFF slides -- each
  `.d` dataset now matches its own `.mis` file

## v1.15.1 (2025-03-08)

### Bug Fixes

- Code quality sweep: mypy fixes, logging improvements, assert cleanup, zarr
  consolidation guards

## v1.15.0 (2025-03-07)

### Features

- Optical image alignment using teaching points from FlexImaging `.mis` files
- Multi-region support with per-pixel `region_number` in `.obs`
- Region metadata (names, spectra counts) stored in `uns["regions"]`
- Streaming converter fixes for large datasets

## v1.14.1 (2025-03-06)

### Bug Fixes

- Resolve all 287 mypy type errors across codebase

## v1.14.0 (2025-03-05)

### Features

- Waters `.raw` MSI reader with MassLynx native library support
- Strategy pattern for instrument detection and continuous mode optimisation
- Intensity threshold filtering moved to reader level

### Bug Fixes

- Handle `ResamplingConfig` dataclass in streaming converter
- Support datasets with >2.1 billion non-zeros in streaming converter
- Lower PCS threshold from 50 GB to 30 GB for memory efficiency
- Fix release workflow to properly detect and publish new versions

## v1.11.0 (2025-02-20)

### Features

- Add `streaming` parameter to `convert_msi()` Python API

## v1.10.0 (2025-02-18)

### Features

- Streaming converter for memory-efficient large dataset conversion
- No-cache CSC streaming for constant memory usage regardless of file size

### Refactoring

- Remove dead code (zero_copy parameter, `_convert_with_scipy`)
- Reduce `_stream_build_coo` complexity

## v1.9.0 (2025-02-10)

### Features

- Unified Bruker MSI folder structure -- BrukerReader accepts parent directory
  containing `.d` folder
- FlexImaging (Rapiflex) reader for Bruker MALDI-TOF data
- Optical image (TIFF) integration for Bruker MSI data
- Teaching point alignment module for FlexImaging optical-MSI registration
- Bruker calibration metadata support with `--interactive-calibration`
- Click-based CLI replacing argparse, with grouped help output
- Configurable sparse matrix format (`--sparse-format csc|csr`)
- `--mass-axis-type` CLI parameter for manual axis type override
- Physics-based mass axis binning with `--resample-width-at-mz`

### Bug Fixes

- Correct MSI-to-image coordinate transformation
- Use consistent pixel spacing to eliminate gaps between MSI pixels
- Prevent duplicate logging by disabling logger propagation
- Make `BrukerReader.close()` idempotent to prevent duplicate cleanup messages
- Resolve pixel size override bug with provenance tracking
- Eliminate multiple reader closures with metadata caching
- Resolve resampling mass axis override bug

### Refactoring

- Reorganise reader package into logical subfolders (bruker/timstof,
  bruker/rapiflex)
- Rename FlexImaging to Rapiflex throughout codebase
- Reduce cyclomatic complexity in alignment and converter modules
- Split SpatialData converter into modular architecture
- Complete package rename from msiconvert to thyra

## v1.8.3 (2025-07-24)

### Bug Fixes

- Remove dry-run mode and fix failing unit tests
- Resolve failing metadata extractor tests

### Refactoring

- Move BaseExtractor to core module for better architecture

## v1.8.2 (2025-07-21)

### Bug Fixes

- Correct pixel size detection metadata for interactive mode

## v1.8.1 (2025-07-21)

### Bug Fixes

- Correct constructor calls in reader implementations
- Resolve batch processing and double progress bar issues

### Refactoring

- Consolidate duplicate base readers and clean up architecture

## v1.8.0 (2025-07-20)

### Features

- Automatic pixel size detection for ImzML and Bruker formats
- Pixel size detection provenance stored in SpatialData metadata

## v1.7.0 (2025-07-20)

### Features

- Add missing reader properties to fix dry-run functionality

## v1.6.0 (2025-07-20)

### Features

- Enhance package metadata and fix GitHub URLs

## v1.5.0 (2025-07-19)

### Features

- Reorganise Bruker reader and fix converter registration

## v1.4.0 (2025-07-17)

### Features

- Structured logging system with configurable levels and file output

## v1.3.0 (2025-07-17)

### Features

- Project planning documents for architecture and roadmap

## v1.2.0 (2025-07-17)

### Features

- Enhanced DLL loading logic for Bruker with cross-platform error handling
- MetadataExtractor class for extracting metadata from MSI readers
- Ontology checking tool with CLI support (`thyra-check-ontology`)

## v1.1.0 (2025-06-16)

### Features

- CV term usage counting and reporting in validator

### Documentation

- Updated documentation for SpatialData structure and average mass spectrum access

## v1.0.0 (2025-06-16)

### Features

- Initial stable release with automated versioning

## v0.1.0 (2025-06-16)

- Initial release
