# CHANGELOG

<!-- version list -->

## v1.8.3 (2025-07-24)

### Bug Fixes

- Remove dry-run mode and fix failing unit tests
  ([`00f5821`](https://github.com/Tomatokeftes/Thyra/commit/00f5821babe23865324631545a86eea69d0ac5ba))

- Resolve failing metadata extractor tests
  ([`b929cff`](https://github.com/Tomatokeftes/Thyra/commit/b929cff296b01c63357d1888cc1716674189a7e0))

### Documentation

- Update TASK.md to reflect completed pixel size metadata work
  ([`5837f3a`](https://github.com/Tomatokeftes/Thyra/commit/5837f3a8f9e5db7c837585f552037ed22b45e419))

### Refactoring

- Move BaseExtractor to core module for better architecture
  ([`c30c407`](https://github.com/Tomatokeftes/Thyra/commit/c30c40709ebe3ec283680e90fdc8543ec62b9bb0))

- Remove redundant and unnecessary comments
  ([`cf86d24`](https://github.com/Tomatokeftes/Thyra/commit/cf86d248445f1979608a19238a6874d3621508c6))


## v1.8.2 (2025-07-21)

### Bug Fixes

- Correct pixel size detection metadata for interactive mode
  ([`5101d66`](https://github.com/Tomatokeftes/Thyra/commit/5101d660ad24191a03a917006a60524b8c3948c7))


## v1.8.1 (2025-07-21)

### Bug Fixes

- Correct constructor calls in reader implementations
  ([`6ac3b4b`](https://github.com/Tomatokeftes/Thyra/commit/6ac3b4b42624f8e4780c4d2e0bd5d4e7afd4f549))

- Resolve batch processing and double progress bar issues
  ([`3e06382`](https://github.com/Tomatokeftes/Thyra/commit/3e063826d5ffec99d0fc975967580c2e0047d504))

- Resolve test failures from constructor changes
  ([`c4f14cf`](https://github.com/Tomatokeftes/Thyra/commit/c4f14cfab24090392a5a28c217741e1241cf2e0f))

- Update test base classes to provide required data_path parameter
  ([`31b04b6`](https://github.com/Tomatokeftes/Thyra/commit/31b04b61748794dbf77c22da88123f7dac6edaeb))

### Refactoring

- Consolidate duplicate base readers and clean up architecture
  ([`a2aaf74`](https://github.com/Tomatokeftes/Thyra/commit/a2aaf7466a742f088a6ba67fef42edfa7863d87f))

### Testing

- Add integration test for real Bruker dataset
  ([`f9170d5`](https://github.com/Tomatokeftes/Thyra/commit/f9170d51cf9cbe00b453e63ac00261fca0b4b245))


## v1.8.0 (2025-07-20)

### Documentation

- Update TASK.md to reflect completed code quality improvements
  ([`3592d7b`](https://github.com/Tomatokeftes/Thyra/commit/3592d7b3fceaae2569d1e21c6a2ab1370edaab49))

### Features

- Add automatic pixel size detection for ImzML and Bruker formats
  ([`e20df80`](https://github.com/Tomatokeftes/Thyra/commit/e20df80179c5ae2d64cee56fbe291b308dd19f1c))

- Add pixel size detection provenance to SpatialData metadata
  ([`1ce7787`](https://github.com/Tomatokeftes/Thyra/commit/1ce778790bb15010b520ee5d7c6861835e6e90e1))

### Refactoring

- Improve code quality with comprehensive flake8 fixes
  ([`afd7f43`](https://github.com/Tomatokeftes/Thyra/commit/afd7f4339680834edf6ec88997eda6492755f22e))


## v1.7.0 (2025-07-20)

### Features

- Add missing reader properties to fix dry-run functionality
  ([`1feee33`](https://github.com/Tomatokeftes/Thyra/commit/1feee3342eb7e4b93cf7d855747e67df2bc78094))


## v1.6.0 (2025-07-20)

### Features

- Enhance package metadata and fix GitHub URLs
  ([`9f5310e`](https://github.com/Tomatokeftes/Thyra/commit/9f5310e2bdcabbc197e6b06c80f320ba10583dec))


## v1.5.1 (2025-07-19)

### Bug Fixes

- Final end-of-file formatting by pre-commit hooks
  ([`a50553f`](https://github.com/Tomatokeftes/Thyra/commit/a50553fc8635ce172ce52cc314ae11215c555071))

### Documentation

- Add commit policy to TASK.md
  ([`2475992`](https://github.com/Tomatokeftes/Thyra/commit/24759924234b59a61387d084fd3e1bceb6080b0b))

- Update TASK.md with accurate completion status
  ([`99578c2`](https://github.com/Tomatokeftes/Thyra/commit/99578c23784b9b4f6a11f6276f1e6a51fae5c195))


## v1.5.0 (2025-07-19)

### Features

- Reorganize Bruker reader and fix converter registration
  ([`49613af`](https://github.com/Tomatokeftes/Thyra/commit/49613afcc3341138269d6db2b74a5638a5ca0f16))


## v1.4.0 (2025-07-17)

### Features

- Implement Logging system
  ([`0fefbb2`](https://github.com/Tomatokeftes/Thyra/commit/0fefbb2f8a009e72394c2342586bde85e17319f7))


## v1.3.0 (2025-07-17)

### Features

- **planning**: Update project planning document with current state assessment and detailed
  architecture
  ([`ccde584`](https://github.com/Tomatokeftes/Thyra/commit/ccde5849c0249e72a05fe4c54e170813f6eabc48))


## v1.2.0 (2025-07-17)

### Documentation

- Create GEMINI.md and PLANNING.md for project context and architecture
  ([`6bdaa77`](https://github.com/Tomatokeftes/Thyra/commit/6bdaa77ea3f3bcd7217a91beeab167c48c17307e))

- Update TASK.md with current development tasks and priorities
  ([`6bdaa77`](https://github.com/Tomatokeftes/Thyra/commit/6bdaa77ea3f3bcd7217a91beeab167c48c17307e))

### Features

- **bruker**: Enhance DLL loading logic and add error handling for unsupported platforms
  ([`827723d`](https://github.com/Tomatokeftes/Thyra/commit/827723d0e2aad2e29fdfab1f31aa3d840db56c86))

- **metadata**: Add MetadataExtractor class for extracting metadata from MSI readers
  ([`6bdaa77`](https://github.com/Tomatokeftes/Thyra/commit/6bdaa77ea3f3bcd7217a91beeab167c48c17307e))

- **tests**: Implement unit tests for MetadataExtractor functionality
  ([`6bdaa77`](https://github.com/Tomatokeftes/Thyra/commit/6bdaa77ea3f3bcd7217a91beeab167c48c17307e))

- **tools**: Add ontology checking tool with CLI support and validation logic
  ([`6bdaa77`](https://github.com/Tomatokeftes/Thyra/commit/6bdaa77ea3f3bcd7217a91beeab167c48c17307e))


## v1.1.0 (2025-06-16)

### Documentation

- Update documentation for SpatialData structure and average mass spectrum access
  ([`fdb7275`](https://github.com/Tomatokeftes/Thyra/commit/fdb72752133b142187b5486aca5c959b484f43af))

### Features

- **validator**: Add CV term usage counting and reporting
  ([`7906830`](https://github.com/Tomatokeftes/Thyra/commit/790683022632df22fec85ad007312da8959aca2a))


## v1.0.0 (2025-06-16)

### Bug Fixes

- Test token permissions for release creation
  ([`839f449`](https://github.com/Tomatokeftes/Thyra/commit/839f449a2f1f5f8740f92d1ef465d56209430e9b))

- Update integration test to match new table structure
  ([`d7e7f66`](https://github.com/Tomatokeftes/Thyra/commit/d7e7f6620654f0552c266fdcf5775d2ef0760abc))

### Features

- Add versioning
  ([`f6cf960`](https://github.com/Tomatokeftes/Thyra/commit/f6cf96087f3785069abc347aecefdf00c1f604cb))


## v0.1.0 (2025-06-16)

- Initial Release
