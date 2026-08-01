# CHANGELOG

<!-- version list -->

## v2.2.4 (2026-08-01)

### Documentation

- **resampling**: Source the SCiLS claims, and stop recommending a forbidden pairing
  ([`90df8cc`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/90df8cc837ed10e0df569a501d9e7550a038ad6f))


## v2.2.3 (2026-08-01)

### Bug Fixes

- **metadata**: Read the declared spectrum representation, don't infer it
  ([`d6c080e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/d6c080e2b8caf1b7cfbbb402fe6809c7e6b83c9d))


## v2.2.2 (2026-08-01)

### Bug Fixes

- **resampling**: Stop handing unknown modalities MALDI-TOF resampling
  ([`f98bc3e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f98bc3e56407b0554e3e1df0842530a7e8db1164))


## v2.2.1 (2026-08-01)

### Bug Fixes

- **imzml**: Parse with ElementTree so unindented CRLF files load
  ([`b9fe017`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/b9fe01726efbcf5270bebe49af10b6dbfdeb5da6))


## v2.2.0 (2026-08-01)

### Features

- **resampling**: Give tic_preserving a gap tolerance, as Cardinal and matter have
  ([`954c95d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/954c95d7f2955adeea03a65d048010a8e25b8e94))


## v2.1.0 (2026-08-01)

### Documentation

- **handouts**: Record the SCiLS baseline for resampling
  ([`35fc7e7`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/35fc7e72aecd8c57494d9eaffae208f486739982))

### Features

- **readers**: Default the raw mass axis cap to SCiLS's 10 million bins
  ([`93eb229`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/93eb229371bd07ba5d43647fc3ea236c59ae3265))


## v2.0.4 (2026-08-01)

### Bug Fixes

- **resampling**: Stop dropping negative bins, and correct the docs
  ([`96938e3`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/96938e3c60b3331106ff14e06f027c6f2fcdc430))


## v2.0.3 (2026-07-31)

### Bug Fixes

- **streaming**: Route on the real bin count, not a guessed one
  ([`3ad7c0a`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/3ad7c0abe21d286c7e49ef488e0a0f419fa3ee62))


## v2.0.2 (2026-07-31)

### Bug Fixes

- **streaming**: Widen CSR column indices and build the string indexes lazily
  ([`edf447a`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/edf447a6a70868640238a9c6a2737f9ebb6c0b4c))


## v2.0.1 (2026-07-31)

### Bug Fixes

- **resampling**: Make tic_preserving actually preserve TIC
  ([`106f15e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/106f15e62c80d1ca107b48103178750c7b75cfb6))

### Documentation

- **handouts**: Mark the TIC finding as fixed
  ([`4c88657`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/4c88657af0caf481c0c0cafba4900474d121b052))

- **handouts**: Report findings for the write-amplification investigation
  ([`6fd1798`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/6fd1798f940a140e7e356c8f07682bc20922e4cc))

### Performance Improvements

- **imzml**: Stream the processed-mode raw mass axis build
  ([`b8562b4`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/b8562b43c38d8bc250c026c5353f8c31e51b5277))

### Refactoring

- **resampling**: Share one TIC-preservation rule between both paths
  ([`e2a2b58`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/e2a2b580445fc739360fe0029c3b04fe33ef7207))


## v2.0.0 (2026-07-31)

### Build System

- Drop Python 3.11 and upgrade the spatialdata cluster
  ([`0ef2263`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/0ef2263e981cd91427bb284563c596ce0f572f75))

### Documentation

- Correct the anndata pandas-3 fix status
  ([`d33460e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/d33460e7d8d05d5facff3a9c5cfbdba86959de28))

- Note that the anndata escape hatches work in combination on pandas 3
  ([`54c0393`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/54c0393556fa39547e7ceaa0cf919bcd4b173e88))

- Record the real anndata 0.13 prerequisites in the pin comment
  ([`951319f`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/951319fa349441334c6de0690346a47c06594f47))

- **handouts**: Refresh handout E after the #1055 rebase
  ([`87209aa`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/87209aaef1bc6021e59807df3199e2244bbbc598))

### Testing

- Make tests a regular package so site-packages cannot shadow it
  ([`f68efee`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f68efee18780b266ed7e2b424c0cff3ea7af3c65))

### Breaking Changes

- Python 3.11 is no longer supported. The minimum is now 3.12, required by spatialdata 0.8.0 and
  anndata 0.13.


## v1.27.3 (2026-07-30)

### Bug Fixes

- **ci**: Mark the complexity monitor's subprocess import as reviewed
  ([`bfe811b`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/bfe811b5a71e169b67a0cdc42a1cf01f3f54706c))

### Build System

- Check out LF line endings on every platform
  ([`04879ab`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/04879ab2dd1bfddae2837f536062635f0804d570))

- Make pyproject.toml the single source of truth for pydocstyle
  ([`eecabb2`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/eecabb2537fcd4bd3c18e739c68fc9ce1f48af42))

### Documentation

- **handouts**: Record how the toolchain-hygiene lane resolved
  ([`41d9c62`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/41d9c62639574b42d16946059585cecb2a191c79))

### Refactoring

- **types**: Drop the two no-any-return errors from a full mypy run
  ([`f207f8d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f207f8d22c1adbb25b41c3582c22190a31273139))


## v1.27.2 (2026-07-30)

### Bug Fixes

- **cli**: Retire the --optimize-chunks pass that never ran
  ([`1d365c0`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/1d365c04c2f4f02ec7bbbd59b48c33c03e649b4e))

### Documentation

- **handouts**: Add handout F and record the branch sweep
  ([`c62afc9`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/c62afc9ad1ee18963b8bf651288d1a21e9832b14))

- **handouts**: Record how the optimize-chunks lane resolved
  ([`a868b41`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a868b41001d3b683f3b53112dc8d02963e5e2270))

### Testing

- **streaming**: Cover the encoding attrs the PCS layout hand-writes
  ([`9786359`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/978635921bee796b7821dd22af565ebb958612d6))


## v1.27.1 (2026-07-30)

### Bug Fixes

- **spatialdata**: Coerce pandas 3 string dtypes to object on write
  ([`86ffb1d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/86ffb1d5dd0af817b02502f4f3ba0bd121098d66))

### Documentation

- Add handouts for the four parallel workstreams
  ([`70be229`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/70be229a0ff7cf989b840da308a6b12cb76f3708))

- **handouts**: Fold in findings from spatialdata PR #1055
  ([`2d0851a`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/2d0851a8c85df7888eec3afdc659c353d761772e))


## v1.27.0 (2026-07-30)

### Bug Fixes

- **cli**: Exit non-zero when conversion fails
  ([`1584f10`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/1584f106af07875b812102af36ef693bf8815c8c))

- **convert**: Write through an extended-length path when the store is deep
  ([`f9c0199`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f9c019938fbeab1133efeb0fdcb6722400036caf))

- **converters**: Derive bin count from axis physics for FT-ICR and Orbitrap
  ([`20e39de`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/20e39defa1b0e87e3278a54d65362a03027544d8))

- **converters**: Reject unknown resampling method and axis type values
  ([`4840b7d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/4840b7da5670b7357738858a60eb1eedcfa2979c))

- **converters**: Retry Zarr metadata renames blocked on Windows
  ([`af19db6`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/af19db6b42a0cc9debf38e82a3f895f3597cb2dc))

- **resampling**: Align the ResamplingConfig reference_mz default with the CLI
  ([`de33301`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/de33301011d4e54af196e79326a728f33f122f7b))

- **resampling**: Emit Orbitrap and FT-ICR axes in ascending m/z order
  ([`40a191e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/40a191e4a1ed14a4ccabb335038b69b30819b0a0))

### Build System

- Make .flake8 the single source of truth for lint settings
  ([`e144f2d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/e144f2d3b626342ddd58d7986def874b5fe01dbb))

### Features

- **cli**: Add --version
  ([`b19f989`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/b19f989ab7419ab1b4416b8c1431e61809fff50b))

### Refactoring

- **resampling**: Drop the dead axis-generator width and bins methods
  ([`c6bf58d`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/c6bf58d5b2035e4576389bd2f3731d66870aa221))

- **resampling**: Drop the unimplemented LINEAR_INTERPOLATION method
  ([`23f8f76`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/23f8f7687273db7ad5be1eeb128a3bab90522894))


## v1.26.0 (2026-07-30)

### Documentation

- Add tutorial and resampling pages, publish the example notebook
  ([`45b3680`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/45b3680549d82a37e918bfc35ba5f924239dc04a))

### Features

- **tools**: Add thyra-example-data synthetic dataset generator
  ([`6566a77`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/6566a7721ce7384dc81b385e9dd06f3587527cdd))


## v1.25.5 (2026-07-09)

### Bug Fixes

- **deps**: Add upper bounds so pip install lands on a working set
  ([#102](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/102),
  [`0210b0e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/0210b0ed55199b81b7947e8313e48d82da933492))


## v1.25.4 (2026-06-22)

### Bug Fixes

- **streaming**: Write PCS images via spatialdata and fail loudly on error
  ([#100](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/100),
  [`5347a76`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/5347a76746406efa9fd838b1e673a3c10a5a68ec))


## v1.25.3 (2026-06-22)

### Bug Fixes

- **convert**: Revert RAM-aware auto-streaming threshold (#98)
  ([#99](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/99),
  [`8474134`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/847413481213c0a264b6ddd3c97ee170aea3bf41))


## v1.25.2 (2026-06-22)

### Bug Fixes

- **convert**: Make auto-streaming threshold RAM-aware
  ([#98](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/98),
  [`4e82aa3`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/4e82aa3a890a211d585a7925723548ec4d392886))


## v1.25.1 (2026-06-22)

### Bug Fixes

- **deps**: Declare defusedxml as runtime dependency for Bruker reader
  ([#97](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/97),
  [`2bc4f86`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/2bc4f86d87ff1d32c7c0f735f1826be511cc0b6b))


## v1.25.0 (2026-06-12)

### Bug Fixes

- **streaming**: Clean up temp storage on every convert() exit path
  ([#96](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/96),
  [`c56f0bb`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/c56f0bbff491b704f64e907e93db49bd1a347a91))

### Documentation

- Credit Nepsis Scriptorium and place logotype on docs landing
  ([`b6e629e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/b6e629e890b9857c7b372982c0eaf244d1bef7e7))

- Switch brand assets from PNG to SVG
  ([`344e152`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/344e1524dafeb78e0cdea5c642173210b98d8b2b))

### Features

- **spatialdata**: Centralize image chunk policy as the sharding seam (foundation)
  ([#95](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/95),
  [`59eb93b`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/59eb93b8bacb2dd5f8bf9feabb25439ee9ea53dc))


## v1.24.0 (2026-06-01)


## v1.23.0 (2026-05-19)

### Bug Fixes

- **resampling**: Permissive timsTOF instrument-name match
  ([`2a200db`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/2a200db4571d8d3cf3d59a6a0fa108ab8906e6b9))

### Features

- **preview**: Add thyra.preview_msi metadata-only shim
  ([`2d06d1f`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/2d06d1f41e2cec63387b9e73c50e713a44c2fb68))

- **preview**: Metadata-only mode for BrukerReader (skip SDK init)
  ([`8f23ef7`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/8f23ef778f0084d4cf99f849a0f97120cd694f7e))

### Performance Improvements

- **preview**: Skip SUM(NumPeaks) scan in metadata_only mode
  ([`ac3f214`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/ac3f21495c26c27138d185e7b08e30527b192f29))


## v1.22.0 (2026-05-06)

### Bug Fixes

- Align global coordinate system across image and shapes
  ([#93](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/93),
  [`36938e7`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/36938e7facec910b9be095f60608931001f3e8a9))

### Documentation

- Document the coordinate-system contract
  ([#93](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/93),
  [`36938e7`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/36938e7facec910b9be095f60608931001f3e8a9))

### Features

- Unified "global" coordinate-system contract for produced zarrs
  ([#93](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/93),
  [`36938e7`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/36938e7facec910b9be095f60608931001f3e8a9))


## v1.21.0 (2026-05-04)

### Bug Fixes

- Bruker region selection, pixel size, and size estimator bugs
  ([#92](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/92),
  [`f4177d5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f4177d5f13d9d01cc585f95f91d80ce73285d655))

- Cast function_types dict keys to string for zarr serialization
  ([#92](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/92),
  [`f4177d5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f4177d5f13d9d01cc585f95f91d80ce73285d655))

- Compute size estimator bins from real resampling axis
  ([#92](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/92),
  [`f4177d5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f4177d5f13d9d01cc585f95f91d80ce73285d655))

- Drop empty pixel rows from obs when polygon != bbox
  ([#92](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/92),
  [`f4177d5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f4177d5f13d9d01cc585f95f91d80ce73285d655))

- Prefer .mis Raster over BeamScanSize for pixel size
  ([#92](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/92),
  [`f4177d5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f4177d5f13d9d01cc585f95f91d80ce73285d655))

### Features

- Accept .mis Area Names on --region; surface mapping
  ([#92](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/92),
  [`f4177d5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f4177d5f13d9d01cc585f95f91d80ce73285d655))

### Testing

- Attach handler directly to module logger for log capture
  ([#92](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/92),
  [`f4177d5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f4177d5f13d9d01cc585f95f91d80ce73285d655))

- Capture log via root logger to fix CI flake
  ([#92](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/92),
  [`f4177d5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/f4177d5f13d9d01cc585f95f91d80ce73285d655))


## v1.20.3 (2026-03-24)

### Performance Improvements

- Vectorise _create_coordinates_dataframe with numpy instead of Python loop
  ([#86](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/86),
  [`7eeaf40`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/7eeaf40a18d8463df678b49106420ff4d47f44eb))

### Refactoring

- Clean up convert.py and align flake8 line length with black
  ([#86](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/86),
  [`7eeaf40`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/7eeaf40a18d8463df678b49106420ff4d47f44eb))

- Code quality improvements ([#86](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/86),
  [`7eeaf40`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/7eeaf40a18d8463df678b49106420ff4d47f44eb))

- Normalise resampling config to ResamplingConfig at init, remove isinstance branching
  ([#86](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/86),
  [`7eeaf40`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/7eeaf40a18d8463df678b49106420ff4d47f44eb))

- Remove duplicate _suppress_upstream_warnings from streaming_converter
  ([#86](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/86),
  [`7eeaf40`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/7eeaf40a18d8463df678b49106420ff4d47f44eb))

- Remove hardcoded path, clean up type ignores, delete dead integration test
  ([#86](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/86),
  [`7eeaf40`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/7eeaf40a18d8463df678b49106420ff4d47f44eb))

### Testing

- Add unit tests for MSIRegistry format detection and registration
  ([#86](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/86),
  [`7eeaf40`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/7eeaf40a18d8463df678b49106420ff4d47f44eb))

- Clean up test suite - remove dead tests, mark integration, add unit tests
  ([#86](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/86),
  [`7eeaf40`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/7eeaf40a18d8463df678b49106420ff4d47f44eb))


## v1.20.2 (2026-03-21)

### Bug Fixes

- Compute bounding box from all points in polygon Area definitions
  ([#85](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/85),
  [`76741a5`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/76741a5bf12c12b24ab7d62e450767078ff3d9db))


## v1.20.1 (2026-03-17)

### Bug Fixes

- Iterate actual frame IDs from MaldiFrameInfo instead of assuming sequential 1..N
  ([#83](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/83),
  [`91642cb`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/91642cbdca502933b3c76f04210f8c53dcc87e36))

- Pass per-region avg spectrum through COO path
  ([#83](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/83),
  [`91642cb`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/91642cbdca502933b3c76f04210f8c53dcc87e36))

- Pass per-region avg spectrum through COO path data structures
  ([#83](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/83),
  [`91642cb`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/91642cbdca502933b3c76f04210f8c53dcc87e36))


## v1.20.0 (2026-03-11)

### Features

- Store thyra_version in essential_metadata
  ([#82](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/82),
  [`7c309e2`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/7c309e215701d5b26bb85bc1282010077cc749c3))


## v1.19.0 (2026-03-10)

### Bug Fixes

- Extract region accumulation helpers to reduce complexity
  ([#81](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/81),
  [`159494e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/159494eb9cbdf518bfa69638142ab24ff130fbd6))

### Features

- Compute and store per-region mean spectrum for multi-region datasets
  ([#81](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/81),
  [`159494e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/159494eb9cbdf518bfa69638142ab24ff130fbd6))

- Per-region mean spectrum for multi-region datasets
  ([#81](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/81),
  [`159494e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/159494eb9cbdf518bfa69638142ab24ff130fbd6))


## v1.18.2 (2026-03-10)

### Bug Fixes

- Update Bruker SDK binaries to support TSF/MALDI on Linux
  ([#80](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/80),
  [`b629c0a`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/b629c0a81806c3e1865d31d7339ba833d7feb037))

### Documentation

- Add Bruker SDK license and attribution for bundled binaries
  ([#80](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/80),
  [`b629c0a`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/b629c0a81806c3e1865d31d7339ba833d7feb037))


## v1.18.1 (2026-03-10)

### Bug Fixes

- Prevent duplicate PyPI publish on concurrent release runs
  ([`40884a6`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/40884a68f6106f43e8640004b8fc44840e0f257a))

### Documentation

- Add prominent documentation links to top of README
  ([`b58f27b`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/b58f27b7bc2f80f1c0b95ae581540acf23e1225e))


## v1.18.0 (2026-03-10)

### Bug Fixes

- Correct changelog dates, dark mode, nav tabs, mypy version
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

- Dark mode admonitions, navigation polish, mypy version
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

- Docs review sweep - accuracy, mobile, and packaging fixes
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

- Final docs sweep - version gap note, obs columns, kwargs type
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

- Pin mkdocs <2 and change streaming default to auto
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

- Resolve tab visibility and layout issues in docs theme
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

### Code Style

- Brand admonition colours for note, tip, info, and warning
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

- Clean tab bar beneath gradient header
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

### Documentation

- Add MkDocs Material documentation site
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

- Complete API reference with metadata, resampling, and registry
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

- Comprehensive documentation overhaul
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))

### Features

- Apply brand colour scheme and typography to docs
  ([#79](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/pull/79),
  [`a2d662c`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/a2d662c59d9dfe67cd713ca095ca1997c3307687))


## v1.17.2 (2026-03-10)

### Bug Fixes

- Add region_number to obs in streaming and 3D converters
  ([`3dea7ca`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/3dea7caf6c07713694ae8daa95dff3d4b414133b))


## v1.17.1 (2026-03-09)

### Bug Fixes

- Store region info as JSON to preserve dict structure in zarr
  ([`49d674e`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/49d674e5e72a23c4700a6efdcb5721dd9ee87cb4))


## v1.17.0 (2026-03-09)

### Features

- Include area names from .mis file in region info output
  ([`2c34961`](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/commit/2c34961a51ecc9f12ab2e783136fc864022611d9))


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
