# Bruker SDK DLL Location

This folder contains the Bruker TDF SDK library files used by Thyra to
read Bruker .d datasets (timsTOF TDF and MALDI TSF formats).

## Bundled Binaries

The binaries shipped here come from
[opentims_bruker_bridge](https://github.com/MatteoLacki/opentims_bruker_bridge)
v1.2.0 and support both `tims_*` (TDF/timsTOF) and `tsf_*` (TSF/MALDI)
functions.

- Windows: `timsdata.dll`
- Linux: `timsdata.so` (rename to `libtimsdata.so` if needed)

## License and Attribution

The SDK binaries are proprietary software provided by Bruker Daltonik GmbH
under the Bruker Software License Agreement. See
[LICENCE-BRUKER.txt](LICENCE-BRUKER.txt) in this directory for the full
license text.

> This software uses TDF SDK software.
> Copyright (c) 2019 by Bruker Daltonik GmbH. All rights reserved.

## Replacing the SDK

To use a different version of the SDK, place the appropriate library file
in this directory:

- Windows: `timsdata.dll`
- Linux: `libtimsdata.so`
- macOS: `libtimsdata.dylib` (limited support)

The library will be automatically detected by Thyra.

### Alternative Locations

If you cannot place the library in this folder, you can also:

1. Set the `BRUKER_SDK_PATH` environment variable to point to the library location
2. Place the library in the same directory as your data (.d folder)
3. Place the library in your current working directory
4. Add the SDK directory to your system PATH

However, placing the library in this repository folder is the most reliable method.

### Where to Get the SDK

1. [opentims_bruker_bridge](https://github.com/MatteoLacki/opentims_bruker_bridge) (bundled here)
2. Bruker Daltonics official channels
3. The Bruker timsTOF software installation
4. Contact your Bruker representative
