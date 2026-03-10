# API Reference

## Main Function

The primary API entry point for converting MSI data:

::: thyra.convert.convert_msi

## Reader Base Class

All format readers inherit from this base class:

::: thyra.core.base_reader.BaseMSIReader
    options:
      members:
        - get_essential_metadata
        - get_comprehensive_metadata
        - get_common_mass_axis
        - get_optical_image_paths
        - iter_spectra
        - get_region_map
        - get_region_info
        - has_shared_mass_axis

## Converter Base Class

All output converters inherit from this base class:

::: thyra.core.base_converter.BaseMSIConverter
    options:
      members:
        - convert
        - pixel_size_um
        - pixel_size_source
        - dataset_id
        - handle_3d
