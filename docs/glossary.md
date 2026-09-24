# Glossary

Words you will meet in these pages, in plain language. The computer terms
come first, then the mass spectrometry ones.

## Computer terms

**Terminal**
: A window where you type commands. On Windows it is called PowerShell; on
  macOS and Linux, Terminal.

**Command**
: A line you type into the terminal and run with Enter, such as
  `thyra --version`.

**Path**
: The address of a file or folder, such as `/Volumes/data/MSI/brain.d` on
  macOS. On Windows, a path starts with a drive letter, such as `D:`.

**Python**
: The programming language Thyra is written in. You need it installed, but
  you do not need to write any.

**pip**
: Python's standard installer for packages such as Thyra.

**uv**
: A fast installer for Python programs. It can also fetch the right version
  of Python for you.

**Environment**
: A separate set of Python packages, so that the packages of one project do
  not interfere with another's. Tools such as conda and uv create them.

**Jupyter notebook**
: A document that mixes text with Python code you can run piece by piece. The
  [Explore the result in Python](explore-the-output.ipynb) page is one.

**SpatialData**
: An open format and Python library for spatial data from biology, such as
  imaging mass spectrometry. Thyra saves your data in this format.

**Zarr**
: The way SpatialData stores data on disk: a folder of many small files that
  programs can read a piece at a time. That is why the result is a folder
  ending in `.zarr`.

**Table (AnnData)**
: The part of the result that holds the spectra: one row per pixel, one
  column per m/z bin. AnnData is the name of the Python format for such
  tables.

**scverse**
: A community of open-source Python tools for single-cell and spatial
  biology. SpatialData, napari-spatialdata, scanpy and squidpy belong to it.

**napari**
: A free image viewer. With the napari-spatialdata plugin, it opens
  SpatialData results.

**Command-line option**
: An extra instruction added to a command, starting with `--`, such as
  `--pixel-size 50`.

## Mass spectrometry terms

**MSI (mass spectrometry imaging)**
: Recording a mass spectrum at every point of a sample, so that each m/z
  value can be shown as an image.

**m/z**
: Mass-to-charge ratio, the x-axis of a mass spectrum.

**Spectrum**
: The intensities measured at one pixel, across the m/z range.

**Pixel size**
: The distance between neighbouring measurement points, in micrometres.

**TIC (total ion current)**
: The sum of all intensities in one spectrum. A TIC image shows it for every
  pixel, and gives a quick overview of the sample.

**Ion image**
: An image of the intensity of one m/z value, or a narrow range around it,
  across all pixels.

**Profile and centroid data**
: Profile data records the full shape of each peak. Centroid data keeps one
  point per peak.

**Resampling and m/z bins**
: Instruments record each spectrum at slightly different m/z values.
  Resampling moves every spectrum onto one shared set of m/z positions,
  called bins, so that the same column means the same m/z in every pixel.

**Region**
: A separately defined area of an acquisition, for example one tissue
  section on a slide that holds several.

**Optical image**
: A photo or microscope scan of the slide, taken before the measurement.

**Teaching points**
: Marks that link positions in the optical image to positions of the
  instrument's stage. Thyra uses them to lay the MSI data over the photo.

**imzML**
: An open file format for imaging mass spectrometry. It comes as two files: an
  `.imzML` file that describes the data and an `.ibd` file that holds the
  numbers.

**Ion mobility (TIMS)**
: A separation of ions by their size and shape, in addition to m/z. Bruker
  timsTOF instruments measure it with trapped ion mobility spectrometry
  (TIMS).

**MS/MS**
: Tandem mass spectrometry: selected ions are broken into fragments, and the
  fragments are measured.
