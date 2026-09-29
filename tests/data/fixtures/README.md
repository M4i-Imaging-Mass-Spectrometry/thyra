# Hand-authored imzML corpus

Eight tiny imzML/`.ibd` pairs whose XML was written out literally, character
by character, and whose binary was packed by hand.

Everything else the test suite reads was produced by pyimzml's own
`ImzMLWriter`. A parser and a writer from the same codebase agree on each
other's mistakes, so that corpus can only ever exercise the shapes the writer
happens to emit — and its emission profile diverges from real vendor files on
nine structural points, every one of which is unreachable from a generated
fixture. These files close some of that gap.

They are small on purpose: 60,885 bytes of imzML and `.ibd` checked out, and
9,211 bytes zlib-compressed — the XML is highly repetitive, so what the
repository actually carries is under 10 KB.

## The files

| Pair | The one thing it carries |
|---|---|
| `iontof_sparse` | ISO-8859-1, CRLF, unindented, no `mzML/@version`, no `IMS:1000052`, three misspelled ontology names, `max dimension` contradicting `pixel size`, and a sparse 6-of-16 acquisition. The whole IONTOF/bellini class in one file. |
| `unit_nanometre` | `IMS:1000046/47` carrying `UO:0000018` nanometre. pyimzml's `imzmldict` drops the unit, so a 4.40625 um pixel would read as 4406.25 um if the extractor did not fetch the unit from the ParamGroup path. |
| `two_scansettings` | Two `<scanSettings>` blocks with different pixel sizes. `__readimzmlmeta` resolves each accession by first-match-anywhere, producing a chimera present in neither block. |
| `two_precision_terms` | One param group declaring both `32-bit float` and `64-bit float`. pyimzml picks by dict insertion order and says nothing. |
| `solarix_fticr` | A profile solariX MRMS export: model term `MS:1001549` inherited through a `CommonInstrumentParams` referenceableParamGroup, `MS:1000079` FT-ICR analyzer on the componentList. The shape the imzML extractor must resolve to `instrument_type = "FT-ICR"` for `FTICRDetector` to ever fire. |
| `timstof_flex_export` | A profile timsTOF fleX export: model term `MS:1003124` directly on the instrumentConfiguration. Outside the native `.d` there is no `GlobalMetadata`, so the surfaced model name is the only way `is_timstof` can recognise the instrument. |
| `mobility_continuous` | A TIMSImaging-style export with a **third binary array**: `mobilityArray` bound to `MS:1003006` (unit `MS:1002814`), continuous mode, one shared feature list with m/z repeated where mobility splits a feature, mobility copied per spectrum as the pyimzML fork writes it. Upstream pyimzml ignores the array; the shared m/z block is not strictly increasing. |
| `mobility_processed` | The same convention in processed mode, TIMSCONVERT-style: a per-pixel (m/z, mobility) point cloud with one m/z at two mobilities. No shared feature axis; the repeats must sum into one bin. |

Deliberately **not** here: zlib compression and continuous mode. Both are
reachable through `ImzMLWriter` — `ZlibCompression()` and
`generate_example_imzml`, which already writes `mode="continuous"` — so a
hand-authored pair would buy nothing.

## Regenerating

```bash
PYTHONPATH=. python tests/data/fixtures/build_fixtures.py
```

Reproducible byte for byte, so a clean run leaves `git status` clean. No test
invokes it; the committed bytes are the fixture. The script is provenance —
read it to find out why a given cvParam is worded the way it is.

## `synthetic_tims.d`: a hand-written Bruker TDF acquisition

Six pixels on a 3 x 2 raster, 240 TIMS scans per frame, three planted ions
per pixel with a Gaussian mobility profile plus seeded noise. Everything in it
is invented; the schema, block layout and calibration model types are those
of a real TDF 3.7 file so that Bruker's own `timsdata` library opens it, reads
every scan, converts indices to m/z and scans to 1/K0, and runs its centroid
extraction on it. `synthetic_tims_expected.json` records what was written,
pair by pair, so `tests/integration/test_bruker_tdf_synthetic.py` can check
the reader against ground truth through the real SDK.

It also carries the calibration record a real MALDI acquisition does, for the
metadata tests: a `CalibrationInfo` table in `analysis.tdf` (the external
calibration the run started with, whose stated standard deviation follows from
its own arrays), and a `calibration.sqlite` beside it holding the online
lock-mass state timsControl writes as a run starts, placeholders included
(`MzStandardDeviationPPM` 0.000000 against one reference peak, a mobility
calibration dated 1999-02-02). The state's per-frame calibrators repeat the
analysis database's coefficients, so every m/z the SDK reports is the same
with the file as without it.

Three more places read it, which is why corrupting it is expensive:
`tests/unit/readers/test_bruker_fragmentation.py` and
`test_bruker_precursor_spectra.py` copy the directory and rewrite the copy's
SQLite to synthesise MS/MS acquisitions, and `docs/explore-the-output.ipynb`
copies it into `example_data/` to build the notebook's PASEF demo. Only the
integration module skips when it cannot open the fixture; the two unit modules
drive the reader with a fake SDK and fail outright.

Built by `build_tdf_fixture.py`, which needs `zstandard`. That is declared in
the `test` dependency group, so `uv sync --group test` installs it — plain
`uv sync`, which is what `docs/contributing.md` tells a contributor to run,
does not. The module docstring documents the frame-block layout and the SDK's
intensity scaling that the fixture sidesteps by declaring a 100 ms
accumulation time. All three files are marked binary by `.gitattributes:39`,
the same guarantee `.gitattributes:15` gives the imzML corpus.

## Pixel size test files

Five more pairs, `pixel_size_*.imzML` and `.ibd`. Each holds the same
synthetic acquisition: 3 by 2 pixels, each 50 um wide, six profile spectra on
one shared m/z axis, declared as a meandering scan. Only the way the header
states the pixel size differs. The five ways are the ones found in the
headers of public imzML files.

| Pair | What the header gives | True pixel size |
|---|---|---|
| `pixel_size_unit_declared` | `IMS:1000046` and `IMS:1000047`, value 50, unit micrometre, and the extent | 50 um |
| `pixel_size_unit_absent` | both terms, value 50, no unit, and the extent | 50 um, the unit is not stated |
| `pixel_size_area_old_name` | `IMS:1000046` alone, named "pixel size", value 2500, no unit, and the extent | 50 um, stated as the area of a pixel |
| `pixel_size_area_old_name_no_extent` | the same without the extent | 50 um, and the file cannot show it |
| `pixel_size_unit_contradiction` | `IMS:1000046` alone, named "pixel size", value 50, unit accession `UO:0000015` (centimetre), unit name "micrometer", no extent | 50 um, and the file cannot show it |

The extent is `IMS:1000044` and `IMS:1000045`, 150 um by 100 um. Until 2017
`IMS:1000046` was named "pixel size" and gave the area of a pixel, see
[imzML Parser Notes](../../../docs/imzml-parser-notes.md#ims1000046-was-an-area-until-2017).

`pixel_size_expected.json` records for each file what a reader or a converter
can conclude from the header alone. `pixel_size_um` is 50 on both axes or
null. `pixel_size_source` says how: `declared`, `unit_assumed` (the value is
there, micrometre is assumed), `derived_from_area` (the square root of the
value, because the square root times the pixel count is the extent), or
`unknown` (the header does not settle it, a person has to give the size).

What Thyra takes from each file is held by
`tests/unit/metadata/extractors/test_pixel_size_fixtures.py`:

| Pair | Read as imzML | Read from an mzPeak archive with the same terms |
|---|---|---|
| `pixel_size_unit_declared` | 50 um | 50 um |
| `pixel_size_unit_absent` | 50 um | 50 um |
| `pixel_size_area_old_name` | none, asks for `--pixel-size` | 50 um, logged as an area |
| `pixel_size_area_old_name_no_extent` | none, asks for `--pixel-size` | none, asks for `--pixel-size` |
| `pixel_size_unit_contradiction` | none, asks for `--pixel-size` | none |

These files are written by pyimzml's `ImzMLWriter`, unlike the corpus above.
The subject here is the geometry terms, so the rest of the file is as plain
as a file can be. The run inside each file is named without the prefix, for
example `area_old_name`. Every line ends in CRLF.

The files are synthetic and carry no measured data. They can be used in the
test suite of any other reader or converter, under the licence of this
repository.

Built by `build_pixel_size_fixtures.py`:

```bash
python tests/data/fixtures/build_pixel_size_fixtures.py
```

Reproducible byte for byte: the script fixes the UUID of each file and sets
the line ends. The same guards as for the corpus above hold the committed
bytes, in `TestCommittedBytes`.

## Three ways these files get destroyed

The first two leave no mark on the worktree file that a reader would notice,
and before `TestCommittedBytes` existed neither failed a single test. The
third is loud instead — and it is the one that actually happened.

1. **Staging.** `.gitattributes:15` is `* text=auto eol=lf`. Without the
   `tests/data/fixtures/*.imzML -text` exemption, staging `iontof_sparse.imzML`
   rewrites all 174 CRLFs to LF — 11,593 bytes become 11,419 — and the CRLF
   property the fixture exists to prove is gone from the repository while the
   worktree still looks right. Only the blob-reading tests can see this;
   measured, removing the exemption fails ten of them and no other test in the
   suite.
2. **Pre-commit.** `mixed-line-ending --fix=lf` does the same thing to the
   worktree file, and `trailing-whitespace`, `end-of-file-fixer` and
   `fix-byte-order-marker` are free to move any byte. All four rewrite in
   place, so all four carry
   `exclude: ^tests/data/fixtures/([^/]+\.(imzML|ibd)|synthetic_tims\.d/.*)$`,
   scoped to the byte-exact fixtures so this README and `build_fixtures.py`
   stay covered. `.pre-commit-config.yaml`'s header names the same four; add a
   fifth in-place rewriter and it needs the same `exclude`.
3. **The build script.** `build_tdf_fixture.py` used to `import zstandard`
   inside `encode_frame`, which the build reaches only after
   `shutil.rmtree(OUT_DIR)` has run and the new `analysis.tdf_bin` is open.
   On a machine without the module — since `xarray-spatial` moved it behind an
   extra, that is any freshly resolved environment — running the script deleted
   `analysis.tdf` and left `analysis.tdf_bin` at its 64-byte header. The 19
   tests that then fail say `no such table: Frames`, naming neither the script
   nor the missing dependency. Closed twice over: the import is at module
   scope, so it fails before a byte is touched, and the build now writes into a
   staging directory that is renamed over the fixture only once complete.

`.gitignore` also needs the negations: `*.imzML` and `*.ibd` ignore every file
with those extensions, so without `!tests/data/fixtures/*.imzML` these files
cannot be added at all. The patterns carry the dot on purpose: the bare
`*imzML` that preceded them also matched the *directory* `thyra/readers/imzml`
on case-insensitive checkouts, which silently kept new modules there out of
every commit.

**If you change a fixture, verify the blob, not the file on disk:**

```bash
git cat-file -p :tests/data/fixtures/iontof_sparse.imzML > /tmp/blob
python -c "d=open('/tmp/blob','rb').read(); print(len(d), d.count(b'\r\n'))"
# 11593 174
```

`test_hand_authored_fixtures.py::TestCommittedBytes` asserts exactly this
against `git cat-file`, so the guard runs in CI as well — but only as long as
the test itself keeps reading the blob.

## Characterisation, not endorsement

Several tests over these fixtures assert behaviour that is currently **wrong**:
a discarded unit, a chimeric metadata dict, an ambiguous precision declaration
that is silently accepted. Every such assertion carries a comment naming the
audit finding it pins and what a fix would change it to. Do not read
`assert ... == 4406.25` as a statement that 4406.25 is correct.
