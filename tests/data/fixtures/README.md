# Hand-authored imzML corpus

Four tiny imzML/`.ibd` pairs whose XML was written out literally, character by
character, and whose binary was packed by hand.

Everything else the test suite reads was produced by pyimzml's own
`ImzMLWriter`. A parser and a writer from the same codebase agree on each
other's mistakes, so that corpus can only ever exercise the shapes the writer
happens to emit — and its emission profile diverges from real vendor files on
nine structural points, every one of which is unreachable from a generated
fixture. These four files close some of that gap.

They are small on purpose: 41,208 bytes of imzML and `.ibd` checked out, and
6,051 bytes zlib-compressed — the XML is highly repetitive, so what the
repository actually carries is under 6 KB.

## The files

| Pair | The one thing it carries |
|---|---|
| `iontof_sparse` | ISO-8859-1, CRLF, unindented, no `mzML/@version`, no `IMS:1000052`, three misspelled ontology names, `max dimension` contradicting `pixel size`, and a sparse 6-of-16 acquisition. The whole IONTOF/bellini class in one file. |
| `unit_nanometre` | `IMS:1000046/47` carrying `UO:0000018` nanometre. pyimzml's `imzmldict` drops the unit, so a 4.40625 um pixel would read as 4406.25 um if the extractor did not fetch the unit from the ParamGroup path. |
| `two_scansettings` | Two `<scanSettings>` blocks with different pixel sizes. `__readimzmlmeta` resolves each accession by first-match-anywhere, producing a chimera present in neither block. |
| `two_precision_terms` | One param group declaring both `32-bit float` and `64-bit float`. pyimzml picks by dict insertion order and says nothing. |

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

## Two ways these files get destroyed

Neither leaves a mark on the worktree file that a reader would notice, and
before `TestCommittedBytes` existed neither failed a single test.

1. **Staging.** `.gitattributes:15` is `* text=auto eol=lf`. Without the
   `tests/data/fixtures/*.imzML -text` exemption, staging `iontof_sparse.imzML`
   rewrites all 174 CRLFs to LF — 11,593 bytes become 11,419 — and the CRLF
   property the fixture exists to prove is gone from the repository while the
   worktree still looks right. Only the blob-reading tests can see this;
   measured, removing the exemption fails ten of them and no other test in the
   suite.
2. **Pre-commit.** `mixed-line-ending --fix=lf` does the same thing to the
   worktree file, and `trailing-whitespace` and `end-of-file-fixer` are free to
   move any byte. All three carry
   `exclude: ^tests/data/fixtures/[^/]+\.(imzML|ibd)$`, scoped to the two
   byte-exact extensions so this README and `build_fixtures.py` stay covered.

`.gitignore` also needs the negations: `*imzML` and `*ibd` are unanchored, so
without `!tests/data/fixtures/*.imzML` these files cannot be added at all.

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
