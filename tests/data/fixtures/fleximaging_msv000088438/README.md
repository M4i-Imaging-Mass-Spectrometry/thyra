# FlexImaging sequence files from MassIVE MSV000088438

Two Bruker timsTOF fleX MALDI imaging runs from the public MassIVE dataset
MSV000088438 (TIMSCONVERT raw data, CC0 1.0), folder
`raw/timsconvert_raw_data/ims/`. Only the flexImaging side files are here: the
sequence file (`.mis`) and the spot list flexImaging wrote for it. The `.d`
folders and the photos are not needed by the tests.

| Run | Raster | Spots | Regions |
|---|---|---|---|
| `20210920_vc_rugose_1mMTCA_gordon` (TSF) | 1000 um | 276 | 4 |
| `20210921_vc_rugose_tims_gordon` (TDF) | 1000 um | 240 | 4 |

The spot list gives every spot's stage position in the teaching-point frame.
Mapped through the three teaching points it says where flexImaging put each
spot on its photo. `tests/unit/alignment/test_raster_lattice.py` checks that
Thyra's lattice fit (D31) lands every spot there.

## Provenance

SHA-256 of the files as MassIVE's download endpoint served them on 2026-10-05:

| File | Bytes | SHA-256 |
|---|---|---|
| `20210920_vc_rugose_1mMTCA_gordon.mis` | 4046 | `0afba0dcab61c37773af0a94d94835aa21cc61715c6d5120806cfc1adbd73bfb` |
| `20210920_vc_rugose_1mMTCA_gordon_spot_list.txt` | 9643 | `7600fda92e62ad6193846890f644df3b2563411e68b1cf94f982083a118fb1af` |
| `20210921_vc_rugose_tims_gordon.mis` | 3983 | `b3080f000bffe15c45b0bd29811fb2d0d954aa3df2b3302545d899eabc5df908` |
| `20210921_vc_rugose_tims_gordon_spot_list.txt` | 8542 | `9f5c4375546f014a609ac46c703d94522da3b9c68b7deea70f71053b68081e06` |

Two changes were made, both outside what the tests read:

- The `.mis` files arrived with CRLF line endings and are stored with LF, as
  the repository's line-ending policy requires.
- Absolute Windows paths were cut to their file names: the method in each
  `<Method>`, the `<OriginalImage>`, and the spot list's own path in its
  header line. The repository's path checks reject absolute paths.
