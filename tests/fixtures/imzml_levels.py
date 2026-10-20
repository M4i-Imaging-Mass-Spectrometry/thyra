"""Write a small imzML whose spectra state MS levels, the two ways files do.

pyimzml's writer states ``ms level`` 0 in the parameter group every spectrum
references, and has no argument for another level. These helpers write the
file with it and then state each spectrum's level one of two ways:

- ``"own"``: a ``MS:1000511`` term on the spectrum itself, the writer's
  group term removed (as the F26 repro of issue #522 does).
- ``"group"``: one parameter group per level, referenced after the writer's
  group, which keeps its level 0.

A level of ``None`` states none for that spectrum.
"""

import re
from pathlib import Path
from typing import Optional, Sequence, Tuple

from pyimzml.ImzMLWriter import ImzMLWriter

#: ``(x, y, z, mzs, intensities)``, 1-based positions as the file states them.
Spectrum = Tuple[int, int, int, Sequence[float], Sequence[float]]

WRITER_LEVEL = '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="0"/>'


def level_term(level: int) -> str:
    """A spectrum's own ``ms level`` term."""
    return (
        f'<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{level}"/>'
    )


def write_imzml_with_levels(
    path: Path,
    spectra: Sequence[Spectrum],
    levels: Sequence[Optional[int]],
    how: str = "own",
    mode: str = "processed",
) -> Path:
    """Write ``spectra`` to ``path`` and state ``levels`` on them ``how``."""
    assert len(levels) == len(spectra)
    with ImzMLWriter(str(path), mode=mode) as writer:
        for x, y, z, mzs, intensities in spectra:
            writer.addSpectrum(list(mzs), list(intensities), (x, y, z))
    text = path.read_text()
    assert WRITER_LEVEL in text
    stated = iter(levels)
    spectrum_start = r"(<spectrum [^>]*>\s*<referenceableParamGroupRef [^>]*/>)"
    if how == "own":
        text = text.replace(WRITER_LEVEL, "")

        def own(match: "re.Match[str]") -> str:
            level = next(stated)
            return match.group(1) + ("" if level is None else level_term(level))

        text = re.sub(spectrum_start, own, text)
    elif how == "group":
        wanted = sorted({level for level in levels if level is not None})
        groups = "".join(
            f'<referenceableParamGroup id="level{level}">{level_term(level)}'
            "</referenceableParamGroup>"
            for level in wanted
        )
        text = re.sub(
            r'<referenceableParamGroupList count="(\d+)">',
            lambda m: f'<referenceableParamGroupList count="{int(m.group(1)) + len(wanted)}">'
            + groups,
            text,
            count=1,
        )

        def grouped(match: "re.Match[str]") -> str:
            level = next(stated)
            if level is None:
                return match.group(1)
            return match.group(1) + f'<referenceableParamGroupRef ref="level{level}"/>'

        text = re.sub(spectrum_start, grouped, text)
    else:
        raise ValueError(how)
    assert next(stated, "done") == "done"
    path.write_text(text)
    return path
