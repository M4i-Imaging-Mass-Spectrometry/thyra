"""Write a small imzML that carries an ion mobility array per spectrum.

pyimzml writes the m/z and intensity arrays; the third array is spliced in
the way the TIMS exporters lay it out, one copy per spectrum unless asked
otherwise. Used by the reader and the converter tests of the shared
(m/z, mobility) axis and of non-finite mobility values.
"""

import re

import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter

_MOBILITY_GROUP = (
    '<referenceableParamGroup id="mobilityArray">'
    '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>'
    '<cvParam cvRef="MS" accession="MS:1003006" '
    'name="mean inverse reduced ion mobility array" value="" unitCvRef="MS" '
    'unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>'
    '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>'
    '<cvParam cvRef="IMS" accession="IMS:1000101" name="external data" '
    'value="true"/>'
    "</referenceableParamGroup>"
)


def _ims_term(accession, name, value):
    return (
        f'<cvParam cvRef="IMS" accession="IMS:1000{accession}" name="{name}" '
        f'value="{value}"/>'
    )


def write_mobility_imzml(
    path,
    mzs,
    mobility,
    *,
    one_mz_block=True,
    one_mobility_block=False,
):
    """An imzML with a mobility array per spectrum, laid out as asked.

    ``mzs`` and ``mobility`` hold one array per spectrum. ``one_mz_block``
    writes the m/z values once for every spectrum, as a continuous file
    does; otherwise each spectrum gets its own block and the file still
    declares continuous mode, as a mislabelled export does.
    ``one_mobility_block`` points every spectrum at the first spectrum's
    mobility block instead of a copy each. Spectra lie on a 12-wide raster.
    """
    n = len(mobility)
    mzs = [np.asarray(m, dtype=np.float64) for m in mzs]
    mode = "continuous" if one_mz_block else "processed"
    with ImzMLWriter(
        str(path), mode=mode, mz_dtype=np.float64, intensity_dtype=np.float64
    ) as writer:
        for i in range(n):
            writer.addSpectrum(
                mzs[i],
                np.arange(1.0, mzs[i].size + 1.0) * (i + 1),
                (i % 12 + 1, i // 12 + 1, 1),
            )
    ibd = bytearray(path.with_suffix(".ibd").read_bytes())
    order = iter(range(n))
    first_offset = []

    def array(_match):
        values = np.asarray(mobility[next(order)], dtype=np.float64)
        if one_mobility_block and first_offset:
            offset = first_offset[0]
        else:
            offset = len(ibd)
            ibd.extend(values.tobytes())
            first_offset.append(offset)
        return (
            '<binaryDataArray encodedLength="0">'
            '<referenceableParamGroupRef ref="mobilityArray"/>'
            + _ims_term(103, "external array length", values.size)
            + _ims_term(104, "external encoded length", values.size * 8)
            + _ims_term(102, "external offset", offset)
            + "</binaryDataArray></binaryDataArrayList>"
        )

    xml = path.read_text("latin-1").replace(
        "</referenceableParamGroupList>",
        _MOBILITY_GROUP + "</referenceableParamGroupList>",
    )
    xml = xml.replace('ParamGroupList count="4"', 'ParamGroupList count="5"')
    xml = xml.replace('ArrayList count="2"', 'ArrayList count="3"')
    if not one_mz_block:
        xml = xml.replace(
            'accession="IMS:1000031" name="processed"',
            'accession="IMS:1000030" name="continuous"',
        )
    path.write_text(re.sub("</binaryDataArrayList>", array, xml), "latin-1")
    path.with_suffix(".ibd").write_bytes(bytes(ibd))
    return path
