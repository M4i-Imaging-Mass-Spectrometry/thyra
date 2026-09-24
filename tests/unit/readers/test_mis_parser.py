"""Tests for the .mis file parser and discovery helpers."""

import logging
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from thyra.errors import ConversionRefused
from thyra.readers.bruker.mis_parser import (
    _extract_areas,
    find_mis_file_for_d_folder,
    parse_mis_file,
)
from thyra.readers.bruker.timstof.timstof_reader import BrukerReader

RECTANGULAR_AREA = (
    '<Area Type="0" Name="01"><Point>10,20</Point><Point>30,40</Point></Area>'
)

# A real FlexImaging polygon, ten points, from the acquisition that produced
# issue #84: its bounding box corners are on different vertices from the first
# two points, so a parser that reads only two points gets it wrong.
POLYGON_AREA = """<Area Type="3" Name="01">
    <Point>24470,4585</Point>
    <Point>24420,3818</Point>
    <Point>24862,3543</Point>
    <Point>25353,3168</Point>
    <Point>26228,3043</Point>
    <Point>26753,4193</Point>
    <Point>26462,5485</Point>
    <Point>25362,5777</Point>
    <Point>24737,4943</Point>
    <Point>24487,4552</Point>
</Area>"""


def _write_mis(
    tmp_path: Path,
    name: str,
    raster: str = "5,5",
    area: str = RECTANGULAR_AREA,
) -> Path:
    mis = tmp_path / name
    mis.write_text(f"""<?xml version="1.0"?>
<ImagingSequence>
<ImageFile>img.tif</ImageFile>
<Raster>{raster}</Raster>
{area}
</ImagingSequence>
""")
    return mis


def test_parse_mis_extracts_raster(tmp_path: Path) -> None:
    mis = _write_mis(tmp_path, "sample.mis", raster="5,5")
    data = parse_mis_file(mis)
    assert data["raster"] == [5, 5]


def test_parse_mis_handles_rectangular_raster(tmp_path: Path) -> None:
    mis = _write_mis(tmp_path, "sample.mis", raster="10,20")
    data = parse_mis_file(mis)
    assert data["raster"] == [10, 20]


def test_find_mis_prefers_matching_stem(tmp_path: Path) -> None:
    """When several .mis files sit in the parent, prefer the one whose stem
    matches the .d folder stem.
    """
    d_folder = tmp_path / "sample_A.d"
    d_folder.mkdir()
    _write_mis(tmp_path, "sample_A.mis", raster="5,5")
    _write_mis(tmp_path, "sample_B.mis", raster="50,50")

    found = find_mis_file_for_d_folder(d_folder)
    assert found is not None
    assert found.name == "sample_A.mis"


def test_find_mis_falls_back_to_any(tmp_path: Path) -> None:
    d_folder = tmp_path / "sample.d"
    d_folder.mkdir()
    other = _write_mis(tmp_path, "different_name.mis", raster="7,7")

    found = find_mis_file_for_d_folder(d_folder)
    assert found == other


def test_find_mis_returns_none_when_missing(tmp_path: Path) -> None:
    d_folder = tmp_path / "sample.d"
    d_folder.mkdir()
    assert find_mis_file_for_d_folder(d_folder) is None


def test_extract_areas_rectangular() -> None:
    """Area extraction with rectangular (Type=0) 2-point areas.

    Moved here from test_rapiflex_reader.py, which called the same assertions
    against RapiflexReader._extract_areas. It passes on both sides of the
    parser merge: the shared parser already carried this logic. What it guards
    is the coverage, not the merge -- it is now asserted on the one code path
    Rapiflex, timsTOF, solariX and BrukerMetadataExtractor share.
    """
    root = ET.fromstring(f"<Root>{RECTANGULAR_AREA}</Root>")
    metadata: dict = {}
    _extract_areas(root, metadata)

    assert len(metadata["areas"]) == 1
    area = metadata["areas"][0]
    assert area["name"] == "01"
    assert area["p1"] == [10, 20]
    assert area["p2"] == [30, 40]


def test_extract_areas_polygon() -> None:
    """Area extraction with polygon (Type=3) N-point areas.

    The regression test for issue #84 / PR #85, which fixed a bounding box
    computed from the first two <Point> elements only. It was written against
    the Rapiflex copy of the parser and stayed there, so the shared parser --
    the one three readers and the metadata extractor use -- carried the fix
    with no test on it. Like the rectangular case it passes before and after
    the merge; it moves so the fix is covered where the code now lives.
    """
    root = ET.fromstring(f"<Root>{POLYGON_AREA}</Root>")
    metadata: dict = {}
    _extract_areas(root, metadata)

    assert len(metadata["areas"]) == 1
    area = metadata["areas"][0]
    assert area["name"] == "01"
    # Bounding box should span ALL points, not just the first two.
    assert area["p1"] == [24420, 3043]
    assert area["p2"] == [26753, 5777]


def test_parse_mis_extracts_polygon_area_bounding_box(tmp_path: Path) -> None:
    """The polygon bounding box survives the public entry point too.

    test_extract_areas_polygon calls the private helper with a tree built in
    memory. This one goes through parse_mis_file from a file on disk, so the
    XML parser, the .//Area search and the helper are exercised together.
    Passes before and after the merge, like the two above it.
    """
    mis = _write_mis(tmp_path, "polygon.mis", area=POLYGON_AREA)
    data = parse_mis_file(mis)

    assert data["areas"] == [{"name": "01", "p1": [24420, 3043], "p2": [26753, 5777]}]


def _write_entity_mis(tmp_path: Path, name: str = "entity.mis") -> Path:
    mis = tmp_path / name
    mis.write_text("""<?xml version="1.0"?>
<!DOCTYPE ImagingSequence [<!ENTITY r "5,5">]>
<ImagingSequence><Raster>&r;</Raster></ImagingSequence>
""")
    return mis


def test_entity_bearing_mis_does_not_expand(tmp_path: Path) -> None:
    """An XML entity in a .mis is refused, not expanded.

    The only test in this file that fails before the defusedxml swap: on the
    unguarded `xml.etree` import the same file parsed happily and yielded
    ``{"raster": [5, 5]}``. defusedxml raises EntitiesForbidden, which is a
    ValueError rather than a ParseError, so the except clause has to name it
    or the exception escapes parse_mis_file into three callers that do not
    catch it.

    No importorskip on defusedxml: it is a hard dependency, and an
    importorskip would turn the one security test in this file into a
    silent pass on exactly the install where the hole is open.
    """
    with pytest.raises(ConversionRefused):
        parse_mis_file(_write_entity_mis(tmp_path))


def test_the_refusal_names_the_file_and_the_reason(tmp_path: Path) -> None:
    """A security refusal must not read like an empty file.

    This started as a warning and an empty dict, which none of the four
    consumers (the Rapiflex, timsTOF and solariX readers, and
    BrukerMetadataExtractor) checks for. The acquisition then ran with no
    areas, no teaching points and no raster step, and the first visible
    symptom was a later ``--region <name>`` failing as "no such region" --
    a message about a region list, pointing away from the file that
    emptied it. So the two things the message has to carry are which file
    and why.
    """
    mis = _write_entity_mis(tmp_path, "brain_section.mis")

    with pytest.raises(ConversionRefused) as excinfo:
        parse_mis_file(mis)

    message = str(excinfo.value)
    assert "brain_section.mis" in message
    assert "entit" in message.lower()


def test_malformed_xml_is_still_only_a_warning(tmp_path: Path, thyra_logs) -> None:
    """Not well-formed is not the same claim as refused.

    A truncated or corrupt .mis costs the optical alignment, which every
    caller already treats as optional -- most acquisitions have no .mis at
    all -- so it stays a warning and an empty result, exactly as before.
    The refusal above is the case where Thyra decided not to read a file it
    could have read, and that decision is the one that has to be audible.

    Not caplog: setup_logging sets propagate=False on the `thyra` logger
    process-globally, so a caplog assertion here would pass alone and fail
    after any test that has invoked the CLI. See the thyra_logs fixture.
    """
    mis = tmp_path / "truncated.mis"
    mis.write_text('<?xml version="1.0"?>\n<ImagingSequence><Raster>5,5')

    with thyra_logs("thyra.readers.bruker.mis_parser", logging.WARNING) as records:
        data = parse_mis_file(mis)

    assert data == {}
    assert any("Failed to parse .mis file" in r.getMessage() for r in records)


def test_the_timstof_reader_does_not_swallow_the_refusal(tmp_path: Path) -> None:
    """The one consumer with a ``except ValueError`` anywhere near it.

    ``BrukerReader._parse_mis_alignment`` wraps
    ``get_teaching_points_file()`` in ``except (ValueError, OSError):
    return {}`` for non-standard folder layouts. ``ConversionRefused`` is a
    ``ValueError``, so widening that try by two lines to cover the
    ``parse_mis_file`` call under it would restore the silence with no
    other visible change -- and this is the consumer where the silence
    hurt most, because ``_parse_mis_alignment`` runs in ``__init__``
    before ``_select_region``, and the areas it fills are what resolves
    ``--region <name>``.

    Called on an uninitialised instance: the method reads nothing but
    ``get_teaching_points_file``, and constructing the reader properly
    would need the vendor library.
    """
    reader = BrukerReader.__new__(BrukerReader)
    mis = _write_entity_mis(tmp_path)
    reader.get_teaching_points_file = lambda: mis  # type: ignore[method-assign]

    with pytest.raises(ConversionRefused, match="entity.mis"):
        reader._parse_mis_alignment()


def test_internal_subset_dtd_still_parses(tmp_path: Path) -> None:
    """A .mis carrying a plain DTD is not collateral damage of the above.

    defusedxml's parse() defaults are forbid_dtd=False, forbid_entities=True,
    forbid_external=True. Only the entity declaration is refused; a document
    type declaration on its own still reads. Passes before and after, and is
    here so a later tightening to forbid_dtd=True cannot pass unnoticed.
    """
    mis = tmp_path / "dtd.mis"
    mis.write_text("""<?xml version="1.0"?>
<!DOCTYPE ImagingSequence [<!ELEMENT Raster (#PCDATA)>]>
<ImagingSequence><Raster>5,5</Raster></ImagingSequence>
""")

    assert parse_mis_file(mis)["raster"] == [5, 5]


class TestTheMisPickIsDeterministic:
    """Issue #303: four locators, three of them order-dependent.

    ``list(glob("*.mis"))[0]`` made the pick depend on directory listing
    order, and the pick decides the pixel pitch -- ``_resolve_pixel_size_um``
    prefers the .mis ``<Raster>`` over ``BeamScanSize`` -- and the
    acquisition areas ``--region`` resolves against. The same dataset
    could convert two ways on two machines.

    The locators also searched different directories, which needed no
    unlucky ordering at all: a .mis inside the .d was visible to the
    areas lookup and invisible to the pitch lookup.
    """

    def test_the_pick_does_not_depend_on_listing_order(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        d_folder = tmp_path / "sample.d"
        d_folder.mkdir()
        _write_mis(tmp_path, "a_other.mis", raster="5,5")
        _write_mis(tmp_path, "z_other.mis", raster="50,50")
        # Two non-matching candidates: refused rather than guessed at.
        assert find_mis_file_for_d_folder(d_folder) is None

        _write_mis(tmp_path, "sample.mis", raster="9,9")
        first = find_mis_file_for_d_folder(d_folder)

        real_glob = Path.glob
        monkeypatch.setattr(
            Path, "glob", lambda self, p: reversed(list(real_glob(self, p)))
        )
        assert find_mis_file_for_d_folder(d_folder) == first
        assert first is not None and first.name == "sample.mis"

    def test_several_non_matching_candidates_are_refused_not_guessed(
        self, tmp_path: Path, thyra_logs
    ) -> None:
        """An arbitrary pick writes a *wrong* pitch into the store.

        ``None`` falls back to ``BeamScanSize``, which is at least the
        instrument's own answer. This is what the solariX reader already
        did.
        """
        d_folder = tmp_path / "sample.d"
        d_folder.mkdir()
        _write_mis(tmp_path, "other_a.mis", raster="5,5")
        _write_mis(tmp_path, "other_b.mis", raster="50,50")

        with thyra_logs("thyra.readers.bruker.mis_parser", logging.WARNING) as records:
            assert find_mis_file_for_d_folder(d_folder) is None
        assert any("refusing to guess" in r.getMessage() for r in records)

    def test_a_lone_non_matching_candidate_is_still_accepted(
        self, tmp_path: Path
    ) -> None:
        """There is nothing to guess between."""
        d_folder = tmp_path / "sample.d"
        d_folder.mkdir()
        other = _write_mis(tmp_path, "different_name.mis", raster="7,7")
        assert find_mis_file_for_d_folder(d_folder) == other

    def test_a_mis_inside_the_d_is_found(self, tmp_path: Path) -> None:
        """Rapiflex writes it there; the pitch lookup could not see it."""
        d_folder = tmp_path / "sample.d"
        d_folder.mkdir()
        inside = _write_mis(d_folder, "sample.mis", raster="5,5")
        assert find_mis_file_for_d_folder(d_folder) == inside

    def test_the_pitch_and_the_areas_resolve_to_the_same_file(
        self, tmp_path: Path
    ) -> None:
        """The two lookups disagreed by construction: one searched the
        data folder first, the other searched only the parent."""
        from thyra.readers.bruker.folder_structure import BrukerFolderStructure

        d_folder = tmp_path / "sample.d"
        d_folder.mkdir()
        (d_folder / "analysis.tsf").write_text("")
        _write_mis(d_folder, "inside.mis", raster="5,5")
        _write_mis(tmp_path, "outside.mis", raster="50,50")

        pitch_pick = find_mis_file_for_d_folder(d_folder)
        areas_pick = BrukerFolderStructure(d_folder)._find_teaching_points_file(
            d_folder
        )
        assert pitch_pick == areas_pick


class TestTheOtherBrukerPicksAreSorted:
    """The same defect in the non-.mis picks of the same folder."""

    def test_which_d_folder_is_converted_does_not_depend_on_order(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        """This pick decides which dataset gets converted, and what
        format it is called."""
        from thyra.readers.bruker.folder_structure import BrukerFolderStructure

        for name in ("b_second.d", "a_first.d"):
            d = tmp_path / name
            d.mkdir()
            (d / "analysis.tsf").write_text("")

        first = BrukerFolderStructure(tmp_path)._detect_format()[1]

        real_glob = Path.glob
        monkeypatch.setattr(
            Path, "glob", lambda self, p: reversed(list(real_glob(self, p)))
        )
        assert BrukerFolderStructure(tmp_path)._detect_format()[1] == first
        assert first.name == "a_first.d"

    def test_strays_in_one_path_do_not_veto_the_stem_match_in_another(
        self, tmp_path: Path
    ) -> None:
        """The stem match is the pick this function is named for.

        Refusing path-by-path let two unrelated .mis inside the .d end
        the search before ``sample.mis`` next to ``sample.d`` was ever
        looked at -- a worse answer than the arbitrary pick it replaced.
        """
        d_folder = tmp_path / "sample.d"
        d_folder.mkdir()
        _write_mis(d_folder, "stray_a.mis", raster="5,5")
        _write_mis(d_folder, "stray_b.mis", raster="6,6")
        expected = _write_mis(tmp_path, "sample.mis", raster="9,9")

        assert find_mis_file_for_d_folder(d_folder) == expected

    def test_the_rapiflex_pick_does_not_reach_into_the_parent(
        self, tmp_path: Path
    ) -> None:
        """Rapiflex writes its .mis inside the data folder, and its pick
        only ever looked there. Adopting a lone slide-level one would put
        a foreign acquisition's teaching points and raster step into this
        store.

        Driven through ``_find_data_files``, not by calling the locator
        with search paths the test supplies itself -- that version passes
        whatever ``rapiflex_reader`` does, including deleting it.
        """
        from thyra.readers.bruker.rapiflex.rapiflex_reader import RapiflexReader

        folder = tmp_path / "run"
        folder.mkdir()
        (folder / "run.dat").write_bytes(b"")
        (folder / "run_info.txt").write_text("")
        (folder / "run_poslog.txt").write_text("")
        _write_mis(tmp_path, "whole_slide.mis", raster="50,50")

        reader = RapiflexReader.__new__(RapiflexReader)
        reader.data_path = folder
        reader._mis_path = None
        reader._find_data_files()

        assert reader._mis_path is None

    def test_the_rapiflex_pick_still_finds_its_own_mis(self, tmp_path: Path) -> None:
        """The other half: the .mis it does own is still picked up."""
        from thyra.readers.bruker.rapiflex.rapiflex_reader import RapiflexReader

        folder = tmp_path / "run"
        folder.mkdir()
        (folder / "run.dat").write_bytes(b"")
        (folder / "run_info.txt").write_text("")
        (folder / "run_poslog.txt").write_text("")
        mine = _write_mis(folder, "run.mis", raster="5,5")
        _write_mis(tmp_path, "whole_slide.mis", raster="50,50")

        reader = RapiflexReader.__new__(RapiflexReader)
        reader.data_path = folder
        reader._mis_path = None
        reader._find_data_files()

        assert reader._mis_path == mine


class TestTheTimstofSequenceReachesTheDocument:
    """The tsf/tdf reader has always parsed the .mis; its metadata now keeps it.

    The rapiflex and solariX extractors hand the parse on as
    ``raw_metadata["mis_metadata"]``; the tsf/tdf one dropped it, so the
    registration of every timsTOF raster onto its optical image reached no
    store and no document. Built on the synthetic TDF, whose tables have the
    real layout, with a sequence file beside it the way flexImaging writes
    one: the three teaching points and the image name of a real one.
    """

    _TEACH_POINTS = (
        "<TeachPoint>4780,784;-26352,26386</TeachPoint>"
        "<TeachPoint>13648,11296;-8793,5388</TeachPoint>"
        "<TeachPoint>32156,724;28003,26505</TeachPoint>"
    )

    def _acquisition_with_a_sequence(self, tmp_path: Path) -> Path:
        import shutil

        fixture = (
            Path(__file__).resolve().parents[2]
            / "data"
            / "fixtures"
            / "synthetic_tims.d"
        )
        d_folder = tmp_path / "synthetic_tims.d"
        shutil.copytree(fixture, d_folder)
        (tmp_path / "synthetic_tims.mis").write_text(
            '<?xml version="1.0"?>\n<ImagingSequence>'
            "<ImageFile>slide_0000.tif</ImageFile>"
            f"{self._TEACH_POINTS}</ImagingSequence>\n"
        )
        return d_folder

    def test_the_document_states_the_image_and_its_teaching_points(
        self, tmp_path: Path
    ) -> None:
        from thyra.metadata.document import read_metadata_document

        document = read_metadata_document(self._acquisition_with_a_sequence(tmp_path))

        alignment = document["alignment"]
        assert alignment["optical_image_file"] == "slide_0000.tif"
        assert alignment["method"] == "teaching points"
        assert [
            (p["image_x_px"], p["image_y_px"], p["stage_x_um"], p["stage_y_um"])
            for p in alignment["teaching_points"]
        ] == [
            (4780.0, 784.0, -26352.0, 26386.0),
            (13648.0, 11296.0, -8793.0, 5388.0),
            (32156.0, 724.0, 28003.0, 26505.0),
        ]

    def test_the_raw_metadata_keeps_the_parse(self, tmp_path: Path) -> None:
        d_folder = self._acquisition_with_a_sequence(tmp_path)
        with BrukerReader(d_folder, metadata_only=True) as reader:
            raw = reader.get_comprehensive_metadata().raw_metadata

        assert raw["mis_metadata"]["ImageFile"] == "slide_0000.tif"
        assert len(raw["mis_metadata"]["teaching_points"]) == 3

    def test_without_a_sequence_there_is_no_section(self, tmp_path: Path) -> None:
        from thyra.metadata.document import read_metadata_document

        d_folder = self._acquisition_with_a_sequence(tmp_path)
        (tmp_path / "synthetic_tims.mis").unlink()

        assert "alignment" not in read_metadata_document(d_folder)
