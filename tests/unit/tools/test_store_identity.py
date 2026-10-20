"""Tests for the store-identity harness, ``tests/tools/store_identity.py``.

The harness is the instrument several design decisions rest on, so it must
never report identity for two stores that differ. Issue #529 found it did,
in each of the ways below; every test builds two small stores that differ
in one of them, or drives ``run`` with a stand-in ``thyra`` package.

These live under ``tests/unit/`` rather than beside the harness because a
test outside ``tests/unit/`` and ``tests/integration/`` belongs to no lane
(``tests/unit/test_lane_markers.py``).
"""

import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import zarr

from thyra.metadata.personal_data import LOCATING_KEYS

REPO_ROOT = Path(__file__).resolve().parents[3]
HARNESS = REPO_ROOT / "tests" / "tools" / "store_identity.py"


def _load_harness():
    spec = importlib.util.spec_from_file_location("store_identity", HARNESS)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


si = _load_harness()

LOCATION_KEYS = ["source_path", *sorted(LOCATING_KEYS)]


def _store(path: Path, **change) -> Path:
    """A small store shaped like a converted one; ``change`` varies one thing.

    Every value a test can vary has a default here, so two calls with the
    same arguments build two stores that must compare equal.
    """
    values = {
        "version": "4.5.1",
        "timestamp": "2026-10-20T12:00:00",
        "offsets": [1, 1, 0],
        "ordered": False,
        "boxes": b"PAR1 pixel boxes PAR1",
        "input": "/data/here/sample.imzML",
        "other_software_version": "2.0",
        "intensities": [1.0, 2.0, 3.0],
        "regions": '{"region": "pixels", "count": 3}',
        "consolidate": True,
    }
    values.update(change)
    root = zarr.open_group(str(path), mode="w")
    root.attrs.update(
        {
            "msi_converter_version": values["version"],
            "conversion_timestamp": values["timestamp"],
            "offsets_px": values["offsets"],
            "coordinate_systems": {
                "global": {
                    "unit": "micrometer",
                    "produced_by": f"thyra/{values['version']}",
                }
            },
        }
    )
    table = root.create_group("tables/t")
    table.attrs["ordered"] = values["ordered"]
    data = table.create_array("X", shape=(3,), dtype="f8")
    data[:] = np.asarray(values["intensities"])
    data.attrs["ordered"] = values["ordered"]

    uns = table.create_group("uns")
    processing = [
        {
            "name": "conversion",
            "software": {"name": "thyra", "version": values["version"]},
        },
        {
            "name": "calibration",
            "software": {"name": "vendor", "version": values["other_software_version"]},
        },
    ]
    for name, text in (
        ("thyra_version", values["version"]),
        ("processing", json.dumps(processing)),
        ("source_path", values["input"]),
        ("regions", values["regions"]),
    ):
        node = uns.create_array(name, shape=(), dtype=str)
        node[...] = text
    vendor = uns.create_group("format_specific")
    for key in sorted(LOCATING_KEYS):
        node = vendor.create_array(key, shape=(), dtype=str)
        node[...] = values["input"]

    boxes = path / "shapes" / "pixels"
    root.create_group("shapes/pixels")  # SpatialData writes the group first
    (boxes / "shapes.parquet").write_bytes(values["boxes"])
    if values["consolidate"]:  # as a converted store is
        zarr.consolidate_metadata(str(path))
    return path


def _differences(tmp_path, mask_keys=None, **change):
    keys = si.DEFAULT_MASK_KEYS if mask_keys is None else mask_keys
    a = si.build_manifest(_store(tmp_path / "a.zarr"), keys)
    b = si.build_manifest(_store(tmp_path / "b.zarr", **change), keys)
    # What `compare` reads: the manifest after a round trip through disk.
    a, b = (json.loads(json.dumps(m, sort_keys=True)) for m in (a, b))
    return si._diff_manifests("demo", a, b)


class TestTwoStoresThatDifferAreNeverIdentical:
    def test_the_same_store_twice_has_no_difference(self, tmp_path):
        """The control: everything below would pass on a harness that
        reports every store as different."""
        assert _differences(tmp_path) == 0

    def test_a_timestamp_alone_is_not_a_difference(self, tmp_path):
        assert _differences(tmp_path, timestamp="2027-01-01T00:00:00") == 0

    def test_the_array_data(self, tmp_path):
        assert _differences(tmp_path, intensities=[1.0, 2.0, 4.0]) == 1

    def test_the_pixel_boxes(self, tmp_path, capsys):
        """F108: ``shapes.parquet`` is not a Zarr array and was never hashed."""
        assert _differences(tmp_path, boxes=b"PAR1 boxes moved 1 um PAR1") == 1
        assert "demo shapes/pixels/shapes.parquet hash" in capsys.readouterr().out

    def test_a_json_type_change_in_attrs(self, tmp_path, capsys):
        """F166: ``1 == 1.0 == True`` in Python, so group attrs holding
        them compared equal; an array's metadata already did not."""
        assert _differences(tmp_path, offsets=[1.0, 1.0, 0.0], ordered=0) == 4
        out = capsys.readouterr().out
        assert "demo <root> attrs" in out
        assert "demo tables/t attrs" in out
        assert "demo tables/t/X metadata" in out
        # The consolidated copy holds the two nested values as well.
        assert "demo <consolidated metadata>" in out

    def test_another_softwares_version(self, tmp_path):
        """Only Thyra's own version is masked in the processing list."""
        assert _differences(tmp_path, other_software_version="2.1") == 1

    def test_a_manifest_from_before_file_hashing_is_not_trusted(self, tmp_path):
        manifest = si.build_manifest(_store(tmp_path / "a.zarr"), si.DEFAULT_MASK_KEYS)
        old = {key: value for key, value in manifest.items() if key != "files"}
        assert si._diff_manifests("demo", old, old) == 1

    def test_key_order_inside_a_stored_json_payload(self, tmp_path):
        """A payload with nothing to mask is hashed as it was written; it
        used to be re-serialised with sorted keys first."""
        assert _differences(tmp_path, regions='{"count": 3, "region": "pixels"}') == 1

    def test_a_store_without_consolidated_metadata(self, tmp_path, capsys):
        assert _differences(tmp_path, consolidate=False) == 1
        assert "demo <consolidated metadata>" in capsys.readouterr().out

    def test_an_attr_written_after_consolidation(self, tmp_path, capsys):
        """The node's own ``zarr.json`` changed and the consolidated copy did
        not. The walk read the copy, so this passed as identical."""
        a = si.build_manifest(_store(tmp_path / "a.zarr"), si.DEFAULT_MASK_KEYS)
        store_b = _store(tmp_path / "b.zarr")
        node = zarr.open_group(str(store_b), mode="r+", use_consolidated=False)
        node["tables/t"].attrs["ordered"] = True
        b = si.build_manifest(store_b, si.DEFAULT_MASK_KEYS)

        assert si._diff_manifests("demo", a, b) == 1
        assert "demo tables/t attrs" in capsys.readouterr().out


class TestWhatIsMasked:
    def test_a_version_only_change_is_no_difference(self, tmp_path):
        """F112: the version sits in root attrs, in ``produced_by``, in the
        ``thyra_version`` arrays and in the processing list's software."""
        assert _differences(tmp_path, version="9.9.9") == 0

    def test_without_the_masks_it_is_one(self, tmp_path):
        assert _differences(tmp_path, mask_keys=set(), version="9.9.9") > 0

    def test_the_location_recipe_masks_every_place_the_input_path_is(self, tmp_path):
        """F111: ``--mask source_path`` alone left five keys holding it."""
        moved = {"input": "/data/there/sample.imzML"}
        only_source_path = si.DEFAULT_MASK_KEYS | {"source_path"}
        recipe = si.DEFAULT_MASK_KEYS | set(LOCATION_KEYS)

        assert _differences(tmp_path / "1", mask_keys=only_source_path, **moved) > 0
        assert _differences(tmp_path / "2", mask_keys=recipe, **moved) == 0

    def test_the_recipe_names_every_locating_key(self):
        """The docstring and the ``--mask`` help name every key that
        ``thyra/metadata/personal_data.py`` keeps as a location."""
        help_text = si._build_run_parser().format_help()
        for key in LOCATION_KEYS:
            assert f"--mask {key}" in (si.__doc__ or "")
            assert key in " ".join(help_text.split())


class TestCompare:
    def _write(self, scratch, label, name, manifest, record=True):
        folder = scratch / label
        folder.mkdir(parents=True, exist_ok=True)
        (folder / f"{name}.hashes.json").write_text(json.dumps(manifest))
        if record:
            si._append_checkpoint(
                folder / "checkpoint.jsonl", {"name": name, "status": "ok"}
            )

    def test_a_label_against_itself_is_refused(self, tmp_path):
        manifest = si.build_manifest(_store(tmp_path / "s.zarr"), si.DEFAULT_MASK_KEYS)
        self._write(tmp_path / "scratch", "main", "demo", manifest)
        assert si.compare_labels("main", "main", tmp_path / "scratch") == 2

    @pytest.mark.parametrize("other", ["main/", "./main"])
    def test_another_spelling_of_one_label_is_refused(self, tmp_path, other):
        manifest = si.build_manifest(_store(tmp_path / "s.zarr"), si.DEFAULT_MASK_KEYS)
        self._write(tmp_path / "scratch", "main", "demo", manifest)
        assert si.compare_labels("main", other, tmp_path / "scratch") == 2

    def test_a_manifest_without_an_ok_record_is_not_compared(self, tmp_path, capsys):
        """A deleted checkpoint left manifests that no record vouched for."""
        manifest = si.build_manifest(_store(tmp_path / "s.zarr"), si.DEFAULT_MASK_KEYS)
        self._write(tmp_path / "scratch", "a", "demo", manifest, record=False)
        self._write(tmp_path / "scratch", "a", "kept", manifest)
        for name in ("demo", "kept"):
            self._write(tmp_path / "scratch", "b", name, manifest)

        assert si.compare_labels("a", "b", tmp_path / "scratch") == 1
        assert "demo <unvouched>" in capsys.readouterr().out

    def test_nothing_to_compare_is_not_identity(self, tmp_path, capsys):
        """F110: no manifest under either label used to exit 0."""
        (tmp_path / "failed_a").mkdir()
        assert si.compare_labels("failed_a", "failed_b", tmp_path) == 2
        out = capsys.readouterr().out
        assert "no manifest" in out and "no such label" in out

    def test_one_mistyped_label_is_a_difference(self, tmp_path):
        manifest = si.build_manifest(_store(tmp_path / "s.zarr"), si.DEFAULT_MASK_KEYS)
        self._write(tmp_path / "scratch", "a", "demo", manifest)
        assert si.compare_labels("a", "aa", tmp_path / "scratch") == 1

    def test_identical_labels_exit_0(self, tmp_path):
        manifest = si.build_manifest(_store(tmp_path / "s.zarr"), si.DEFAULT_MASK_KEYS)
        for label in ("a", "b"):
            self._write(tmp_path / "scratch", label, "demo", manifest)
        assert si.compare_labels("a", "b", tmp_path / "scratch") == 0


FAKE_THYRA = """\
import json
import sys
from pathlib import Path

source, out, options = sys.argv[1], sys.argv[2], sys.argv[3:]
stem = Path(source).stem
if stem == "asks":
    if not sys.stdin.readline():
        print("t - thyra.cli - ERROR - no answer on stdin")
        sys.exit(1)
if stem == "refused":
    print("t - thyra.convert - ERROR - refused on purpose")
    print("t - thyra.cli - ERROR - Conversion failed.")
    print("Traceback (most recent call last)", file=sys.stderr)
    sys.exit(1)
# A one-group Zarr v3 store, written by hand: importing zarr in every
# child would cost most of each test's time.
Path(out).mkdir(parents=True)
group = {"zarr_format": 3, "node_type": "group", "attributes": {"options": options}}
(Path(out) / "zarr.json").write_text(json.dumps(group))
"""


@pytest.fixture
def fake_source(tmp_path, monkeypatch):
    """A checkout whose ``thyra`` writes a store recording its options.

    ``run`` puts the source first on ``PYTHONPATH`` and runs in it, so
    ``python -m thyra`` finds this package, not the real one. The
    datasets live in a root on ``THYRA_DATASET_ROOTS``, and the home
    folder is empty, so no machine registry is read.
    """
    source = tmp_path / "checkout"
    (source / "thyra").mkdir(parents=True)
    (source / "thyra" / "__init__.py").write_text("")
    (source / "thyra" / "__main__.py").write_text(FAKE_THYRA)
    data = tmp_path / "data"
    data.mkdir()
    for stem in ("good", "other", "refused", "asks"):
        (data / f"{stem}.imzML").write_text("")
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("THYRA_DATASET_ROOTS", str(data))
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    return source


def _run(
    source, scratch, datasets=("good",), options=(), masks=(), force=False, label="x"
):
    return si.run_datasets(
        label=label,
        datasets=list(datasets),
        source=source,
        interpreter=sys.executable,
        scratch=scratch,
        thyra_options=list(options),
        mask_keys=si.DEFAULT_MASK_KEYS | set(masks),
        force=force,
    )


class TestCompareAfterRun:
    def test_a_label_partly_redone_under_new_code_is_refused(
        self, fake_source, tmp_path, capsys
    ):
        """Redoing one name of a label left the others' stores from the old
        code, and ``compare`` let them answer for the new code."""
        scratch = tmp_path / "scratch"
        _run(fake_source, scratch, datasets=["good", "other"], label="br")
        main = fake_source / "thyra" / "__main__.py"
        main.write_text(main.read_text() + "# an edit\n")
        _run(fake_source, scratch, datasets=["good"], label="br")
        _run(fake_source, scratch, datasets=["good", "other"], label="main")
        capsys.readouterr()

        assert si.compare_labels("main", "br", scratch) == 1
        assert "different settings" in capsys.readouterr().out

    def test_a_name_that_failed_under_both_labels_is_reported(
        self, fake_source, tmp_path, capsys
    ):
        scratch = tmp_path / "scratch"
        for label in ("a", "b"):
            _run(fake_source, scratch, datasets=["good", "refused"], label=label)
        capsys.readouterr()

        assert si.compare_labels("a", "b", scratch) == 1
        assert "refused <failed>" in capsys.readouterr().out

    def test_a_name_skipped_under_both_labels_is_reported(
        self, fake_source, tmp_path, capsys
    ):
        scratch = tmp_path / "scratch"
        for label in ("a", "b"):
            _run(fake_source, scratch, datasets=["good", "nosuch"], label=label)
        capsys.readouterr()

        assert si.compare_labels("a", "b", scratch) == 1
        assert "nosuch <skipped>" in capsys.readouterr().out

    def test_two_runs_of_the_same_code_compare_equal(self, fake_source, tmp_path):
        """The control for the three above."""
        scratch = tmp_path / "scratch"
        for label in ("a", "b"):
            _run(fake_source, scratch, datasets=["good", "other"], label=label)
        assert si.compare_labels("a", "b", scratch) == 0


class TestRun:
    def test_a_rerun_under_other_options_converts_again(
        self, fake_source, tmp_path, capsys
    ):
        """F109: the cached store from the first options answered for both."""
        scratch = tmp_path / "scratch"
        _run(fake_source, scratch)
        _run(fake_source, scratch)
        assert "[skip-cached] good" in capsys.readouterr().out

        _run(fake_source, scratch, options=["--resample-bins", "50"])
        assert "[rerun] good" in capsys.readouterr().out
        store = zarr.open_group(str(scratch / "x" / "good.zarr"), mode="r")
        assert store.attrs["options"] == ["--resample-bins", "50"]

    def test_an_uncommitted_code_edit_converts_again(
        self, fake_source, tmp_path, capsys
    ):
        """The commit alone keyed the cache, but the child runs the files on
        disk: a second edit at one commit reused the first edit's store."""
        scratch = tmp_path / "scratch"
        _run(fake_source, scratch)
        main = fake_source / "thyra" / "__main__.py"
        main.write_text(main.read_text() + "# an edit\n")

        _run(fake_source, scratch)
        assert "[rerun] good: label 'x' was made with other code" in (
            capsys.readouterr().out
        )

    def test_a_rewritten_input_converts_again(self, fake_source, tmp_path, capsys):
        """The path alone keyed the cache: a dataset re-exported to the same
        place reused the old store. The .ibd beside an imzML counts too."""
        scratch = tmp_path / "scratch"
        _run(fake_source, scratch)
        (tmp_path / "data" / "good.ibd").write_bytes(b"new binary data")

        _run(fake_source, scratch)
        assert "[rerun] good: label 'x' was made with other input_state" in (
            capsys.readouterr().out
        )

    def test_a_rewritten_file_named_in_an_option_converts_again(
        self, fake_source, tmp_path, capsys
    ):
        scratch = tmp_path / "scratch"
        named = tmp_path / "optical.tif"
        named.write_bytes(b"one")
        _run(fake_source, scratch, options=["--optical", str(named)])
        named.write_bytes(b"other image")

        _run(fake_source, scratch, options=["--optical", str(named)])
        assert "[rerun] good" in capsys.readouterr().out

    def test_another_sdk_path_converts_again(
        self, fake_source, tmp_path, capsys, monkeypatch
    ):
        scratch = tmp_path / "scratch"
        _run(fake_source, scratch)
        monkeypatch.setenv("BRUKER_SDK_PATH", str(tmp_path / "another_sdk"))

        _run(fake_source, scratch)
        assert "[rerun] good: label 'x' was made with other environment" in (
            capsys.readouterr().out
        )

    def test_a_source_without_thyra_is_refused(self, fake_source, tmp_path):
        """The child imported the installed thyra, another checkout's."""
        with pytest.raises(si.SourceNotUsed, match="holds no thyra package"):
            _run(fake_source / "thyra", tmp_path / "scratch")
        with pytest.raises(si.SourceNotUsed, match="holds no thyra package"):
            _run(fake_source.parent, tmp_path / "scratch")

    def test_a_source_whose_thyra_is_not_the_one_imported_is_refused(
        self, fake_source, monkeypatch
    ):
        elsewhere = subprocess.CompletedProcess(
            [], 0, stdout="/another/checkout/thyra/__init__.py\n", stderr=""
        )
        monkeypatch.setattr(si.subprocess, "run", lambda cmd, **kwargs: elsewhere)
        with pytest.raises(si.SourceNotUsed, match="not from"):
            si._check_source(fake_source, sys.executable)

    def test_an_environment_that_cannot_be_listed_matches_nothing(self, tmp_path):
        """``None == None`` let two unknown environments share a store."""
        missing = str(tmp_path / "no-such-python")
        assert si._environment_hash(missing) != si._environment_hash(missing)

    def test_the_fingerprint_covers_what_a_conversion_reads(self, tmp_path):
        """Files beside the input (a ``.mis`` next to a ``.d``) and files
        named as ``--flag=value``; not the log file the child writes."""
        run = tmp_path / "slide" / "run.d"
        run.mkdir(parents=True)
        (run / "analysis.tsf").write_bytes(b"frames")
        mis = tmp_path / "slide" / "run.mis"
        mis.write_text("<teaching points/>")
        plan = tmp_path / "plan.json"
        plan.write_text("{}")
        log = tmp_path / "run.log"
        log.write_text("")
        options = [f"--mass-axis-plan={plan}", "--log-file", str(log)]

        def fingerprint():
            return si._input_fingerprint(run, options, tmp_path)

        before = fingerprint()
        log.write_text("the child wrote its log")
        assert fingerprint() == before
        mis.write_text("<teaching points moved/>")
        after_mis = fingerprint()
        assert after_mis != before
        plan.write_text('{"bins": 50}')
        assert fingerprint() != after_mis

    def test_code_changed_during_the_run_fails_the_name(
        self, fake_source, tmp_path, monkeypatch
    ):
        hashes = iter(["before", "after", "after"])
        monkeypatch.setattr(si, "_code_hash", lambda source: next(hashes))

        results = _run(fake_source, tmp_path / "scratch")
        assert results[0]["status"] == "failed"
        assert "changed during the run" in results[0]["error_lines"][0]
        assert not (tmp_path / "scratch" / "x" / "good.hashes.json").exists()

    def test_a_rerun_under_other_masks_converts_again(
        self, fake_source, tmp_path, capsys
    ):
        scratch = tmp_path / "scratch"
        _run(fake_source, scratch)
        _run(fake_source, scratch, masks=["source_path"])
        assert "[rerun] good" in capsys.readouterr().out

    def test_a_record_from_before_the_settings_were_kept_is_redone(
        self, fake_source, tmp_path, capsys
    ):
        scratch = tmp_path / "scratch"
        _run(fake_source, scratch)
        checkpoint = scratch / "x" / "checkpoint.jsonl"
        old = {"name": "good", "options": [], "status": "ok"}
        checkpoint.write_text(json.dumps(old) + "\n")

        _run(fake_source, scratch)
        assert "[rerun] good" in capsys.readouterr().out

    def test_a_skipped_name_drops_its_manifest_from_an_earlier_run(
        self, fake_source, tmp_path, monkeypatch
    ):
        scratch = tmp_path / "scratch"
        _run(fake_source, scratch)
        assert (scratch / "x" / "good.hashes.json").exists()

        monkeypatch.setenv("THYRA_DATASET_ROOTS", str(tmp_path / "offline"))
        (fake_source.parent / "data" / "good.imzML").unlink()
        _run(fake_source, scratch)  # same settings: still not this run's store
        assert not (scratch / "x" / "good.hashes.json").exists()

    def test_a_failed_record_keeps_the_error_lines_from_stdout(
        self, fake_source, tmp_path, capsys
    ):
        """F113: Thyra logs to stdout, and the record kept only stderr."""
        results = _run(fake_source, tmp_path / "scratch", datasets=["refused"])

        record = results[0]
        assert record["status"] == "failed"
        assert record["error_lines"][0].endswith("refused on purpose")
        assert "refused on purpose" in capsys.readouterr().out

    def test_the_child_cannot_wait_for_an_answer(self, fake_source, tmp_path):
        """F49: with stdin open, a thyra that asks a question waits forever."""
        results = _run(fake_source, tmp_path / "scratch", datasets=["asks"])
        assert results[0]["status"] == "failed"
        assert results[0]["error_lines"][0].endswith("no answer on stdin")

    def test_the_child_runs_with_stdin_closed(self, fake_source, tmp_path, monkeypatch):
        seen = {}

        def fake_run(cmd, **kwargs):
            seen.update(kwargs)
            return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="")

        monkeypatch.setattr(si.subprocess, "run", fake_run)
        monkeypatch.setattr(si, "_check_source", lambda source, interpreter: None)
        _run(fake_source, tmp_path / "scratch")
        assert seen["stdin"] is subprocess.DEVNULL


class TestDatasetResolution:
    def _registry(self, tmp_path, monkeypatch, folder):
        home = tmp_path / "home"
        home.mkdir()
        (home / ".msi-datasets.toml").write_text(
            f"[local]\nslide = '{folder.as_posix()}'\n"
        )
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setenv("USERPROFILE", str(home))

    def test_a_registered_folder_of_several_runs_is_skipped(
        self, tmp_path, monkeypatch
    ):
        """F49: passed as it was, thyra asked on stdin which run to take."""
        folder = tmp_path / "slide"
        (folder / "a.d").mkdir(parents=True)
        (folder / "b.d").mkdir()
        self._registry(tmp_path, monkeypatch, folder)
        with pytest.raises(si.AmbiguousDatasetName, match="2 .d runs"):
            si.resolve_dataset("slide")

    def test_runs_are_counted_as_thyra_counts_them(self, tmp_path, monkeypatch):
        """Folders only, the suffix in any letter case."""
        folder = tmp_path / "slide"
        (folder / "a.d").mkdir(parents=True)
        (folder / "B.D").mkdir()
        (folder / "notes.d").write_text("a file, not a run")
        self._registry(tmp_path, monkeypatch, folder)
        with pytest.raises(si.AmbiguousDatasetName, match="2 .d runs"):
            si.resolve_dataset("slide")

    def test_a_registered_folder_of_one_run_resolves_to_it(self, tmp_path, monkeypatch):
        folder = tmp_path / "slide"
        (folder / "a.d").mkdir(parents=True)
        self._registry(tmp_path, monkeypatch, folder)
        assert si.resolve_dataset("slide") == folder / "a.d"

    def test_a_registered_folder_with_no_run_is_passed_as_it_is(
        self, tmp_path, monkeypatch
    ):
        """A Waters .raw or a rapifleX folder holds no .d and is the input."""
        folder = tmp_path / "sample.raw"
        folder.mkdir()
        self._registry(tmp_path, monkeypatch, folder)
        assert si.resolve_dataset("slide") == folder

    def test_every_default_dataset_resolves_to_a_fixture(self, tmp_path, monkeypatch):
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.setenv("USERPROFILE", str(tmp_path))
        monkeypatch.delenv("THYRA_DATASET_ROOTS", raising=False)
        fixtures = REPO_ROOT / "tests" / "data" / "fixtures"
        for name in si.DEFAULT_DATASETS:
            resolved = si.resolve_dataset(name)
            assert resolved is not None and resolved.parent == fixtures, name

    def test_every_fixture_is_in_the_default_set_or_named_as_left_out(self):
        """F114: two fixtures that convert were missing from the set."""
        source = HARNESS.read_text(encoding="utf-8")
        comment = source.split("DEFAULT_DATASETS = (")[0].rsplit("\n\n", 1)[-1]
        named = set(re.findall(r"[A-Za-z0-9_]+", comment))
        fixtures = REPO_ROOT / "tests" / "data" / "fixtures"
        for path in [*fixtures.glob("*.imzML"), *fixtures.glob("*.d")]:
            assert path.stem in si.DEFAULT_DATASETS or path.stem in named, path.name
