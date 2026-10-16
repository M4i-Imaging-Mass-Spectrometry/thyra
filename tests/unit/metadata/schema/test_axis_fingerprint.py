"""The mass axis fingerprint: the hash rule, the model and the store check."""

import hashlib
import json

import numpy as np
import pytest
import zarr
from pydantic import ValidationError

from thyra.metadata.schema import (
    MassAxisFingerprint,
    Provenance,
    check_store_var_conventions,
    mass_axis_sha256,
)

AXIS = np.array([100.0, 100.5, 101.25, 250.125, 999.999], dtype=np.float64)


def _reference(values) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


class TestHashRule:
    def test_equals_hashlib_over_little_endian_float64_bytes(self):
        assert mass_axis_sha256(AXIS) == _reference(AXIS)
        assert len(mass_axis_sha256(AXIS)) == 64

    def test_float32_input_hashes_as_its_float64_values(self):
        narrow = AXIS.astype(np.float32)
        assert mass_axis_sha256(narrow) == _reference(narrow.astype(np.float64))

    def test_big_endian_input_hashes_like_little_endian(self):
        assert mass_axis_sha256(AXIS.astype(">f8")) == mass_axis_sha256(AXIS)

    def test_non_contiguous_view_hashes_like_its_copy(self):
        wide = np.stack([AXIS, AXIS * 2], axis=1)
        view = wide[:, 0]
        assert not view.flags["C_CONTIGUOUS"]
        assert mass_axis_sha256(view) == mass_axis_sha256(view.copy())

    def test_one_ulp_changes_the_hash(self):
        nudged = AXIS.copy()
        nudged[2] = np.nextafter(nudged[2], np.inf)
        assert mass_axis_sha256(nudged) != mass_axis_sha256(AXIS)


class TestModel:
    def test_accepts_a_real_digest(self):
        fingerprint = MassAxisFingerprint(sha256=mass_axis_sha256(AXIS), n_bins=5)
        assert fingerprint.n_bins == 5

    @pytest.mark.parametrize("digest", ["A" * 64, "a" * 63, "g" * 64])
    def test_rejects_a_malformed_digest(self, digest):
        with pytest.raises(ValidationError):
            MassAxisFingerprint(sha256=digest, n_bins=5)

    def test_rejects_zero_bins(self):
        with pytest.raises(ValidationError):
            MassAxisFingerprint(sha256="a" * 64, n_bins=0)

    def test_provenance_without_it_serialises_without_the_key(self):
        dumped = json.loads(
            Provenance(thyra_version="1").model_dump_json(exclude_none=True)
        )
        assert "mass_axis" not in dumped

    def test_provenance_with_it_keeps_it(self):
        provenance = Provenance(
            thyra_version="1",
            mass_axis=MassAxisFingerprint(sha256="b" * 64, n_bins=3),
        )
        dumped = provenance.model_dump(exclude_none=True)
        assert dumped["mass_axis"] == {"sha256": "b" * 64, "n_bins": 3}


def _store_with_block(tmp_path, mz, recorded):
    """A minimal store whose table carries a block recording ``recorded``."""
    import anndata as ad

    root = zarr.open_group(str(tmp_path / "store.zarr"), mode="a")
    table = root.create_group("tables").create_group("t")
    table.create_group("var").create_array("mz", data=np.asarray(mz))
    uns = table.create_group("uns")
    block = {"provenance": {"thyra_version": "1", "mass_axis": recorded}}
    ad.io.write_elem(uns, "msi_metadata", block)
    return tmp_path / "store.zarr"


class TestStoreCheck:
    def test_matching_fingerprint_passes(self, tmp_path):
        recorded = {"sha256": mass_axis_sha256(AXIS), "n_bins": AXIS.size}
        store = _store_with_block(tmp_path, AXIS, recorded)
        assert check_store_var_conventions(store) == {"t": []}

    def test_tampered_hash_is_an_error_at_the_field(self, tmp_path):
        recorded = {"sha256": "0" * 64, "n_bins": AXIS.size}
        store = _store_with_block(tmp_path, AXIS, recorded)
        issues = check_store_var_conventions(store)["t"]
        assert [(i.severity, i.location) for i in issues] == [
            ("error", "provenance.mass_axis")
        ]
        assert "does not match" in issues[0].message

    def test_wrong_bin_count_is_an_error(self, tmp_path):
        recorded = {"sha256": mass_axis_sha256(AXIS), "n_bins": AXIS.size + 1}
        store = _store_with_block(tmp_path, AXIS, recorded)
        issues = check_store_var_conventions(store)["t"]
        assert [i.location for i in issues] == ["provenance.mass_axis"]

    def test_a_field_that_is_not_a_mapping_is_reported_not_raised(self, tmp_path):
        store = _store_with_block(tmp_path, AXIS, "abc")
        issues = check_store_var_conventions(store)["t"]
        assert [(i.severity, i.location) for i in issues] == [
            ("error", "provenance.mass_axis")
        ]
        assert "not a mapping" in issues[0].message

    @pytest.mark.parametrize("n_bins", ["many", [1, 2], 3.5, True, None])
    def test_a_bin_count_that_is_not_an_integer_is_reported_not_raised(
        self, tmp_path, n_bins
    ):
        """Issue #494 (D1): int() raised on text and lists, losing the verdict."""
        recorded = {"sha256": mass_axis_sha256(AXIS)}
        if n_bins is not None:
            recorded["n_bins"] = n_bins
        store = _store_with_block(tmp_path, AXIS, recorded)
        issues = check_store_var_conventions(store)["t"]
        assert [(i.severity, i.location) for i in issues] == [
            ("error", "provenance.mass_axis")
        ]
        assert "no integer n_bins" in issues[0].message

    @pytest.mark.parametrize("kind", [np.int64, float])
    def test_a_whole_bin_count_of_another_type_is_accepted(self, tmp_path, kind):
        recorded = {"sha256": mass_axis_sha256(AXIS), "n_bins": kind(AXIS.size)}
        store = _store_with_block(tmp_path, AXIS, recorded)
        assert check_store_var_conventions(store) == {"t": []}
