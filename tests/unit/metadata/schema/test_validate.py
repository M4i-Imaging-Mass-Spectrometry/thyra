"""Document validation: version rules, structure, ontology checks."""

from thyra.metadata.schema import (
    MSI_METADATA_SCHEMA_VERSION,
    build_msi_metadata,
    validate_document,
)


def _valid_doc() -> dict:
    return build_msi_metadata(None, pixel_size_um=(20.0, 20.0)).to_uns_dict()


def _errors(issues):
    return [i for i in issues if i.severity == "error"]


def _warnings(issues):
    return [i for i in issues if i.severity == "warning"]


class TestSchemaVersion:
    def test_valid_document_has_no_issues(self):
        meta, issues = validate_document(_valid_doc())
        assert meta is not None
        assert issues == []

    def test_missing_schema_version_is_an_error(self):
        doc = _valid_doc()
        del doc["schema_version"]
        meta, issues = validate_document(doc)
        assert meta is None
        assert any(i.location == "schema_version" for i in _errors(issues))

    def test_different_major_version_is_an_error(self):
        doc = _valid_doc()
        doc["schema_version"] = "99.0.0"
        meta, issues = validate_document(doc)
        assert meta is None
        assert _errors(issues)

    def test_newer_minor_version_warns_but_validates(self):
        doc = _valid_doc()
        major, minor, patch = (int(p) for p in MSI_METADATA_SCHEMA_VERSION.split("."))
        doc["schema_version"] = f"{major}.{minor + 1}.{patch}"
        meta, issues = validate_document(doc)
        assert meta is not None
        assert _warnings(issues)
        assert not _errors(issues)

    def test_garbage_version_is_an_error(self):
        doc = _valid_doc()
        doc["schema_version"] = "not-a-version"
        meta, issues = validate_document(doc)
        assert meta is None


class TestStructure:
    def test_non_mapping_document_is_an_error(self):
        meta, issues = validate_document(["not", "a", "mapping"])
        assert meta is None
        assert _errors(issues)

    def test_structural_errors_carry_their_location(self):
        doc = _valid_doc()
        doc["ms_analysis"]["pixel_size_um"]["x"] = -1.0
        meta, issues = validate_document(doc)
        assert meta is None
        assert any(i.location == "ms_analysis.pixel_size_um.x" for i in _errors(issues))


class TestOntologyChecks:
    def test_unknown_ms_accession_is_an_error(self):
        doc = _valid_doc()
        doc["ms_analysis"]["analyzer_term"] = {
            "accession": "MS:9999999",
            "name": "made up",
        }
        _, issues = validate_document(doc)
        assert any(i.location == "ms_analysis.analyzer_term" for i in _errors(issues))

    def test_mismatched_term_name_is_a_warning(self):
        doc = _valid_doc()
        doc["ms_analysis"]["analyzer_term"] = {
            "accession": "MS:1000084",
            "name": "definitely not time-of-flight",
        }
        meta, issues = validate_document(doc)
        assert meta is not None
        assert any(i.location == "ms_analysis.analyzer_term" for i in _warnings(issues))
        assert not _errors(issues)

    def test_wrong_prefix_for_designated_ontology_is_a_warning(self):
        doc = _valid_doc()
        doc["sample"] = {
            "organism": "Mus musculus",
            "organism_term": {"accession": "CHEBI:12345", "name": "mouse"},
        }
        meta, issues = validate_document(doc)
        assert meta is not None
        assert any(i.location == "sample.organism_term" for i in _warnings(issues))

    def test_foreign_prefixes_are_not_checked_against_local_tables(self):
        # NCBITaxon is not shipped locally; a plausible term must pass
        # without an existence check.
        doc = _valid_doc()
        doc["sample"] = {
            "organism": "Mus musculus",
            "organism_term": {
                "accession": "NCBITaxon:10090",
                "name": "Mus musculus",
            },
        }
        meta, issues = validate_document(doc)
        assert meta is not None
        assert issues == []
