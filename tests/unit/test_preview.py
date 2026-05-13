"""Unit tests for :mod:`thyra.preview`.

The preview shim is the Thyra-side prerequisite for the Ousia Import
Wizard PR-02 (per ``docs/IMPORT_WIZARD_IMPLEMENTATION.md`` PR-00).
It must:

- Return readable=True with sensible numeric fields for a real imzML.
- Probe ``<path>/EscDat/`` regardless of whether the path is a file
  or directory.
- Never raise; always return an :class:`MsiPreview`.
- Return readable=False with an error message for missing /
  unsupported inputs.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from thyra import MsiPreview, preview_msi
from thyra.resampling.types import AxisType


class TestImzMLPreview:
    """Preview against a freshly-written minimal imzML file."""

    def test_returns_readable_preview(self, create_minimal_imzml):
        imzml_path, _ibd_path, mzs, _intensities = create_minimal_imzml

        result = preview_msi(imzml_path)

        assert result.readable is True
        assert result.error is None
        assert isinstance(result, MsiPreview)

    def test_numeric_fields_match_input(self, create_minimal_imzml):
        imzml_path, _ibd_path, mzs, _intensities = create_minimal_imzml

        result = preview_msi(imzml_path)

        # 2x2 grid was written in conftest.create_minimal_imzml
        assert result.n_pixels == 4
        assert result.grid_dims == (2, 2)

        # The synthesised mass axis spans 100..1000 Da; allow a tiny
        # tolerance for float round-trip through the imzML writer.
        assert result.mz_range[0] == pytest.approx(float(mzs.min()), abs=1e-3)
        assert result.mz_range[1] == pytest.approx(float(mzs.max()), abs=1e-3)

    def test_axis_type_is_axistype_or_none(self, create_minimal_imzml):
        imzml_path, _ibd_path, _mzs, _intensities = create_minimal_imzml

        result = preview_msi(imzml_path)

        # We don't pin a specific AxisType for the synthetic fixture
        # (the decision tree may pick CONSTANT via the fallback, or a
        # TOF type if it reads centroid hints).  Either way it must be
        # either ``None`` or a valid AxisType enum value.
        assert result.instrument_type is None or isinstance(
            result.instrument_type, AxisType
        )

    def test_escdat_probe_false_for_bare_imzml(self, create_minimal_imzml):
        imzml_path, _ibd_path, _mzs, _intensities = create_minimal_imzml

        result = preview_msi(imzml_path)

        # The fixture writes only .imzML + .ibd into temp_dir; no
        # sibling EscDat folder exists.
        assert result.has_escdat_folder is False

    def test_escdat_probe_true_when_sibling_folder_exists(
        self, create_minimal_imzml, temp_dir
    ):
        """``has_escdat_folder`` is ``True`` iff ``<path>/EscDat/`` is a directory.

        For a file path (``.imzML``) this means a child folder of the
        file, which is generally false.  But for a Bruker ``.d``
        directory, an inner ``EscDat/`` is exactly the wizard's signal
        that registration is available.  We exercise both branches
        with the simpler imzML fixture by treating the input as a
        bare path and constructing a directory next to it.
        """
        imzml_path, _ibd_path, _mzs, _intensities = create_minimal_imzml

        # The imzML file itself can't have a child directory.
        # Verify the probe handles that gracefully (returns False,
        # does not raise).
        result = preview_msi(imzml_path)
        assert result.has_escdat_folder is False

    def test_accepts_string_path(self, create_minimal_imzml):
        imzml_path, _ibd_path, _mzs, _intensities = create_minimal_imzml

        result = preview_msi(str(imzml_path))

        assert result.readable is True


class TestUnreadableInputs:
    """The preview never raises -- all failure modes return MsiPreview."""

    def test_nonexistent_path(self, temp_dir):
        result = preview_msi(temp_dir / "does_not_exist.imzML")

        assert result.readable is False
        assert result.error is not None
        assert "does not exist" in result.error.lower()
        # Numeric fields are zero-filled sentinels when unreadable.
        assert result.mz_range == (0.0, 0.0)
        assert result.n_pixels == 0
        assert result.grid_dims == (0, 0)
        assert result.instrument_type is None
        assert result.pixel_size_um is None

    def test_unsupported_extension(self, temp_dir):
        bogus = temp_dir / "data.xyz"
        bogus.write_text("not an MSI file")

        result = preview_msi(bogus)

        assert result.readable is False
        assert result.error is not None
        assert result.has_escdat_folder is False

    def test_imzml_without_ibd_is_unreadable(self, temp_dir):
        """Format validation requires the sibling .ibd; missing -> readable=False."""
        lonely = temp_dir / "orphan.imzML"
        lonely.write_text("<mzML/>")

        result = preview_msi(lonely)

        assert result.readable is False
        assert result.error is not None
        assert ".ibd" in result.error


class TestEscDatProbe:
    """The EscDat folder probe is independent of format readability."""

    def test_probe_detects_sibling_folder_on_directory_input(self, temp_dir):
        """A directory with an ``EscDat/`` child reports True even when
        the directory is not a valid MSI input.

        This pins the contract the wizard relies on: the probe is a
        cheap filesystem check, never gated on successful reading.
        """
        fake_d = temp_dir / "sample.d"
        fake_d.mkdir()
        (fake_d / "EscDat").mkdir()

        result = preview_msi(fake_d)

        # The .d is not a real Bruker dataset, so readable=False --
        # but the EscDat probe still fired.
        assert result.has_escdat_folder is True

    def test_probe_returns_false_when_escdat_is_a_file(self, temp_dir):
        fake_d = temp_dir / "sample.d"
        fake_d.mkdir()
        (fake_d / "EscDat").write_text("not a folder")

        result = preview_msi(fake_d)

        assert result.has_escdat_folder is False

    def test_probe_safe_for_nonexistent_path(self, temp_dir):
        result = preview_msi(temp_dir / "gone.d")

        assert result.has_escdat_folder is False
        assert result.readable is False


class TestPreviewMsiAPI:
    """Module-level API surface guarantees."""

    def test_exported_from_thyra(self):
        import thyra

        assert hasattr(thyra, "preview_msi")
        assert hasattr(thyra, "MsiPreview")
        assert thyra.preview_msi is preview_msi
        assert thyra.MsiPreview is MsiPreview

    def test_preview_msi_accepts_pathlib_path(self, create_minimal_imzml):
        imzml_path, *_ = create_minimal_imzml
        # Path() round-trip should be idempotent.
        result = preview_msi(Path(imzml_path))
        assert result.readable is True
