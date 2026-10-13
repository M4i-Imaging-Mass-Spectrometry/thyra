"""``plan_mass_axis`` reads sources' headers and plans from them (issue #449, D43)."""

from __future__ import annotations

from pathlib import Path

import pytest

import thyra
from thyra.errors import ConversionRefused
from thyra.mass_axis_plan import _source_names, plan_mass_axis
from thyra.resampling.plan import MassAxisPlan

_FIXTURES = Path(__file__).resolve().parents[1] / "data" / "fixtures"
_IMZML = _FIXTURES / "mobility_continuous.imzML"


def test_the_names_are_public_and_lazy():
    assert thyra.plan_mass_axis is plan_mass_axis
    assert thyra.MassAxisPlan is MassAxisPlan
    assert {"plan_mass_axis", "MassAxisPlan"} <= set(thyra.__all__)


def test_a_source_that_does_not_exist_is_refused_by_name(tmp_path):
    with pytest.raises(ConversionRefused, match="missing.imzML"):
        plan_mass_axis([_IMZML, tmp_path / "missing.imzML"])


def test_a_source_that_cannot_be_read_is_refused_by_name(tmp_path):
    broken = tmp_path / "broken.imzML"
    broken.write_text("not xml", encoding="utf-8")

    with pytest.raises(ConversionRefused, match="broken.imzML"):
        plan_mass_axis([broken])


def test_a_file_that_states_no_range_needs_both_bounds():
    """This fixture records no observed m/z, and a plan will not decode spectra."""
    with pytest.raises(ConversionRefused, match="mobility_continuous.imzML"):
        plan_mass_axis(_IMZML)

    plan = plan_mass_axis(_IMZML, min_mz=100.0, max_mz=900.0)

    assert plan.sources == (_IMZML.name,)
    assert plan.n_bins > 0 and len(plan.sha256) == 64


def test_names_carry_the_folder_only_where_two_would_collide(tmp_path):
    paths = [tmp_path / "x" / "a.imzML", tmp_path / "y" / "a.imzML", tmp_path / "b"]

    assert _source_names(paths) == ["x/a.imzML", "y/a.imzML", "b"]


def test_duplicates_apart_in_the_list_are_numbered_in_turn(tmp_path):
    a, b = tmp_path / "x" / "a.imzML", tmp_path / "y" / "b.imzML"

    assert _source_names([a, b, a]) == ["x/a.imzML (1)", "b.imzML", "x/a.imzML (2)"]


@pytest.mark.parametrize("copies", [2, 3])
def test_the_same_path_given_more_than_once_gets_a_number_each(tmp_path, copies):
    paths = [tmp_path / "x" / "a.imzML"] * copies

    names = _source_names(paths)

    assert len(set(names)) == copies
    assert names == [f"x/a.imzML ({i + 1})" for i in range(copies)]
