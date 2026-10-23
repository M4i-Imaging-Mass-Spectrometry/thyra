"""``--waters-lock-mass`` is a transport for ``WatersReader(lock_mass=...)`` (D84).

Omitted, nothing is forwarded and the store keeps what the run carries.
Given, it is a reference m/z or ``off``, checked before anything is read.
"""

from __future__ import annotations

import click
import pytest

from thyra.__main__ import _build_reader_options, _parse_lock_mass


class TestParsing:
    @pytest.mark.parametrize(
        "text, value", [("526.5545", 526.5545), ("off", "off"), ("OFF", "off")]
    )
    def test_a_reference_or_off_is_accepted(self, text, value):
        assert _parse_lock_mass(None, None, text) == value

    def test_omitted_stays_none(self):
        assert _parse_lock_mass(None, None, None) is None

    @pytest.mark.parametrize("text", ["auto", "-1", "0", "nan", "inf"])
    def test_anything_else_is_refused(self, text):
        with pytest.raises(click.BadParameter):
            _parse_lock_mass(None, None, text)


class TestBuildReaderOptions:
    def test_omitted_is_not_forwarded(self):
        assert "lock_mass" not in _build_reader_options(True, None)

    @pytest.mark.parametrize("value", [526.5545, "off"])
    def test_a_value_is_forwarded_as_the_reader_keyword(self, value):
        options = _build_reader_options(True, None, waters_lock_mass=value)
        assert options["lock_mass"] == value
