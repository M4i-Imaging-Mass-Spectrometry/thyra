# tests/unit/test_cli_help_grouping.py
"""Tests that ``thyra --help`` really is organised by category.

``GroupedCommand.format_help`` sweeps anything missing from
``GroupedCommand.GROUPS`` into a trailing ``Options:`` section. That makes a
forgotten option invisible as a bug -- the help still renders, the option is
still listed, and only a reader comparing it against docs/cli.md's promise of
"all options organised by category" would notice. Two options sat there for
two releases exactly that way.

These tests make the sweep bucket assert-on-non-empty, so classifying a new
option is part of adding it rather than a thing to remember.
"""

from __future__ import annotations

import click
import pytest
from click.testing import CliRunner

from thyra.__main__ import GroupedCommand, main


@pytest.fixture
def ctx():
    with click.Context(main) as context:
        yield context


def _visible_params(ctx):
    """Options that actually render in --help.

    Hidden options (``--optimize-chunks``) return no help record and are
    excluded, which is correct: they are not shown, so they need no section.
    """
    return [
        p for p in ctx.command.get_params(ctx) if p.get_help_record(ctx) is not None
    ]


def _spellings(param):
    """Every name a param answers to, primary and secondary."""
    return set(getattr(param, "opts", [])) | set(getattr(param, "secondary_opts", []))


def _groups_for(param):
    """The GROUPS sections that claim this param, by any of its spellings."""
    names = _spellings(param)
    return {
        group
        for group, members in GroupedCommand.GROUPS.items()
        if names & set(members)
    }


class TestEveryOptionIsClassified:
    """The mapping between options and sections must be total and unique."""

    def test_every_visible_option_is_in_exactly_one_group(self, ctx):
        misfiled = {
            sorted(_spellings(p))[0]: sorted(_groups_for(p))
            for p in _visible_params(ctx)
            if len(_groups_for(p)) != 1
        }

        assert not misfiled, (
            "every option shown in --help must belong to exactly one "
            "GroupedCommand.GROUPS section. Offenders map to the sections "
            f"that claim them (an empty list means none): {misfiled}"
        )

    def test_groups_lists_no_option_that_does_not_exist(self, ctx):
        real = set()
        for param in _visible_params(ctx):
            real |= _spellings(param)

        phantom = {
            group: [name for name in members if name not in real]
            for group, members in GroupedCommand.GROUPS.items()
            if any(name not in real for name in members)
        }

        assert not phantom, (
            "GroupedCommand.GROUPS names options the command does not have. "
            "A renamed or removed option leaves its old spelling behind here, "
            "where it silently stops matching anything: " + repr(phantom)
        )


class TestRenderedHelp:
    """The end-to-end assertion: what a user actually sees."""

    def test_help_has_no_ungrouped_options_section(self):
        result = CliRunner().invoke(main, ["--help"])

        assert result.exit_code == 0, result.output
        assert "\nOptions:\n" not in result.output, (
            "--help ended with an ungrouped 'Options:' section, so at least "
            "one option is missing from GroupedCommand.GROUPS. Rendered "
            "help:\n" + result.output
        )

    def test_help_renders_every_declared_section(self):
        result = CliRunner().invoke(main, ["--help"])

        assert result.exit_code == 0, result.output
        missing = [
            group
            for group in GroupedCommand.GROUPS
            if f"\n{group}:\n" not in result.output
        ]

        assert not missing, (
            "sections declared in GROUPS did not render, which means every "
            f"option in them has gone: {missing}"
        )

    @pytest.mark.parametrize(
        "option, group",
        [
            ("--resample-gap-tolerance", "Resampling (advanced)"),
            ("--spectrum-type", "imzML-specific"),
            ("--version", "General"),
            ("--help", "General"),
        ],
    )
    def test_previously_ungrouped_options_are_classified(self, option, group):
        """The four that were in the sweep bucket, pinned to their sections."""
        assert option in GroupedCommand.GROUPS[group]
