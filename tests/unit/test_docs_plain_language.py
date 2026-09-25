# tests/unit/test_docs_plain_language.py
"""The beginner pages stay short and plain.

The documentation has two layers (see "Writing the documentation" in
docs/contributing.md): beginner pages under Home, "Get started" and
"Guides", written for readers who have never used a terminal, and a
technical reference with every detail. In September 2026 the pages grew from
about 24,000 to 71,000 words, mostly because each change appended its own
history, reasons and measurements to the page a newcomer reads first.

These tests hold the beginner pages to the rules that restructure set:

- no sentence runs past ``MAX_SENTENCE_WORDS``, and a page averages at most
  ``MAX_AVERAGE_WORDS`` a sentence, collapsed "Advanced:" boxes included;
- the main text of each page, outside those boxes, stays within its word
  budget, so a page that grows does so by a deliberate edit to
  ``WORD_BUDGETS`` rather than by accretion;
- no version history ("since v3.24", "used to") and no issue numbers: history
  belongs in the changelog, reasons in docs/design-decisions.md.

The pages are read from the nav in mkdocs.yml, so a page added to one of those
sections is checked without touching this file -- and fails until it is given
a budget.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, List, Tuple

import pytest
import yaml

_REPO_ROOT = Path(__file__).resolve().parents[2]
_DOCS = _REPO_ROOT / "docs"
_MKDOCS = _REPO_ROOT / "mkdocs.yml"

#: Nav sections, besides the home page, that hold beginner pages.
BEGINNER_SECTIONS = ("Get started", "Guides")

#: The contributing guide's rule: "split any over 25".
MAX_SENTENCE_WORDS = 25
MAX_AVERAGE_WORDS = 20

#: Words of main text (outside "Advanced:" boxes, code and tables) a page may
#: hold. Raising one is allowed; it should be a decision, taken in review.
#: Set on 2026-09-25 to each page's length then plus about 10%, rounded up to
#: the next 50.
WORD_BUDGETS: Dict[str, int] = {
    "index.md": 350,
    "install.md": 400,
    "getting-started.md": 500,
    "look-at-the-result.md": 350,
    "tutorial.md": 1150,
    "troubleshooting.md": 600,
    "glossary.md": 700,
    "which-files.md": 700,
    "settings.md": 550,
    "describe-your-data.md": 550,
}

_FENCE = re.compile(r"^\s*(```|~~~)")
_ADMONITION = re.compile(r"^(\s*)(?:!!!|\?\?\?\+?)\s+([\w-]+)")
_STRUCTURE = re.compile(
    r"^\s*(?:#{1,6}\s|\||===\s|!\[|</?(?:p|div|img|br)\b|---\s*$|\*\[)"
)
_ITEM = re.compile(r"^(\s*)(?:[-*+]\s+|\d+\.\s+|:\s+)")
_INLINE_CODE = re.compile(r"`[^`]*`")
_IMAGE = re.compile(r"!\[[^\]]*\]\([^)]*\)")
_LINK = re.compile(r"\[([^\]]*)\]\([^)]*\)")
_TAG = re.compile(r"<[^>]+>")
_EMPHASIS = re.compile(r"(\*\*|__|\*)")
_SENTENCE_END = re.compile(r"(?<=[.!?;:])\s+")
_WORD = re.compile(r"\w")

_HISTORY = (
    re.compile(r"\b(?:since|until|before|as of|up to)\s+(?:Thyra\s+)?v?\d+\.\d+", re.I),
    re.compile(r"\bv\d+\.\d+(?:\.\d+)?\b"),
    re.compile(r"\bused to\b", re.I),
    re.compile(r"\bissue\s+#?\d+", re.I),
    re.compile(r"(?<![\w/&])#\d{2,5}\b"),
)


@dataclass
class Unit:
    """One paragraph or list item of prose, and whether it sits in a box."""

    text: str
    advanced: bool


@dataclass
class _Walk:
    """Line-walk state: open fence, open box, and the unit being built."""

    fence: bool = False
    box_indent: int = -1
    box_advanced: bool = False
    lines: List[str] = field(default_factory=list)
    units: List[Unit] = field(default_factory=list)

    def flush(self) -> None:
        """Close the unit being built, if any."""
        if self.lines:
            self.units.append(Unit(" ".join(self.lines), self.box_advanced))
            self.lines = []


def _indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _leave_box_if_dedented(state: _Walk, line: str) -> None:
    if state.box_indent >= 0 and line.strip() and _indent(line) <= state.box_indent:
        state.flush()
        state.box_indent, state.box_advanced = -1, False


def _step(state: _Walk, line: str) -> None:
    """Feed one line of Markdown into the walk."""
    if _FENCE.match(line):
        state.flush()
        state.fence = not state.fence
        return
    if state.fence:
        return
    _leave_box_if_dedented(state, line)
    box = _ADMONITION.match(line)
    if box:
        state.flush()
        state.box_indent = len(box.group(1))
        state.box_advanced = box.group(2) == "advanced"
        return
    if not line.strip() or _STRUCTURE.match(line):
        state.flush()
        return
    item = _ITEM.match(line)
    if item:
        state.flush()
        line = line[item.end() :]
    state.lines.append(line.strip())


def prose_units(markdown: str) -> List[Unit]:
    """The prose of a page: paragraphs and list items, without code or tables."""
    state = _Walk()
    for line in markdown.splitlines():
        _step(state, line)
    state.flush()
    return state.units


def plain(text: str) -> str:
    """Prose as a reader sees it: code as one word, links as their text."""
    text = _INLINE_CODE.sub("code", text)
    text = _IMAGE.sub("", text)
    text = _LINK.sub(r"\1", text)
    text = _TAG.sub("", text)
    return _EMPHASIS.sub("", text)


def sentences(unit: Unit) -> Iterator[Tuple[int, str]]:
    """Each sentence of a unit with its word count."""
    for sentence in _SENTENCE_END.split(plain(unit.text)):
        words = [w for w in sentence.split() if _WORD.search(w)]
        if words:
            yield len(words), sentence


def main_text_words(units: List[Unit]) -> int:
    """Words outside the collapsed "Advanced:" boxes."""
    return sum(n for unit in units if not unit.advanced for n, _ in sentences(unit))


def _nav_pages(entry: Any) -> Iterator[str]:
    if isinstance(entry, str):
        yield entry
    elif isinstance(entry, list):
        for item in entry:
            yield from _nav_pages(item)
    elif isinstance(entry, dict):
        for value in entry.values():
            yield from _nav_pages(value)


def beginner_pages() -> List[str]:
    """The Markdown pages under Home and the beginner nav sections."""
    nav = yaml.safe_load(_MKDOCS.read_text(encoding="utf-8"))["nav"]
    pages: List[str] = []
    for entry in nav:
        for title, value in entry.items():
            if title == "Home" or title in BEGINNER_SECTIONS:
                pages.extend(p for p in _nav_pages(value) if p.endswith(".md"))
    return pages


def _units(page: str) -> List[Unit]:
    return prose_units((_DOCS / page).read_text(encoding="utf-8"))


def test_the_nav_still_has_beginner_pages() -> None:
    """Guard the discovery itself: renaming a section must not empty it."""
    pages = beginner_pages()
    assert "index.md" in pages and "install.md" in pages, pages


@pytest.mark.parametrize("page", beginner_pages())
def test_no_sentence_runs_past_the_limit(page: str) -> None:
    """Split a long sentence in two; one idea each reads faster."""
    long_ones = [
        f"{n} words: {text}"
        for unit in _units(page)
        for n, text in sentences(unit)
        if n > MAX_SENTENCE_WORDS
    ]
    assert not long_ones, "\n".join(long_ones)


@pytest.mark.parametrize("page", beginner_pages())
def test_sentences_stay_short_on_average(page: str) -> None:
    """The plain-language target is 15 to 20 words a sentence."""
    counts = [n for unit in _units(page) for n, _ in sentences(unit)]
    average = sum(counts) / len(counts)
    assert average <= MAX_AVERAGE_WORDS, f"{page}: {average:.1f} words a sentence"


@pytest.mark.parametrize("page", beginner_pages())
def test_main_text_stays_within_its_budget(page: str) -> None:
    """Growth is a decision: raise the budget here, in review, or cut."""
    assert page in WORD_BUDGETS, f"{page} has no entry in WORD_BUDGETS"
    words = main_text_words(_units(page))
    assert words <= WORD_BUDGETS[page], (
        f"{page}: {words} words of main text, budget {WORD_BUDGETS[page]}. "
        "Move detail into an 'Advanced:' box or the technical reference."
    )


@pytest.mark.parametrize("page", beginner_pages())
def test_no_version_history_or_issue_numbers(page: str) -> None:
    """History goes to the changelog; reasons to design-decisions.md."""
    found = [
        f"{pattern.pattern!r}: {plain(unit.text)}"
        for unit in _units(page)
        for pattern in _HISTORY
        if pattern.search(plain(unit.text))
    ]
    assert not found, "\n".join(found)
