"""Extract explicit key-value floral traits from structured species pages.

Many regional floras and specialist sites render a characteristic name on one
line and its value on the next.  The sentence-rule extractor intentionally does
not join arbitrary neighbouring text, so this module handles only a small set of
explicit, auditable field labels.  Each result retains the exact label and value
as the supporting excerpt.
"""

from __future__ import annotations

import re
from collections.abc import Sequence

from island_v2.open_web_evidence import Page


COLOUR_MAP: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"\b(?:white|cream)\b|白", re.IGNORECASE), "white"),
    (re.compile(r"\b(?:yellow|orange)\b|黄|橙", re.IGNORECASE), "yellow_orange"),
    (re.compile(r"\b(?:red|pink|rose)\b|赤|紅|桃", re.IGNORECASE), "red_pink"),
    (re.compile(r"\b(?:blue|purple|violet)\b|青|紫", re.IGNORECASE), "blue_purple"),
    (
        re.compile(r"\b(?:green|brown|inconspicuous)\b|緑|褐|目立たない", re.IGNORECASE),
        "green_brown_inconspicuous",
    ),
)

LABELS: dict[str, tuple[re.Pattern[str], ...]] = {
    "flower_primary_color": (
        re.compile(r"^(?:flower\s+)?petal\s+colou?r$", re.IGNORECASE),
        re.compile(r"^flower\s+colou?r$", re.IGNORECASE),
        re.compile(r"^corolla\s+colou?r$", re.IGNORECASE),
        re.compile(r"^花(?:冠|弁)?(?:の)?色$"),
    ),
    "floral_symmetry": (
        re.compile(r"^flower\s+symmetry$", re.IGNORECASE),
        re.compile(r"^corolla\s+symmetry$", re.IGNORECASE),
        re.compile(r"^花(?:冠)?(?:の)?相称$"),
    ),
    "flower_size_class": (
        re.compile(r"^flower\s+(?:diameter|length|size)$", re.IGNORECASE),
        re.compile(r"^corolla\s+(?:diameter|length|size)$", re.IGNORECASE),
        re.compile(r"^花(?:冠)?(?:の)?(?:直径|長さ|大きさ)$"),
    ),
    "cleistogamy": (
        re.compile(r"^cleistogamous\s+flowers?$", re.IGNORECASE),
        re.compile(r"^cleistogamy$", re.IGNORECASE),
        re.compile(r"^閉鎖花$"),
    ),
}

MEASUREMENT = re.compile(
    r"(?P<low>\d+(?:\.\d+)?)\s*(?:[-–—~〜～]\s*(?P<high>\d+(?:\.\d+)?))?\s*"
    r"(?P<unit>mm|㎜|cm|㎝)",
    re.IGNORECASE,
)


def _text(value: object) -> str:
    return " ".join(str(value or "").replace("\xa0", " ").split())


def _pairs(text: str) -> Sequence[tuple[str, str]]:
    lines = [_text(line) for line in text.splitlines() if _text(line)]
    return tuple(zip(lines, lines[1:]))


def _matches_label(label: str, trait_name: str) -> bool:
    return any(pattern.fullmatch(label) for pattern in LABELS.get(trait_name, ()))


def _colour(value: str) -> str:
    states = [state for pattern, state in COLOUR_MAP if pattern.search(value)]
    states = list(dict.fromkeys(states))
    if len(states) == 1:
        return states[0]
    if len(states) > 1:
        return "multicolored_variable"
    return ""


def _symmetry(value: str) -> str:
    if re.search(
        r"radially\s+symmetr|actinomorph|two\s+or\s+more\s+ways\s+to\s+evenly\s+divide|"
        r"放射相称|辐射对称",
        value,
        re.IGNORECASE,
    ):
        return "actinomorphic"
    if re.search(
        r"bilaterally\s+symmetr|zygomorph|only\s+one\s+way\s+to\s+evenly\s+divide|"
        r"左右相称|两侧对称",
        value,
        re.IGNORECASE,
    ):
        return "zygomorphic"
    return ""


def _size(value: str) -> str:
    match = MEASUREMENT.search(value)
    if not match:
        return ""
    low = float(match.group("low"))
    high = float(match.group("high") or match.group("low"))
    midpoint = (low + high) / 2
    if match.group("unit").casefold() in {"cm", "㎝"}:
        midpoint *= 10
    if midpoint <= 5:
        return "very_small"
    if midpoint <= 15:
        return "small"
    if midpoint <= 30:
        return "medium"
    if midpoint <= 60:
        return "large"
    return "very_large"


def _cleistogamy(value: str) -> str:
    if re.search(
        r"\b(?:no|none|absent|not)\b.{0,30}cleistogam|"
        r"there\s+are\s+no\s+cleistogamous\s+flowers|閉鎖花.{0,10}(?:ない|無し|なし)",
        value,
        re.IGNORECASE,
    ):
        return "absent"
    if re.search(r"facultative|chasmogamous\s+and\s+cleistogamous|開放花と閉鎖花", value, re.IGNORECASE):
        return "facultative"
    if re.search(r"obligate|only\s+cleistogamous|閉鎖花のみ", value, re.IGNORECASE):
        return "obligate"
    return ""


def extract_structured_characteristics(
    page: Page,
    *,
    trait_name: str,
) -> list[tuple[str, str, str]]:
    """Return raw value, normalized value, and exact label-value excerpt."""

    extractors = {
        "flower_primary_color": _colour,
        "floral_symmetry": _symmetry,
        "flower_size_class": _size,
        "cleistogamy": _cleistogamy,
    }
    extractor = extractors.get(trait_name)
    if extractor is None:
        return []
    rows: list[tuple[str, str, str]] = []
    for label, value in _pairs(page.text):
        if not _matches_label(label, trait_name):
            continue
        normalized = extractor(value)
        if not normalized:
            continue
        excerpt = f"{label}: {value}"
        rows.append((value, normalized, excerpt))
    return list(dict.fromkeys(rows))
