"""Extract explicit key-value floral traits from structured species pages.

Regional floras and specialist sites often render a characteristic label followed
by one or more values.  This module reads the complete value block rather than
an arbitrary sentence.  Multiple reported states remain explicit instead of
being collapsed to whichever option happens to occur last on the page.
"""

from __future__ import annotations

import re

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


def _matches_label(label: str, trait_name: str) -> bool:
    return any(pattern.fullmatch(label) for pattern in LABELS.get(trait_name, ()))


def _known_label(line: str) -> bool:
    return any(
        pattern.fullmatch(line)
        for patterns in LABELS.values()
        for pattern in patterns
    )


def _looks_like_next_field(line: str) -> bool:
    """Recognize a short characteristic heading, not a value sentence."""

    if _known_label(line):
        return True
    if line.casefold() in {"na", "n/a", "no", "yes", "unknown"}:
        return False
    if len(line) > 80 or len(line.split()) > 9:
        return False
    if re.search(r"[.!?;:]$", line):
        return False
    return bool(re.fullmatch(r"[A-Z][A-Za-z ()/,'-]{2,79}", line))


def _field_blocks(text: str, trait_name: str) -> list[tuple[str, list[str]]]:
    lines = [_text(line) for line in text.splitlines() if _text(line)]
    blocks: list[tuple[str, list[str]]] = []
    for index, label in enumerate(lines):
        if not _matches_label(label, trait_name):
            continue
        values: list[str] = []
        for value in lines[index + 1 : index + 7]:
            if _looks_like_next_field(value):
                break
            values.append(value)
        if values:
            blocks.append((label, values))
    return blocks


def _colour(value: str) -> str:
    states = [state for pattern, state in COLOUR_MAP if pattern.search(value)]
    states = list(dict.fromkeys(states))
    if len(states) == 1:
        return states[0]
    if len(states) > 1:
        return "multicolored_variable"
    return ""


def _symmetry(value: str) -> str:
    radial = bool(
        re.search(
            r"radially\s+symmetr|actinomorph|two\s+or\s+more\s+ways\s+to\s+evenly\s+divide|"
            r"放射相称|辐射对称",
            value,
            re.IGNORECASE,
        )
    )
    bilateral = bool(
        re.search(
            r"bilaterally\s+symmetr|zygomorph|only\s+one\s+way\s+to\s+evenly\s+divide|"
            r"左右相称|两侧对称",
            value,
            re.IGNORECASE,
        )
    )
    if radial and bilateral:
        return "multistate_variable"
    if radial:
        return "actinomorphic"
    if bilateral:
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
    absent = bool(
        re.search(
            r"\b(?:no|none|absent|not)\b.{0,40}cleistogam|"
            r"there\s+are\s+no\s+cleistogamous\s+flowers|"
            r"閉鎖花.{0,10}(?:ない|無し|なし)",
            value,
            re.IGNORECASE,
        )
    )
    present = bool(
        re.search(
            r"(?:has|have|with|some|both).{0,35}cleistogam|"
            r"cleistogamous\s+flowers?.{0,20}(?:present|occur)|"
            r"閉鎖花.{0,10}(?:ある|あり|つける)",
            value,
            re.IGNORECASE,
        )
    )
    obligate = bool(
        re.search(
            r"obligate|only\s+cleistogamous|閉鎖花のみ",
            value,
            re.IGNORECASE,
        )
    )
    facultative = bool(
        re.search(
            r"facultative|chasmogamous\s+and\s+cleistogamous|"
            r"some\s+cleistogamous|開放花と閉鎖花",
            value,
            re.IGNORECASE,
        )
    )
    if present and absent:
        return "multistate_variable"
    if obligate:
        return "obligate"
    if facultative or present:
        return "facultative"
    if absent:
        return "absent"
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
    for label, values in _field_blocks(page.text, trait_name):
        value = "; ".join(values)
        normalized = extractor(value)
        if not normalized:
            continue
        excerpt = f"{label}: {value}"
        rows.append((value, normalized, excerpt))
    return list(dict.fromkeys(rows))
