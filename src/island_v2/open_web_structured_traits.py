"""Extract explicit key-value floral traits from structured species pages.

Regional floras and specialist sites often render a characteristic label followed
by one or more values.  This module reads the complete value block rather than
an arbitrary sentence.  Multiple reported states remain explicit instead of
being collapsed to whichever option happens to occur last on the page.
"""

from __future__ import annotations

import re

from island_v2.open_web_evidence import Page, primary_treatment_text

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
    "floral_form": (
        re.compile(r"^flower\s+shape$", re.IGNORECASE),
        re.compile(r"^corolla\s+(?:form|shape)$", re.IGNORECASE),
    ),
    "flower_size_class": (
        re.compile(r"^flower\s+(?:diameter|length|size)$", re.IGNORECASE),
        re.compile(r"^corolla\s+(?:diameter|length|size)$", re.IGNORECASE),
        re.compile(r"^花(?:冠)?(?:の)?(?:直径|長さ|大きさ)$"),
    ),
    "inflorescence_display": (
        re.compile(r"^flower\s+inflorescence$", re.IGNORECASE),
        re.compile(r"^inflorescence(?:\s+type)?$", re.IGNORECASE),
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

NARRATIVE_FRAGMENT_SPLIT = re.compile(
    r"(?<=[;!?])\s+|(?<=[a-z\]\)])\.\s+(?=[A-Z])|\n+"
)
PRIMARY_FLORAL_ORGAN = re.compile(
    r"\b(?:flowers?|corollas?|petals?|tepals?|perianths?)\b", re.IGNORECASE
)
NON_PRIMARY_FLORAL_ORGAN = re.compile(
    r"\b(?:caly(?:x|ces)|sepals?|bracts?|anthers?|stamens?|styles?|stigmas?|"
    r"ovaries|ovary)\b",
    re.IGNORECASE,
)
FLOWER_SIZE_BRIDGE_BLOCKER = re.compile(
    r"\b(?:inflorescences?|racemes?|spikes?|panicles?|umbels?|corymbs?|"
    r"capitula|heads?|peduncles?|pedicels?|caly(?:x|ces)|sepals?|bracts?|lobes?)\b",
    re.IGNORECASE,
)
INFLORESCENCE_SUBJECT = re.compile(
    r"\binflorescences?\b|\bflowers?\b.{0,100}"
    r"\b(?:solitary|few[- ]flowered|raceme|spike|panicle|umbel|corymb|"
    r"capitulum|flower head)\b",
    re.IGNORECASE,
)
TUBE_SUBJECT = re.compile(r"\b(?:corolla\s+|floral\s+)?tube\b", re.IGNORECASE)
DESCRIPTION_END = re.compile(
    r"\n(?:Flowering|Distribution(?: and occurrence)?|Habitat|Phenology|Notes?)\s*:",
    re.IGNORECASE,
)


def _text(value: object) -> str:
    return " ".join(str(value or "").replace("\xa0", " ").split())


def _matches_label(label: str, trait_name: str) -> bool:
    return any(pattern.fullmatch(label.rstrip(":")) for pattern in LABELS.get(trait_name, ()))


def _known_label(line: str) -> bool:
    return any(
        pattern.fullmatch(line.rstrip(":"))
        for patterns in LABELS.values()
        for pattern in patterns
    )


def _looks_like_next_field(line: str) -> bool:
    """Recognize a short characteristic heading, not a value sentence."""

    if _known_label(line):
        return True
    if line.endswith(":"):
        return True
    if line.casefold() in {"na", "n/a", "no", "yes", "unknown"}:
        return False
    return False


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


def _form(value: str) -> str:
    mappings = (
        (r"\bbell\b|campanulat", "bell_campanulate"),
        (r"\btubular\b", "tubular"),
        (r"salver", "salverform"),
        (r"\bfunnel\b|\btrumpet\b", "funnel_trumpet"),
        (r"\burn\b|urceolat", "urn_urceolate"),
        (r"papilion|\bpea(?:-like)?\b", "papilionaceous"),
        (r"\blipped\b|bilabiat|two[- ]lipped", "bilabiate"),
        (r"\bspur(?:red)?\b", "spurred"),
        (r"\bcup\b|\bsaucer\b|\bstar\b|rotate", "open_radial"),
    )
    states = [state for pattern, state in mappings if re.search(pattern, value, re.IGNORECASE)]
    states = list(dict.fromkeys(states))
    return states[0] if len(states) == 1 else "other_described" if states else ""


def _inflorescence(value: str) -> str:
    mappings = (
        (r"\bsolitary\b|single flower", "solitary"),
        (r"few[- ]flower", "few_flowered"),
        (r"\braceme\b|\bspike\b|\bpanicle\b|\bcatkin\b", "raceme_spike_panicle"),
        (r"\bumbel\b|\bcorymb\b|\bcyme\b", "umbel_corymb"),
        (r"\bhead\b|capitulum", "composite_display"),
        (r"brush|puff", "brush_puff_display"),
    )
    states = [state for pattern, state in mappings if re.search(pattern, value, re.IGNORECASE)]
    states = list(dict.fromkeys(states))
    return states[0] if len(states) == 1 else "other_described" if states else ""


def _size(value: str) -> str:
    if re.search(r"<\s*\d", value):
        return ""
    inch = re.search(
        r"(?P<low>\d+(?:\.\d+)?)\s*(?:[-–—]\s*(?P<high>\d+(?:\.\d+)?))?\s*"
        r"(?:inches|inch|in\.)\b",
        value,
        re.IGNORECASE,
    )
    match = inch or MEASUREMENT.search(value)
    if not match:
        return ""
    low = float(match.group("low"))
    high = float(match.group("high") or match.group("low"))
    midpoint = (low + high) / 2
    if inch:
        midpoint *= 25.4
        if midpoint <= 5:
            return "very_small"
        if midpoint <= 15:
            return "small"
        if midpoint <= 30:
            return "medium"
        if midpoint <= 60:
            return "large"
        return "very_large"
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


def _narrative_fragments(
    page: Page, *, treatment_end_pattern: str = ""
) -> list[str]:
    text = primary_treatment_text(
        page,
        treatment_end_pattern=treatment_end_pattern,
    )
    description = re.search(r"(?:^|\n)Description\s*:\s*", text, re.IGNORECASE)
    if description:
        text = text[description.end() :]
    end = DESCRIPTION_END.search(text)
    if end:
        text = text[: end.start()]
    return [
        _text(fragment)
        for fragment in NARRATIVE_FRAGMENT_SPLIT.split(text)
        if len(_text(fragment)) >= 8
    ]


def _tube_depth(value: str) -> str:
    subject = TUBE_SUBJECT.search(value)
    if not subject:
        return ""
    match = MEASUREMENT.search(value, subject.end())
    if not match or match.start() - subject.end() > 100:
        return ""
    low = float(match.group("low"))
    high = float(match.group("high") or match.group("low"))
    midpoint = (low + high) / 2
    if match.group("unit").casefold() in {"cm", "㎝"}:
        midpoint *= 10
    return "shallow" if midpoint <= 5 else "intermediate" if midpoint <= 20 else "deep"


def _primary_organ_subject(fragment: str) -> re.Match[str] | None:
    subject = PRIMARY_FLORAL_ORGAN.search(fragment)
    if subject and not NON_PRIMARY_FLORAL_ORGAN.search(fragment[: subject.start()]):
        return subject
    return None


def extract_botanical_description(
    page: Page,
    *,
    trait_name: str,
    treatment_end_pattern: str = "",
) -> list[tuple[str, str, str]]:
    """Extract organ-scoped traits from narrative flora descriptions.

    Botanical abbreviations such as ``c.`` must not split the organ from its
    measurement or colour.  Every returned value therefore comes from one
    verbatim description fragment with an explicit floral subject.
    """

    rows: list[tuple[str, str, str]] = []
    for fragment in _narrative_fragments(
        page, treatment_end_pattern=treatment_end_pattern
    ):
        normalized = ""
        raw_value = fragment
        if trait_name == "flower_primary_color":
            if _primary_organ_subject(fragment) and not NON_PRIMARY_FLORAL_ORGAN.search(
                fragment
            ):
                normalized = _colour(fragment)
        elif trait_name == "floral_form":
            if _primary_organ_subject(fragment):
                normalized = _form(fragment)
        elif trait_name == "floral_symmetry":
            if _primary_organ_subject(fragment):
                normalized = _symmetry(fragment)
        elif trait_name == "flower_size_class":
            subject = _primary_organ_subject(fragment)
            if subject:
                measurement = MEASUREMENT.search(fragment, subject.end())
                bridge = fragment[subject.end() : measurement.start()] if measurement else ""
                if (
                    measurement
                    and measurement.start() - subject.end() <= 120
                    and not FLOWER_SIZE_BRIDGE_BLOCKER.search(bridge)
                ):
                    normalized = _size(fragment[subject.start() :])
                    raw_value = measurement.group(0)
        elif trait_name == "tube_depth_class":
            normalized = _tube_depth(fragment)
        elif trait_name == "inflorescence_display":
            if INFLORESCENCE_SUBJECT.search(fragment):
                normalized = _inflorescence(fragment)
        elif trait_name == "cleistogamy" and re.search(
            r"\bcleistogam", fragment, re.IGNORECASE
        ):
            normalized = _cleistogamy(fragment)
        if normalized and normalized != "other_described":
            rows.append((raw_value, normalized, fragment))
    return list(dict.fromkeys(rows))


def extract_structured_characteristics(
    page: Page,
    *,
    trait_name: str,
) -> list[tuple[str, str, str]]:
    """Return raw value, normalized value, and exact label-value excerpt."""

    extractors = {
        "flower_primary_color": _colour,
        "floral_symmetry": _symmetry,
        "floral_form": _form,
        "flower_size_class": _size,
        "inflorescence_display": _inflorescence,
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
