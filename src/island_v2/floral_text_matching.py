"""Match floral colour statements in text, identically in every lane.

Two acquisition lanes mine floral colour from free text: protologues scanned by
BHL, and herbarium specimen labels carried in GBIF occurrences. They read very
different prose but they answer the same question, feed the same candidate
ledger and are read by the same reviewer, so a colour word that is evidence in
one lane and invisible in the other is not two policies -- it is one policy
applied inconsistently.

That is measured, not hypothetical. The label lane was written first with thirty
English colour words and plain substring matching; the protologue lane later
grew word boundaries, sentence walls, Unicode folding and a 255-term
nine-language vocabulary. By the time both were wired to the same endemic queue
they disagreed on ordinary sentences, and the label lane was the one producing
wrong answers:

    "Flowers on reduced peduncles"                    -> red_pink, matching
                                                        "red" inside "reduced"
    "Leaves coriaceous, beneath brown. Flowers seen." -> green_brown, a leaf
                                                        clause reaching across
                                                        a full stop

Neither is missing evidence. Both are wrong evidence, which is worse, because a
reviewer reads a populated cell as an observation. This module is the single
implementation both lanes call, so the rules cannot drift again. What stays in
each lane's own config is what genuinely differs: how wide the organ window is,
and what voids a statement.

The vocabulary itself lives in ``config/floral_vocabulary.yml``, which both lane
configs inherit through ``extends``.
"""

from __future__ import annotations

import re
import unicodedata
from pathlib import Path
from typing import Any, Mapping

import typer
import yaml

ACCEPTED = "candidate"
REJECT_NO_TEXT = "no_source_text"
REJECT_NEGATED = "statement_negated_or_uncertain"
REJECT_NO_COLOUR = "no_colour_term"
REJECT_NO_ORGAN = "colour_not_anchored_to_floral_organ"
REJECT_COMPETING = "colour_belongs_to_non_floral_organ"

# Japanese voicing marks. Stripping these along with the accents would fold
# distinct kana together, so they survive folding and are recomposed.
VOICING_MARKS = ("゙", "゚")


def load_config(path: Path, required: set[str]) -> dict[str, Any]:
    """Load a lane config, resolving the shared vocabulary it inherits.

    A lane config names its base in ``extends``, resolved relative to its own
    directory. Keys the lane declares itself win, which is how the protologue
    lane keeps a 60-character organ window against the label lane's 40 without
    either of them owning a private copy of the vocabulary.
    """
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise typer.BadParameter(f"{path} is not a mapping")

    inherited = config.pop("extends", None)
    if inherited:
        base_path = path.parent / str(inherited)
        base = yaml.safe_load(base_path.read_text(encoding="utf-8"))
        if not isinstance(base, dict):
            raise typer.BadParameter(f"{base_path} is not a mapping")
        base.update(config)
        config = base

    if not required.issubset(config):
        raise typer.BadParameter(f"config must contain {sorted(required)}")
    return config


def fold(text: str) -> tuple[str, list[int]]:
    """Fold text for matching and keep a map back to the original offsets.

    Folding lowercases, strips combining marks so ``blanchâtre`` matches
    ``blanchatre``, and expands eszett so ``weiß`` matches ``weiss``. Because
    eszett expands to two characters the offsets shift, so each folded character
    carries the index of the original character it came from -- that is what
    lets an accepted match be quoted verbatim rather than in folded form.

    Japanese voicing is the one mark that is not an accent. Dakuten and
    handakuten do not decorate a kana, they select a different one, so
    discarding them turns ピンク into ヒンク and collapses ハラ, バラ and パラ onto
    a single form. The marks are therefore kept and recomposed, which also folds
    halfwidth ﾋﾟ onto precomposed ピ.
    """
    folded: list[str] = []
    origin: list[int] = []
    for index, char in enumerate(text):
        if char in ("ß", "ẞ"):  # eszett, capital eszett
            expanded = "ss"
        else:
            decomposed = unicodedata.normalize("NFKD", char)
            expanded = "".join(
                c
                for c in decomposed
                if c in VOICING_MARKS or not unicodedata.combining(c)
            ).lower()
        for produced in expanded:
            if produced in VOICING_MARKS and folded:
                recomposed = unicodedata.normalize("NFC", folded[-1] + produced)
                if len(recomposed) == 1:
                    folded[-1] = recomposed
                    continue
            folded.append(produced)
            origin.append(index)
    return "".join(folded), origin


# Scripts that write words without spaces. A term in one of these is matched as
# a plain substring, because there is no boundary to anchor to.
SPACELESS_SCRIPT = re.compile(
    r"[぀-ヿ㐀-䶿一-鿿豈-﫿가-힯]"
)


def boundary_spans(text: str, term: str) -> list[tuple[int, int]]:
    """Every whole-word occurrence of ``term`` in already-folded ``text``.

    The guard is chosen by the script the term is written in, because the two
    families need opposite treatment.

    An alphabetic term must not match inside a longer word -- that rule is what
    stops ``rotundifolia`` being read as German ``rot`` and ``reduced`` as
    English ``red``. The guard has to cover every letter, not just ASCII: with a
    Latin-only class, Russian ``белый`` matches inside ``белыйцветок`` and
    Cyrillic loses the protection Latin has.

    A Han, Kana or Hangul term needs the opposite. Those scripts do not separate
    words, so ``白色`` legitimately sits flush against its neighbours in
    ``花は白色`` and any boundary guard would reject every real match.
    """
    if not term:
        return []
    if SPACELESS_SCRIPT.search(term):
        pattern = re.compile(re.escape(term))
    else:
        pattern = re.compile(rf"(?<![^\W_]){re.escape(term)}(?![^\W_])", re.UNICODE)
    return [(m.start(), m.end()) for m in pattern.finditer(text)]


SENTENCE_DELIMITERS = ".;\n。！？"

# A period inside one of these does not end a sentence. Written folded, and
# without the trailing period, which is what gets masked.
DEFAULT_ABBREVIATIONS = ("fl", "fr", "var", "subsp", "ssp", "cf", "aff", "sp")

# Reviewers quote with an ellipsis, which is elision rather than a full stop.
_ELLIPSIS = re.compile(r"\.{2,}|…")


def mask_non_terminal_periods(text: str, abbreviations: tuple[str, ...]) -> str:
    """Blank the periods that elide or abbreviate, preserving every offset.

    Both replacements are length-preserving so that the spans this function
    feeds still index the caller's text.
    """
    masked = _ELLIPSIS.sub(lambda m: " " * len(m.group()), text)
    for abbreviation in abbreviations:
        masked = re.sub(
            rf"(?<![^\W_])({re.escape(abbreviation)})\.",
            r"\1 ",
            masked,
            flags=re.UNICODE,
        )
    return masked


def sentence_spans(
    text: str, abbreviations: tuple[str, ...] = DEFAULT_ABBREVIATIONS
) -> list[tuple[int, int]]:
    """Split text into the segments an organ description may not escape.

    Botanical descriptions are organised one organ per sentence or
    semicolon-delimited clause: "Folia coriacea, subtus fusca. Flores albi."
    Searching for the nearest organ across that boundary attaches the leaf's
    brown to the flower, which is the exact failure the WFO ledger measured at
    38.6% of its rejections. A period between them is a wall.

    Two kinds of period are not that wall, and both were measured costing real
    evidence in the reviewed ledger rather than merely being possible:

    An ellipsis is elision. A reviewer quoting ``petala ... albida`` means one
    statement, but split on the dots the colour loses its organ and the row is
    rejected as unanchored.

    An abbreviation's period is part of the word. ``fl. pink`` is a floral
    colour statement; split on the period, ``pink`` has no organ in its segment.
    Specimen labels are written almost entirely in such abbreviations, so this
    matters more on a label than in protologue prose.
    """
    prepared = mask_non_terminal_periods(text, abbreviations)
    spans: list[tuple[int, int]] = []
    start = 0
    for index, char in enumerate(prepared):
        if char in SENTENCE_DELIMITERS:
            if index > start:
                spans.append((start, index))
            start = index + 1
    if start < len(prepared):
        spans.append((start, len(prepared)))
    return spans


def containing_span(spans: list[tuple[int, int]], position: int) -> tuple[int, int]:
    for start, end in spans:
        if start <= position < end:
            return start, end
    return position, position


def nearest_organ(
    span: tuple[int, int],
    organs: list[tuple[int, int, str]],
    bounds: tuple[int, int],
    window: int,
) -> str:
    """Return "floral", "competing" or "" for the organ closest to a colour word.

    Only organ terms inside ``bounds`` -- the colour word's own sentence -- are
    considered. Within it, distance is measured to the nearest whole-word organ
    on either side, and when nothing separates them ("flores rubri fructus", one
    character each way) the competing organ wins: the conservative reading is
    the one the WFO rejection ledger says matters.
    """
    start, end = span
    low, high = bounds
    best_kind = ""
    best_distance = window + 1
    for organ_start, organ_end, kind in organs:
        if organ_start < low or organ_end > high:
            continue
        if organ_start <= start and end <= organ_end:
            distance = 0
        else:
            distance = min(abs(start - organ_end), abs(organ_start - end))
        if distance <= window and (
            distance < best_distance
            or (distance == best_distance and kind == "competing")
        ):
            best_distance = distance
            best_kind = kind
    return best_kind


def sentence_around(text: str, start: int, end: int, limit: int = 500) -> str:
    """The verbatim sentence containing an accepted match, capped at ``limit``."""
    marks = tuple(SENTENCE_DELIMITERS)
    left = max((text.rfind(mark, 0, start) for mark in marks), default=-1)
    right_candidates = [
        position
        for position in (text.find(mark, end) for mark in marks)
        if position != -1
    ]
    right = min(right_candidates) + 1 if right_candidates else len(text)
    return " ".join(text[left + 1 : right].split())[:limit]


# Between two colour words, whitespace, hyphens and the declared modifier words
# mean they qualify one another. Anything else -- a comma, "demum", "and"
# -- separates statements.
_COMPOUND_PUNCTUATION = re.compile(r"[\s\-]+")


def is_compound_gap(gap: str, gap_terms: tuple[str, ...]) -> bool:
    """Whether what sits between two colour words merely qualifies them."""
    remainder = gap
    for term in sorted(gap_terms, key=len, reverse=True):
        remainder = remainder.replace(term, " ")
    return not _COMPOUND_PUNCTUATION.sub("", remainder)


def expand_colour_terms(config: Mapping[str, Any]) -> dict[str, str]:
    """The full folded colour vocabulary, with Latin declensions applied.

    Latin adjectives agree with the organ noun, so which form appears is fixed
    by the sentence -- "corolla alba" but "flores albi" but "floribus albis" --
    and protologues use the plural at least as often as the dictionary form.
    Enumerating only the singular would drop most real statements, so stems are
    declined here instead. A generated non-word never matches, which makes
    over-generation free and under-generation the only real cost.
    """
    terms: dict[str, str] = {}
    for stems_key, endings_key in (
        ("latin_adjective_stems", "latin_adjective_endings"),
        ("latin_third_declension_stems", "latin_third_declension_endings"),
    ):
        endings = [str(e) for e in (config.get(endings_key) or [])]
        for stem, value in (config.get(stems_key) or {}).items():
            for ending in endings:
                terms[fold(f"{stem}{ending}")[0]] = str(value)
    # Literal entries win over generated ones: an irregular listed by hand is
    # there precisely because the stem rule cannot produce it correctly.
    for term, value in config["colour_terms"].items():
        terms[fold(str(term))[0]] = str(value)
    return terms


def plain_colour_values(config: Mapping[str, Any]) -> set[str]:
    """The ontology values the analysis counts as plain.

    Declared in configuration rather than imported, so that acquisition does not
    pull in the analysis dependency chain. A test asserts the two agree.
    """
    return {str(value) for value in config["plain_colour_values"]}


def binary_plain_class(matched_terms: str, config: Mapping[str, Any]) -> str:
    """Return "plain", "nonplain" or "unresolved" for the terms that matched.

    This is the level the model actually consumes, and it survives cases the
    five-value ontology cannot decide. The reviewer's own frozen call on
    "Corolla pale reddish purple" was that red-pink versus blue-purple is not
    uniquely resolvable while non-plain is certain, and the same split resolves
    it here: both candidate values sit on the non-plain side of the line.

    A description straddling the line -- "greenish yellow", plain and non-plain
    at once -- stays unresolved, which is what the reviewer recorded for it.
    """
    if not matched_terms:
        return "unresolved"
    colours = expand_colour_terms(config)
    plain = plain_colour_values(config)
    sides = {
        ("plain" if colours[term] in plain else "nonplain")
        for term in matched_terms.split("+")
        if term in colours
    }
    return sides.pop() if len(sides) == 1 else "unresolved"


def _organ_positions(text: str, config: Mapping[str, Any]) -> list[tuple[int, int, str]]:
    """Every whole-word organ occurrence, found once for the whole document.

    An OCR'd page against two hundred organ terms and three hundred colour forms
    is otherwise a five-figure number of regex scans.
    """
    organs: list[tuple[int, int, str]] = []
    for kind, key in (
        ("competing", "competing_organ_terms"),
        ("floral", "floral_organ_terms"),
    ):
        for raw in config[key]:
            for start, end in boundary_spans(text, fold(str(raw))[0]):
                organs.append((start, end, kind))
    return organs


def extract_floral_colour(
    description: str,
    config: Mapping[str, Any],
    empty_text_outcome: str = REJECT_NO_TEXT,
) -> tuple[str, str, str, str]:
    """Return (outcome, normalized_value, matched_terms, verbatim_quote).

    ``empty_text_outcome`` is the one thing the lanes name differently, because
    a protologue with no OCR text and an occurrence with no label text are
    different failures worth counting apart.
    """
    if not description.strip():
        return empty_text_outcome, "", "", ""

    text, origin = fold(description)

    for marker in config["negation_markers"]:
        folded_marker, _ = fold(str(marker))
        if folded_marker and folded_marker in text:
            return REJECT_NEGATED, "", "", ""

    colours = expand_colour_terms(config)
    window = int(config["organ_proximity_chars"])
    spans = sentence_spans(
        text, tuple(config.get("abbreviations") or DEFAULT_ABBREVIATIONS)
    )
    organs = _organ_positions(text, config)

    accepted: list[tuple[tuple[int, int], str, str]] = []
    quote = ""
    saw_colour = False
    saw_competing = False

    for term in sorted(colours, key=len, reverse=True):
        for span in boundary_spans(text, term):
            saw_colour = True
            kind = nearest_organ(span, organs, containing_span(spans, span[0]), window)
            if kind == "competing":
                saw_competing = True
                continue
            if kind != "floral":
                continue
            # A shorter term inside one already accepted is the same word:
            # "yellow" inside "yellowish" is not a second colour statement.
            if any(
                start <= span[0] and span[1] <= end for (start, end), _, _ in accepted
            ):
                continue
            accepted.append((span, colours[term], term))
            if not quote:
                quote = sentence_around(
                    description, origin[span[0]], origin[span[1] - 1] + 1
                )

    if not saw_colour:
        return REJECT_NO_COLOUR, "", "", ""
    if not accepted:
        return (REJECT_COMPETING if saw_competing else REJECT_NO_ORGAN), "", "", ""

    accepted.sort(key=lambda hit: hit[0])
    matched = "+".join(sorted(term for _, _, term in accepted))
    values = list(dict.fromkeys(value for _, value, _ in accepted))
    if len(values) == 1:
        return ACCEPTED, values[0], matched, quote

    # Several ontology values. Whether that is one hue or several colours turns
    # on what sits between the words: "pale reddish purple" is a single hue
    # written with a modifier, while "alba demum rosea" is a corolla that
    # changes. Collapsing both to "variable" asserts a variability the first one
    # does not have, so they are separated here.
    gap_terms = tuple(
        fold(str(term))[0] for term in (config.get("compound_gap_terms") or ())
    )
    compound = all(
        is_compound_gap(text[left[0][1] : right[0][0]], gap_terms)
        for left, right in zip(accepted, accepted[1:])
    )
    if compound:
        return (
            ACCEPTED,
            str(config.get("compound_hue_result", "ontology_unresolved")),
            matched,
            quote,
        )
    return (
        ACCEPTED,
        str(config.get("multiple_value_result", "multicolored_variable")),
        matched,
        quote,
    )
