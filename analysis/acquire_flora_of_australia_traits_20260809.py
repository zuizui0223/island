"""Acquire six visible floral traits from the official Flora of Australia.

The acquisition unit is the complete Flora of Australia species-profile index,
not one search query per island species.  The public collection index is frozen
once, exact accepted species are selected locally, and each matching profile is
cached independently so interrupted runs resume without refetching completed
treatments.  Only the species-level ``Description`` field is interpreted.

Promotion remains separate: this module writes common-schema candidate evidence
and a deterministic 200-row audit queue.  A reviewed audit must pass the shared
source-package gate before any row enters formal coverage.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import re
import time
import unicodedata
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.trait_measurements import derive_class, load_rules

OPUS_ID = "0ded7a77-9efb-4684-8df0-48cbb1933684"
OPUS_SHORT_NAME = "foa"
INDEX_URL = "https://profiles.ala.org.au/profile/search/taxon/name"
PROFILE_URL = "https://profiles.ala.org.au/opus/{opus_id}/profile/{profile_id}/json"
SOURCE_LANDING_URL = "https://profiles.ala.org.au/opus/foa"
SOURCE_RUN_ID = "flora-of-australia-source-scale-20260809"
SOURCE_ARTIFACT = "flora-of-australia-public-profile-api"
SOURCE_FILE = "flora_of_australia_evidence.csv.gz"
ACCEPTANCE_CONTRACT = "official_flora_exact_species_description_quote_audit_v1"
USER_AGENT = "island-trait-research/1.0 (github.com/zuizui0223/island)"
AUDIT_SIZE = 200
AUDIT_SEED = "flora-of-australia-audit-20260809"

AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_symmetry": "floral_structural_complexity",
    "floral_form": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
}

COLOUR_PATTERNS = (
    (r"\b(?:white|whitish|cream(?:y)?|ivory)\b", "white"),
    (
        (
            r"\b(?:yellow(?:ish)?|orange(?:ish)?|gold(?:en)?|ochre|stramineous|"
            r"straw[- ]colou?red)\b"
        ),
        "yellow_orange",
    ),
    (
        (
            r"\b(?:red(?:dish)?|pink(?:ish)?|rose(?:ate)?|scarlet|crimson|"
            r"magenta|maroon|salmon|coral)\b"
        ),
        "red_pink",
    ),
    (
        (
            r"\b(?:blue|bluish|purple|purplish|violet|violaceous|lilac|mauve|"
            r"lavender)\b"
        ),
        "blue_purple",
    ),
    (
        (
            r"\b(?:green(?:ish)?|brown(?:ish)?|grey(?:ish)?|gray(?:ish)?|fuscous|"
            r"ferruginous|ferrugineous|inconspicuous)\b"
        ),
        "green_brown_inconspicuous",
    ),
    (r"\bblack(?:ish)?\b", "other_described"),
)

FORM_PATTERNS = (
    (r"\bcampanulat\w*|\bbell[- ]shaped\b", "bell_campanulate"),
    (r"\btubular\b|\btubulose\b", "tubular"),
    (r"\bsalverform\w*|\bhypocrateriform\w*", "salverform"),
    (r"\binfundibuliform\w*|\bfunnel[- ]shaped\b|\btrumpet[- ]shaped\b", "funnel_trumpet"),
    (r"\burn[- ]shaped\b|\burceolat\w*", "urn_urceolate"),
    (r"\bpapilion\w*", "papilionaceous"),
    (r"\bbilabiat\w*|\btwo[- ]lipped\b", "bilabiate"),
    (r"\bspurred\b|\bcalcarate\b", "spurred"),
    (r"\brotate\b|\bstar[- ]shaped\b|\bsaucer[- ]shaped\b", "open_radial"),
)

INFLORESCENCE_PATTERNS = (
    (r"\bsolitary\b", "solitary"),
    (r"\bfew[- ]flowered\b|\bflowers?\s+(?:in\s+)?pairs?\b", "few_flowered"),
    (r"\b(?:racemes?|spikes?|panicles?|thyrses?|catkins?)\b", "raceme_spike_panicle"),
    (r"\b(?:umbels?|corymbs?|cymes?|cincinni)\b", "umbel_corymb"),
    (r"\b(?:capitula|capitulum|flower[- ]heads?)\b", "composite_display"),
    (r"\b(?:brush|puff)[- ](?:like|shaped)\b", "brush_puff_display"),
)

TARGET_SUBJECT = re.compile(
    r"\b(?:flowers?|corollas?|petals?|perianths?|tepals?|labellums?|labella)\b",
    re.IGNORECASE,
)
WHOLE_FLOWER_SUBJECT = re.compile(
    r"\b(?:flowers?|corollas?|perianths?)\b",
    re.IGNORECASE,
)
NON_TARGET_SUBJECT = re.compile(
    r"\b(?:plants?|fruits?|berries|leaves|leaf|foliage|stems?|branches?|bracts?|"
    r"bracteoles?|calyces|calyx|sepals?|anthers?|filaments?|styles?|stigmas?|"
    r"ovaries|ovary|hypanthia|hypanthium|capsules?|rays?|seeds?|hairs?|"
    r"indumentum|pubescence)\b",
    re.IGNORECASE,
)
SYMMETRY_NON_TARGET_SUBJECT = re.compile(
    rf"(?:{NON_TARGET_SUBJECT.pattern}|\b(?:lobes?|segments?)\b)",
    re.IGNORECASE,
)
MEASUREMENT = re.compile(
    r"(?P<raw>(?:c\.?|ca\.?)?\s*\d+(?:\.\d+)?(?:\s*[\-–—]\s*\d+(?:\.\d+)?)?)"
    r"(?:\s*\(\s*[+\-–—]?\s*\d+(?:\.\d+)?\s*\))?"
    r"\s*(?P<unit>mm|cm|m)\b",
    re.IGNORECASE,
)
CULTIVAR = re.compile(
    r"\b(?:cultivars?|cultivat(?:ed|ion)|hybrids?|grex|horticultural selection)\b|[×✕]",
    re.IGNORECASE,
)


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _normalise_text(value: Any) -> str:
    return unicodedata.normalize("NFKC", _text(value)).replace("—", "-").replace("–", "-")


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _json_gzip_write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        json.dump(value, handle, ensure_ascii=False, sort_keys=True)


def _json_gzip_read(path: Path) -> Any:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def _request_json(
    session: requests.Session,
    url: str,
    *,
    params: dict[str, Any] | None = None,
    attempts: int = 5,
) -> Any:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            response = session.get(url, params=params, timeout=180)
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(f"failed to fetch {url}") from last_error


def acquire_index(cache_dir: Path) -> list[dict[str, Any]]:
    """Return the frozen complete species index, using one public API request."""

    path = cache_dir / "flora_of_australia_species_index.json.gz"
    if path.exists():
        rows = _json_gzip_read(path)
    else:
        session = requests.Session()
        session.headers.update({"User-Agent": USER_AGENT})
        rows = _request_json(
            session,
            INDEX_URL,
            params={
                "opusId": OPUS_ID,
                "taxon": "species",
                "scientificName": "*",
                "max": 25_000,
                "offset": 0,
                "countChildren": "false",
                "sortBy": "name",
                "immediateChildrenOnly": "false",
            },
        )
        if not isinstance(rows, list) or len(rows) < 19_000:
            raise RuntimeError(f"unexpected Flora of Australia species index size: {len(rows)}")
        _json_gzip_write(path, rows)
    profile_ids = [_text(row.get("profileId")) for row in rows]
    names = [_text(row.get("scientificName")) for row in rows]
    if any(not value for value in profile_ids) or len(profile_ids) != len(set(profile_ids)):
        raise RuntimeError("Flora of Australia index has missing or duplicate profile IDs")
    if any(not value for value in names) or len(names) != len(set(names)):
        raise RuntimeError("Flora of Australia index has missing or duplicate species names")
    return rows


def _fetch_profile(profile_id: str, cache_dir: Path) -> tuple[str, str]:
    path = cache_dir / "profiles" / f"{profile_id}.json.gz"
    if path.exists():
        return profile_id, "cached"
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    payload = _request_json(
        session,
        PROFILE_URL.format(opus_id=OPUS_ID, profile_id=profile_id),
        # The compact response still contains the complete taxonomic chain but
        # omits child-profile expansion.  On the public service this is roughly
        # an order of magnitude faster for a source-scale acquisition.
        params={"fullClassification": "false"},
    )
    profile = payload.get("profile", {}) if isinstance(payload, dict) else {}
    if _text(profile.get("uuid")) != profile_id:
        raise RuntimeError(f"profile identity mismatch for {profile_id}")
    _json_gzip_write(path, payload)
    return profile_id, "fetched"


def acquire_profiles(
    selected_index: pd.DataFrame,
    cache_dir: Path,
    *,
    workers: int,
) -> dict[str, int]:
    """Fetch only exact master matches; completed profile files are never repeated."""

    ids = sorted(set(selected_index["profileId"].map(_text)))
    counts = {"requested": len(ids), "fetched": 0, "cached": 0, "failed": 0}
    failures: list[tuple[str, str]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(_fetch_profile, profile_id, cache_dir): profile_id
            for profile_id in ids
        }
        for completed, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            profile_id = futures[future]
            try:
                _, status = future.result()
                counts[status] += 1
            except (OSError, RuntimeError, ValueError) as exc:  # pragma: no cover
                counts["failed"] += 1
                failures.append((profile_id, repr(exc)))
            if completed % 250 == 0 or completed == len(futures):
                print(
                    f"profiles {completed}/{len(futures)} "
                    f"fetched={counts['fetched']} cached={counts['cached']} "
                    f"failed={counts['failed']}",
                    flush=True,
                )
    if failures:
        raise RuntimeError(f"{len(failures)} profile fetches failed; first={failures[0]}")
    return counts


def _attribute(profile: dict[str, Any], title: str) -> str:
    for item in profile.get("attributes") or []:
        if _text(item.get("title")).casefold() == title.casefold():
            return _normalise_text(item.get("plainText"))
    return ""


def _family(profile: dict[str, Any]) -> str:
    for item in profile.get("classification") or []:
        if _text(item.get("rank")).casefold() == "family":
            return _text(item.get("name"))
    return ""


def _sentences(description: str) -> list[str]:
    text = _normalise_text(description)
    # Preserve the exact normalized source span.  Semicolons remain inside a
    # sentence because botanical subjects commonly carry across them.  Requiring
    # an uppercase next token avoids cutting the ubiquitous ``c. 5 mm`` and
    # ``ca. 2 cm`` measurement abbreviations away from their floral subject.
    return [
        value.strip()
        for value in re.split(
            r"(?<!\bc\.)(?<!\bca\.)(?<=[.!?])\s+(?=[A-Z0-9])",
            text,
        )
        if value.strip()
    ]


def _nearest_subject(
    sentence: str,
    position: int,
    target_pattern: re.Pattern[str] = TARGET_SUBJECT,
    non_target_pattern: re.Pattern[str] = NON_TARGET_SUBJECT,
) -> str:
    # A floral noun can be the object of a comparison rather than the owner of
    # the following description: ``bracts overtopping flowers, pinkish green``.
    # Treating that occurrence as a floral subject would transfer bract colour
    # to the flower.
    def target_is_object(match: re.Match[str]) -> bool:
        prefix = sentence[max(0, match.start() - 40) : match.start()]
        suffix = sentence[match.end() : match.end() + 12]
        return bool(
            re.search(
                r"\b(?:overtopping|surrounding|subtending|exceeding|than)\s+$",
                prefix,
                re.IGNORECASE,
            )
            or re.search(r"\bwithout\b[^.;]{0,35}$", prefix, re.IGNORECASE)
            or re.match(r"\s+heads?\b", suffix, re.IGNORECASE)
        )

    subjects: list[tuple[int, str]] = []
    target_matches = [
        match
        for match in target_pattern.finditer(sentence)
        if match.start() <= position and not target_is_object(match)
    ]
    for match in target_matches:
        if match.start() <= position:
            subjects.append((match.end(), "target"))
    for match in non_target_pattern.finditer(sentence):
        if match.start() <= position:
            subjects.append((match.end(), "non_target"))

    # In noun phrases such as ``white hairs`` or ``yellow anthers``, the noun
    # immediately following the state owns it even when an earlier sentence
    # subject was ``flowers``.  This also positively handles ``white petals``.
    following: list[tuple[int, str]] = []
    for match in target_pattern.finditer(sentence, position):
        if not target_is_object(match):
            following.append((match.start(), "target"))
    for match in non_target_pattern.finditer(sentence, position):
        following.append((match.start(), "non_target"))
    if following:
        following_start, following_type = min(following)
        intervening = sentence[position:following_start]
        if following_start - position <= 24 and not re.search(r"[,;:.]", intervening):
            return following_type

    if subjects and position - max(subjects)[0] <= 140:
        nearest = max(subjects)
        if nearest[1] == "target":
            return "target"
        # Botanical descriptions frequently keep the leading floral subject
        # through a comparison: ``Petals shorter than calyx, yellow`` and
        # ``Corolla tube longer than sepals, urceolate``.  The comparison object
        # is not the owner of the post-comma colour or form.
        if target_matches:
            target = target_matches[-1]
            bridge = sentence[target.end() : position]
            if len(bridge) <= 170 and re.search(
                r"\b(?:shorter|longer|equal|as\s+long\s+as|exceed(?:s|ing)?)\b"
                r"[^,;]{0,45}\b(?:calyx|calyces|sepals?|bracts?)\b\s*,",
                bridge,
                re.IGNORECASE,
            ):
                return "target"
        return "non_target"
    if following and min(following)[0] - position <= 35:
        return min(following)[1]
    return ""


def _extract_pattern_states(
    description: str,
    patterns: Iterable[tuple[str, str]],
    *,
    require_floral_subject: bool,
    target_pattern: re.Pattern[str] = TARGET_SUBJECT,
    non_target_pattern: re.Pattern[str] = NON_TARGET_SUBJECT,
) -> tuple[set[str], list[str]]:
    states: set[str] = set()
    quotes: list[str] = []
    for sentence in _sentences(description):
        sentence_states: set[str] = set()
        for pattern, state in patterns:
            for match in re.finditer(pattern, sentence, re.IGNORECASE):
                prefix = sentence[max(0, match.start() - 35) : match.start()]
                if re.search(
                    r"(?:\bnot\b|\bnever\b|\bwithout\b|\bnon[- ])"
                    r"(?:\s+\w+){0,2}\s*$",
                    prefix,
                    re.IGNORECASE,
                ):
                    continue
                if (
                    require_floral_subject
                    and _nearest_subject(
                        sentence,
                        match.start(),
                        target_pattern,
                        non_target_pattern,
                    )
                    != "target"
                ):
                    continue
                sentence_states.add(state)
        if sentence_states:
            states.update(sentence_states)
            quotes.append(sentence)
    return states, quotes


def _extract_symmetry(description: str) -> tuple[set[str], list[str]]:
    patterns = (
        (r"\bactinomorphic\b", "actinomorphic"),
        (r"\bzygomorphic\b", "zygomorphic"),
        (r"\basymmetric\b", "asymmetric"),
    )
    return _extract_pattern_states(
        description,
        patterns,
        require_floral_subject=True,
        target_pattern=WHOLE_FLOWER_SUBJECT,
        non_target_pattern=SYMMETRY_NON_TARGET_SUBJECT,
    )


def _extract_form(description: str) -> tuple[set[str], list[str]]:
    whole_patterns = tuple(item for item in FORM_PATTERNS if item[1] != "spurred")
    whole_states, whole_quotes = _extract_pattern_states(
        description,
        whole_patterns,
        require_floral_subject=True,
        target_pattern=WHOLE_FLOWER_SUBJECT,
    )
    spur_states, spur_quotes = _extract_pattern_states(
        description,
        tuple(item for item in FORM_PATTERNS if item[1] == "spurred"),
        require_floral_subject=True,
        target_pattern=TARGET_SUBJECT,
    )
    return whole_states | spur_states, list(dict.fromkeys(whole_quotes + spur_quotes))


def _extract_inflorescence(description: str) -> tuple[set[str], list[str]]:
    states: set[str] = set()
    quotes: list[str] = []
    for sentence in _sentences(description):
        if re.search(r"\binfructescen", sentence, re.IGNORECASE):
            continue
        has_display_subject = bool(
            re.search(
                r"\b(?:inflorescences?|flowers?|flower[- ]heads?|capitula|capitulum|"
                r"racemes?|spikes?|panicles?|thyrses?|catkins?|umbels?|corymbs?|"
                r"cymes?|cincinni)\b",
                sentence,
                re.IGNORECASE,
            )
        )
        if not has_display_subject:
            continue
        sentence_states = {
            state
            for pattern, state in INFLORESCENCE_PATTERNS
            if re.search(pattern, sentence, re.IGNORECASE)
        }
        # A bare "solitary" must qualify a flower, not a leaf or plant.
        if "solitary" in sentence_states and not re.search(
            r"\b(?:flowers?\s+(?:usually\s+)?solitary|solitary\s+flowers?)\b",
            sentence,
            re.IGNORECASE,
        ):
            sentence_states.remove("solitary")
        if "solitary" in sentence_states and re.search(
            r"\bflowers?\s+solitary\s+within\s+(?:each|a)\s+bract\b|"
            r"\binflorescences?\b.*\bflowered\b.*\bflowers?\s+solitary\b",
            sentence,
            re.IGNORECASE,
        ):
            sentence_states.remove("solitary")
        if sentence_states:
            states.update(sentence_states)
            quotes.append(sentence)
    return states, quotes


@dataclass(frozen=True)
class MeasurementResult:
    states: frozenset[str]
    quotes: tuple[str, ...]
    raw_measurements: tuple[str, ...]


def _extract_measurements(
    description: str,
    trait: str,
    rules: dict[str, Any],
) -> MeasurementResult:
    states: set[str] = set()
    quotes: list[str] = []
    raw_values: list[str] = []
    for sentence in _sentences(description):
        normalized = sentence.replace("–", "-").replace("—", "-")
        if trait == "tube_depth_class":
            subject_pattern = re.compile(
                r"\b(?:(?:corolla|perianth|floral)[- ]+tube|"
                r"tube\s+of\s+the\s+corolla|"
                r"(?:corolla|perianth)\b[^.;]{0,80};\s*tube)\b",
                re.IGNORECASE,
            )
        else:
            subject_pattern = WHOLE_FLOWER_SUBJECT
        for subject in subject_pattern.finditer(normalized):
            prefix = normalized[max(0, subject.start() - 90) : subject.start()]
            if trait == "flower_size_class" and (
                (subject.start() > 0 and normalized[subject.start() - 1] == "-")
                or re.search(r"\b(?:in|during)\s+$", prefix, re.IGNORECASE)
                or re.search(
                    r"\b(?:with|bearing)\s+(?:\d+(?:\s*-\s*\d+)?\s+)?$",
                    prefix,
                    re.IGNORECASE,
                )
                or re.search(
                    r"\b(?:lip|limb|lobe|segment|part|portion)\s+of\s+(?:the\s+)?$",
                    prefix,
                    re.IGNORECASE,
                )
                or re.search(
                    r"\b(?:calyx|tube|limb|lobe|segment|bract|pedicel|cluster|"
                    r"bracteole|bud|cyme|inflorescence|panicle|raceme|spike|"
                    r"head|axis|fruit)s?\b[^.;:]{0,70}$",
                    prefix,
                    re.IGNORECASE,
                )
                or re.search(
                    r"\b(?:fruit|bract|sepal|calyx)\b[^.;:]{0,55}"
                    r"\b(?:exceeding|surrounding|subtending)\s+(?:the\s+)?$",
                    prefix,
                    re.IGNORECASE,
                )
                or re.search(r"\b(?:to|per)\s+(?:a\s+)?$", prefix, re.IGNORECASE)
            ):
                continue
            if re.search(r"\bwithout\b[^.;]{0,35}$", prefix, re.IGNORECASE):
                continue
            if trait == "flower_size_class" and re.search(
                r"\b(?:pedicels?|peduncles?|stalks?)\s+(?:of\s+)?(?:open\s+)?$",
                prefix,
                re.IGNORECASE,
            ):
                continue
            if re.search(
                r"\b(?:pedicels?|sepals?|lobes?|filaments?|bracts?|bracteoles?|segments?)\b"
                r"[^.;]{0,55}\b(?:shorter|longer|equal|emerging)\b"
                r"[^.;]{0,25}(?:than|to|from)\s+(?:the\s+)?$|"
                r"\bpedicels?\s*\+\s*$",
                prefix,
                re.IGNORECASE,
            ):
                continue
            if re.search(
                r"\b(?:scales?|bracts?|bracteoles?|lobes?|segments?)\b"
                r"[^.;]{0,70}\b(?:equal\s+to\s+or\s+)?exceed(?:s|ing)?\s+"
                r"(?:the\s+)?$",
                prefix,
                re.IGNORECASE,
            ):
                continue
            if trait == "flower_size_class" and re.match(
                r"(?:[- ]+tube|\s+(?:stage|heads?|limbs?|lobes?|segments?|parts?|"
                r"remains?|claws?))\b",
                normalized[subject.end() : subject.end() + 12],
                re.IGNORECASE,
            ):
                continue
            window = normalized[subject.end() : subject.end() + 110]
            measurement = MEASUREMENT.search(window)
            if not measurement:
                continue
            bridge = window[: measurement.start()]
            bridge_without_comparison = re.sub(
                r"\b(?:shorter|longer)\s+than\s+(?:the\s+)?"
                r"(?:petals?|tepals?|sepals?|calyx|bracts?)\b",
                "",
                bridge,
                flags=re.IGNORECASE,
            )
            if re.search(
                r"\b(?:petals?|tepals?|sepals?|calyx|bracts?|bracteoles?|pedicels?|"
                r"peduncles?|lips?|lobes?|segments?|standards?|wings?|keels?|"
                r"vexilla?|anthers?|stamens?|staminodes?|gynophores?|styles?|ovary|"
                r"ovaries|androecium|gynoecium|disc|scales?|tubes?|hypanthium|"
                r"operculum|filaments?|spathes?|limbs?|glumes?|lemmas?|paleas?|awns?|"
                r"bristles?|buds?|"
                r"fruit|fruiting|capsules?|leaves|leaf|inflorescences?|"
                r"panicles?|pseudoracemes?|racemes?|spikes?|clusters?|fascicles?|"
                r"flower[- ]heads?|heads?|fertile\s+portion|"
                r"stalks?|scapes?|rachis|axes?|axis)\b",
                bridge_without_comparison,
                re.IGNORECASE,
            ):
                continue
            if re.search(r"\bflowers?\s+not\s+seen\b", normalized, re.IGNORECASE):
                continue
            after = window[measurement.end() : measurement.end() + 25]
            allowed_dimensions = (
                r"\s*(?:long|deep)\b"
                if trait == "tube_depth_class"
                else r"\s*(?:long|wide|diam(?:eter)?|across|high|deep)\b"
            )
            if not re.match(allowed_dimensions, after, re.IGNORECASE):
                continue
            if trait == "flower_size_class" and re.match(
                r"\s*(?:long|wide|diam(?:eter)?|across|high|deep)\s*,\s*"
                r"(?:\w+[ -]+){0,2}(?:the\s+)?"
                r"(?:pedicels?|peduncles?|stalks?|bracts?|fruits?)\b",
                after,
                re.IGNORECASE,
            ):
                continue
            if trait == "flower_size_class" and re.match(
                r"\s*(?:long|wide|diam(?:eter)?|across|high|deep)\s+"
                r"(?:\w+[ -]+){0,2}(?:on\s+)?(?:the\s+)?"
                r"(?:pedicels?|peduncles?|stalks?|axes?|axis)\b",
                after,
                re.IGNORECASE,
            ):
                continue
            raw = measurement.group("raw").strip().replace("-", "-")
            unit = measurement.group("unit").casefold()
            state = derive_class(trait, raw, unit, rules)
            if state:
                states.add(state)
                quotes.append(sentence)
                raw_values.append(f"{raw} {unit}")
    return MeasurementResult(
        states=frozenset(states),
        quotes=tuple(dict.fromkeys(quotes)),
        raw_measurements=tuple(dict.fromkeys(raw_values)),
    )


def extract_description(description: str, rules: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Return normalized values and exact description spans for all six traits."""

    extracted: dict[str, dict[str, Any]] = {}
    color_states, color_quotes = _extract_pattern_states(
        description,
        COLOUR_PATTERNS,
        require_floral_subject=True,
        target_pattern=TARGET_SUBJECT,
    )
    form_states, form_quotes = _extract_form(description)
    symmetry_states, symmetry_quotes = _extract_symmetry(description)
    display_states, display_quotes = _extract_inflorescence(description)
    values = (
        ("flower_primary_color", color_states, color_quotes, ()),
        ("floral_form", form_states, form_quotes, ()),
        ("floral_symmetry", symmetry_states, symmetry_quotes, ()),
        ("inflorescence_display", display_states, display_quotes, ()),
    )
    for trait, states, quotes, measurements in values:
        if states:
            extracted[trait] = {
                "value": "|".join(sorted(states)),
                "quote": " ".join(dict.fromkeys(quotes)),
                "raw_measurements": measurements,
            }
    for trait in ("flower_size_class", "tube_depth_class"):
        result = _extract_measurements(description, trait, rules)
        if result.states:
            extracted[trait] = {
                "value": "|".join(sorted(result.states)),
                "quote": " ".join(result.quotes),
                "raw_measurements": result.raw_measurements,
            }
    return extracted


def _source_lineage(source: str, profile_id: str) -> tuple[str, str]:
    normalized = re.sub(r"[^a-z0-9]+", " ", source.casefold()).strip()
    if normalized:
        return f"citation:{_sha(normalized)[:32]}", "normalized_source_attribute_fingerprint"
    return f"flora-of-australia-profile:{profile_id}", "official_profile_uuid_fallback"


def _evidence_row(
    *,
    species: str,
    profile_id: str,
    trait: str,
    value: str,
    quote: str,
    source: str,
) -> dict[str, str]:
    lineage, lineage_method = _source_lineage(source, profile_id)
    url = f"https://profiles.ala.org.au/opus/{OPUS_SHORT_NAME}/profile/{profile_id}"
    record_key = _sha(f"{profile_id}|{trait}|{value}|{quote}")[:24]
    return {
        "accepted_species": species,
        "axis": AXIS[trait],
        "trait_name": trait,
        "normalized_value": value,
        "quality": "high",
        "source_group": "latest_public_web",
        "source_provider": "flora_of_australia_official",
        "source_url": url,
        "source_record_id": f"foa-{record_key}",
        "source_citation": source or "Flora of Australia official species profile",
        "source_excerpt": quote,
        "evidence_scope": "species_direct",
        "name_match_method": "accepted_name_exact_plus_family_exact",
        "source_lineage": lineage,
        "lineage_method": lineage_method,
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "source_file": SOURCE_FILE,
        "acceptance_contract": ACCEPTANCE_CONTRACT,
    }


def build_evidence(
    selected_index: pd.DataFrame,
    master: pd.DataFrame,
    cache_dir: Path,
    baseline_pairs: set[tuple[str, str]],
    rules: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    master_family = dict(zip(master["accepted_species"], master["family"], strict=True))
    rows: list[dict[str, str]] = []
    treatments: list[dict[str, Any]] = []
    for item in selected_index.sort_values("scientificName").itertuples(index=False):
        payload = _json_gzip_read(cache_dir / "profiles" / f"{item.profileId}.json.gz")
        profile = payload.get("profile", {})
        species = _text(profile.get("scientificName"))
        family = _family(profile)
        expected_family = _text(master_family.get(species))
        rank = _text(profile.get("rank")).casefold()
        description = _attribute(profile, "Description")
        source = _attribute(profile, "Source")
        identity_ok = bool(
            species == _text(item.scientificName)
            and species in master_family
            and rank == "species"
            and family
            and family == expected_family
        )
        cultivar_excluded = bool(CULTIVAR.search(f"{species} {description}"))
        extracted = (
            extract_description(description, rules)
            if identity_ok and description and not cultivar_excluded
            else {}
        )
        novel_traits = {
            trait: value
            for trait, value in extracted.items()
            if (species, trait) not in baseline_pairs
        }
        for trait, value in novel_traits.items():
            rows.append(
                _evidence_row(
                    species=species,
                    profile_id=_text(item.profileId),
                    trait=trait,
                    value=_text(value["value"]),
                    quote=_text(value["quote"]),
                    source=source,
                )
            )
        treatments.append(
            {
                "accepted_species": _text(item.scientificName),
                "profile_id": _text(item.profileId),
                "profile_scientific_name": species,
                "profile_rank": rank,
                "profile_family": family,
                "master_family": expected_family,
                "identity_gate_passed": identity_ok,
                "cultivar_excluded": cultivar_excluded,
                "description_chars": len(description),
                "all_extracted_traits": "|".join(sorted(extracted)),
                "novel_extracted_traits": "|".join(sorted(novel_traits)),
                "source_attribute": source,
                "source_url": (
                    f"https://profiles.ala.org.au/opus/{OPUS_SHORT_NAME}/profile/"
                    f"{_text(item.profileId)}"
                ),
            }
        )
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
    if not evidence.empty:
        evidence = evidence.drop_duplicates("source_record_id").sort_values(
            ["trait_name", "accepted_species", "source_record_id"], kind="stable"
        )
    return evidence.reset_index(drop=True), pd.DataFrame(treatments)


def deterministic_audit_queue(evidence: pd.DataFrame, size: int = AUDIT_SIZE) -> pd.DataFrame:
    if len(evidence) < size:
        raise ValueError(f"only {len(evidence)} candidates available for {size}-row audit")
    pool = evidence.copy()
    pool["_audit_hash"] = pool["source_record_id"].map(
        lambda value: _sha(f"{AUDIT_SEED}|{value}")
    )
    selected: list[pd.DataFrame] = []
    # Six traits x 20 rows guarantees every scope can be independently approved.
    for _, group in pool.groupby("trait_name", sort=True):
        selected.append(group.sort_values("_audit_hash").head(min(20, len(group))))
    audit = pd.concat(selected, ignore_index=False).drop_duplicates("source_record_id")
    remaining = pool.loc[~pool.index.isin(audit.index)].sort_values("_audit_hash")
    if len(audit) > size:
        audit = audit.sort_values("_audit_hash").head(size)
    elif len(audit) < size:
        audit = pd.concat([audit, remaining.head(size - len(audit))])
    audit = audit.sort_values(["trait_name", "_audit_hash"], kind="stable").rename(
        columns={
            "source_record_id": "candidate_id",
            "source_provider": "provider",
            "source_excerpt": "exact_supporting_quote",
        }
    )
    audit["accepted_correct"] = ""
    audit["cultivar_status"] = "wild_or_naturalised_flora_treatment"
    audit["reviewer"] = ""
    audit["reviewed_at_utc"] = ""
    audit["audit_reason"] = ""
    return audit[
        [
            "candidate_id",
            "accepted_species",
            "trait_name",
            "normalized_value",
            "provider",
            "source_url",
            "source_lineage",
            "source_citation",
            "exact_supporting_quote",
            "name_match_method",
            "cultivar_status",
            "accepted_correct",
            "reviewer",
            "reviewed_at_utc",
            "audit_reason",
        ]
    ].reset_index(drop=True)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(
    *,
    master_csv: Path,
    universe_coverage_csv: Path,
    baseline_evidence_csv: Path,
    cache_dir: Path,
    output_dir: Path,
    measurement_config: Path,
    workers: int,
) -> dict[str, Any]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    universe = pd.read_csv(
        universe_coverage_csv,
        dtype=str,
        usecols=["accepted_species"],
    ).fillna("")
    universe_names = set(universe["accepted_species"])
    if len(universe_names) != 106_295:
        raise ValueError("coverage universe must contain exactly 106,295 accepted species")
    master = master.loc[master["accepted_species"].isin(universe_names)].copy()
    if len(master) != 106_295 or master["accepted_species"].nunique() != 106_295:
        raise ValueError("master does not uniquely cover the 106,295-species universe")
    baseline = pd.read_csv(baseline_evidence_csv, dtype=str).fillna("")
    baseline_pairs = set(zip(baseline["accepted_species"], baseline["trait_name"], strict=True))
    index_rows = acquire_index(cache_dir)
    index = pd.DataFrame(index_rows)
    index["scientificName"] = index["scientificName"].map(_text)
    selected = index.loc[index["scientificName"].isin(set(master["accepted_species"]))].copy()
    selected = selected.drop_duplicates(["scientificName", "profileId"])
    fetch_counts = acquire_profiles(selected, cache_dir, workers=workers)
    rules = load_rules(measurement_config)
    evidence, treatments = build_evidence(
        selected,
        master,
        cache_dir,
        baseline_pairs,
        rules,
    )
    audit = deterministic_audit_queue(evidence)
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / SOURCE_FILE
    treatment_path = output_dir / "flora_of_australia_treatment_inventory.csv.gz"
    audit_path = output_dir / "flora_of_australia_audit_queue_200.csv"
    index_path = output_dir / "flora_of_australia_species_index.csv.gz"
    evidence.to_csv(evidence_path, index=False, compression="gzip")
    treatments.to_csv(treatment_path, index=False, compression="gzip")
    audit.to_csv(audit_path, index=False)
    index.to_csv(index_path, index=False, compression="gzip")
    summary = {
        "contract": "flora_of_australia_source_scale_acquisition_v1",
        "source_run_id": SOURCE_RUN_ID,
        "source_landing_url": SOURCE_LANDING_URL,
        "collection_species_profiles": len(index),
        "exact_master_profiles": len(selected),
        "identity_gate_passed": int(treatments["identity_gate_passed"].sum()),
        "profiles_with_description": int(treatments["description_chars"].gt(0).sum()),
        "profile_fetch": fetch_counts,
        "search_queries": 0,
        "collection_index_requests": 0 if fetch_counts["fetched"] == 0 else 1,
        "profile_requests_this_run": fetch_counts["fetched"],
        "api_cost_usd": 0.0,
        "baseline_species_trait_pairs": len(baseline_pairs),
        "novel_candidate_rows": len(evidence),
        "novel_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "novel_species_axis": int(
            evidence[["accepted_species", "axis"]].drop_duplicates().shape[0]
        ),
        "unique_species_with_novel_evidence": int(evidence["accepted_species"].nunique()),
        "by_trait": (
            evidence.groupby("trait_name")
            .agg(rows=("accepted_species", "size"), species=("accepted_species", "nunique"))
            .reset_index()
            .to_dict("records")
        ),
        "audit_queue_rows": len(audit),
        "audit_seed": AUDIT_SEED,
        "measurement_classification_rule_version": rules.get("classification_rule_version"),
        "guardrails": {
            "species_identity": "accepted_name_exact_plus_family_exact",
            "rank": "species_only",
            "description_field": "Description_only",
            "cultivar_hybrid_excluded": True,
            "snippet_as_evidence": False,
            "family_inference": False,
            "global_fallback": False,
            "axis_level_extraction": False,
        },
    }
    summary_path = output_dir / "flora_of_australia_source_scale_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    manifest = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "file_manifest.sha256":
            manifest.append(f"{_sha256_file(path)}  {path.name}")
    (output_dir / "file_manifest.sha256").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--master-csv",
        type=Path,
        default=Path("data/v2/staging/gbif/collected/island_taxa.csv"),
    )
    parser.add_argument("--baseline-evidence-csv", type=Path, required=True)
    parser.add_argument("--universe-coverage-csv", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--measurement-config",
        type=Path,
        default=Path("config/measurement_classification.yml"),
    )
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    if not 1 <= args.workers <= 32:
        raise ValueError("workers must be between 1 and 32")
    print(
        json.dumps(
            run(
                master_csv=args.master_csv,
                universe_coverage_csv=args.universe_coverage_csv,
                baseline_evidence_csv=args.baseline_evidence_csv,
                cache_dir=args.cache_dir,
                output_dir=args.output_dir,
                measurement_config=args.measurement_config,
                workers=args.workers,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
