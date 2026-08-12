"""Build reviewed Medium flower-colour evidence from IOSPE treatments.

The Internet Orchid Species Photo Encyclopedia (IOSPE) is treated as a
professional specialist Web flora, not as an image classifier.  Only text from
the species-treatment description is eligible.  Search/index text, common
names, photographs, cultivar captions, synonym redirects, genus inference and
non-floral colours are excluded.

The parser is deliberately source-specific only at the page-structure seam.
Its evidence contract is the repository-wide species-direct ledger: an exact
accepted-species heading, an exact supporting quote, immutable URL and page
hash, source lineage, and a complete multistate colour set for the quoted
eligible floral organs.
"""

from __future__ import annotations

import hashlib
import json
import re
import unicodedata
from pathlib import Path
from typing import Annotated, Any
from urllib.parse import urlsplit, urlunsplit

import pandas as pd
import typer
from bs4 import BeautifulSoup

from island_v2.integrated_trait_coverage import TRAIT_TO_AXIS

app = typer.Typer(add_completion=False, no_args_is_help=True)

TRAIT = "flower_primary_color"
RUN_ID = "iospe-orchid-colour-20260812"
ARTIFACT = "iospe-orchid-colour-reviewed-source-package-20260812"
RETRIEVAL_DATE = "2026-08-12"
AUDIT_SEED = "iospe-orchid-colour-source-holdout-audit-20260812-v2"
AUDIT_SIZE = 200

COLOUR_TERMS: dict[str, tuple[str, ...]] = {
    "white": (
        "white",
        "whitish",
        "whiter",
        "withe",
        "cream",
        "creamy",
        "creamish",
        "ivory",
    ),
    "yellow_orange": (
        "yellow",
        "yellowish",
        "orange",
        "apricot",
        "golden",
        "gold",
        "ochre",
        "ochreous",
        "ocher",
        "ocherous",
        "terracotta",
        "sulphur",
        "sulfur",
    ),
    "red_pink": (
        "red",
        "reddish",
        "pink",
        "pinkish",
        "scarlet",
        "crimson",
        "rose",
        "reose",
        "rosy",
        "magenta",
        "carmine",
        "salmon",
        "coral",
        "wine-coloured",
        "wine-colored",
    ),
    "blue_purple": (
        "blue",
        "bluish",
        "purple",
        "purplish",
        "purplsh",
        "pirple",
        "violet",
        "lavender",
        "lilac",
        "mauve",
    ),
    "green_brown_inconspicuous": (
        "green",
        "greenish",
        "brown",
        "brownish",
        "maroon",
        "black",
        "blackish",
        "buff",
        "fawn",
        "chocolate",
        "bronze",
        "bronzy",
        "gereenish",
    ),
}
TERM_TO_STATE = {
    term: state for state, terms in COLOUR_TERMS.items() for term in terms
}
COLOUR_PATTERN = re.compile(
    "|".join(
        rf"(?<![A-Za-z]){re.escape(term)}(?![A-Za-z])"
        for term in sorted(TERM_TO_STATE, key=len, reverse=True)
    ),
    re.IGNORECASE,
)

# Sepals are part of the orchid display and are therefore eligible here.
# Small coloured markings on an eligible organ are retained so multistate
# composition is not flattened.  Explicitly non-target organs still win when
# they are closer to a colour word than an eligible floral organ.
ELIGIBLE_ORGAN = re.compile(
    r"\b(?:flowers?|florets?|corollas?|petals?|sepals?|perianths?|tepals?|"
    r"lips?|labellums?|labella|spurs?|pouches?)\b",
    re.IGNORECASE,
)
EXCLUDED_ORGAN = re.compile(
    r"\b(?:leaves?|leaflets?|foliage|fruits?|berries|seeds?|stems?|bark|"
    r"branches?|branchlets?|roots?|tubers?|pseudobulbs?|rhizomes?|rachis|rachillae?|"
    r"peduncles?|pedicels?|bracts?|involucres?|paleae?|caly(?:x|ces)|"
    r"hypanth(?:ium|ia)|ovaries|ovary|carpels?|anthers?|stamens?|filaments?|"
    r"styles?|stigmas?|columns?|column|glands?|"
    r"hairs?|trichomes?|indumentum|pubescence|spines?|prickles?|glumes?|"
    r"lemmas?|spikelets?|sheaths?|tomentum|inflorescences?|drupes?|buds?|"
    r"appendages?|feet|foot)\b",
    re.IGNORECASE,
)
COMPARATIVE_OR_UNCERTAIN = re.compile(
    r"\b(?:unlike|similar to|similar flowers|differs? from|close to|"
    r"compared to|compared with|related to|distinguished from|"
    r"separated morphologically|other species|can be confused|differs? in|"
    r"all (?:two|three|four|five|\d+) have|(?:which|that) (?:has|have)|"
    r"may be|maybe|perhaps|misidentified|incorrect|no way of verifying|"
    r"no way to verify|photo is suspect)\b",
    re.IGNORECASE,
)
CULTIVAR_OR_HYBRID = re.compile(
    r"\b(?:cultivar|garden hybrid|hybrid cultivar|clone|variety|cv\.)\b",
    re.IGNORECASE,
)
DESCRIPTION_START = re.compile(
    r"\bfound\s+(?:in|only\s+in|on)\b", re.IGNORECASE
)
SYNONYM_REDIRECT = re.compile(
    r"\bsee\s*-\s*[A-Z][a-z]+\s+[a-z][a-z-]+", re.IGNORECASE
)
WEAK_FRAGMENT = re.compile(
    r"^\s*(?:,|from\s+[A-Z][^.;]{0,100},\s+has\b)", re.IGNORECASE
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _id(*parts: object, length: int = 24) -> str:
    payload = "\n".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_url(value: object) -> str:
    parts = urlsplit(_text(value))
    return urlunsplit((parts.scheme.casefold(), parts.netloc.casefold(), parts.path, "", ""))


def _normalized_excerpt(value: object) -> str:
    value = unicodedata.normalize("NFKC", _text(value)).casefold()
    return " ".join(re.findall(r"[\w]+", value, flags=re.UNICODE))


def _distance(match: re.Match[str], position: int) -> int:
    if match.end() <= position:
        return position - match.end()
    return match.start() - position


def _eligible_windows(text: str) -> list[tuple[int, int]]:
    """Return narrow display-description windows around eligible organs."""

    windows: list[tuple[int, int]] = []
    for organ in ELIGIBLE_ORGAN.finditer(text):
        left_sentence = max(text.rfind(".", 0, organ.start()), text.rfind(";", 0, organ.start())) + 1
        right_candidates = [
            position
            for position in (text.find(".", organ.end()), text.find(";", organ.end()))
            if position >= 0
        ]
        right_sentence = min(right_candidates) if right_candidates else len(text)
        if organ.group(0).casefold().startswith("flower"):
            carrying = max(
                text.rfind("carrying", left_sentence, organ.start()),
                text.rfind("carries", left_sentence, organ.start()),
            )
            if carrying >= 0 and organ.start() - carrying <= 240:
                left = carrying
            else:
                left = max(left_sentence, organ.start() - 100)
            right = min(right_sentence, organ.end() + 320)
        else:
            left = max(left_sentence, organ.start() - 55)
            right = min(right_sentence, organ.end() + 140)
        windows.append((left, right))
    return windows


def observed_orchid_colour_states(value: object) -> set[str]:
    """Return colours locally attached to eligible orchid floral organs."""

    text = "".join(
        character
        for character in unicodedata.normalize("NFKD", _text(value))
        if not unicodedata.combining(character)
    )
    windows = _eligible_windows(text)
    states: set[str] = set()
    for colour in COLOUR_PATTERN.finditer(text):
        if not any(start <= colour.start() < end for start, end in windows):
            continue
        after = text[colour.end() : colour.end() + 70]
        before = text[max(0, colour.start() - 70) : colour.start()]
        # A colour immediately modifying a floral bract is not petaloid.  The
        # ``carrying <colour> flowers`` form remains eligible because the last
        # carrying cue occurs after the bract description.
        bract_after = re.match(
            r"^[^.;]{0,42}\bfloral bracts?\b", after, re.IGNORECASE
        )
        if bract_after and not re.search(r"\bcarrying\b", bract_after.group(0), re.IGNORECASE):
            continue
        if re.search(
            r"\b(?:not|without|does not have|do not have)\b[^.;]{0,35}$",
            before,
            re.IGNORECASE,
        ):
            continue
        # Colours immediately before ``flowers`` are predicate complements
        # even when ``inflorescence`` or ``bracts`` is closer in raw character
        # distance (e.g. ``carrying white flowers``).
        flower_after = re.search(r"^[^.;]{0,60}\bflowers?\b", after, re.IGNORECASE)
        excluded_between = (
            EXCLUDED_ORGAN.search(flower_after.group(0)) if flower_after else None
        )
        forced_flower_predicate = bool(flower_after and not excluded_between)
        local_start = max(0, colour.start() - 55)
        local_end = min(len(text), colour.end() + 55)
        local = text[local_start:local_end]
        local_colour = colour.start() - local_start
        eligible = list(ELIGIBLE_ORGAN.finditer(local))
        excluded = list(EXCLUDED_ORGAN.finditer(local))
        nearest_eligible = min(
            (_distance(item, local_colour) for item in eligible), default=10_000
        )
        nearest_excluded = min(
            (_distance(item, local_colour) for item in excluded), default=10_000
        )
        if nearest_eligible > 45 and not forced_flower_predicate:
            continue
        # A directly named non-target organ wins, except when an eligible
        # organ is also explicitly coordinated in the same short phrase
        # (for example, ``lip and column are white``).
        coordinated_display_and_excluded = re.search(
            r"\b(?:lips?|labellums?|labella|spurs?|pouches?|petals?|sepals?|"
            r"tepals?)\b[^,.;]{0,30}\b(?:and|or)\b[^,.;]{0,12}"
            r"\b(?:columns?|anthers?)\b[^,.;]{0,18}$",
            before,
            re.IGNORECASE,
        )
        if (
            nearest_excluded < nearest_eligible
            and nearest_excluded <= 22
            and not coordinated_display_and_excluded
            and not (
                forced_flower_predicate
                and re.search(r"\bcarrying\b[^.;]{0,24}$", before, re.IGNORECASE)
            )
        ):
            continue
        scent_context = text[max(0, colour.start() - 12) : colour.end() + 28]
        if re.search(
            r"\brose\s+to\s+(?:spicy|sweet|floral)\s+scented\b",
            scent_context,
            re.IGNORECASE,
        ):
            continue
        states.add(TERM_TO_STATE[colour.group(0).casefold()])
    return states


def _heading_is_exact(strings: list[str], species: str) -> tuple[bool, str]:
    pattern = re.compile(rf"^\s*{re.escape(species)}(?:\s|$)", re.IGNORECASE)
    for item in strings[:15]:
        if pattern.search(item):
            return True, item
        if species.casefold() in item.casefold():
            return False, item
    return False, ""


def _description_nodes(strings: list[str]) -> list[str]:
    start = next(
        (index for index, item in enumerate(strings) if DESCRIPTION_START.search(item)),
        None,
    )
    if start is None:
        return []
    nodes: list[str] = []
    for item in strings[start:]:
        if item.casefold() in {"synonyms", "references"}:
            break
        nodes.append(item)
    return nodes


def extract_page_candidates(
    *, species: str, source_url: str, page_sha256: str, html: bytes
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Extract fail-closed source-text candidates from one fetched page."""

    soup = BeautifulSoup(html, "html.parser")
    text_root = soup.body or soup
    strings = [_text(item) for item in text_root.stripped_strings if _text(item)]
    identity_exact, heading = _heading_is_exact(strings, species)
    redirect = any(SYNONYM_REDIRECT.search(item) for item in strings[:15])
    description = _description_nodes(strings)
    selected: list[dict[str, str]] = []
    if identity_exact and not redirect:
        for node_index, quote in enumerate(description):
            if (
                COMPARATIVE_OR_UNCERTAIN.search(quote)
                or CULTIVAR_OR_HYBRID.search(quote)
                or WEAK_FRAGMENT.search(quote)
            ):
                continue
            states = observed_orchid_colour_states(quote)
            if not states:
                continue
            fingerprint = hashlib.sha256(
                _normalized_excerpt(quote).encode("utf-8")
            ).hexdigest()
            lineage = f"url:{_canonical_url(source_url)}"
            value = "|".join(sorted(states))
            # Candidate identity is stable across ontology corrections to the
            # same immutable source quote; value is intentionally not hashed.
            record_id = "iospe-colour:" + _id(species, TRAIT, lineage, fingerprint)
            selected.append(
                {
                    "accepted_species": species,
                    "axis": TRAIT_TO_AXIS[TRAIT],
                    "trait_name": TRAIT,
                    "normalized_value": value,
                    "quality": "medium",
                    "source_group": "iospe_orchid_colour",
                    "source_provider": (
                        "Internet Orchid Species Photo Encyclopedia (IOSPE)"
                    ),
                    "source_url": _canonical_url(source_url),
                    "source_record_id": record_id,
                    "source_citation": (
                        "Pfahl, J. Internet Orchid Species Photo Encyclopedia, "
                        f"species treatment: {species}"
                    ),
                    "source_excerpt": quote,
                    "evidence_scope": "species_direct",
                    "name_match_method": "exact_accepted_species_page_heading",
                    "source_lineage": lineage,
                    "lineage_method": (
                        "canonical_species_treatment_url_plus_content_fingerprint"
                    ),
                    "source_run_id": RUN_ID,
                    "source_artifact": ARTIFACT,
                    "source_file": f"live-page-sha256:{page_sha256}",
                    "retrieval_date": RETRIEVAL_DATE,
                    "acceptance_contract": (
                        "iospe_exact_heading_description_text_orchid_flower_colour_v1"
                    ),
                    "content_fingerprint": fingerprint,
                    "page_sha256": page_sha256,
                    "page_title": (
                        _text(soup.title.get_text(" ", strip=True))
                        if soup.title
                        else species
                    ),
                    "description_node_index": str(node_index),
                }
            )
    audit = {
        "accepted_species": species,
        "source_url": _canonical_url(source_url),
        "identity_exact": identity_exact,
        "heading_text": heading,
        "synonym_redirect": redirect,
        "description_nodes": len(description),
        "candidate_nodes": len(selected),
        "selected": bool(selected),
    }
    return selected, audit


def deterministic_audit_template(
    evidence: pd.DataFrame,
    *,
    audit_size: int = AUDIT_SIZE,
    excluded_species_quotes: frozenset[tuple[str, str]] = frozenset(),
) -> pd.DataFrame:
    """Sample pages reproducibly, covering every represented genus first."""

    representatives = evidence.copy()
    if excluded_species_quotes:
        retained = [
            (_text(species), _text(quote)) not in excluded_species_quotes
            for species, quote in zip(
                representatives["accepted_species"],
                representatives["source_excerpt"],
                strict=True,
            )
        ]
        representatives = representatives.loc[retained].copy()
    representatives["genus"] = representatives["accepted_species"].str.split().str[0]
    representatives["n_states"] = representatives["normalized_value"].str.count(r"\|") + 1
    representatives["quote_len"] = representatives["source_excerpt"].str.len()
    representatives = (
        representatives.sort_values(
            ["accepted_species", "n_states", "quote_len", "source_record_id"],
            ascending=[True, False, True, True],
            kind="stable",
        )
        .drop_duplicates("accepted_species")
        .copy()
    )
    representatives["audit_hash"] = representatives["source_record_id"].map(
        lambda value: _id(AUDIT_SEED, value, length=64)
    )
    first_per_genus = (
        representatives.sort_values(["genus", "audit_hash"], kind="stable")
        .groupby("genus", as_index=False)
        .head(1)
    )
    if len(first_per_genus) > audit_size:
        first_per_genus = first_per_genus.sort_values("audit_hash").head(audit_size)
    remaining = representatives.loc[
        ~representatives["source_record_id"].isin(first_per_genus["source_record_id"])
    ].sort_values("audit_hash", kind="stable")
    sample = pd.concat(
        [first_per_genus, remaining.head(audit_size - len(first_per_genus))],
        ignore_index=True,
    ).sort_values("audit_hash", kind="stable")
    if len(sample) != audit_size:
        raise ValueError(
            f"IOSPE audit requires {audit_size} unique treatment pages, found {len(sample)}"
        )
    return pd.DataFrame(
        {
            "candidate_id": sample["source_record_id"],
            "accepted_species": sample["accepted_species"],
            "trait_name": sample["trait_name"],
            "normalized_value": sample["normalized_value"],
            "source_url": sample["source_url"],
            "page_title": sample["page_title"],
            "exact_supporting_quote": sample["source_excerpt"],
            "content_fingerprint": sample["content_fingerprint"],
            "accepted_correct": "",
            "cultivar_status": "",
            "reviewer": "",
            "reviewed_at_utc": "",
            "audit_reason": "",
            "source_lineage": sample["source_lineage"],
            "source_provider": sample["source_provider"],
            "source_citation": sample["source_citation"],
            "source_excerpt": sample["source_excerpt"],
            "audit_hash": sample["audit_hash"],
        }
    ).reset_index(drop=True)


def build_checkpoint(
    fetch_manifest: pd.DataFrame,
    development_audit: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required = {
        "accepted_species",
        "source_url",
        "fetch_status",
        "content_sha256",
        "cache_file",
    }
    if missing := required.difference(fetch_manifest.columns):
        raise ValueError(f"IOSPE fetch manifest missing columns: {sorted(missing)}")
    evidence_rows: list[dict[str, str]] = []
    selection_rows: list[dict[str, Any]] = []
    for record in fetch_manifest.fillna("").to_dict("records"):
        if _text(record["fetch_status"]) != "http:200":
            selection_rows.append(
                {
                    "accepted_species": _text(record["accepted_species"]),
                    "source_url": _canonical_url(record["source_url"]),
                    "identity_exact": False,
                    "heading_text": "",
                    "synonym_redirect": False,
                    "description_nodes": 0,
                    "candidate_nodes": 0,
                    "selected": False,
                    "selection_reason": "fetch_failed",
                }
            )
            continue
        path = Path(_text(record["cache_file"]))
        if not path.exists():
            raise ValueError(f"IOSPE cache file is missing: {path}")
        page_hash = _sha256(path)
        if page_hash != _text(record["content_sha256"]):
            raise ValueError(f"IOSPE page hash mismatch: {path}")
        rows, audit = extract_page_candidates(
            species=_text(record["accepted_species"]),
            source_url=_text(record["source_url"]),
            page_sha256=page_hash,
            html=path.read_bytes(),
        )
        evidence_rows.extend(rows)
        audit["selection_reason"] = (
            "selected_exact_description_colour"
            if rows
            else "no_eligible_exact_description_colour"
        )
        selection_rows.append(audit)
    evidence = pd.DataFrame(evidence_rows)
    if evidence.empty:
        raise ValueError("IOSPE checkpoint produced no evidence")
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("IOSPE checkpoint produced duplicate source record IDs")
    evidence = evidence.sort_values(
        ["accepted_species", "source_record_id"], kind="stable"
    ).reset_index(drop=True)
    selection = pd.DataFrame(selection_rows).sort_values(
        ["selected", "accepted_species"], ascending=[False, True], kind="stable"
    )
    excluded_species_quotes: frozenset[tuple[str, str]] = frozenset()
    if development_audit is not None:
        required_development = {"accepted_species", "exact_supporting_quote"}
        if missing := required_development.difference(development_audit.columns):
            raise ValueError(
                f"development audit missing columns: {sorted(missing)}"
            )
        excluded_species_quotes = frozenset(
            zip(
                development_audit["accepted_species"].map(_text),
                development_audit["exact_supporting_quote"].map(_text),
                strict=True,
            )
        )
    audit = deterministic_audit_template(
        evidence, excluded_species_quotes=excluded_species_quotes
    )
    return evidence, audit, selection.reset_index(drop=True)


@app.command("build")
def build(
    fetch_manifest_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    development_audit_csv: Annotated[
        Path | None, typer.Option(exists=True, dir_okay=False)
    ] = None,
) -> None:
    manifest = pd.read_csv(fetch_manifest_csv, dtype=str).fillna("")
    development_audit = (
        pd.read_csv(development_audit_csv, dtype=str).fillna("")
        if development_audit_csv
        else None
    )
    evidence, audit, selection = build_checkpoint(manifest, development_audit)
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "iospe_orchid_colour_evidence_20260812.csv.gz"
    audit_path = output_dir / "iospe_orchid_colour_audit_200_20260812.csv"
    selection_path = output_dir / "iospe_orchid_colour_selection_20260812.csv.gz"
    evidence.to_csv(
        evidence_path, index=False, compression={"method": "gzip", "mtime": 0}
    )
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    selection.to_csv(
        selection_path, index=False, compression={"method": "gzip", "mtime": 0}
    )
    summary = {
        "contract": "iospe_orchid_colour_checkpoint_v1",
        "source_run_id": RUN_ID,
        "source_artifact": ARTIFACT,
        "fetch_manifest_sha256": _sha256(fetch_manifest_csv),
        "development_audit_sha256": (
            _sha256(development_audit_csv) if development_audit_csv else ""
        ),
        "fetched_pages": int(manifest["fetch_status"].eq("http:200").sum()),
        "failed_pages": int(manifest["fetch_status"].ne("http:200").sum()),
        "selected_evidence_rows": len(evidence),
        "selected_species_trait": len(
            evidence[["accepted_species", "trait_name"]].drop_duplicates()
        ),
        "selected_genera": int(
            evidence["accepted_species"].str.split().str[0].nunique()
        ),
        "audit_rows": len(audit),
        "guardrails": {
            "source_text_only": True,
            "exact_species_heading": True,
            "exact_quote": True,
            "cultivar_text_rejected": True,
            "family_inference": False,
            "global_fallback": False,
            "axis_only_join": False,
            "cross_trait_substitution": False,
        },
    }
    (output_dir / "iospe_orchid_colour_manifest_20260812.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


@app.command("record-review")
def record_review(
    audit_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_csv: Annotated[Path, typer.Option(dir_okay=False)],
    reviewer: Annotated[str, typer.Option()],
    reviewed_at_utc: Annotated[str, typer.Option()],
    rejected_id: Annotated[list[str] | None, typer.Option()] = None,
) -> None:
    """Record decisions after direct inspection of the deterministic sample."""

    audit = pd.read_csv(audit_csv, dtype=str).fillna("")
    rejected = set(rejected_id or [])
    if unknown := rejected.difference(audit["candidate_id"]):
        raise ValueError(f"unknown rejected candidate IDs: {sorted(unknown)}")
    accepted = ~audit["candidate_id"].isin(rejected)
    audit["accepted_correct"] = accepted.map({True: "true", False: "false"})
    audit["cultivar_status"] = "wild_or_unspecified_species_treatment"
    audit["reviewer"] = reviewer
    audit["reviewed_at_utc"] = reviewed_at_utc
    audit["audit_reason"] = accepted.map(
        {
            True: (
                "accepted_exact_species_heading_and_explicit_orchid_flower_colour_quote"
            ),
            False: "rejected_direct_treatment_review",
        }
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_csv, index=False, lineterminator="\n")


if __name__ == "__main__":
    app()
