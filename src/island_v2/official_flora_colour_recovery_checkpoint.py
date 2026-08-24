"""Recover direct flower-colour evidence from acquired official flora text.

The broad flora acquisition retained lexical colour hits so that later gates
could be improved without refetching pages.  This checkpoint is the strict,
organ-aware re-audit of that cache.  A colour is retained only when its local
grammatical subject is a flower, floret, corolla, petal, perianth, tepal, or
ligule.  Colours of hairs, glands, bracts, calyces, fruits, foliage, and other
non-target structures are rejected.  The complete observed colour-state set
must equal the stored state set, so partially flattened multicolour records
fail closed.

No genus, family, or global inference is performed here.  The output uses the
common reviewed source-package contract and a deterministic 200-page audit.
"""

from __future__ import annotations

import hashlib
import json
import re
import unicodedata
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.all_evidence_trait_audit import canonical_trait_name
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS, TRAIT_TO_AXIS

app = typer.Typer(add_completion=False, no_args_is_help=True)

PACKAGE_RUN_ID = "official-flora-colour-recovery-20260811"
PACKAGE_ARTIFACT = "official-flora-colour-recovery-reviewed-source-package-20260811"
PRIOR_FORMAL_RUN_ID = "31445281207"
PRIOR_FORMAL_ARTIFACT = "reviewed-open-web-evidence-31445281207"
REVIEWER = "OpenAI Codex source-text holdout review"
REVIEWED_AT_UTC = "2026-08-11T00:00:00Z"
AUDIT_SEED = "official-flora-colour-recovery-holdout-20260811-v1"
AUDIT_SIZE = 200
TRAIT = "flower_primary_color"

TARGET_PROVIDERS = frozenset(
    {
        "brazilian_flora_2020",
        "flora_of_australia_official",
        "flora_of_panama",
        "kew_flora_tropical_east_africa",
        "kew_flora_west_tropical_africa",
        "kew_flora_zambesiaca_current",
        "nybg_floraneotropica",
        "nybg_memoirs",
        "plantnet_nsw_flora_source_scale",
    }
)
DIRECT_SCOPES = frozenset({"species_direct", "synonym_direct"})
STRICT_NAME_METHODS = frozenset(
    {
        "accepted_name_exact_plus_family_exact",
        "wfo_accepted_usage_from_exact_fixed_synonym",
        "wfo_june_2026_exact_accepted_species_family_consistent",
        "wfo_june_2026_exact_synonym_species_family_consistent",
        "wfo_plant_list_exact_accepted_or_synonym_species_family_consistent",
        "wfo_synonym_cluster_from_exact_fixed_name",
    }
)

COLOUR_TERMS: dict[str, tuple[str, ...]] = {
    "white": (
        "white",
        "whitish",
        "cream",
        "creamy",
        "creamish",
        "ivory",
    ),
    "yellow_orange": (
        "yellow",
        "yellowish",
        "orange",
        "golden",
        "ochre",
        "ochreous",
    ),
    "red_pink": (
        "red",
        "reddish",
        "pink",
        "pinkish",
        "scarlet",
        "crimson",
        "rose",
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
        "inconspicuous",
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

# Only petaloid/whole-flower subjects define the strict flower-colour trait.
# Sepals and calyces are intentionally excluded; perianth and tepals remain
# eligible for taxa without differentiated petals.
POSITIVE_ORGAN = re.compile(
    r"\b(?:flowers?(?![- ]heads?)|florets?|corollas?|petals?|perianths?|tepals?|ligules?|"
    r"standards?|wings?|keel[- ]petals?|crests?)\b",
    re.IGNORECASE,
)
NEGATIVE_ORGAN = re.compile(
    r"\b(?:leaves?|leaflets?|foliage|fruits?|berries|seeds?|stems?|bark|"
    r"branches?|branchlets?|rachis|rachillae?|peduncles?|pedicels?|bracts?|"
    r"involucres?|paleae?|sepals?|caly(?:x|ces)|hypanth(?:ium|ia)|ovaries|"
    r"ovary|carpels?|anthers?|stamens?|filaments?|styles?|stigmas?|glands?|hairs?|"
    r"trichomes?|indumentum|pubescence|spines?|prickles?|glumes?|lemmas?|"
    r"spikelets?|sheaths?|tomentum|inflorescences?|infructescences?|drupes?|"
    r"buds?|aerial parts?|subterranean parts?|pulp)\b",
    re.IGNORECASE,
)
NEGATIVE_AFTER = re.compile(
    r"^[\s,;:/-]{0,8}(?:[A-Za-z-]+\s+){0,2}"
    r"(?:glands?|hairs?|trichomes?|tomentum|pubescence|pubescent|pilose|"
    r"puberulous|puberulent|scurfy|scales?|tomentose|velvety|hairy|ciliate|spines?|prickles?|"
    r"bracts?|sepals?|caly(?:x|ces)|anthers?|glumes?|lemmas?|fruits?|berries)\b",
    re.IGNORECASE,
)
NEGATIVE_BEFORE = re.compile(
    r"\b(?:glands?|hairs?|trichomes?|tomentum|pubescence|pubescent|pilose|"
    r"puberulous|puberulent|scurfy|scales?|tomentose|velvety|hairy|ciliate|spines?|prickles?|"
    r"bracts?|sepals?|caly(?:x|ces)|anthers?|glumes?|lemmas?|fruits?|berries)"
    r"(?:\s+[A-Za-z-]+){0,2}[\s,;:/-]{0,8}$",
    re.IGNORECASE,
)
NEGATED = re.compile(r"\b(?:not|never|without|non)[ -]*$", re.IGNORECASE)
CULTIVAR_OR_HYBRID = re.compile(
    r"\b(?:cultivar|horticultural|garden hybrid|hybrid cultivar|cv\.)\b",
    re.IGNORECASE,
)
COMPARATIVE_TAXON = re.compile(
    r"\b(?:unlike|as in|other species|related species|allied species|"
    r"distinguished from)\b",
    re.IGNORECASE,
)
BUD_CONTEXT = re.compile(r"\b(?:in bud|buds?)\b", re.IGNORECASE)
DRY_CONTEXT = re.compile(
    r"\b(?:when dry|when dried|drying|dried|in sicco)\b", re.IGNORECASE
)
INVALID_POSITIVE_PREFIX = re.compile(
    r"\b(?:fruiting|with(?:\s+the)?|among|between|subtending|beneath|"
    r"producing|as\s+long\s+as(?:\s+the)?|longer\s+than(?:\s+the)?|per)\s+$",
    re.IGNORECASE,
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _id(*parts: object, length: int = 24) -> str:
    payload = "\n".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]


def _normalized_excerpt(value: object) -> str:
    text = unicodedata.normalize("NFKC", _text(value)).casefold()
    return " ".join(re.findall(r"[\w]+", text, flags=re.UNICODE))


def _state_set(value: object) -> set[str]:
    return {token for token in _text(value).split("|") if token}


def _repair_colour_hyphenation(text: str) -> str:
    """Repair OCR/PDF line-break hyphens inside known colour words only."""

    repaired = text
    for term in sorted(TERM_TO_STATE, key=len, reverse=True):
        if "-" in term or " " in term:
            continue
        for split in range(2, len(term) - 1):
            pattern = re.compile(
                rf"(?<![A-Za-z]){re.escape(term[:split])}-\s+"
                rf"{re.escape(term[split:])}(?![A-Za-z])",
                re.IGNORECASE,
            )
            repaired = pattern.sub(term, repaired)
    return repaired


def _clause_bounds(text: str, position: int) -> tuple[int, int]:
    """Return a flora-description clause without splitting on abbreviations."""

    left = max(text.rfind(";", 0, position), text.rfind(".", 0, position)) + 1
    semicolon = text.find(";", position)
    period = text.find(".", position)
    ends = [value for value in (semicolon, period) if value >= 0]
    right = min(ends) if ends else len(text)
    return left, right


def _organ_distance(match: re.Match[str], position: int) -> int:
    if match.end() <= position:
        return position - match.end()
    return match.start() - position


def _linked_to_positive_organ(text: str, colour: re.Match[str]) -> bool:
    start, end = _clause_bounds(text, colour.start())
    clause = text[start:end]
    local_position = colour.start() - start
    before = clause[max(0, local_position - 60) : local_position]
    after = clause[colour.end() - start : colour.end() - start + 60]
    if NEGATED.search(before):
        return False
    if NEGATIVE_BEFORE.search(before) or NEGATIVE_AFTER.search(after):
        return False

    positive: list[re.Match[str]] = []
    for item in POSITIVE_ORGAN.finditer(clause):
        prefix = clause[max(0, item.start() - 30) : item.start()]
        if INVALID_POSITIVE_PREFIX.search(prefix):
            continue
        if re.search(r"\bof(?:\s+the)?\s+$", prefix, re.IGNORECASE):
            prior = clause[: item.start()]
            if NEGATIVE_ORGAN.search(prior):
                continue
        positive.append(item)
    negative = list(NEGATIVE_ORGAN.finditer(clause))
    if not positive:
        return False
    nearest_positive = min(_organ_distance(item, local_position) for item in positive)
    nearest_negative = min(
        (_organ_distance(item, local_position) for item in negative),
        default=10_000,
    )
    if nearest_positive > 180 or nearest_negative < nearest_positive:
        return False

    # ``petals ... mauve in bud`` describes a transient bud state, not the
    # open flower.  Other colours in the same clause can still be retained.
    bud = BUD_CONTEXT.search(clause)
    if bud and abs(bud.start() - local_position) <= 45:
        return False
    dry = DRY_CONTEXT.search(clause)
    return not (dry and abs(dry.start() - local_position) <= 55)


def observed_colour_states(excerpt: object) -> set[str]:
    """Return complete floral colour states linked to eligible floral organs."""

    quote = _repair_colour_hyphenation(_text(excerpt))
    observed: set[str] = set()
    for match in COLOUR_PATTERN.finditer(quote):
        if _linked_to_positive_organ(quote, match):
            observed.add(TERM_TO_STATE[match.group(0).casefold()])
    return observed


def high_precision_candidate(record: dict[str, Any]) -> tuple[bool, str]:
    trait = canonical_trait_name(_text(record.get("trait_name")))
    quote = _text(record.get("source_excerpt"))
    expected = _state_set(record.get("normalized_value"))
    if _text(record.get("source_provider")) not in TARGET_PROVIDERS:
        return False, "provider_outside_checkpoint"
    if trait != TRAIT:
        return False, "trait_outside_checkpoint"
    if _text(record.get("evidence_scope")).casefold() not in DIRECT_SCOPES:
        return False, "not_species_or_synonym_direct"
    if _text(record.get("name_match_method")) not in STRICT_NAME_METHODS:
        return False, "name_match_not_strict_exact"
    if not _text(record.get("accepted_species")) or not quote or not expected:
        return False, "missing_identity_quote_or_value"
    if CULTIVAR_OR_HYBRID.search(quote):
        return False, "cultivar_or_horticultural_hybrid_context"
    if COMPARATIVE_TAXON.search(quote):
        return False, "comparative_taxon_context"
    observed = observed_colour_states(quote)
    if not observed:
        return False, "no_organ_linked_flower_colour"
    if observed != expected:
        return False, "explicit_state_set_mismatch"
    return True, "selected_for_manual_audit_gate"


def _completed_pairs(*frames: pd.DataFrame) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for frame in frames:
        if missing := {"accepted_species", "trait_name"}.difference(frame.columns):
            raise ValueError(f"baseline evidence missing columns: {sorted(missing)}")
        work = frame.copy().fillna("")
        if "resolution_status" in work:
            work = work.loc[work["resolution_status"].eq("resolved")]
        work["trait_name"] = work["trait_name"].map(canonical_trait_name)
        pairs.update(zip(work["accepted_species"], work["trait_name"], strict=True))
    return pairs


def _content_lineage(species: str, excerpt: str) -> tuple[str, str]:
    fingerprint = _id(species, _normalized_excerpt(excerpt), length=40)
    return f"official-flora-treatment-content:{fingerprint}", fingerprint


def _audit_template(evidence: pd.DataFrame) -> pd.DataFrame:
    sample = evidence.copy()
    sample["audit_hash"] = sample["source_record_id"].map(
        lambda value: _id(AUDIT_SEED, value, length=64)
    )
    sample = (
        sample.sort_values("audit_hash", kind="stable")
        .drop_duplicates("source_url")
        .head(AUDIT_SIZE)
    )
    if len(sample) != AUDIT_SIZE:
        raise ValueError(
            f"flower colour audit requires {AUDIT_SIZE} unique pages, found {len(sample)}"
        )
    return pd.DataFrame(
        {
            "candidate_id": sample["source_record_id"],
            "trait_name": sample["trait_name"],
            "accepted_species": sample["accepted_species"],
            "normalized_value": sample["normalized_value"],
            "source_provider": sample["source_provider"],
            "source_url": sample["source_url"],
            "source_citation": sample["source_citation"],
            "source_excerpt": sample["source_excerpt"],
            "content_fingerprint": sample["content_fingerprint"],
            "accepted_correct": "",
            "cultivar_status": "",
            "reviewer": "",
            "reviewed_at_utc": "",
            "audit_reason": "",
            "audit_hash": sample["audit_hash"],
        }
    ).reset_index(drop=True)


def apply_reviewed_audit(template: pd.DataFrame, reviewed: pd.DataFrame) -> pd.DataFrame:
    immutable = [
        "candidate_id",
        "trait_name",
        "accepted_species",
        "normalized_value",
        "source_provider",
        "source_url",
        "source_citation",
        "source_excerpt",
        "content_fingerprint",
        "audit_hash",
    ]
    decisions = [
        "accepted_correct",
        "cultivar_status",
        "reviewer",
        "reviewed_at_utc",
        "audit_reason",
    ]
    if missing := set(immutable + decisions).difference(reviewed.columns):
        raise ValueError(f"reviewed audit missing columns: {sorted(missing)}")
    if reviewed["candidate_id"].map(_text).duplicated().any():
        raise ValueError("reviewed audit has duplicate candidate_id values")
    if set(template["candidate_id"]) != set(reviewed["candidate_id"]):
        raise ValueError("reviewed audit candidate IDs do not match deterministic sample")
    aligned = reviewed.set_index("candidate_id").loc[template["candidate_id"]].reset_index()
    for column in immutable:
        if not template[column].map(_text).eq(aligned[column].map(_text)).all():
            raise ValueError(f"reviewed audit changed immutable column: {column}")
    if aligned[decisions].apply(lambda column: column.map(_text).eq("").any()).any():
        raise ValueError("reviewed audit has incomplete decisions or provenance")
    return aligned.reindex(columns=template.columns)


def build_checkpoint(
    candidates: pd.DataFrame,
    baseline_direct: pd.DataFrame,
    prior_public_web: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return evidence, deterministic audit template, and selection audit."""

    if missing := set(EVIDENCE_COLUMNS).difference(candidates.columns):
        raise ValueError(f"candidate cache missing columns: {sorted(missing)}")
    completed = _completed_pairs(baseline_direct, prior_public_web)
    selection_rows: list[dict[str, str]] = []
    evidence_rows: list[dict[str, str]] = []
    for record in candidates.fillna("").to_dict("records"):
        trait = canonical_trait_name(_text(record.get("trait_name")))
        species = _text(record.get("accepted_species"))
        selected, reason = high_precision_candidate(record)
        if selected and (species, trait) in completed:
            selected, reason = False, "already_acquired_direct_species_trait"
        selection_rows.append(
            {
                "original_source_record_id": _text(record.get("source_record_id")),
                "accepted_species": species,
                "trait_name": trait,
                "normalized_value": _text(record.get("normalized_value")),
                "source_provider": _text(record.get("source_provider")),
                "selected": str(selected).lower(),
                "selection_reason": reason,
            }
        )
        if not selected:
            continue
        quote = _text(record.get("source_excerpt"))
        lineage, fingerprint = _content_lineage(species, quote)
        value = "|".join(sorted(_state_set(record.get("normalized_value"))))
        source_record_id = "official-flora-colour-recovery:" + _id(
            species, trait, value, lineage
        )
        evidence_rows.append(
            {
                "accepted_species": species,
                "axis": TRAIT_TO_AXIS[trait],
                "trait_name": trait,
                "normalized_value": value,
                "quality": "high",
                "source_group": "official_flora_colour_recovery",
                "source_provider": _text(record.get("source_provider")),
                "source_url": _text(record.get("source_url")),
                "source_record_id": source_record_id,
                "source_citation": _text(record.get("source_citation")),
                "source_excerpt": quote,
                "evidence_scope": _text(record.get("evidence_scope")),
                "name_match_method": _text(record.get("name_match_method")),
                "source_lineage": lineage,
                "lineage_method": (
                    "accepted_species_plus_normalized_excerpt_content_fingerprint"
                ),
                "source_run_id": _text(record.get("source_run_id")),
                "source_artifact": _text(record.get("source_artifact")),
                "source_file": _text(record.get("source_file")),
                "acceptance_contract": (
                    "official_flora_organ_linked_complete_colour_state_reaudit_v1"
                ),
                "content_fingerprint": fingerprint,
            }
        )

    evidence = pd.DataFrame(evidence_rows)
    if evidence.empty:
        raise ValueError("official flora colour recovery selected no evidence")
    evidence = evidence.sort_values(
        ["accepted_species", "source_provider", "source_record_id"], kind="stable"
    ).drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"],
        keep="first",
    )
    evidence = evidence.reset_index(drop=True)
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("official flora recovery generated duplicate source_record_id")
    audit = _audit_template(evidence)
    selection = pd.DataFrame(selection_rows).sort_values(
        ["selected", "source_provider", "accepted_species", "original_source_record_id"],
        kind="stable",
    )
    return evidence, audit, selection.reset_index(drop=True)


@app.command("build")
def build(
    candidates_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    baseline_direct_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    prior_public_web_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    reviewed_audit_csv: Annotated[
        Path | None, typer.Option(exists=True, dir_okay=False)
    ] = None,
) -> None:
    candidates = pd.read_csv(candidates_csv, dtype=str).fillna("")
    baseline = pd.read_csv(baseline_direct_csv, dtype=str).fillna("")
    prior = pd.read_csv(prior_public_web_csv, dtype=str).fillna("")
    evidence, audit, selection = build_checkpoint(candidates, baseline, prior)
    if reviewed_audit_csv is not None:
        reviewed = pd.read_csv(reviewed_audit_csv, dtype=str).fillna("")
        audit = apply_reviewed_audit(audit, reviewed)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "official_flora_colour_evidence_20260811.csv.gz"
    audit_path = output_dir / "official_flora_colour_audit_200_20260811.csv"
    selection_path = output_dir / "official_flora_colour_selection_20260811.csv.gz"
    evidence.reindex(columns=[*EVIDENCE_COLUMNS, "content_fingerprint"]).to_csv(
        evidence_path,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    selection.to_csv(
        selection_path,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )
    manifest: dict[str, Any] = {
        "contract": "official_flora_colour_recovery_checkpoint_v1",
        "package_run_id": PACKAGE_RUN_ID,
        "package_artifact": PACKAGE_ARTIFACT,
        "prior_formal_run_id": PRIOR_FORMAL_RUN_ID,
        "prior_formal_artifact": PRIOR_FORMAL_ARTIFACT,
        "providers": sorted(evidence["source_provider"].unique()),
        "source_runs": sorted(evidence["source_run_id"].unique()),
        "source_artifacts": sorted(evidence["source_artifact"].unique()),
        "traits": [TRAIT],
        "quality": "high",
        "family_inference": False,
        "global_fallback": False,
        "trait_substitution": False,
        "candidate_rows": len(evidence),
        "species_traits": len(
            evidence[["accepted_species", "trait_name"]].drop_duplicates()
        ),
        "species_axes": len(evidence[["accepted_species", "axis"]].drop_duplicates()),
        "species": evidence["accepted_species"].nunique(),
        "provider_rows": evidence["source_provider"].value_counts().sort_index().to_dict(),
        "audit_rows": len(audit),
        "audit_unique_urls": audit["source_url"].nunique(),
        "audit_review_complete": bool(
            audit[["accepted_correct", "reviewer", "reviewed_at_utc"]]
            .apply(lambda column: column.map(_text).ne("").all())
            .all()
        ),
        "inputs": {
            str(candidates_csv): _sha256(candidates_csv),
            str(baseline_direct_csv): _sha256(baseline_direct_csv),
            str(prior_public_web_csv): _sha256(prior_public_web_csv),
        },
    }
    manifest["files"] = [
        {
            "path": path.name,
            "sha256": _sha256(path),
            "size_bytes": path.stat().st_size,
        }
        for path in (evidence_path, audit_path, selection_path)
    ]
    (output_dir / "official_flora_colour_manifest_20260811.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    app()
