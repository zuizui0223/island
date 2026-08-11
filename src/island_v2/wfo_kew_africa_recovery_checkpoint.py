"""Recover fail-closed direct floral evidence from the acquired Kew WFO flora wave.

The upstream acquisition deliberately retained broad lexical hits.  This
checkpoint does not fetch new pages or infer traits.  It re-audits the exact
species-description excerpts from Flora of Tropical East Africa, Flora of
West Tropical Africa, and Flora Zambesiaca, excludes already acquired
``accepted_species x trait_name`` cells, and emits the common source-package
schema plus a deterministic 200-record manual audit.

Form terms are accepted only when the nearest governing organ is a flower,
floret, corolla, or perianth.  This prevents shapes of calyces, hypanthia,
labella, and receptacle tubes from being transferred to the whole flower.
Every explicit supported state must equal the stored state set; incomplete
multistates fail closed.
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

from analysis.prepare_wfo_high_yield_source_package_20260810 import (
    DISPLAY_STATE_TERMS,
    ambiguous_display_quote,
)
from island_v2.all_evidence_trait_audit import canonical_trait_name
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS, TRAIT_TO_AXIS

app = typer.Typer(add_completion=False, no_args_is_help=True)

SOURCE_RUN_ID = "wfo-kew-africa-three-floras-20260810"
SOURCE_ARTIFACT = "wfo-kew-africa-official-dwca-20260810"
PACKAGE_RUN_ID = "wfo-kew-africa-recovery-20260811"
PACKAGE_ARTIFACT = "wfo-kew-africa-recovery-reviewed-source-package-20260811"
PRIOR_FORMAL_RUN_ID = "31440559943"
PRIOR_FORMAL_ARTIFACT = "reviewed-open-web-evidence-31440559943"
REVIEWER = "OpenAI Codex source-text holdout review"
REVIEWED_AT_UTC = "2026-08-11T00:00:00Z"
AUDIT_SEED = "wfo-kew-africa-recovery-holdout-20260811-v3"
AUDIT_QUOTA = {"floral_form": 50, "inflorescence_display": 150}
TARGET_PROVIDERS = frozenset(
    {
        "kew_flora_tropical_east_africa",
        "kew_flora_west_tropical_africa",
        "kew_flora_zambesiaca_current",
    }
)
TARGET_TRAITS = frozenset(AUDIT_QUOTA)
DIRECT_SCOPES = frozenset({"species_direct", "synonym_direct"})
STRICT_NAME_METHODS = frozenset(
    {
        "wfo_plant_list_exact_accepted_or_synonym_species_family_consistent",
        "wfo_june_2026_exact_accepted_species_family_consistent",
        "wfo_synonym_cluster_from_exact_fixed_name",
    }
)

FORM_TERMS = {
    "campanulate": "bell_campanulate",
    "bell-shaped": "bell_campanulate",
    "bell shaped": "bell_campanulate",
    "tubular": "tubular",
    "tubiform": "tubular",
    "tube-shaped": "tubular",
    "tube shaped": "tubular",
    "funnelform": "funnel_trumpet",
    "funnel-shaped": "funnel_trumpet",
    "funnel shaped": "funnel_trumpet",
    "trumpet-shaped": "funnel_trumpet",
    "trumpet shaped": "funnel_trumpet",
    "infundibuliform": "funnel_trumpet",
    "rotate": "open_radial",
    "stelliform": "open_radial",
    "star-shaped": "open_radial",
    "star shaped": "open_radial",
    "salverform": "salverform",
    "salver-shaped": "salverform",
    "hypocrateriform": "salverform",
    "urceolate": "urn_urceolate",
    "urn-shaped": "urn_urceolate",
    "urn shaped": "urn_urceolate",
    "bilabiate": "bilabiate",
    "two-lipped": "bilabiate",
    "2-lipped": "bilabiate",
    "two-labiate": "bilabiate",
    "2-labiate": "bilabiate",
    "papilionaceous": "papilionaceous",
    "spurred": "spurred",
    "brush-like": "brush_puff",
    "brushlike": "brush_puff",
    "puffball": "brush_puff",
    "powder-puff": "brush_puff",
    "powder puff": "brush_puff",
}
FORM_PATTERN = re.compile(
    "|".join(
        rf"(?<![A-Za-z]){re.escape(term)}(?![A-Za-z])"
        for term in sorted(FORM_TERMS, key=len, reverse=True)
    ),
    re.IGNORECASE,
)
POSITIVE_FORM_ORGAN = re.compile(
    r"\b(?:corollas?|perianths?|flowers?|florets?)\b", re.IGNORECASE
)
NEGATIVE_FORM_ORGAN = re.compile(
    r"\b(?:receptacle(?:[- ]tubes?)?|caly(?:x|ces)|hypanth(?:ium|ia)|"
    r"labellum|petals?|sepals?|bracts?|involucres?)\b",
    re.IGNORECASE,
)
INCOMPLETE_FORM_CONTEXT = re.compile(
    r"\bin bud\b|\bcup[- ]shaped\b", re.IGNORECASE
)
UNRESOLVED_RAY_DISK_FORM_CONTEXT = re.compile(
    r"\bray[- ]florets?\b", re.IGNORECASE
)
RAY_FLORETS_ABSENT = re.compile(
    r"\bray[- ]florets?\s+(?:0|absent|none)\b", re.IGNORECASE
)
EXTRA_DISPLAY_STATE_TERMS = {
    "solitary": re.compile(
        r"(?:\bflowers?\b[^;]{0,120}\bsolitary\b|"
        r"\bsolitary\b[^;]{0,60}\bflowers?\b)",
        re.IGNORECASE,
    ),
    "few_flowered": re.compile(
        r"(?:\b[1-5]\s*(?:-|\u2013|\u2014|to|or)\s*[1-5]"
        r"\s*[- ]flowered\b|\b[1-5]\s*[- ]flowered\b|"
        r"\(\s*[1-5]\s*-\s*\)\s*[1-5]\s*\(\s*-\s*[1-5]\s*\)"
        r"\s*[- ]flowered\b|"
        r"\b[1-5]\s*\(\s*[1-5]\s*\)\s*[- ]flowered\b|"
        r"\b(?:with|of)\s+[1-5]\s*(?:-|\u2013|\u2014|to|or)\s*[1-5]"
        r"\s+flowers?\b|\b(?:with|of)\s+[1-5]\s+flowers?\b|"
        r"\bflowers?\s+[1-5]\s*(?:-|\u2013|\u2014|to|or)\s*[1-5]\b|"
        r"\bflowers?\s+few\b)",
        re.IGNORECASE,
    ),
    "composite_display": re.compile(r"\bheads?\b", re.IGNORECASE),
    "umbel_corymb": re.compile(
        r"\b(?:subumbell\w*|subcorymb\w*)\b", re.IGNORECASE
    ),
}
UNMAPPED_SOLITARY_ALTERNATE = re.compile(
    r"(?:\bflowers?\b[^;]{0,100}\bsolitary\b[^;]{0,80}\b(?:or|and|rarely)\b"
    r"[^;]{0,80}\b(?:fascic\w*|cluster\w*|groups?|pairs?|paired|several|"
    r"ternate|geminate|[2-9](?:\s*(?:-|to|or)\s*[2-9])?)\b|"
    r"\b(?:fascic\w*|cluster\w*|groups?|pairs?|paired|several|ternate|geminate|"
    r"[2-9](?:\s*(?:-|to|or)\s*[2-9])?)\b[^;]{0,80}"
    r"\b(?:or|and)\b"
    r"[^;]{0,80}\bsolitary\b)",
    re.IGNORECASE,
)
UNMAPPED_FASCICLE_ALTERNATE = re.compile(
    r"(?:\bfascic\w*\b[^;]{0,80}\bor\b[^;]{0,80}"
    r"\b(?:racem\w*|panicul\w*|cym\w*|umbel\w*|spikes?)\b|"
    r"\b(?:racem\w*|panicul\w*|cym\w*|umbel\w*|spikes?)\b"
    r"[^;]{0,80}\bor\b[^;]{0,80}\bfascic\w*\b)",
    re.IGNORECASE,
)
INCOMPLETE_FLOWER_COUNT_VARIATION = re.compile(
    r"\bfew[- ](?:many|several)[- ]flowered\b|"
    r"\bfew\s*-?\s*(?:to|or)\s*many[- ]flowered\b|"
    r"\bfew-\s+and\b[^;]{0,25}\bflowered\b|"
    r"\b(?:sometimes|rarely|occasionally)\b[^;]{0,60}\bfew[- ]flowered\b",
    re.IGNORECASE,
)
CONTRADICTORY_FLOWER_COUNT_CONTEXT = re.compile(
    r"(?=.*\bfew[- ]flowered\b)(?=.*\bmany[- ]flowered\b)",
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


def observed_form_states(excerpt: object) -> set[str]:
    """Return only form states governed by a whole-flower organ.

    A 300-character window is intentional: flora descriptions contain unit
    abbreviations such as ``cm.`` that look like sentence boundaries.  The
    nearest positive versus negative organ within that window is more robust
    than splitting on punctuation, while remaining fail closed.
    """

    quote = _text(excerpt)
    observed: set[str] = set()
    for match in FORM_PATTERN.finditer(quote):
        prefix = quote[max(0, match.start() - 300) : match.start()]
        positive = [item.start() for item in POSITIVE_FORM_ORGAN.finditer(prefix)]
        negative = [item.start() for item in NEGATIVE_FORM_ORGAN.finditer(prefix)]
        if not positive or (negative and max(negative) > max(positive)):
            continue
        immediate = quote[max(0, match.start() - 20) : match.start()]
        if re.search(r"\b(?:not|never|non)[ -]*$", immediate, re.IGNORECASE):
            continue
        observed.add(FORM_TERMS[match.group(0).casefold()])
    return observed


def observed_display_states(excerpt: object) -> set[str]:
    quote = _text(excerpt)
    observed = {
        state for state, pattern in DISPLAY_STATE_TERMS.items() if pattern.search(quote)
    }
    observed.update(
        state
        for state, pattern in EXTRA_DISPLAY_STATE_TERMS.items()
        if pattern.search(quote)
    )
    return observed


def high_precision_candidate(record: dict[str, Any]) -> tuple[bool, str]:
    trait = canonical_trait_name(_text(record.get("trait_name")))
    quote = _text(record.get("source_excerpt"))
    expected = _state_set(record.get("normalized_value"))
    if _text(record.get("source_provider")) not in TARGET_PROVIDERS:
        return False, "provider_outside_checkpoint"
    if trait not in TARGET_TRAITS:
        return False, "trait_outside_checkpoint"
    if _text(record.get("source_run_id")) != SOURCE_RUN_ID:
        return False, "source_run_outside_checkpoint"
    if _text(record.get("source_artifact")) != SOURCE_ARTIFACT:
        return False, "source_artifact_outside_checkpoint"
    if _text(record.get("evidence_scope")).casefold() not in DIRECT_SCOPES:
        return False, "not_species_or_synonym_direct"
    if _text(record.get("name_match_method")) not in STRICT_NAME_METHODS:
        return False, "name_match_not_strict_wfo_exact"
    if not _text(record.get("accepted_species")) or not quote or not expected:
        return False, "missing_identity_quote_or_value"
    if trait == "floral_form":
        if INCOMPLETE_FORM_CONTEXT.search(quote):
            return False, "bud_or_unsupported_cup_form_context"
        if (
            UNRESOLVED_RAY_DISK_FORM_CONTEXT.search(quote)
            and not RAY_FLORETS_ABSENT.search(quote)
        ):
            return False, "disk_form_with_unresolved_ray_floret_form"
        observed = observed_form_states(quote)
    else:
        if (
            ambiguous_display_quote(quote)
            or UNMAPPED_SOLITARY_ALTERNATE.search(quote)
            or UNMAPPED_FASCICLE_ALTERNATE.search(quote)
            or INCOMPLETE_FLOWER_COUNT_VARIATION.search(quote)
            or CONTRADICTORY_FLOWER_COUNT_CONTEXT.search(quote)
        ):
            return False, "sex_specific_or_partial_inflorescence_context"
        observed = observed_display_states(quote)
    if not observed:
        return False, "no_safe_trait_grammar"
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
    return f"wfo-flora-treatment-content:{fingerprint}", fingerprint


def _audit_template(evidence: pd.DataFrame) -> pd.DataFrame:
    chosen: list[pd.DataFrame] = []
    used_urls: set[str] = set()
    for trait, quota in AUDIT_QUOTA.items():
        group = evidence.loc[evidence["trait_name"].eq(trait)].copy()
        group["audit_hash"] = group["source_record_id"].map(
            lambda value, name=trait: _id(AUDIT_SEED, name, value, length=64)
        )
        group = group.loc[~group["source_url"].isin(used_urls)]
        sample = (
            group.sort_values("audit_hash", kind="stable")
            .drop_duplicates("source_url")
            .head(quota)
        )
        if len(sample) != quota:
            raise ValueError(
                f"{trait} audit requires {quota} unique pages, found {len(sample)}"
            )
        used_urls.update(sample["source_url"])
        chosen.append(sample)
    sample = pd.concat(chosen, ignore_index=True).sort_values(
        ["trait_name", "audit_hash"], kind="stable"
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
        source_record_id = "wfo-africa-recovery:" + _id(
            species, trait, value, lineage
        )
        evidence_rows.append(
            {
                "accepted_species": species,
                "axis": TRAIT_TO_AXIS[trait],
                "trait_name": trait,
                "normalized_value": value,
                "quality": "high",
                "source_group": "wfo_kew_official_flora_recovery",
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
                "source_run_id": SOURCE_RUN_ID,
                "source_artifact": SOURCE_ARTIFACT,
                "source_file": "wfo_kew_africa_evidence_candidates.csv.gz",
                "acceptance_contract": (
                    "wfo_official_flora_exact_quote_recovery_reaudit_v1"
                ),
                "content_fingerprint": fingerprint,
            }
        )

    evidence = pd.DataFrame(evidence_rows)
    if evidence.empty:
        raise ValueError("WFO Africa recovery selected no evidence")
    evidence = evidence.sort_values(
        ["trait_name", "accepted_species", "source_provider", "source_record_id"],
        kind="stable",
    ).drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"],
        keep="first",
    )
    evidence = evidence.reset_index(drop=True)
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("WFO Africa recovery generated duplicate source_record_id values")
    audit = _audit_template(evidence)
    selection = pd.DataFrame(selection_rows).sort_values(
        ["selected", "trait_name", "accepted_species", "original_source_record_id"],
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
    evidence_path = output_dir / "wfo_kew_africa_recovery_evidence_20260811.csv.gz"
    audit_path = output_dir / "wfo_kew_africa_recovery_audit_200_20260811.csv"
    selection_path = output_dir / "wfo_kew_africa_recovery_selection_20260811.csv.gz"
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
        "contract": "wfo_kew_africa_recovery_checkpoint_v1",
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "package_run_id": PACKAGE_RUN_ID,
        "package_artifact": PACKAGE_ARTIFACT,
        "prior_formal_run_id": PRIOR_FORMAL_RUN_ID,
        "prior_formal_artifact": PRIOR_FORMAL_ARTIFACT,
        "providers": sorted(TARGET_PROVIDERS),
        "traits": sorted(TARGET_TRAITS),
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
        "trait_rows": evidence["trait_name"].value_counts().sort_index().to_dict(),
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
    (output_dir / "wfo_kew_africa_recovery_manifest_20260811.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    app()
