"""Build a fail-closed Medium checkpoint from an acquired eFloras artifact.

The upstream eFloras extractor intentionally keeps broad rule hits as pending.
This module narrows those already-acquired rows to grammatical contexts that
refer to the flower/corolla/perianth itself, excludes completed direct cells,
and emits the common source-package schema plus a deterministic 200-row audit.
It never fetches pages and never promotes a pending row without the separate
manual audit gate in :mod:`island_v2.open_web_common`.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.angiosperm_scope import classify_scope, load_config
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS, TRAIT_TO_AXIS

app = typer.Typer(add_completion=False)

SOURCE_RUN_ID = "30688236057"
SOURCE_ARTIFACT = "efloras-all-floras-combined"
SOURCE_ARTIFACT_ID = "8818398786"
TARGET_TRAITS = {
    "floral_form",
    "flower_size_class",
    "inflorescence_display",
    "tube_depth_class",
}
AUDIT_ALLOCATION = {
    "flower_size_class": 69,
    "inflorescence_display": 60,
    "floral_form": 31,
    "tube_depth_class": 40,
}

_CULTIVAR = re.compile(
    r"\b(?:cultivar|cultivars|cv\.|horticultural|double[- ]flowered|garden hybrid)\b",
    re.IGNORECASE,
)
_SUBORDINATE_ORGAN = re.compile(
    r"\b(?:tubes?|lobes?|petals?|sepals?|calyces|calyx|pedicels?|peduncles?|"
    r"bracts?|anthers?|filaments?|stamens?|styles?|ovaries|ovary|spikelets?|"
    r"discs?|florets?|tepals?|hypanthia|hypanthium|inflorescences?|nectaries|nectary|"
    r"corona|appendices?|mouths?|neuter)\b",
    re.IGNORECASE,
)
_TUBE_SUBJECT_SWITCH = re.compile(
    r"\b(?:stamens?|filaments?|scales?|lobes?|inserted|bridged|fruiting|in fruit|"
    r"wide|diameter)\b",
    re.IGNORECASE,
)
_FORM_TERMS = {
    "campanulate": "bell_campanulate",
    "bell-shaped": "bell_campanulate",
    "bell shaped": "bell_campanulate",
    "tubular": "tubular",
    "tube-shaped": "tubular",
    "tube shaped": "tubular",
    "funnelform": "funnel_trumpet",
    "funnel-shaped": "funnel_trumpet",
    "trumpet-shaped": "funnel_trumpet",
    "trumpet shaped": "funnel_trumpet",
    "infundibuliform": "funnel_trumpet",
    "rotate": "open_radial",
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
    "papilionaceous": "papilionaceous",
    "spurred": "spurred",
    "brush-like": "brush_puff",
    "brushlike": "brush_puff",
    "puffball": "brush_puff",
    "powder-puff": "brush_puff",
    "powder puff": "brush_puff",
}
_FORM_PATTERN = "|".join(
    re.escape(term) for term in sorted(_FORM_TERMS, key=len, reverse=True)
)
_FORM_VALUE_ALIASES = {
    "campanulate": "bell_campanulate",
    "funnelform": "funnel_trumpet",
    "rotate": "open_radial",
    "urceolate": "urn_urceolate",
}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _measurement_class(low: str, high: str, unit: str, trait: str) -> str:
    midpoint = (float(low) + float(high or low)) / 2.0
    mm = midpoint * (10.0 if unit.casefold() == "cm" else 1.0)
    thresholds = (
        (
            (2.0, "very_small"),
            (10.0, "small"),
            (30.0, "medium"),
            (100.0, "large"),
            (float("inf"), "very_large"),
        )
        if trait == "flower_size_class"
        else (
            (1.0, "absent_or_open"),
            (5.0, "shallow"),
            (15.0, "intermediate"),
            (float("inf"), "deep"),
        )
    )
    return next(label for maximum, label in thresholds if mm <= maximum)


def _safe_size_values(excerpt: str) -> set[str]:
    values: set[str] = set()
    two_dimensional = re.compile(
        r"\b(flowers?|corollas?|perianths?)\b"
        r"(?P<between>[^;.!?]{0,100}?)"
        r"(?P<low>\d+(?:\.\d+)?)\s*"
        r"(?:(?:-|\u2013|\u2014)\s*(?P<high>\d+(?:\.\d+)?))?\s*"
        r"(?:x|\u00d7)\s*\d+(?:\.\d+)?\s*"
        r"(?:(?:-|\u2013|\u2014)\s*\d+(?:\.\d+)?)?\s*"
        r"(?P<unit>mm|cm)\b",
        re.IGNORECASE,
    )
    for match in two_dimensional.finditer(excerpt):
        sentence_prefix = re.split(r"[;.!?]", excerpt[: match.start()])[-1]
        if re.search(
            r"\b(?:from|between|in|neuter)\s+$", sentence_prefix, re.IGNORECASE
        ):
            continue
        if re.search(
            r"\b(?:inflorescences?|spikes?|racemes?|panicles?|cymes?|thyrses?)\b",
            sentence_prefix,
            re.IGNORECASE,
        ):
            continue
        if _SUBORDINATE_ORGAN.search(match.group("between")):
            continue
        if re.match(
            r"\s*(?:pedicels?|peduncles?|tepals?|petals?|sepals?|lobes?|spurs?)\b",
            excerpt[match.end() : match.end() + 30],
            re.IGNORECASE,
        ):
            continue
        values.add(
            _measurement_class(
                match.group("low"),
                match.group("high") or "",
                match.group("unit"),
                "flower_size_class",
            )
        )
    if values:
        return values
    pattern = re.compile(
        r"\b(flowers?|corollas?|perianths?)\b"
        r"(?P<between>[^;.!?]{0,100}?)"
        r"(?P<low>\d+(?:\.\d+)?)\s*"
        r"(?:[-–—]\s*(?P<high>\d+(?:\.\d+)?))?\s*"
        r"(?P<unit>mm|cm)\b",
        re.IGNORECASE,
    )
    for match in pattern.finditer(excerpt):
        sentence_prefix = re.split(r"[;.!?]", excerpt[: match.start()])[-1]
        if re.search(
            r"\b(?:from|between|in|neuter)\s+$", sentence_prefix, re.IGNORECASE
        ):
            continue
        if re.search(
            r"\b(?:inflorescences?|spikes?|racemes?|panicles?|cymes?|thyrses?)\b",
            sentence_prefix,
            re.IGNORECASE,
        ):
            continue
        if _SUBORDINATE_ORGAN.search(match.group("between")):
            continue
        if re.match(
            r"\s*(?:pedicels?|peduncles?|tepals?|petals?|sepals?|lobes?|spurs?)\b",
            excerpt[match.end() : match.end() + 30],
            re.IGNORECASE,
        ):
            continue
        values.add(
            _measurement_class(
                match.group("low"),
                match.group("high") or "",
                match.group("unit"),
                "flower_size_class",
            )
        )
    return values


def _safe_tube_values(excerpt: str) -> set[str]:
    values: set[str] = set()
    pattern = re.compile(
        r"\b(?:corolla|perianth)\s+tube\b"
        r"(?P<between>[^;.!?]{0,45}?)"
        r"(?P<low>\d+(?:\.\d+)?)\s*"
        r"(?:[-–—]\s*(?P<high>\d+(?:\.\d+)?))?\s*"
        r"(?P<unit>mm|cm)\b",
        re.IGNORECASE,
    )
    for match in pattern.finditer(excerpt):
        before = excerpt[max(0, match.start() - 100) : match.start()]
        after = excerpt[match.end() : match.end() + 15]
        if _TUBE_SUBJECT_SWITCH.search(before + " " + match.group("between")):
            continue
        if re.search(r"\b(?:wide|diam(?:eter)?\.?)\b", after, re.IGNORECASE):
            continue
        values.add(
            _measurement_class(
                match.group("low"),
                match.group("high") or "",
                match.group("unit"),
                "tube_depth_class",
            )
        )
    return values


def _safe_form_values(excerpt: str) -> set[str]:
    values: set[str] = set()
    pattern = re.compile(
        r"\b(?P<organ>corolla|perianth|flowers?)\b"
        r"(?P<between>[^;.!?]{0,80}?)"
        rf"(?P<form>{_FORM_PATTERN})\b",
        re.IGNORECASE,
    )
    for match in pattern.finditer(excerpt):
        between = match.group("between")
        before_organ = excerpt[max(0, match.start() - 25) : match.start()]
        if re.search(r"\b(?:than|in)\s+$", before_organ, re.IGNORECASE):
            continue
        if match.group("organ").casefold().startswith("flower") and re.search(
            r"\b(?:calyces|calyx|hypanthia|hypanthium|involucres?|bracts?)\b",
            between,
            re.IGNORECASE,
        ):
            continue
        prefix = excerpt[max(0, match.start("form") - 20) : match.start("form")]
        if re.search(r"\b(?:not|never|non)[ -]*$", prefix, re.IGNORECASE):
            continue
        values.add(_FORM_TERMS[match.group("form").casefold()])
    return values


def _normalised_candidate_values(trait_name: str, value: str) -> set[str]:
    values = {_text(token) for token in str(value).split("|") if _text(token)}
    if trait_name == "floral_form":
        return {_FORM_VALUE_ALIASES.get(token, token) for token in values}
    return values


def high_precision_candidate(record: dict[str, Any]) -> tuple[bool, str]:
    """Return whether one acquired pending row passes the narrow grammar gate."""

    trait = _text(record.get("trait_name"))
    excerpt = _text(record.get("source_excerpt"))
    value = _text(record.get("candidate_value"))
    if trait not in TARGET_TRAITS:
        return False, "trait_outside_checkpoint"
    if _CULTIVAR.search(excerpt):
        return False, "cultivar_or_horticultural_context"
    if re.search(
        r"\b(?:treated\b[^.]{0,80}\bas a species|does not agree with|"
        r"circumscription of|which some authors have included in|is similar but)\b",
        excerpt,
        re.IGNORECASE,
    ):
        return False, "taxonomic_comparative_note"
    expected = _normalised_candidate_values(trait, value)
    if trait == "flower_size_class":
        if re.search(r"\bor much (?:smaller|larger)\b", excerpt, re.IGNORECASE):
            return False, "incomplete_multistate_context"
        if re.search(r"\bflowers?\b[^;.!?]{0,50}\bin bud\b", excerpt, re.IGNORECASE):
            return False, "bud_measurement_not_open_flower"
        observed = _safe_size_values(excerpt)
    elif trait == "tube_depth_class":
        observed = _safe_tube_values(excerpt)
    elif trait == "floral_form":
        observed = _safe_form_values(excerpt)
    else:
        if re.search(r"\b(?:pollen|seed)[ -]?cones?\b", excerpt, re.IGNORECASE):
            return False, "non_flower_cone_context"
        if re.search(
            r"\b(?:in contrast|unlike|compared with)\b", excerpt, re.IGNORECASE
        ):
            return False, "comparative_other_taxon_context"
        direct_start = re.match(
            r"^(?:[^.;:]{0,25}:\s*)?"
            r"(?:inflorescences?|flowers?|capitula|synflorescences?|umbels?|"
            r"racemes?|panicles?|spikes?|male inflorescences?|female inflorescences?|"
            r"pistillate flowers?|staminate flowers?)\b",
            excerpt,
            re.IGNORECASE,
        )
        if not direct_start:
            return False, "not_direct_morphology_fragment"
        if "solitary" in expected and re.search(
            r"(?:\bsolitary\b[^;]{0,45}\bor\b[^;]{0,45}"
            r"(?:pairs?|paired|together|clusters?|fascicles?|cymes?|racemes?|"
            r"panicles?|spikes?|umbels?|verticils?|branched|\d)|"
            r"(?:pairs?|paired|together|clusters?|fascicles?|cymes?|racemes?|"
            r"panicles?|spikes?|umbels?|verticils?|branched|\d)"
            r"[^;]{0,45}\bor\b[^;]{0,45}\bsolitary\b)",
            excerpt,
            re.IGNORECASE,
        ):
            return False, "incomplete_multistate_context"
        if "solitary" in expected and re.search(
            r"\b(?:rarely|occasionally)\s+(?:in\s+)?(?:pairs?|clusters?|fascicles?|"
            r"cymes?|racemes?|panicles?|spikes?|umbels?|verticils?)\b|"
            r"\bcyath(?:ium|ia)\b",
            excerpt,
            re.IGNORECASE,
        ):
            return False, "incomplete_multistate_context"
        observed = expected
    if not observed:
        return False, "no_safe_trait_grammar"
    if observed != expected:
        return False, "safe_grammar_value_mismatch"
    return True, "selected_for_manual_audit_gate"


def _candidate_id(record: dict[str, Any]) -> str:
    payload = "|".join(
        _text(record.get(column))
        for column in (
            "flora_id",
            "taxon_id",
            "accepted_species",
            "trait_name",
            "candidate_value",
            "source_excerpt",
        )
    )
    return "efloras-medium:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def build_checkpoint(
    candidates: pd.DataFrame,
    baseline_direct: pd.DataFrame,
    master: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    direct_pairs = set(zip(baseline_direct["accepted_species"], baseline_direct["trait_name"]))
    scope = classify_scope(master, load_config(Path("config/angiosperm_scope.yml")))
    eligible_species = set(
        scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"]
    )
    selection_rows: list[dict[str, str]] = []
    evidence_rows: list[dict[str, str]] = []
    for record in candidates.fillna("").to_dict("records"):
        candidate_id = _candidate_id(record)
        species = _text(record.get("accepted_species"))
        trait = _text(record.get("trait_name"))
        selected, reason = high_precision_candidate(record)
        if _text(record.get("name_match_method")) != "exact_accepted_name_match":
            selected, reason = False, "name_match_not_exact_accepted"
        elif _text(record.get("evidence_scope")) != "species_direct":
            selected, reason = False, "not_species_direct"
        elif species not in eligible_species:
            selected, reason = False, "outside_angiosperm_scope"
        elif (species, trait) in direct_pairs:
            selected, reason = False, "already_resolved_direct_species_trait"
        selection_rows.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": species,
                "trait_name": trait,
                "candidate_value": _text(record.get("candidate_value")),
                "selected": str(selected),
                "selection_reason": reason,
            }
        )
        if not selected:
            continue
        normalised = sorted(
            _normalised_candidate_values(trait, _text(record.get("candidate_value")))
        )
        flora_id = _text(record.get("flora_id"))
        taxon_id = _text(record.get("taxon_id"))
        evidence_rows.append(
            {
                "accepted_species": species,
                "axis": TRAIT_TO_AXIS[trait],
                "trait_name": trait,
                "normalized_value": "|".join(normalised),
                "quality": "medium",
                "source_group": "efloras",
                "source_provider": f"eFloras:{_text(record.get('flora_name'))}",
                "source_url": _text(record.get("source_url")),
                "source_record_id": candidate_id,
                "source_citation": (
                    f"{_text(record.get('flora_name'))}; treatment "
                    f"{flora_id}:{taxon_id}"
                ),
                "source_excerpt": _text(record.get("source_excerpt")),
                "evidence_scope": "species_direct",
                "name_match_method": "accepted_name_exact",
                "source_lineage": f"efloras-treatment:{flora_id}:{taxon_id}",
                "lineage_method": "provider_flora_taxon_treatment",
                "source_run_id": SOURCE_RUN_ID,
                "source_artifact": SOURCE_ARTIFACT,
                "source_file": "all_trait_candidates.csv",
                "acceptance_contract": "efloras_medium_reaudit_v1",
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS).sort_values(
        ["trait_name", "accepted_species", "source_record_id"]
    ).reset_index(drop=True)
    selection = pd.DataFrame(selection_rows).sort_values("candidate_id").reset_index(drop=True)
    audit_parts: list[pd.DataFrame] = []
    reviewed_urls: set[str] = set()
    for trait, count in AUDIT_ALLOCATION.items():
        group = evidence.loc[evidence["trait_name"].eq(trait)].copy()
        group["audit_hash"] = group["source_record_id"].map(
            lambda value: hashlib.sha256(
                f"efloras-medium-audit-v1|{value}".encode()
            ).hexdigest()
        )
        group = group.loc[~group["source_url"].isin(reviewed_urls)]
        chosen = group.sort_values("audit_hash").drop_duplicates("source_url").head(count)
        reviewed_urls.update(chosen["source_url"])
        audit_parts.append(chosen)
    sample = pd.concat(audit_parts, ignore_index=True).sort_values(
        ["trait_name", "audit_hash"]
    )
    audit = pd.DataFrame(
        {
            "candidate_id": sample["source_record_id"],
            "trait_name": sample["trait_name"],
            "accepted_species": sample["accepted_species"],
            "normalized_value": sample["normalized_value"],
            "source_url": sample["source_url"],
            "source_excerpt": sample["source_excerpt"],
            "accepted_correct": "",
            "cultivar_status": "",
            "reviewer": "",
            "reviewed_at_utc": "",
            "audit_reason": "",
            "audit_hash": sample["audit_hash"],
        }
    ).reset_index(drop=True)
    return evidence, audit, selection


def apply_reviewed_audit(template: pd.DataFrame, reviewed: pd.DataFrame) -> pd.DataFrame:
    """Attach review decisions without allowing sampled evidence to change."""

    immutable = [
        "candidate_id",
        "trait_name",
        "accepted_species",
        "normalized_value",
        "source_url",
        "source_excerpt",
        "audit_hash",
    ]
    decisions = [
        "accepted_correct",
        "cultivar_status",
        "reviewer",
        "reviewed_at_utc",
        "audit_reason",
    ]
    missing = set(immutable + decisions).difference(reviewed.columns)
    if missing:
        raise ValueError(f"reviewed audit is missing columns: {sorted(missing)}")
    if reviewed["candidate_id"].duplicated().any():
        raise ValueError("reviewed audit has duplicate candidate_id values")
    if set(template["candidate_id"]) != set(reviewed["candidate_id"]):
        raise ValueError("reviewed audit candidate IDs do not match deterministic sample")
    aligned = reviewed.set_index("candidate_id").loc[template["candidate_id"]].reset_index()
    for column in immutable:
        if not template[column].map(_text).eq(aligned[column].map(_text)).all():
            raise ValueError(f"reviewed audit changed immutable column: {column}")
    if aligned[decisions].apply(lambda column: column.map(_text).eq("").any()).any():
        raise ValueError("reviewed audit has incomplete decisions or provenance")
    accepted = aligned["accepted_correct"].map(_text).str.casefold()
    if not accepted.isin({"true", "false", "1", "0", "yes", "no"}).all():
        raise ValueError("reviewed audit has invalid accepted_correct values")
    result = template.copy()
    result[decisions] = aligned[decisions].to_numpy()
    return result


@app.command("build")
def build(
    candidates_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    baseline_direct_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    master_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    reviewed_audit_csv: Annotated[
        Path | None, typer.Option(exists=True, dir_okay=False)
    ] = None,
) -> None:
    candidates = pd.read_csv(candidates_csv, dtype=str).fillna("")
    baseline = pd.read_csv(baseline_direct_csv, dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    evidence, audit, selection = build_checkpoint(candidates, baseline, master)
    if len(audit) != sum(AUDIT_ALLOCATION.values()):
        raise ValueError(f"audit has {len(audit)} rows, expected 200")
    if audit["source_url"].duplicated().any():
        raise ValueError("audit must contain 200 distinct source pages")
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("evidence source_record_id values are not unique")
    reviewed_input_hash = ""
    if reviewed_audit_csv is not None:
        reviewed = pd.read_csv(reviewed_audit_csv, dtype=str).fillna("")
        audit = apply_reviewed_audit(audit, reviewed)
        reviewed_input_hash = _sha256(reviewed_audit_csv)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "efloras_medium_evidence_20260811.csv.gz"
    audit_path = output_dir / "efloras_medium_manual_audit_200_20260811.csv"
    selection_path = output_dir / "efloras_medium_selection_audit_20260811.csv.gz"
    evidence.to_csv(
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
    manifest = {
        "contract": "efloras_medium_checkpoint_v1",
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "source_artifact_id": SOURCE_ARTIFACT_ID,
        "source_run_conclusion": "failure_with_successful_final_collect_artifact",
        "input_hashes": {
            "candidates_csv_sha256": _sha256(candidates_csv),
            "baseline_direct_csv_sha256": _sha256(baseline_direct_csv),
            "master_csv_sha256": _sha256(master_csv),
            "reviewed_audit_csv_sha256": reviewed_input_hash,
        },
        "counts": {
            "candidate_rows": len(candidates),
            "selected_evidence_rows": len(evidence),
            "selected_species_trait": int(
                evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "audit_rows": len(audit),
            "selection_reasons": {
                str(key): int(value)
                for key, value in selection["selection_reason"].value_counts().items()
            },
            "selected_by_trait": {
                str(key): int(value)
                for key, value in evidence["trait_name"].value_counts().items()
            },
        },
        "output_hashes": {
            evidence_path.name: _sha256(evidence_path),
            audit_path.name: _sha256(audit_path),
            selection_path.name: _sha256(selection_path),
        },
        "manual_audit_required_before_promotion": reviewed_audit_csv is None,
        "manual_audit": {
            "reviewed": int(audit["reviewer"].map(_text).ne("").sum()),
            "accepted_correct": int(
                audit["accepted_correct"]
                .map(_text)
                .str.casefold()
                .isin({"true", "1", "yes"})
                .sum()
            ),
            "unique_source_pages": int(audit["source_url"].nunique()),
        },
        "family_inference": False,
        "global_fallback": False,
    }
    manifest_path = output_dir / "efloras_medium_checkpoint_manifest_20260811.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    typer.echo(json.dumps(manifest, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
