"""Audit committed trait staging files for unintegrated strict reproductive evidence.

The public recovery baseline is Wave53-derived, while later acquisition work left
many reviewed files under ``data/v2/staging/traits``.  This audit walks committed
CSV/CSV.GZ/JSON records, compares them to the current cumulative direct ledger,
and reports only exact fixed-universe reproductive candidates whose species-axis
cell is still empty.

Nothing is promoted automatically.  Rows are split into:

* ``strict_ready_candidates``: High/Medium, ontology-valid, species/synonym-direct
  and provenance-complete, with an accepted review status;
* ``review_needed_candidates``: potentially useful rows that fail one or more
  strict admission checks.

This makes cached evidence recovery auditable without treating stale candidate,
probabilistic, inferred, machine-only or cross-trait rows as direct evidence.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml

AXIS = "reproductive_assurance"
STRICT_TRAITS = {
    "self_incompatibility",
    "autonomous_selfing_capacity",
    "mating_system",
    "cleistogamy",
}
DIRECT_SCOPES = {"species_direct", "synonym_direct"}
ACCEPTED_REVIEW_STATUSES = {
    "accepted_direct_statement",
    "accepted_correct",
    "accepted_source_backed_direct",
    "source_methodology_reviewed_reference_backed_strict_direct",
}
SKIP_PATH_TOKENS = {
    "/source_batches/",
}
SKIP_FILENAMES = {
    "restart_orchid_pending_20260907.csv",
}


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def boolish_false(value: object) -> bool:
    return text(value).casefold() in {"", "false", "0", "no", "none", "nan"}


def load_table(path: Path) -> pd.DataFrame | None:
    suffixes = [s.casefold() for s in path.suffixes]
    try:
        if suffixes[-2:] == [".csv", ".gz"]:
            return pd.read_csv(path, dtype=str).fillna("")
        if path.suffix.casefold() == ".csv":
            return pd.read_csv(path, dtype=str).fillna("")
        if path.suffix.casefold() == ".json":
            payload = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(payload, list):
                return pd.DataFrame(payload).fillna("")
            if isinstance(payload, dict):
                for key in ("records", "rows", "evidence", "data"):
                    value = payload.get(key)
                    if isinstance(value, list):
                        return pd.DataFrame(value).fillna("")
    except Exception:
        return None
    return None


def first(row: dict[str, object], *keys: str) -> str:
    for key in keys:
        if key in row:
            value = text(row.get(key))
            if value:
                return value
    return ""


def source_provenance(row: dict[str, object]) -> tuple[str, str, str]:
    url = first(row, "source_url", "url", "evidence_url", "source_record_url")
    lineage = first(row, "source_lineage", "provider_lineage", "lineage", "source_citation")
    excerpt = first(row, "source_excerpt", "excerpt", "evidence_quote", "quote", "source_raw_value")
    return url, lineage, excerpt


def iter_candidate_files(root: Path) -> Iterable[Path]:
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = "/" + str(path.as_posix()).lstrip("/")
        if any(token in rel for token in SKIP_PATH_TOKENS) or path.name in SKIP_FILENAMES:
            continue
        if path.suffix.casefold() == ".json" or path.suffix.casefold() == ".csv" or path.name.casefold().endswith(".csv.gz"):
            yield path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--staging-root", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--direct-ledger", type=Path, required=True)
    p.add_argument("--ontology", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    ontology = yaml.safe_load(args.ontology.read_text(encoding="utf-8"))
    allowed = {
        trait: set(ontology["traits"][trait]["allowed_values"]) - {"unresolved"}
        for trait in STRICT_TRAITS
    }

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106_295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected fixed reproductive universe, got {len(rep)}")
    universe = set(rep["accepted_species"])
    unresolved = set(rep.loc[rep["quality"].eq(""), "accepted_species"])

    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    current_pairs = set(zip(direct["accepted_species"], direct["trait_name"], strict=False))

    accepted: list[dict[str, object]] = []
    review: list[dict[str, object]] = []
    files_scanned = 0
    rows_scanned = 0
    reproductive_rows_seen = 0
    parseable_files = 0

    for path in iter_candidate_files(args.staging_root):
        files_scanned += 1
        frame = load_table(path)
        if frame is None or frame.empty:
            continue
        parseable_files += 1
        rows_scanned += len(frame)
        if "accepted_species" not in frame.columns or "trait_name" not in frame.columns:
            continue
        rel = str(path.as_posix())
        for raw in frame.to_dict("records"):
            trait = text(raw.get("trait_name"))
            if trait not in STRICT_TRAITS:
                continue
            reproductive_rows_seen += 1
            species = text(raw.get("accepted_species"))
            if species not in universe or species not in unresolved:
                continue
            if (species, trait) in current_pairs:
                continue

            value = first(raw, "normalized_value", "value")
            states_raw = raw.get("state_set", "")
            states: list[str] = []
            if isinstance(states_raw, list):
                states = [text(x) for x in states_raw if text(x)]
            elif text(states_raw):
                s = text(states_raw)
                try:
                    parsed = json.loads(s)
                    if isinstance(parsed, list):
                        states = [text(x) for x in parsed if text(x)]
                except Exception:
                    states = [x for x in s.split("|") if x]
            if not states and value:
                states = [x for x in value.split("|") if x]
            states = sorted(set(states))

            quality = first(raw, "quality", "selected_quality", "evidence_quality").casefold()
            scope = first(raw, "evidence_scope", "scope").casefold()
            review_status = first(raw, "review_status", "status")
            name_match = first(raw, "name_match_method", "match_method")
            url, lineage, excerpt = source_provenance(raw)
            genus_training = raw.get("genus_rule_training_allowed", False)

            checks = {
                "quality_high_medium": quality in {"high", "medium"},
                "ontology_valid": bool(states) and set(states).issubset(allowed[trait]),
                "direct_scope": scope in DIRECT_SCOPES,
                "accepted_review_status": review_status in ACCEPTED_REVIEW_STATUSES,
                "provenance_complete": bool(url and lineage and excerpt),
                "name_match_present": bool(name_match),
                "genus_training_not_required": boolish_false(genus_training),
            }
            out = {
                "accepted_species": species,
                "axis": AXIS,
                "trait_name": trait,
                "normalized_value": "|".join(states) if states else value,
                "quality": quality,
                "evidence_scope": scope,
                "review_status": review_status,
                "name_match_method": name_match,
                "source_url": url,
                "source_lineage": lineage,
                "source_excerpt": excerpt,
                "source_file": rel,
                "source_row_index": raw.get("source_row_index", raw.get("source_excel_row", "")),
                **checks,
                "strict_ready": all(checks.values()),
            }
            (accepted if out["strict_ready"] else review).append(out)

    ready = pd.DataFrame(accepted)
    needs = pd.DataFrame(review)
    if not ready.empty:
        ready = ready.sort_values(["source_file", "accepted_species", "trait_name"]).drop_duplicates(
            ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
        )
    if not needs.empty:
        needs = needs.sort_values(["source_file", "accepted_species", "trait_name"]).drop_duplicates(
            ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
        )
    ready.to_csv(args.output / "strict_ready_candidates.csv", index=False)
    needs.to_csv(args.output / "review_needed_candidates.csv", index=False)

    if len(ready):
        by_source = (
            ready.groupby("source_file", dropna=False)
            .agg(rows=("accepted_species", "size"), species=("accepted_species", "nunique"))
            .reset_index()
            .sort_values(["species", "rows"], ascending=False)
        )
    else:
        by_source = pd.DataFrame(columns=["source_file", "rows", "species"])
    by_source.to_csv(args.output / "strict_ready_by_source.csv", index=False)

    summary = {
        "contract": "committed_reproductive_residual_audit_v1",
        "current_reproductive_unresolved": len(unresolved),
        "files_scanned": files_scanned,
        "parseable_files": parseable_files,
        "rows_scanned": rows_scanned,
        "reproductive_rows_seen": reproductive_rows_seen,
        "strict_ready_rows": int(len(ready)),
        "strict_ready_species": int(ready["accepted_species"].nunique()) if len(ready) else 0,
        "strict_ready_species_trait_pairs": int(ready[["accepted_species", "trait_name"]].drop_duplicates().shape[0]) if len(ready) else 0,
        "review_needed_rows": int(len(needs)),
        "review_needed_species": int(needs["accepted_species"].nunique()) if len(needs) else 0,
        "strict_ready_source_files": int(ready["source_file"].nunique()) if len(ready) else 0,
        "formal_gain": 0,
        "promotion_allowed": False,
        "claim_limit": "Audit only. Strict-ready rows still require collision/source-lineage review before any batch integration.",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
