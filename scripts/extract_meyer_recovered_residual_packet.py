from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

EXPECTED = {"Eulobus californicus", "Syzygium nervosum"}
SOURCE_URL = "https://doi.org/10.5061/dryad.cc2fqz6hr"
SOURCE_LINEAGE = "dataset:dryad.cc2fqz6hr:v4"
REVIEW_STATUS = "source_methodology_reviewed_reference_backed_strict_direct"


def text(v: object) -> str:
    if v is None or pd.isna(v):
        return ""
    return " ".join(str(v).strip().split())


def truth(v: object) -> bool:
    return text(v).casefold() in {"true", "1", "yes"}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--manual-audit", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--direct-ledger", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    evidence = pd.read_csv(args.evidence, dtype=str).fillna("")
    audit = pd.read_csv(args.manual_audit, dtype=str).fillna("")
    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")

    rep = coverage.loc[coverage["axis"].eq("reproductive_assurance")]
    if len(rep) != 106_295 or rep["accepted_species"].nunique() != 106_295:
        raise ValueError("current checkpoint is not the fixed 106,295-species reproductive universe")
    unresolved = set(rep.loc[rep["quality"].eq(""), "accepted_species"])
    direct_pairs = set(zip(direct["accepted_species"], direct["trait_name"], strict=True))

    approved_audit = audit.loc[
        audit["decision"].str.casefold().eq("accept")
        & audit["species_identity_correct"].map(truth)
        & audit["value_correct"].map(truth)
        & audit["provenance_complete"].map(truth)
        & ~audit["cultivar_contamination"].map(truth)
    ].copy()
    approved_ids = set(approved_audit["candidate_id"])

    src = evidence.loc[
        evidence["candidate_id"].isin(approved_ids)
        & evidence["source_provider"].eq("Meyer, Galloway & Eckert 2026 Dryad v4")
        & evidence["source_url"].eq(SOURCE_URL)
        & evidence["source_lineage"].eq(SOURCE_LINEAGE)
        & evidence["trait_name"].eq("self_incompatibility")
        & evidence["normalized_value"].isin({"SI", "SC"})
        & evidence["evidence_quality"].str.casefold().isin({"high", "medium"})
        & evidence["evidence_scope"].eq("species_direct")
        & evidence["name_match_method"].eq("accepted_name_exact")
        & evidence["evidence_status"].eq("accepted_individual_structured_source_audit")
    ].copy()
    src = src.loc[src["accepted_species"].isin(unresolved)].copy()
    src = src.loc[
        [
            (species, "self_incompatibility") not in direct_pairs
            for species in src["accepted_species"]
        ]
    ].copy()
    src = src.sort_values("accepted_species", kind="stable").drop_duplicates(
        ["accepted_species", "trait_name"], keep="first"
    )

    observed = set(src["accepted_species"])
    if observed != EXPECTED or len(src) != 2:
        raise ValueError(f"expected exactly {sorted(EXPECTED)}, observed {sorted(observed)}")

    out = pd.DataFrame(
        {
            "accepted_species": src["accepted_species"],
            "axis": "reproductive_assurance",
            "trait_name": "self_incompatibility",
            "normalized_value": src["normalized_value"],
            "quality": src["evidence_quality"].str.casefold(),
            "source_url": src["source_url"],
            "source_lineage": src["source_lineage"],
            "source_reference_raw": src.apply(
                lambda r: f"{text(r['source_citation'])}; {text(r['source_record_id'])}", axis=1
            ),
            "source_column": "mating_system",
            "source_raw_value": src["raw_value"],
            "name_match_method": src["name_match_method"],
            "review_status": REVIEW_STATUS,
            "promotion_allowed": "false",
            "genus_rule_training_allowed": "false",
            "source_excel_row": src["source_record_id"],
            "source_article_url": SOURCE_URL,
        }
    )
    out.to_csv(
        args.output / "meyer_recovered_residual_batch.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    audit_subset = approved_audit.loc[approved_audit["candidate_id"].isin(set(src["candidate_id"]))].copy()
    audit_subset.to_csv(args.output / "meyer_recovered_residual_manual_audit.csv", index=False)

    summary = {
        "contract": "meyer_recovered_public_residual_packet_v1",
        "source": "Meyer, Galloway & Eckert 2026 Dryad v4",
        "source_url": SOURCE_URL,
        "source_lineage": SOURCE_LINEAGE,
        "source_mapping": "mating_system=sc/si -> self_incompatibility=SC/SI only",
        "reviewed_current_gap_rows": int(len(out)),
        "reviewed_current_gap_species": sorted(out["accepted_species"].tolist()),
        "state_counts": {str(k): int(v) for k, v in out["normalized_value"].value_counts().items()},
        "manual_audit_required_fields": {
            "decision": "accept",
            "species_identity_correct": True,
            "value_correct": True,
            "provenance_complete": True,
            "cultivar_contamination": False,
        },
        "formal_gain": 0,
        "promotion_allowed": False,
        "genus_rule_training_allowed": False,
        "integration_deferred_to_source_batch": True,
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
