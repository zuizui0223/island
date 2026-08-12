"""Freeze the fully reviewed IISc India Flora Online bulk morphology acquisition.

The input is the deterministic candidate table extracted from the downloaded
IISc species treatments.  Every candidate is reviewed here: false transfers
from fruits, stems, bracts or spathes are rejected, horticultural hybrid
identity is rejected, and legitimate multistate flower colours are retained.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd

REJECTED: dict[tuple[str, str], tuple[str, bool]] = {
    ("Blumea densiflora", "flower_primary_color"): (
        "red describes the pappus, not the corolla or another flower-colour organ",
        False,
    ),
    ("Cordia grandis", "flower_primary_color"): (
        "white describes the drupe after a semicolon; the flower colour is absent",
        False,
    ),
    ("Croton tiglium", "flower_primary_color"): (
        "whitish and greenish brown describe fruits; flower colour is absent",
        False,
    ),
    ("Dinochloa andamanica", "flower_primary_color"): (
        "dark green describes internodes; flower colour is absent",
        False,
    ),
    ("Euphorbia griffithii", "flower_primary_color"): (
        "red is stated only for the involucre and is not mapped to strict flower colour",
        False,
    ),
    ("Milicia excelsa", "flower_primary_color"): (
        "white and green describe catkins and spikes without an explicit flower-colour statement",
        False,
    ),
    ("Nicotiana sanderae", "flower_primary_color"): (
        "the treatment is horticultural hybrid material and fails the wild-species identity gate",
        True,
    ),
    ("Quercus semecarpifolia", "flower_primary_color"): (
        "brown describes nuts; flower colour is absent",
        False,
    ),
    ("Quercus semiserrata", "flower_primary_color"): (
        "brown describes nuts; flower colour is absent",
        False,
    ),
    ("Spiraea bumalda", "flower_primary_color"): (
        "the landscaping treatment is horticultural hybrid material and fails the identity gate",
        True,
    ),
    ("Swintonia floribunda", "flower_primary_color"): (
        "reddish describes drupes; flower colour is absent",
        False,
    ),
    ("Licuala peltata", "floral_form"): (
        "tubular describes spathes rather than the flowers or corolla",
        False,
    ),
}


OVERRIDES = {
    ("Castanopsis indica", "flower_primary_color"): "white",
    ("Cleisostoma discolor", "flower_primary_color"): (
        "green_brown_inconspicuous|white|yellow_orange"
    ),
    ("Clerodendrum schmidtii", "flower_primary_color"): "white",
    ("Clerodendrum urticifolium", "flower_primary_color"): "red_pink",
    ("Crotalaria hirta", "flower_primary_color"): "yellow_orange",
    ("Croton joufra", "flower_primary_color"): "yellow_orange",
    ("Derris scandens", "flower_primary_color"): "red_pink|white",
    ("Endiandra firma", "flower_primary_color"): "yellow_orange",
    ("Eurya acuminata", "flower_primary_color"): "white",
    ("Grewia eriocarpa", "flower_primary_color"): "yellow_orange",
    ("Ilex odorata", "flower_primary_color"): "white",
    ("Miliusa roxburghiana", "flower_primary_color"): "red_pink|white",
    ("Mussaenda philippica", "flower_primary_color"): (
        "green_brown_inconspicuous|white|yellow_orange"
    ),
    ("Rhopalocnemis phalloides", "flower_primary_color"): "white",
    ("Sterculia hamiltonii", "flower_primary_color"): "yellow_orange",
    ("Syzygium laetum", "flower_primary_color"): "red_pink|white",
    ("Trema orientale", "flower_primary_color"): "green_brown_inconspicuous",
}


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def curate(candidate_csv: Path, output_dir: Path) -> dict[str, object]:
    candidates = pd.read_csv(candidate_csv, dtype=str).fillna("")
    required = {
        "accepted_species",
        "page_id",
        "axis",
        "trait_name",
        "normalized_value",
        "supporting_excerpt",
        "comments",
        "source_url",
        "page_content_sha256",
    }
    missing = required - set(candidates.columns)
    if missing:
        raise ValueError(f"candidate table missing columns: {sorted(missing)}")
    if len(candidates) != 262:
        raise ValueError(f"expected 262 extracted candidates, found {len(candidates)}")
    keys = set(zip(candidates["accepted_species"], candidates["trait_name"]))
    missing_decisions = sorted((set(REJECTED) | set(OVERRIDES)) - keys)
    if missing_decisions:
        raise ValueError(f"review decisions absent from candidates: {missing_decisions}")

    reviewed: list[dict[str, str]] = []
    accepted: list[dict[str, str]] = []
    reviewed_at = "2026-08-12T13:18:00Z"
    for row in candidates.to_dict("records"):
        key = (row["accepted_species"], row["trait_name"])
        candidate_id = _sha(
            "|".join(
                [
                    "iisc-india-flora-online",
                    row["page_id"],
                    row["accepted_species"],
                    row["trait_name"],
                    row["supporting_excerpt"],
                ]
            )
        )[:24]
        rejection = REJECTED.get(key)
        decision = "reject" if rejection else "accept"
        value = OVERRIDES.get(key, row["normalized_value"])
        cultivar = bool(rejection and rejection[1])
        reason = rejection[0] if rejection else (
            "Accepted after exact master-name identity, complete species-page quote, "
            "trait-ontology and non-cultivar review; multistate values retained."
        )
        reviewed.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": row["accepted_species"],
                "trait_name": row["trait_name"],
                "normalized_value": value,
                "source_url": row["source_url"],
                "supporting_excerpt": row["supporting_excerpt"],
                "decision": decision,
                "species_identity_correct": str(not cultivar).lower(),
                "value_correct": str(not rejection).lower(),
                "provenance_complete": "true",
                "cultivar_contamination": str(cultivar).lower(),
                "false_positive_reason": rejection[0] if rejection else "",
                "decision_reason": reason,
                "reviewer": "Codex full IISc India Flora Online candidate audit",
                "reviewed_at_utc": reviewed_at,
            }
        )
        if rejection:
            continue
        cultivated = bool(
            pd.Series([row["comments"]]).str.contains(
                r"cultivat|ornamental|planted|garden|plantation",
                case=False,
                regex=True,
            ).iloc[0]
        )
        accepted.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": row["accepted_species"],
                "page_id": row["page_id"],
                "axis": row["axis"],
                "trait_name": row["trait_name"],
                "normalized_value": value,
                "supporting_excerpt": row["supporting_excerpt"],
                "source_url": row["source_url"],
                "page_content_sha256": row["page_content_sha256"],
                "cultivar_status": (
                    "species_level_cultivated_record_not_cultivar_limited"
                    if cultivated
                    else "wild_or_species_level_statement_not_cultivar_limited"
                ),
            }
        )

    evidence = pd.DataFrame(accepted).sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    )
    audit = pd.DataFrame(reviewed).sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    )
    if len(evidence) != 250 or len(audit) != 262:
        raise ValueError(f"unexpected accepted/reviewed counts: {len(evidence)}/{len(audit)}")
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(
        output_dir / "india_flora_online_reviewed_records_20260812.csv",
        index=False,
        lineterminator="\n",
    )
    audit.to_csv(
        output_dir / "india_flora_online_full_candidate_audit_20260812.csv",
        index=False,
        lineterminator="\n",
    )
    return {
        "reviewed": len(audit),
        "accepted_correct": len(evidence),
        "precision": len(evidence) / len(audit),
        "cultivar_contamination_rate": int(
            audit["cultivar_contamination"].eq("true").sum()
        )
        / len(audit),
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "accepted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(curate(args.candidate_csv, args.output_dir))


if __name__ == "__main__":
    main()
