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

SYNONYM_REJECTED: dict[tuple[str, str], str] = {
    ("Pachygone ovata", "flower_primary_color"): (
        "red describes the ripe drupe; the treatment does not state a flower colour"
    ),
}

SYNONYM_OVERRIDES = {
    ("Archidendron bigeminum", "flower_primary_color"): "white",
    ("Mappianthus hookerianus", "flower_primary_color"): "white",
    ("Ziziphus xylopyrus", "flower_primary_color"): (
        "green_brown_inconspicuous|yellow_orange"
    ),
}

PROTA_REJECTED: dict[tuple[str, str], str] = {
    ("Croton sakamaliensis", "flower_primary_color"): (
        "yellowish describes ovary pubescence, not the primary floral display"
    ),
    ("Euphorbia antso", "flower_primary_color"): (
        "purple describes only the ovary and is not a primary flower-colour statement"
    ),
    ("Euphorbia candelabrum", "flower_primary_color"): (
        "golden-yellow describes the cyathial involucre, not the flowers"
    ),
    ("Euphorbia cooperi", "flower_primary_color"): (
        "golden-yellow describes the cyathial involucre, not the flowers"
    ),
    ("Euphorbia granulata", "flower_primary_color"): (
        "pink or white describes cyathial glands, not the flowers"
    ),
    ("Euphorbia schimperiana", "flower_primary_color"): (
        "green to brownish red describes cyathial appendages, not the flowers"
    ),
    ("Ficus bussei", "flower_primary_color"): (
        "green or yellow describes mature figs, not enclosed flower colour"
    ),
    ("Ficus glumosa", "flower_primary_color"): (
        "orange to red describes figs at the fruiting stage"
    ),
    ("Ficus natalensis", "flower_primary_color"): (
        "reddish-orange to brown describes mature figs"
    ),
    ("Ficus politoria", "flower_primary_color"): (
        "yellow to red-brown describes mature figs"
    ),
    ("Ficus sur", "flower_primary_color"): (
        "red to dark orange describes mature figs"
    ),
    ("Ficus tremula", "flower_primary_color"): (
        "greenish to brown describes mature figs"
    ),
    ("Ficus vogeliana", "flower_primary_color"): (
        "red, orange and pale spots describe mature figs"
    ),
    ("Milicia regia", "flower_primary_color"): (
        "white describes hairs on the catkin; individual flower colour is absent"
    ),
    ("Musanga cecropioides", "flower_primary_color"): (
        "greenish white explicitly describes the compound fruit"
    ),
    ("Uapaca guineensis", "flower_primary_color"): (
        "bright yellow describes involucral bracts enclosing the flower buds"
    ),
    ("Hyphaene thebaica", "floral_form"): (
        "tubular describes only the calyx base; overall floral form is not stated"
    ),
    ("Telfairia occidentalis", "floral_form"): (
        "campanulate describes the receptacle, not the flower or corolla form"
    ),
    ("Telfairia pedata", "floral_form"): (
        "campanulate describes the receptacle, not the flower or corolla form"
    ),
    ("Catharanthus trichophyllus", "autonomous_selfing_capacity"): (
        "the match occurs only in a cited reference title about variant periwinkle strains"
    ),
    ("Catharanthus trichophyllus", "self_incompatibility"): (
        "the self-incompatible state is explicitly restricted to Catharanthus roseus strains"
    ),
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


def curate_synonym_candidates(candidate_csv: Path, output_dir: Path) -> dict[str, object]:
    """Review IISc treatments resolved by exact WFO+GBIF synonym agreement."""

    candidates = pd.read_csv(candidate_csv, dtype=str).fillna("")
    required = {
        "accepted_species",
        "searched_name",
        "page_id",
        "family",
        "wfo_taxon_id",
        "wfo_accepted_usage_id",
        "gbif_usage_key",
        "gbif_accepted_usage_key",
        "axis",
        "trait_name",
        "normalized_value",
        "supporting_excerpt",
        "source_url",
        "page_content_sha256",
        "name_match_method",
    }
    missing = required - set(candidates.columns)
    if missing:
        raise ValueError(f"synonym candidate table missing columns: {sorted(missing)}")
    if len(candidates) != 71:
        raise ValueError(f"expected 71 synonym candidates, found {len(candidates)}")
    if not candidates["name_match_method"].eq(
        "strict_exact_two_backbone_wfo_gbif_synonym"
    ).all():
        raise ValueError("synonym candidates must have exact WFO+GBIF agreement")
    if not candidates["page_content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("synonym candidates contain invalid page hashes")

    reviewed: list[dict[str, str]] = []
    accepted: list[dict[str, str]] = []
    reviewed_at = "2026-08-12T14:29:00Z"
    for row in candidates.to_dict("records"):
        key = (row["accepted_species"], row["trait_name"])
        rejection = SYNONYM_REJECTED.get(key, "")
        value = SYNONYM_OVERRIDES.get(key, row["normalized_value"])
        candidate_id = _sha(
            "|".join(
                [
                    row["accepted_species"],
                    row["trait_name"],
                    value,
                    f"provider_treatment:iisc-india-flora-online:{row['page_id']}",
                ]
            )
        )[:24]
        decision = "reject" if rejection else "accept"
        reason = rejection or (
            "Accepted after exact species-rank synonym and family agreement in WFO "
            "June 2026 and GBIF, complete treatment quote, ontology and non-cultivar "
            "review; multistate floral values retained."
        )
        reviewed.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": row["accepted_species"],
                "searched_name": row["searched_name"],
                "trait_name": row["trait_name"],
                "normalized_value": value,
                "source_url": row["source_url"],
                "supporting_excerpt": row["supporting_excerpt"],
                "decision": decision,
                "species_identity_correct": "true",
                "value_correct": str(not rejection).lower(),
                "provenance_complete": "true",
                "cultivar_contamination": "false",
                "false_positive_reason": rejection,
                "decision_reason": reason,
                "reviewer": "Codex IISc two-backbone synonym full candidate audit",
                "reviewed_at_utc": reviewed_at,
            }
        )
        if rejection:
            continue
        accepted.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": row["accepted_species"],
                "searched_name": row["searched_name"],
                "page_id": row["page_id"],
                "family": row["family"],
                "wfo_taxon_id": row["wfo_taxon_id"],
                "wfo_accepted_usage_id": row["wfo_accepted_usage_id"],
                "gbif_usage_key": row["gbif_usage_key"],
                "gbif_accepted_usage_key": row["gbif_accepted_usage_key"],
                "axis": row["axis"],
                "trait_name": row["trait_name"],
                "normalized_value": value,
                "supporting_excerpt": row["supporting_excerpt"],
                "source_url": row["source_url"],
                "page_content_sha256": row["page_content_sha256"],
                "name_match_method": "exact_synonym",
                "cultivar_status": "wild_or_species_level_statement_not_cultivar_limited",
            }
        )

    evidence = pd.DataFrame(accepted).sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    )
    audit = pd.DataFrame(reviewed).sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    )
    if len(evidence) != 70 or len(audit) != 71:
        raise ValueError(
            f"unexpected synonym accepted/reviewed counts: {len(evidence)}/{len(audit)}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(
        output_dir / "india_flora_online_synonym_reviewed_records_20260812.csv",
        index=False,
        lineterminator="\n",
    )
    audit.to_csv(
        output_dir / "india_flora_online_synonym_full_candidate_audit_20260812.csv",
        index=False,
        lineterminator="\n",
    )
    return {
        "reviewed": len(audit),
        "accepted_correct": len(evidence),
        "precision": len(evidence) / len(audit),
        "cultivar_contamination_rate": 0.0,
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "accepted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
    }


def curate_prota_candidates(candidate_csv: Path, output_dir: Path) -> dict[str, object]:
    """Freeze the full source-backed audit of exact-name PROTA treatments."""

    candidates = pd.read_csv(candidate_csv, dtype=str).fillna("")
    required = {
        "accepted_species",
        "page_id",
        "axis",
        "trait_name",
        "normalized_value",
        "supporting_excerpt",
        "page_title",
        "source_url",
        "page_content_sha256",
        "revision_id",
        "revision_timestamp",
        "citation",
        "family",
    }
    missing = required - set(candidates.columns)
    if missing:
        raise ValueError(f"PROTA candidate table missing columns: {sorted(missing)}")
    if len(candidates) != 469:
        raise ValueError(f"expected 469 PROTA candidates, found {len(candidates)}")
    if candidates[["accepted_species", "trait_name", "supporting_excerpt"]].duplicated().any():
        raise ValueError("PROTA candidates are not unique at the reviewed source statement")
    if not candidates["page_content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("PROTA candidates contain invalid revision-content hashes")
    keys = set(zip(candidates["accepted_species"], candidates["trait_name"]))
    missing_decisions = sorted(set(PROTA_REJECTED) - keys)
    if missing_decisions:
        raise ValueError(f"PROTA review decisions absent from candidates: {missing_decisions}")

    reviewed: list[dict[str, str]] = []
    accepted: list[dict[str, str]] = []
    reviewed_at = "2026-08-12T14:35:00Z"
    for row in candidates.to_dict("records"):
        key = (row["accepted_species"], row["trait_name"])
        candidate_id = _sha(
            "|".join(
                [
                    "plantuse-prota",
                    row["revision_id"],
                    row["accepted_species"],
                    row["trait_name"],
                    row["supporting_excerpt"],
                ]
            )
        )[:24]
        rejection = PROTA_REJECTED.get(key)
        decision = "reject" if rejection else "accept"
        reason = rejection or (
            "Accepted after exact master-name and family identity, complete PROTA "
            "revision, species-specific description, trait ontology, non-cultivar "
            "and quoted-value review; multistate floral colours are retained."
        )
        reviewed.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": row["accepted_species"],
                "trait_name": row["trait_name"],
                "normalized_value": row["normalized_value"],
                "source_url": row["source_url"],
                "supporting_excerpt": row["supporting_excerpt"],
                "decision": decision,
                "species_identity_correct": "true",
                "value_correct": str(not rejection).lower(),
                "provenance_complete": "true",
                "cultivar_contamination": "false",
                "false_positive_reason": rejection or "",
                "decision_reason": reason,
                "reviewer": "Codex full PROTA monograph candidate audit",
                "reviewed_at_utc": reviewed_at,
            }
        )
        if rejection:
            continue
        accepted.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": row["accepted_species"],
                "family": row["family"],
                "page_id": row["page_id"],
                "revision_id": row["revision_id"],
                "revision_timestamp": row["revision_timestamp"],
                "axis": row["axis"],
                "trait_name": row["trait_name"],
                "normalized_value": row["normalized_value"],
                "supporting_excerpt": row["supporting_excerpt"],
                "page_title": row["page_title"],
                "source_url": row["source_url"],
                "citation": row["citation"],
                "page_content_sha256": row["page_content_sha256"],
                "cultivar_status": (
                    "wild_or_species_level_monograph_not_cultivar_limited"
                ),
            }
        )

    evidence = pd.DataFrame(accepted).sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    )
    audit = pd.DataFrame(reviewed).sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    )
    if len(evidence) != 448 or len(audit) != 469:
        raise ValueError(f"unexpected PROTA accepted/reviewed counts: {len(evidence)}/{len(audit)}")
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(
        output_dir / "plantuse_prota_reviewed_records_20260812.csv",
        index=False,
        lineterminator="\n",
    )
    audit.to_csv(
        output_dir / "plantuse_prota_full_candidate_audit_20260812.csv",
        index=False,
        lineterminator="\n",
    )
    return {
        "reviewed": len(audit),
        "accepted_correct": len(evidence),
        "precision": len(evidence) / len(audit),
        "cultivar_contamination_rate": 0.0,
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "accepted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "accepted_species_axis": int(
            evidence[["accepted_species", "axis"]].drop_duplicates().shape[0]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--synonym-candidate-csv", type=Path)
    parser.add_argument("--prota-candidate-csv", type=Path)
    args = parser.parse_args()
    result = {"accepted_name": curate(args.candidate_csv, args.output_dir)}
    if args.synonym_candidate_csv is not None:
        result["two_backbone_synonym"] = curate_synonym_candidates(
            args.synonym_candidate_csv, args.output_dir
        )
    if args.prota_candidate_csv is not None:
        result["prota_monograph"] = curate_prota_candidates(
            args.prota_candidate_csv, args.output_dir
        )
    print(result)


if __name__ == "__main__":
    main()
