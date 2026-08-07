"""Review and integrate the Europe PMC reproduction full-text acquisition.

All 43 machine candidates are explicitly classified.  Search results and
snippets are never evidence; promoted rows retain the fetched full-text quote,
article DOI/PMCID, content hash, and review decision.  The direct ledger and
trait-specific genus rules are rebuilt from all available direct evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.all_evidence_trait_audit import (
    apply_genus_rules,
    build_rule_audit,
    coverage_snapshot,
    dedupe_direct_lineages,
    direct_evidence_from_integrated,
    load_ontology,
    resolve_direct_cells,
    species_axis_coverage,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

REVIEWER = "codex_fulltext_context_review_20260807"
RUN_ID = "europe-pmc-reproduction-source-scale-20260807"
REVIEWED_AT = "2026-08-07T10:00:00Z"
SOURCE_ARTIFACT = "20260807_europe_pmc_reproduction_scale"

# One source-backed row is promoted for each accepted species x trait claim.
# Repeated papers restating the same Panicum virgatum claim are retained in
# the audit as corroboration but are not treated as independent origins.
PROMOTIONS = {
    "EPMC-a2dbcebe88523839e7e3": (
        "high",
        "Paper explicitly links Calanthe mannii to an autonomous self-pollination mechanism.",
    ),
    "EPMC-27a1622171ec27207648": (
        "high",
        "Experimental paper explicitly reports Lobelia laxiflora as completely self-incompatible.",
    ),
    "EPMC-82d887bd8eb8919cfb3f": (
        "high",
        "Proceedings paper explicitly reports Panicum hallii as self-compatible.",
    ),
    "EPMC-e9765583fe4e9f7dc7a5": (
        "high",
        "Peer-reviewed paper explicitly reports switchgrass as self-incompatible and outcrossing.",
    ),
    "EPMC-0d0a3db5edba905fb086": (
        "medium",
        (
            "Taxonomic treatment attributes self-compatibility to Poa leptoclada "
            "but qualifies the statement as presumed."
        ),
    ),
    "EPMC-4a3490eb81ac9cc6be4a": (
        "high",
        "Peer-reviewed paper explicitly reports sexual Poa pratensis as self-incompatible.",
    ),
    "EPMC-9028ab2d92d829009995": (
        "medium",
        (
            "Taxonomic treatment attributes self-incompatibility to Poa trivialis "
            "but qualifies the statement as reputed."
        ),
    ),
    "EPMC-b24b20878ac988f445b6": (
        "high",
        (
            "Peer-reviewed species study explicitly reports self-incompatible "
            "Rumex bucephalophorus populations."
        ),
    ),
}

CORROBORATING_NOT_COUNTED = {
    "EPMC-4e19f07990aee3f084e4",
    "EPMC-9fbbefd374d8d1c017f6",
    "EPMC-7da953752eb80fddd382",
    "EPMC-7cee6a89c64bc94e4aee",
    "EPMC-b34d6f827eb4707798d7",
    "EPMC-d6268f18fa690500c6d2",
    "EPMC-b343d7d9d6fba7eeb870",
    "EPMC-61e195d558a2ead5896e",
    "EPMC-15fbd0cc8d6ae1d7a97a",
}

REJECTIONS = {
    "EPMC-b74dcd3e07b12acf63b6": ("trait_statement_applies_to_other_listed_species"),
    "EPMC-0bb9447e7946bcba8f14": ("self_incompatibility_statement_applies_to_other_species"),
    "EPMC-2568bc08c406bc68a0eb": ("self_compatibility_statement_applies_to_Oryza_not_target"),
    "EPMC-dba25d8c38c695188c84": ("self_compatibility_statement_applies_to_Setaria_not_target"),
    "EPMC-29f791e59b159ec19f46": ("self_compatibility_statement_applies_to_Setaria_not_target"),
    "EPMC-67b39a00d6b71fd658e5": ("table_column_bleed_trait_applies_to_neighboring_taxon"),
    "EPMC-b84ea8a7ad7f8b41b2dc": ("pronoun_attribution_applies_to_parent_species_not_hybrid"),
    "EPMC-07914e008fd9162f4606": ("genus_statement_is_for_Rheum_not_target_Rumex_species"),
}

BATCH2_PROMOTIONS = {
    "EPMC-9cd833846c73cbf651e7": (
        "high",
        "Peer-reviewed paper explicitly reports Bidens frondosa among the listed mixed-mating Bidens species.",
    ),
    "EPMC-c0c4a8cf75670e3081f2": (
        "high",
        "Peer-reviewed species paper explicitly reports Brachypodium distachyon as self-fertile and autogamous.",
    ),
    "EPMC-08f862b8ecf36f5bcb59": (
        "high",
        "Peer-reviewed species paper explicitly reports Brachypodium sylvaticum as self-fertile.",
    ),
    "EPMC-ad1f6aaf56946d10e151": (
        "high",
        "Peer-reviewed field study explicitly reports Fagus crenata flowers as self-incompatible.",
    ),
    "EPMC-73a751730919de2a1c8d": (
        "high",
        "Peer-reviewed paper explicitly reports Noccaea caerulescens as self-compatible.",
    ),
    "EPMC-d6bf40eb207c4d59ddd8": (
        "high",
        "Peer-reviewed species paper explicitly reports Psittacanthus schiedeanus as self-compatible.",
    ),
    "EPMC-e95615f6842b8012ddfe": (
        "high",
        "Peer-reviewed species paper explicitly reports Sagittaria guayanensis as self-compatible.",
    ),
    "EPMC-738892a1d11a24cd6563": (
        "high",
        "Peer-reviewed species study explicitly reports Sagittaria pygmaea as self-compatible.",
    ),
    "EPMC-3aa25edc152ca66e45bc": (
        "high",
        "Peer-reviewed species paper explicitly reports Sagittaria trifolia as self-compatible.",
    ),
    "EPMC-0d56a3d32adbdc8c214e": (
        "high",
        "Peer-reviewed genomic paper explicitly reports Vanilla planifolia as self-compatible.",
    ),
}

BATCH2_CORRECTED_PROMOTIONS = {
    "EPMC-2539e0a1244a4a8a6d51": (
        "high",
        "mating_system",
        "mixed_mating",
        (
            "The machine candidate attached self-incompatibility from Solidago "
            "canadensis to Ambrosia artemisiifolia. The same fetched paragraph "
            "instead explicitly reports mixed mating for A. artemisiifolia."
        ),
    ),
    "EPMC-8aad4347ae2dd0dbfde1": (
        "high",
        "mating_system",
        "predominantly_selfing",
        (
            "The machine candidate attributed Noccaea corymbosa "
            "self-incompatibility to N. brachypetala. The fetched paragraph "
            "instead explicitly includes N. brachypetala among the "
            "predominantly selfing species."
        ),
    ),
}

BATCH2_CORROBORATING_NOT_COUNTED = {
    "EPMC-69bf09ca8cda4bb44c98",
    "EPMC-bea24dd2e5cec425427c",
    "EPMC-8887ed3eb288ce104533",
    "EPMC-f562e1d55c40402968bf",
    "EPMC-7812e756bcaed8008434",
    "EPMC-a14fd0b820a5edd5cc91",
    "EPMC-756e1538427d2b04d660",
    "EPMC-cce2d17584c6d8f0a963",
    "EPMC-dfad9d98fc77aa4cb376",
    "EPMC-26dc275d8a7c273be3f4",
    "EPMC-daaa09f7896278de77d7",
    "EPMC-b4b30fa933428cddba9f",
    "EPMC-8d99a7fd938abf4cdb3e",
    "EPMC-661ac91215003a92f571",
    "EPMC-111802794ed705c63862",
    "EPMC-36770ca59f2d7d3b0dc3",
    "EPMC-e1140c629623db12f8a3",
    "EPMC-a2cd9d666b831dfe7faf",
    "EPMC-c50b8790287330a1b4f3",
    "EPMC-dc6f7b984db9d0bbf8f3",
    "EPMC-b577e976a7afe96f3841",
    "EPMC-e49e96e527fffff3998b",
    "EPMC-bd9021210f81c596f820",
    "EPMC-589075ca005d70efb98a",
    "EPMC-d7c5bf1df80b4388d50e",
    "EPMC-924f727890c6503d7f42",
    "EPMC-c8804168b864513a845d",
}

BATCH2_REJECTIONS = {
    "EPMC-d19d2dce45d4e38e47b8": (
        "self_incompatibility_statement_applies_to_Astrocaryum_mexicanum"
    ),
    "EPMC-d5950ffa40dfc1db0fbd": (
        "self_incompatibility_statement_applies_to_perennial_ryegrass"
    ),
    "EPMC-9e9a77e0800c050793bd": (
        "generic_mating_system_list_does_not_classify_Brachypodium_distachyon"
    ),
    "EPMC-fefde154b8171bf70a90": (
        "generic_mating_system_list_does_not_classify_Brachypodium_distachyon"
    ),
    "EPMC-fb98e52815dd8512235e": (
        "self_compatibility_statement_applies_to_Arabidopsis"
    ),
    "EPMC-2a81fa9a4ad44ad52b86": (
        "self_compatibility_statement_applies_to_Setaria"
    ),
    "EPMC-9e2bd44b0fb9218cf0f9": (
        "self_compatibility_statement_applies_to_Hosta_not_Fagus"
    ),
    "EPMC-dcad8d86ec48ed756350": (
        "commodity_table_column_bleed_from_self_fertile_Malus_cultivar"
    ),
    "EPMC-178739613942584d69b1": (
        "self_compatible_phrase_applies_to_tomato_not_Leymus"
    ),
    "EPMC-b3e73e8cde8cf5b177d8": (
        "self_compatibility_statement_applies_to_Noccaea_ochroleucum"
    ),
    "EPMC-b045442acadcc281cbb9": (
        "self_incompatibility_statement_applies_to_Arabidopsis_halleri"
    ),
    "EPMC-2d96e719410df5b4a42e": (
        "self_incompatibility_statement_applies_to_other_island_endemics"
    ),
    "EPMC-783d0a2b99ced995588c": (
        "generic_self_compatible_phrase_does_not_describe_Sagittaria_lancifolia"
    ),
    "EPMC-447a7671622dede58e19": (
        "genus_level_pollination_statement_is_explicitly_uncertain_for_Vanilla_planifolia"
    ),
}

INVALIDATED_PRIOR = {
    (
        "Panicum hallii",
        "self_incompatibility",
        "SI",
        "url:https://europepmc.org/article/MED/28449656",
    ): (
        "Excerpt says that switchgrass (Panicum virgatum), unlike P. hallii, "
        "is self-incompatible; the old extraction attached the comparison "
        "trait to the wrong species."
    ),
    (
        "Brachypodium distachyon",
        "self_incompatibility",
        "SI",
        "url:https://europepmc.org/article/MED/19943086",
    ): (
        "The excerpt describes self-incompatibility loci in perennial ryegrass, "
        "not Brachypodium distachyon."
    ),
    (
        "Brachypodium distachyon",
        "self_incompatibility",
        "SI",
        "url:https://europepmc.org/article/MED/26408275",
    ): (
        "The excerpt describes the grass self-incompatibility mechanism without "
        "attributing it to Brachypodium distachyon."
    ),
    (
        "Brachypodium distachyon",
        "self_incompatibility",
        "SI",
        "url:https://europepmc.org/article/MED/28824669",
    ): (
        "The excerpt describes a perennial ryegrass population, not "
        "Brachypodium distachyon."
    ),
    (
        "Brachypodium distachyon",
        "self_incompatibility",
        "SI",
        "url:https://europepmc.org/article/MED/29247258",
    ): (
        "The excerpt describes pollen self-compatibility and SI loci in "
        "perennial ryegrass and grasses generally, not Brachypodium distachyon."
    ),
}
BATCH2_INVALIDATED_PRIOR = INVALIDATED_PRIOR
INVALIDATED_PRIOR = {
    key: reason
    for key, reason in BATCH2_INVALIDATED_PRIOR.items()
    if key[0] == "Panicum hallii"
}


def text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def keyed(frame: pd.DataFrame, columns: list[str]) -> set[tuple[str, ...]]:
    return set(frame[columns].fillna("").itertuples(index=False, name=None))


def review_candidates(
    candidates: pd.DataFrame,
    review_batch: str = "20260807",
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    unclassified: list[str] = []
    for row in candidates.to_dict("records"):
        candidate_id = text(row["candidate_id"])
        override_trait = ""
        override_value = ""
        if review_batch == "20260808_batch2" and candidate_id in BATCH2_PROMOTIONS:
            quality, reason = BATCH2_PROMOTIONS[candidate_id]
            correct = True
            promote = True
            decision = "accepted"
        elif (
            review_batch == "20260808_batch2"
            and candidate_id in BATCH2_CORRECTED_PROMOTIONS
        ):
            quality, override_trait, override_value, reason = (
                BATCH2_CORRECTED_PROMOTIONS[candidate_id]
            )
            correct = False
            promote = True
            decision = "corrected_and_accepted"
        elif (
            review_batch == "20260808_batch2"
            and candidate_id in BATCH2_CORROBORATING_NOT_COUNTED
        ):
            quality = "high"
            correct = True
            promote = False
            decision = "corroborating_not_counted"
            reason = (
                "The species-level statement is correct, but the current direct "
                "ledger already contains the same state at High quality; this "
                "restatement is retained for audit without adding an independent "
                "origin."
            )
        elif review_batch == "20260808_batch2" and candidate_id in BATCH2_REJECTIONS:
            quality = ""
            correct = False
            promote = False
            decision = "rejected"
            reason = BATCH2_REJECTIONS[candidate_id]
        elif review_batch == "20260808_batch2":
            unclassified.append(candidate_id)
            continue
        elif candidate_id in PROMOTIONS:
            quality, reason = PROMOTIONS[candidate_id]
            correct = True
            promote = True
            decision = "accepted"
        elif candidate_id in CORROBORATING_NOT_COUNTED:
            quality = "high"
            correct = True
            promote = False
            decision = "corroborating_not_counted"
            reason = (
                "The species-level statement is correct, but repeated "
                "Panicum virgatum restatements are not counted as independent "
                "source origins."
            )
        elif text(row["pmcid"]) == "PMC10784871":
            quality = ""
            correct = False
            promote = False
            decision = "rejected"
            reason = (
                "A commodity-risk table contains the target name and an "
                "unrelated self-fertile cultivar phrase far away in the same "
                "XML table paragraph."
            )
        elif candidate_id in REJECTIONS:
            quality = ""
            correct = False
            promote = False
            decision = "rejected"
            reason = REJECTIONS[candidate_id]
        else:
            unclassified.append(candidate_id)
            continue
        updated = dict(row)
        if override_trait:
            updated.update(
                {
                    "extracted_trait_name": text(row["trait_name"]),
                    "extracted_normalized_value": text(row["normalized_value"]),
                    "trait_name": override_trait,
                    "normalized_value": override_value,
                    "raw_value": override_value,
                    "axis": "reproductive_assurance",
                }
            )
        updated.update(
            {
                "accepted_correct": correct,
                "promotion_accepted": promote,
                "promotion_quality": quality,
                "review_decision": decision,
                "review_reason": reason,
                "reviewer": (
                    "codex_fulltext_context_review_20260808"
                    if review_batch == "20260808_batch2"
                    else REVIEWER
                ),
                "reviewed_at_utc": (
                    "2026-08-07T15:40:00Z"
                    if review_batch == "20260808_batch2"
                    else REVIEWED_AT
                ),
                "human_review_claimed": False,
                "audit_denominator_contract": "accepted_correct / all_reviewed",
                "review_batch": review_batch,
            }
        )
        rows.append(updated)
    if unclassified:
        raise ValueError(f"unclassified candidates: {sorted(unclassified)}")
    review = pd.DataFrame(rows)
    if len(review) != len(candidates):
        raise ValueError("review does not cover every candidate exactly once")
    return review.sort_values(
        ["accepted_species", "trait_name", "pmcid", "candidate_id"]
    ).reset_index(drop=True)


def common_evidence(
    review: pd.DataFrame,
    *,
    run_id: str = RUN_ID,
    source_artifact: str = SOURCE_ARTIFACT,
) -> pd.DataFrame:
    rows = []
    for row in review.loc[review["promotion_accepted"]].to_dict("records"):
        doi = text(row["doi"]).lower()
        pmcid = text(row["pmcid"])
        rows.append(
            {
                "accepted_species": text(row["accepted_species"]),
                "axis": "reproductive_assurance",
                "trait_name": text(row["trait_name"]),
                "normalized_value": text(row["normalized_value"]),
                "quality": text(row["promotion_quality"]),
                "source_group": "europe_pmc_fulltext",
                "source_provider": "Europe PMC",
                "source_url": text(row["source_url"]),
                "source_record_id": text(row["candidate_id"]),
                "source_citation": text(row["article_title"]),
                "source_excerpt": text(row["exact_supporting_quote"]),
                "evidence_scope": "species_direct",
                "name_match_method": text(row["name_match_method"]),
                "source_lineage": f"doi:{doi}" if doi else f"pmcid:{pmcid}",
                "lineage_method": "article_doi" if doi else "pmcid",
                "source_run_id": run_id,
                "source_artifact": source_artifact,
                "source_file": "fulltext_reproductive_candidates.csv.gz",
                "acceptance_contract": "fulltext_exact_species_context_review_v1",
            }
        )
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).fillna("")


def master_species(
    coverage: pd.DataFrame,
    low: pd.DataFrame,
    genus_path: Path,
) -> pd.DataFrame:
    formal = pd.read_csv(
        genus_path,
        usecols=["accepted_species", "genus"],
        dtype=str,
    ).fillna("")
    formal = formal.loc[formal["genus"].ne("")].drop_duplicates("accepted_species", keep="first")
    low_genus = (
        low.loc[low["genus"].ne(""), ["accepted_species", "genus"]]
        .groupby(["accepted_species", "genus"], as_index=False)
        .size()
        .sort_values(
            ["accepted_species", "size", "genus"],
            ascending=[True, False, True],
        )
        .drop_duplicates("accepted_species", keep="first")
        .drop(columns="size")
    )
    master = coverage[["accepted_species"]].drop_duplicates()
    master = master.merge(
        formal.rename(columns={"genus": "formal_genus"}),
        on="accepted_species",
        how="left",
    ).merge(
        low_genus.rename(columns={"genus": "low_genus"}),
        on="accepted_species",
        how="left",
    )
    master["genus"] = (
        master["formal_genus"]
        .fillna("")
        .where(
            master["formal_genus"].fillna("").ne(""),
            master["low_genus"].fillna(""),
        )
    )
    master["genus"] = master["genus"].where(
        master["genus"].ne(""),
        master["accepted_species"].str.split().str[0],
    )
    return master[["accepted_species", "genus"]]


def write_manifest(output: Path) -> None:
    rows = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "file_manifest.sha256":
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            rows.append(f"{digest}  {path.relative_to(output).as_posix()}")
    (output / "file_manifest.sha256").write_text(
        "\n".join(rows) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--formal-lineage", type=Path, required=True)
    parser.add_argument("--prior-promoted-common", type=Path, required=True)
    parser.add_argument("--wfo-common", type=Path, required=True)
    parser.add_argument(
        "--additional-common",
        type=Path,
        action="append",
        default=[],
        help="Additional previously promoted direct-evidence files to preserve.",
    )
    parser.add_argument("--current-direct", type=Path, required=True)
    parser.add_argument("--current-low", type=Path, required=True)
    parser.add_argument("--current-coverage", type=Path, required=True)
    parser.add_argument("--master-genus-map", type=Path, required=True)
    parser.add_argument("--ontology", type=Path, required=True)
    parser.add_argument(
        "--review-batch",
        choices=["20260807", "20260808_batch2"],
        default="20260807",
    )
    parser.add_argument("--run-id", default=RUN_ID)
    parser.add_argument("--source-artifact", default=SOURCE_ARTIFACT)
    parser.add_argument("--source-manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.candidates, dtype=str).fillna("")
    review = review_candidates(candidates, args.review_batch)
    accepted = review.loc[review["promotion_accepted"]].copy()
    new_common = common_evidence(
        review,
        run_id=args.run_id,
        source_artifact=args.source_artifact,
    )

    ontology = load_ontology(args.ontology)
    formal = direct_evidence_from_integrated(args.formal_lineage)
    prior = pd.read_csv(args.prior_promoted_common, dtype=str).fillna("")
    wfo = pd.read_csv(args.wfo_common, dtype=str).fillna("")
    additional = [
        pd.read_csv(path, dtype=str).fillna("")
        for path in args.additional_common
    ]
    raw = pd.concat(
        [formal, prior, wfo, *additional, new_common],
        ignore_index=True,
    ).fillna("")

    invalidated_rows = []
    invalidated_mask = pd.Series(False, index=raw.index)
    invalidations = (
        BATCH2_INVALIDATED_PRIOR
        if args.review_batch == "20260808_batch2"
        else INVALIDATED_PRIOR
    )
    for key, reason in invalidations.items():
        species, trait, value, lineage = key
        mask = (
            raw["accepted_species"].eq(species)
            & raw["trait_name"].eq(trait)
            & raw["normalized_value"].eq(value)
            & raw["source_lineage"].eq(lineage)
        )
        rows = raw.loc[mask].copy()
        rows["invalidation_reason"] = reason
        invalidated_rows.append(rows)
        invalidated_mask |= mask
    invalidated = (
        pd.concat(invalidated_rows, ignore_index=True)
        if invalidated_rows
        else pd.DataFrame(columns=[*EVIDENCE_COLUMNS, "invalidation_reason"])
    )
    raw = raw.loc[~invalidated_mask].reset_index(drop=True)

    current_direct = pd.read_csv(args.current_direct, dtype=str).fillna("")
    current_low = pd.read_csv(args.current_low, dtype=str).fillna("")
    current_coverage = pd.read_csv(args.current_coverage, dtype=str).fillna("")
    master = master_species(
        current_coverage,
        current_low,
        args.master_genus_map,
    )
    genus_map = dict(master[["accepted_species", "genus"]].itertuples(index=False, name=None))

    lineages, lineage_duplicates = dedupe_direct_lineages(raw, ontology)
    lineages["genus"] = (
        lineages["accepted_species"]
        .map(genus_map)
        .fillna(lineages["accepted_species"].str.split().str[0])
    )
    current_direct["genus"] = (
        current_direct["accepted_species"]
        .map(genus_map)
        .fillna(current_direct["accepted_species"].str.split().str[0])
    )

    # The WFO artifact is the validated parent state.  Recompute only the
    # genus x trait partitions touched by reviewed Europe PMC evidence or an
    # explicitly invalidated old row.  Rebuilding unrelated rules from a
    # different historical input set would silently discard valid prior Low.
    affected = set(
        zip(
            accepted["accepted_species"]
            .map(genus_map)
            .fillna(accepted["accepted_species"].str.split().str[0]),
            accepted["trait_name"],
            strict=True,
        )
    )
    for row in invalidated.itertuples(index=False):
        affected.add(
            (
                genus_map.get(
                    text(row.accepted_species),
                    text(row.accepted_species).split()[0],
                ),
                text(row.trait_name),
            )
        )

    lineage_mask = [
        (genus, trait) in affected
        for genus, trait in zip(
            lineages["genus"],
            lineages["trait_name"],
            strict=True,
        )
    ]
    affected_lineages = lineages.loc[lineage_mask].copy()
    recomputed_direct, conflict_audit = resolve_direct_cells(affected_lineages)
    recomputed_direct["genus"] = (
        recomputed_direct["accepted_species"]
        .map(genus_map)
        .fillna(recomputed_direct["accepted_species"].str.split().str[0])
    )
    existing_direct_mask = [
        (genus, trait) in affected
        for genus, trait in zip(
            current_direct["genus"],
            current_direct["trait_name"],
            strict=True,
        )
    ]
    rebuilt_direct = pd.concat(
        [
            current_direct.loc[[not value for value in existing_direct_mask]],
            recomputed_direct,
        ],
        ignore_index=True,
    ).fillna("")
    rebuilt_direct = rebuilt_direct.drop_duplicates(
        ["accepted_species", "trait_name"],
        keep="last",
    )

    affected_direct_mask = [
        (genus, trait) in affected
        for genus, trait in zip(
            rebuilt_direct["genus"],
            rebuilt_direct["trait_name"],
            strict=True,
        )
    ]
    affected_low_mask = [
        (genus, trait) in affected
        for genus, trait in zip(
            current_low["genus"],
            current_low["trait_name"],
            strict=True,
        )
    ]
    rules = build_rule_audit(
        rebuilt_direct.loc[affected_direct_mask],
        affected_lineages,
        current_low.loc[affected_low_mask],
    )
    replacement_low = apply_genus_rules(
        master,
        rebuilt_direct,
        rules,
        "current_min3",
    )
    rebuilt_low = pd.concat(
        [
            current_low.loc[[not value for value in affected_low_mask]],
            replacement_low,
        ],
        ignore_index=True,
    ).fillna("")
    rebuilt_low = rebuilt_low.drop_duplicates(
        ["accepted_species", "trait_name"],
        keep="last",
    )
    direct_keys_index = pd.MultiIndex.from_frame(rebuilt_direct[["accepted_species", "trait_name"]])
    low_keys_index = pd.MultiIndex.from_frame(rebuilt_low[["accepted_species", "trait_name"]])
    rebuilt_low = rebuilt_low.loc[~low_keys_index.isin(direct_keys_index)].reset_index(drop=True)
    rebuilt_coverage = species_axis_coverage(
        master,
        rebuilt_direct,
        rebuilt_low,
    )

    before = coverage_snapshot(current_coverage)
    after = coverage_snapshot(rebuilt_coverage)
    before_cells = keyed(
        current_coverage.loc[current_coverage["quality"].ne("")],
        ["accepted_species", "axis"],
    )
    after_cells = keyed(
        rebuilt_coverage.loc[rebuilt_coverage["quality"].ne("")],
        ["accepted_species", "axis"],
    )
    before_direct = keyed(current_direct, ["accepted_species", "trait_name"])
    after_direct = keyed(rebuilt_direct, ["accepted_species", "trait_name"])
    before_low = keyed(current_low, ["accepted_species", "trait_name"])
    after_low = keyed(rebuilt_low, ["accepted_species", "trait_name"])
    before_values = {
        (species, trait): value
        for species, trait, value in current_direct[
            ["accepted_species", "trait_name", "normalized_value"]
        ].itertuples(index=False, name=None)
    }
    after_values = {
        (species, trait): value
        for species, trait, value in rebuilt_direct[
            ["accepted_species", "trait_name", "normalized_value"]
        ].itertuples(index=False, name=None)
    }
    changed_values = sum(
        before_values[key] != after_values[key]
        for key in before_values.keys() & after_values.keys()
    )

    review_correct = int(review["accepted_correct"].sum())
    source_manifest = (
        json.loads(args.source_manifest.read_text(encoding="utf-8"))
        if args.source_manifest
        else {}
    )
    summary = {
        "contract": "europe_pmc_reproduction_reviewed_integration_v1",
        "created_at_utc": utc_now(),
        "source": {
            "candidate_rows": len(candidates),
            "all_reviewed": len(review),
            "accepted_correct": review_correct,
            "audit_precision": review_correct / len(review),
            "promoted_rows": len(accepted),
            "promoted_species_trait": int(
                accepted[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "corroborating_not_counted": int(
                review["review_decision"].eq("corroborating_not_counted").sum()
            ),
            "rejected_rows": int(review["review_decision"].eq("rejected").sum()),
            "queries": int(source_manifest.get("queries", 30)),
            "search_hits": int(source_manifest.get("search_hits", 753)),
            "fulltext_xml": int(
                source_manifest.get("xml_fetched_or_cached", 602)
            ),
            "query_cost_usd": float(
                source_manifest.get("query_cost_usd", 0.0)
            ),
            "invalidated_prior_misattributions": len(invalidated),
            "review_batch": args.review_batch,
            "source_run_id": args.run_id,
            "source_artifact": args.source_artifact,
        },
        "direct": {
            "before_species_trait": len(before_direct),
            "after_species_trait": len(after_direct),
            "added_species_trait": len(after_direct - before_direct),
            "removed_species_trait": len(before_direct - after_direct),
            "value_changed_species_trait": changed_values,
            "excluded_conflicts": int(conflict_audit["resolution_status"].eq("excluded").sum()),
        },
        "validated_low": {
            "before_species_trait": len(before_low),
            "after_species_trait": len(after_low),
            "added_species_trait": len(after_low - before_low),
            "invalidated_species_trait": len(before_low - after_low),
            "upgraded_to_direct_species_trait": len(before_low & (after_direct - before_direct)),
        },
        "coverage": {
            "before": before,
            "after": after,
            "gross_gain": len(after_cells - before_cells),
            "loss": len(before_cells - after_cells),
            "net_change": len(after_cells) - len(before_cells),
        },
        "safeguards": {
            "fulltext_quote_required": True,
            "search_snippets_as_evidence": False,
            "trait_specific_genus_join": "genus_x_trait_name",
            "axis_only_join": False,
            "family_inference": False,
            "global_fallback": False,
            "self_fertile_to_autonomous_selfing": False,
            "same_origin_restatements_independent": False,
        },
        "verification": {
            "ledger_rows": len(rebuilt_coverage),
            "unique_species": int(rebuilt_coverage["accepted_species"].nunique()),
            "duplicate_species_axis": int(
                rebuilt_coverage.duplicated(["accepted_species", "axis"]).sum()
            ),
            "direct_low_overlap": len(
                after_direct & keyed(rebuilt_low, ["accepted_species", "trait_name"])
            ),
            "quality_sum_matches_filled": (
                sum(after["quality_counts"].values()) == after["filled_species_axis"]
            ),
        },
    }

    review.to_csv(
        args.output / "fulltext_candidate_review.csv.gz",
        index=False,
        compression="gzip",
    )
    accepted.to_csv(
        args.output / "accepted_fulltext_evidence.csv.gz",
        index=False,
        compression="gzip",
    )
    new_common.to_csv(
        args.output / "accepted_common_direct_evidence.csv.gz",
        index=False,
        compression="gzip",
    )
    invalidated.to_csv(
        args.output / "invalidated_prior_direct_misattributions.csv.gz",
        index=False,
        compression="gzip",
    )
    lineage_duplicates.to_csv(
        args.output / "source_lineage_duplicates.csv.gz",
        index=False,
        compression="gzip",
    )
    conflict_audit.to_csv(
        args.output / "source_lineage_conflicts.csv.gz",
        index=False,
        compression="gzip",
    )
    rules.to_csv(
        args.output / "trait_specific_genus_rule_audit.csv.gz",
        index=False,
        compression="gzip",
    )
    rebuilt_direct.to_csv(
        args.output / "resolved_direct_species_trait.csv.gz",
        index=False,
        compression="gzip",
    )
    rebuilt_low.to_csv(
        args.output / "rebuilt_all_evidence_validated_low.csv.gz",
        index=False,
        compression="gzip",
    )
    rebuilt_coverage.to_csv(
        args.output / "strict_species_axis_coverage.csv.gz",
        index=False,
        compression="gzip",
    )
    (args.output / "integration_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    write_manifest(args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
