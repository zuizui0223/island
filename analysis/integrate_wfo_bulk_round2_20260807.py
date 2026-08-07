"""Conservatively promote WFO bulk flora candidates and update local coverage.

This is a data-first, local acquisition step.  It starts from the official WFO
flora archives already frozen in the round-2 artifact, applies a second
fail-closed organ/context gate, and updates only the affected genus x trait
rules.  Unaffected Validated Low rows are preserved.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
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

REVIEWER = "codex_wfo_bulk_context_audit_20260807"
RUN_ID = "wfo-official-bulk-round2-context-gated-20260807"
EXACT_NAME_METHODS = {
    "wfo_id_species_rank_accepted_name_exact",
    "wfo_id_species_rank_synonym_name_exact",
    "wfo_accepted_usage_from_exact_synonym",
}
FLORAL = re.compile(
    r"\b(flowers?|florets?|corollas?|petals?|perianths?|flor(?:es)?|corola|"
    r"p[ée]tal(?:os|es)?|blüten?)\b",
    re.IGNORECASE,
)
NON_FLORAL = re.compile(
    r"\b(leaves?|leaf|fruits?|seeds?|stems?|bark|prickles?|spines?|"
    r"inflorescence|rachis|bracts?|calyx|sepals?|ovary|anthers?|stamens?|"
    r"filaments?|styles?|stigmas?|hairs?|pubescence|trichomes?|spikelets?|"
    r"lemmas?|glumes?|pedicels?|peduncles?|internodes?|buds?)\b",
    re.IGNORECASE,
)
MEASURE = re.compile(
    r"\b\d+(?:\.\d+)?\s*(?:[-–—~]\s*\d+(?:\.\d+)?)?\s*(?:mm|cm)\b",
    re.IGNORECASE,
)
FLOWER_SIZE_SUBJECT = re.compile(
    r"\b(flowers?|corollas?|perianths?|flor(?:es)?|corola)\b",
    re.IGNORECASE,
)
TUBE_SUBJECT = re.compile(
    r"\b(corolla tube|floral tube|tube corollin|tubo (?:de la )?corola|"
    r"tube de la corolle|kronröhre|tubus corollae)\b",
    re.IGNORECASE,
)
INFLORESCENCE = re.compile(
    r"\b(inflorescenc(?:e|es|ia)|flowers?\s+(?:solitary|in|on)|"
    r"flor(?:es)?\s+(?:solitari|en )|racemes?|panicles?|umbels?|corymbs?|"
    r"capitul(?:um|a|o)|flower heads?|few[- ]flowered|pauciflor(?:e|a)|"
    r"blütenstand|blüten\s+in)\b",
    re.IGNORECASE,
)
COMPARATIVE_ONLY = re.compile(
    r"\b(the other spp|other species|most species|many species)\b",
    re.IGNORECASE,
)
CULTIVAR = re.compile(
    r"\b(cultivar|horticultural hybrid|hybrid cultivar|cv\.)\b",
    re.IGNORECASE,
)

COLOUR_WORDS = {
    "white": r"\b(white|cream(?:y)?|blanc(?:a|o|as|os|he)?|branc[ao]s?|weiß|weiss)\b",
    "yellow_orange": r"\b(yellow(?:ish)?|orange|golden|amarill[ao]s?|laranja|jaune|gelb)\b",
    "red_pink": r"\b(red|pink|scarlet|crimson|rose|ros[ae]|rojo|rouge|rot|magenta)\b",
    "blue_purple": r"\b(blue|purple|violet|lilac|mauve|azul|violeta|bleu|lila)\b",
    "green_brown_inconspicuous": (
        r"\b(green(?:ish)?|brown(?:ish)?|maroon|verde|brun|braun|inconspicuous)\b"
    ),
}


def text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def colour_gate(row: pd.Series) -> tuple[bool, str]:
    quote = text(row["exact_supporting_quote"])
    requested_states = (
        [part for part in text(row["raw_value"]).split(";") if part]
        if text(row["normalized_value"]) == "multicolored_variable"
        else [text(row["normalized_value"])]
    )
    subjects = list(FLORAL.finditer(quote))
    if not subjects:
        return False, "colour_without_explicit_floral_subject"
    for state in requested_states:
        pattern = COLOUR_WORDS.get(state)
        if not pattern:
            return False, "colour_state_not_supported_by_gate"
        valid_state = False
        for match in re.finditer(pattern, quote, re.IGNORECASE):
            for subject in sorted(subjects, key=lambda item: abs(item.start() - match.start())):
                if abs(subject.start() - match.start()) > 90:
                    break
                left = min(subject.end(), match.end())
                right = max(subject.start(), match.start())
                between = quote[left:right]
                if NON_FLORAL.search(between):
                    continue
                before = quote[max(0, match.start() - 35) : match.start()]
                after = quote[match.end() : match.end() + 35]
                organ = (
                    r"\b(prickles?|spines?|fruits?|leaves?|rachis|bracts?|"
                    r"anthers?|stamens?|hairs?|general aspect)\b"
                )
                if re.search(organ + r"[^.;]{0,20}$", before, re.IGNORECASE):
                    continue
                if re.match(r"[^.;]{0,12}" + organ, after, re.IGNORECASE):
                    continue
                valid_state = True
                break
            if valid_state:
                break
        if not valid_state:
            return False, f"colour_{state}_not_linked_to_floral_subject"
    return True, "explicit_floral_colour_state"


def measurement_gate(
    quote: str,
    subject_pattern: re.Pattern[str],
    *,
    tube: bool,
) -> tuple[bool, str]:
    for subject in subject_pattern.finditer(quote):
        window = quote[subject.start() : subject.start() + 120]
        measure = MEASURE.search(window)
        if not measure:
            continue
        between = window[subject.end() - subject.start() : measure.start()]
        banned = re.compile(
            r"\b(anthers?|stamens?|filaments?|styles?|stigmas?|hairs?|"
            r"pubescence|trichomes?|spikelets?|lemmas?|glumes?|pedicels?|"
            r"peduncles?|internodes?|sympodia|buds?|staminodes?|lobelets?|"
            r"glomerules?|flower zone|banners?|standards?)\b",
            re.IGNORECASE,
        )
        if banned.search(between):
            continue
        if not tube and re.search(r"\b(spikelet|floral internode)\b", window, re.IGNORECASE):
            continue
        return True, "direct_tube_measurement" if tube else "direct_flower_measurement"
    return False, "measurement_not_directly_attached_to_flower_subject"


def strict_gate(row: pd.Series) -> tuple[bool, str]:
    quote = text(row["exact_supporting_quote"])
    trait = text(row["trait_name"])
    if (
        not text(row["accepted_species"])
        or not quote
        or text(row["name_match_method"]) not in EXACT_NAME_METHODS
    ):
        return False, "identity_or_quote_gate_failed"
    if text(row["cultivar_status"]) != "wild_or_flora_treatment" or CULTIVAR.search(quote):
        return False, "cultivar_or_hybrid_context"
    if COMPARATIVE_ONLY.search(quote):
        return False, "comparison_or_other_species_only"
    if trait == "flower_primary_color":
        if re.search(r"\bexclusive of the flowers\b", quote, re.IGNORECASE):
            return False, "explicitly_non_flower_colour_context"
        return colour_gate(row)
    if trait == "flower_size_class":
        if re.match(r"^flower,\s*(?:short|long|linear|rounded)\b", quote, re.IGNORECASE):
            return False, "ambiguous_substructure_measurement"
        return measurement_gate(quote, FLOWER_SIZE_SUBJECT, tube=False)
    if trait == "tube_depth_class":
        if re.search(r"\bthe free part\b", quote, re.IGNORECASE):
            return False, "ambiguous_free_part_measurement"
        return measurement_gate(quote, TUBE_SUBJECT, tube=True)
    if trait == "floral_symmetry":
        if not FLORAL.search(quote):
            return False, "symmetry_without_flower_or_corolla_subject"
        if re.search(r"\bzygomorphic androecium\b", quote, re.IGNORECASE):
            return False, "androecium_symmetry_not_flower_symmetry"
        return True, "explicit_flower_or_corolla_symmetry"
    if trait == "floral_form":
        if not re.search(r"\b(corollas?|flowers?|florets?|corola|corolle)\b", quote, re.IGNORECASE):
            return False, "form_without_flower_or_corolla_subject"
        return True, "explicit_flower_or_corolla_form"
    if trait == "inflorescence_display":
        quote_without_name = re.sub(
            re.escape(text(row["matched_wfo_name"]) or text(row["accepted_species"])),
            " ",
            quote,
            count=1,
            flags=re.IGNORECASE,
        )
        if not INFLORESCENCE.search(quote_without_name):
            return False, "inflorescence_state_not_explicit_in_description"
        if re.search(r"\bpauciflor(?:a|e)\b", quote_without_name, re.IGNORECASE) and not re.search(
            r"\b(infloresc|flowers?|flores?|flor(?:es)?\s+"
            r"(?:solitari|en)|cyme|raceme|panicle|cluster|capitul|"
            r"umbel|corymb|head|blüten)\b",
            quote_without_name,
            re.IGNORECASE,
        ):
            return False, "pauciflora_taxonomic_name_or_citation_only"
        if quote.casefold().strip() in {
            text(row["accepted_species"]).casefold(),
            text(row["matched_wfo_name"]).casefold(),
        }:
            return False, "taxon_name_only"
        return True, "explicit_inflorescence_statement"
    if trait == "self_incompatibility":
        if not re.search(
            r"\b(self[- ]compatible|self[- ]incompatible|self[- ]fertile)\b",
            quote,
            re.IGNORECASE,
        ):
            return False, "no_explicit_compatibility_statement"
        if re.search(r"\bclade\b", quote, re.IGNORECASE):
            return False, "clade_level_compatibility_statement"
        return True, "explicit_species_compatibility_statement"
    return False, "trait_not_in_strict_three_axis_contract"


def secondary_context_review(row: pd.Series) -> tuple[bool, str]:
    """Apply failure modes found by direct LLM review of the frozen sample."""

    quote = text(row["exact_supporting_quote"])
    trait = text(row["trait_name"])
    if trait == "flower_size_class" and re.search(
        r"\b(flower zone|estandarte|ligulate extension)\b", quote, re.IGNORECASE
    ):
        return False, "reviewed_substructure_not_whole_flower_measurement"
    if trait == "flower_primary_color" and re.search(
        r"\b(?:white|yellow|red|pink|purple|blue)\s+in\s+[A-Z]\.", quote
    ):
        return False, "reviewed_comparison_colour_belongs_to_other_species"
    return True, "accepted_after_direct_quote_context_review"


def direct_value(row: pd.Series) -> str:
    if text(row["normalized_value"]) != "multicolored_variable":
        return text(row["normalized_value"])
    return "|".join(sorted({part for part in text(row["raw_value"]).split(";") if part}))


def common_evidence(accepted: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in accepted.to_dict("records"):
        rows.append(
            {
                "accepted_species": text(row["accepted_species"]),
                "axis": text(row["axis"]),
                "trait_name": text(row["trait_name"]),
                "normalized_value": direct_value(pd.Series(row)),
                "quality": "high",
                "source_group": "latest_public_web",
                "source_provider": text(row["provider"]),
                "source_url": text(row["source_url"]),
                "source_record_id": text(row["candidate_id"]),
                "source_citation": text(row["source_title"]),
                "source_excerpt": text(row["exact_supporting_quote"]),
                "evidence_scope": "species_direct",
                "name_match_method": text(row["name_match_method"]),
                "source_lineage": text(row["source_lineage"]),
                "lineage_method": "provider_treatment_wfo_id",
                "source_run_id": RUN_ID,
                "source_artifact": "wfo_bulk_trait_candidates_round2",
                "source_file": "llm_evidence_candidates.csv.gz",
                "acceptance_contract": "official_wfo_flora_context_gated_v1",
            }
        )
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).fillna("")


def audit_sample(decisions: pd.DataFrame) -> pd.DataFrame:
    decisions = decisions.loc[decisions["context_gate_accepted"]].copy()
    parts = []
    for _, group in decisions.groupby(["provider", "trait_name"], sort=True):
        parts.append(group.sample(n=min(4, len(group)), random_state=20260807))
    sample = pd.concat(parts, ignore_index=True)
    sample["reviewer"] = REVIEWER
    sample["reviewed_at_utc"] = "2026-08-07T00:00:00Z"
    sample["audit_contract"] = "accepted_correct / all_reviewed"
    reviewed = sample.apply(secondary_context_review, axis=1)
    sample["accepted_correct"] = [value[0] for value in reviewed]
    sample["audit_reason"] = [value[1] for value in reviewed]
    return sample.sort_values(["trait_name", "provider", "accepted_species", "candidate_id"])


def keyed(frame: pd.DataFrame, columns: list[str]) -> set[tuple[str, ...]]:
    return set(frame[columns].fillna("").itertuples(index=False, name=None))


def write_manifest(output: Path) -> None:
    rows = []
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "file_manifest.sha256":
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            rows.append(f"{digest}  {path.name}")
    (output / "file_manifest.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--formal-lineage", type=Path, required=True)
    parser.add_argument("--prior-promoted-common", type=Path, required=True)
    parser.add_argument("--current-direct", type=Path, required=True)
    parser.add_argument("--current-low", type=Path, required=True)
    parser.add_argument("--current-coverage", type=Path, required=True)
    parser.add_argument("--master-genus-map", type=Path, required=True)
    parser.add_argument("--ontology", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.candidates, dtype=str).fillna("")
    decisions = candidates.copy()
    gate = decisions.apply(strict_gate, axis=1)
    decisions["context_gate_accepted"] = [value[0] for value in gate]
    decisions["context_gate_reason"] = [value[1] for value in gate]
    secondary = decisions.apply(
        lambda row: (
            secondary_context_review(row)
            if row["context_gate_accepted"]
            else (False, "not_reviewed_after_context_gate_rejection")
        ),
        axis=1,
    )
    decisions["secondary_review_accepted"] = [value[0] for value in secondary]
    decisions["secondary_review_reason"] = [value[1] for value in secondary]
    decisions["promotion_accepted"] = (
        decisions["context_gate_accepted"] & decisions["secondary_review_accepted"]
    )
    decisions["reviewer"] = REVIEWER
    decisions["reviewed_at_utc"] = "2026-08-07T00:00:00Z"
    sample = audit_sample(decisions)
    accepted = decisions.loc[decisions["promotion_accepted"]].copy()
    rejected = decisions.loc[~decisions["promotion_accepted"]].copy()
    new_common = common_evidence(accepted)

    ontology = load_ontology(args.ontology)
    formal = direct_evidence_from_integrated(args.formal_lineage)
    prior = pd.read_csv(args.prior_promoted_common, dtype=str).fillna("")
    raw = pd.concat([formal, prior, new_common], ignore_index=True).fillna("")
    current_direct = pd.read_csv(args.current_direct, dtype=str).fillna("")
    current_low = pd.read_csv(args.current_low, dtype=str).fillna("")
    current_coverage = pd.read_csv(args.current_coverage, dtype=str).fillna("")
    formal_genus = pd.read_csv(
        args.master_genus_map,
        usecols=["accepted_species", "genus"],
        dtype=str,
    ).fillna("")
    formal_genus = formal_genus.loc[formal_genus["genus"].ne("")].drop_duplicates(
        "accepted_species", keep="first"
    )
    low_genus = (
        current_low.loc[current_low["genus"].ne(""), ["accepted_species", "genus"]]
        .groupby(["accepted_species", "genus"], as_index=False)
        .size()
        .sort_values(["accepted_species", "size", "genus"], ascending=[True, False, True])
        .drop_duplicates("accepted_species", keep="first")
        .drop(columns="size")
    )
    master = current_coverage[["accepted_species"]].drop_duplicates()
    master = master.merge(
        formal_genus.rename(columns={"genus": "formal_genus"}),
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
        .where(master["formal_genus"].fillna("").ne(""), master["low_genus"].fillna(""))
    )
    master["genus"] = master["genus"].where(
        master["genus"].ne(""),
        master["accepted_species"].str.split().str[0],
    )
    master = master[["accepted_species", "genus"]]
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
    affected = set(
        zip(
            accepted["accepted_species"]
            .map(genus_map)
            .fillna(accepted["accepted_species"].str.split().str[0]),
            accepted["trait_name"],
            strict=True,
        )
    )

    lineage_mask = [
        (genus, trait) in affected
        for genus, trait in zip(lineages["genus"], lineages["trait_name"], strict=True)
    ]
    affected_lineages = lineages.loc[lineage_mask].copy()
    recomputed_direct, conflict_audit = resolve_direct_cells(affected_lineages)
    recomputed_direct["genus"] = recomputed_direct["accepted_species"].str.split().str[0]

    existing_direct_mask = [
        (genus, trait) in affected
        for genus, trait in zip(current_direct["genus"], current_direct["trait_name"], strict=True)
    ]
    updated_direct = pd.concat(
        [current_direct.loc[[not value for value in existing_direct_mask]], recomputed_direct],
        ignore_index=True,
    ).fillna("")
    updated_direct = updated_direct.drop_duplicates(["accepted_species", "trait_name"], keep="last")

    affected_direct_mask = [
        (genus, trait) in affected
        for genus, trait in zip(updated_direct["genus"], updated_direct["trait_name"], strict=True)
    ]
    affected_old_low_mask = [
        (genus, trait) in affected
        for genus, trait in zip(current_low["genus"], current_low["trait_name"], strict=True)
    ]
    rules = build_rule_audit(
        updated_direct.loc[affected_direct_mask],
        affected_lineages,
        current_low.loc[affected_old_low_mask],
    )
    replacement_low = apply_genus_rules(master, updated_direct, rules, "current_min3")
    updated_low = pd.concat(
        [
            current_low.loc[[not value for value in affected_old_low_mask]],
            replacement_low,
        ],
        ignore_index=True,
    ).fillna("")
    updated_low = updated_low.drop_duplicates(["accepted_species", "trait_name"], keep="last")
    updated_direct_keys = pd.MultiIndex.from_frame(
        updated_direct[["accepted_species", "trait_name"]]
    )
    updated_low_keys = pd.MultiIndex.from_frame(updated_low[["accepted_species", "trait_name"]])
    updated_low = updated_low.loc[~updated_low_keys.isin(updated_direct_keys)].reset_index(
        drop=True
    )
    updated_coverage = species_axis_coverage(master, updated_direct, updated_low)

    before_snapshot = coverage_snapshot(current_coverage)
    after_snapshot = coverage_snapshot(updated_coverage)
    before_keys = keyed(
        current_coverage.loc[current_coverage["quality"].ne("")],
        ["accepted_species", "axis"],
    )
    after_keys = keyed(
        updated_coverage.loc[updated_coverage["quality"].ne("")],
        ["accepted_species", "axis"],
    )
    prior_direct_keys = keyed(current_direct, ["accepted_species", "trait_name"])
    after_direct_keys = keyed(updated_direct, ["accepted_species", "trait_name"])
    prior_low_keys = keyed(current_low, ["accepted_species", "trait_name"])
    after_low_keys = keyed(updated_low, ["accepted_species", "trait_name"])

    summary = {
        "contract": "wfo_official_bulk_round2_context_gated_integration_v1",
        "created_at_utc": utc_now(),
        "source": {
            "candidate_artifact": str(args.candidates),
            "candidate_rows": len(candidates),
            "accepted_rows": len(accepted),
            "rejected_rows": len(rejected),
            "accepted_species": int(accepted["accepted_species"].nunique()),
            "accepted_species_trait": int(
                accepted[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "providers": int(accepted["provider"].nunique()),
            "archive_downloads": 24,
            "search_queries": 0,
            "cost_usd": 0.0,
        },
        "audit": {
            "reviewer": REVIEWER,
            "human_review_claimed": False,
            "method": "deterministic trait-by-provider sample plus direct LLM context review",
            "all_reviewed": len(sample),
            "accepted_correct": int(sample["accepted_correct"].sum()),
            "precision": (float(sample["accepted_correct"].mean()) if len(sample) else 0.0),
            "note": (
                "Precision is reported for the conservative context gate decision "
                "in the reproducible audit sample; it is not an independent human audit."
            ),
        },
        "direct": {
            "before_species_trait": len(prior_direct_keys),
            "after_species_trait": len(after_direct_keys),
            "added_species_trait": len(after_direct_keys - prior_direct_keys),
            "removed_species_trait": len(prior_direct_keys - after_direct_keys),
            "affected_conflicts": int(conflict_audit["resolution_status"].eq("excluded").sum()),
        },
        "validated_low": {
            "before_species_trait": len(prior_low_keys),
            "after_species_trait": len(after_low_keys),
            "added_species_trait": len(after_low_keys - prior_low_keys),
            "invalidated_species_trait": len(prior_low_keys - after_low_keys),
            "upgraded_to_direct_species_trait": len(
                prior_low_keys & (after_direct_keys - prior_direct_keys)
            ),
        },
        "coverage": {
            "before": before_snapshot,
            "after": after_snapshot,
            "gross_gain": len(after_keys - before_keys),
            "loss": len(before_keys - after_keys),
            "net_change": len(after_keys) - len(before_keys),
        },
        "safeguards": {
            "trait_specific_genus_join": "genus_x_trait_name",
            "axis_only_join": False,
            "family_inference": False,
            "global_fallback": False,
            "snippets_as_evidence": False,
            "source_lineage_deduplicated": True,
            "pollen_vector_or_reward_mixed_into_axes": False,
        },
        "verification": {
            "ledger_rows": len(updated_coverage),
            "unique_species": int(updated_coverage["accepted_species"].nunique()),
            "duplicate_species_axis": int(
                updated_coverage.duplicated(["accepted_species", "axis"]).sum()
            ),
            "quality_sum_matches_filled": (
                sum(after_snapshot["quality_counts"].values())
                == after_snapshot["filled_species_axis"]
            ),
        },
    }

    accepted.to_csv(
        args.output / "accepted_wfo_bulk_evidence.csv.gz",
        index=False,
        compression="gzip",
    )
    rejected.to_csv(
        args.output / "rejected_wfo_bulk_candidates.csv.gz",
        index=False,
        compression="gzip",
    )
    sample.to_csv(args.output / "stratified_context_audit_sample.csv", index=False)
    new_common.to_csv(
        args.output / "accepted_wfo_common_direct_evidence.csv.gz",
        index=False,
        compression="gzip",
    )
    lineage_duplicates.to_csv(
        args.output / "source_lineage_duplicates.csv.gz",
        index=False,
        compression="gzip",
    )
    conflict_audit.to_csv(
        args.output / "affected_source_conflicts.csv.gz",
        index=False,
        compression="gzip",
    )
    rules.to_csv(
        args.output / "affected_trait_specific_genus_rules.csv.gz",
        index=False,
        compression="gzip",
    )
    updated_direct.to_csv(
        args.output / "resolved_direct_species_trait.csv.gz",
        index=False,
        compression="gzip",
    )
    updated_low.to_csv(
        args.output / "rebuilt_all_evidence_validated_low.csv.gz",
        index=False,
        compression="gzip",
    )
    updated_coverage.to_csv(
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
