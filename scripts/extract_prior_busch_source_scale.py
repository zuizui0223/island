"""Source-scale audit and reviewed packet builder for Prior & Busch (2021).

The Meyer et al. public workbook retains the Prior & Busch population rows with
`quant = tm`, where tm is the multilocus outcrossing rate.  Prior & Busch define
realized selfing as 1 - tm and calculate species-level summaries across the
population estimates.  The repository already uses the same multilocus
outcrossing-rate mating-system contract for the Moeller source:

    tm < 0.2   -> predominantly_selfing
    0.2-0.8    -> mixed_mating
    tm > 0.8   -> predominantly_outcrossing

This adapter preserves every population row and its range/heterogeneity, then
uses the arithmetic species mean tm as the species-level aggregation contract.
Only exact fixed-universe species that are currently unresolved on the
reproductive-assurance axis enter the reviewed packet.  The packet remains
unpromoted until the formal cumulative source-batch integrator consumes it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
AXIS = "reproductive_assurance"
TRAIT = "mating_system"
SOURCE_REPO = "elenacicada/interesting_flowers"
SOURCE_COMMIT = "3d57315fe77097fe582025b884917550139801e8"
SOURCE_BLOB = "6c1f7de0f3e5a6c62d1fca2647d7adf33a19c380"
SOURCE_ARTIFACT = "Input_Data_MS.xlsx"
MEYER_DOI = "10.5061/dryad.cc2fqz6hr"
PRIOR_BUSCH_DATA_DOI = "10.5061/dryad.3j9kd51jw"
PRIOR_BUSCH_ARTICLE_DOI = "10.1002/ajb2.1766"
PRIOR_BUSCH_DATA_URL = f"https://doi.org/{PRIOR_BUSCH_DATA_DOI}"
PRIOR_BUSCH_ARTICLE_URL = f"https://doi.org/{PRIOR_BUSCH_ARTICLE_DOI}"
REVIEW_STATUS = "source_methodology_reviewed_reference_backed_strict_direct"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def locate_prior_busch_sheet(path: Path) -> tuple[str, pd.DataFrame]:
    book = pd.ExcelFile(path)
    hits: list[tuple[str, pd.DataFrame]] = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet, dtype=object).fillna("")
        lookup = {str(c).strip().casefold(): c for c in frame.columns}
        if "ref" not in lookup or "genus_species" not in lookup:
            continue
        ref_col = lookup["ref"]
        subset = frame.loc[
            frame[ref_col].map(text).str.casefold().str.contains("prior and busch", regex=False)
        ].copy()
        if not subset.empty:
            hits.append((sheet, subset))
    if len(hits) != 1:
        raise ValueError(
            f"expected exactly one Prior & Busch source sheet, found {[h[0] for h in hits]}"
        )
    return hits[0]


def mating_system_from_tm(value: float) -> str:
    """Apply the repository's established multilocus-outcrossing-rate contract."""
    if not 0 <= value <= 1:
        return ""
    if value < 0.2:
        return "predominantly_selfing"
    if value > 0.8:
        return "predominantly_outcrossing"
    return "mixed_mating"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-xlsx", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    sheet, source = locate_prior_busch_sheet(args.input_xlsx)
    col = {str(c).strip().casefold(): c for c in source.columns}
    required_source = {"genus_species", "quant", "ref", "notes"}
    missing_source = required_source.difference(col)
    if missing_source:
        raise ValueError(f"Prior & Busch block lacks columns: {sorted(missing_source)}")

    species_col = col["genus_species"]
    quant_col = col["quant"]
    ref_col = col["ref"]
    notes_col = col["notes"]

    source = source.copy()
    source["source_sheet"] = sheet
    source["source_row_in_sheet"] = source.index + 2
    source["accepted_species_candidate"] = (
        source[species_col].map(text).str.replace("_", " ", regex=False)
    )
    source["quant_raw"] = source[quant_col].map(text)
    source["tm"] = pd.to_numeric(source[quant_col], errors="coerce")
    source["notes_raw"] = source[notes_col].map(text)
    source["notes_confirm_quant_is_tm"] = source["notes_raw"].str.casefold().str.contains(
        r"quant\s*=\s*tm", regex=True
    )
    if not source["notes_confirm_quant_is_tm"].all():
        bad = source.loc[~source["notes_confirm_quant_is_tm"], [species_col, notes_col]].head(10)
        raise ValueError(f"Prior & Busch quant semantics are not uniformly tm:\n{bad}")
    if source["tm"].isna().any() or not source["tm"].between(0, 1, inclusive="both").all():
        raise ValueError("Prior & Busch tm must be numeric and in [0,1] for every retained row")

    source["selfing_rate_1_minus_tm"] = 1.0 - source["tm"]
    source["population_mating_system_from_tm"] = source["tm"].map(mating_system_from_tm)
    source["exact_binomial"] = source["accepted_species_candidate"].map(
        lambda x: bool(BINOMIAL.fullmatch(x))
    )
    source["trait_name"] = TRAIT
    source["axis"] = AXIS
    source["source_lineage"] = f"doi:{PRIOR_BUSCH_DATA_DOI}"
    source["evidence_scope"] = "source_scale_population_multilocus_outcrossing_rate"
    source["genus_rule_training_allowed"] = False

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[
        coverage["axis"].eq(AXIS), ["accepted_species", "quality"]
    ].drop_duplicates()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected fixed 106295-species reproductive coverage, got {len(rep)}")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    universe = set(quality)
    source["in_fixed_universe_exact"] = source["accepted_species_candidate"].isin(universe)
    source["current_reproductive_quality"] = (
        source["accepted_species_candidate"].map(quality).fillna("")
    )
    source["currently_unresolved_reproductive"] = (
        source["in_fixed_universe_exact"] & source["current_reproductive_quality"].eq("")
    )

    eligible_rows = source.loc[
        source["currently_unresolved_reproductive"] & source["exact_binomial"]
    ].copy()

    species_rows: list[dict[str, object]] = []
    packet_rows: list[dict[str, object]] = []
    source_reference = (
        f"Prior & Busch (2021), article DOI {PRIOR_BUSCH_ARTICLE_DOI}; "
        f"data DOI {PRIOR_BUSCH_DATA_DOI}; Meyer mirror DOI {MEYER_DOI}"
    )
    for species, group in eligible_rows.groupby("accepted_species_candidate", sort=True):
        tm = group["tm"].astype(float)
        mean_tm = float(tm.mean())
        state = mating_system_from_tm(mean_tm)
        population_states = sorted(set(group["population_mating_system_from_tm"]) - {""})
        species_rows.append(
            {
                "accepted_species": species,
                "axis": AXIS,
                "trait_name": TRAIT,
                "source_rows": int(len(group)),
                "tm_min": float(tm.min()),
                "tm_max": float(tm.max()),
                "tm_mean": mean_tm,
                "selfing_rate_mean_1_minus_tm": float((1.0 - tm).mean()),
                "population_state_set": "|".join(population_states),
                "population_state_variable": len(population_states) > 1,
                "species_mean_mating_system": state,
                "aggregation_contract": "arithmetic mean of population multilocus outcrossing rate tm",
                "mapping_contract": "tm<0.2 predominantly_selfing; 0.2<=tm<=0.8 mixed_mating; tm>0.8 predominantly_outcrossing",
                "review_status": REVIEW_STATUS,
                "packet_eligible": bool(state),
            }
        )
        if not state:
            continue
        raw = (
            f"mean_tm={mean_tm:.6f}; n_population_rows={len(group)}; "
            f"tm_range={float(tm.min()):.6f}-{float(tm.max()):.6f}; "
            f"population_states={'|'.join(population_states)}"
        )
        packet_rows.append(
            {
                "accepted_species": species,
                "source_species_name": species,
                "axis": AXIS,
                "source_reference_raw": source_reference,
                "source_url": PRIOR_BUSCH_DATA_URL,
                "source_article_url": PRIOR_BUSCH_ARTICLE_URL,
                "source_lineage": f"doi:{PRIOR_BUSCH_DATA_DOI}",
                "name_match_method": "exact_accepted_species_in_fixed_universe",
                "unresolved_reproductive_before_source": True,
                "promotion_allowed": False,
                "genus_rule_training_allowed": False,
                "trait_name": TRAIT,
                "normalized_value": state,
                "quality": "medium",
                "source_column": "species_mean_multilocus_outcrossing_rate_tm",
                "source_raw_value": raw,
                "review_status": REVIEW_STATUS,
                "acceptance_basis": (
                    "Meyer source notes explicitly define quant=tm. Prior & Busch define tm as "
                    "multilocus outcrossing rate, use 1-tm for realized selfing, and calculate "
                    "species-level means. The repository already accepts the same tm thresholds "
                    "for mating_system. Population heterogeneity is retained in the audit and is "
                    "not erased by the species-mean classification."
                ),
            }
        )

    species = pd.DataFrame(species_rows)
    packet = pd.DataFrame(packet_rows)
    packet_columns = [
        "accepted_species", "source_species_name", "axis", "source_reference_raw",
        "source_url", "source_article_url", "source_lineage", "name_match_method",
        "unresolved_reproductive_before_source", "promotion_allowed",
        "genus_rule_training_allowed", "trait_name", "normalized_value", "quality",
        "source_column", "source_raw_value", "review_status", "acceptance_basis",
    ]
    if packet.empty:
        packet = pd.DataFrame(columns=packet_columns)
    else:
        packet = packet[packet_columns].sort_values("accepted_species").reset_index(drop=True)

    source.to_csv(
        args.output / "prior_busch_full_source_rows.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    species.to_csv(args.output / "prior_busch_unresolved_species_candidates.csv", index=False)
    packet.to_csv(
        args.output / "prior_busch_unresolved_mating_batch.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    numeric = source["tm"].dropna()
    summary = {
        "contract": "prior_busch_2021_source_scale_tm_review_v2",
        "source_repo": SOURCE_REPO,
        "source_commit": SOURCE_COMMIT,
        "source_blob": SOURCE_BLOB,
        "source_artifact": SOURCE_ARTIFACT,
        "source_sha256": sha256(args.input_xlsx),
        "meyer_source_doi": MEYER_DOI,
        "prior_busch_dataset_doi": PRIOR_BUSCH_DATA_DOI,
        "prior_busch_article_doi": PRIOR_BUSCH_ARTICLE_DOI,
        "source_sheet": sheet,
        "prior_busch_rows": int(len(source)),
        "raw_unique_taxon_labels": int(source["accepted_species_candidate"].nunique()),
        "exact_binomial_rows": int(source["exact_binomial"].sum()),
        "exact_fixed_universe_species": int(
            source.loc[source["in_fixed_universe_exact"], "accepted_species_candidate"].nunique()
        ),
        "exact_current_unresolved_reproductive_species": int(
            source.loc[
                source["currently_unresolved_reproductive"], "accepted_species_candidate"
            ].nunique()
        ),
        "numeric_tm_rows": int(numeric.size),
        "tm_min": float(numeric.min()) if len(numeric) else None,
        "tm_max": float(numeric.max()) if len(numeric) else None,
        "quant_semantics_verified_from_source_notes": bool(
            source["notes_confirm_quant_is_tm"].all()
        ),
        "reviewed_strict_packet_rows": int(len(packet)),
        "reviewed_strict_packet_state_counts": (
            packet["normalized_value"].value_counts().sort_index().astype(int).to_dict()
            if len(packet)
            else {}
        ),
        "formal_gain": 0,
        "formal_gain_pending_batch_integration": int(len(packet)),
        "promotion_allowed": False,
        "semantics_gate": "resolved: quant=tm; species mean tm mapped with established repository mating-system thresholds",
        "population_variation_policy": "retain population tm range/state set in audit; classify only the source-defined species mean",
    }
    (args.output / "prior_busch_source_scale_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
