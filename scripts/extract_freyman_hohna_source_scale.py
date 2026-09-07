"""Extract the complete Freyman & Höhna Onagraceae SC/SI dataset at source scale.

The published analysis repository exposes a complete species-state table with
`sc=0si=1`. This adapter reads the entire file, keeps only explicit 0/1 states,
requires exact binomial membership in the fixed 106,295-species analysis
universe, and measures novelty against the current reproductive-assurance
coverage. Unknown `?` rows are retained only in the exclusion audit.

Important provenance rule: this is one compiled source dataset. Species rows are
species-direct candidates, but all rows share one compilation source lineage and
`genus_rule_training_allowed=false`; the compilation must not masquerade as
hundreds of independent studies for Validated Low.

No formal coverage is changed by this script. Promotion remains false until the
completed source packet is reviewed and included in a later manual batch
integration.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
TRAIT = "self_incompatibility"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
SOURCE_REPO = "wf8/onagraceae"
SOURCE_BRANCH = "master"
SOURCE_BLOB = "5ac0a3d897f4c11f70e84fda7007597b921ec7e4"
SOURCE_FILE = "sse_analyses/data/selfing_data_complete.csv"
SOURCE_URL = f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_BRANCH}/{SOURCE_FILE}"
DATASET_DOI = "10.5061/dryad.9jj428f"
ARTICLE_DOI = "10.1093/sysbio/syy078"
COMPILATION_LINEAGE = f"doi:{ARTICLE_DOI}"


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-csv", type=Path, required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    source = pd.read_csv(args.source_csv, dtype=str).fillna("")
    required = {"taxa", "sc=0si=1"}
    if missing := required.difference(source.columns):
        raise ValueError(f"Freyman source lacks columns: {sorted(missing)}")

    source = source.copy()
    source["source_row"] = source.index + 2
    source["source_species_name"] = source["taxa"].map(text).str.replace("_", " ", regex=False)
    source["source_state"] = source["sc=0si=1"].map(text)
    source["exact_binomial"] = source["source_species_name"].map(lambda x: bool(BINOMIAL.fullmatch(x)))
    source["normalized_value"] = source["source_state"].map({"0": "SC", "1": "SI"}).fillna("")
    source["explicit_state"] = source["normalized_value"].ne("")

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected unique 106295-species reproductive coverage, got {len(rep)}")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    universe = set(quality)

    source["accepted_species"] = source["source_species_name"]
    source["in_fixed_universe_exact"] = source["accepted_species"].isin(universe)
    source["current_reproductive_quality"] = source["accepted_species"].map(quality).fillna("")
    source["currently_unresolved_reproductive"] = source["in_fixed_universe_exact"] & source["current_reproductive_quality"].eq("")

    duplicate_species = source.loc[
        source["explicit_state"] & source["accepted_species"].duplicated(keep=False),
        ["accepted_species", "normalized_value", "source_row"],
    ].copy()
    conflicts = (
        source.loc[source["explicit_state"]]
        .groupby("accepted_species")["normalized_value"]
        .nunique()
    )
    conflicting_species = set(conflicts.loc[conflicts.gt(1)].index)

    selected = source.loc[
        source["explicit_state"]
        & source["exact_binomial"]
        & source["in_fixed_universe_exact"]
        & source["currently_unresolved_reproductive"]
        & ~source["accepted_species"].isin(conflicting_species)
    ].copy()
    # One compiled state per species. Duplicate same-state rows collapse to one
    # candidate while the full source audit preserves every row.
    selected = selected.sort_values("source_row").drop_duplicates("accepted_species", keep="first")

    candidates = pd.DataFrame(
        {
            "accepted_species": selected["accepted_species"],
            "source_species_name": selected["source_species_name"],
            "axis": AXIS,
            "trait_name": TRAIT,
            "normalized_value": selected["normalized_value"],
            "quality": "high",
            "evidence_scope": "species_direct_published_compilation",
            "source_provider": "Freyman & Hohna Onagraceae mating-system dataset",
            "source_url": SOURCE_URL,
            "source_doi": DATASET_DOI,
            "source_article_doi": ARTICLE_DOI,
            "source_record_id": selected["source_row"].map(lambda n: f"freyman-hohna-row:{int(n)}"),
            "source_excerpt": selected.apply(
                lambda r: f"taxa={r['taxa']}; sc=0si=1={r['source_state']}; 0=SC, 1=SI",
                axis=1,
            ),
            "name_match_method": "exact_frozen_analysis_species",
            "source_lineage": COMPILATION_LINEAGE,
            "lineage_method": "single_published_compilation_not_independent_underlying_studies",
            "review_status": "source_scale_extracted_pending_batch_review",
            "promotion_allowed": False,
            "genus_rule_training_allowed": False,
            "current_public_cell_status": "unresolved",
        }
    )

    source["exclusion_reason"] = ""
    source.loc[~source["explicit_state"], "exclusion_reason"] = "unknown_or_nonbinary_source_state"
    source.loc[source["explicit_state"] & ~source["exact_binomial"], "exclusion_reason"] = "not_exact_binomial"
    source.loc[source["explicit_state"] & source["exact_binomial"] & ~source["in_fixed_universe_exact"], "exclusion_reason"] = "not_exact_fixed_universe_name"
    source.loc[source["explicit_state"] & source["in_fixed_universe_exact"] & ~source["currently_unresolved_reproductive"], "exclusion_reason"] = "already_resolved_reproductive_axis"
    source.loc[source["accepted_species"].isin(conflicting_species), "exclusion_reason"] = "within_source_species_state_conflict"

    source.to_csv(args.output / "freyman_hohna_full_source_audit.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    candidates.to_csv(args.output / "freyman_hohna_unresolved_direct_candidates.csv", index=False)
    duplicate_species.to_csv(args.output / "freyman_hohna_duplicate_species_rows.csv", index=False)

    state_counts = {str(k): int(v) for k, v in source["source_state"].value_counts().sort_index().items()}
    selected_counts = {str(k): int(v) for k, v in candidates["normalized_value"].value_counts().sort_index().items()} if len(candidates) else {}
    summary = {
        "contract": "freyman_hohna_source_scale_reproductive_v1",
        "source_repo": SOURCE_REPO,
        "source_branch": SOURCE_BRANCH,
        "source_blob": SOURCE_BLOB,
        "source_file": SOURCE_FILE,
        "source_sha256": sha256(args.source_csv),
        "dataset_doi": DATASET_DOI,
        "article_doi": ARTICLE_DOI,
        "source_rows": int(len(source)),
        "raw_unique_taxa": int(source["source_species_name"].nunique()),
        "source_state_counts": state_counts,
        "explicit_sc_si_rows": int(source["explicit_state"].sum()),
        "explicit_sc_si_species": int(source.loc[source["explicit_state"], "accepted_species"].nunique()),
        "exact_fixed_universe_sc_si_species": int(source.loc[source["explicit_state"] & source["in_fixed_universe_exact"], "accepted_species"].nunique()),
        "exact_current_unresolved_reproductive_species": int(source.loc[source["explicit_state"] & source["currently_unresolved_reproductive"], "accepted_species"].nunique()),
        "within_source_conflicting_species": int(len(conflicting_species)),
        "selected_unresolved_direct_candidates": int(candidates["accepted_species"].nunique()),
        "selected_by_value": selected_counts,
        "shared_compilation_source_lineage": COMPILATION_LINEAGE,
        "genus_rule_training_allowed": False,
        "formal_gain": 0,
        "promotion_allowed": False,
        "next_gate": "review complete source packet, taxonomy mismatches, and provenance; batch integrate later only",
    }
    (args.output / "freyman_hohna_source_scale_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
