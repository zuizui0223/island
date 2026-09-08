from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

SPECIES = "Hauya heydeana"
AXIS = "reproductive_assurance"
TRAIT = "self_incompatibility"
EXPECTED_SOURCE_SHA256 = "76d2357323d3cc5ad2a55370eb2a780c062feb567d5328bf6af2e4c74377b489"
FREYMAN_COMMIT = "6099e27bb238a5e85fb0d8249553e780302e3793"
FREYMAN_URL = (
    "https://raw.githubusercontent.com/wf8/onagraceae/"
    + FREYMAN_COMMIT
    + "/sse_analyses/data/selfing_data_complete.csv"
)
RAVEN_DOI = "10.1080/0028825X.1979.10432572"
RAVEN_URL = "https://doi.org/10.1080/0028825X.1979.10432572"
MONOGRAPH_URL = "https://repository.si.edu/bitstreams/a5edce1e-b430-4663-ad88-bd8166e496d0/download"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--source-csv", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--direct-ledger", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    observed_sha = sha256(args.source_csv)
    if observed_sha != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"Freyman source hash changed: {observed_sha}")
    source = pd.read_csv(args.source_csv, dtype=str).fillna("")
    if list(source.columns) != ["taxa", "sc=0si=1"] or len(source) != 357:
        raise ValueError("unexpected Freyman source schema/row count")
    hit = source.loc[source["taxa"].eq("Hauya_heydeana")].copy()
    if len(hit) != 1 or hit.iloc[0]["sc=0si=1"] != "1":
        raise ValueError("Hauya source state is not the pinned SI=1 row")
    source_row = int(hit.index[0]) + 2

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    cell = coverage.loc[
        coverage["accepted_species"].eq(SPECIES) & coverage["axis"].eq(AXIS)
    ]
    if len(cell) != 1 or cell.iloc[0]["quality"] != "":
        raise ValueError("Hauya reproductive cell is not currently unresolved")
    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    if (
        direct["accepted_species"].eq(SPECIES)
        & direct["trait_name"].eq(TRAIT)
    ).any():
        raise ValueError("Hauya self-incompatibility already exists in direct ledger")

    refs = (
        "Raven 1979, A survey of reproductive biology in Onagraceae, "
        f"doi:{RAVEN_DOI}; Wagner, Hoch & Raven 2007, Revised Classification of the "
        "Onagraceae, Systematic Botany Monographs 83; Freyman & Höhna 2019, "
        "doi:10.1093/sysbio/syy078"
    )
    packet = pd.DataFrame(
        [
            {
                "accepted_species": SPECIES,
                "axis": AXIS,
                "trait_name": TRAIT,
                "normalized_value": "SI",
                "quality": "medium",
                "source_url": RAVEN_URL,
                "source_lineage": f"doi:{RAVEN_DOI}",
                "source_reference_raw": refs,
                "source_column": "reproductive_biology_and_binary_state",
                "source_raw_value": (
                    "Raven: both species of Hauya self-incompatible; "
                    "Freyman matrix: Hauya_heydeana sc=0si=1 -> 1"
                ),
                "name_match_method": "exact_frozen_analysis_species_and_named_taxon",
                "review_status": "source_methodology_reviewed_reference_backed_strict_direct",
                "promotion_allowed": "false",
                "genus_rule_training_allowed": "false",
                "source_excel_row": f"freyman-hohna-row:{source_row}",
                "source_article_url": RAVEN_URL,
            }
        ]
    )
    packet.to_csv(
        args.output / "hauya_unresolved_si_batch.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    summary = {
        "contract": "hauya_si_reviewed_source_packet_v1",
        "accepted_species": SPECIES,
        "source_state": "SI",
        "freyman_source_commit": FREYMAN_COMMIT,
        "freyman_source_sha256": observed_sha,
        "freyman_source_row": source_row,
        "raven_1979_doi": RAVEN_DOI,
        "raven_support": "both species of Hauya are self-incompatible",
        "wagner_hoch_raven_2007_confirmation": "both species self-incompatible",
        "monograph_url": MONOGRAPH_URL,
        "current_cell_unresolved": True,
        "formal_gain": 0,
        "promotion_allowed": False,
        "genus_rule_training_allowed": False,
        "integration_deferred_to_source_batch": True,
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
