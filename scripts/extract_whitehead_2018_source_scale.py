"""Build the remaining strict Whitehead et al. (2018) mating-system packet.

Whitehead et al. provide a population-level tm table plus a source-defined
species summary table.  This repository already reviewed and hash-pinned that
supplement on 2026-08-11, accepted 23 exact-name species as high-quality direct
mating-system evidence, and fixed the mapping contract:

    mean_tm <= 0.2 -> predominantly_selfing
    0.2 < mean_tm < 0.8 -> mixed_mating
    mean_tm >= 0.8 -> predominantly_outcrossing

This adapter reuses that existing contract only for exact fixed-universe species
that remain unresolved today.  It does not infer from genus/family, does not
recompute a new species aggregation from population rows, and does not let the
new packet train genus rules.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

AXIS = "reproductive_assurance"
TRAIT = "mating_system"
TARGET = "Xanthorrhoea johnsonii"
ARTICLE_DOI = "10.3389/fevo.2018.00038"
DATA_DOI = "10.5061/dryad.3c340"
ARTICLE_URL = f"https://doi.org/{ARTICLE_DOI}"
DATA_URL = f"https://doi.org/{DATA_DOI}"
EXPECTED_ZIP_SHA256 = "7067064bf2bc37ade116b33dd3dc710f2938a36fc069265b0628499f887658f8"
EXPECTED_MEAN_TM = 0.95125
EXPECTED_POPS = 4
REVIEW_STATUS = "source_methodology_reviewed_reference_backed_strict_direct"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def classify(mean_tm: float) -> str:
    if not 0 <= mean_tm <= 1:
        raise ValueError(f"mean_tm outside [0,1]: {mean_tm}")
    if mean_tm <= 0.2:
        return "predominantly_selfing"
    if mean_tm < 0.8:
        return "mixed_mating"
    return "predominantly_outcrossing"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--source-zip", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--direct-ledger", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    observed_hash = sha256(args.source_zip)
    if observed_hash != EXPECTED_ZIP_SHA256:
        raise ValueError(f"Whitehead supplement hash mismatch: {observed_hash}")

    with ZipFile(args.source_zip) as archive:
        with archive.open("SpeciesTmData_Whitehead_etal.csv") as handle:
            species = pd.read_csv(handle, dtype=str).fillna("")
        with archive.open("PopTmData_Whitehead_etal.csv") as handle:
            populations = pd.read_csv(handle, dtype=str).fillna("")

    required_species = {"species", "mean_tm", "pops"}
    required_pop = {"species_name", "first_author", "year_published"}
    if not required_species.issubset(species.columns):
        raise ValueError(f"Whitehead species table missing {sorted(required_species - set(species.columns))}")
    if not required_pop.issubset(populations.columns):
        raise ValueError(f"Whitehead population table missing {sorted(required_pop - set(populations.columns))}")

    species = species.copy()
    species["accepted_species_candidate"] = species["species"].str.replace("_", " ", regex=False).map(text)
    species["mean_tm_numeric"] = pd.to_numeric(species["mean_tm"], errors="coerce")
    species["pops_numeric"] = pd.to_numeric(species["pops"], errors="coerce")

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected fixed 106295-species reproductive coverage, got {len(rep)}")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    universe = set(quality)
    unresolved = {name for name, q in quality.items() if not q}

    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    existing_target = direct.loc[
        direct["accepted_species"].eq(TARGET)
        & direct["trait_name"].isin(["mating_system", "self_incompatibility"])
    ]
    if len(existing_target):
        raise ValueError(f"{TARGET} already has direct reproductive evidence")

    eligible = species.loc[
        species["accepted_species_candidate"].isin(universe)
        & species["accepted_species_candidate"].isin(unresolved)
        & species["mean_tm_numeric"].notna()
        & species["pops_numeric"].ge(3)
    ].copy()

    if set(eligible["accepted_species_candidate"]) != {TARGET}:
        raise ValueError(
            "expected exactly the current Whitehead residual target; found "
            + repr(sorted(eligible["accepted_species_candidate"].unique()))
        )
    target = eligible.loc[eligible["accepted_species_candidate"].eq(TARGET)]
    if len(target) != 1:
        raise ValueError(f"expected one source-defined species row for {TARGET}, got {len(target)}")
    row = target.iloc[0]
    mean_tm = float(row["mean_tm_numeric"])
    pops = int(float(row["pops_numeric"]))
    if abs(mean_tm - EXPECTED_MEAN_TM) > 1e-9 or pops != EXPECTED_POPS:
        raise ValueError(f"{TARGET} source summary changed: mean_tm={mean_tm}, pops={pops}")

    source_name = text(row["species"])
    contributing = populations.loc[populations["species_name"].map(text).eq(source_name)].copy()
    if len(contributing) != pops:
        raise ValueError(f"{TARGET} contributing population rows {len(contributing)} != source pops {pops}")
    reference_keys = sorted(
        {
            f"{text(r.first_author)}:{int(float(r.year_published))}"
            for r in contributing.itertuples()
            if text(r.first_author) and text(r.year_published)
        }
    )
    if not reference_keys:
        raise ValueError(f"{TARGET} has no retained underlying citation keys")
    source_lineage = "citation-set:" + hashlib.sha256(
        "|".join(reference_keys).encode("utf-8")
    ).hexdigest()[:20]
    state = classify(mean_tm)
    if state != "predominantly_outcrossing":
        raise ValueError(f"unexpected Whitehead state for {TARGET}: {state}")

    source_row_number = int(target.index[0] + 2)
    full_audit = species.copy()
    full_audit["in_fixed_universe_exact"] = full_audit["accepted_species_candidate"].isin(universe)
    full_audit["currently_unresolved_reproductive"] = full_audit["accepted_species_candidate"].isin(unresolved)
    full_audit.to_csv(
        args.output / "whitehead_2018_species_source_audit.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    contributing.to_csv(args.output / "whitehead_2018_target_population_rows.csv", index=False)

    refs = (
        "Whitehead MR, Lanfear R, Mitchell RJ, Karron JD (2018), "
        f"Plant Mating Systems Often Vary Widely Among Populations, DOI {ARTICLE_DOI}; "
        f"Dryad dataset DOI {DATA_DOI}; underlying citation keys " + "|".join(reference_keys)
    )
    packet = pd.DataFrame([
        {
            "accepted_species": TARGET,
            "source_species_name": source_name.replace("_", " "),
            "axis": AXIS,
            "source_reference_raw": refs,
            "source_url": DATA_URL,
            "source_article_url": ARTICLE_URL,
            "source_lineage": source_lineage,
            "name_match_method": "exact_accepted_species_in_fixed_universe",
            "unresolved_reproductive_before_source": True,
            "promotion_allowed": False,
            "genus_rule_training_allowed": False,
            "trait_name": TRAIT,
            "normalized_value": state,
            "quality": "high",
            "source_column": "SpeciesTmData_Whitehead_etal.csv:mean_tm",
            "source_raw_value": f"mean_tm={mean_tm:.6f}; pops={pops}; source_row={source_row_number}",
            "review_status": REVIEW_STATUS,
            "acceptance_basis": (
                "Reuses the repository's hash-pinned 2026-08-11 Whitehead contract. "
                "The original supplement provides a source-defined species mean_tm and pops count; "
                "the repository's existing reviewed mapping classifies mean_tm>=0.8 as "
                "predominantly_outcrossing. No new population aggregation, proxy mapping, or "
                "genus/family inference is introduced."
            ),
        }
    ])
    packet.to_csv(
        args.output / "whitehead_2018_unresolved_mating_batch.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    summary = {
        "contract": "whitehead_2018_hash_pinned_source_defined_species_mean_tm_v1",
        "source_zip_sha256": observed_hash,
        "article_doi": ARTICLE_DOI,
        "dataset_doi": DATA_DOI,
        "source_species_rows": int(len(species)),
        "source_population_rows": int(len(populations)),
        "exact_current_unresolved_species_with_mean_tm_pops_ge3": 1,
        "reviewed_strict_packet_rows": 1,
        "reviewed_strict_packet_species": [TARGET],
        "target_mean_tm": mean_tm,
        "target_pops": pops,
        "target_state": state,
        "target_source_lineage": source_lineage,
        "underlying_reference_keys": reference_keys,
        "mapping_contract": "mean_tm<=0.2 predominantly_selfing; 0.2<mean_tm<0.8 mixed_mating; mean_tm>=0.8 predominantly_outcrossing",
        "mapping_contract_origin": "goal_reproduction_checkpoint_manifest_20260811.json plus tests/test_goal_reproduction_checkpoint.py",
        "formal_gain": 0,
        "formal_gain_pending_batch_integration": 1,
        "promotion_allowed": False,
        "genus_rule_training_allowed": False,
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
