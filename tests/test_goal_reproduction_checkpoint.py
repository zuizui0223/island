from __future__ import annotations

import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence

ROOT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "goal_reproduction_checkpoint_20260811"
)


def _canonical_bytes(path: Path) -> bytes:
    payload = path.read_bytes()
    if path.suffix.lower() in {".csv", ".json", ".md", ".txt"}:
        payload = payload.replace(b"\r\n", b"\n")
    return payload


def test_checkpoint_is_hash_pinned_and_passes_the_formal_gate() -> None:
    manifest = json.loads(
        (ROOT / "goal_reproduction_checkpoint_manifest_20260811.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = pd.read_csv(
        ROOT / "combined_curated_evidence_20260811.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        ROOT / "combined_curated_manual_audit_20260811.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence,
        audit,
        master,
        Path("config/trait_ontology.yml"),
    )

    assert len(evidence) == len(audit) == len(accepted) == 248
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["source_lineage"].ne("").all()
    assert audit["decision"].eq("accept").all()
    assert audit["cultivar_contamination"].str.casefold().eq("false").all()
    assert manifest["incremental_reviewed_rows"] == 78
    assert manifest["policy"]["trait_specific_genus_join"] is True
    assert manifest["policy"]["minimum_genus_species"] == 3
    assert manifest["policy"]["family_inference"] is False
    assert manifest["policy"]["global_fallback"] is False
    assert manifest["policy"]["n2_formal_adoption"] is False
    assert manifest["policy"]["cross_trait_substitution"] is False

    for filename, metadata in manifest["files"].items():
        payload = _canonical_bytes(ROOT / filename)
        assert len(payload) == metadata["size_bytes"]
        assert hashlib.sha256(payload).hexdigest() == metadata["sha256"]


def test_increment_has_expected_sources_and_trait_specific_targets() -> None:
    evidence = pd.read_csv(
        ROOT / "incremental_reviewed_evidence_20260811.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        ROOT / "incremental_reviewed_audit_20260811.csv", dtype=str
    ).fillna("")

    assert len(evidence) == len(audit) == 78
    assert evidence["source_group"].value_counts().to_dict() == {
        "nparks_near_rule_morphology_checkpoint": 49,
        "whitehead_2018_population_mating_system": 23,
        "silene_primary_reproduction_checkpoint": 2,
        "goal_reproduction_near_rule_checkpoint": 1,
        "goal_reproduction_primary_mating_checkpoint": 3,
    }
    targets = evidence.loc[
        evidence["source_group"].eq("goal_reproduction_near_rule_checkpoint"),
        ["accepted_species", "trait_name", "normalized_value", "evidence_quality"],
    ]
    assert set(map(tuple, targets.to_numpy())) == {
        ("Magnolia salicifolia", "mating_system", "mixed_mating", "medium"),
    }
    primary_mating = evidence.loc[
        evidence["source_group"].eq("goal_reproduction_primary_mating_checkpoint"),
        ["accepted_species", "trait_name", "normalized_value", "evidence_quality"],
    ]
    assert set(map(tuple, primary_mating.to_numpy())) == {
        ("Ficus arpazusa", "mating_system", "predominantly_outcrossing", "high"),
        ("Ficus petiolaris", "mating_system", "predominantly_outcrossing", "high"),
        ("Dioscorea alata", "mating_system", "predominantly_outcrossing", "high"),
    }
    assert "Hieracium umbellatum" not in set(evidence["accepted_species"])
    assert not set(evidence["trait_name"]).intersection(
        {"reward_type", "pollen_vector_mode"}
    )


def test_whitehead_rows_are_backed_by_the_committed_source_table() -> None:
    evidence = pd.read_csv(
        ROOT / "incremental_reviewed_evidence_20260811.csv", dtype=str
    ).fillna("")
    whitehead = evidence.loc[
        evidence["source_group"].eq("whitehead_2018_population_mating_system")
    ].copy()

    with ZipFile(ROOT / "whitehead_2018_supplementary_data.zip") as archive:
        with archive.open("SpeciesTmData_Whitehead_etal.csv") as source_file:
            source = pd.read_csv(source_file, dtype=str).fillna("")
        with archive.open("PopTmData_Whitehead_etal.csv") as population_file:
            populations = pd.read_csv(population_file, dtype=str).fillna("")

    source_species = set(source["species"].str.replace("_", " ", regex=False))
    assert len(whitehead) == 23
    assert set(whitehead["accepted_species"]).issubset(source_species)
    assert whitehead["trait_name"].eq("mating_system").all()
    assert whitehead["source_excerpt"].str.contains("mean_tm").all()

    for record in whitehead.to_dict("records"):
        source_row_number = int(record["source_record_id"].rsplit(":", 1)[1])
        source_record = source.iloc[source_row_number - 2]
        source_name = source_record["species"].replace("_", " ")
        assert source_name == record["accepted_species"]
        assert abs(float(source_record["mean_tm"]) - float(record["raw_value"])) < 1e-6
        assert int(float(source_record["pops"])) >= 3

        population_name = source_record["species"]
        contributing = populations.loc[populations["species_name"].eq(population_name)]
        reference_keys = sorted(
            {
                f"{row.first_author}:{int(float(row.year_published))}"
                for row in contributing.itertuples()
                if row.first_author and row.year_published
            }
        )
        expected_lineage = "citation-set:" + hashlib.sha256(
            "|".join(reference_keys).encode("utf-8")
        ).hexdigest()[:20]
        assert record["source_lineage"] == expected_lineage
