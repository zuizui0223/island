import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.goodwillie_2005_mixed_mating_checkpoint import (
    SOURCE_CSV_SHA256,
    SOURCE_README_SHA256,
    SOURCE_XLS_SHA256,
)
from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/goodwillie_2005_mixed_mating_checkpoint_20260821"
)
SOURCE_DIR = Path("data/v2/external/traits/goodwillie_2005")
SOURCE_GROUP = "goodwillie_2005_mixed_mating_checkpoint_20260821"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(CHECKPOINT / name, dtype=str).fillna("")


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _text_sha256(path: Path) -> str:
    canonical = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def test_goodwillie_source_snapshot_is_complete_and_immutable() -> None:
    assert (
        _file_sha256(SOURCE_DIR / "mixed_mating_outcrossing_rate_database.xls") == SOURCE_XLS_SHA256
    )
    assert _text_sha256(SOURCE_DIR / "README.txt") == SOURCE_README_SHA256
    normalized = SOURCE_DIR / "mixed_mating_outcrossing_rate_database_normalized.csv"
    assert _text_sha256(normalized) == SOURCE_CSV_SHA256
    source = pd.read_csv(normalized, dtype=str).fillna("")
    assert len(source) == 469
    assert not source["source_row_number"].duplicated().any()


def test_goodwillie_accepts_only_explicit_natural_population_modes() -> None:
    evidence = _read("goodwillie_2005_mode_of_selfing_evidence_20260821.csv")
    audit = _read("goodwillie_2005_mode_of_selfing_manual_audit_20260821.csv")

    assert len(evidence) == len(audit) == 9
    assert evidence["accepted_species"].nunique() == 8
    assert evidence["source_lineage"].nunique() == 8
    assert not evidence.duplicated(["accepted_species", "trait_name"]).any()
    assert set(evidence["evidence_scope"]) == {"species_direct"}
    assert set(evidence["evidence_quality"]) == {"high"}
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["wild_cultivated_cultivar_status"].eq("wild_natural_population").all()
    assert evidence["source_excerpt"].str.contains("PopType=0", regex=False).all()
    assert set(evidence["trait_name"]) == {
        "autonomous_selfing_capacity",
        "cleistogamy",
    }
    assert not evidence["trait_name"].eq("mating_system").any()
    autonomy = evidence["trait_name"].eq("autonomous_selfing_capacity")
    cleistogamy = evidence["trait_name"].eq("cleistogamy")
    assert autonomy.sum() == 7
    assert cleistogamy.sum() == 2
    assert evidence.loc[autonomy, "raw_value"].str.contains("auton", regex=False).all()
    assert evidence.loc[cleistogamy, "raw_value"].str.contains("cleis", regex=False).all()
    assert set(evidence.loc[autonomy, "normalized_value"]) == {"autonomous"}
    assert set(evidence.loc[cleistogamy, "normalized_value"]) == {"facultative"}
    assert audit["decision"].eq("accept").all()
    assert audit["species_identity_correct"].eq("true").all()
    assert audit["value_correct"].eq("true").all()
    assert audit["provenance_complete"].eq("true").all()
    assert audit["cultivar_contamination"].eq("false").all()


def test_goodwillie_lineage_is_underlying_study_not_dryad_provider() -> None:
    evidence = _read("goodwillie_2005_mode_of_selfing_evidence_20260821.csv")
    begonia = evidence.loc[evidence["accepted_species"].str.startswith("Begonia ")]
    assert len(begonia) == 2
    assert begonia["source_lineage"].nunique() == 1
    assert begonia["source_lineage"].iloc[0].startswith("citation:")
    assert (
        evidence["lineage_method"]
        .eq("normalized_underlying_study_citation_not_dryad_redistributor")
        .all()
    )


def test_goodwillie_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = _read("combined_curated_evidence_20260821.csv")
    audit = _read("combined_curated_manual_audit_20260821.csv")
    master = pd.read_csv("data/v2/staging/gbif/collected/island_taxa.csv", dtype=str).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_944
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 9


def test_goodwillie_manifest_records_dry_run_and_fail_closed_mapping() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_goodwillie_2005.json").read_text(
            encoding="utf-8"
        )
    )

    assert manifest["baseline_formal_run_id"] == 32445414491
    assert manifest["source_rows"] == 469
    assert manifest["exact_master_mating_system_rows"] == 243
    assert manifest["novel_mating_system_rows"] == 0
    assert manifest["accepted_evidence_rows"] == 9
    assert manifest["accepted_species"] == 8
    assert manifest["accepted_source_lineages"] == 8
    assert manifest["formal_search_api_queries"] == 0
    assert manifest["search_cost_usd"] == 0.0
    assert manifest["local_formal_simulation"]["net_change"] == 660
    assert manifest["local_formal_simulation"]["reproductive_assurance_after"] == 31058
    assert manifest["guardrails"]["none_mapped_to_absent_autonomous_selfing"] is False
    assert manifest["guardrails"]["geitonogamy_mapped_to_autonomous_selfing"] is False
    assert manifest["guardrails"]["genus_axis_only_join"] is False
