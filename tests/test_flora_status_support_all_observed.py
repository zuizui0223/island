import json
from pathlib import Path

import pandas as pd

from island_v2.flora_status_support import STRATA, run, stratum_mask


def test_all_observed_retains_every_flora_row_regardless_of_status():
    frame = pd.DataFrame(
        {
            "origin_status": ["native", "introduced", "unresolved", "native"],
            "floristic_status": [
                "native_nonendemic",
                "introduced",
                "unresolved",
                "endemic",
            ],
        }
    )
    mask = stratum_mask(frame, "all_observed")
    assert mask.dtype == bool
    assert mask.tolist() == [True, True, True, True]


def test_all_observed_is_explicitly_part_of_supported_strata():
    assert STRATA[0] == "all_observed"
    assert set(STRATA) >= {
        "all_observed",
        "all_native",
        "native_nonendemic",
        "endemic",
        "introduced",
        "unresolved",
    }


def test_manifest_keeps_v2_contract_and_declares_all_observed_extension(tmp_path: Path):
    island_species = tmp_path / "island_species.csv"
    status = tmp_path / "status.csv"
    evidence = tmp_path / "evidence.csv"
    covariates = tmp_path / "covariates.csv"
    config = tmp_path / "config.yml"
    output = tmp_path / "out"

    pd.DataFrame({"island_id": ["i1"], "accepted_species": ["Species one"]}).to_csv(
        island_species, index=False
    )
    pd.DataFrame(
        {
            "island_id": ["i1"],
            "accepted_species": ["Species one"],
            "origin_status": ["native"],
            "endemic_status": ["nonendemic"],
        }
    ).to_csv(status, index=False)
    pd.DataFrame(
        {
            "accepted_species": ["Species one"],
            "trait_name": ["trait_a"],
            "evidence_scope": ["species_direct"],
            "resolution_status": ["resolved"],
        }
    ).to_csv(evidence, index=False)
    pd.DataFrame({"island_id": ["i1"], "x": [1.0]}).to_csv(covariates, index=False)
    config.write_text(
        '{"outcomes":{"outcome_a":{"trait_name":"trait_a"}},'
        '"support":{"min_direct_species_per_island":1}}\n',
        encoding="utf-8",
    )

    run(
        island_species_csv=island_species,
        status_ledger_csv=status,
        direct_evidence_csv=evidence,
        covariates_csv=covariates,
        output_dir=output,
        config_path=config,
    )

    manifest = json.loads((output / "flora_status_support_manifest.json").read_text())
    assert manifest["contract"] == "flora_status_support_v2"
    assert manifest["extensions"] == ["all_observed_stratum_v1"]
    assert manifest["strata"][0] == "all_observed"
