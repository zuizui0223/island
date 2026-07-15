from pathlib import Path

import pandas as pd

from island_v2.izu_public_data import build


def test_build_izu_public_data(tmp_path: Path) -> None:
    species = pd.DataFrame(
        [
            {"island_id": "gshhg_2.3.7_h_3f48d601bf60ade348ae", "species": "Alpha beta", "n_records": "2", "n_unique_gbif_ids": "2"},
            {"island_id": "gshhg_2.3.7_h_80bd057071c9043cebcc", "species": "Alpha beta", "n_records": "1", "n_unique_gbif_ids": "1"},
            {"island_id": "gshhg_2.3.7_h_80bd057071c9043cebcc", "species": "Gamma delta", "n_records": "3", "n_unique_gbif_ids": "3"},
            {"island_id": "outside", "species": "Ignore me", "n_records": "9", "n_unique_gbif_ids": "9"},
        ]
    )
    effort = pd.DataFrame(
        [
            {"island_id": "gshhg_2.3.7_h_3f48d601bf60ade348ae", "n_records": "2", "n_species": "1", "n_datasets": "1"},
            {"island_id": "gshhg_2.3.7_h_80bd057071c9043cebcc", "n_records": "4", "n_species": "2", "n_datasets": "2"},
        ]
    )
    species_path = tmp_path / "species.csv.gz"
    effort_path = tmp_path / "effort.csv"
    out = tmp_path / "out"
    species.to_csv(species_path, index=False, compression="gzip")
    effort.to_csv(effort_path, index=False)

    summary = build(species_path, effort_path, out)

    assert summary["n_island_species_pairs"] == 3
    assert summary["n_unique_species_labels"] == 2
    assert summary["n_total_records"] == 6
    assert (out / "izu_island_species.csv.gz").exists()
    incidence = pd.read_csv(out / "izu_species_incidence.csv")
    alpha = incidence.loc[incidence["species"] == "Alpha beta"].iloc[0]
    assert alpha["n_islands"] == 2
