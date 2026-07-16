import json

import pandas as pd

from island_v2 import acquisition_rate_report as report


def _write_master(path, species: list[str]) -> None:
    pd.DataFrame({"accepted_species": species, "genus": [s.split()[0] for s in species]}).to_csv(
        path, index=False
    )


def _write_candidate_index(path, rows: list[dict]) -> None:
    columns = [
        "accepted_species",
        "trait_name",
        "candidate_value",
        "source_url",
        "source_citation",
        "source_excerpt",
        "source_type",
        "evidence_scope",
        "confidence",
        "campaign_task",
        "campaign_phase",
    ]
    frame = pd.DataFrame(rows)
    for column in columns:
        if column not in frame.columns:
            frame[column] = "x"
    frame[columns].to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def _config(tmp_path, master_path) -> dict:
    return {
        "master_taxa_csv": str(master_path),
        "master_species_column": "accepted_species",
        "core_traits": {
            "colour": ["flower_primary_color"],
            "form": ["floral_form", "floral_symmetry"],
        },
        "reported_traits": ["flower_primary_color", "floral_form", "sex_system"],
        "sources": {
            "machine": {
                "adapter": "candidate_index",
                "path": str(tmp_path / "machine_candidate_index.csv.gz"),
                "counts_as_broad": True,
                "track": "machine",
            },
            "llm": {
                "adapter": "search_enabled_llm",
                "glob": str(tmp_path / "batch_*.jsonl"),
                "counts_as_broad": False,
                "track": "llm",
            },
        },
    }


def test_broad_rate_counts_source_backed_and_excludes_inference(tmp_path):
    master = ["Aaa bbb", "Ccc ddd", "Eee fff", "Ggg hhh"]
    _write_master(tmp_path / "master.csv", master)
    _write_candidate_index(
        tmp_path / "machine_candidate_index.csv.gz",
        [
            {
                "accepted_species": "Aaa bbb",
                "trait_name": "flower_primary_color",
                "candidate_value": "white",
            },
            {"accepted_species": "Aaa bbb", "trait_name": "floral_form", "candidate_value": "open"},
            # duplicate species-trait with a second value must not double count
            {
                "accepted_species": "Aaa bbb",
                "trait_name": "flower_primary_color",
                "candidate_value": "cream",
            },
            {
                "accepted_species": "Ccc ddd",
                "trait_name": "flower_primary_color",
                "candidate_value": "red",
            },
            # empty value must not count as acquired
            {
                "accepted_species": "Eee fff",
                "trait_name": "flower_primary_color",
                "candidate_value": "unknown",
            },
            # species outside the master must be dropped
            {
                "accepted_species": "Zzz zzz",
                "trait_name": "flower_primary_color",
                "candidate_value": "blue",
            },
        ],
    )
    llm_result = "Ggg hhh,yellow,star-shaped,bees,notes,mixed_mating,likely_SI,inference,low"
    (tmp_path / "batch_0.jsonl").write_text(
        json.dumps({"result": llm_result}) + "\n", encoding="utf-8"
    )

    config = _config(tmp_path, tmp_path / "master.csv")
    master_list = report.load_master_species(tmp_path / "master.csv", "accepted_species")
    assert master_list == master

    contributions, off_master = report.collect_contributions(config, set(master_list))
    assert off_master["machine"] == 1  # Zzz zzz dropped

    summary = report.build_summary(contributions, config, len(master_list))
    # colour: Aaa bbb + Ccc ddd (Eee fff 'unknown' excluded, Ggg hhh is inference-only)
    assert summary["n_colour_broad"] == 2
    # core minimum needs colour AND form: only Aaa bbb
    assert summary["n_core_minimum_broad"] == 1
    # any-trait broad excludes the inference-only LLM species
    assert summary["n_any_trait_broad"] == 2
    # any-channel includes the LLM species
    assert summary["n_any_trait_any_channel"] == 3

    by_trait = report.build_by_trait(contributions, config, len(master_list)).set_index(
        "trait_name"
    )
    assert by_trait.loc["flower_primary_color", "n_species_broad"] == 2
    # LLM colour for Ggg hhh shows only in the any-channel count
    assert by_trait.loc["flower_primary_color", "n_species_any_channel"] == 3

    by_track = report.build_by_track(contributions, len(master_list)).set_index("track")
    assert bool(by_track.loc["machine", "counts_as_broad"]) is True
    assert bool(by_track.loc["llm", "counts_as_broad"]) is False


def test_missing_source_files_contribute_nothing(tmp_path):
    master = ["Aaa bbb", "Ccc ddd"]
    _write_master(tmp_path / "master.csv", master)
    config = _config(tmp_path, tmp_path / "master.csv")
    # no candidate index or batch files written
    contributions, _ = report.collect_contributions(config, set(master))
    assert contributions.empty
    summary = report.build_summary(contributions, config, len(master))
    assert summary["n_any_trait_broad"] == 0
    assert summary["any_trait_broad_rate"] == 0.0


def test_master_wide_ledgers_keep_every_species_and_trait_cell(tmp_path):
    master = ["Aaa bbb", "Ccc ddd", "GenusOnly"]
    _write_master(tmp_path / "master.csv", master)
    _write_candidate_index(
        tmp_path / "machine_candidate_index.csv.gz",
        [
            {
                "accepted_species": "Aaa bbb",
                "trait_name": "flower_primary_color",
                "candidate_value": "white",
            }
        ],
    )
    config = _config(tmp_path, tmp_path / "master.csv")
    contributions, _ = report.collect_contributions(config, set(master))

    cell_ledger = report.build_species_trait_ledger(
        master,
        contributions,
        config["reported_traits"],
    )
    species_ledger = report.build_species_ledger(
        master,
        cell_ledger,
        config["reported_traits"],
    ).set_index("accepted_species")

    assert len(cell_ledger) == len(master) * len(config["reported_traits"])
    assert set(cell_ledger["accepted_species"]) == set(master)
    acquired = cell_ledger.loc[
        cell_ledger["accepted_species"].eq("Aaa bbb")
        & cell_ledger["trait_name"].eq("flower_primary_color")
    ].iloc[0]
    assert acquired["acquisition_status"] == "source_backed_candidate"
    assert acquired["unresolved_reason"] == ""
    missing = cell_ledger.loc[
        cell_ledger["accepted_species"].eq("Ccc ddd") & cell_ledger["trait_name"].eq("sex_system")
    ].iloc[0]
    assert missing["acquisition_status"] == "no_candidate_yet"
    assert missing["unresolved_reason"] == "not_found_in_configured_sources_yet"
    assert species_ledger.loc["Aaa bbb", "acquisition_status"] == "partial_source_backed"
    assert species_ledger.loc["GenusOnly", "acquisition_status"] == "no_candidate_yet"


def test_candidate_long_reads_bulk_standardized_value(tmp_path):
    path = tmp_path / "trait_candidates.csv.gz"
    pd.DataFrame(
        [
            {
                "accepted_species": "Aaa bbb",
                "trait_name": "flower_primary_color",
                "standardized_value": "white",
            }
        ]
    ).to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})

    pairs = report.adapter_candidate_long({"glob": str(path)})

    assert pairs.to_dict("records") == [
        {"accepted_species": "Aaa bbb", "trait_name": "flower_primary_color"}
    ]


def test_taxonomically_inapplicable_species_remain_in_both_ledgers():
    master = ["Aaa bbb", "Fern ccc"]
    traits = ["flower_primary_color", "sex_system"]
    applicability = pd.DataFrame(
        [
            {
                "accepted_species": "Aaa bbb",
                "taxonomic_group": "angiosperm",
                "angiosperm_analysis_eligible": True,
            },
            {
                "accepted_species": "Fern ccc",
                "taxonomic_group": "non_seed_vascular",
                "angiosperm_analysis_eligible": False,
            },
        ]
    )
    empty = pd.DataFrame(columns=report.CONTRIBUTION_COLUMNS)

    cells = report.build_species_trait_ledger(
        master,
        empty,
        traits,
        applicability,
        set(traits),
    )
    species = report.build_species_ledger(master, cells, traits).set_index("accepted_species")

    fern_cells = cells.loc[cells["accepted_species"].eq("Fern ccc")]
    assert set(fern_cells["acquisition_status"]) == {"not_applicable_taxonomic_scope"}
    assert species.loc["Fern ccc", "acquisition_status"] == ("not_applicable_taxonomic_scope")
    assert species.loc["Fern ccc", "n_applicable_traits"] == 0
    assert species.loc["Aaa bbb", "n_applicable_traits"] == 2
