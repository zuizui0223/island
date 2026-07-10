from pathlib import Path

import pandas as pd

from island_v2 import interaction_guild_machine as guild_machine, pollinator_evidence


def test_classify_partner_guild_prefers_bombus_over_generic_bees():
    assert guild_machine.classify_partner_guild("Bombus pascuorum")[0] == "bumblebees"
    assert guild_machine.classify_partner_guild("Apis mellifera")[0] == "other_bees"
    assert guild_machine.classify_partner_guild("Calliphoridae")[0] == "flies"
    assert guild_machine.classify_partner_guild("Cantharidae")[0] == "beetles_wasps_ants"
    assert guild_machine.classify_partner_guild("Unmapped taxon")[0] == ""


def test_build_machine_guild_candidates_splits_candidates_and_holdouts():
    globi = pd.DataFrame(
        [
            {
                "evidence_id": "g1",
                "query_taxon": "Plant A",
                "interaction_type": "pollinates",
                "interaction_claim_class": "pollination_claim",
                "partner_taxon_name": "Bombus terrestris",
                "study_citation": "Study A",
            },
            {
                "evidence_id": "g2",
                "query_taxon": "Plant B",
                "interaction_type": "pollinates",
                "interaction_claim_class": "pollination_claim",
                "partner_taxon_name": "",
            },
        ]
    )

    candidates, holdouts = guild_machine.build_machine_guild_candidates(globi)

    assert len(candidates) == 1
    assert len(holdouts) == 1
    assert candidates.loc[0, "accepted_species"] == "Plant A"
    assert candidates.loc[0, "pollinator_guild"] == "bumblebees"
    assert candidates.loc[0, "effectiveness_class"] == "unknown"
    assert holdouts.loc[0, "accepted_species"] == "Plant B"


def test_long_format_output_validates_against_pollinator_schema():
    candidates = pd.DataFrame(
        [
            {
                "accepted_species": "Plant A",
                "pollinator_guild": "flies",
                "evidence_type": "floral_visitor",
                "effectiveness_class": "unknown",
                "source_citation": "Study",
                "source_url": "https://example.test",
                "study_location": "",
                "event_date": "2026-07-10",
                "wild_or_cultivated": "unknown",
            }
        ]
    )

    long_format = guild_machine.to_long_pollinator_evidence(candidates)
    schema = pollinator_evidence.load_schema(Path("config/pollinator_evidence_schema.yml"))
    issues = pollinator_evidence.validate_long_format(long_format, schema)

    assert issues == []
    assert long_format.loc[0, "review_status"] == "unreviewed"
    assert long_format.loc[0, "season"] == "07"


def test_machine_guild_index_flags_multi_guild_signal():
    candidates = pd.DataFrame(
        [
            {"accepted_species": "Plant A", "pollinator_guild": "bumblebees"},
            {"accepted_species": "Plant A", "pollinator_guild": "flies"},
            {"accepted_species": "Plant B", "pollinator_guild": "other_bees"},
        ]
    )

    index = guild_machine.machine_guild_index(candidates)

    plant_a = index.loc[index["accepted_species"].eq("Plant A")].iloc[0]
    plant_b = index.loc[index["accepted_species"].eq("Plant B")].iloc[0]
    assert plant_a["machine_functional_replacement_signal"]
    assert not plant_b["machine_functional_replacement_signal"]
