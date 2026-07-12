import pandas as pd

from island_v2 import angiosperm_scope as scope


def _config() -> dict:
    return {
        "family_column": "family",
        "eligible_group": "angiosperm",
        "unresolved_group": "family_unresolved",
        "non_angiosperm_families": {
            "gymnosperm": ["Pinaceae", "Cycadaceae"],
            "non_seed_vascular": ["Dryopteridaceae", "Polypodiaceae"],
            "fossil_non_angiosperm": ["Lyginopteridaceae"],
        },
    }


def test_classify_scope_groups_and_eligibility():
    master = pd.DataFrame(
        {
            "accepted_species": [
                "Rosa canina",        # angiosperm family -> eligible
                "Pinus sylvestris",   # gymnosperm -> out
                "Dryopteris filix",   # fern -> out
                "Lyginopteris oldhamia",  # fossil -> out
                "Mystery plant",      # empty family -> unresolved
                "Cardiopteris moluccana",  # fern-sounding but angiosperm family -> eligible
            ],
            "family": ["Rosaceae", "Pinaceae", "Dryopteridaceae", "Lyginopteridaceae", "", "Cardiopteridaceae"],
        }
    )
    result = scope.classify_scope(master, _config()).set_index("accepted_species")

    assert result.loc["Rosa canina", "taxonomic_group"] == "angiosperm"
    assert bool(result.loc["Rosa canina", "angiosperm_analysis_eligible"]) is True
    assert result.loc["Pinus sylvestris", "taxonomic_group"] == "gymnosperm"
    assert bool(result.loc["Pinus sylvestris", "angiosperm_analysis_eligible"]) is False
    assert result.loc["Dryopteris filix", "taxonomic_group"] == "non_seed_vascular"
    assert result.loc["Lyginopteris oldhamia", "taxonomic_group"] == "fossil_non_angiosperm"
    assert result.loc["Mystery plant", "taxonomic_group"] == "family_unresolved"
    assert bool(result.loc["Mystery plant", "angiosperm_analysis_eligible"]) is False
    # fern-sounding family that is actually an angiosperm stays eligible
    assert bool(result.loc["Cardiopteris moluccana", "angiosperm_analysis_eligible"]) is True


def test_scope_summary_counts():
    master = pd.DataFrame(
        {
            "accepted_species": ["A one", "B two", "C three"],
            "family": ["Rosaceae", "Pinaceae", "Polypodiaceae"],
        }
    )
    summary = scope.scope_summary(scope.classify_scope(master, _config()))
    assert summary["n_master_species"] == 3
    assert summary["n_angiosperm_eligible"] == 1
    assert summary["n_out_of_scope"] == 2
    assert summary["by_taxonomic_group"]["gymnosperm"] == 1
    assert summary["by_taxonomic_group"]["non_seed_vascular"] == 1
