from __future__ import annotations

import pandas as pd

from island_v2.taxon_name_expansion import expand_names


def test_expand_names_keeps_aliases_attached_to_accepted_species() -> None:
    def getter(url: str, params: dict[str, object]) -> dict[str, object]:
        if url.endswith("/species/match"):
            return {
                "usageKey": 10,
                "acceptedUsageKey": 20,
                "canonicalName": "Plantus alba",
                "status": "ACCEPTED",
            }
        if url.endswith("/species/20/synonyms"):
            return {
                "results": [
                    {
                        "key": 30,
                        "scientificName": "Oldplantus albus Author",
                        "canonicalName": "Oldplantus albus",
                        "taxonomicStatus": "SYNONYM",
                    },
                    {
                        "key": 31,
                        "scientificName": "Oldplantus albus Author",
                        "canonicalName": "Oldplantus albus",
                        "taxonomicStatus": "SYNONYM",
                    },
                ]
            }
        raise AssertionError(url)

    frame, report = expand_names(
        pd.DataFrame({"accepted_species": ["Plantus alba"]}),
        getter,
        max_taxa=10,
        max_synonyms_per_species=50,
    )

    assert frame["accepted_species"].tolist() == ["Plantus alba", "Plantus alba"]
    assert frame["name_role"].tolist() == ["accepted_input", "gbif_synonym"]
    assert frame.iloc[1]["search_name"] == "Oldplantus albus Author"
    assert frame.iloc[1]["gbif_accepted_key"] == "20"
    assert report["n_species_with_synonyms"] == 1
    assert report["n_synonym_names"] == 1
