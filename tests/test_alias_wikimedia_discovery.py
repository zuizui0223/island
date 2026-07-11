from __future__ import annotations

import pandas as pd

from island_v2.alias_wikimedia_discovery import discover_alias_pages


def test_alias_only_qid_is_verified_against_alias_p225() -> None:
    aliases = pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alba",
                "search_name": "Plantus alba",
                "name_role": "accepted_input",
                "gbif_canonical_name": "Plantus alba",
            },
            {
                "accepted_species": "Plantus alba",
                "search_name": "Oldplantus albus Author",
                "name_role": "gbif_synonym",
                "gbif_canonical_name": "Oldplantus albus",
            },
        ]
    )

    def getter(url: str, params: dict[str, object]) -> dict[str, object]:
        if "wikipedia.org" in url:
            titles = str(params["titles"]).split("|")
            pages = {}
            for index, title in enumerate(titles, start=1):
                if title == "Plantus alba":
                    pages[str(index)] = {
                        "title": title,
                        "pageprops": {"wikibase_item": "Q1"},
                    }
                elif title in {"Oldplantus albus Author", "Oldplantus albus"}:
                    pages[str(index)] = {
                        "title": title,
                        "pageprops": {"wikibase_item": "Q2"},
                    }
            return {"query": {"pages": pages}}
        if "wikidata.org" in url:
            return {
                "entities": {
                    "Q1": {
                        "claims": {
                            "P225": [
                                {"mainsnak": {"datavalue": {"value": "Plantus alba"}}}
                            ]
                        }
                    },
                    "Q2": {
                        "claims": {
                            "P225": [
                                {"mainsnak": {"datavalue": {"value": "Oldplantus albus"}}}
                            ]
                        }
                    },
                }
            }
        raise AssertionError(url)

    frame, report = discover_alias_pages(aliases, getter)

    alias_rows = frame.loc[frame["qid"].eq("Q2")]
    assert alias_rows["verification_status"].eq("p225_verified").all()
    assert alias_rows["is_alias_only_qid"].eq("true").all()
    assert report["n_species_with_alias_only_verified_qid"] == 1
    assert report["alias_only_species"] == ["Plantus alba"]
