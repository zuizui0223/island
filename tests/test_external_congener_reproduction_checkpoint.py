from __future__ import annotations

import pandas as pd

from island_v2 import external_congener_reproduction_checkpoint as checkpoint


def response(
    source_name: str,
    family: str,
    *,
    accepted_species: str | None = None,
    gbif_family: str | None = None,
) -> dict[str, object]:
    accepted = accepted_species or source_name
    genus, epithet = accepted.split()
    return {
        "source_name": source_name,
        "retrieved_at_utc": "2026-08-28T00:00:00Z",
        "wfo_status": 200,
        "wfo": {
            "parsedName": {"rank": "species", "canonical_form": source_name},
            "params": {
                "fuzzyNameParts": 0,
                "checkHomonyms": True,
                "checkRank": True,
                "classificationVersion": "2026-06",
            },
            "match": {
                "placement": f"Code/Plantae/Angiosperms/{family}/{genus}/{epithet}",
                "wfo_id": f"wfo-{genus}-{epithet}",
            },
        },
        "gbif_status": 200,
        "gbif": {
            "matchType": "EXACT",
            "rank": "SPECIES",
            "kingdom": "Plantae",
            "species": accepted,
            "family": gbif_family or family,
            "usageKey": 123,
        },
    }


def test_mapping_gate_excludes_target_species_and_family_conflict() -> None:
    records = [
        response("Alpha outside", "Testaceae"),
        response("Alpha target", "Testaceae"),
        response("Beta outside", "Testaceae", gbif_family="Otheraceae"),
    ]

    audit = checkpoint.build_mapping_audit(
        records=records,
        expected_families={
            "Alpha outside": "Testaceae",
            "Alpha target": "Testaceae",
            "Beta outside": "Testaceae",
        },
        target_species={"Alpha target"},
        provider="fixture",
    )

    assert dict(zip(audit["source_name"], audit["mapping_reason"], strict=True)) == {
        "Alpha outside": "accepted_strict_two_backbone",
        "Alpha target": "mapped_into_fixed_target_universe",
        "Beta outside": "family_conflict",
    }


def test_checkpoint_builds_external_support_without_target_direct_rows() -> None:
    master = pd.DataFrame(
        [{"accepted_species": "Alpha target", "family": "Testaceae"}]
    )
    bsdb = pd.DataFrame(
        [
            {
                "Infrasp": "",
                "BreedingSys": "SC",
                "ISI_value": "0.3",
                "ISItype": "ISIfruits",
                "bs.Source": "Study A",
                "tnrs_family": "Testaceae",
                "tnrs_Sciname": "Alpha outside",
                "tnrs_infrasp": "",
            },
            {
                "Infrasp": "",
                "BreedingSys": "SI",
                "ISI_value": "0.9",
                "ISItype": "ISIfruits",
                "bs.Source": "Study target",
                "tnrs_family": "Testaceae",
                "tnrs_Sciname": "Alpha target",
                "tnrs_infrasp": "",
            },
        ]
    )
    rodger = pd.DataFrame(
        [
            {
                "source_row": "2",
                "source.dataset": "fixture",
                "study": "Study B",
                "genus.species": "Beta_outside",
                "TPL.family": "Testaceae",
                "auto.fs.x": "0.2",
                "auto.spfr.x": "",
                "auto.spofl.x": "",
                "auto.spofr.x": "",
                "auto.spfl.x": "",
                "auto.spp.x": "",
            }
        ]
    )

    evidence, mapping, selection = checkpoint.build_checkpoint(
        bsdb_source=bsdb,
        bsdb_records=[response("Alpha outside", "Testaceae")],
        rodger_source=rodger,
        rodger_records=[response("Beta outside", "Testaceae")],
        master=master,
    )

    assert len(evidence) == 2
    assert set(evidence["accepted_species"]) == {"Alpha outside", "Beta outside"}
    assert set(evidence["trait_name"]) == {
        "self_incompatibility",
        "autonomous_selfing_capacity",
    }
    assert set(evidence["evidence_scope"]) == {"external_congener_species_direct"}
    assert set(evidence["name_match_method"]) == {
        "strict_wfo_gbif_two_backbone"
    }
    assert "Alpha target" not in set(evidence["accepted_species"])
    assert len(mapping) == 2
    assert set(selection["selection_reason"]) == {
        "selected",
        "not_strict_external_two_backbone",
    }
