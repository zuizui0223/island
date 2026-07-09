import pandas as pd
import pytest
import typer

from island_v2 import taxon_identity_b_lite as identity


def test_identity_audit_matches_exact_and_reports_name_change():
    species = pd.DataFrame({"accepted_species": ["Olda olda", "Missinga missinga"]})
    authority = pd.DataFrame(
        [
            {
                "scientificName": "Olda olda",
                "acceptedNameUsage": "Newa newa",
                "family": "Fabaceae",
                "taxonomicStatus": "synonym",
                "taxonID": "urn:lsid:test:1",
                "source_url": "https://example.org/Olda_olda",
            }
        ]
    )

    frame, report = identity.audit_identity(species, authority, "test_authority")

    olda = frame.set_index("submitted_name").loc["Olda olda"]
    assert olda["authority_accepted_name"] == "Newa newa"
    assert olda["authority_family"] == "Fabaceae"
    assert olda["match_status"] == "matched_exact_name"
    missing = frame.set_index("submitted_name").loc["Missinga missinga"]
    assert missing["match_status"] == "unmatched"
    assert report["n_matched_exact_name"] == 1
    assert report["n_unmatched"] == 1
    assert report["n_name_changes_suggested"] == 1


def test_identity_audit_accepts_name_as_minimal_authority_column():
    species = pd.DataFrame({"accepted_species": ["Aaa aaa"]})
    authority = pd.DataFrame({"name": ["Aaa aaa"]})

    frame, _ = identity.audit_identity(species, authority)

    assert frame.loc[0, "authority_accepted_name"] == "Aaa aaa"
    assert frame.loc[0, "match_status"] == "matched_exact_name"


def test_identity_audit_requires_species_and_authority_name_columns():
    with pytest.raises(typer.BadParameter, match="accepted_species"):
        identity.audit_identity(pd.DataFrame({"species": ["Aaa aaa"]}), pd.DataFrame({"name": ["Aaa aaa"]}))

    with pytest.raises(typer.BadParameter, match="authority CSV"):
        identity.audit_identity(pd.DataFrame({"accepted_species": ["Aaa aaa"]}), pd.DataFrame({"family": ["Fabaceae"]}))
