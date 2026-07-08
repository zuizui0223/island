import importlib

import pandas as pd
import pytest
import typer

guarded = importlib.import_module("island_v2.trait_source_discovery_guarded")


def _staged(island_id: str = "island_core") -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": "Example species",
                "island_id": island_id,
                "genus": "Example",
                "family": "Exampleaceae",
                "trait_candidate_status": "staged_not_curated",
                "release_gate": "accepted taxon scope + establishment review",
            }
        ]
    )


def test_guard_requires_all_staged_rows_to_match_requested_island():
    result = guarded.validate_staged_taxa_context(_staged(), "island_core")

    assert result.iloc[0]["island_id"] == "island_core"

    with pytest.raises(typer.BadParameter, match="exactly island_id"):
        guarded.validate_staged_taxa_context(_staged("wrong_island"), "island_core")


def test_guard_requires_curation_gate_columns_and_rejects_finalized_values():
    missing_gate = _staged().drop(columns=["release_gate"])
    with pytest.raises(typer.BadParameter, match="release_gate"):
        guarded.validate_staged_taxa_context(missing_gate, "island_core")

    finalized = _staged().assign(trait_value="white")
    with pytest.raises(typer.BadParameter, match="prohibited"):
        guarded.validate_staged_taxa_context(finalized, "island_core")


def test_release_gate_is_attached_to_every_unreviewed_source_lead():
    leads = pd.DataFrame({"query_taxon": ["Example species", "Unmatched species"]})
    result = guarded.attach_release_gate(leads, _staged())

    assert result.loc[0, "release_gate"] == "accepted taxon scope + establishment review"
    assert result.loc[1, "release_gate"] == ""


def _pilot() -> pd.DataFrame:
    # As written by scripts/build_validation_pilot.py: identity + seeded-unknown axes.
    return pd.DataFrame(
        [
            {"accepted_species": "Acampe ochracea", "family": "Orchidaceae",
             "native_status": "unknown", "growth_form": "unknown",
             "flower_conspicuousness": "unknown", "n_islands": "2", "n_records": "2"},
            {"accepted_species": "Antigonon leptopus", "family": "Polygonaceae",
             "native_status": "unknown", "growth_form": "unknown",
             "flower_conspicuousness": "unknown", "n_islands": "209", "n_records": "3836"},
            {"accepted_species": "Acampe ochracea", "family": "Orchidaceae",
             "native_status": "unknown", "growth_form": "unknown",
             "flower_conspicuousness": "unknown", "n_islands": "2", "n_records": "2"},
        ]
    )


def test_prepare_pilot_taxa_drops_prohibited_axes_and_dedups():
    taxa = guarded.prepare_pilot_taxa(_pilot())
    # Placeholder trait axes (incl. the prohibited native_status) are gone.
    assert set(taxa.columns) == {"accepted_species", "family"}
    assert not guarded.discovery.PROHIBITED_OUTPUT_COLUMNS.intersection(taxa.columns)
    # Duplicate species collapsed.
    assert list(taxa["accepted_species"]) == ["Acampe ochracea", "Antigonon leptopus"]


def test_prepare_pilot_taxa_requires_accepted_species():
    with pytest.raises(typer.BadParameter, match="accepted_species"):
        guarded.prepare_pilot_taxa(pd.DataFrame({"family": ["Orchidaceae"]}))


def test_prepare_pilot_taxa_output_passes_discovery_output_guard():
    taxa = guarded.prepare_pilot_taxa(_pilot())
    # discover_sources guards its input; the prepared frame must clear that guard.
    guarded.discovery._guard_output(taxa)  # must not raise
