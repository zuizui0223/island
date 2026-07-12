"""Production patch for complete GBIF Bombus species facet discovery.

GBIF occurrence search uses the camelCase ``facetLimit`` query parameter.
The earlier snake_case spelling was ignored by the API, silently returning the
default ten facet values.  This module patches only species discovery and
reuses the existing acquisition, environmental extraction, and CLI contract.
"""

from __future__ import annotations

import pandas as pd
import typer

from island_v2 import bombus_niche_inputs as base


def discover_bombus_species(
    min_occurrence_records: int,
    max_species: int,
) -> pd.DataFrame:
    """Discover all eligible Bombus species using GBIF's documented facetLimit."""
    requested_facets = max(max_species * 3, 500)
    response = base.httpx.get(
        f"{base.GBIF}/occurrence/search",
        params={
            "taxon_key": base._gbif_genus_key(),
            "has_coordinate": "true",
            "occurrence_status": "PRESENT",
            "limit": 0,
            "facet": "speciesKey",
            "facetLimit": requested_facets,
        },
        timeout=120.0,
    )
    response.raise_for_status()
    facets = response.json().get("facets", [])
    counts = facets[0].get("counts", []) if facets else []
    if len(counts) <= 10 and requested_facets > 10:
        raise typer.BadParameter(
            "GBIF returned at most ten Bombus species facets despite requesting "
            f"{requested_facets}; refusing to treat a truncated registry as global"
        )

    rows: list[dict[str, object]] = []
    for item in counts:
        count = int(item.get("count", 0))
        if count < min_occurrence_records:
            continue
        species_key = int(item["name"])
        species_response = base.httpx.get(
            f"{base.GBIF}/species/{species_key}", timeout=60.0
        )
        species_response.raise_for_status()
        payload = species_response.json()
        scientific_name = str(
            payload.get("canonicalName") or payload.get("scientificName") or ""
        ).strip()
        rank = str(payload.get("rank") or "").upper()
        if rank != "SPECIES" or not scientific_name.startswith("Bombus "):
            continue
        rows.append(
            {
                "species_key": species_key,
                "bombus_species": scientific_name,
                "gbif_coordinate_record_count": count,
            }
        )

    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return (
        result.drop_duplicates("species_key")
        .sort_values(
            ["gbif_coordinate_record_count", "bombus_species"],
            ascending=[False, True],
        )
        .head(max_species)
        .reset_index(drop=True)
    )


# The existing acquire_global command resolves this name from its defining
# module at runtime, so patching the module global keeps the CLI stable.
base.discover_bombus_species = discover_bombus_species
app = base.app


if __name__ == "__main__":
    app()
