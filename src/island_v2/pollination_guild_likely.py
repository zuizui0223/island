"""Low-confidence family-level priors for main pollination guild.

These priors are intended only for expanded-sensitivity analyses when species-level
reported pollination evidence is unavailable. They are never promoted to reported data.
"""

from __future__ import annotations

import pandas as pd

FAMILY_GUILD_PRIORS = {
    "Orchidaceae": "mixed",
    "Asteraceae": "mixed",
    "Rubiaceae": "mixed",
    "Fabaceae": "bees",
    "Poaceae": "wind",
    "Cyperaceae": "wind",
    "Myrtaceae": "mixed",
    "Euphorbiaceae": "mixed",
    "Rosaceae": "mixed",
    "Apocynaceae": "mixed",
    "Lamiaceae": "bees",
    "Melastomataceae": "bees",
    "Malvaceae": "bees",
    "Ericaceae": "bees",
    "Arecaceae": "mixed",
    "Acanthaceae": "bees",
    "Lauraceae": "mixed",
    "Primulaceae": "mixed",
    "Gesneriaceae": "mixed",
    "Phyllanthaceae": "mixed",
    "Araceae": "flies",
    "Ranunculaceae": "bees",
    "Annonaceae": "mixed",
    "Brassicaceae": "bees",
    "Sapindaceae": "mixed",
}


def generate_pollination_guild_priors(species_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for record in species_table.fillna("").to_dict("records"):
        species = str(record.get("accepted_species", "")).strip()
        family = str(record.get("family", "")).strip()
        guild = FAMILY_GUILD_PRIORS.get(family)
        if not species or not guild:
            continue
        rows.append(
            {
                "accepted_species": species,
                "trait_name": "pollination_functional_guild",
                "provisional_candidate_value": {
                    "bees": "other_bees",
                    "flies": "flies",
                    "mixed": "mixed_or_generalist",
                    "wind": "abiotic_wind",
                }[guild],
                "candidate_class": "proxy",
                "inference_status": "likely",
                "confidence": "low",
                "evidence_quote": f"family={family}; family-level dominant pollination functional prior",
                "inference_rule": "declared_family_pollination_guild_prior",
                "searched_name": species,
                "search_query": "",
                "search_engine": "",
                "search_rank": "",
                "source_url": "",
                "source_title": "",
                "source_type": "declared_family_pollination_prior",
                "evidence_scope": "proxy_not_reported",
                "review_status": "unreviewed_proxy",
            }
        )
    return pd.DataFrame(rows)
