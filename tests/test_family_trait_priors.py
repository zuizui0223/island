from __future__ import annotations

import pandas as pd

from island_v2.family_trait_priors import generate_family_trait_priors


def test_family_priors_fill_only_reported_trait_gaps() -> None:
    species = pd.DataFrame(
        [
            {"accepted_species": "Daisy alpha", "family": "Asteraceae"},
            {"accepted_species": "Mint beta", "family": "Lamiaceae"},
        ]
    )
    existing = pd.DataFrame(
        [
            {
                "accepted_species": "Mint beta",
                "trait_name": "floral_form",
                "provisional_candidate_value": "bell_campanulate",
                "candidate_class": "reported",
            }
        ]
    )

    frame, report = generate_family_trait_priors(species, existing)

    got = {
        (row.accepted_species, row.trait_name, row.provisional_candidate_value)
        for row in frame.itertuples(index=False)
    }
    assert ("Daisy alpha", "floral_form", "composite_head") in got
    assert ("Mint beta", "floral_symmetry", "zygomorphic") in got
    assert ("Mint beta", "floral_form", "tubular") not in got
    assert set(frame["candidate_class"]) == {"proxy"}
    assert set(frame["inference_status"]) == {"likely"}
    assert set(frame["confidence"]) == {"low"}
    assert report["n_species_with_prior"] == 2
