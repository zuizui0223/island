import pandas as pd

from island_v2 import trait_proxy_generator as proxy


def _rules():
    return {
        "floral_syndrome_proxy": {
            "wind_like_families": ["Poaceae"],
            "large_bee_or_Bombus_like": {
                "colours": ["blue", "purple", "white", "yellow"],
                "symmetries": ["zygomorphic"],
            },
            "butterfly_or_moth_like": {"colours": ["red", "orange", "pink"]},
            "open_or_generalist_insect_like": {
                "colours": ["white", "cream", "yellow", "green", "brown"],
                "symmetries": ["actinomorphic"],
            },
        },
        "compatibility_system_proxy": {
            "likely_self_compatible_families": ["Brassicaceae"],
            "likely_self_incompatible_families": [],
        },
        "reproductive_assurance_proxy": {
            "reproductive_assurance_like_families": ["Plantaginaceae"],
        },
    }


def test_large_bee_proxy_from_reported_phenotype_but_not_pollinator_identity():
    species = pd.DataFrame(
        {"accepted_species": ["Aaa aaa"], "family": ["Fabaceae"]}
    )
    reported = pd.DataFrame(
        [
            {
                "accepted_species": "Aaa aaa",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "purple",
            },
            {
                "accepted_species": "Aaa aaa",
                "trait_name": "floral_symmetry",
                "provisional_candidate_value": "zygomorphic",
            },
        ]
    )

    frame, report = proxy.generate_proxies(species, reported, _rules())

    assert report["n_proxy_candidates"] == 1
    row = frame.iloc[0]
    assert row["trait_name"] == "floral_syndrome_proxy"
    assert row["proxy_value"] == "large_bee_or_Bombus_like_floral_phenotype_proxy"
    assert row["candidate_class"] == "proxy"
    assert "flower_primary_color=purple" in row["basis_traits"]
    assert row["evidence_scope"] == "proxy_not_reported"


def test_reported_pollination_vector_blocks_floral_proxy():
    species = pd.DataFrame({"accepted_species": ["Beea beea"], "family": ["Fabaceae"]})
    reported = pd.DataFrame(
        [
            {
                "accepted_species": "Beea beea",
                "trait_name": "pollination_vector_reported",
                "provisional_candidate_value": "bee",
            },
            {
                "accepted_species": "Beea beea",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "purple",
            },
            {
                "accepted_species": "Beea beea",
                "trait_name": "floral_symmetry",
                "provisional_candidate_value": "zygomorphic",
            },
        ]
    )

    frame, _ = proxy.generate_proxies(species, reported, _rules())

    assert "floral_syndrome_proxy" not in set(frame["trait_name"])


def test_family_prior_can_generate_compatibility_proxy_only_when_no_reported_value():
    species = pd.DataFrame(
        {
            "accepted_species": ["Selfa selfa", "Knowna knowna"],
            "family": ["Brassicaceae", "Brassicaceae"],
        }
    )
    reported = pd.DataFrame(
        [
            {
                "accepted_species": "Knowna knowna",
                "trait_name": "self_compatibility_reported",
                "provisional_candidate_value": "self_compatible",
            }
        ]
    )

    frame, _ = proxy.generate_proxies(species, reported, _rules())
    compat = frame.loc[frame["trait_name"] == "compatibility_system_proxy"]

    assert list(compat["accepted_species"]) == ["Selfa selfa"]
    assert compat.iloc[0]["proxy_value"] == "likely_self_compatible_proxy"


def test_gbif_descriptive_table_is_normalised_as_reported(tmp_path):
    path = tmp_path / "m0.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Opena opena",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "yellow",
                "evidence_excerpt": "Flowers yellow.",
            }
        ]
    ).to_csv(path, index=False)

    loaded = proxy.load_reported([path])

    assert loaded.loc[0, "candidate_class"] == "reported"
    assert loaded.loc[0, "raw_description"] == "Flowers yellow."
