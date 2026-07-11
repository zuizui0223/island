from __future__ import annotations

from island_v2.trait_recovery_benchmark import benchmark_lanes


def test_benchmark_reports_union_and_marginal_gain() -> None:
    species = ["A", "B", "C"]
    lanes = {
        "fixed": {("A", "flower_primary_color"), ("B", "floral_form")},
        "synonym": {("A", "flower_primary_color"), ("C", "flower_primary_color")},
        "image": {("B", "flower_primary_color"), ("C", "floral_form")},
    }

    coverage, marginal, summary = benchmark_lanes(
        species=species,
        lanes=lanes,
        order=["fixed", "synonym", "image"],
    )

    union_any = coverage[
        (coverage["lane"] == "__union_all__") & (coverage["trait_name"] == "__any_trait__")
    ].iloc[0]
    assert int(union_any["n_species_recovered"]) == 3
    assert marginal["new_species_trait_pairs"].tolist() == [2, 1, 2]
    assert summary["union_species_trait_pairs"] == 5
    assert summary["union_species_any_trait"] == 3
