# Acquisition trigger: resume from completed seed run 29742199817.
import pandas as pd

from island_v2.three_pass_trait_campaign import build_tasks


def test_three_pass_routing_is_species_axis_quality_ordered():
    master = pd.DataFrame({"accepted_species": ["Alpha one", "Alpha two", "Beta one"]})
    evidence = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "flower_primary_color",
                "evidence_quality": "high",
                "evidence_scope": "species_direct",
            },
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
                "evidence_quality": "medium",
                "evidence_scope": "species_direct",
            },
            {
                "accepted_species": "Alpha two",
                "trait_name": "mating_system",
                "evidence_quality": "low",
                "evidence_scope": "genus_inference",
            },
            {
                "accepted_species": "Beta one",
                "trait_name": "flower_primary_color",
                "evidence_quality": "low",
                "evidence_scope": "family_inference",
            },
        ]
    )

    tasks, summary = build_tasks(master, evidence)

    # 3 species x 3 axes form the immutable denominator.
    assert summary["species_axis_tasks"] == 9
    assert summary["family_inference_allowed"] is False

    # Existing high evidence removes only that species-axis from pass 1.
    p1 = tasks.loc[tasks["pass"].eq(1)]
    assert not ((p1.accepted_species == "Alpha one") & (p1.axis == "flower_colour")).any()

    # Existing medium evidence is retained but may still be upgraded in pass 1;
    # it is not searched again in pass 2.
    assert ((p1.accepted_species == "Alpha one") & (p1.axis == "floral_structural_complexity")).any()
    p2 = tasks.loc[tasks["pass"].eq(2)]
    assert not ((p2.accepted_species == "Alpha one") & (p2.axis == "floral_structural_complexity")).any()

    # Existing genus evidence prevents another pass-3 inference attempt.
    p3 = tasks.loc[tasks["pass"].eq(3)]
    assert not ((p3.accepted_species == "Alpha two") & (p3.axis == "reproductive_assurance")).any()

    # Family inference is ignored, so the task remains unresolved and eligible.
    assert ((p3.accepted_species == "Beta one") & (p3.axis == "flower_colour")).any()
