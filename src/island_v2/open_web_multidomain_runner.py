"""Run the multi-domain pilot with the manuscript's strict axis contract.

Reward type and pollen-vector mode are extracted as separate traits but never
counted as floral-structural-complexity coverage.
"""

from island_v2 import open_web_multidomain as pilot

pilot.TRAIT_AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_form": "floral_structural_complexity",
    "floral_symmetry": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
    "reward_type": None,
    "pollen_vector_mode": None,
    "self_incompatibility": "reproductive_assurance",
    "autonomous_selfing_capacity": "reproductive_assurance",
    "mating_system": "reproductive_assurance",
    "cleistogamy": "reproductive_assurance",
}

if __name__ == "__main__":
    pilot.app()
