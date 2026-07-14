import numpy as np
import pandas as pd

from island_v2.bombus_niche_hypervolume import score_niche_hypervolumes
from island_v2.bombus_occurrence_quality import quality_filter_and_thin_occurrences


def test_occurrence_quality_filter_and_spatial_thinning():
    d = pd.DataFrame({
        "bombus_species": ["Bombus test"] * 5,
        "decimal_latitude": [35.01, 35.02, 36.0, 37.0, 38.0],
        "decimal_longitude": [139.01, 139.02, 140.0, 141.0, 142.0],
        "year": [2020, 2021, 2020, 1900, 2020],
        "coordinate_uncertainty_m": [100, 200, 100, 100, 50000],
        "basis_of_record": ["HUMAN_OBSERVATION"] * 5,
    })
    out = quality_filter_and_thin_occurrences(
        d, grid_degrees=0.25, max_coordinate_uncertainty_m=20000, min_year=1950
    )
    # first two records occupy the same 0.25-degree cell; old and imprecise records drop
    assert len(out) == 2
    assert out["year"].astype(int).min() >= 1950


def test_correlated_predictors_are_reduced_before_ellipsoidal_fit():
    rng = np.random.default_rng(42)
    n = 100
    x = rng.normal(size=n)
    occ = pd.DataFrame({
        "bombus_species": ["Bombus test"] * n,
        "bio1": x,
        "bio5": x * 2 + rng.normal(scale=0.01, size=n),
        "bio12": rng.normal(size=n),
    })
    target = pd.DataFrame({
        "island_id": ["i1"],
        "bombus_species": ["Bombus test"],
        "bio1": [0.0],
        "bio5": [0.0],
        "bio12": [0.0],
    })
    result = score_niche_hypervolumes(
        occ, target, ["bio1", "bio5", "bio12"], min_occurrences=50
    )
    row = result.iloc[0]
    assert row["model_status"] == "scored"
    assert row["model_type"] == "regularized_ellipsoidal_environmental_niche"
    assert int(row["n_pca_dimensions"]) < 3
    assert 0 <= float(row["environmental_compatibility"]) <= 1


def test_score_is_environmental_compatibility_not_absence_probability():
    rng = np.random.default_rng(7)
    occ = pd.DataFrame({
        "bombus_species": ["Bombus test"] * 60,
        "bio1": rng.normal(size=60),
        "bio12": rng.normal(size=60),
    })
    target = pd.DataFrame({
        "island_id": ["near", "far"],
        "bombus_species": ["Bombus test", "Bombus test"],
        "bio1": [0.0, 8.0],
        "bio12": [0.0, 8.0],
    })
    result = score_niche_hypervolumes(
        occ, target, ["bio1", "bio12"], min_occurrences=50
    ).set_index("island_id")
    assert result.loc["near", "environmental_compatibility"] > result.loc["far", "environmental_compatibility"]
    assert set(result["model_type"]) == {"regularized_ellipsoidal_environmental_niche"}
