"""Hand-computable checks for the floral-profile geometry, not biological results."""

import importlib

import numpy as np
import pytest
from scipy.spatial.distance import pdist, squareform


def geometry():
    return importlib.import_module("island_v2.floral_space_geometry")


def test_multistate_distance_averages_traits_without_cartesian_expansion():
    profiles = [
        (("red", "white"), ("tubular",)),
        (("white",), ("tubular",)),
        (("red",), ("open",)),
    ]
    distance = geometry().profile_distances(profiles)
    np.testing.assert_allclose(
        distance, [[0.0, 0.25, 0.75], [0.25, 0.0, 1.0], [0.75, 1.0, 0.0]]
    )


def test_state_order_and_repeated_labels_do_not_create_extra_profiles():
    profiles = [
        (("red", "white"), ("open", "tubular")),
        (("white", "red", "red"), ("tubular", "open")),
    ]
    np.testing.assert_array_equal(geometry().profile_distances(profiles), np.zeros((2, 2)))


@pytest.mark.parametrize(
    "profiles",
    [
        [],
        [((), ("tubular",))],
        [(("red",), ("tubular",)), (("white",),)],
        [("red", ("tubular",))],
    ],
)
def test_distance_rejects_missing_or_misaligned_traits(profiles):
    with pytest.raises(ValueError):
        geometry().profile_distances(profiles)


def test_pcoa_reconstructs_euclidean_distances_and_reports_projection_loss():
    distance = np.array([[0.0, 3.0, 4.0], [3.0, 0.0, 5.0], [4.0, 5.0, 0.0]])
    result = geometry().pcoa(distance)
    assert result.coordinates.shape == (3, 2)
    assert result.lingoes_constant == 0.0
    assert result.raw_negative_inertia_fraction == 0.0
    np.testing.assert_allclose(squareform(pdist(result.coordinates)), distance, atol=1e-12)
    assert result.diagnostics[0]["raw_rmse"] > 0.1
    assert result.diagnostics[-1]["raw_rmse"] < 1e-12


def test_lingoes_corrects_negative_axes_but_retains_original_distance_error():
    # Four combinations of two binary nominal traits: D = mean mismatch.
    distance = np.array(
        [[0.0, 0.5, 0.5, 1.0], [0.5, 0.0, 1.0, 0.5],
         [0.5, 1.0, 0.0, 0.5], [1.0, 0.5, 0.5, 0.0]]
    )
    result = geometry().pcoa(distance)
    assert result.lingoes_constant == pytest.approx(0.25)
    assert result.raw_negative_inertia_fraction == pytest.approx(0.2)
    assert result.corrected_negative_inertia_fraction == 0.0
    expected = np.array(
        [[0.0, np.sqrt(0.75), np.sqrt(0.75), np.sqrt(1.5)],
         [np.sqrt(0.75), 0.0, np.sqrt(1.5), np.sqrt(0.75)],
         [np.sqrt(0.75), np.sqrt(1.5), 0.0, np.sqrt(0.75)],
         [np.sqrt(1.5), np.sqrt(0.75), np.sqrt(0.75), 0.0]]
    )
    np.testing.assert_allclose(result.corrected_distance, expected)
    np.testing.assert_allclose(squareform(pdist(result.coordinates)), expected, atol=1e-12)
    assert result.diagnostics[-1]["corrected_rmse"] < 1e-12
    assert result.diagnostics[-1]["raw_rmse"] > 0.1


def test_sqrt_transformation_has_its_own_reconstruction_target():
    distance = np.array(
        [[0.0, 0.5, 0.5, 1.0], [0.5, 0.0, 1.0, 0.5],
         [0.5, 1.0, 0.0, 0.5], [1.0, 0.5, 0.5, 0.0]]
    )
    result = geometry().pcoa(distance, correction="sqrt")
    assert result.correction == "sqrt"
    assert result.lingoes_constant == 0.0
    assert result.corrected_negative_inertia_fraction == 0.0
    np.testing.assert_allclose(
        squareform(pdist(result.coordinates)), np.sqrt(distance), atol=1e-12
    )
    assert result.diagnostics[-1]["raw_rmse"] > 0.1


def test_pcoa_refuses_duplicate_profiles_before_correction_can_separate_them():
    with pytest.raises(ValueError, match="unique profiles"):
        geometry().pcoa(np.zeros((2, 2)))


@pytest.mark.parametrize(
    "distance",
    [
        np.array([[0.0, 1.0], [0.5, 0.0]]),
        np.array([[0.0, np.nan], [np.nan, 0.0]]),
        np.array([[1.0, 1.0], [1.0, 0.0]]),
        np.array([[0.0, -1.0], [-1.0, 0.0]]),
        np.zeros((2, 3)),
    ],
)
def test_pcoa_rejects_invalid_distance_matrices(distance):
    with pytest.raises(ValueError):
        geometry().pcoa(distance)


def test_pcoa_single_profile_has_no_artificial_axis():
    result = geometry().pcoa(np.zeros((1, 1)))
    assert result.coordinates.shape == (1, 0)
    assert result.diagnostics == []


def test_profile_metrics_separate_number_of_types_from_species_per_type():
    metrics = geometry().profile_metrics(["a", "a", "a", "b", "c"])
    assert metrics["n_species"] == 5
    assert metrics["n_profiles"] == 3
    assert metrics["singleton_profiles"] == 2
    assert metrics["singleton_fraction"] == pytest.approx(2 / 3)
    assert metrics["singleton_species_fraction"] == pytest.approx(2 / 5)
    assert metrics["mean_redundancy"] == pytest.approx(5 / 3)
    assert metrics["effective_profiles_shannon"] == pytest.approx(2.586409289858933)
    assert metrics["effective_profiles_simpson"] == pytest.approx(25 / 11)
    assert metrics["largest_profile_share"] == pytest.approx(0.6)


def test_repeated_identical_species_increase_redundancy_only():
    metrics = geometry().profile_metrics([7] * 10)
    assert metrics["n_profiles"] == 1
    assert metrics["mean_redundancy"] == 10
    assert metrics["singleton_profiles"] == 0
    assert metrics["effective_profiles_shannon"] == 1


def test_empty_profile_metrics_do_not_report_observed_singleton_fraction_zero():
    metrics = geometry().profile_metrics([])
    assert metrics["n_species"] == 0
    assert metrics["n_profiles"] == 0
    assert np.isnan(metrics["singleton_fraction"])
    assert np.isnan(metrics["mean_redundancy"])


@pytest.mark.parametrize("ids", [[None], [np.nan]])
def test_missing_profile_ids_are_not_counted_as_biological_types(ids):
    with pytest.raises(ValueError):
        geometry().profile_metrics(ids)


def test_exact_rarefaction_uses_sampling_without_replacement():
    ids = ["a", "a", "a", "b", "c"]
    # For two draws, P(a observed)=9/10 and P(each singleton observed)=2/5.
    assert geometry().rarefied_profile_richness(ids, 2) == pytest.approx(1.7)
    assert geometry().rarefied_profile_richness(ids, 3) == pytest.approx(2.2)
    assert geometry().rarefied_profile_richness(ids, 0) == 0.0
    assert geometry().rarefied_profile_richness(ids, 5) == 3.0


@pytest.mark.parametrize("sample_size", [-1, 6, 1.5])
def test_rarefaction_never_silently_extrapolates(sample_size):
    with pytest.raises(ValueError):
        geometry().rarefied_profile_richness([0, 0, 0, 1, 2], sample_size)
