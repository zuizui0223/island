"""Geometry and diversity of recorded, possibly multistate floral profiles.

One profile contains one state *set* per trait. A multistate record stays a
single profile; its Cartesian product is never interpreted as observed flowers.
All traits passed to ``profile_distances`` are nominal and receive equal weight.
An ordinal depth/size analysis needs an explicitly different distance contract.

PCoA is fitted once to unique profiles, then callers map species to those fixed
coordinates. In particular, Lingoes correction must not separate duplicate
profiles. Geometric diagnostics retain error against the original dissimilarity
as well as against any transformed distance.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Hashable, Iterable, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln


def profile_distances(profiles: Sequence[Sequence[Iterable[str]]]) -> np.ndarray:
    """Mean within-trait set-Jaccard dissimilarity between complete profiles.

    For single-state traits this is the usual equally weighted nominal Gower
    mismatch. Repeated states and state order have no effect. Missing traits are
    rejected: pairwise deletion would silently give different trait weights to
    different species pairs.
    """
    if not len(profiles) or not len(profiles[0]):
        raise ValueError("At least one complete profile and one trait are required")
    n_traits = len(profiles[0])
    normalized = []
    for profile in profiles:
        if len(profile) != n_traits:
            raise ValueError("Profiles must have the same ordered traits")
        row = []
        for states in profile:
            if isinstance(states, str):
                raise ValueError("Each trait must be a nonempty state collection, not a string")
            values = frozenset(states)
            if not values or any(not isinstance(value, str) or not value for value in values):
                raise ValueError("Missing traits and empty/non-string states are not distances")
            row.append(values)
        normalized.append(row)

    n_profiles = len(normalized)
    distance = np.zeros((n_profiles, n_profiles), dtype=float)
    for trait_index in range(n_traits):
        states = sorted(set().union(*(profile[trait_index] for profile in normalized)))
        indicator = np.array(
            [[state in profile[trait_index] for state in states] for profile in normalized],
            dtype=np.int64,
        )
        intersections = indicator @ indicator.T
        sizes = indicator.sum(axis=1)
        unions = sizes[:, None] + sizes[None, :] - intersections
        distance += 1.0 - intersections / unions
    distance /= n_traits
    np.fill_diagonal(distance, 0.0)
    return distance


@dataclass(frozen=True)
class PCoAResult:
    """A common coordinate system, with all fidelity diagnostics explicit.

    ``eigenvalues`` corresponds to the columns of ``coordinates`` and contains
    only positive axes. Raw/corrected eigenvalue vectors include negative and
    numerical-zero axes. Negative inertia is sum(abs(negative eigenvalues)) /
    sum(abs(all eigenvalues)), after a scale-relative numerical tolerance.

    ``raw_relative_rmse`` uses the mean original off-diagonal distance as its
    denominator. No variance percentage substitutes for original-distance error.
    """

    coordinates: np.ndarray
    eigenvalues: np.ndarray
    raw_eigenvalues: np.ndarray
    corrected_eigenvalues: np.ndarray
    correction: str
    raw_negative_inertia_fraction: float
    corrected_negative_inertia_fraction: float
    lingoes_constant: float
    corrected_distance: np.ndarray
    diagnostics: list[dict[str, float | int]]


def _distance_eigendecomposition(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    squared = np.square(distance)
    means = squared.mean(axis=1)
    gram = -0.5 * (squared - means[:, None] - means[None, :] + squared.mean())
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    eigenvalues = eigenvalues[::-1]
    eigenvectors = eigenvectors[:, ::-1]
    scale = float(np.max(np.abs(eigenvalues)))
    tolerance = 100 * np.finfo(float).eps * max(1, len(distance)) * scale
    eigenvalues[np.abs(eigenvalues) <= tolerance] = 0.0
    return eigenvalues, eigenvectors, tolerance


def _negative_inertia(eigenvalues: np.ndarray) -> float:
    absolute_sum = float(np.abs(eigenvalues).sum())
    if absolute_sum == 0:
        return 0.0
    return float(-eigenvalues[eigenvalues < 0].sum() / absolute_sum)


def pcoa(distance: np.ndarray, *, correction: str = "lingoes") -> PCoAResult:
    """Classical PCoA on unique profiles, optionally correcting non-Euclidean D.

    ``lingoes`` adds twice the magnitude of the most negative raw Gram
    eigenvalue to squared off-diagonal distances, before taking square roots.
    ``sqrt`` instead embeds sqrt(D); ``none`` embeds the positive raw axes.
    The latter two are not assumed Euclidean: their residual negative inertia
    is reported. All modes report projection error against untransformed D.

    Zero-distance duplicate profiles must be collapsed before fitting. A single
    profile returns zero axes; it does not receive an artificial dimension.
    """
    original = np.asarray(distance, dtype=float)
    if original.ndim != 2 or original.shape[0] != original.shape[1] or not len(original):
        raise ValueError("Distance must be a nonempty square matrix")
    if not np.isfinite(original).all() or np.any(original < 0):
        raise ValueError("Distances must be finite and nonnegative")
    if not np.allclose(original, original.T, rtol=1e-10, atol=1e-12):
        raise ValueError("Distance must be symmetric")
    if not np.allclose(np.diag(original), 0, rtol=0, atol=1e-12):
        raise ValueError("Distance must have a zero diagonal")
    pairs = np.triu_indices(len(original), k=1)
    if np.any(original[pairs] == 0):
        raise ValueError("Fit PCoA on unique profiles; map duplicate species to those coordinates")
    if correction not in {"lingoes", "sqrt", "none"}:
        raise ValueError("Correction must be 'lingoes', 'sqrt', or 'none'")
    original = (original + original.T) / 2
    np.fill_diagonal(original, 0.0)

    raw_values, raw_vectors, _ = _distance_eigendecomposition(original)
    lingoes_constant = 0.0
    target = original.copy()
    if correction == "lingoes" and raw_values[-1] < 0:
        lingoes_constant = float(-raw_values[-1])
        target = np.sqrt(np.square(original) + 2 * lingoes_constant)
        np.fill_diagonal(target, 0.0)
    elif correction == "sqrt":
        target = np.sqrt(original)

    if np.array_equal(original, target):
        target_values, target_vectors = raw_values, raw_vectors
    else:
        target_values, target_vectors, _ = _distance_eigendecomposition(target)
    positive = target_values > 0
    positive_values = target_values[positive]
    coordinates = target_vectors[:, positive] * np.sqrt(positive_values)
    # Fix arbitrary eigenvector signs for stable plots without asserting a
    # particular orientation inside tied-eigenvalue subspaces.
    for dimension in range(coordinates.shape[1]):
        anchor = int(np.argmax(np.abs(coordinates[:, dimension])))
        if coordinates[anchor, dimension] < 0:
            coordinates[:, dimension] *= -1

    diagnostics = []
    if len(pairs[0]) and len(positive_values):
        raw_pairs, target_pairs = original[pairs], target[pairs]
        mean_raw_distance = float(raw_pairs.mean())
        squared_pairs = np.zeros(len(raw_pairs))
        cumulative_inertia = np.cumsum(positive_values) / positive_values.sum()
        for axis in range(coordinates.shape[1]):
            squared_pairs += np.square(
                coordinates[pairs[0], axis] - coordinates[pairs[1], axis]
            )
            fitted_pairs = np.sqrt(squared_pairs)
            raw_error = fitted_pairs - raw_pairs
            corrected_error = fitted_pairs - target_pairs
            raw_rmse = float(np.sqrt(np.mean(np.square(raw_error))))
            diagnostics.append(
                {
                    "dimension": axis + 1,
                    "positive_inertia_fraction": float(cumulative_inertia[axis]),
                    "raw_mae": float(np.mean(np.abs(raw_error))),
                    "raw_rmse": raw_rmse,
                    "raw_relative_rmse": raw_rmse / mean_raw_distance,
                    "corrected_mae": float(np.mean(np.abs(corrected_error))),
                    "corrected_rmse": float(np.sqrt(np.mean(np.square(corrected_error)))),
                }
            )
    return PCoAResult(
        coordinates=coordinates,
        eigenvalues=positive_values,
        raw_eigenvalues=raw_values,
        corrected_eigenvalues=target_values,
        correction=correction,
        raw_negative_inertia_fraction=_negative_inertia(raw_values),
        corrected_negative_inertia_fraction=_negative_inertia(target_values),
        lingoes_constant=lingoes_constant,
        corrected_distance=target,
        diagnostics=diagnostics,
    )


def _profile_counts(profile_ids: Iterable[Hashable]) -> np.ndarray:
    ids = list(profile_ids)
    if any(
        value is None or (isinstance(value, (float, np.floating)) and np.isnan(value))
        for value in ids
    ):
        raise ValueError("Missing profile IDs cannot be counted as biological types")
    return np.asarray(list(Counter(ids).values()), dtype=np.int64)


def profile_metrics(profile_ids: Iterable[Hashable]) -> dict[str, int | float]:
    """Diversity and species support, taking one profile ID per unique species.

    Singleton fraction has profiles as denominator; singleton_species_fraction
    has species as denominator. Empty data leave undefined ratios as NaN rather
    than treating unobserved assemblages as evidence for zero singleton share.
    """
    counts = _profile_counts(profile_ids)
    n_species, n_profiles = int(counts.sum()), len(counts)
    singleton_profiles = int(np.count_nonzero(counts == 1))
    if n_species:
        proportions = counts / n_species
        effective_shannon = float(np.exp(-np.sum(proportions * np.log(proportions))))
        effective_simpson = float(1 / np.sum(np.square(proportions)))
        largest_share = float(proportions.max())
    else:
        effective_shannon = effective_simpson = 0.0
        largest_share = float("nan")
    return {
        "n_species": n_species,
        "n_profiles": n_profiles,
        "singleton_profiles": singleton_profiles,
        "singleton_fraction": singleton_profiles / n_profiles if n_profiles else float("nan"),
        "singleton_species_fraction": singleton_profiles / n_species if n_species else float("nan"),
        "mean_redundancy": n_species / n_profiles if n_profiles else float("nan"),
        "effective_profiles_shannon": effective_shannon,
        "effective_profiles_simpson": effective_simpson,
        "largest_profile_share": largest_share,
    }


def rarefied_profile_richness(profile_ids: Iterable[Hashable], sample_size: int) -> float:
    """Exact expected number of profiles in a uniform sample of species.

    Sampling is without replacement. For a profile represented by c of N
    species, its observation probability is 1 - choose(N-c,m)/choose(N,m).
    There is no extrapolation beyond the observed number of species.
    """
    counts = _profile_counts(profile_ids)
    n_species = int(counts.sum())
    if (
        isinstance(sample_size, (bool, np.bool_))
        or not isinstance(sample_size, (int, np.integer))
        or sample_size < 0
        or sample_size > n_species
    ):
        raise ValueError("Sample size must be an integer between zero and observed species count")
    if sample_size == 0:
        return 0.0
    if sample_size == n_species:
        return float(len(counts))
    possible_absence = n_species - counts >= sample_size
    log_absence = np.full(len(counts), -np.inf)
    remaining = n_species - counts[possible_absence]
    log_absence[possible_absence] = (
        gammaln(remaining + 1)
        - gammaln(remaining - sample_size + 1)
        - gammaln(n_species + 1)
        + gammaln(n_species - sample_size + 1)
    )
    # expm1 retains precision for rare types whose absence probability is near 1.
    return float(np.sum(-np.expm1(np.minimum(log_absence, 0))))
