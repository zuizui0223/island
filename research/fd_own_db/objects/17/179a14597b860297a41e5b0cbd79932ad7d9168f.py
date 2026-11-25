"""Distance-based floral phenotype diversity on categorical state sets.

The data are recorded flower-colour and floral-form state sets, not direct
measurements of pollinator performance. Do not expand multistate sets into
unobserved flower combinations. Calculations use the original profile
dissimilarities, never a low-dimensional PCoA projection.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Mapping

import numpy as np

Profile = tuple[frozenset[str], frozenset[str]]


def csv_rows(path: str | Path):
    """Read regular or gzip-compressed UTF-8 CSV by actual compression suffix."""
    filename = Path(path)
    opener = gzip.open if filename.suffix == ".gz" else open
    with opener(filename, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def _states(raw: str) -> frozenset[str]:
    values = raw.split("|") if isinstance(raw, str) else []
    if not values or any(not s or s != s.strip() for s in values):
        raise ValueError("Missing, empty or untrimmed floral state")
    if len(values) != len(set(values)):
        raise ValueError("Duplicated state in a recorded floral state set")
    return frozenset(values)


def load_profiles(path: str | Path) -> dict[int, Profile]:
    """Load the 496-profile common reference from island/all/profile_coordinates.csv.

    The PCoA columns in that CSV are intentionally ignored.
    """
    result: dict[int, Profile] = {}
    for row in csv_rows(path):
        key = int(row["profile_id"])
        if key in result:
            raise ValueError(f"Duplicate profile ID: {key}")
        result[key] = (_states(row["flower_primary_color"]), _states(row["floral_form"]))
    if not result:
        raise ValueError("No reference profiles")
    if len(set(result.values())) != len(result):
        raise ValueError("Different IDs encode duplicate state-set profiles")
    return result


def distance(a: Profile, b: Profile) -> float:
    """Equal-weight mean Jaccard distance between colour and form state sets."""
    return sum(1.0 - len(x & y) / len(x | y) for x, y in zip(a, b)) / 2.0


class ProfileSpace:
    """Common original-distance matrix reused for every scope and island."""

    def __init__(self, profiles: Mapping[int, Profile]) -> None:
        self.ids = sorted(profiles)
        self.lookup = {key: idx for idx, key in enumerate(self.ids)}
        n = len(self.ids)
        self.matrix = np.zeros((n, n), dtype=np.float64)
        for i in range(n):
            for j in range(i + 1, n):
                d = distance(profiles[self.ids[i]], profiles[self.ids[j]])
                self.matrix[i, j] = self.matrix[j, i] = d
        if not np.isfinite(self.matrix).all() or np.any(self.matrix < 0):
            raise ValueError("Invalid profile dissimilarity")

    def vector(self, carriers: Mapping[int, int]) -> tuple[np.ndarray, np.ndarray, int]:
        if not carriers:
            raise ValueError("Empty island is not a zero-diversity observation")
        missing = set(carriers) - set(self.lookup)
        if missing:
            raise ValueError(f"Unknown profile IDs: {sorted(missing)[:5]}")
        keys = sorted(carriers)
        counts = np.asarray([carriers[k] for k in keys], dtype=np.float64)
        if np.any(counts <= 0) or not np.isfinite(counts).all() or np.any(counts != np.floor(counts)):
            raise ValueError("Carrier counts must be positive integers")
        indices = np.asarray([self.lookup[k] for k in keys], dtype=np.intp)
        return indices, counts / counts.sum(), int(counts.sum())

    def island(self, carriers: Mapping[int, int]) -> dict[str, float | int]:
        ids, weights, n = self.vector(carriers)
        sub = self.matrix[np.ix_(ids, ids)]
        q = float(weights @ sub @ weights)
        equal = np.full(len(ids), 1.0 / len(ids))
        return {
            "n_names": n,
            "n_profiles": len(ids),
            "rao_q": q,
            # Finite-name correction is a descriptive counterpart, not a
            # correction for nonrandom trait missingness.
            "rao_q_without_replacement": q * n / (n - 1) if n > 1 else float("nan"),
            "rao_q_profile_equal": float(equal @ sub @ equal),
        }

    def pair(self, a: Mapping[int, int], b: Mapping[int, int]) -> dict[str, float | int]:
        """Whole-island set/weight comparisons; shared taxa are NOT removed.

        Weighted nearest-neighbour distance = mean of the two directional
        nearest-profile distances, with carriers weighting the source profiles.
        This measures functional separation of recorded profile supports; it is
        not an exact taxon-excluded or genus-conditioned beta-diversity test.
        """
        ia, pa, na = self.vector(a)
        ib, pb, nb = self.vector(b)
        cross = self.matrix[np.ix_(ia, ib)]
        sa, sb = set(a), set(b)
        return {
            "n_names_a": na,
            "n_names_b": nb,
            "n_profiles_a": len(ia),
            "n_profiles_b": len(ib),
            "profile_jaccard": len(sa & sb) / len(sa | sb),
            "nn_ab": float(pa @ cross.min(axis=1)),
            "nn_ba": float(pb @ cross.min(axis=0)),
            "functional_nn_dissimilarity": float(
                (pa @ cross.min(axis=1) + pb @ cross.min(axis=0)) / 2.0
            ),
            "mean_cross_profile_distance": float(pa @ cross @ pb),
        }


def load_carriers(path: str | Path) -> dict[str, dict[int, int]]:
    """Public aggregated island x profile x number-of-name-carriers input."""
    result: dict[str, dict[int, int]] = {}
    for row in csv_rows(path):
        island = row["island_id"]
        if not island or island != island.strip():
            raise ValueError("Blank or invalid island ID")
        profile = int(row["profile_id"])
        count = int(row["n_species"])
        if count < 1:
            raise ValueError("Each observed profile must have >=1 carrier")
        group = result.setdefault(island, {})
        if profile in group:
            raise ValueError(f"Duplicate island/profile: {island}/{profile}")
        group[profile] = count
    return result


def load_island_metadata(path: str | Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in csv_rows(path):
        island = row["island_id"]
        if not island or island in result:
            raise ValueError("Blank or duplicate island in metadata")
        result[island] = row
    return result
