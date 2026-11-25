"""Marginal versus joint floral set-state similarity kernels.

The original D is exactly the equally weighted mean of Jaccard distances on
COLOUR SETS and FORM SETS, so original K = 1-D is ADDITIVE in these traits.
The additive K-score cannot distinguish colour/form pairings when both
marginals are identical. A Schur-product Jaccard kernel can probe joint
colour/form association without multiplying recorded multistate descriptions
into fabricated profiles.

All kernels are exploratory ecological proxy scores, not pollinator fitness.
All profile identifiers and state sets here come from the PUBLIC 496-row
reference. This module never reads or emits private island/taxon data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping

import numpy as np

from .core import Profile, ProfileSpace, load_profiles
from .pipeline import file_sha256
from .proper_scoring import inspect_kernel

KERNEL_NAMES = ("colour", "form", "additive", "joint_product")


def _set_jaccard(a: frozenset[str], b: frozenset[str]) -> float:
    if not a or not b:
        raise ValueError("Reference trait-state set cannot be empty")
    return len(a & b) / len(a | b)


def phenotype_kernels(profiles: Mapping[int, Profile]) -> dict[str, np.ndarray]:
    """Return K_colour, K_form, original K_additive and K_colour * K_form.

    A multistate value is one unchanged recorded SET (not the Cartesian
    product of species-level phenotype combinations).
    """
    if len(profiles) < 2:
        raise ValueError("Need at least two distinct reference profiles")
    ids = sorted(profiles)
    n = len(ids)
    colour = np.eye(n, dtype=float)
    form = np.eye(n, dtype=float)
    for i, id_i in enumerate(ids):
        c_i, f_i = profiles[id_i]
        for j in range(i + 1, n):
            c_j, f_j = profiles[ids[j]]
            a = _set_jaccard(c_i, c_j)
            b = _set_jaccard(f_i, f_j)
            colour[i, j] = colour[j, i] = a
            form[i, j] = form[j, i] = b
    additive = (colour + form) * .5
    joint = colour * form  # Schur product of positive-semidefinite set kernels
    return {"colour": colour, "form": form,
            "additive": additive, "joint_product": joint}


def kernel_quality(profiles: Mapping[int, Profile],
                   *, verify_psd: bool = True) -> tuple[dict, dict]:
    """Validate additive D identity and strictly distinguishable joint K."""
    kernels = phenotype_kernels(profiles)
    info = {}
    for name, kernel in kernels.items():
        if not np.isfinite(kernel).all() or not np.allclose(
                np.diag(kernel), 1., atol=1e-12):
            raise ValueError(f"Invalid reference state-set kernel: {name}")
        if not np.allclose(kernel, kernel.T, atol=1e-12):
            raise ValueError(f"Asymmetric reference state-set kernel: {name}")
        if verify_psd:
            audit = inspect_kernel(1 - kernel)
            info[name] = {
                "rank_positive_1e8": audit["rank_positive_1e8"],
                "minimum_eigenvalue": audit["minimum_eigenvalue"],
                "n_reference_profiles": len(kernel),
            }
    return kernels, info


def quadratic_kernel_score(prob: np.ndarray, observed_id: int,
                           kernel: np.ndarray) -> float:
    """Proper quadratic kernel score for one categorical heldout observation."""
    p = np.asarray(prob, dtype=float)
    k = np.asarray(kernel, dtype=float)
    if p.ndim != 1 or k.shape != (len(p), len(p)):
        raise ValueError("Predicted probability and kernel dimensions differ")
    if not (0 <= observed_id < len(p)) or np.any(p < 0) or (
            not np.isfinite(p).all()) or not np.isclose(p.sum(), 1):
        raise ValueError("Invalid categorical prediction/target")
    return float(p @ k @ p - 2 * (k @ p)[observed_id] + k[observed_id, observed_id])



def heldout_axis_source_gains(part, model, island: str,
                              kernels: Mapping[str, np.ndarray]) -> dict[str, float]:
    """Post-hoc donor gain split into colour, form, additive and joint scores.

    Scores use the SAME held-out test names, donor model and target island
    for all kernels. Each profile remains its original intact colour/form
    set-state identifier. Joint product tests joint phenotype similarity,
    but a joint gain is NOT by itself a pure interaction effect.
    """
    from .proper_scoring import heldout_probability_scores

    if set(kernels) != set(KERNEL_NAMES):
        raise ValueError("Incomplete trait-kernel score family")
    result = {}
    for name in KERNEL_NAMES:
        scores = heldout_probability_scores(part, model, island, kernels[name])
        result[f"genus_{name}"] = float(
            scores["global_prior_kernel"] - scores["global_genus_kernel"])
        result[f"source_{name}"] = float(
            scores["global_genus_kernel"] - scores["local_genus_kernel"])
    for prefix in ("genus_", "source_"):
        expected = .5 * (result[prefix + "colour"] +
                         result[prefix + "form"])
        if not np.isclose(result[prefix + "additive"], expected, atol=1e-9):
            raise AssertionError("Additive trait-kernel score must equal marginal mean")
    return result

def audit_reference(reference: Path, output: Path) -> dict:
    """Public-only CLI run; log only kernel diagnostics, never taxon records."""
    reference, output = Path(reference), Path(output)
    profiles = load_profiles(reference)
    if len(profiles) != 496 or sorted(profiles) != list(range(496)):
        raise ValueError("Frozen original 496-profile reference changed")
    kernels, info = kernel_quality(profiles)
    # Compare against the INDEPENDENT original distance implementation,
    # not a tautological average of the just-constructed component matrices.
    original_similarity = 1. - ProfileSpace(profiles).matrix
    additive_error = float(np.max(np.abs(
        kernels["additive"] - original_similarity
    )))
    if additive_error > 1e-12:
        raise AssertionError("Trait-score kernel disagrees with original categorical D")
    expected_ranks = {
        "colour": 31, "form": 46, "additive": 76, "joint_product": 496,
    }
    for name, expected in expected_ranks.items():
        if info[name]["rank_positive_1e8"] != expected:
            raise AssertionError(
                f"Frozen public floral kernel rank drift: {name} "
                f"{info[name]['rank_positive_1e8']} != {expected}")
    if info["joint_product"]["minimum_eigenvalue"] <= 1e-8:
        raise AssertionError("Joint colour/form source kernel ceased to be strictly PD")
    # Kadd = 1 - 0.5*(Jaccard colour distance + Jaccard form distance).
    result = {
        "schema": "fd_public_original_floral_trait_kernel_audit_v1",
        "status": "public_reference_mathematical_identity_only",
        "reference_sha256": file_sha256(reference),
        "reference_profiles": len(profiles),
        "kernel_quality": info,
        "additive_formula": "Kadd=(K_colour+K_form)/2=1-original_D",
        "joint_formula": "Kjoint=K_colour * K_form (Hadamard/Schur)",
        "original_distance_max_identity_error": additive_error,
        "frozen_expected_numeric_ranks": expected_ranks,
        "joint_min_eigenvalue_must_be_positive": True,
        "scope": "trait scoring algebra only; no island-level effect, no source prediction, no adaptation claim",
        "interaction_limit": "product kernel can detect joint colour/form association, but a joint K-score gain is not uniquely pure interaction after conditioning on changing marginals",
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "trait_kernel_audit.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"reference_profiles": result["reference_profiles"],
                      "kernels": result["kernel_quality"],
                      "original_distance_max_identity_error": additive_error}))
    return result


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    audit_reference(a.reference, a.output)


if __name__ == "__main__":
    main()
