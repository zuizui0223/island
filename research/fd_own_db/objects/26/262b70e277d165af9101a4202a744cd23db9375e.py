"""Marginal-preserving test of colour x form association in floral source predictions.

Rake each local-source categorical probability vector to EXACTLY match the
global-genus model's predicted COLOUR-SET and FORM-SET marginal probabilities.
Comparing the raked local model with global genus via the strictly proper
joint-profile Brier and full-rank product-Jaccard kernel scores tests only
conditional colour/form association beyond the *same fixed marginals*.

All colour/form pairings remain the 496 ORIGINAL recorded profiles; no new
Cartesian-product categories are fabricated. Selection of source/placebo
controls remains outcome-blind in n2_source_assignment_placebo.

This is a private-input-only observational prediction test, never evidence
of selected pollinator syndromes or actual ancestral mainland sources.
"""
from __future__ import annotations

from typing import Mapping

import numpy as np
import pandas as pd

from .core import Profile
from .trait_kernels import phenotype_kernels


def profile_axis_codes(profiles: Mapping[int, Profile]) -> tuple[np.ndarray, np.ndarray]:
    """Stable IDs for complete multistate colour and form SET categories."""
    ids = sorted(profiles)
    colours = {v[0] for v in profiles.values()}
    forms = {v[1] for v in profiles.values()}
    c_id = {c: i for i, c in enumerate(sorted(colours, key=lambda x: tuple(sorted(x))))}
    f_id = {f: i for i, f in enumerate(sorted(forms, key=lambda x: tuple(sorted(x))))}
    return (np.asarray([c_id[profiles[i][0]] for i in ids], dtype=int),
            np.asarray([f_id[profiles[i][1]] for i in ids], dtype=int))


def margins(p: np.ndarray, colour: np.ndarray,
            form: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(p, dtype=float)
    c, f = np.asarray(colour), np.asarray(form)
    if x.ndim != 1 or c.shape != x.shape or f.shape != x.shape:
        raise ValueError("Marginal dimensions are inconsistent")
    return (np.bincount(c, weights=x, minlength=int(c.max()) + 1),
            np.bincount(f, weights=x, minlength=int(f.max()) + 1))


def match_marginals(
    global_genus: np.ndarray, local_genus: np.ndarray,
    colour: np.ndarray, form: np.ndarray, *,
    tol: float = 2e-7, max_iter: int = 3000,
) -> tuple[np.ndarray, dict]:
    """KL-minimal IPF projection of local predictions onto global marginals.

    Both distributions must be strictly positive on the SAME actual
    496-profile support. Multiplicative row/column scaling preserves
    local colour/form association odds ratios on every supported 2x2
    subtable while forcing both predicted marginal distributions to be
    identical. Return convergence metadata to prevent silent bias.
    """
    g, w = np.asarray(global_genus, dtype=float), np.asarray(local_genus, dtype=float).copy()
    c, f = np.asarray(colour, dtype=int), np.asarray(form, dtype=int)
    if (g.shape != w.shape or g.ndim != 1 or c.shape != g.shape or f.shape != g.shape
            or np.any(c < 0) or np.any(f < 0)
            or np.any(g <= 0) or np.any(w <= 0)
            or not np.isfinite(g).all() or not np.isfinite(w).all()
            or not np.isclose(g.sum(), 1., atol=1e-10)
            or not np.isclose(w.sum(), 1., atol=1e-10)):
        raise ValueError("IPF requires strictly positive, normalized matched-support models")
    if max_iter < 1 or tol <= 0:
        raise ValueError("Invalid IPF tolerance/iteration budget")
    target_c, target_f = margins(g, c, f)
    converged = False
    max_diff = float("inf")
    for j in range(max_iter):
        current_c = np.bincount(c, weights=w, minlength=len(target_c))
        if np.any(current_c <= 0):
            raise ValueError("IPF lost colour support")
        w *= target_c[c] / current_c[c]
        current_f = np.bincount(f, weights=w, minlength=len(target_f))
        if np.any(current_f <= 0):
            raise ValueError("IPF lost form support")
        w *= target_f[f] / current_f[f]
        if (j + 1) % 5 == 0 or j == max_iter - 1:
            now_c, now_f = margins(w, c, f)
            max_diff = max(float(np.max(np.abs(now_c - target_c))),
                           float(np.max(np.abs(now_f - target_f))))
            if max_diff <= tol:
                converged = True
                break
    if not converged:
        raise RuntimeError(f"IPF failed to match marginals, max_err={max_diff:g}")
    if not np.isclose(w.sum(), 1., atol=1e-9):
        raise AssertionError("Raked probability mass does not sum to one")
    return w, {"iterations": j + 1, "max_abs_marginal_error": max_diff}


def brier_and_joint_kernel_loss(
    probabilities: np.ndarray, observations: np.ndarray,
    joint_kernel: np.ndarray,
) -> tuple[float, float]:
    """Strictly proper Brier and joint positive-definite kernel scores."""
    p = np.asarray(probabilities, dtype=float)
    obs = np.asarray(observations, dtype=int)
    k = np.asarray(joint_kernel, dtype=float)
    if (p.ndim != 1 or k.shape != (len(p), len(p)) or obs.ndim != 1
            or len(obs) == 0 or np.any(obs < 0) or np.any(obs >= len(p))):
        raise ValueError("Invalid full-category prediction/observation")
    kp = k @ p
    brier = float(1 + p @ p - 2 * p[obs].mean())
    joint = float(1 + p @ kp - 2 * kp[obs].mean())
    return brier, joint


def marginal_fixed_source_gain(
    observed: pd.DataFrame, model, focal_source: str,
    reference_profiles: Mapping[int, Profile], *,
    colour_codes: np.ndarray | None = None,
    form_codes: np.ndarray | None = None,
    joint_kernel: np.ndarray | None = None,
) -> dict[str, float]:
    """One island vs source assignment, holding genus predictions' marginals.

    The global genus model is the same for any source assignment, and is
    the fixed marginal target. This returns only mean *score* diagnostics;
    never serialize or print individual island, genus or taxon identities.
    """
    if observed.empty or observed.accepted_species.duplicated().any():
        raise ValueError("Invalid scored island name units")
    if colour_codes is None or form_codes is None:
        colour_codes, form_codes = profile_axis_codes(reference_profiles)
    if joint_kernel is None:
        joint_kernel = phenotype_kernels(reference_profiles)["joint_product"]
    lookup = {key: index for index, key in enumerate(sorted(reference_profiles))}
    gains = np.zeros(2, dtype=float)
    total = 0
    iterations = 0
    largest_error = 0.
    n_genus = 0
    for genus, group in observed.groupby("genus", sort=False):
        pred = np.asarray(model.probabilities(focal_source, genus)[[2, 3]], dtype=float)
        global_p, local_p = pred
        matched, info = match_marginals(global_p, local_p, colour_codes, form_codes)
        obs = np.asarray([lookup[int(x)] for x in group.profile_id], dtype=int)
        base = brier_and_joint_kernel_loss(global_p, obs, joint_kernel)
        candidate = brier_and_joint_kernel_loss(matched, obs, joint_kernel)
        gains += len(group) * (np.asarray(base) - np.asarray(candidate))
        total += len(group)
        iterations += info["iterations"]
        largest_error = max(largest_error, info["max_abs_marginal_error"])
        n_genus += 1
    return {
        "gain_marginal_fixed_brier": float(gains[0]/total),
        "gain_marginal_fixed_joint": float(gains[1]/total),
        "mean_ipf_iterations": iterations/n_genus,
        "maximum_ipf_marginal_error": largest_error,
        "n_evaluated_name_carriers": total,
    }


def marginal_fixed_placebo_gain(
    focal: pd.DataFrame, model, actual: str, controls: list[str],
    profiles: Mapping[int, Profile], *,
    colour_codes: np.ndarray | None = None,
    form_codes: np.ndarray | None = None,
    joint_kernel: np.ndarray | None = None,
) -> dict:
    """Actual vs exactly same two controls, with global genus marginals fixed.

    Precompute reference kernels ONCE per run for large island cohorts:
    constructing the 496x496 kernel separately per target is costly.
    """
    if len(controls) != 2 or actual in controls or controls[0] == controls[1]:
        raise ValueError("Require two distinct alternative source assignments")
    if colour_codes is None or form_codes is None:
        colour, form = profile_axis_codes(profiles)
    else:
        colour, form = colour_codes, form_codes
    kernel = (phenotype_kernels(profiles)["joint_product"] if
              joint_kernel is None else joint_kernel)
    actual_gain = marginal_fixed_source_gain(
        focal, model, actual, profiles,
        colour_codes=colour, form_codes=form, joint_kernel=kernel)
    alternative = [marginal_fixed_source_gain(
        focal, model, candidate, profiles,
        colour_codes=colour, form_codes=form, joint_kernel=kernel)
        for candidate in controls]
    return {
        "actual_minus_placebo_marginal_fixed_joint": float(
            actual_gain["gain_marginal_fixed_joint"] -
            np.mean([a["gain_marginal_fixed_joint"] for a in alternative])),
        "actual_minus_placebo_marginal_fixed_brier": float(
            actual_gain["gain_marginal_fixed_brier"] -
            np.mean([a["gain_marginal_fixed_brier"] for a in alternative])),
        "maximum_ipf_error": float(max(
            [actual_gain["maximum_ipf_marginal_error"]]
            + [a["maximum_ipf_marginal_error"] for a in alternative])),
        "mean_ipf_iterations": float(np.mean(
            [actual_gain["mean_ipf_iterations"]]
            + [a["mean_ipf_iterations"] for a in alternative])),
    }
