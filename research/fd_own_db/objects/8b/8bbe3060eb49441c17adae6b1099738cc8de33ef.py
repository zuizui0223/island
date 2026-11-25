"""Why additive mixed-trait Gower diversity cannot detect floral trait coupling.

On genuine same-accepted-species joint flower categories (colour, form),
the usual equally-weighted additive nominal Gower distance is

    D((a,b),(a',b')) = w_a [a!=a'] + w_b [b!=b'].

Its Rao Q is w_a * Q(marginal A) + w_b * Q(marginal B), and therefore
is INVARIANT to changing how colour and form are paired as long as their
marginal category distributions are held fixed.

Nominal Rao Q of the FULL observed joint category (a,b), by contrast,
depends on actual pairing. Its excess over the counterfactual independent
source marginals is a descriptive joint-composition index; neither the
index nor Gower distance identifies pollinator selection or evolution.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Hashable, Sequence

import numpy as np


def nominal_two_axis_diversity(
    source_probability: Sequence[float],
    category_pairs: Sequence[tuple[Hashable, Hashable]],
    *,
    weight_first: float = 0.5,
    weight_second: float = 0.5,
) -> dict[str, float | bool | int]:
    """Exact joint, additive nominal-Gower and independent-marginal Rao Q.

    Category pairs must represent two real trait axes attached to THE SAME
    biological taxon; this function does not fabricate joint occurrence
    records from separate species. Zero-probability combinations may exist
    in the ontology. A missing biological state is NOT a zero-probability
    state and requires separate documented handling.
    """
    p = np.asarray(source_probability, dtype=float)
    if (p.ndim != 1 or len(p) < 2 or not np.isfinite(p).all()
            or (p < 0).any() or not np.isclose(p.sum(), 1, atol=1e-10)):
        raise ValueError("Joint categorical probabilities must form a finite simplex")
    if len(category_pairs) != len(p):
        raise ValueError("The joint category ontology must align with probabilities")
    try:
        if (len(set(category_pairs)) != len(category_pairs)
                or not all(isinstance(pair, tuple) and len(pair) == 2
                           for pair in category_pairs)):
            raise ValueError("All pairs must be unique actual two-axis category labels")
    except TypeError as error:
        raise ValueError("Category labels must be hashable for source axis grouping") from error
    weights = np.asarray([weight_first, weight_second], dtype=float)
    if (not np.isfinite(weights).all() or (weights < 0).any()
            or not np.isclose(weights.sum(), 1, atol=1e-10)):
        raise ValueError("Axis-distance weights must be nonnegative and sum to one")

    first: defaultdict[Hashable, float] = defaultdict(float)
    second: defaultdict[Hashable, float] = defaultdict(float)
    for (a, b), mass in zip(category_pairs, p):
        first[a] += float(mass)
        second[b] += float(mass)

    qa = float(1 - sum(m*m for m in first.values()))
    qb = float(1 - sum(m*m for m in second.values()))
    q_joint = float(1 - p @ p)
    q_additive_gower = float(weight_first * qa + weight_second * qb)
    q_joint_independent = float(
        1 - (1-qa)*(1-qb)
    )
    dependence_excess = float(q_joint - q_joint_independent)
    if not np.isclose(q_joint_independent, qa+qb-qa*qb, atol=1e-12):
        raise AssertionError("Independent marginal state repertoire identity fails")
    return {
        "nominal_Gower_additive_Rao_Q": q_additive_gower,
        "nominal_first_axis_Rao_Q": qa,
        "nominal_second_axis_Rao_Q": qb,
        "true_joint_category_Rao_Q": q_joint,
        "joint_Q_if_source_axes_independent": q_joint_independent,
        "true_joint_Q_minus_independent_marginal_Q": dependence_excess,
        "n_distinct_joint_ontology_states": len(p),
        "n_first_axis_categories": len(first),
        "n_second_axis_categories": len(second),
        "additive_Gower_cannot_detect_cross_axis_association": True,
        "joint_excess_is_descriptive_not_evolutionary_causality": True,
    }
