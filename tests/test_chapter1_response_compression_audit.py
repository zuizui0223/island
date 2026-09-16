import numpy as np

from island_v2.chapter1_response_compression_audit import _contrast_matrices


def test_contrast_matrices_partition_two_three_member_families():
    outcomes = ["a1", "a2", "a3", "r1", "r2", "r3"]
    config = {
        "response_families": {
            "accessibility_generalization": ["a1", "a2", "a3"],
            "reproductive_assurance": ["r1", "r2", "r3"],
        }
    }
    family, within = _contrast_matrices(outcomes, config)
    assert family.shape == (2, 6)
    assert within.shape == (4, 6)
    np.testing.assert_allclose(family.sum(axis=1), [1.0, 1.0])
    # Within-family contrasts are orthogonal to both family means.
    np.testing.assert_allclose(family @ within.T, np.zeros((2, 4)))
