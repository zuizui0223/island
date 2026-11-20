"""Synthetic-only original GIFT HIGH donor source predictor matched target audit."""
import numpy as np
import pytest

from fd_original_stage2_name_disjoint_flower_forecast import (
    HIGH_GRADE_PILOT_ENDPOINTS, HIGH_GRADE_GAINS,
    high_grade_source_gain,high_grade_pair_positive,
)
from fd.two_stage_floral_forecast import proper_losses


def test_high_grade_donor_and_old_grade_source_score_exact_same_target_codes():
    # Exactly the same held-out accepted species real Direct category indices;
    # only the SOURCE OTHER species evidence-quality restriction changes.
    codes=[0,0,0,0,1,0,0,1,0,0]
    p_full=np.array([.6,.4])
    p_high=np.array([.8,.2])
    diag=high_grade_source_gain(codes,p_high,p_full)
    assert set(diag)==set(HIGH_GRADE_GAINS)
    for j,loss in enumerate(("Brier","logloss")):
        actual=proper_losses(p_full,codes)[j]-proper_losses(p_high,codes)[j]
        assert diag[f"high_grade_external_congener_gain_vs_high_medium_{loss}"]==pytest.approx(actual)
    assert diag[HIGH_GRADE_GAINS[0]]>0
    assert diag[HIGH_GRADE_GAINS[1]]>0


def test_high_grade_comparison_does_not_change_biological_target_cohort():
    assert HIGH_GRADE_PILOT_ENDPOINTS==("floral_symmetry","flower_primary_color")
    codes=[0,1,0,1,0]
    source=[.5,.5]
    improved=[.7,.3]
    a=high_grade_source_gain(codes,improved,source)
    b=high_grade_source_gain(list(reversed(codes)),improved,source)
    assert set(a)==set(b)
    # If all plant states are used, permutation of target plant order cannot
    # change a proper-score comparison.
    assert a==pytest.approx(b,abs=1e-12)


def test_high_grade_dual_loss_and_two_original_geographic_group_gates():
    def summar(lo_b,lo_l,qual=True):
        return {
            "supported_for_inference":qual,
            HIGH_GRADE_GAINS[0]:{"bootstrap_95":[lo_b,.03]},
            HIGH_GRADE_GAINS[1]:{"bootstrap_95":[lo_l,.04]},
        }
    assert high_grade_pair_positive(
        summar(.002,.003),summar(.001,.001))
    assert not high_grade_pair_positive(
        summar(-.002,.003),summar(.001,.001))
    assert not high_grade_pair_positive(
        summar(.002,.003),summar(.001,-.004))
    assert not high_grade_pair_positive(
        summar(.002,.003),summar(.001,.001,qual=False))


def test_original_high_model_source_probabilities_must_be_finite_and_same_ontology():
    with pytest.raises(ValueError):
        high_grade_source_gain([0,1,0,0,0],[.1,.1,.8],[.6,.4])
    with pytest.raises(ValueError):
        high_grade_source_gain([0,1,0,0,0],[.8,.2],[-.1,1.1])
