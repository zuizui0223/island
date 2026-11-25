"""Tests for coarse fold diagnostic, not a new inferential p-value."""
import numpy as np
import pandas as pd
import pytest

from fd.fold_coupling_sensitivity import fold_sign_sensitivity, ENDPOINTS


def test_five_fold_positive_differences_have_exact_discrete_tail():
    rows=[]
    for scope in ("all","direct"):
        for fold in range(5):
            rows.append(dict(
                scope=scope,fold=fold,n_islands=fold+1,
                mean_joint_before=.001,
                mean_joint_marginal_fixed=.0002,
                mean_brier_marginal_fixed=.0003,
            ))
    out=fold_sign_sensitivity(pd.DataFrame(rows))
    assert len(out)==6
    assert all(row["n_positive_fold_means"]==5 for row in out)
    assert all(row["sign_configurations"]==32 for row in out)
    assert all(row["two_sided_sign_tail_if_exchangeable"]==pytest.approx(2/32)
               for row in out)
    assert all(row["leave_one_fold_out_min"]>0 for row in out)


def test_one_negative_fold_cannot_claim_five_independent_confirmations():
    data=pd.DataFrame([dict(
        scope="all",fold=i,n_islands=10,
        mean_joint_before=.002,
        mean_joint_marginal_fixed=.01 if i<4 else -.002,
        mean_brier_marginal_fixed=.003,
    ) for i in range(5)])
    results=fold_sign_sensitivity(data)
    joint=next(x for x in results if x["endpoint"]=="mean_joint_marginal_fixed")
    assert joint["n_positive_fold_means"]==4
    assert joint["two_sided_sign_tail_if_exchangeable"]>=2/32
    assert joint["leave_one_fold_out_min"]<=joint["leave_one_fold_out_max"]


def test_missing_folds_fail_closed():
    rows=[dict(scope="all",fold=i,n_islands=20,
        **{m:.01 for m in ENDPOINTS}) for i in range(4)]
    with pytest.raises(ValueError,match="Exactly five"):
        fold_sign_sensitivity(pd.DataFrame(rows))
    rows.append(dict(rows[0]))
    with pytest.raises(ValueError):
        fold_sign_sensitivity(pd.DataFrame(rows))
