"""Synthetic privacy + robust block-bootstrap tests for source coupling."""
import numpy as np
import pandas as pd
import pytest

from fd.coupling_heterogeneity import (
    heterogeneity_summary, trim_winsor, OUTPUT_FIELDS,
)


def test_ten_percent_trim_winsor_discards_outlying_island_changes():
    values=np.concatenate([np.zeros(18),[100.,-100.]])
    a,b,median=trim_winsor(values)
    assert a==pytest.approx(0.)
    assert b==pytest.approx(0.)
    assert median==0.
    with pytest.raises(ValueError):
        trim_winsor([.1,.2,.3])


def test_block_bootstrap_heterogeneity_and_no_private_island_metadata():
    rows=[]
    for block in range(6):
        for i in range(6):
            rows.append({
                "block":f"private_block{block}",
                "island_id":f"private_island{block}_{i}",
                "accepted_species":f"secret_taxon{block}_{i}",
                "gain_joint_marginal_fixed":
                    (.0005 if i%2==0 else -.0001) + block*.00001,
            })
    frame=pd.DataFrame(rows)
    a=heterogeneity_summary(frame,scope="all",draws=299,seed=17)
    b=heterogeneity_summary(frame,scope="all",draws=299,seed=17)
    assert a==b
    assert a["n_islands"]==36 and a["n_geographic_blocks"]==6
    assert a["fraction_positive_islands"]==.5
    assert a["trimmed10_mean"]>0
    assert a["trimmed10_block_q025"]>0
    assert 0<=a["share_top10pct_of_positive_gains"]<=1
    assert set(a)==set(OUTPUT_FIELDS)
    assert not {"block","island_id","accepted_species"} & a.keys()


def test_nonindependent_single_block_and_bad_score_fail_closed():
    part=pd.DataFrame({"block":["same"]*12,
         "gain_joint_marginal_fixed":[.1]*12})
    with pytest.raises(ValueError):
        heterogeneity_summary(part,scope="all",draws=299)
    part["block"]=[f"B{i}" for i in range(12)]
    part.loc[0,"gain_joint_marginal_fixed"]=np.nan
    with pytest.raises(ValueError):
        heterogeneity_summary(part,scope="all",draws=299)
