"""No private island IDs, datasets or live GBIF requests in unit tests."""
import numpy as np
import pandas as pd
import pytest

from fd.gbif_temporal_opportunity import (
    REGIMES, select_outcome_blind_islands, aggregate_repeated_dataset_year_evidence
)

def _geography():
    return pd.DataFrame([
        {"island_id":f"{r}-{j}", "analysis_regime":r, "spatial_block":f"{r}-block-{j%12}"}
        for r in REGIMES for j in range(25)
    ])


def _records_and_scan(geography):
    rows=[]
    scan=[]
    for r in geography.itertuples(index=False):
        for period,yrs in (("early",(2012,2014)),("recent",(2021,2023))):
            scan.append(dict(island_id=r.island_id,period=period,
                complete=True,truncated=False,api_error="",
                n_raw_examined=6,n_usable_records=6))
            for yr in yrs:
                for k in range(3):
                    rows.append(dict(island_id=r.island_id,period=period,
                         dataset_key="persistent_dataset",year=yr,
                         gbif_id=f"{r.island_id}:{yr}:{k}",
                         genus="Bombus" if k==0 and period=="early" else "Apis"))
    return pd.DataFrame(rows),pd.DataFrame(scan)


def test_outcome_blind_selects_20_each_regime_and_is_stable():
    g=_geography()
    allowed=set(g.island_id)
    selected=select_outcome_blind_islands(g,allowed,allowed)
    assert len(selected)==80 and len(set(selected))==80
    assert selected==select_outcome_blind_islands(g.sample(frac=1,random_state=3),allowed,allowed)
    frame=g[g.island_id.isin(selected)]
    assert frame.analysis_regime.value_counts().to_dict()=={r:20 for r in REGIMES}
    assert frame.groupby("analysis_regime").spatial_block.nunique().min()>=8


def test_repeated_dataset_year_is_not_population_decline():
    g=_geography()
    selected=select_outcome_blind_islands(g,set(g.island_id),set(g.island_id))
    g=g[g.island_id.isin(selected)]
    raw,scans=_records_and_scan(g)
    result=aggregate_repeated_dataset_year_evidence(raw,scans,g,bootstraps=999)
    assert result["n_islands_with_persistent_dataset_year_evidence"]==80
    assert result["n_regimes_reaching_predeclared_20_island_8_block_gate"]==4
    assert result["no_population_decline_or_shift_inferred"] is True
    assert result["no_original_island_dataset_gbif_species_coordinates_exported"] is True
    assert all(x["block_equal_mean_documented_change"] == pytest.approx(-1.) for x in result["by_regime"])


def test_scan_incomplete_cannot_establish_persistent_dataset_comparison():
    g=_geography()
    selected=select_outcome_blind_islands(g,set(g.island_id),set(g.island_id))
    g=g[g.island_id.isin(selected)]
    raw,scans=_records_and_scan(g)
    scans.loc[scans.period.eq("early"),"complete"]=False
    result=aggregate_repeated_dataset_year_evidence(raw,scans,g,bootstraps=999)
    assert result["n_islands_with_persistent_dataset_year_evidence"]==0
    assert result["n_regimes_reaching_predeclared_20_island_8_block_gate"]==0
    assert all(x["status"]=="TEMPORAL_DOCUMENTATION_CHANGE_NOT_IDENTIFIABLE" for x in result["by_regime"])


def test_truncated_window_never_marked_complete():
    g=_geography()
    selected=select_outcome_blind_islands(g,set(g.island_id),set(g.island_id))
    g=g[g.island_id.isin(selected)]
    raw,scans=_records_and_scan(g)
    scans.loc[0,"truncated"]=True
    with pytest.raises(ValueError,match="Truncated"):
        aggregate_repeated_dataset_year_evidence(raw,scans,g,bootstraps=999)
