"""Geographically disjoint geology-match controls: synthetic outcome-independent pairs."""
import pandas as pd
import numpy as np
from fd.geological_disjoint_matched import METRICS, matched_oceanic_continental_source_audit


def _example():
    rows=[]
    for i in range(11):
        for origin in ("continental","oceanic"):
            q=dict(
                island_id=f"{origin}_{i}", origin=origin,
                regime="tropical",scope="all_recorded",
                spatial_block=f"{origin}_block_{i}",
                archipelago=f"{origin}_arch_{i}",
                area_km2=95+i*4,
                distance_to_continent_km=60+i*6,
                n_original_direct_species=40+i,
                n_source_genera=50+i,
                climate_pc1=0.1*i,
                climate_pc2=0.2*i,
                climate_pc3=0.3*i,
                climate_pc4=0.4*i,
            )
            for metric in METRICS:
                q[metric]= 0.01*i + (0.004 if origin=="oceanic" else 0.0)
            rows.append(q)
    return pd.DataFrame(rows)


def test_compatible_11pairs_no_outcome_driven_eligibility():
    data=_example()
    out=matched_oceanic_continental_source_audit(data)
    a=[x for x in out if x["scope"]=="all_recorded" and x["region"]=="tropical"][0]
    assert a["n_spatial_and_verified_archipelago_disjoint_pairs"]>=6
    assert a["n_unknown_archipelago_endpoints"]==0
    key="pair_oceanic_minus_continental_source_relative_marginal_rao"
    assert a[key+"_n"]>=6
    assert a[key+"_exact_two_sided_p"]>=0
    assert a[key+"_exact_two_sided_p"]<=1
    assert "island_id" not in a and "archipelago" not in a
    # Changing floral outcome values must never affect chosen geological pair count.
    data2=data.copy()
    for met in METRICS:
        data2[met]*=-1
    b=[x for x in matched_oceanic_continental_source_audit(data2)
       if x["scope"]=="all_recorded" and x["region"]=="tropical"][0]
    assert a["n_spatial_and_verified_archipelago_disjoint_pairs"]==b[
        "n_spatial_and_verified_archipelago_disjoint_pairs"]


def test_missing_covariate_does_not_become_zero_or_imputed():
    data=_example()
    data.loc[data.index[:3],"climate_pc2"]=np.nan
    data.loc[data.index[5],"distance_to_continent_km"]=-1
    out=matched_oceanic_continental_source_audit(data)
    t=[x for x in out if x["scope"]=="all_recorded" and x["region"]=="tropical"][0]
    assert t["n_source_eligible_islands_before_geography_complete_case"]==22
    assert t["n_source_eligible_islands_complete_geo_covariates"]==18
