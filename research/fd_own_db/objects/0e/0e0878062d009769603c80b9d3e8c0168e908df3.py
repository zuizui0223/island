"""Privacy-safe PRE-FROZEN 80-island GBIF Apidae/Bombus two-window pilot.

No survey, effective pollination, or historical population decline is inferred.
Original island IDs, GBIF record IDs, dataset names and georeferences never leave
the private process; only grouped regional counts are returned.
"""
from __future__ import annotations
import hashlib
import numpy as np
import pandas as pd

PERIODS={"early":(2010,2017),"recent":(2018,2025)}
REGIMES=("northern_high_latitude","northern_midlatitude",
         "southern_extratropical","tropical")
REQUIRED_GEO=("island_id","analysis_regime","spatial_block")
REQUIRED_ROW=("island_id","period","dataset_key","year","gbif_id","genus")
REQUIRED_META=("island_id","period","complete","truncated","api_error",
               "n_raw_examined","n_usable_records")


def _stable(item):
    return hashlib.sha256(("fd-bombus-temporal-v1:"+str(item)).encode("utf-8")).hexdigest()


def select_outcome_blind_islands(geo: pd.DataFrame, eligible_islands: set,
                                matched_polygons: set, *, per_regime: int=20):
    """Select exactly 20 per context with geography/block dispersion first.

    No pollinator, functional FD, potential source, or study outcome is read.
    Selection is only from original Direct eligibility and matching original WKB.
    """
    if per_regime!=20 or not set(REQUIRED_GEO).issubset(geo):
        raise ValueError("Frozen pilot sample/geo contract changed")
    if geo[list(REQUIRED_GEO)].isna().any().any() or geo.island_id.duplicated().any():
        raise ValueError("Invalid original island geography")
    x=geo.loc[geo.island_id.isin(eligible_islands & matched_polygons)].copy()
    selected=[]
    for regime in REGIMES:
        sub=x.loc[x.analysis_regime.eq(regime)]
        if len(sub)<per_regime:
            raise ValueError("Outcome-blind regional support below 20 matched islands")
        blocks={}
        for block,g in sub.groupby("spatial_block"):
            blocks[str(block)]=sorted(g.island_id.astype(str),key=_stable)
        keys=sorted(blocks,key=_stable)
        n=0
        while n<per_regime:
            progressed=False
            for block in keys:
                if blocks[block] and n<per_regime:
                    selected.append(blocks[block].pop(0))
                    n+=1
                    progressed=True
            if not progressed: raise ValueError("Original block-stratified region exhausted")
    if len(selected)!=80 or len(set(selected))!=80:
        raise ValueError("Nonunique or incomplete 80-island pilot")
    return selected


def aggregate_repeated_dataset_year_evidence(
        original_rows: pd.DataFrame, scan_meta: pd.DataFrame,
        selected_geography: pd.DataFrame, *,
        bootstraps: int=1999, seed: int=20261113,
):
    """Status and regional response, never individual/island-level output.

    Not even a positive recent-minus-early result can be called decline/gain
    without comparable independent survey effort; datasetKey merely offers a
    minimal common source for a prospective ecological validation.
    """
    if bootstraps<999:raise ValueError("Insufficient frozen block bootstrap")
    for name,tab,fields in (("records",original_rows,REQUIRED_ROW),
                             ("scans",scan_meta,REQUIRED_META),
                             ("geo",selected_geography,REQUIRED_GEO)):
        if not set(fields).issubset(tab):
            raise ValueError(name+" missing frozen original fields")
    geo=selected_geography[list(REQUIRED_GEO)].copy()
    if geo.duplicated("island_id").any() or geo.isna().any().any():
        raise ValueError("Original pilot geography ambiguous")
    if len(geo)!=80 or geo.groupby("analysis_regime").size().to_dict()!={
       regime:20 for regime in REGIMES}:
        raise ValueError("Original frozen 20×4 pilot drifted")
    scans=scan_meta[list(REQUIRED_META)].copy()
    if scans.duplicated(["island_id","period"]).any() or len(scans)!=160:
        raise ValueError("Missing original 80×2 period scan")
    if set(scans.island_id)!=set(geo.island_id) or set(scans.period)!=set(PERIODS):
        raise ValueError("Period coverage does not match predefined islands")
    scans["complete"]=scans.complete.astype(bool)
    scans["truncated"]=scans.truncated.astype(bool)
    scans["api_error"]=scans.api_error.astype(str)
    if (scans.complete & (scans.truncated|scans.api_error.ne(""))).any():
        raise ValueError("Truncated or failed scan cannot be exhaustive")
    for field in ("n_raw_examined","n_usable_records"):
        scans[field]=pd.to_numeric(scans[field],errors="coerce")
        if scans[field].isna().any() or scans[field].lt(0).any() or scans[field].gt(900).any():
            raise ValueError("Invalid per-window GBIF record count / cap")
    rows=original_rows[list(REQUIRED_ROW)].copy()
    if not rows.empty:
        if rows[list(REQUIRED_ROW)].isna().any().any():
            raise ValueError("Missing year/GBIF ID/provider on qualifying record")
        if rows.duplicated(["island_id","period","gbif_id"]).any():
            raise ValueError("Duplicate GBIF ID within island × period")
        if not set(rows.island_id).issubset(set(geo.island_id)) or not set(rows.period).issubset(PERIODS):
            raise ValueError("Record outside original 80-island period frame")
        rows["year"]=pd.to_numeric(rows.year,errors="coerce")
        if rows.year.isna().any():raise ValueError("Non-numeric year")
        for window,(a,b) in PERIODS.items():
            if not rows.loc[rows.period.eq(window),"year"].between(a,b).all():
                raise ValueError("GBIF record outside window")
        rows["bombus"]=rows.genus.astype(str).str.strip().eq("Bombus")
        # The GBIF API is a source of records, not per-visit surveys.
        yearly=rows.groupby(["island_id","dataset_key","period","year"],as_index=False).agg(
           background_records=("gbif_id","nunique"),target_present=("bombus","max"))
        # Require >=3 independent Apidae GBIF records per year×dataset cell.
        yearly=yearly.loc[yearly.background_records.ge(3)].copy()
        grouped=yearly.groupby(["island_id","dataset_key","period"],as_index=False).agg(
           n_observed_years=("year","nunique"),
           bombus_year_fraction=("target_present","mean"),
           n_background_records=("background_records","sum"))
    else:
        grouped=pd.DataFrame(columns=["island_id","dataset_key","period",
                                      "n_observed_years","bombus_year_fraction","n_background_records"])
    ready=scans.set_index(["island_id","period"])["complete"]
    comparisons=[]
    for (island,provider),g in grouped.groupby(["island_id","dataset_key"]):
        if not (ready.get((island,"early"),False) and ready.get((island,"recent"),False)):
            continue
        if set(g.period)!={"early","recent"}:
            continue
        a=g.set_index("period")
        if not (a.loc["early","n_observed_years"]>=2 and a.loc["recent","n_observed_years"]>=2):
            continue
        comparisons.append(dict(island_id=island,dataset_key=provider,
          recent_minus_early=float(a.loc["recent","bombus_year_fraction"]-
                                   a.loc["early","bombus_year_fraction"])))
    cmp=pd.DataFrame(comparisons)
    islands=(cmp.groupby("island_id").recent_minus_early.mean() if not cmp.empty
             else pd.Series(dtype=float))
    private=geo.merge(islands.rename("recorded_delta"),left_on="island_id",right_index=True,how="left")
    complete=scans.pivot(index="island_id",columns="period",values="complete")
    complete_both=complete.index[complete["early"] & complete["recent"]]
    raw=scans.pivot(index="island_id",columns="period",values="n_raw_examined")
    out=[]
    for regime,group in private.groupby("analysis_regime"):
        ids=set(group.island_id)
        available=group.dropna(subset=["recorded_delta"])
        rec=dict(regime=regime,n_outcome_blind_selected=len(group),
            n_original_spatial_blocks=group.spatial_block.nunique(),
            n_both_windows_exhaustive=len(ids & set(complete_both)),
            n_persistent_dataset_islands=len(available),
            n_persistent_dataset_spatial_blocks=available.spatial_block.nunique(),
            n_qualifying_island_dataset_pairs=int(
                (cmp.island_id.isin(ids)).sum()) if not cmp.empty else 0,
            n_early_raw_records_scanned=int(scans.loc[
                scans.island_id.isin(ids)&scans.period.eq("early"),"n_raw_examined"].sum()),
            n_recent_raw_records_scanned=int(scans.loc[
                scans.island_id.isin(ids)&scans.period.eq("recent"),"n_raw_examined"].sum()),
            change_from_GBIF_documentation_not_population_decline=True)
        # No inferential statistic if fewer than 20 matched islands or 8 blocks
        # in a regional sample of 20. Return only support by default.
        if len(available)>=20 and available.spatial_block.nunique()>=8:
            means=available.groupby("spatial_block").recorded_delta.mean().to_numpy(float)
            rng=np.random.default_rng(seed+REGIMES.index(regime))
            draws=means[rng.integers(len(means),size=(bootstraps,len(means)))].mean(axis=1)
            rec.update(status="DOCUMENTATION_CHANGE_ESTIMABLE_NOT_ECOLOGICAL_DECLINE",
              block_equal_mean_documented_change=float(means.mean()),
              conditional_95_low=float(np.quantile(draws,.025)),
              conditional_95_high=float(np.quantile(draws,.975)))
        else:
            rec["status"]="TEMPORAL_DOCUMENTATION_CHANGE_NOT_IDENTIFIABLE"
        out.append(rec)
    return dict(
       status="ORIGINAL_GBIF_TEMPORAL_PILOT_RECORDING_ONLY",
       data_type="opportunistic_GBIF_occurrence_not_standardized_surveys",
       fixed_islands=80, windows={k:list(v) for k,v in PERIODS.items()},
       max_records_examined_per_island_period=900,
       n_islands_with_persistent_dataset_year_evidence=len(islands),
       n_regimes_reaching_predeclared_20_island_8_block_gate=sum(
           x["status"].startswith("DOCUMENTATION_CHANGE_ESTIMABLE") for x in out),
       by_regime=out,
       target="Bombus_genus_occurrences_within_Apidae_background",
       year_dataset_observability_min_apidae_records=3,
       minimum_years_per_window_per_persistent_dataset=2,
       no_population_decline_or_shift_inferred=True,
       no_longitudinal_flora_or_fitness_measurement=True,
       missing_GBIF_target_or_truncated_window_never_biological_absence=True,
       no_original_island_dataset_gbif_species_coordinates_exported=True,
    )
