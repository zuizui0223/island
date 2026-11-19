"""GIFT original 4505 island GENUS COMMUNITY prediction using geography only.

Source-trait-free first stage towards predicting floral functional repertoire.
Predict the normalized genus composition of an island's accepted botanical
flora from OTHER original GIFT island floras whose geographic 10-degree block
AND first candidate GIFT mainland-source entity differ from the focal. The
focal botanical genus roster/plant species are read ONLY for proper loss
scoring, never for selecting genera or geographic neighbors.

Freeze source geography: original GSHHG latitude/longitude, climate_pc1..4,
log1p(area_km2) and log1p(distance_to_continent_km), and GIFT source geo_k5
candidate entity labels (as availability proxies, NOT independently observed
ancestral floras). Candidate models:
  equal region source island genus profile (baseline);
  geographically closest 40 source island floras (spatial benchmark);
  40 most similar source climate PCs (climate benchmark);
  40 source islands similar in climate+area+isolation (environment);
  40 source islands selected by environment + original SOURCE candidate
  ranked-entity Jaccard (environment and source opportunity).
All predictions share the EXACT SAME source islands and a fixed shrinkage
toward the same region-profile baseline. Genus names absent from every
eligible source flora are pooled as a SINGLE unknown genus class with the
same p_unknown=.02 across models; this avoids a target taxon-name oracle.

Scoring uses the target accepted whole flora's source genus membership only
AFTER fitting every predictive distribution. Checklist nonlisting is NOT
ecological true species/genus absence; this is an observed floristic taxon
composition forecast, NOT a causal colonization establishment probability.

All output is anonymous aggregates, no accepted names, genus names, source
entity identifiers, target island coords or source species rosters.
"""
from __future__ import annotations
from collections import Counter,defaultdict
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from fd_interregional_genus_name_exclusive_transfer import genus,sha

K=40
DELTA_UNKNOWN=.02
SHRINK=.20
SOURCE_WEIGHT=.75
NBOOT=399
MIN_SOURCE_ISLANDS=30
MIN_FOCAL_SPECIES=10
MIN_FOCAL_GENERA=5
MIN_EVAL_ISLANDS=30
MIN_CLUSTERS=10
REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude")
SCOPES=("all_recorded","original_native")
FEATURES=("climate_pc1","climate_pc2","climate_pc3","climate_pc4",
          "area_km2","distance_to_continent_km")
MODELS=("geo_nearest","climate_nearest","environment_nearest",
        "environment_source_nearest")
OUTPUT_COLS=tuple(f"{m}_minus_region_equal_Brier" for m in MODELS)+tuple(
    f"{m}_minus_region_equal_logloss" for m in MODELS)
GROUPS=("spatial_block","first_source")


def _floats(g:pd.DataFrame)->np.ndarray:
    values=g.loc[:,FEATURES].apply(pd.to_numeric,errors="coerce").to_numpy(dtype=float,copy=True)
    if (values[:,4:6]<0).any():
        raise ValueError("Original GSHHG island geographic area/distance negative")
    values[:,4:6]=np.log1p(values[:,4:6])
    valid=np.isfinite(values).all(axis=1)
    if not valid.any():
        raise ValueError("Original GIFT GSHHG all climate PCs missing")
    # Geography is PUBLIC and knows no botanical response; fixed source
    # normalization may include target GEOGRAPHICAL covariates, not taxa.
    med=np.nanmedian(values[valid],axis=0)
    q75=np.nanquantile(values[valid],.75,axis=0)
    q25=np.nanquantile(values[valid],.25,axis=0)
    scale=np.maximum(q75-q25,1e-6)
    standardized=np.clip((values-med)/scale,-10.,10.)
    return standardized


def _xyz(lat,lon)->np.ndarray:
    a,b=np.radians(np.asarray(lat,float)),np.radians(np.asarray(lon,float))
    xyz=np.column_stack((np.cos(a)*np.cos(b),
        np.cos(a)*np.sin(b),np.sin(a)))
    # Sphere embedding avoids spurious distance across the 180-degree meridian.
    return xyz


def genus_source_matrix(floras:pd.DataFrame,ids:list[str],scope:str):
    if not {"island_id","accepted_species","origin_status"}<=set(floras):
        raise ValueError("Original GIFT full accepted species/name/status fields missing")
    if floras.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate original botanical accepted island×species")
    if scope not in SCOPES:raise ValueError("Original native/all-recorded scope")
    sel=(floras if scope=="all_recorded" else floras.loc[
        floras.origin_status.eq("native")])
    n_all=floras.island_id.nunique()
    if n_all!=4505:
        raise ValueError("Original all 4505 island botanical parent is mandatory")
    loc={name:i for i,name in enumerate(ids)}
    if set(floras.island_id)!=set(ids):
        raise ValueError("All GIFT island names and GIS id alignment failed")
    tmp=sel[["island_id","accepted_species"]].copy()
    tmp["g"]=tmp.accepted_species.map(genus)
    tmp["rid"]=tmp.island_id.map(loc)
    genders=sorted(tmp.g.unique())
    lookup={g:j for j,g in enumerate(genders)}
    count=tmp.groupby(["rid","g"]).size().reset_index(name="species")
    nrows=len(ids)
    n_species=sel.groupby("island_id").size().reindex(ids,fill_value=0).to_numpy(int)
    ngen=count.groupby("rid").size().reindex(range(nrows),fill_value=0).to_numpy(int)
    row=count.rid.to_numpy(int)
    col=count.g.map(lookup).to_numpy(int)
    denom=np.maximum(n_species[row],1)
    matrix=csr_matrix((count.species.to_numpy(float)/denom,(row,col)),
        shape=(nrows,len(genders)),dtype=float)
    return matrix,n_species,ngen,len(genders)


def _jaccard(a:set[str],b:set[str])->float:
    union=a|b
    return float(len(a&b)/len(union)) if union else 0.


def _nearest_probabilities(genera:csr_matrix,valid_source:np.ndarray,
                           source_score:np.ndarray,p_region:np.ndarray,
                           *,k:int=K)->np.ndarray:
    if k<2 or len(valid_source)<2:raise ValueError("Need ≥2 donor islands")
    scores=np.asarray(source_score,float)
    if len(scores)!=len(valid_source) or not np.isfinite(scores).all():
        raise ValueError("Unusable geographic source neighbor scores")
    take=np.argsort(scores,kind="stable")[:min(k,len(scores))]
    d=scores[take]
    half=max(float(np.median(d)),.01)
    w=np.exp(-.5*(d/half)**2)
    w/=w.sum()
    local=np.asarray(genera[valid_source[take]].multiply(
         w[:,None]).sum(axis=0)).ravel()
    p=SHRINK*p_region+(1.-SHRINK)*local
    if (not np.isfinite(p).all() or not np.isclose(p.sum(),1.)):
        raise AssertionError("Model source genus probability failed to normalize")
    return p


def _proper_genus_score(p:np.ndarray,truth:csr_matrix,
                        source_known:np.ndarray)->tuple[float,float]:
    """Source genus category vocabulary + ONE unknown class, no target labels."""
    if truth.shape[0]!=1 or truth.shape[1]!=len(p):
        raise ValueError("Original botanical target category comparison mismatch")
    if (not np.isfinite(p).all() or (p<0).any()
            or not np.isclose(p.sum(),1.)):
        raise ValueError("Source normalized whole-genus distribution invalid")
    q=truth.toarray().ravel()
    if not np.isclose(q.sum(),1.):
        raise ValueError("Original accepted target plant genus composition not 1")
    known=q*source_known
    extra=float(q.sum()-known.sum())
    # All source models use the same eligible donor genus dictionary; no
    # unobserved target genus gets its own learned probability mass.
    prob=(1.-DELTA_UNKNOWN)*p
    hits=float(known@prob+extra*DELTA_UNKNOWN)
    brier=float(1.+np.dot(prob,prob)+DELTA_UNKNOWN**2-2.*hits)
    logloss=float(-np.dot(known,np.log(np.maximum(prob,1e-15))) -
                  extra*np.log(DELTA_UNKNOWN))
    return brier,logloss


def _group_summary(frame:pd.DataFrame,key:str,*,seed:int,nboot:int=NBOOT)->dict:
    filt=frame.dropna(subset=[key,*OUTPUT_COLS])
    n=len(filt);nc=filt[key].nunique()
    if n<MIN_EVAL_ISLANDS or nc<MIN_CLUSTERS:
        return dict(status="NOT_IDENTIFIABLE_SOURCE_GEOGRAPHIC_GROUP_SUPPORT",
                    n_original_source_scored_floras=n,n_source_groups=nc)
    means=filt.groupby(key,sort=True)[list(OUTPUT_COLS)].mean().to_numpy()
    rng=np.random.default_rng(seed)
    draw=means[rng.integers(0,nc,size=(nboot,nc))].mean(axis=1)
    out=dict(status="ORIGINAL_GIFT_TAXONOMIC_GENUS_SOURCE_COMPOSITION_GROUP_PREDICTION",
             n_original_source_scored_floras=n,
             n_source_groups=nc,
             n_original_outcome_blind_group_resamples=nboot)
    for i,m in enumerate(OUTPUT_COLS):
        out[m+"_group_equal_mean"]=float(means[:,i].mean())
        out[m+"_group_resample_95"]=[
            float(np.quantile(draw[:,i],.025)),
            float(np.quantile(draw[:,i],.975))]
    return out


def analyze(floras:pd.DataFrame,geo:pd.DataFrame,assign:pd.DataFrame,
            *,scopes=SCOPES,nboot:int=NBOOT)->dict:
    if not {"island_id","analysis_regime","spatial_block",
      "island_latitude","island_longitude",*FEATURES}<=set(geo):
        raise ValueError("Original full GSHHG climate/area/distance missing")
    if geo.island_id.duplicated().any():raise ValueError("GSHHG duplicate island")
    ids=sorted(floras.island_id.unique())
    if len(ids)!=4505 or floras.duplicated([
        "island_id","accepted_species"]).any():
        raise ValueError("Original complete GIFT 4505 accepted botanical flora changed")
    geo=geo.set_index("island_id").reindex(ids)
    z=_floats(geo)
    lat=pd.to_numeric(geo.island_latitude,errors="coerce").to_numpy(float)
    lon=pd.to_numeric(geo.island_longitude,errors="coerce").to_numpy(float)
    coord_valid=(np.isfinite(lat)&np.isfinite(lon)&
                 (np.abs(lat)<=90)&(np.abs(lon)<=180))
    xyz=_xyz(np.where(coord_valid,lat,0),np.where(coord_valid,lon,0))
    valid=np.isfinite(z).all(axis=1)&coord_valid
    regions=geo.analysis_regime.fillna("geography_unresolved").astype(str).to_numpy()
    blocks=geo.spatial_block.fillna("").astype(str).to_numpy()
    source=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    need={"island_id","source_mode","source_rank","entity_ID"}
    if not need<=set(assign) or source.empty:
        raise ValueError("Original GIFT source-region entity labels not present")
    source.source_rank=pd.to_numeric(source.source_rank,errors="raise")
    if source.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original source geo_k5 entity ranks duplicate")
    first=(source.loc[source.source_rank.eq(1)]
           .drop_duplicates("island_id").set_index("island_id").entity_ID.to_dict())
    fsource=np.array([str(first.get(x) or "") for x in ids])
    fifth=(source.loc[source.source_rank.le(5)]
           .groupby("island_id").entity_ID.agg(set).to_dict())
    candidate_sets=[fifth.get(x,set()) for x in ids]
    valid=valid&(blocks!="")&(fsource!="")&np.isin(regions,REGIONS)
    all_dims=len(ids);orig_reg=Counter(regions)
    if sum(orig_reg.values())!=all_dims:
        raise AssertionError("Original source all island region sum drifted")

    results=[]
    for scope in scopes:
        Q,npl,ngenus,Kall=genus_source_matrix(floras,ids,scope)
        scope_valid=valid&(npl>=MIN_FOCAL_SPECIES)&(ngenus>=MIN_FOCAL_GENERA)
        for reg in REGIONS:
            region_ids=np.where(scope_valid&(regions==reg))[0]
            if len(region_ids)==0:
                for_empty=dict(scope=scope,region=reg,
                    n_original_region_botanical_floras=int(orig_reg[reg]),
                    n_source_independently_qualified_focal_floras=0,
                    n_original_genus_source_predictions=0,
                    status="NOT_IDENTIFIABLE_NO_SOURCE_GIFT_GENUS_FLORAS")
                results.append(for_empty);continue
            scored=[]
            unsupported=Counter()
            # The SAME independent source genus universe & candidates for
            # environment, source identity and geo benchmark on each island.
            for focal in region_ids:
                donor=region_ids[
                    (blocks[region_ids]!=blocks[focal])&
                    (fsource[region_ids]!=fsource[focal])]
                if len(donor)<MIN_SOURCE_ISLANDS:
                    unsupported["n_focal_islands_without_30_distant_source_islands"]+=1
                    continue
                dsource=fsource[focal]
                b=blocks[focal]
                if (len(set(blocks[donor]))<MIN_CLUSTERS or
                    len(set(fsource[donor]))<MIN_CLUSTERS):
                    unsupported["n_focal_islands_lacking_10_donor_source_groups"]+=1
                    continue
                uniform=np.asarray(Q[donor].mean(axis=0)).ravel()
                if not np.isclose(uniform.sum(),1):
                    raise AssertionError("Original source regional genus weights do not sum to 1")
                known=uniform>0
                d_env=np.linalg.norm(z[donor]-z[focal],axis=1)
                d_climate=np.linalg.norm(z[donor,:4]-z[focal,:4],axis=1)
                d_geo=np.linalg.norm(xyz[donor]-xyz[focal],axis=1)
                target_set=candidate_sets[focal]
                d_withsrc=d_env+SOURCE_WEIGHT*np.asarray(
                    [1.-_jaccard(target_set,candidate_sets[int(i)])
                     for i in donor],dtype=float)
                models={
                    "geo_nearest":d_geo,
                    "climate_nearest":d_climate,
                    "environment_nearest":d_env,
                    "environment_source_nearest":d_withsrc,
                }
                truth=Q.getrow(int(focal))
                if npl[focal]<MIN_FOCAL_SPECIES:raise AssertionError("Missing original target")
                base=_proper_genus_score(uniform,truth,known)
                rec=dict(spatial_block=blocks[focal],first_source=fsource[focal],
                    n_source_donor_island_floras=int(len(donor)),
                    n_original_botanical_flora_species=int(npl[focal]),
                    n_original_botanical_flora_genera=int(ngenus[focal]),
                    original_genus_not_seen_by_ANY_source_target_species_fraction=float(
                        truth.toarray().ravel()[~known].sum()))
                for model,dist in models.items():
                    p=_nearest_probabilities(Q,donor,dist,uniform,k=K)
                    compared=_proper_genus_score(p,truth,known)
                    rec[model+"_minus_region_equal_Brier"]=base[0]-compared[0]
                    rec[model+"_minus_region_equal_logloss"]=base[1]-compared[1]
                scored.append(rec)
            frame=pd.DataFrame(scored)
            row=dict(scope=scope,region=reg,
                n_original_region_botanical_floras=int(orig_reg[reg]),
                n_original_source_genus_names_global_vocab=Kall,
                n_original_source_islands_with_full_geo_and_genus_baseline=int(len(region_ids)),
                n_original_floras_scored_for_genus_community_environmental_prediction=len(scored),
                n_unscored_donor_coverage_by_reason=dict(unsupported),
                n_original_focal_islands_no_complete_geographic_or_full_genus_record=int(
                    orig_reg[reg]-len(region_ids)),
                n_genus_overlap_heldout_group_source_PSEUDOABSENCE_not_inferred=True)
            if scored:
                row["original_source_unknown_genus_species_fraction_mean"]=float(
                    frame.original_genus_not_seen_by_ANY_source_target_species_fraction.mean())
                row["original_mean_target_accepted_species"]=float(
                    frame.n_original_botanical_flora_species.mean())
                row["original_mean_target_genera"]=float(
                    frame.n_original_botanical_flora_genera.mean())
                public_groups={
                    "spatial_block":"physical_10degree_block",
                    "first_source":"first_candidate_source_entity_group",
                }
                for key,outkey in public_groups.items():
                    row[outkey]=_group_summary(frame,key,
                        seed=int.from_bytes(hashlib.sha256(
                            (scope+reg+key).encode()).digest()[:8],"big"),
                        nboot=nboot)
                row["status"]=(
                    "ORIGINAL_GIFT_GENUS_OCCUPANCY_COMPOSITION_ENVIRONMENTAL_TEST_BOTH_SOURCE_GROUPS"
                    if all(row[outkey]["status"]==
                    "ORIGINAL_GIFT_TAXONOMIC_GENUS_SOURCE_COMPOSITION_GROUP_PREDICTION"
                    for outkey in public_groups.values()) else
                    "NOT_IDENTIFIABLE_INDEPENDENT_SOURCE_GEOGRAPHIC_GROUPING")
            else:
                row["status"]="NOT_IDENTIFIABLE_NO_SOURCE_DONOR_MATCH"
            results.append(row)
    return dict(
      status="ORIGINAL_GIFT_4505_NO_FOCAL_GENUS_SOURCE_ISLAND_ENVIRONMENTAL_GENUS_ASSEMBLY",
      original_GIFT_parent_island_floras=4505,
      n_full_source_geo_valid_islands=int(valid.sum()),
      source_geo_valid_by_original_region={
        r:int(sum(valid&(regions==r))) for r in REGIONS},
      source_region_parent_island_counts=dict(orig_reg),
      source_genus_models=list(MODELS),
      fixed_exposure_features=list(FEATURES),
      source_original_40_nearest_donor_island_count=K,
      source_original_environmental_feature_global_IQR_transductive_without_taxa=True,
      original_no_oracle_original_island_genus_composition_source_predictions=results,
      controls=dict(
        original_4505_island_botanical_flora_population_retained=True,
        original_focal_true_genus_community_not_used_to_predict=True,
        original_true_accepted_species_genus_used_ONLY_FOR_SCORING=True,
        original_training_GIFT_islands_in_focal_10degree_block_EXCLUDED=True,
        original_training_GIFT_islands_with_same_first_GIFT_candidate_mainland_source_EXCLUDED=True,
        GIFT_source_geo_k5_ranked_entity_overlap_is_NOT_measured_seed_dispersal=True,
        GIFT_source_donor_island_taxonomy_NOT_independent_mainland_ancestor=True,
        regional_genus_frequency_baseline_uses_same_eligible_source_islands=True,
        source_species_normalized_by_island_and_UNSEEN_GENUS_held_in_common_unknown_category=True,
        source_taxon_genus_nonlisting_is_NOT_confirmed_ecological_absence=True,
        no_true_island_genus_colonization_probability_OR_independent_trait_adaptation=True,
        no_pollinator_guild_selection_or_seedset_causality_claim=True,
        source_labels_and_true_focal_island_species_names_NEVER_EXPORTED=True,
        source_islands_not_independent_evolutionary_replicates=True,
        source_group_bootstrap_exploratory_not_multiplicity_adjusted=True,
      ),
    )


def main():
    p=argparse.ArgumentParser()
    for key in ("original-root","geography","output"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--nboot",type=int,default=NBOOT)
    a=p.parse_args()
    paths=dict(floras=a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
      assign=a.original_root/"fixed/source/result/island_source_assignments.csv.gz")
    expected=dict(floras="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
      assign="8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44")
    for k,q in paths.items():
        if sha(q)!=expected[k]:
            raise ValueError("Original SHA-locked GIFT 4505 or geography donor opportunity changed: "+k)
    flora=pd.read_csv(paths["floras"],dtype=str,
        usecols=["island_id","accepted_species","origin_status"])
    geo=pd.read_csv(a.geography,dtype=str,
        usecols=["island_id","analysis_regime","spatial_block",
            "island_latitude","island_longitude",*FEATURES])
    if geo.island_id.nunique()!=8264:
        raise ValueError("Independent original GSHHG geography source universe changed")
    source=pd.read_csv(paths["assign"],dtype=str,
        usecols=["island_id","source_mode","source_rank","entity_ID"])
    result=analyze(flora,geo,source,nboot=a.nboot)
    raw=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    for private in (
        "island_id","accepted_species","entity_ID",
        "genus","genus_name","source_genus","original_genus",
        "latitude","longitude","island_latitude","island_longitude",
        "spatial_block","source_cluster","source_region_id"):
        if '"'+private+'"' in raw:
            raise ValueError("Original source GIFT island accepted taxonomy or geography leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    for row in result["original_no_oracle_original_island_genus_composition_source_predictions"]:
        print("ANONYMOUS_GIFT_ENV_GENUS_OCCUPANCY="+json.dumps({
          k:row.get(k) for k in ("scope","region","status",
            "n_original_region_botanical_floras",
            "n_original_source_islands_with_full_geo_and_genus_baseline",
            "n_original_floras_scored_for_genus_community_environmental_prediction",
            "n_original_focal_islands_no_complete_geographic_or_full_genus_record",
            "original_source_unknown_genus_species_fraction_mean",
            "physical_10degree_block","first_candidate_source_entity_group")
        },sort_keys=True))
    print("ANONYMOUS_GIFT4505_ENV_GENUS_SUMMARY="+json.dumps(dict(
       status=result["status"],parent_original_GIFT_floras=4505,
       n_original_geo_complete_GIFT_islands=result["n_full_source_geo_valid_islands"],
       first_source_candidate_group_taxonomic_NOT_ancestral=True,
       no_local_genus_roster_predictor=True,
       model_labels=list(MODELS)),sort_keys=True))
if __name__=="__main__":
    main()
