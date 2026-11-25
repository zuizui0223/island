"""DoPI v2 exact recorded study × geographic site × year/month flower visits.

Original GIFT accepted island plant names and original DoPI site keys remain
INSIDE protected island CI. Two guilds must be positively documented at the
SAME published location coordinate/grid and observation year (optionally
month). Missing site/time NEVER becomes an absence or a shared site.

This is a stringent observation-support preflight and external European
visitor/flower morphology sensitivity. It cannot establish local tropical
island guild shifts or plant reproductive fitness.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from hashlib import sha256
import math
import numpy as np

from .dopi_gift_guild_motif_bridge import (
    _study, _s, GUILDS,PAIRS,botanical_genus,
)

RESOLUTIONS=("same_study_site_year","same_study_site_year_month")
EVIDENCE=("positive_visitation_123","pollen_supported_12")
MIN_STUDIES=10
MIN_MOTIF_STUDIES=10
MIN_PLANTS_PER_GUILD=2
MIN_COMMON_GENERA=2


def _site(row:dict)->tuple|None:
    """Published documented site locator, never fabricated from country alone."""
    country=_s(row.get("country")).casefold()
    gridletter=_s(row.get("grid_letter")).upper()
    gridcode=_s(row.get("grid_code"))
    if gridletter and gridcode and gridcode.isdecimal():
        # National Grid code may span different sizes; same code is NOT a
        # verified individual field site. Source claims only shared grid cell.
        return ("published_british_grid",country,gridletter,gridcode)
    lat=_s(row.get("latitude"))
    lon=_s(row.get("longitude"))
    try:
        a,b=float(lat),float(lon)
    except (ValueError,TypeError):
        return None
    if not (math.isfinite(a) and math.isfinite(b)
            and -90<=a<=90 and -180<=b<=180):
        return None
    if not (lat.count(".")==1 and lon.count(".")==1
            and len(lat.split(".")[1])>=2
            and len(lon.split(".")[1])>=2):
        # Reject low precision country/centroid coordinates as "same site".
        return None
    return ("published_exact_coordinate_string",country,lat,lon)


def _time(row:dict)->tuple[str,str|None]|None:
    y=_s(row.get("year"))
    if not y.isdecimal() or len(y)!=4 or not 1800<=int(y)<=2026:
        return None
    m=_s(row.get("month"))
    month=str(int(m)) if m.isdecimal() and 1<=int(m)<=12 else None
    return (y,month)


def site_time_positive_index(rows, *,fields:set[str],island_names:set[str],
                             motifs:dict,classifier)->dict:
    """Stream original 2026 DoPI CSV; retain only private positive visitor units.

    Site keys contain published coordinates/grid and study IDs, and MUST NOT
    be serialized in final output. Only site_period_audit is source safe.
    """
    required={"plant_species","pollinator_species","pollinator_family",
              "pollinator_order","doi","publication","pollination",
              "country","grid_letter","grid_code","latitude","longitude",
              "year","month"}
    if missing:=required-fields:
        raise ValueError("DoPI site/year source fields unavailable: "+repr(sorted(missing)))
    if not island_names or not motifs:
        raise ValueError("original GIFT floristic source missing")
    documented=set().union(*(set(m["known"]) for m in motifs.values()))
    units=defaultdict(set)
    metrics=Counter()
    source_pub=set()
    for record in rows:
        metrics["n_original_source_rows"]+=1
        plant=_s(record.get("plant_species"))
        if plant not in island_names or plant not in documented:
            continue
        metrics["n_original_island_Direct_plant_rows"]+=1
        paper=_study(record)
        site=_site(record)
        when=_time(record)
        if not paper or site is None or when is None:
            metrics["n_study_site_or_year_unresolved"]+=1
            continue
        animal=classifier(record.get("pollinator_order"),
                          record.get("pollinator_family"),
                          record.get("pollinator_species"))
        if animal not in GUILDS:
            continue
        grade=_s(record.get("pollination"))
        if grade not in ("1","2","3"):
            continue
        source_pub.add(paper)
        metrics["n_true_positive_site_year_rows"]+=1
        year,month=when
        key=(paper,site,year)
        units[(RESOLUTIONS[0],EVIDENCE[0],key,animal)].add(plant)
        if grade in ("1","2"):
            units[(RESOLUTIONS[0],EVIDENCE[1],key,animal)].add(plant)
        if month is not None:
            monthkey=(paper,site,year,month)
            units[(RESOLUTIONS[1],EVIDENCE[0],monthkey,animal)].add(plant)
            if grade in ("1","2"):
                units[(RESOLUTIONS[1],EVIDENCE[1],monthkey,animal)].add(plant)
    return dict(private_units=units,
                n_original_source_publication_keys=len(source_pub),
                n_original_site_guild_positive_units=len(units),
                source_scan_counts=dict(metrics),
                true_guild_absence_not_inferred=True)


def _study_bootstrap(study:dict[str,list[float]], positives:set[str],
                     seed:str,draws:int)->dict:
    n=len(study)
    out=dict(n_source_publication_keys=n,
             n_source_publications_with_true_motif_positive=len(positives))
    if n<MIN_STUDIES:
        out["status"]="NOT_IDENTIFIABLE_COMMON_SITE_INDEPENDENT_STUDIES"
        return out
    if len(positives)<MIN_MOTIF_STUDIES:
        out["status"]="NOT_IDENTIFIABLE_SITE_MOTIF_POSITIVE_SUPPORT"
        return out
    # One vote per distinct source publication regardless of number of sites,
    # local positive records, flowers or observation sessions.
    v=np.asarray([float(np.mean(study[k])) for k in sorted(study)],float)
    rng=np.random.default_rng(int.from_bytes(sha256(seed.encode()).digest()[:8],"big"))
    samples=v[rng.integers(len(v),size=(draws,len(v)))].mean(axis=1)
    out.update(status="DESCRIPTIVE_POSITIVE_GUILD_COMMON_SITE_PERIOD",
               guild_a_minus_guild_b_motif_frequency=float(v.mean()),
               bootstrap95=[float(z) for z in np.quantile(samples,[.025,.975])])
    return out


def site_period_audit(index:dict, *,motifs:dict, draws:int=999)->dict:
    if draws<199:raise ValueError("insufficient study bootstrap draws")
    if not index.get("true_guild_absence_not_inferred"):
        raise ValueError("source-positive site-period index required")
    units=index["private_units"]
    groups={(r,e):sorted({key for a,b,key,g in units if (a,b)==(r,e)})
            for r in RESOLUTIONS for e in EVIDENCE}
    result=[]
    for resolution in RESOLUTIONS:
        for tier in EVIDENCE:
            keys=groups[(resolution,tier)]
            for ga,gb in PAIRS:
                common=[key for key in keys
                        if (resolution,tier,key,ga) in units
                        and (resolution,tier,key,gb) in units]
                for name,profile in motifs.items():
                    known=set(profile["known"])
                    positive=set(profile["joint"])
                    raw_studies=defaultdict(list)
                    matched_studies=defaultdict(list)
                    pos_raw=set();pos_matched=set()
                    n_site_raw=n_site_genus=0
                    for key in common:
                        a=units[(resolution,tier,key,ga)]&known
                        b=units[(resolution,tier,key,gb)]&known
                        if min(len(a),len(b))<MIN_PLANTS_PER_GUILD:
                            continue
                        n_site_raw+=1
                        paper=key[0]
                        raw_studies[paper].append(
                            len(a&positive)/len(a)-len(b&positive)/len(b))
                        if (a|b)&positive:pos_raw.add(paper)
                        ga_by=defaultdict(set);gb_by=defaultdict(set)
                        for p in a:ga_by[botanical_genus(p)].add(p)
                        for p in b:gb_by[botanical_genus(p)].add(p)
                        gen=set(ga_by)&set(gb_by)
                        if len(gen)<MIN_COMMON_GENERA:continue
                        n_site_genus+=1
                        vals=[len(ga_by[g]&positive)/len(ga_by[g])-
                              len(gb_by[g]&positive)/len(gb_by[g])
                              for g in gen]
                        matched_studies[paper].append(float(np.mean(vals)))
                        if any((ga_by[g]|gb_by[g])&positive for g in gen):
                            pos_matched.add(paper)
                    item=dict(resolution=resolution,
                              evidence_tier=tier,guild_a=ga,guild_b=gb,
                              floral_motif=name,
                              n_recorded_flower_axes=profile.get("axis_count",
                                 profile.get("n_axes",len(profile.get("trait_axes",())))),
                              n_study_site_period_units_with_both_guilds=len(common),
                              n_study_site_period_units_with_phenotype_data=n_site_raw,
                              n_study_site_period_units_with_two_shared_genera=n_site_genus)
                    item["same_site_positive_visitor"]=_study_bootstrap(
                        raw_studies,pos_raw,
                        seed=f"{resolution}|{tier}|{ga}|{gb}|{name}|raw",
                        draws=draws)
                    item["same_site_equal_genus"]=_study_bootstrap(
                        matched_studies,pos_matched,
                        seed=f"{resolution}|{tier}|{ga}|{gb}|{name}|equal",
                        draws=draws)
                    result.append(item)
    return dict(
        status="DOPI_V2_ORIGINAL_GIFT_SAME_PUBLISHED_LOCATION_PERIOD_VISITOR_SUPPORT",
        n_original_global_GIFT_botanical_islands=4505,
        n_original_DOPI_source_rows=index["source_scan_counts"].get("n_original_source_rows",0),
        n_original_source_publication_keys_with_geographic_year_evidence=index[
            "n_original_source_publication_keys"],
        n_original_positive_site_guild_units=index["n_original_site_guild_positive_units"],
        source_scan_counts=index["source_scan_counts"],
        n_visitor_pairs=len(PAIRS),
        n_direct_flower_profiles=len(motifs),
        n_site_temporal_scopes=len(RESOLUTIONS)*len(EVIDENCE),
        n_reporting_panels=len(result),
        panels=result,
        limitations=dict(
            source_geography_and_period_not_independently_verified=True,
            exact_grid_cell_not_verified_same_field_site=True,
            year_month_not_guaranteed_same_observation_session=True,
            unknown_site_or_time_not_imputed=True,
            no_undocumented_guild_absence=True,
            source_positive_records_not_local_island_visitation=True,
            same_site_source_studies_not_local_guild_preference=True,
            pollen_effectiveness_not_measured_by_grade3=True,
            no_exact_source_study_or_site_or_botanical_identifiers_released=True))
