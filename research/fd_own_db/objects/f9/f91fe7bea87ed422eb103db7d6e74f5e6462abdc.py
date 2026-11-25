#!/usr/bin/env python3
"""Downscaled native bat ranges x historical public floral state-set aggregates.
Exploratory source pilot; no plant-genus null, realized pollination, or absence outside range product.
"""
import argparse, hashlib, io, json, pathlib, re, subprocess, tempfile, time
import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Polygon, box
from shapely.affinity import translate
from tropical_pollinator_public_flora_bridge import BASE, GEO_SHA, gshhg
from tropical_pollinator_source_preflight import get

RDS_URL="https://zenodo.org/records/10974868/files/Chiroptera.rds?download=1"
RDS_MD5="493d2e09b0810a5291326c55892526d3"
GSHHG_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
DIET_URL="https://figshare.com/ndownloader/files/5631084"
PINS={"direct":"6c834cb58105993b3bda10b6816198a6023a25ff4ce407eebccee5b93d5d1c7f",
      "direct_single_state":"ebb615091b971606a6626a488c44322340e6f0526c330b86b88ffee287eb9bc2"}
OUTCOMES=("white","white_open_brush","red_pink","tubular")
SEED=20261125

def sha(b):return hashlib.sha256(b).hexdigest()
def norm(x):return " ".join(str(x).replace("_"," ").split()).casefold()
def download(url,cap):
    for trial in range(3):
        try:return get(url,maxbytes=cap)
        except Exception:
            if trial==2:raise
            time.sleep(5)

def diet():
    raw=download(DIET_URL,5_000_000)
    d=pd.read_csv(io.BytesIO(raw),sep="\t",encoding="latin-1",low_memory=False)
    cols={re.sub("[^a-z0-9]","",x.lower()):x for x in d.columns}
    assert {"scientific","dietnect","mswfamilylatin"}<=set(cols)
    t=pd.DataFrame({"name":d[cols["scientific"]].map(norm),
                    "nectar":pd.to_numeric(d[cols["dietnect"]],errors="coerce"),
                    "family":d[cols["mswfamilylatin"]]})
    t=t[t.nectar.between(0,100)&t.name.ne("nan")&t.name.ne("")].copy()
    bad=t.groupby("name").nectar.nunique()
    t=t[~t.name.isin(bad[bad>1].index)].drop_duplicates("name")
    assert len(t)>4000
    return t.set_index("name"),{"sha256":sha(raw),"valid_species":len(t),"conflicting_names":int((bad>1).sum())}

def polygons(g):
    if g is None or g.is_empty:return []
    if g.geom_type=="Polygon":return [g]
    if g.geom_type in ("MultiPolygon","GeometryCollection"):
        return [q for x in g.geoms for q in polygons(x)]
    return []

def geographic_parts(poly):
    """Split genuinely dateline-crossing rings before planar projection.

    A wide bounding box alone is not evidence of a crossing. Only a longitude
    jump >180 degrees triggers shortest-arc unwrapping. Align holes with the
    exterior and conserve the repaired polygon's planar degree-area when
    clipping. Nonclosing/winding rings fail rather than invent coverage.
    """
    rings=[np.asarray(poly.exterior.coords,float)[:,:2]]+[np.asarray(h.coords,float)[:,:2] for h in poly.interiors]
    for a in rings:
        if not np.isfinite(a).all() or np.any(np.abs(a[:,0])>180+1e-7) or np.any(np.abs(a[:,1])>90+1e-7):
            raise ValueError("Out-of-domain geographic coordinates")
    crossing=any(np.any(np.abs(np.diff(a[:,0]))>180+1e-8) for a in rings)
    wide=poly.bounds[2]-poly.bounds[0]>180
    if not crossing:
        return polygons(shapely.make_valid(poly)),False,wide
    def unwrap(a):
        b=a.copy()
        b[:,0]=np.rad2deg(np.unwrap(np.deg2rad(a[:,0])))
        if not np.allclose(b[0],b[-1],rtol=0,atol=1e-7):
            raise ValueError("Nonclosing longitude ring: ambiguous global winding")
        b[-1]=b[0]
        if np.ptp(b[:,0])>=360-1e-7:
            raise ValueError("Ambiguous global longitude extent")
        return b
    ext=unwrap(rings[0])
    center=(ext[:,0].min()+ext[:,0].max())/2
    holes=[]
    for a in rings[1:]:
        h=unwrap(a)
        h[:,0]+=360*round((center-(h[:,0].min()+h[:,0].max())/2)/360)
        holes.append(h)
    shell=Polygon(ext)
    if any(not shell.covers(Polygon(h).representative_point()) for h in holes):
        raise ValueError("Dateline hole could not be aligned")
    original=polygons(shapely.make_valid(Polygon(ext,holes)))
    pieces=[]
    for g in original:
        lo=int(np.floor((g.bounds[0]+180)/360));hi=int(np.floor((g.bounds[2]+180)/360))
        for band in range(lo,hi+1):
            clipped=g.intersection(box(-180+360*band,-90,180+360*band,90))
            pieces.extend(translate(q,xoff=-360*band) for q in polygons(clipped) if q.area>0)
    before=sum(q.area for q in original);after=sum(q.area for q in pieces)
    if before<=0 or not np.isclose(before,after,rtol=1e-8,atol=1e-10):
        raise ValueError("Dateline clipping failed area conservation")
    if any(q.bounds[0]<-180-1e-7 or q.bounds[2]>180+1e-7 for q in pieces):
        raise ValueError("Dateline clipping escaped longitude domain")
    return pieces,True,wide

def geometry_self_test():
    # Narrow crossing rectangle should not cover Greenwich. Preserve a hole
    # crossing the same seam, and distinguish genuinely wide no-jump rings.
    p=Polygon([(179,-2),(-179,-2),(-179,2),(179,2),(179,-2)],
              [[(179.5,-1),(-179.5,-1),(-179.5,1),(179.5,1),(179.5,-1)]])
    g,c,w=geographic_parts(p)
    assert c and w and len(g)==2 and np.isclose(sum(x.area for x in g),6.)
    assert all(not x.intersects(box(-10,-3,10,3)) for x in g)
    assert all(x.is_valid for x in g)
    wide=Polygon([(-110,-2),(0,-2),(110,-2),(110,2),(0,2),(-110,2),(-110,-2)])
    g,c,w=geographic_parts(wide)
    assert w and not c and np.isclose(sum(x.area for x in g),880.)
    ordinary=box(30,-2,32,2)
    g,c,w=geographic_parts(ordinary)
    assert not c and not w and np.isclose(sum(x.area for x in g),8.)
    print("BAT_DATELINE_GEOMETRY_TEST_PASS",flush=True)

def ranges(tmp):
    raw=download(RDS_URL,100_000_000)
    assert hashlib.md5(raw).hexdigest()==RDS_MD5
    digest=sha(raw)
    (tmp/"bat.rds").write_bytes(raw);del raw
    subprocess.run(["Rscript",str(pathlib.Path(__file__).with_name("export_mdd_bat_ranges.R")),
                    str(tmp/"bat.rds"),str(tmp/"bat.gpkg")],check=True)
    g=gpd.read_file(tmp/"bat.gpkg")
    assert g.crs is not None and "sciname" in g and g.sciname.notna().all()
    g=g.to_crs(4326)
    names=g.sciname.map(norm)
    assert names.nunique()>500 and names.ne("").all()
    # Inspect raw rings before make_valid, which could otherwise repair a
    # crossing ring into a perfectly valid polygon on the wrong hemisphere.
    pieces=[];crossing=wide_no_jump=0;crossing_species=set()
    for s,geom in zip(names,g.geometry):
        for q in polygons(geom):
            try: parts,changed,wide=geographic_parts(q)
            except ValueError as e:
                raise RuntimeError("Ambiguous range geometry; source not scored: "+str(e)) from e
            crossing+=int(changed);wide_no_jump+=int(wide and not changed)
            if changed:crossing_species.add(s)
            pieces.extend((s,z) for z in parts if z.area>0)
    print("BAT_GEOMETRY_AUDIT",json.dumps({"crossing_polygons":crossing,
        "crossing_species":len(crossing_species),"wide_no_jump":wide_no_jump}),flush=True)
    p=gpd.GeoDataFrame(pieces,columns=["name","geometry"],geometry="geometry",crs=4326).to_crs(6933)
    assert p.name.nunique()==names.nunique()
    p.geometry=shapely.make_valid(p.geometry.to_numpy())
    assert np.isfinite(shapely.area(p.geometry.to_numpy())).all()
    return p,{"rds_md5":RDS_MD5,"rds_sha256":digest,"map_species":int(names.nunique()),"polygon_pieces":len(p),"crossing_polygons":crossing,"crossing_species":len(crossing_species),"wide_no_jump":wide_no_jump}

def incidence(units,p,d,strict_family_map=None,exact_diet=None):
    u=units.to_crs(6933).reset_index(drop=True)
    arr=shapely.make_valid(u.geometry.to_numpy())
    area=shapely.area(arr)
    assert (area>0).all() and np.isfinite(area).all()
    tree=shapely.STRtree(arr)
    # columns: non-nectar, partial, specialist, unmatched; distinct map species.
    counts={t:np.zeros((len(u),4),dtype=int) for t in (.1,.5)}
    from bat_family_context import allocate, add_counts, attach_counts
    families={t:allocate(len(u)) for t in counts}
    from bat_family_local_completeness import allocate as local_allocate, add as local_add, attach as local_attach
    local=None
    if strict_family_map is not None:
        assert exact_diet is not None
        local={t:{z:local_allocate(len(u)) for z in ('local_exact','local_mdd')} for t in counts}
        exact_unknown={t:np.zeros(len(u),dtype=int) for t in counts}
    occurrence={t:set() for t in counts}
    matched=set(p.name)&set(d.index)
    for s,g in p.groupby("name",sort=True):
        q=shapely.union_all(g.geometry.to_numpy())
        hit=tree.query(q,predicate="intersects")
        if not len(hit):continue
        cover=np.clip(shapely.area(shapely.intersection(arr[hit],q))/area[hit],0,1)
        k=3 if s not in matched else 2 if d.loc[s,"nectar"]>=50 else 1 if d.loc[s,"nectar"]>0 else 0
        for t,c in counts.items():
            selected=hit[cover>=t-1e-12]
            c[selected,k]+=1
            if k!=3:add_counts(families[t],d.loc[s,"family"],k,selected)
            if local is not None:
                ke=3 if s not in exact_diet.index else 2 if exact_diet.loc[s,'nectar']>=50 else 1 if exact_diet.loc[s,'nectar']>0 else 0
                family=strict_family_map.get(s)
                local_add(local[t]['local_exact'],family,ke,selected)
                local_add(local[t]['local_mdd'],family,k,selected)
                if ke==3:exact_unknown[t][selected]+=1
            if len(selected):occurrence[t].add(s)
    tables={};summary={}
    for t,c in counts.items():
        v=pd.DataFrame(c,columns=["other","partial","specialist","unknown"])
        v["island_id"]=u.island_id.to_numpy()
        v["range_area"]=area/1e6
        v["total"]=c.sum(1)
        v["nectar"]=v.partial+v.specialist
        attach_counts(v,families[t])
        if local is not None:
            for prefix,a in local[t].items():local_attach(v,a,prefix)
            v['local_exact_global_unknown']=exact_unknown[t]
            v['local_mdd_global_unknown']=v.unknown
        # Unmatched taxa remain unknown, never assigned zero nectar.
        v["complete"]=v.unknown.eq(0)&v.total.gt(0)
        v["match_share"]=np.divide((v.total-v.unknown).to_numpy(),v.total.to_numpy(),out=np.full(len(v),np.nan),where=v.total.gt(0).to_numpy())
        tables[t]=v
        summary[str(t)]={"map_species_with_island_overlap":len(occurrence[t]),
                         "matched_species_with_overlap":len(occurrence[t]&matched),
                         "potential_bat_positive_islands":int(v.total.gt(0).sum()),
                         "potential_nectar_positive_lower_bound":int(v.nectar.gt(0).sum()),
                         "complete_bat_diet_positive_islands":int(v.complete.sum()),
                         "coverage_ge95_bat_positive":int(v.match_share.ge(.95).sum()),
                         "unknown_species_island_pairs":int(v.unknown.sum())}
    return tables,summary|{"exact_map_species_diet_matches":len(matched),"map_species_unmatched":int(p.name.nunique()-len(matched))}

def flora():
    raw=download(BASE+"floral_space_20261106/all/profile_coordinates.csv",2_000_000)
    ref=pd.read_csv(io.BytesIO(raw))
    assert ref.profile_id.is_unique
    color=ref.flower_primary_color.str.split("|",regex=False).map(set)
    form=ref.floral_form.str.split("|",regex=False).map(set)
    assert any("white" in x for x in color)
    ref["white"]=color.map(lambda x:"white" in x)
    ref["white_open_brush"]=ref.white&form.map(lambda x:bool(x&{"open_radial","brush_puff"}))
    ref["red_pink"]=color.map(lambda x:"red_pink" in x)
    ref["tubular"]=form.map(lambda x:"tubular" in x)
    out={}
    for scope,pin in PINS.items():
        b=download(BASE+"floral_space_20261106/"+scope+"/island_profile_carriers.csv.gz",2_000_000)
        assert sha(b)==pin
        c=pd.read_csv(io.BytesIO(b),compression="gzip")
        assert not c.duplicated(["island_id","profile_id"]).any() and c.n_species.gt(0).all()
        j=c.merge(ref[["profile_id",*OUTCOMES]],on="profile_id",validate="many_to_one")
        assert len(j)==len(c)
        for o in OUTCOMES:j[o]=j[o].astype(int)*j.n_species
        v=j.groupby("island_id")[["n_species",*OUTCOMES]].sum().reset_index()
        for o in OUTCOMES:v[o]/=v.n_species
        out[scope]=v
    return out,sha(raw)

def estimate(d,features,controls,draws=1999):
    cols=features+controls
    d=d.dropna(subset=["block",*cols,*OUTCOMES]).copy()
    blocks=[]
    for _,v in d.groupby("block",sort=True):
        if len(v)<2 or v[features[0]].nunique()<2:continue
        x=v[cols].to_numpy(float); y=v[list(OUTCOMES)].to_numpy(float)
        if not(np.isfinite(x).all() and np.isfinite(y).all()):continue
        blocks.append(((x-x.mean(0))/np.sqrt(len(x)),(y-y.mean(0))/np.sqrt(len(x)),len(x)))
    k=len(blocks);n=sum(b[2] for b in blocks)
    base={"eligible":len(d),"informative_blocks":k,"informative_islands":n,"features":features,"controls":controls,
          "nectar_positive_eligible":int(d[features[0]].gt(0).sum())}
    def fit(w):
        active=np.flatnonzero(w)
        if not len(active) or sum(blocks[i][2] for i in active)-len(active)-len(cols)<4:return None
        x=np.vstack([blocks[i][0]*np.sqrt(w[i]) for i in active])
        y=np.vstack([blocks[i][1]*np.sqrt(w[i]) for i in active])
        if np.linalg.matrix_rank(x)<len(cols):return None
        b=np.linalg.lstsq(x,y,rcond=None)[0]
        return b[0]-b[1]
    v=fit(np.ones(k,dtype=int))
    if v is None:return base|{"status":"not_estimable"}
    rng=np.random.default_rng(SEED)
    boot=[];leave=[]
    for _ in range(draws):
        z=fit(np.bincount(rng.integers(0,k,k),minlength=k))
        if z is not None:boot.append(z)
    for i in range(k):
        w=np.ones(k,dtype=int);w[i]=0;z=fit(w)
        if z is not None:leave.append(z)
    gate=k>=30 and n>=200 and len(boot)>=.9*draws and len(leave)==k
    res={}
    for j,o in enumerate(OUTCOMES):
        lv=np.asarray(leave)[:,j] if leave else []
        res[o]={"contrast":float(v[j]),
                "group_bootstrap_interval":np.quantile(np.asarray(boot)[:,j],[.025,.975]).tolist() if len(boot)>=.9*draws else None,
                "leave_block_range":[float(min(lv)),float(max(lv))] if len(lv) else None}
    primary=res["white"]
    return base|{"status":"exploratory_support" if gate and primary["group_bootstrap_interval"][0]>0 and primary["leave_block_range"][0]>0 else "exploratory_not_supported" if gate else "insufficient_support",
                 "valid_bootstrap":len(boot),"valid_leave_block":len(leave),"endpoints":res}

def panels(f,tab,geo):
    out={}
    for scope,v in f.items():
        for overlap,c in tab.items():
            d=v.merge(c,on="island_id",validate="one_to_one").merge(geo,on="island_id",validate="one_to_one")
            # Historical public aggregate; tropical, bat-positive, fully mapped diets only.
            d=d[d.analysis_regime.eq("tropical")&d.complete&d.n_species.ge(5)].copy()
            d=d[d.distance_to_continent_km.ge(0)&d.range_area.gt(0)].copy()
            d["block"]=(np.floor(d.island_latitude/10).astype("Int64").astype(str)+"_"+np.floor(d.island_longitude/10).astype("Int64").astype(str))
            d["log_area"]=np.log2(d.range_area)
            d["log_distance"]=np.log2(1+d.distance_to_continent_km)
            d["log_n"]=np.log2(d.n_species)
            d["abs_lat"]=d.island_latitude.abs()
            for x in ("nectar","other","partial","specialist"):d["log_"+x]=np.log1p(d[x])
            controls=["log_area","log_distance","log_n","abs_lat",*[f"climate_pc{i}" for i in range(1,5)]]
            for amin in (20,100):
                for nmin in (5,50):
                    sub=d[d.range_area.ge(amin)&d.n_species.ge(nmin)]
                    for label,features in (("broad",["log_nectar","log_other"]),("split",["log_specialist","log_partial","log_other"])):
                        key=f"{scope}_cover{overlap}_area{amin}_n{nmin}_{label}"
                        out[key]=estimate(sub,features,controls)
    return out

def self_test():
    geometry_self_test()
    from bat_family_context import self_test as family_self_test
    family_self_test()
    from bat_family_local_completeness import self_test as local_self_test
    local_self_test()
    rng=np.random.default_rng(SEED);n=240
    d=pd.DataFrame({"block":np.repeat(np.arange(30),8),"a":np.tile(np.arange(8)/7,30),"b":rng.normal(size=n),"c":rng.normal(size=n)})
    for o in OUTCOMES:d[o]=.5*d.a-.1*d.b+rng.normal(0,.001,n)
    x=estimate(d,["a","b"],["c"],draws=199)
    assert x["status"]=="exploratory_support" and abs(x["endpoints"]["white"]["contrast"]-.6)<.01
    d["a"]=0
    assert estimate(d,["a","b"],["c"],draws=199)["status"]=="not_estimable"
    print("BAT_SYNTHETIC_ESTIMATOR_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");p.add_argument("--shared-support",action="store_true");p.add_argument("--synonym-sensitivity",action="store_true");p.add_argument("--unknown-diet-sensitivity",action="store_true")
    a=p.parse_args()
    if a.self_test:self_test();return
    if a.out is None:p.error("--out required")
    a.out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        tmp=pathlib.Path(td)
        d,dm=diet();p,pm=ranges(tmp)
        raw=download(BASE+"geography_20260924/corrected_geography_covariates.csv",8_000_000)
        assert sha(raw)==GEO_SHA
        geo=pd.read_csv(io.BytesIO(raw))
        assert len(geo)==8264 and geo.island_id.is_unique
        units,gh=gshhg(tmp,set(geo.island_id))
        assert gh==GSHHG_SHA
        tab,counts=incidence(units[units.area_km2.ge(20)],p,d)
        source={"range":pm,"diet":dm,"geo_exact_ids":len(units),"source_islands_ge20":int(units.area_km2.ge(20).sum()),"coverage":counts}
        print("BAT_SOURCE",json.dumps(source),flush=True)
        (a.out/"source_summary.json").write_text(json.dumps(source,indent=2,allow_nan=False))
        f,fs=flora()
        res=panels(f,tab,geo)
        out={"source":source,"profile_sha256":fs,"primary_panel":"direct_cover0.5_area20_n5_broad",
             "primary_endpoint":"white","panels":res,
             "warning":"Downscaled expert native range overlap, not observed or confirmed current absence; 1 km simplification particularly limits small-island inference. Strict exact diet matches in bat-positive islands only. Public historical species-name category sets, not true flower dimensions, nocturnal traits, or plant genus-adjusted effects. All panels exploratory, 10-degree physical blocks are not archipelagos. No new causal claim."}
        (a.out/"bat_public_flora.json").write_text(json.dumps(out,indent=2,allow_nan=False))
        for k,v in res.items():print("BAT_PANEL",k,json.dumps(v),flush=True)
        if a.shared_support:
            from multiguild_shared_support import same_cohort_support
            shared=same_cohort_support(tab,geo)
            (a.out/"multiguild_shared_support.json").write_text(json.dumps(shared,indent=2,allow_nan=False))
            print("MULTIGUILD_SHARED_SUPPORT",json.dumps(shared,allow_nan=False),flush=True)
        if a.unknown_diet_sensitivity:
            from bat_unknown_diet_sensitivity import missingness_panels
            m0=missingness_panels(f,tab,geo)
            (a.out/"bat_unknown_diet_sensitivity.json").write_text(json.dumps(m0,indent=2,allow_nan=False))
            for k,v in m0["panels"].items():print("BAT_UNKNOWN_PANEL",k,json.dumps(v),flush=True)
        if a.synonym_sensitivity:
            from bat_mdd_synonym_audit import build
            aug,audit,family_map,family_audit=build(set(p.name),d,return_taxonomy=True)
            print("BAT_MDD_CROSSWALK",json.dumps(audit,allow_nan=False),flush=True)
            (a.out/"mdd_crosswalk_audit.json").write_text(json.dumps(audit,indent=2,allow_nan=False))
            tab2,counts2=incidence(units[units.area_km2.ge(20)],p,aug,strict_family_map=family_map,exact_diet=d)
            print("BAT_MDD_COVERAGE",json.dumps(counts2,allow_nan=False),flush=True)
            for t in tab:
                old=tab[t].set_index("island_id");new=tab2[t].set_index("island_id").loc[old.index]
                assert np.array_equal(old.total,new.total)
                assert np.array_equal(old.unknown,new.local_exact_global_unknown)
                assert np.array_equal(new.unknown,new.local_mdd_global_unknown)
                assert (new.unknown<=old.unknown).all() and (new.nectar>=old.nectar).all()
                assert (new.loc[old.complete,["other","partial","specialist","unknown"]]==old.loc[old.complete,["other","partial","specialist","unknown"]]).all().all()
                assert new.loc[old.complete,"complete"].all()
            res2=panels(f,tab2,geo)
            from bat_family_context import clade_panels
            family_results={"exact":clade_panels(f,tab,geo,estimate),"mdd_singleton":clade_panels(f,tab2,geo,estimate)}
            (a.out/"bat_family_context.json").write_text(json.dumps(family_results,indent=2,allow_nan=False))
            print("BAT_FAMILY_CONTEXT",json.dumps(family_results,allow_nan=False),flush=True)
            from bat_family_local_completeness import panels as local_panels
            local_results=local_panels(f,tab2,geo,estimate)
            local_results['taxonomy']=family_audit
            (a.out/'bat_family_local_completeness.json').write_text(json.dumps(local_results,indent=2,allow_nan=False))
            print('BAT_FAMILY_LOCAL_COMPLETENESS',json.dumps(local_results,allow_nan=False),flush=True)
            from bat_synonym_cohort_diagnostic import diagnostic
            cohort=diagnostic(f,tab,tab2,geo,estimate)
            for t in (.5,.1):
                for ck,pr in (("original",res),("expanded_all",res2)):
                    expected=pr[f"direct_cover{t}_area20_n5_broad"]
                    found=cohort["panels"][str(t)]["fits"][ck]
                    assert found["eligible"]==expected["eligible"]
                    assert found["informative_blocks"]==expected["informative_blocks"]
                    assert found["informative_islands"]==expected["informative_islands"]
                    assert found["status"]==expected["status"]
                    if "endpoints" in expected:
                        for endpoint in expected["endpoints"]:
                            assert np.isclose(found["endpoints"][endpoint]["contrast"],expected["endpoints"][endpoint]["contrast"],rtol=1e-9,atol=1e-10)
            (a.out/"bat_synonym_cohort_diagnostic.json").write_text(json.dumps(cohort,indent=2,allow_nan=False))
            print("BAT_SYNONYM_COHORT_DIAGNOSTIC",json.dumps(cohort,allow_nan=False),flush=True)
            # Same original fully mapped cohort has exactly the same counts.
            # Expanded complete cohorts are a sensitivity, not a fresh test.
            summary={"crosswalk":audit,"coverage":counts2,"panels":res2,
                     "warning":"Post-pilot nomenclatural sensitivity; original exact-match primary unchanged. Expanded complete cohort is not independent confirmation. Singleton synonym linkage is not proof of equal historical circumscription."}
            (a.out/"bat_mdd_synonym_sensitivity.json").write_text(json.dumps(summary,indent=2,allow_nan=False))
            for k,v in res2.items():print("BAT_MDD_PANEL",k,json.dumps(v),flush=True)
            if a.shared_support:
                m2=same_cohort_support(tab2,geo)
                (a.out/"multiguild_mdd_sensitivity.json").write_text(json.dumps(m2,indent=2,allow_nan=False))
                print("MULTIGUILD_MDD_SUPPORT",json.dumps(m2,allow_nan=False),flush=True)
if __name__=="__main__":main()
