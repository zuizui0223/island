#!/usr/bin/env python3
"""Exact public geography bridge from the 305-bee study to frozen public floral aggregates.
No private GIFT species ledger or claimed radial x deep endpoint; no effect fit in this stage.
"""
import argparse,hashlib,io,json,math,os,pathlib,subprocess,sys,tempfile,zipfile,urllib.request
from collections import defaultdict
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely import make_valid
from shapely.ops import unary_union
from shapely.geometry import Polygon,MultiPolygon
from tropical_pollinator_source_preflight import get
from tropical_three_guild_candidate_audit import fetch
from tropical_bee_bird_crosswalk import decode_source,BEE_MD5,MODEL_SHA
from fd.core import ProfileSpace,load_profiles,load_carriers,load_island_metadata
BASE="https://raw.githubusercontent.com/zuizui0223/island/d4a5dafab5f0737ea6c1c1c0263fc20decf6035c/results/"
GEO_URL="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
GEO_SHA="08218617b5880153986012e79817de65a74dd4e44b5a65e8a537dab8cb2cd0a7"

def raw(url,cap):
    return get(url,maxbytes=cap)
def extract_members(z,prefix,out):
    for m in z.infolist():
        if pathlib.PurePosixPath(m.filename).name.startswith(prefix+"."):
            (out/pathlib.PurePosixPath(m.filename).name).write_bytes(z.read(m))
def polygonal(geom):
    geom=make_valid(geom)
    if isinstance(geom,(Polygon,MultiPolygon)):return geom
    parts=[x for x in getattr(geom,"geoms",[]) if isinstance(x,Polygon)]
    return MultiPolygon(parts) if len(parts)>1 else parts[0] if parts else None

def gshhg(tmp,reference):
    p=tmp/"gshhg";p.mkdir()
    # Match original upstream geometry construction and hash definition, no approximate ID fallback.
    data=raw(GEO_URL,180_000_000)
    digest=hashlib.sha256(data).hexdigest()
    with zipfile.ZipFile(io.BytesIO(data)) as z:extract_members(z,"GSHHS_h_L1",p)
    del data
    coast=gpd.read_file(p/"GSHHS_h_L1.shp")
    if coast.crs is None:coast=coast.set_crs(4326)
    sibling=next((c for c in coast.columns if c.lower()=="sibling_id"),None)
    id_name=sibling or next((c for c in coast.columns if c.lower()=="id"),None)
    if id_name is None:raise ValueError("missing original GSHHG feature ID")
    if sibling:
        idn=pd.to_numeric(coast[sibling],errors="raise")
        assert np.isfinite(idn).all() and (idn>=0).all() and (idn%1==0).all()
        coast=coast.copy();coast[sibling]=idn.astype("int64")
    # Singleton groups still pass original unary_union/make_valid, preserving its hash semantics.
    records=[]
    for identifier,group in coast.groupby(id_name,dropna=False,sort=True):
        geom=make_valid(unary_union(group.geometry.tolist()))
        if not geom.is_empty:records.append({"source_feature_id":str(identifier),"geometry":geom})
    units=gpd.GeoDataFrame(records,geometry="geometry",crs=coast.crs).to_crs(4326)
    del coast,records
    units["geometry"]=units.geometry.map(polygonal)
    units=units.loc[units.geometry.notna() & ~units.geometry.is_empty].copy()
    units["area_km2"]=units.to_crs(6933).area/1e6
    units=units[units.area_km2.between(5,7000000)].copy()
    units["island_id"]=["gshhg_2.3.7_h_"+hashlib.sha256(g.wkb).hexdigest()[:20] for g in units.geometry]
    assert units.island_id.is_unique
    matched=units[units.island_id.isin(reference)].copy()
    print("GSHHG_REPLAY",json.dumps({"zip_sha256":digest,"reconstructed":len(units),"exact_reference_ids":len(matched),"reference":len(reference),"missing_reference":len(reference-set(matched.island_id)),"geos":shapely.geos_version_string}),flush=True)
    if len(matched)<.9*len(reference):raise RuntimeError("GSHHG replay insufficient exact reference hashes; no approximate IDs used")
    return matched,digest

def bee_geometry(tmp):
    p=tmp/"bee";p.mkdir()
    meta=json.loads(raw("https://zenodo.org/api/records/23222636",20_000_000))
    f=next(f for f in meta["files"] if f["key"]=="Data.zip")
    data=raw(f["links"].get("content") or f["links"]["self"],300_000_000)
    assert len(data)==262554015 and hashlib.md5(data).hexdigest()==BEE_MD5
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        model_raw=z.read("Data/model_table.csv")
        assert hashlib.sha256(model_raw).hexdigest()==MODEL_SHA
        model=pd.read_csv(io.StringIO(decode_source(model_raw)),dtype={"OBJECTID_1":str})
        extract_members(z,"islands",p)
    del data
    geo=gpd.read_file(p/"islands.shp",encoding="CP1252")
    key=next(c for c in geo.columns if c.lower()=="objectid_1")
    geo["bee_id"]=geo[key].astype(str).str.replace(r"\.0$","",regex=True)
    assert len(geo)==309 and len(model)==305
    ambiguous=sorted(geo.loc[geo.bee_id.duplicated(False),"bee_id"].unique())
    geo=geo[~geo.bee_id.isin(ambiguous)].copy()
    assert len(geo)==302 and geo.bee_id.is_unique
    return geo[["bee_id","geometry"]],model,ambiguous

def spatial_join(bees,units):
    # Equal-area reference in the original bee study CRS. Area intersection not centroid nearest.
    # Exclude unsupported / invalid geometry; do not force one-to-many or approximate match.
    proj=units.to_crs(bees.crs).reset_index(drop=True)
    gs=proj.geometry.to_numpy()
    good=shapely.is_valid(gs)&~shapely.is_empty(gs)&np.isfinite(shapely.area(gs))
    proj=proj[good].reset_index(drop=True);gs=proj.geometry.to_numpy()
    tree=shapely.STRtree(gs);rows=[]
    for r in bees.itertuples():
        b=r.geometry
        rec={"bee_id":r.bee_id,"status":"unresolved","best_gshhg_id":None,"bee_overlap":None,"gshhg_overlap":None}
        if b is None or b.is_empty or not b.is_valid or not math.isfinite(b.area) or b.area<=0:
            rec["status"]="invalid_bee_geometry";rows.append(rec);continue
        ids=tree.query(b)
        cand=[]
        for j in ids:
            g=gs[j]
            if g.area<=0:continue
            if min(g.area,b.area)/max(g.area,b.area)<.90:continue
            try:
                a=b.intersection(g).area
            except shapely.GEOSException:continue
            if a>0:cand.append((a/b.area,a/g.area,proj.iloc[j].island_id))
        cand.sort(key=lambda v:min(v[0],v[1]),reverse=True)
        if cand:
            u,v,id=cand[0];rec.update(best_gshhg_id=id,bee_overlap=u,gshhg_overlap=v)
            strict=[c for c in cand if min(c[:2])>=.90]
            if len(strict)==1:rec["status"]="matched"
            elif len(strict)>1:rec["status"]="ambiguous"
        rows.append(rec)
    a=pd.DataFrame(rows)
    dup=a.loc[a.status.eq("matched"),"best_gshhg_id"].value_counts()
    a.loc[a.status.eq("matched")&a.best_gshhg_id.isin(dup[dup>1].index),"status"]="reused_gshhg_id"
    return a
def flowers(tmp):
    refs=tmp/"profiles.csv"
    refs.write_bytes(raw(BASE+"floral_space_20261106/all/profile_coordinates.csv",2_000_000))
    p=ProfileSpace(load_profiles(refs));out={}
    for scope in ("direct","direct_single_state"):
        stem=BASE+"floral_space_20261106/"+scope+"/"
        ca=tmp/(scope+"_carriers.csv.gz");md=tmp/(scope+"_metrics.csv")
        ca.write_bytes(raw(stem+"island_profile_carriers.csv.gz",2_000_000))
        md.write_bytes(raw(stem+"island_metrics.csv",3_000_000))
        carriers=load_carriers(ca);meta=load_island_metadata(md)
        assert set(carriers)==set(meta)
        rows=[]
        for island,c in carriers.items():
            v=p.island(c);m=meta[island]
            assert v["n_names"]==int(m["n_species"]) and v["n_profiles"]==int(m["n_profiles"])
            rows.append({"island_id":island,"analysis_regime":m["analysis_regime"],"joint_coverage":m.get("joint_coverage",""),**v})
        out[scope]=pd.DataFrame(rows)
        print("PUBLIC_FLORAL_SOURCE",json.dumps({"scope":scope,"n":len(rows),"n50":sum(r["n_names"]>=50 for r in rows),"carrier_sha256":hashlib.sha256(ca.read_bytes()).hexdigest(),"metric_sha256":hashlib.sha256(md.read_bytes()).hexdigest()}),flush=True)
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    with tempfile.TemporaryDirectory() as td:
        tmp=pathlib.Path(td)
        ref=raw(BASE+"geography_20260924/corrected_geography_covariates.csv",8_000_000)
        assert hashlib.sha256(ref).hexdigest()==GEO_SHA,"wrong frozen public geography"
        r=pd.read_csv(io.BytesIO(ref));reference=set(r.island_id)
        assert len(reference)==8264
        units,digest=gshhg(tmp,reference)
        b,model,amb=bee_geometry(tmp);cross=spatial_join(b,units)
        cross.to_csv(a.out/"bee_gshhg_geometry_crosswalk.csv",index=False)
        m=cross[cross.status.eq("matched")].rename(columns={"best_gshhg_id":"island_id"})
        assert m.island_id.is_unique
        # Join only rights-safe public aggregate metrics; no new 4505-Direct flower claim.
        bird=fetch("bee_bird");but=fetch("butterfly")
        assert len(bird)==52 and len(but)==98
        fl=flowers(tmp);results={}
        for scope,f in fl.items():
            join=m.merge(f,on="island_id",validate="one_to_one")
            join.to_csv(a.out/(scope+"_bee_public_flower_bridge.csv"),index=False)
            panels=[]
            for label,v in [("all_matched",join),("bird_candidate",join[join.bee_id.isin(bird.bee_id)]),("butterfly_positive",join[join.bee_id.isin(but.loc[but.all_records.gt(0),"bee_id"])]),("three_guild_candidate",join[join.bee_id.isin(set(bird.bee_id)&set(but.loc[but.all_records.gt(0),"bee_id"]))])]:
                panels.append({"scope":label,"island_units":len(v),"n5":int(v.n_names.ge(5).sum()),"n50":int(v.n_names.ge(50).sum()),"regime_counts":v.analysis_regime.value_counts().to_dict(),"n50_regime_counts":v.loc[v.n_names.ge(50),"analysis_regime"].value_counts().to_dict()})
            results[scope]=panels
        out={"reference_count":len(reference),"exact_replayed_count":len(units),"gshhg_archive_sha256":digest,"bee_geometry_ambiguous_ids":amb,"crosswalk_status":cross.status.value_counts().to_dict(),"public_aggregate_flower_support":results,"warning":"Exact geography overlap gate, not bird geometry certification. Public historical PR321 colour/form name-unit aggregate only; not original 4505 radial x deep/symmetry flora. No statistical ecological effect fitted."}
        (a.out/"summary.json").write_text(json.dumps(out,indent=2))
        print("PUBLIC_FLORA_GEOGRAPHY_BRIDGE",json.dumps(out),flush=True)
if __name__=="__main__":main()
