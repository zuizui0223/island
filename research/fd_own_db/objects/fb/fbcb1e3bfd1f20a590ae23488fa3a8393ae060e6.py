#!/usr/bin/env python3
"""Generic protected range-map incidence adapter, separate from frozen PR66 bird CLI.

Inputs must be locally authorized. No acquisition, upload or public export of ranges,
island IDs, species lists or incidence matrices. All generated outputs stay protected.
A provider-specific, source-pinned codebook JSON supplies exact field names and codes.
"""
from __future__ import annotations
import argparse, hashlib, json, pathlib, re
import numpy as np
import pandas as pd
from scipy import sparse
from fd.range_map_nectarivore import island_incidence_streaming, normalise_binomial, crosswalk_names, build_species_table
from range_map_nectarivore_incidence import read_islands, projector, sha256

REQUIRED_FIELDS=("scientific","presence","origin","seasonal")

def validate_config(cfg):
    for k in ("provider","release","codebook_url","authorized_use_note","hashes","fields","keep","diet_mode"):
        if k not in cfg:raise ValueError("missing provider configuration: "+k)
    for k in ("provider","release","codebook_url","authorized_use_note"):
        if not isinstance(cfg[k],str) or not cfg[k].strip():raise ValueError("empty metadata")
    if cfg["diet_mode"] not in {"nectar","fruit"}:raise ValueError("only vetted diet-defined range sources")
    if set(cfg["fields"])!=set(REQUIRED_FIELDS) or len(set(cfg["fields"].values()))!=4:raise ValueError("explicit distinct field map required")
    for v in cfg["fields"].values():
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*",v):raise ValueError("unsafe field identifier")
    if set(cfg["keep"])!={"presence","origin","seasonal"}:raise ValueError("explicit code sets required")
    for codes in cfg["keep"].values():
        if not isinstance(codes,list) or not codes or any(type(c)!=int or c<0 for c in codes):raise ValueError("invalid code set")
    if not {"ranges","islands","traits"}<=set(cfg["hashes"]):raise ValueError("source pins required")
    for digest in cfg["hashes"].values():
        if not isinstance(digest,str) or not re.fullmatch("[0-9a-f]{64}",digest):raise ValueError("invalid SHA256")
    return cfg

def chunks(path,layer,cfg,chunk):
    import pyogrio
    import shapely
    info=pyogrio.read_info(path,layer=layer)
    f=cfg["fields"]
    if not set(f.values())<=set(info["fields"]):raise ValueError("source schema differs from audited codebook")
    if not info["crs"]:raise ValueError("missing CRS")
    where=" AND ".join(f"{f[k]} IN ({','.join(map(str,cfg['keep'][k]))})" for k in ("presence","origin","seasonal"))
    project=projector(info["crs"]);offset=0
    while True:
        _,_,wkb,values=pyogrio.raw.read(path,layer=layer,columns=[f["scientific"]],where=where,skip_features=offset,max_features=chunk)
        if len(wkb)==0:return
        if any(pd.isna(v) for v in values[0]):raise ValueError("unknown taxon name in retained ranges")
        geom=shapely.make_valid(project(shapely.from_wkb(wkb)))
        yield list(values[0]),geom
        offset+=len(wkb)

def traits_for_mode(path,mode):
    t=pd.read_csv(path,dtype={"scientific":str,"family":str,"clade":str})
    need={"scientific","diet_fraction","family","clade"}
    if not need<=set(t.columns) or t[list(need)].isna().any().any():raise ValueError("complete audited taxon/trait/clade rows required")
    t["key"]=t.scientific.map(normalise_binomial)
    if t.key.duplicated().any():raise ValueError("duplicate normalized taxon")
    if t.family.str.strip().eq("").any() or t.clade.str.strip().eq("").any():raise ValueError("empty taxonomy")
    frac=pd.to_numeric(t.diet_fraction,errors="raise")
    if not np.isfinite(frac).all() or not frac.between(0,100).all():raise ValueError("invalid diet fraction")
    t["diet_nect"]=np.where(frac>0,100.,0.) if mode=="fruit" else frac
    t["genus"]=t.key.str.split().str[0]
    if t.groupby("genus").clade.nunique().max()>1:raise ValueError("clade assignment inconsistent within genus")
    return t

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--config",type=pathlib.Path,required=True)
    p.add_argument("--ranges",type=pathlib.Path,required=True);p.add_argument("--layer")
    p.add_argument("--islands",type=pathlib.Path,required=True);p.add_argument("--islands-layer")
    p.add_argument("--island-id-col",required=True)
    p.add_argument("--traits",type=pathlib.Path,required=True);p.add_argument("--synonyms",type=pathlib.Path)
    p.add_argument("--out",type=pathlib.Path,required=True);p.add_argument("--chunk",type=int,default=2000)
    a=p.parse_args();cfg=validate_config(json.loads(a.config.read_text()))
    paths={"ranges":a.ranges,"islands":a.islands,"traits":a.traits}
    if a.synonyms:paths["synonyms"]=a.synonyms
    if set(paths)!=set(cfg["hashes"]):raise ValueError("each input needs an exact pin")
    for k,v in paths.items():
        if sha256(v)!=cfg["hashes"][k]:raise ValueError("source hash mismatch "+k)
    ids,geoms=read_islands(a.islands,a.islands_layer,a.island_id_col)
    import shapely
    area=shapely.area(geoms)/1e6
    inc=island_incidence_streaming(ids,geoms,chunks(a.ranges,a.layer,cfg,a.chunk),min_cover=0.,island_buffer=0.)
    t=traits_for_mode(a.traits,cfg["diet_mode"])
    syn=pd.read_csv(a.synonyms,dtype=str) if a.synonyms else None
    if syn is not None and ("source" not in syn or syn.source.fillna("").str.strip().eq("").any()):raise ValueError("each synonym needs source")
    cw=crosswalk_names(inc.species,t.scientific,syn)
    # Multiple retained range taxon names mapping to one diet taxon can double count.
    # Require explicit upstream species reconciliation instead of silently merging.
    if len(set(cw.mapping.values()))!=len(cw.mapping):raise ValueError("many-to-one range/trait taxonomy crosswalk")
    sub,sp,report=build_species_table(inc,t,cw,max_unmatched_share=.05,island_mask=area>=1)
    clades=t.set_index("key").clade
    # Primary map intersection is 50%; keep fractional coverage for prespecified sensitivity.
    a.out.mkdir(parents=True,exist_ok=True)
    sparse.save_npz(a.out/"cover_primary.npz",sub.cover)
    pd.DataFrame({"island":ids,"area_km2":area}).to_csv(a.out/"islands.csv",index=False)
    pd.DataFrame({"species":sp.species,"diet_class":sp.diet,"genus_block":sp.genus,"family_block":sp.family,
      "clade":[clades.loc[cw.mapping[s]] for s in sp.species]}).to_csv(a.out/"species_primary.csv",index=False)
    report.update(provider=cfg["provider"],release=cfg["release"],codebook_url=cfg["codebook_url"],hashes=cfg["hashes"],
      map_zero_meaning="retained range noncoverage, not absence",islands=len(ids),qualifying_area=int((area>=1).sum()),
      coverage_pairs_ge50=int((sub.cover.data>=.5).sum()),config_sha256=hashlib.sha256(a.config.read_bytes()).hexdigest())
    (a.out/"incidence_report.json").write_text(json.dumps(report,indent=2,allow_nan=False))
    # No species/island identifiers in stdout.
    print("PARTNER_RANGE_SUMMARY",json.dumps(report),flush=True)

if __name__=="__main__":main()
