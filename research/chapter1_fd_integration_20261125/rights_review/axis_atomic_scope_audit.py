#!/usr/bin/env python3
"""Exact frozen axis-to-trait provenance audit. Aggregate only; no rights promotion.

The historical Direct snapshot expands a species-axis cell's quality and lineage
set to all its component traits. This audit distinguishes that original scope
from the selected same-species same-state atomic Direct ledger. It does not
assert that absent atomic rows are biologically wrong, or fit new models.
"""
from __future__ import annotations
import argparse, hashlib, json, zipfile
from pathlib import Path
import pandas as pd
PINS={
 "analysis":("10058653212-chapter1-progressive-analysis-34232450884.zip","b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999"),
 "ledger":("10042377593-source-scale-batch-integration-34191508045.zip","5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd"),
 "rights":("10186293306-chapter1-database-v1-cell-rights-34566271397.zip","6c5f2577e85e95e841b8659427e27c72ce631f59a3da5a335f54f505de5ea32a"),
}
CORE={"self_incompatibility","mating_system","autonomous_selfing_capacity"}
STATES={
 "self_incompatibility":{"SC","SI"},
 "mating_system":{"predominantly_outcrossing","predominantly_selfing","obligate_selfing"},
 "autonomous_selfing_capacity":{"autonomous","delayed","absent"},
 "flower_primary_color":{"white","green_brown_inconspicuous","yellow_orange","red_pink","blue_purple"},
 "floral_form":{"open_radial","brush_puff","composite_head","tubular","bell_campanulate","funnel_trumpet","urn_urceolate","bilabiate","salverform","papilionaceous","spurred"},
 "floral_symmetry":{"actinomorphic","zygomorphic"},
 "tube_depth_class":{"absent_or_open","shallow","deep"},
}
SNAP_AXIS="snapshot/chapter1_trait_snapshot_species_axis.csv.gz"
SNAP_TRAIT="snapshot/chapter1_trait_ledger_direct_only.csv.gz"
SNAP_STATE="snapshot/chapter1_trait_state_audit_direct.csv.gz"
PINNED_INNER={
 "snapshot_axis":"a443cda3025f0cadd2ccb58a0529bfbe30459f445cf8909093b6b88124d51aea",
 "snapshot_trait":"d7d6a0372ce7f3ebe958047a222ae0504f01564e364f93bab360273235c6661b",
 "snapshot_state":"4890da7df53b2e240c4b1493ef93c0814a68435a5e1670e3c0e7c3e6717b90e7",
 "v1_axis":"a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a",
 "v1_atomic":"7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577",
}
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def read(z,n,pin=None):
 if pin is not None:assert hashlib.sha256(z.read(n)).hexdigest()==PINNED_INNER[pin],pin
 return pd.read_csv(z.open(n),compression="gzip",dtype=str,keep_default_na=False)
def counts(x):return {str(k):int(v) for k,v in x.value_counts().items()}
def tokens(x):return {s for s in x.split("|") if s}
def run(folder):
 z={}
 for k,(n,h) in PINS.items():
  p=folder/n;assert sha(p)==h,k;z[k]=zipfile.ZipFile(p)
 key=["accepted_species","axis"];tkey=["accepted_species","trait_name"]
 v=read(z["ledger"],"species_axis_coverage.csv.gz","v1_axis").sort_values(key).reset_index(drop=True)
 s=read(z["analysis"],SNAP_AXIS,"snapshot_axis").sort_values(key).reset_index(drop=True)
 assert len(v)==318885 and not v.duplicated(key).any() and s.equals(v)
 m=json.loads(z["analysis"].read("snapshot/chapter1_trait_snapshot_manifest.json"))
 assert m["source_species_axis_sha256"]==PINNED_INNER["v1_axis"]
 d=read(z["ledger"],"direct_species_trait_ledger.csv.gz","v1_atomic")
 assert hashlib.sha256(z["analysis"].read("trait-source/direct_species_trait_ledger.csv.gz")).hexdigest()==PINNED_INNER["v1_atomic"]
 assert not d.duplicated(tkey).any()
 a=read(z["analysis"],SNAP_STATE,"snapshot_state")
 t=read(z["analysis"],SNAP_TRAIT,"snapshot_trait")
 assert not a.duplicated(tkey).any() and not t.duplicated(tkey).any()
 a=a[a.trait_name.isin(STATES)&a.resolved_for_primary.eq("True")]
 a=a[a.apply(lambda x:x.canonical_signature in STATES[x.trait_name],axis=1)]
 x=a.merge(t,on=tkey,validate="one_to_one")
 assert len(x)==len(a) and x.canonical_signature.eq(x.normalized_value).all()
 x=x.merge(v[key+["source_lineages","quality"]].rename(columns={"source_lineages":"axis_lineages","quality":"axis_quality"}),on=key,validate="many_to_one")
 assert x.source_lineages.eq(x.axis_lineages).all() and x.quality.eq(x.axis_quality).all()
 r=read(z["rights"],"CELL_RELEASE_RIGHTS.csv.gz")
 assert not r.duplicated(key).any()
 x=x.merge(r[key+["release_status","blocked_source_families"]],on=key,validate="many_to_one")
 assert x.release_status.isin(["redistributable","rights_review_required"]).all()
 cols=["state_set","quality","resolution_status","source_lineages"]
 rename={"quality":"atomic_quality","source_lineages":"atomic_lineages"}
 x=x.merge(d[tkey+cols].rename(columns=rename),on=tkey,how="left",indicator=True,validate="one_to_one")
 acols=["state_set","atomic_quality","resolution_status","atomic_lineages"]
 x[acols]=x[acols].fillna("")
 x["atomic_state"]=x.state_set.map(lambda y:"|".join(sorted(json.loads(y))) if y else "")
 x["exact_atomic"]=x.atomic_state.eq(x.canonical_signature)&x.resolution_status.eq("resolved")&x.atomic_quality.isin(["high","medium"])
 x["atomic_category"]="state_or_quality_mismatch"
 x.loc[x._merge.eq("left_only"),"atomic_category"]="absent_atomic_row"
 x.loc[x.exact_atomic,"atomic_category"]="exact_high_medium_state"
 tables={}
 for name,g in x.groupby("trait_name"):
  tables[name]={"axis_direct_rows":len(g),"atomic_categories":counts(g.atomic_category),"axis_rights_status":counts(g.release_status)}
 c=x[x.trait_name.isin(CORE)]
 exact=c[c.exact_atomic]
 assert len(c)==12401 and len(exact)==10654
 assert (c.atomic_category=="absent_atomic_row").sum()==1747
 subset=sum(tokens(y.atomic_lineages)<=tokens(y.source_lineages) for y in exact.itertuples())
 equal=sum(tokens(y.atomic_lineages)==tokens(y.source_lineages) for y in exact.itertuples())
 assert subset==10654 and equal==4527
 cc=c.drop_duplicates(key)
 targets=["database:floraweb_biolflor","database:pladias","domain:europepmc.org","ken-fern-useful-plants:species","domain:gobotany.nativeplanttrust.org"]
 families={}
 for f in targets:
  g=c[c.blocked_source_families.map(lambda s:f in tokens(s))]
  sole=g[g.blocked_source_families.eq(f)]
  families[f]={"matched_core_axis_rows_affected":len(g),"sole_family_axis_rows":len(sole),"sole_family_rows_by_trait":counts(sole.trait_name)}
 complete=lambda a,n:int(a.groupby("accepted_species").size().eq(n).sum())
 pair=x[x.trait_name.isin(["floral_symmetry","tube_depth_class"])]
 return {"source_archive_sha256":{k:h for k,(_,h) in PINS.items()},"inner_pins":PINNED_INNER,
  "all_species_axis_rows_identical":len(v),"all_axis_columns_identical":True,"axis_metadata_inherited_for_every_informative_H1_row":True,
  "trait_tables":tables,"core":{"records":len(c),"unique_species_axis_cells":len(cc),"axis_rights_status":counts(c.release_status),"unique_axis_rights_status":counts(cc.release_status),
  "exact_atomic":len(exact),"absent_atomic":int(c.atomic_category.eq("absent_atomic_row").sum()),"atomic_state_or_quality_mismatches":int(c.atomic_category.eq("state_or_quality_mismatch").sum()),"atomic_lineage_sets_subset_of_axis":subset,"atomic_lineage_sets_equal_axis":equal,"axis_blocker_families":families},
  "complete_three_binary_roster":{"axis_direct":complete(c,3),"exact_atomic":complete(exact,3)},
  "informative_single_state_binary_symmetry_depth_roster":{"axis_direct":complete(pair,2),"exact_atomic":complete(pair[pair.exact_atomic],2)},
  "interpretation":"Axis snapshot identity supports exact axis-cell rights crosswalk for all historical H1 informative rows. Axis source_lineages are inherited by component traits, not trait-specific provenance. Exact atomic matches are a stricter data-source roster, not a license grant. Absent atomic rows are not proven erroneous. Binary symmetry-depth roster here uses final H1 binary states, not the full paired-trait analysis cohort. No new model fit, rights promotion, frozen manuscript or selector changes."}
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--source-dir",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
 a=p.parse_args();res=run(a.source_dir)
 a.out.write_text(json.dumps(res,indent=2,sort_keys=True,allow_nan=False)+"\n")
 print("AXIS_ATOMIC_SCOPE_RESULT "+json.dumps(res,sort_keys=True,allow_nan=False),flush=True)
