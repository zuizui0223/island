#!/usr/bin/env python3
"""Check comprehensive review sources and export complete anonymous repertoire summaries.
No raw input restoration, biological refit, permission or publication promotion.
"""
import argparse,csv,hashlib,json,math,re,shutil
from pathlib import Path
R=Path(__file__).resolve().parent
ROOT=R.parents[2]
A=ROOT/"research/fd_own_db/current/functional_reassembly_20261128/results/20261127_source_repertoire_turnover"
M="MANUSCRIPT_BLINDED_REVIEW_COMPREHENSIVE_20261128.md"
S="SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW_COMPREHENSIVE_20261128.md"
INV=("Table_S12b_migrated_evidence_inventory.csv","Table_S12c_recent_evidence_inventory.csv")
PINS={
 "panels_from_verified_job_log.json":("88dbbb149f34e55ad5831b0faabd956648bd4418",45953256799,12),
 "baseline_partition_from_verified_job_log.json":("54c22654366f3406aae4dc21807a107ab388316b",45957671153,14),
 "directional_heterogeneity_partition.json":("5b16c954ca7a9ecd093eede58c42bd51b5c38499",45970178262,7),
}
BANNED=re.compile(r"(?i)zuizui0223|zhang ruiqi|@gmail\.com|@tohoku|10\.5281/zenodo\.22704973")
def read(p):return p.read_text(encoding="utf-8")
def blob(p):
 b=p.read_bytes()
 return hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()
def near(a,b,tol=1e-9):assert math.isclose(a,b,abs_tol=tol,rel_tol=0),(a,b)
def load(n):
 p=A/n;h,run,count=PINS[n]
 assert blob(p)==h,n
 z=json.loads(read(p))
 assert int(z["source_run"])==run and len(z["panels"])==count,n
 return z["panels"]
def rows(p):
 with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def write(path,data):
 assert data
 fields=list(dict.fromkeys(k for row in data for k in row))
 with path.open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)
 assert not BANNED.search(read(path)),path
 assert len(rows(path))==len(data)
def metric_rows(panels):
 out=[]
 for name,p in panels.items():
  gs=p.get("group_weightings",{})
  if not gs:
   out.append({"panel":name,"panel_status":p.get("status",""),"k":p.get("k","")});continue
  for g,rs in gs.items():
   for region,d in rs.items():
    base={"panel":name,"panel_status":p.get("status",""),"k":p.get("k",""),"grouping":g,"region":region,"gate":d.get("gate")}
    support=p.get("support",{}).get(region,{})
    for k in ("islands","blocks","source_groups","movable_species","movable_genera","median_target","median_donor"):
     if k in support:base[k]=support[k]
    if not d.get("metrics"):
     out.append(base);continue
    for metric,v in d["metrics"].items():
     r=dict(base,metric=metric)
     for k in ("observed","null_mean","residual","coherent_p","primary_holm_p","primary_supported"):
      if k in v:r[k]=v[k]
     for k in ("observed_ci","residual_ci","residual_leave_range"):
      if isinstance(v.get(k),list) and len(v[k])==2:r[k+"_low"],r[k+"_high"]=v[k]
     if all(isinstance(v.get(k),(int,float)) for k in ("observed","null_mean","residual")):near(v["observed"]-v["null_mean"],v["residual"])
     out.append(r)
 return out
def direction_rows(panels):
 out=[]
 for name,gs in panels.items():
  for g,rs in gs.items():
   for region,v in rs.items():
    near(v["total"],v["coherent"]+v["heterogeneity"])
    out.append(dict(panel=name,grouping=g,region=region,systematic_distance=v["total"],directional_mean=v["coherent"],heterogeneity=v["heterogeneity"]))
 return out
def check(out):
 m=read(R/M);s=read(R/S)
 old=read(R/"MANUSCRIPT_BLINDED_REVIEW_POLISHED_20261128.md")
 assert blob(R/"MANUSCRIPT_BLINDED_REVIEW_POLISHED_20261128.md")=="ba18179f87f11fba40ff2708f4121aefd7b89707"
 assert blob(R/"SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW_TRANSPORT_20261128.md")=="f61e852a9a8d37e5c74cfa63211afee877b159fd"
 assert s.startswith(read(R/"SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW_TRANSPORT_20261128.md").rstrip()+"\n")
 def result(t):return t.split("## Results\n",1)[1].split("## Discussion\n",1)[0]
 assert all(p in result(m) for p in result(old).strip().split("\n\n") if p.strip())
 a=m.split("## Abstract",1)[1].split("**Keywords:**",1)[0]
 b="## Introduction"+m.split("## Introduction",1)[1].split("## References",1)[0]
 assert len(a.split())==260 and len(b.split())==4983
 for h in ("Aim","Location","Time period","Major taxa studied","Methods","Results","Main conclusions"):
  assert a.count("**"+h+":**")==1
 assert not BANNED.search(m+s)
 assert "PENDING" in m and "176,414" in m
 assert all(f"# Appendix S{i}." in s for i in range(1,13))
 for n in ("0.06981","0.06057","10.99","0.02571","0.03669","0.02003","0.01108","0.13458","0.12907","0.01296"):
  assert n in m,n
 registry=json.loads(read(R/"FD_MIGRATION_PUBLICATION_COVERAGE_20261128.json"))
 assert blob(R/"FD_MIGRATION_PUBLICATION_COVERAGE_20261128.json")=="013cfaf081ac26b96eab3a2687c4eff262160216"
 inv,late=[rows(R/n) for n in INV]
 assert len(inv)==53 and len(late)==13
 assert {r["record_id"] for r in inv}=={r["id"] for r in registry["records"]}
 assert len({r["record_id"] for r in inv})==53
 excluded=[r for r in inv if r["disposition"]=="not used as empirical evidence"]
 assert len(excluded)==3 and all(not r["reported_result"] for r in excluded)
 assert all(not BANNED.search(read(R/n)) and n in s for n in INV)
 p,q,d=[load(n) for n in PINS]
 key="floral_symmetry__tube_depth_class__all_recorded__geo_k5__disjoint__n5"
 ck="flower_primary_color__floral_form__all_recorded__geo_k5__disjoint__n5"
 def metric(z,k,g,name):return z[k]["group_weightings"][g]["tropical_minus_north"]["metrics"][name]
 for k,expected,holm in ((key,.129073,.088),(ck,.134580,.102)):
  v=metric(p,k,"block","balanced_exclusive_cells_per_draw")
  near(v["observed"],expected,1e-6);near(v["primary_holm_p"],holm)
 v=metric(q,key,"block","systematic_occupancy_distance")
 near(v["observed"],.038840,1e-6);near(v["residual"],.012962,1e-6);near(v["primary_holm_p"],.028)
 near(d[key]["block"]["tropical_minus_north"]["total"],v["observed"])
 out.mkdir(parents=True,exist_ok=True)
 exports=(("Table_S11d_repertoire_all_panels.csv",metric_rows(p)),("Table_S11e_sampling_baseline_all_panels.csv",metric_rows(q)),("Table_S11f_directional_heterogeneity_all_panels.csv",direction_rows(d)))
 for n,rs in exports:
  assert n in s
  write(out/n,rs)
 for n in INV:shutil.copy2(R/n,out/n)
 audit={"status":"GEB_COMPREHENSIVE_SOURCE_PASS","abstract_words":len(a.split()),"body_words":len(b.split()),"original_results_preserved":True,"old_si_preserved":True,"inventory_records":len(inv),"excluded_records":len(excluded),"later_families":len(late),"panel_counts":[len(p),len(q),len(d)],"export_rows":{n:len(rs) for n,rs in exports},"source_scope":"report-level inventory plus exact three pinned recent aggregate receipts; not full verification of 73 historical files","ready_to_submit":False,"new_biological_fit":False}
 print(json.dumps(audit,allow_nan=False),flush=True)
 return audit
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path)
 check(p.parse_args().out)
