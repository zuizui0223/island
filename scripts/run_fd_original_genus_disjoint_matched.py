"""SHA-locked GIFT source: within-focal matched genus-disjoint floral beta."""
import argparse,hashlib,json
from pathlib import Path
import pandas as pd
from fd_island_functional_beta import strict_traits
from fd_island_genus_disjoint_matched_beta import analyze_matched
PINS={
"traits":("fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745"),
"flora":("fixed/canonical/input/chapter1_status_flora.csv.gz","6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a")}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();paths={};pins={}
    for k,(path,expect) in PINS.items():
        f=a.original_root/path;sha=hashlib.sha256()
        with f.open("rb") as r:
            for chunk in iter(lambda:r.read(1<<20),b""):
                sha.update(chunk)
        if sha.hexdigest()!=expect:raise ValueError("GIFT source checksum drift: "+k)
        paths[k]=f;pins[k]=expect
    ledger=pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    flora=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    if flora.island_id.nunique()!=4505:raise ValueError("original GIFT N changed")
    result=analyze_matched(strict_traits(ledger),flora,
        pd.read_csv(a.geology,dtype=str),calipers=(2.,4.),boot=499)
    if result["eligible_direct_islands"]!={"all_recorded":2324,"original_native":327}:
        raise ValueError("Direct masks drifted")
    result["source_sha256"]=pins
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for label in ("island_id","accepted_species","spatial_block",
                  "island_latitude","island_longitude","genus_name",
                  "source_region_id","source_signature"):
        if '"'+label+'"' in body:raise ValueError("Original private data export blocked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    for row in result["results"]:
        if row["scope"]=="all_recorded" and row["covariate_factor"]==2.:
            print("ANONYMOUS_G0_MATCHED="+json.dumps(row,sort_keys=True))
    print(json.dumps(dict(status=result["status"],
        n_original=result["n_source_accepted_island_floras"],
        eligible_direct_islands=result["eligible_direct_islands"],
        n_qualified_panels=len(result["results"]),
        no_causal_guild_inference=True)))
if __name__=="__main__":main()
