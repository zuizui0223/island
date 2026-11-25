"""Inspect EuPPollNet source-original flower abundance join to visits.

Public source only, pinned 2025 Git blobs. Does not access any protected
GIFT plant/island metadata or imply ecological visit negatives or flower
availability in a given island community.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from urllib.request import Request,urlopen
from fd.euppollnet_public_preflight import SOURCE_REPO,SOURCE_COMMIT
from fd.euppollnet_flower_availability import audit_public_flower_source

def get_data(name:str,maximum:int)->bytes:
    if name not in ("Flower_counts.csv","Interaction_data.csv.gz"):
        raise ValueError("Unapproved public original flower/visitor source")
    url=("https://raw.githubusercontent.com/"+SOURCE_REPO+"/"+
         SOURCE_COMMIT+"/Data/3_Final_data/"+name)
    request=Request(url,headers={"User-Agent":"fd-flower-source-linkage-quality-1.0"})
    with urlopen(request,timeout=150) as f:
        blob=f.read(maximum+1)
    if not blob or len(blob)>maximum:
        raise ValueError("Original public source too large or unavailable")
    return blob

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    flower=get_data("Flower_counts.csv",25000000)
    interactions=get_data("Interaction_data.csv.gz",4000000)
    result=audit_public_flower_source(flower,interactions)
    body=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    for field in ("island_id","accepted_species","site_id",
                  "source_region_id","network_id","reference_id"):
        if '"'+field+'"' in body:raise ValueError("Private source identifer in public preflight")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    print("ANONYMOUS_EUPPOLLNET_FLOWER_LINK="+json.dumps({
       "n_source_flower_rows":result["n_original_flowers_raw_rows"],
       "n_source_interaction_rows":result["n_original_interaction_raw_rows"],
       "n_flower_join_keys":result["n_distinct_flower_join_keys"],
       "n_visit_join_keys":result["n_distinct_interaction_join_keys"],
       "n_shared_join_keys":result["n_original_distinct_join_keys_shared"],
       "n_shared_and_unique_flower_join_keys":
          result["n_original_join_keys_shared_and_unique_flower_count_rows"],
       "n_interaction_rows_with_unique_flower_link":
          result["n_original_interaction_rows_with_unique_flower_key"],
       "columns":result["flower_source_column_names"],
       "candidate_flower_exposure_fields":
          result["original_quantitative_flower_abundance_column_candidates"],
       "no_effect_estimate":True,
    },sort_keys=True))
if __name__=="__main__":main()
