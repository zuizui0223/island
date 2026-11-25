"""Download and audit pinned public EuPPollNet field-site metadata (no GIFT IDs).

An original source study is NOT proved independent of DoPI/GloBI/WOL.
This is a data-access/qualifying sample audit, not a biological result.
"""
import argparse,json
from pathlib import Path
from urllib.request import Request,urlopen

from fd.euppollnet_public_preflight import (
 SOURCE_REPO,SOURCE_COMMIT,SOURCE_FILES,summarize_public_euppollnet
)

def get(name:str)->bytes:
    if name not in SOURCE_FILES:
        raise ValueError("Unpinned independent literature source file")
    url=("https://raw.githubusercontent.com/"+SOURCE_REPO+"/"+
         SOURCE_COMMIT+"/Data/3_Final_data/"+name)
    req=Request(url,headers={"User-Agent":"fd-source-study-metadata-preflight/1.0"})
    with urlopen(req,timeout=90) as res:
        data=res.read(10000000)
    if len(data)>10000000:
        raise ValueError("Public source file exceeds expected budget")
    return data

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    i=get("Interaction_data.csv.gz")
    m=get("Metadata.csv")
    out=summarize_public_euppollnet(i,m)
    raw=json.dumps(out,indent=2,sort_keys=True,allow_nan=False)
    for key in ("accepted_species","island_id","original_GIFT_island",
                "source_region_id","reference_id","provider_id"):
        if '"'+key+'"' in raw:
            raise ValueError("Private original GIFT identifier in public source preflight")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    print("ANONYMOUS_EUPPOLLNET_SOURCE="+json.dumps({
        "status":out["status"],
        "n_rows":out["n_original_source_interaction_rows"],
        "n_studies":out["n_distinct_original_study_ids"],
        "n_networks":out["n_distinct_original_study_network_combinations"],
        "n_source_plant_taxa":out["n_original_source_distinct_named_plants"],
        "n_source_pollinator_taxa":out["n_original_source_distinct_named_pollinators"],
        "n_unique_documented_GPS_pairs":out["n_distinct_geographic_coordinate_pairs"],
        "n_geocoded_rows":out["n_original_georeferenced_rows"],
        "available_fields":out["original_field_contract"],
        "not_independent_from_previous_data_sources":True,
    },sort_keys=True))
if __name__=="__main__":
    main()
