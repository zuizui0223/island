"""Count how widely GBIF-recorded taxa occur among sampled islands.

This is a *sampling-distribution audit*, not true island endemism or a
functional-trait finding. The only output is an aggregate JSON file.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def summarize(path: Path, status: Path | None = None) -> dict:
    category = {'1':0,'2':0,'3_5':0,'6_10':0,'11_plus':0}
    names: set[str] = set()
    counts: list[int] = []
    total_pairs = 0
    with path.open('r',encoding='utf-8',newline='') as handle:
        reader = csv.DictReader(handle)
        columns = {'accepted_species','genus','n_islands','n_records'}
        if not columns.issubset(reader.fieldnames or []):
            raise ValueError(f'Missing required source columns: {columns-set(reader.fieldnames or [])}')
        for i,row in enumerate(reader, start=2):
            name = row['accepted_species'].strip()
            if not name or name in names:
                raise ValueError(f'Blank/duplicated accepted_species at CSV record {i}')
            names.add(name)
            try:
                islands, records = int(row['n_islands']), int(row['n_records'])
            except (ValueError,TypeError) as e:
                raise ValueError(f'Invalid counts at record {i}') from e
            if islands<1 or records<islands:
                raise ValueError(f'Impossible count at record {i}')
            counts.append(islands)
            total_pairs+=islands
            key = '1' if islands==1 else '2' if islands==2 else '3_5' if islands<=5 else '6_10' if islands<=10 else '11_plus'
            category[key]+=1
    if not counts:
        raise ValueError('No records')
    counts.sort()
    n=len(counts)
    result={
       'scientific_scope':'GBIF recorded island frequencies, not true endemism or functional traits',
       'n_accepted_taxon_names':n,
       'n_island_taxon_occurrences':total_pairs,
       'taxa_recorded_on_one_island':category['1'],
       'fraction_recorded_on_one_island':category['1']/n,
       'island_count_by_taxon_recording_frequency':category,
       'median_recorded_islands_per_taxon':counts[(n-1)//2],
       'p90_recorded_islands_per_taxon':counts[int(.9*(n-1))],
       'p99_recorded_islands_per_taxon':counts[int(.99*(n-1))],
       'max_recorded_islands_per_taxon':counts[-1],
    }
    if status is not None:
        state=json.loads(status.read_text(encoding='utf-8'))
        if n!=state['n_accepted_species'] or total_pairs!=state['n_island_species_pairs']:
            raise ValueError('Taxon file and GBIF campaign status disagree')
        result['n_islands_with_any_record']=int(state['n_islands_with_records'])
        result['n_source_blocks']=int(state['n_collected_succeeded_blocks'])
    return result


def main(argv:list[str]|None=None)->None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--taxa',type=Path,required=True)
    parser.add_argument('--status',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    summary=summarize(args.taxa,args.status)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
