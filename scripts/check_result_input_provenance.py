"""Validate DECLARED result input ancestry; not a substitute for inspecting actual runners."""
import json
from pathlib import Path

def classify(sources, ids):
    found=set()
    def visit(key, stack=()):
        if key in stack or key not in sources:
            raise ValueError(f"unknown/cyclic source: {key}")
        s=sources[key]
        kind=s["type"]
        if kind=="derived":
            parents=s.get("parents",[])
            if not parents:
                raise ValueError("derived source missing parents")
            for p in parents: visit(p, stack+(key,))
        elif kind in {"own","external","neutral","unknown"}:
            found.add(kind)
        else:
            raise ValueError(f"unknown source type: {kind}")
    for key in ids: visit(key)
    if "own" in found:
        return "island", ("mixed" if "external" in found else "own_dependent_incomplete" if "unknown" in found else "own")
    if "unknown" in found: return "hold","unverified"
    if "external" in found: return "fd","external_only"
    return "hold","no_biological_input"

def validate(registry):
    ids=set()
    for r in registry["results"]:
        if r["id"] in ids: raise ValueError("duplicate result id")
        ids.add(r["id"])
        if not r["inputs"] or not r["source_url"] or not r["audit"]:
            raise ValueError("missing provenance")
        if classify(registry["sources"],r["inputs"]) != (r["owner"],r["kind"]):
            raise ValueError(f"incorrect owner: {r['id']}")
    return len(ids)

if __name__=="__main__":
    root=Path(__file__).resolve().parents[1]
    reg=json.loads((root/"research/result_input_provenance_20261125.json").read_text())
    print(f"Validated {validate(reg)} declared result input lineages")
