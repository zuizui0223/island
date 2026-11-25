#!/usr/bin/env python3
"""Check review-candidate numerical tables against saved original-run receipts."""
import csv,json
from pathlib import Path
from hashlib import sha1
from math import isclose
R=Path(__file__).resolve().parent
M=json.loads((R/"REVIEW_PACKAGE_MANIFEST.json").read_text())
def blob(p):
    b=p.read_bytes()
    return sha1(f"blob {len(b)}\0".encode()+b).hexdigest()
for x in M["receipts"]:
    assert blob(R/x["path"])==x["git_blob_sha"],x["path"]
for x in M["tables"]:
    with (R/x["path"]).open(newline="") as f:r=list(csv.DictReader(f))
    assert len(r)==x["rows"],x["path"]
h1=json.loads((R/"FLORAL_AXIS_DECOUPLING_RESULTS.json").read_text())
h2=json.loads((R/"FLORAL_AXIS_H2_RESULT.json").read_text())
w=json.loads((R/"FLORAL_AXIS_H2_WEIGHTS_RESULT.json").read_text())
c=json.loads((R/"FLORAL_AXIS_COMPONENT_RESULT.json").read_text())
h4=json.loads((R/"FLORAL_AXIS_H4_RESULT.json").read_text())
assert len(h1["panels"])==60 and len(h2["panels"])==len(w["panels"])==24
assert sum(x.get("positive_parent_and_joint_sensitivity",False) for x in w["panels"])==4
assert all(not x["support"]["pass"] and x["status"]=="insufficient_support" for x in h4["panels"])
for p in c["panels"]:
    x=p["metrics"]
    assert isclose(x["radial_given_deep"]["beta"]-x["radial_given_shallow_open"]["beta"],x["deep_minus_shallow"]["beta"],abs_tol=1e-10)
def rows(name):
    with (R/"supplement"/name).open(newline="") as f:return list(csv.DictReader(f))
for a,b in zip(rows("Table_S8a_flower_pair_isolation.csv"),h1["panels"],strict=True):
    assert a["scope"]==b["scope"] and a["endpoint"]==b["endpoint"] and a["region"]==b["region"]
    assert int(a["islands"])==b["n"] and a["status"]==b["status"]
    if b["status"]=="fit":assert isclose(float(a["beta"]),b["beta"],abs_tol=1e-14)
for a,b in zip(rows("Table_S8b_reproductive_conditioning.csv"),h2["panels"],strict=True):
    assert (a["scope"],a["roster"],a["score"])==(b["scope"],b["roster"],b["score"])
    if b["status"]=="fit":assert isclose(float(a["adjusted_beta"]),b["adjusted_same_islands"]["beta"],abs_tol=1e-14)
for a,(p,k) in zip(rows("Table_S8c_depth_specific_components.csv"),[(p,k) for p in c["panels"] for k in c["metric_order"]],strict=True):
    assert a["metric"]==k and isclose(float(a["beta"]),p["metrics"][k]["beta"],abs_tol=1e-14)
for a,b in zip(rows("Table_S8e_GloPL_joint_support.csv"),[c for p in h4["panels"] for c in p["support"]["four_cells"]],strict=True):
    assert int(a["cell_species"])==b["species"]
assert all((R/p).exists() for p in ["MANUSCRIPT_REVIEW.md","SUPPLEMENTARY_INFORMATION_REVIEW.md","FLORAL_AXIS_COMPONENT_FIGURE.svg","PRIOR_ART_FLORAL_AXIS.md"])
print("REVIEW_PACKAGE_PASS: 6 source receipts; 5 tables; 117 table rows; exact component identities; no H4 fit")
