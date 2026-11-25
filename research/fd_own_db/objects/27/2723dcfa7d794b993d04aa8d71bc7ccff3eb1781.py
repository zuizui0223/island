"""Multi-resolution TRUE species-Direct flower phenotype battery for DoPI visitors.

One-axis, genuine same accepted-plant two-axis, and three/four-axis motifs
are separate sensitivity panels. No global "triple complete" gate and no
product of marginals observed on different accepted species. Every recorded
motif is nonexclusive, NOT a named pollination syndrome.
"""
from __future__ import annotations
from collections import Counter
from .real_joint_syndrome_access import JOINT_SPECS, make_joint_species
from .syndrome_ensemble_source import MOTIFS, make_motifs

SINGLE_PROFILES={
 "deep_tube":("tube_depth_class",frozenset({"deep"})),
 "tubular_bell":("floral_form",frozenset({"tubular","bell_campanulate"})),
 "zygomorphic":("floral_symmetry",frozenset({"zygomorphic"})),
 "open_radial":("floral_form",frozenset({"open_radial"})),
}


def make_real_multiresolution_motifs(states:dict[str,dict[str,str]])->dict:
    """27 true Direct floral phenotypes; 4 single + 11 joint + 12 triple/quad."""
    multi=make_motifs(states)
    pairs=make_joint_species(states)
    if len(multi)!=len(MOTIFS) or len(pairs)!=len(JOINT_SPECS):
        raise ValueError("original registered real profiles missing")
    out={}
    out.update(multi)
    out.update(pairs)
    for name,(axis,values) in SINGLE_PROFILES.items():
        if axis not in states:
            raise ValueError("original Direct single flower category absent")
        known=frozenset(states[axis])
        true=frozenset(s for s in known if states[axis][s] in values)
        out[name]=dict(known=known,joint=true,n_axes=1,axes=(axis,))
    if len(out)!=len(MOTIFS)+len(JOINT_SPECS)+len(SINGLE_PROFILES):
        raise AssertionError("multi-resolution name collision")
    for name,m in out.items():
        if not set(m["joint"])<=set(m["known"]):
            raise ValueError(f"named phenotype positive outside documented direct species {name}")
        if not m["known"]:
            raise ValueError(f"original GIFT named Direct motif without source observations: {name}")
        if name in multi:
            size=m["n_axes"]
        elif name in pairs:
            size=len(m["trait_axes"])
        else:
            size=1
        if size<1 or size>4:
            raise ValueError("invalid nonexclusive phenotype axis count")
        m["axis_count"]=size
    return out


def botanical_axis_coverage(motifs:dict)->dict:
    return dict(sorted(Counter(x["axis_count"] for x in motifs.values()).items()))
