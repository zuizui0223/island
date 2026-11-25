"""Exact algebra kernel for the FD depth-conditional architecture contract.

No I/O, identifiers, external access or protected workflow. Input is one row
per distinct scored species, tuple (symmetry, depth, genus, exact grade).
The custodian must enforce source pins, roster identity, native scope,
eligibility and geographic/taxon resampling. Fractions preserve identities.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import product

SYMMETRY=("actinomorphic","zygomorphic","asymmetric")
DEPTH=("absent_or_open","shallow","intermediate","deep")
CELLS=tuple(product(SYMMETRY,DEPTH))
VECTORS=("total","common_within_stratum","common_weight_change","target_only","source_only")

def freq(counts,keys):
    n=sum(counts.values())
    if not n:raise ValueError("empty carrier set")
    return tuple(F(counts[k],n) for k in keys)

def joint_parts(rows):
    """Population-frequency identity, with no finite-n scaling."""
    rows=list(rows)
    n=len(rows)
    if not n:raise ValueError("empty roster")
    c=Counter((s,d) for s,d,g,q in rows)
    if set(c)-set(CELLS):raise ValueError("unknown ontology state")
    a=Counter(s for s,d,g,q in rows);b=Counter(d for s,d,g,q in rows)
    groups=defaultdict(list)
    for s,d,g,q in rows:
        if not g or not q:raise ValueError("missing stratum")
        groups[(g,q)].append((s,d))
    observed=freq(c,CELLS)
    margin=tuple(F(a[s]*b[d],n*n) for s,d in CELLS)
    expected=[F(0)]*len(CELLS)
    for group in groups.values():
        h=len(group);ac=Counter(s for s,d in group);bc=Counter(d for s,d in group)
        for j,(s,d) in enumerate(CELLS):
            expected[j]+=F(ac[s]*bc[d],n*h)
    between=tuple(x-y for x,y in zip(expected,margin))
    within=tuple(x-y for x,y in zip(observed,expected))
    total=tuple(x-y for x,y in zip(observed,margin))
    assert all(t==x+y for t,x,y in zip(total,between,within))
    assert all(sum(x)==0 for x in (between,within,total))
    return {"n":n,"joint":observed,"symmetry_marginal":freq(a,SYMMETRY),
            "depth_marginal":freq(b,DEPTH),"independent_joint":margin,
            "genus_grade_expected_joint":tuple(expected),
            "covariance":total,"between":between,"within":within}

def subset(rows,depth):
    if depth=="nondeep":return [r for r in rows if r[1]!="deep"]
    if depth not in DEPTH:raise ValueError("unknown depth")
    return [r for r in rows if r[1]==depth]

def conditional_parts(target,source,depth):
    """Exact symmetric decomposition of conditional symmetry frequencies.

    One-sided stratum contributions retain their sign; common composition
    and common within-stratum terms are symmetric, not causal effects.
    """
    def summarize(rows):
        by=defaultdict(Counter)
        for s,d,g,q in subset(rows,depth):
            by[(g,q)][s]+=1
        return by,sum(sum(c.values()) for c in by.values())
    t,nt=summarize(target);u,nu=summarize(source)
    if not nt or not nu:return {"status":"undefined","target_n":nt,"source_n":nu}
    common=set(t)&set(u);only_t=set(t)-set(u);only_u=set(u)-set(t)
    out={k:[F(0)]*3 for k in VECTORS}
    tq=[F(0)]*3;sq=[F(0)]*3
    for k in sorted(set(t)|set(u)):
        a=t.get(k,Counter());b=u.get(k,Counter())
        n=sum(a.values());m=sum(b.values())
        w=F(n,nt);v=F(m,nu)
        qa=freq(a,SYMMETRY) if n else (F(0),)*3
        qb=freq(b,SYMMETRY) if m else (F(0),)*3
        for i,(x,y) in enumerate(zip(qa,qb)):
            tq[i]+=w*x;sq[i]+=v*y
            if k in common:
                out["common_within_stratum"][i]+=(w+v)*(x-y)/2
                out["common_weight_change"][i]+=(w-v)*(x+y)/2
            elif k in only_t:out["target_only"][i]+=w*x
            else:out["source_only"][i]-=v*y
    out["total"]=[x-y for x,y in zip(tq,sq)]
    assert sum(tq)==sum(sq)==1
    for i in range(3):assert out["total"][i]==sum(out[k][i] for k in VECTORS[1:])
    assert sum(out["total"])==sum(out["common_within_stratum"])==0
    return {"status":"defined","target_n":nt,"source_n":nu,
            "target_conditional":tuple(tq),"source_conditional":tuple(sq),
            **{k:tuple(v) for k,v in out.items()},
            "common_strata":len(common),"target_only_strata":len(only_t),
            "source_only_strata":len(only_u),
            "common_target_mass":F(sum(sum(t[k].values()) for k in common),nt),
            "common_source_mass":F(sum(sum(u[k].values()) for k in common),nu),
            "common_variable_strata":sum(len(set(t[k])|set(u[k]))>1 for k in common),
            "common_changed_distribution_strata":sum(freq(t[k],SYMMETRY)!=freq(u[k],SYMMETRY) for k in common)}

def pair(target,source,min_each=10):
    target=list(target);source=list(source)
    t=joint_parts(target);s=joint_parts(source)
    c={d:conditional_parts(target,source,d) for d in (*DEPTH,"nondeep")}
    deep,other=c["deep"],c["nondeep"]
    eligible=all(v["status"]=="defined" and min(v["target_n"],v["source_n"])>=min_each for v in (deep,other))
    delta={k:tuple(x-y for x,y in zip(t[k],s[k])) for k in
           ("joint","symmetry_marginal","depth_marginal","independent_joint","covariance","between","within")}
    assert all(j==m+x+y for j,m,x,y in zip(delta["joint"],delta["independent_joint"],delta["between"],delta["within"]))
    contrast=tuple(x-y for x,y in zip(deep["total"],other["total"])) if eligible else None
    return {"eligible":eligible,"delta":delta,"conditional":c,
            "deep_vs_nondeep_difference_of_differences":contrast}

def self_test():
    from random import Random
    a,b,z=SYMMETRY;d=DEPTH[-1];e=DEPTH[0]
    t=[(a,d,"A","high|high"),(b,d,"B","medium|high"),(a,e,"A","high|high")]
    s=[(b,d,"A","high|high"),(a,d,"C","high|medium"),(b,e,"A","high|high")]
    c=conditional_parts(t,s,d)
    assert c["total"]==(F(0),F(0),F(0))
    assert c["common_within_stratum"]==(F(1,2),F(-1,2),F(0))
    assert c["target_only"]==(F(0),F(1,2),F(0))
    assert c["source_only"]==(F(-1,2),F(0),F(0))
    assert not pair(t,s)["eligible"] and pair(t,s,min_each=1)["eligible"]
    assert pair(t,s,min_each=1)["deep_vs_nondeep_difference_of_differences"]==(F(-1),F(1),F(0))
    rng=Random(20261125)
    for _ in range(100):
        def make(n):return [(rng.choice(SYMMETRY),rng.choice(DEPTH),str(rng.randrange(5)),str(rng.randrange(2))) for _ in range(n)]
        t=make(45);s=make(60)
        x=pair(t,s,min_each=1);y=pair(s,t,min_each=1)
        for k in x["delta"]:assert x["delta"][k]==tuple(-v for v in y["delta"][k])
        for d in (*DEPTH,"nondeep"):
            p=x["conditional"][d];q=y["conditional"][d]
            if p["status"]=="defined":
                assert p["total"]==tuple(-v for v in q["total"])
                assert p["target_only"]==tuple(-v for v in q["source_only"])
    print("FD_DEPTH_CONDITIONAL_SYNTHETIC_PASS: exact identities, 100 random pairs, no empirical input")

if __name__=="__main__":self_test()
