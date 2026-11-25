"""Frozen plant syndrome proxies. No animal-guild assignment or inferred missing traits."""
from itertools import combinations
import numpy as np

CORE=("flower_primary_color","floral_form","floral_symmetry","tube_depth_class")
SIZE="flower_size_class"
VALID={
 CORE[0]:set("white green_brown_inconspicuous yellow_orange red_pink blue_purple multicolored_variable other_described".split()),
 CORE[1]:set("open_radial bell_campanulate tubular salverform funnel_trumpet urn_urceolate brush_puff composite_head papilionaceous bilabiate spurred reduced_wind other_described".split()),
 CORE[2]:set("actinomorphic zygomorphic asymmetric".split()),
 CORE[3]:set("absent_or_open shallow intermediate deep".split()),
 SIZE:set("very_small small medium large very_large".split())}
TEMPLATES={
 "specialist_bee":{CORE[0]:"blue_purple yellow_orange",CORE[1]:"bilabiate papilionaceous",CORE[2]:"zygomorphic",CORE[3]:"intermediate",SIZE:"medium"},
 "bird":{CORE[0]:"red_pink yellow_orange",CORE[1]:"tubular funnel_trumpet",CORE[3]:"deep",SIZE:"large very_large"},
 "hawkmoth":{CORE[0]:"white",CORE[1]:"salverform tubular",CORE[3]:"deep"},
 "bat":{CORE[0]:"white green_brown_inconspicuous",CORE[1]:"bell_campanulate brush_puff",SIZE:"large very_large"},
 "open_generalist":{CORE[0]:"white yellow_orange",CORE[1]:"open_radial composite_head",CORE[2]:"actinomorphic",CORE[3]:"absent_or_open shallow",SIZE:"very_small small"}}

def encode(frame,axes=CORE,red_only=False,multistate=False):
    """Return compatibility indicators, fixed pair weights, and per-species scores.
    A state-set is accepted ONLY in the separately labelled public pilot.
    """
    xs=[];pairs=[];weights=[];labels=[];scores=[];strict=[]
    for k,spec in TEMPLATES.items():
        positions=[]
        for axis in axes:
            if axis not in spec:continue
            allowed=set(spec[axis].split())
            if red_only and k=="bird" and axis==CORE[0]:allowed={"red_pink"}
            v=frame[axis].astype(str)
            if multistate:
                vals=v.str.split("|",regex=False)
                assert vals.map(lambda a:bool(a) and set(a)<=VALID[axis]).all()
                z=vals.map(lambda a:bool(set(a)&allowed)).to_numpy(float)
            else:
                assert v.isin(VALID[axis]).all()
                z=v.isin(allowed).to_numpy(float)
            positions.append(len(xs));xs.append(z)
        assert len(positions)>=2
        pp=list(combinations(positions,2));pairs.extend(pp)
        weights.extend([1/(len(TEMPLATES)*len(pp))]*len(pp));labels.extend([k]*len(pp))
        mat=np.column_stack([xs[j] for j in positions])
        scores.append(mat.mean(1));strict.append(mat.min(1))
    return np.column_stack(xs),np.asarray(pairs,int),np.asarray(weights),np.asarray(labels),np.column_stack(scores),np.column_stack(strict)

def components(X,pairs,group=None,counts=None):
    """Signed excess joint compatibility; finite-n corrected exact decomposition."""
    counts=np.ones(len(X)) if counts is None else np.asarray(counts,float)
    assert np.isfinite(counts).all() and (counts>0).all() and np.all(counts==np.floor(counts))
    n=float(counts.sum());assert n>1
    w=counts/n
    mean=w@X;xy=w@(X[:,pairs[:,0]]*X[:,pairs[:,1]])
    raw=xy-mean[pairs[:,0]]*mean[pairs[:,1]]
    scale=n/(n-1)
    var=mean*(1-mean)
    denom=np.sqrt(var[pairs[:,0]]*var[pairs[:,1]])
    phi=np.divide(raw,denom,out=np.full(len(raw),np.nan),where=denom>1e-12)
    out={"total":scale*raw,"phi":phi}
    if group is not None:
        _,g=np.unique(group,return_inverse=True)
        mass=np.bincount(g,weights=counts);weights_g=mass/n
        mu=np.column_stack([np.bincount(g,weights=counts*X[:,j])/mass for j in range(X.shape[1])])
        between=weights_g@(mu[:,pairs[:,0]]*mu[:,pairs[:,1]])-mean[pairs[:,0]]*mean[pairs[:,1]]
        out["between"]=scale*between;out["within"]=scale*(raw-between)
        assert np.allclose(out["within"]+out["between"],out["total"],atol=1e-12)
    return out

def common_phi(a,b,labels):
    # Keep the same pairs within each template and equal template weights.
    ok=np.isfinite(a)&np.isfinite(b)
    v=[]
    for k in TEMPLATES:
        mask=ok&(labels==k)
        if not mask.any():return None
        v.append(float((a[mask]-b[mask]).mean()))
    return float(np.mean(v))

def self_test():
    # Equal margins, opposite pairing, and exact within/between identities.
    p=np.array([[0,1]]);x=np.array([[0,0],[0,0],[1,1],[1,1]],float)
    y=np.array([[0,0],[0,1],[1,0],[1,1]],float)
    g=np.array([0,0,1,1])
    a=components(x,p,g);b=components(y,p,g)
    assert np.isclose(a["total"][0],1/3) and np.isclose(b["total"][0],0)
    assert np.isclose(a["within"][0],0) and np.isclose(a["between"][0],1/3)
    assert np.isclose(components(x,p,counts=[2,2,2,2])["total"][0],2/7)
    assert np.isclose(a["phi"][0],1)
    print("PLANT_SYNDROME_CORE_TEST_PASS",flush=True)
if __name__=="__main__":self_test()
