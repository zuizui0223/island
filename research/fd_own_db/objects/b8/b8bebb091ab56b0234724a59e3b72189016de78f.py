"""Patristic distances on a dated ORIGINAL phylogenetic backbone tree.

This helper never makes taxonomy-derived branches or imputes missing tips.
Uses source original phylo edge lists with root-to-tip distances and binary
lifting for LCA. A compact distance matrix is computed only for source-
verified exact original Direct species tips, held in private CI memory.

Branch lengths are those published by pinned GBOTB.extended.*.rda, not
re-estimated from user floral traits or pollen interactions.
"""
from __future__ import annotations
import numpy as np


def true_tip_patristic_distance(
    edges:np.ndarray, edge_lengths:np.ndarray,
    tip_nodes:np.ndarray, *,
    n_original_tree_tips:int
)->np.ndarray:
    edge=np.asarray(edges,dtype=np.int64)
    w=np.asarray(edge_lengths,dtype=np.float64)
    tips=np.asarray(tip_nodes,dtype=np.int64)
    if edge.ndim!=2 or edge.shape[1]!=2 or len(edge)!=len(w):
        raise ValueError("Original phylo source edge or lengths invalid")
    if len(edge)<4 or n_original_tree_tips<2:
        raise ValueError("Insufficient dated original backbone")
    if np.isnan(w).any() or not np.isfinite(w).all() or (w<=0).any():
        raise ValueError("Nonpositive or nonfinite original tree branch length")
    n_nodes=int(edge.max())+1
    if int(edge.min())<1 or n_nodes>400000 or len(edge)!=n_nodes-2:
        raise ValueError("Original rooted tree is not a single node-indexed tree")
    if np.unique(tips).size!=len(tips) or len(tips)<2:
        raise ValueError("Original exact-backed tips not distinct")
    if tips.min()<1 or tips.max()>n_original_tree_tips:
        raise ValueError("Requested taxon is not a measured original tree TIP")
    parents=np.zeros(n_nodes,dtype=np.int32)
    parents[edge[:,1]]=edge[:,0].astype(np.int32)
    if np.any(np.bincount(edge[:,1],minlength=n_nodes)[1:]>1):
        raise ValueError("Source is not a valid rooted tree: multiple parents")
    roots=np.setdiff1d(np.unique(edge[:,0]),np.unique(edge[:,1]))
    if len(roots)!=1:
        raise ValueError("Root not unique in original dated source tree")
    root=int(roots[0])
    child_groups=[[] for _ in range(n_nodes)]
    for p,c in edge:
        child_groups[int(p)].append(int(c))
    weights=np.zeros(n_nodes,dtype=float)
    weights[edge[:,1]]=w
    height=np.zeros(n_nodes,dtype=np.int32)
    rootdist=np.zeros(n_nodes,dtype=np.float64)
    queue=[root]
    ix=0
    while ix<len(queue):
        node=queue[ix]
        ix+=1
        for child in child_groups[node]:
            height[child]=height[node]+1
            rootdist[child]=rootdist[node]+weights[child]
            queue.append(child)
    if len(queue)!=n_nodes-1:
        raise ValueError("Original source tree disconnected, unrooted or cyclic")
    K=int(np.ceil(np.log2(max(1,height.max()))))+1
    ancestors=np.empty((K,n_nodes),dtype=np.int32)
    ancestors[0]=parents
    for bit in range(1,K):
        ancestors[bit]=ancestors[bit-1,ancestors[bit-1]]
    D=np.zeros((len(tips),len(tips)),dtype=np.float64)
    for i in range(len(tips)-1):
        others=tips[i+1:].copy().astype(np.int32)
        a=np.full(len(others),int(tips[i]),dtype=np.int32)
        b=others.copy()
        deeper=height[a]<height[b]
        swap=a[deeper].copy()
        a[deeper]=b[deeper]
        b[deeper]=swap
        diff=height[a]-height[b]
        for k in range(K):
            which=((diff>>k)&1).astype(bool)
            a[which]=ancestors[k,a[which]]
        neq=a!=b
        for k in range(K-1,-1,-1):
            different=(ancestors[k,a]!=ancestors[k,b]) & neq
            a[different]=ancestors[k,a[different]]
            b[different]=ancestors[k,b[different]]
        lca=np.where(neq,parents[a],a)
        if np.any(lca<=0):
            raise ValueError("Source tip LCA missing")
        vals=rootdist[tips[i]]+rootdist[others]-2*rootdist[lca]
        if np.any(vals<(-1e-5)):
            raise ValueError("Published branch lengths not a valid positive tree metric")
        vals=np.maximum(0.,vals)
        D[i,i+1:]=vals
        D[i+1:,i]=vals
    if not np.allclose(D,D.T) or not np.allclose(np.diag(D),0):
        raise ValueError("Source patristic metric is not symmetric")
    return D
