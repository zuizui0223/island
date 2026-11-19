#!/usr/bin/env Rscript
# Export original published GBOTB topology/branch lengths only to PRIVATE CI.
# Source plant acceptance mapping is confidential; no names in public logs.
args=commandArgs(trailingOnly=TRUE)
if(length(args)!=5)stop("Rscript exporter backbone_data source_species.tsv tree_version edge.tsv tip.tsv")
root=args[[1]]
sourcefile=args[[2]]
version=args[[3]]
outEdge=args[[4]]
outTips=args[[5]]
if(!version %in% c("TPL","LCVP","WP"))stop("Unknown tree source")
env=new.env(parent=baseenv())
nm=load(file.path(root,paste0("GBOTB.extended.",version,".rda")),envir=env)
tree=env[[nm[[1]]]]
if(anyNA(tree$edge.length)||any(tree$edge.length<=0))stop("Original branch lengths invalid")
source=read.delim(sourcefile,quote="",stringsAsFactors=FALSE,check.names=FALSE)
if(!"accepted_species" %in% names(source)||nrow(source)!=2234)
  stop("Original strict 2234 Direct source plant roster missing")
key=gsub(" ","_",trimws(source$accepted_species),fixed=TRUE)
tips=match(key,tree$tip.label)
valid=!is.na(tips)
result=data.frame(accepted_species=source$accepted_species[valid],
                  tip_node_id=as.integer(tips[valid]),stringsAsFactors=FALSE)
if(anyDuplicated(result$accepted_species)||anyDuplicated(result$tip_node_id))
  stop("Source accepted tip identity not unique")
edges=data.frame(parent_node_id=tree$edge[,1],
                 child_node_id=tree$edge[,2],
                 original_branch_length=tree$edge.length)
write.table(result,outTips,quote=FALSE,sep="\t",row.names=FALSE)
write.table(edges,outEdge,quote=FALSE,sep="\t",row.names=FALSE)
cat("ORIGINAL_TREE",version,"N_TIPS",length(tree$tip.label),
    "N_EDGES",nrow(edges),"ORIGINAL_DIRECT_TIP_CANDIDATES",nrow(result),"\n")
