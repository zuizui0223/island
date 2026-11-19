#!/usr/bin/env Rscript
# Private original Direct plant species × source V.PhyloMaker2 backbone tip census.
# Data columns containing original accepted names NEVER go to CI artifacts/logs.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("Usage: Rscript private_tip_match.R backbone_data input.tsv output.tsv")
root <- args[[1]]
inp <- args[[2]]
out <- args[[3]]
required <- c("accepted_species","source_family","source_genus","rarity_tier")
source <- read.delim(inp,quote="",stringsAsFactors=FALSE,
                     check.names=FALSE,fileEncoding="UTF-8")
if (!all(required %in% names(source)) || anyDuplicated(source$accepted_species))
  stop("Frozen original Direct input schema/accepted taxa invalid")
source$normalized_tip <- gsub(" ","_",trimws(source$accepted_species),fixed=TRUE)
if (anyDuplicated(source$normalized_tip))
  stop("Original normalized species tips not unique")
if (nrow(source)!=2234) stop("2234 source species Direct input drifted")

all_result <- list()
for (label in c("TPL","LCVP","WP")) {
  e <- new.env(parent=baseenv())
  loaded_tree <- load(file.path(root,paste0("GBOTB.extended.",label,".rda")),envir=e)
  loaded_tips <- load(file.path(root,paste0("tips.info.",label,".rda")),envir=e)
  tree <- e[[loaded_tree[1]]]
  info <- e[[loaded_tips[1]]]
  if (!all(c("species","genus","family") %in% names(info)))
    stop("Original dated backbone tips.info format changed")
  if (length(tree$tip.label)!=nrow(info) || anyNA(tree$edge.length) ||
      any(tree$edge.length<=0))
    stop("Tree backbone tip/branch-length provenance or coherence invalid")
  if (anyDuplicated(tree$tip.label)) stop("Duplicated backbone tip labels")
  named_species <- as.character(info$species)
  # tips.info species expected standardized canonical underscored names.
  normalized_info <- gsub(" ","_",trimws(named_species),fixed=TRUE)
  if (length(unique(normalized_info))!=length(normalized_info))
    stop("Ambiguous original backbone info tip list")
  # Source tree and genus names must be internally congruent even if a
  # vendor's extra taxonomic synonyms occur; never accept unverified matches.
  tip_ids <- match(source$normalized_tip,tree$tip.label)
  info_idx <- match(source$normalized_tip,normalized_info)
  exact <- !is.na(tip_ids) & !is.na(info_idx)
  family_ok <- exact & !is.na(source$source_family) &
     nzchar(source$source_family) &
     tolower(trimws(source$source_family)) ==
       tolower(trimws(as.character(info$family)[info_idx]))
  genus_ok <- exact & !is.na(source$source_genus) &
     nzchar(source$source_genus) &
     tolower(trimws(source$source_genus)) ==
       tolower(trimws(as.character(info$genus)[info_idx]))
  exact_ok <- family_ok & genus_ok
  # A taxonomy candidate in this dataset is unreviewed, so require both
  # independently curated source family and backbone tip agreement.
  known <- !is.na(source$source_family) & nzchar(source$source_family)
  all_genera <- unique(paste(tolower(as.character(info$genus)),
                       tolower(as.character(info$family)),sep="||"))
  all_families <- unique(tolower(as.character(info$family)))
  potential_genus <- !exact_ok & known & paste(
      tolower(source$source_genus),
      tolower(source$source_family),sep="||") %in% all_genera
  potential_family <- !exact_ok & known & !potential_genus &
      tolower(source$source_family) %in% all_families
  possible <- ifelse(exact_ok,"exact_backbone_tip",
     ifelse(potential_genus,"genus_backbone_insertion_only",
       ifelse(potential_family,"family_backbone_insertion_only",
         "unresolved_no_valid_source_family")))
  # Separate exact *species tree name* from conflicting family assignment.
  # Never silently call conflict exact phylogenetic placement.
  all_result[[label]] <- data.frame(
      accepted_species=source$accepted_species,
      rarity_tier=source$rarity_tier,
      tree_version=label,
      backbone_tip_name_exact=!is.na(tip_ids),
      backbone_tip_and_info_names_exact=exact,
      exact_source_genus_and_family_agree=exact_ok,
      exact_tree_tip_family_disagrees=exact & known & !exact_ok,
      supported_species_placement=possible,
      stringsAsFactors=FALSE)
  cat("Original source",label,"n_tips",length(tree$tip.label),
      "original_exact_backbone_names",sum(exact),
      "strict_family_genus_matched_tips",sum(exact_ok),
      "genus_only_unresolved_tips",sum(potential_genus),
      "family_only_unresolved_tips",sum(potential_family),"\n")
}
result <- do.call(rbind,all_result)
write.table(result,file=out,quote=FALSE,sep="\t",row.names=FALSE)
