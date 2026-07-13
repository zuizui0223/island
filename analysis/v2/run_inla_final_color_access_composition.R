#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(INLA); library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) stop('usage: run_inla_final_color_access_composition.R <input.csv> <output_dir>')
infile <- args[[1]]; outdir <- args[[2]]; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(20260714)
d <- fread(infile)

colors <- c('white','blue','yellow','red','green','cream','rare_other')
access <- c('open_accessible','intermediate_bell_funnel','restricted_tubular','specialized_restricted','multi_flower_display','other')
endpoints <- list(
  color=list(role='primary', all_prefix='all_fine_color_', all_trials='all_fine_color_trials', animal_prefix='fine_color_', animal_trials='fine_color_trials', categories=colors),
  access=list(role='secondary_exploratory', all_prefix='all_access_', all_trials='all_access_trials', animal_prefix='animal_access_', animal_trials='animal_access_trials', categories=access)
)
regime_levels <- c('northern_midlatitude','northern_high_latitude','tropical','southern_extratropical')
d[, analysis_regime := factor(analysis_regime, levels=regime_levels)]
d[, wind_share := fifelse(animal_status_observed>0, n_wind_species/animal_status_observed, NA_real_)]
d[, mixed_share := fifelse(animal_status_observed>0, n_mixed_species/animal_status_observed, NA_real_)]

base_required <- c('island_id','analysis_regime','spatial_block','log_distance_to_continent_km','log_island_area_km2','climate_pc1','climate_pc2','climate_pc3','climate_pc4','bombus_deficit','sc_share','animal_status_observed','n_wind_species','n_mixed_species','sc_trials','sc_successes')
required <- unique(c(base_required, unlist(lapply(endpoints, function(x) c(x$all_trials,x$animal_trials,paste0(x$all_prefix,x$categories),paste0(x$animal_prefix,x$categories))))))
missing <- setdiff(required,names(d)); if(length(missing)) stop('missing final-model columns: ',paste(missing,collapse=', '))
for(ep in names(endpoints)){
  x <- endpoints[[ep]]; ac <- paste0(x$all_prefix,x$categories); nc <- paste0(x$animal_prefix,x$categories)
  if(any(rowSums(d[,..ac],na.rm=TRUE)!=fifelse(is.na(d[[x$all_trials]]),0,d[[x$all_trials]]))) stop(ep,' all-flora categories do not sum to trials')
  if(any(rowSums(d[,..nc],na.rm=TRUE)!=fifelse(is.na(d[[x$animal_trials]]),0,d[[x$animal_trials]]))) stop(ep,' animal-flora categories do not sum to trials')
}

for(nm in c('log_distance_to_continent_km','log_island_area_km2','climate_pc1','climate_pc2','climate_pc3','climate_pc4')){
  mu <- mean(d[[nm]],na.rm=TRUE); sig <- sd(d[[nm]],na.rm=TRUE); if(!is.finite(sig)||sig<=0) stop('invalid scale for ',nm)
  d[[paste0('z_',nm)]] <- (d[[nm]]-mu)/sig
}
set(d,j='z_isolation_sq',value=d$z_log_distance_to_continent_km^2)

north <- d[analysis_regime=='northern_midlatitude']
mech_candidates <- c('bombus_deficit','sc_share','wind_share','mixed_share')
mech_audit <- rbindlist(lapply(mech_candidates,function(nm){x<-north[[nm]]; n<-sum(is.finite(x)); u<-uniqueN(x[is.finite(x)]); s<-sd(x,na.rm=TRUE); ok<-n>=2&&u>=2&&is.finite(s)&&s>0; data.table(variable=nm,n_nonmissing=n,n_unique=u,sd=s,estimable=ok,reason=if(ok)'estimable' else 'zero_or_insufficient_variation')}))
fwrite(mech_audit,file.path(outdir,'final_mechanism_variable_audit.csv'))
for(nm in mech_audit[estimable==TRUE,variable]){mu<-mean(north[[nm]],na.rm=TRUE); s<-sd(north[[nm]],na.rm=TRUE); d[[paste0('z_',nm)]]<-(d[[nm]]-mu)/s}
if(!'z_bombus_deficit'%in%names(d)) stop('bombus_deficit not estimable in northern mechanism population')
if(!'z_sc_share'%in%names(d)) stop('sc_share not estimable in northern mechanism population')
poll_controls <- mech_audit[variable%in%c('wind_share','mixed_share') & estimable==TRUE,variable]

compute <- list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE)
fixed_rows<-list(); score_rows<-list(); support_rows<-list()
fit_comp <- function(data,endpoint,role,model_name,prefix,trials,categories,rhs){
  cols<-paste0(prefix,categories); ids<-intersect(c('island_id','spatial_block','analysis_regime',trials,'z_log_distance_to_continent_km','z_isolation_sq','z_log_island_area_km2','z_climate_pc1','z_climate_pc2','z_climate_pc3','z_climate_pc4','z_bombus_deficit','z_sc_share','z_wind_share','z_mixed_share'),names(data))
  long<-melt(data,id.vars=ids,measure.vars=cols,variable.name='category',value.name='count'); long[,category:=factor(sub(paste0('^',prefix),'',category),levels=categories)]; long[,island_id_re:=as.integer(factor(island_id))]; long[,spatial_id:=as.integer(factor(spatial_block))]; long[,obs_id:=seq_len(.N)]; long[,log_trials:=log(get(trials))]
  fit<-inla(as.formula(paste('count ~',rhs,'+ offset(log_trials)')),family='poisson',data=long,control.compute=compute,control.predictor=list(compute=TRUE),verbose=FALSE)
  fx<-as.data.table(fit$summary.fixed,keep.rownames='parameter'); fx[,`:=`(endpoint=endpoint,role=role,model=model_name)]; fixed_rows[[length(fixed_rows)+1L]]<<-fx
  cpo<-fit$cpo$cpo; cpo<-cpo[is.finite(cpo)&cpo>0]; score_rows[[length(score_rows)+1L]]<<-data.table(endpoint=endpoint,role=role,model=model_name,n_islands=uniqueN(data$island_id),n_rows=nrow(long),log_cpo_sum=if(length(cpo))sum(log(cpo)) else NA_real_,waic=fit$waic$waic,dic=fit$dic$dic)
}
geo_terms <- 'category:z_log_island_area_km2 + category:z_climate_pc1 + category:z_climate_pc2 + category:z_climate_pc3 + category:z_climate_pc4'
random_terms <- "+ f(island_id_re, model='iid') + f(spatial_id, model='iid') + f(obs_id, model='iid')"

for(ep in names(endpoints)){
  x<-endpoints[[ep]]
  global_support<-d[!is.na(analysis_regime)&get(x$all_trials)>0&complete.cases(d[,.(z_log_distance_to_continent_km,z_isolation_sq,z_log_island_area_km2,z_climate_pc1,z_climate_pc2,z_climate_pc3,z_climate_pc4)])]
  global_rhs<-paste('0 + category + category:analysis_regime + category:analysis_regime:z_log_distance_to_continent_km + category:analysis_regime:z_isolation_sq +',geo_terms,random_terms)
  fit_comp(global_support,paste0('all_flora_',ep),x$role,'G_global_isolation_by_regime_no_bombus',x$all_prefix,x$all_trials,x$categories,global_rhs)
  support_rows[[length(support_rows)+1L]]<-data.table(endpoint=paste0('all_flora_',ep),role=x$role,layer='global_all_flora_no_bombus',n_islands=nrow(global_support))

  base_vars<-c('z_log_distance_to_continent_km','z_isolation_sq','z_log_island_area_km2','z_climate_pc1','z_climate_pc2','z_climate_pc3','z_climate_pc4')
  base_rhs<-paste('0 + category + category:z_log_distance_to_continent_km + category:z_isolation_sq +',geo_terms,random_terms)

  if(ep=='color'){
    mech_vars<-c(base_vars,'z_bombus_deficit','z_sc_share',paste0('z_',poll_controls))
    ns<-d[analysis_regime=='northern_midlatitude'&get(x$animal_trials)>0&complete.cases(d[,..mech_vars])]
    models<-list(
      N0_isolation=base_rhs,
      N1_bombus=paste(base_rhs,'+ category:z_bombus_deficit'),
      N2_bombus_sc=paste(base_rhs,'+ category:z_bombus_deficit + category:z_sc_share')
    )
    extra<-c('category:z_bombus_deficit','category:z_sc_share',paste0('category:z_',poll_controls)); models$N3_full<-paste(base_rhs,'+',paste(extra,collapse=' + '))
    for(nm in names(models)) fit_comp(ns,'animal_flora_color','primary',nm,x$animal_prefix,x$animal_trials,x$categories,models[[nm]])
    support_rows[[length(support_rows)+1L]]<-data.table(endpoint='animal_flora_color',role='primary',layer='northern_bombus_color_mechanism',n_islands=nrow(ns))
  } else {
    # Shape/access structure is secondary: diagnose isolation and pollination-mode composition only.
    sec_vars<-c(base_vars,paste0('z_',poll_controls))
    ns<-d[analysis_regime=='northern_midlatitude'&get(x$animal_trials)>0&complete.cases(d[,..sec_vars])]
    models<-list(S0_isolation=base_rhs)
    if(length(poll_controls)) models$S1_pollination_mode_adjusted<-paste(base_rhs,'+',paste(paste0('category:z_',poll_controls),collapse=' + '))
    for(nm in names(models)) fit_comp(ns,'animal_flora_access','secondary_exploratory',nm,x$animal_prefix,x$animal_trials,x$categories,models[[nm]])
    support_rows[[length(support_rows)+1L]]<-data.table(endpoint='animal_flora_access',role='secondary_exploratory',layer='northern_access_pollination_mode_diagnostic',n_islands=nrow(ns))
  }
}

# Explicit SC equation only in the northern Bombus-relevant regime.
sc_vars<-c('z_log_distance_to_continent_km','z_isolation_sq','z_log_island_area_km2','z_climate_pc1','z_climate_pc2','z_climate_pc3','z_climate_pc4','z_bombus_deficit',paste0('z_',poll_controls))
scs<-d[analysis_regime=='northern_midlatitude'&sc_trials>0&!is.na(sc_successes)&complete.cases(d[,..sc_vars])]; scs[,spatial_id_sc:=as.integer(factor(spatial_block))]
sc_rhs<-c('z_log_distance_to_continent_km','z_isolation_sq','z_bombus_deficit',paste0('z_',poll_controls),'z_log_island_area_km2','z_climate_pc1','z_climate_pc2','z_climate_pc3','z_climate_pc4',"f(spatial_id_sc, model='iid')")
sc_fit<-inla(as.formula(paste('sc_successes ~',paste(sc_rhs,collapse=' + '))),family='betabinomial',data=scs,Ntrials=sc_trials,control.compute=compute,control.predictor=list(compute=TRUE),verbose=FALSE)
sc_fx<-as.data.table(sc_fit$summary.fixed,keep.rownames='parameter'); sc_fx[,`:=`(endpoint='self_compatibility',role='mechanism',model='SC_bombus_reproductive_assurance_north_only')]; fixed_rows[[length(fixed_rows)+1L]]<-sc_fx

fixed<-rbindlist(fixed_rows,fill=TRUE); fixed[,excludes_zero:=(`0.025quant`>0&`0.975quant`>0)|(`0.025quant`<0&`0.975quant`<0)]
scores<-rbindlist(score_rows,fill=TRUE); support<-rbindlist(support_rows,fill=TRUE)
fwrite(fixed,file.path(outdir,'final_color_access_fixed_effects.csv')); fwrite(scores,file.path(outdir,'final_color_access_model_comparisons.csv')); fwrite(support,file.path(outdir,'final_color_access_support.csv')); fwrite(sc_fx,file.path(outdir,'final_sc_pathway.csv'))
key<-fixed[grepl('z_log_distance_to_continent_km|z_isolation_sq|z_bombus_deficit|z_sc_share|z_wind_share|z_mixed_share',parameter)]; fwrite(key,file.path(outdir,'final_color_access_key_effects.csv'))
fwrite(key[role=='primary'],file.path(outdir,'primary_color_key_effects.csv'))
fwrite(key[role=='secondary_exploratory'],file.path(outdir,'secondary_access_key_effects.csv'))
regime_key<-fixed[model=='G_global_isolation_by_regime_no_bombus' & grepl('analysis_regime.*z_(log_distance_to_continent_km|isolation_sq)$',parameter)]
fwrite(regime_key,file.path(outdir,'global_regime_isolation_effects_no_bombus.csv'))
fwrite(regime_key[endpoint=='all_flora_color'],file.path(outdir,'primary_global_color_regime_effects.csv'))
fwrite(regime_key[endpoint=='all_flora_access'],file.path(outdir,'secondary_global_access_regime_effects.csv'))
meta<-list(
  canonical_primary_endpoint='fine flower-colour composition',
  primary_color_categories=colors,
  primary_global_model='all trait-resolved flora; isolation by regime; no Bombus outside the northern mechanism model',
  primary_northern_mechanism='animal-pollinated colour composition: isolation + Bombus deficit + SC + estimable pollination-mode controls',
  self_compatibility_pathway='explicit northern SC equation retained as a parallel reproductive-assurance mechanism test',
  floral_access_role='secondary exploratory diagnostic only',
  floral_access_interpretation='not used as evidence for Bombus-mediated selection; evaluated against isolation and estimable wind/mixed composition because form may track pollination-mode/floristic composition',
  tropical_southern_interpretation='global colour isolation effects estimated without Bombus terms because Bombus is not the regional baseline channel',
  joint_colour_form='sensitivity only'
)
write_json(meta,file.path(outdir,'final_color_access_metadata.json'),pretty=TRUE,auto_unbox=TRUE)
