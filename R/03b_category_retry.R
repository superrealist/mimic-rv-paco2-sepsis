source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
# Conditional recovery: requires exactly one missing categorical chain.
p<-file.path(CODE,'R/03_categorical_mi.R')
ll<-readLines(p);cut<-grep('^args<-',ll)[1]-1L
eval(parse(text=paste(ll[seq_len(cut)],collapse='\n')),envir=.GlobalEnv)
missing<-which(!file.exists(file.path(OUT,'checkpoints',sprintf('imp_%03d.rds',1:100))))
stopifnot(length(missing)==1L)
inp<-setup();source(FUN('numerical_envelope.R'))
i<-missing;seed<-as.integer(ct$seed_base+7919L*i)
z<-prag_impute_envelope(inp,m=1L,numit=40L,seed=seed,rjlimit=100000L,maxit=1000L)
stopifnot(isTRUE(z$exact_envelope_diagnostics$all_returned_glms_converged),length(z$rejection_failure_context)==0L)
prag_validate(z$impDatasets[[1]],inp)
z$metadata$source_hashes<-HASHES;z$metadata$imputation_id<-i
z$metadata$numerical_repair<-list(reason='Pzero bias-reduced nuisance regression did not meet convergence at1000 iterations.',
 same_seed=seed,same_formula=TRUE,same_covariate_conditional_model=TRUE,
 deterministic_retry='For a nonconverged brglmFit only: same target,start at previous finite coefficients,slowit0.25,maxit10000; no extra RNG calls.',
 repair_script_sha256=digest::digest(file=file.path(CODE,'R/03b_category_retry.R'),algo='sha256'),
 envelope_sha256=digest::digest(file=FUN('numerical_envelope.R'),algo='sha256'))
saveRDS(z,file.path(OUT,'checkpoints',sprintf('imp_%03d.rds',i)))
ids<-1:100;pieces<-lapply(ids,function(k)readRDS(file.path(OUT,'checkpoints',sprintf('imp_%03d.rds',k))))
for(a in pieces)stopifnot(identical(a$metadata$source_hashes,HASHES),isTRUE(a$exact_envelope_diagnostics$all_returned_glms_converged),length(a$rejection_failure_context)==0L)
combined<-pieces[[1]];combined$impDatasets<-lapply(pieces,function(a)a$impDatasets[[1]])
combined$smCoefIter<-abind::abind(lapply(pieces,`[[`,'smCoefIter'),along=1L)
combined$warnings<-unlist(lapply(pieces,`[[`,'warnings'))
combined$metadata$M<-100L;combined$metadata$ids<-ids
combined$metadata$seeds<-vapply(pieces,`[[`,integer(1),'seed')
combined$metadata$numerical_repair<-z$metadata$numerical_repair
invisible(lapply(combined$impDatasets,prag_validate,inp=inp))
path<-file.path(OUT,'imputations.rds');stopifnot(!file.exists(path));saveRDS(combined,path)
jsonlite::write_json(list(status='CATEGORICAL_MI_COMPLETE',M=100,N=nrow(inp$data),iterations=40,
 uniform_numerical_fallback_target_unchanged=TRUE,repaired_id=i,repair=z$metadata$numerical_repair,
 all_final_glms_converged=TRUE,warning_count=length(combined$warnings),
 sha256=digest::digest(file=path,algo='sha256')),file.path(OUT,'logs/MI_complete.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('ALL100 CATEGORY IMPUTATIONS COMPLETE; numerical repair documented\n')
