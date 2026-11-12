source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
OUT<-file.path(OUT,'category')
make_dirs(OUT,c('checkpoints','tables','logs'))
ct<-list(categorical_formula=FORM_CATEGORY_IMPUTATION,seed_base=202611350L,M=M_CATEGORY,workers=WORKERS)
MODULES<-c(INPUT,file.path(CODE,'R/config.R'),FUN('imputation.R'),FUN('exact_envelope.R'),FUN('category_setup.R'))
HASHES<-vapply(MODULES,function(p)digest::digest(file=p,algo='sha256'),character(1))
setup<-function(){
 source(FUN('imputation.R'))
 source(FUN('exact_envelope.R'))
 source(FUN('category_setup.R'))
 inp<-cat_read_input(INPUT,KN);prag_formula<<-ct$categorical_formula
 cat_setup();inp
}
job<-function(i){
 inp<-setup();seed<-as.integer(ct$seed_base+7919L*i)
 folder<-file.path(OUT,'checkpoints');dir.create(folder,showWarnings=FALSE)
 path<-file.path(folder,sprintf('imp_%03d.rds',i))
 if(file.exists(path)){
  z<-readRDS(path);stopifnot(identical(z$metadata$source_hashes,HASHES),z$seed==seed)
 }else{
  z<-prag_impute_envelope(inp,m=1L,numit=40L,seed=seed,rjlimit=100000L,maxit=1000L)
  z$metadata$source_hashes<-HASHES;z$metadata$imputation_id<-i
 }
 e<-z$exact_envelope_diagnostics;w<-z$warning_context
 expected<-length(w)>0L&&all(vapply(w,function(a)identical(a$variable,'plow')&&
  grepl('brglmFit: algorithm did not converge|brglmFit: fitted probabilities numerically 0 or 1',a$message),logical(1)))
 checked<-isTRUE(e$all_returned_glms_converged)&&!is.null(e$constant_binary_intercept_checks)&&
  e$constant_binary_intercept_checks>0&&e$constant_binary_max_analytic_error<1e-6
 fatal<-grepl('not converge|not converged|failed to converge|rejection sampling',z$warnings,ignore.case=TRUE)
 stopifnot(!any(fatal)||(expected&&checked&&length(z$rejection_failure_context)==0L))
 z$metadata$constant_initialization_warning_adjudicated<-isTRUE(expected&&checked)
 stopifnot(isTRUE(e$all_returned_glms_converged),e$max_scaled_probability<=1+1e-8,
  e$max_offset_difference<1e-9,e$max_inactive_matrix_error<1e-9,length(z$rejection_failure_context)==0L)
 prag_validate(z$impDatasets[[1]],inp)
 for(v in c('co2lo','co2hi'))stopifnot(identical(z$impDatasets[[1]][[v]],inp$data[[v]]))
 saveRDS(z,path);cat('MI',i,'seconds',round(z$seconds,1),'warnings',length(z$warnings),'PASS\n');flush.console();TRUE
}
args<-commandArgs(TRUE);ids<-if('--pilot'%in%args)1L else seq_len(ct$M)
cl<-parallel::makePSOCKcluster(if('--pilot'%in%args)1L else ct$workers,
 outfile=file.path(OUT,'logs/MI_workers.txt'))
parallel::clusterExport(cl,c('CODE','OUT','INPUT','FUN','ct','KN','HASHES','setup','job'),envir=environment())
invisible(parallel::clusterEvalQ(cl,{Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');requireNamespace('smcfcs',quietly=TRUE);NULL}))
started<-Sys.time();ans<-tryCatch(parallel::parLapplyLB(cl,ids,job,chunk.size=1L),finally=parallel::stopCluster(cl))
stopifnot(all(unlist(ans)))
if('--pilot'%in%args){cat('CATEGORICAL PILOT PASS\n');quit(save='no')}
inp<-setup();pieces<-lapply(ids,function(i)readRDS(file.path(OUT,'checkpoints',sprintf('imp_%03d.rds',i))))
z<-pieces[[1]];z$impDatasets<-lapply(pieces,function(p)p$impDatasets[[1]])
z$smCoefIter<-abind::abind(lapply(pieces,`[[`,'smCoefIter'),along=1)
z$metadata$M<-length(ids);z$metadata$ids<-ids;z$metadata$seeds<-vapply(pieces,`[[`,integer(1),'seed')
z$warnings<-unlist(lapply(pieces,`[[`,'warnings'));invisible(lapply(z$impDatasets,prag_validate,inp=inp))
path<-file.path(OUT,'imputations.rds');stopifnot(!file.exists(path));saveRDS(z,path)
writeLines(capture.output(sessionInfo()),file.path(OUT,'logs/R_session_info.txt'))
jsonlite::write_json(list(status='CATEGORICAL_MI_COMPLETE',M=length(ids),N=nrow(inp$data),iterations=40,
  elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),sha256=digest::digest(file=path,algo='sha256')),
 file.path(OUT,'logs/MI_complete.json'),auto_unbox=TRUE,pretty=TRUE)
cat('CATEGORICAL MI COMPLETE\n')
