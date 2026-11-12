source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
DEST <- file.path(OUT,'primary')
make_dirs(DEST,c('checkpoints','fits','tables','quality','logs'))
source(FUN('primary_setup.R'))
contract <- list(primary_formula=FORM_PRIMARY)
make_input <- function(){
 z<-setup();stopifnot(nrow(z$data)==2051L,sum(z$data$y)==998L,all(z$data$hist%in%c(0,1)));z
}
input<-make_input()
source_paths<-c(INPUT,file.path(CODE,'R/config.R'),FUN('imputation.R'),FUN('exact_envelope.R'),FUN('primary_setup.R'))
hashes<-vapply(source_paths,function(p)digest::digest(file=p,algo='sha256'),character(1))
signature<-digest::digest(hashes,algo='sha256')
jsonlite::write_json(list(N=2051,deaths=998,M=M_PRIMARY,iterations=ITERATIONS,formula=FORM_PRIMARY,
 history_retained_as_imputation_auxiliary=TRUE,PaCO2_knots=c(30,40,55),clinical_knots=KN,
 signature=signature),file.path(DEST,'quality/model_contract.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
job <- function(i){
  path <- file.path(DEST,'checkpoints',sprintf('imp_%03d.rds',i))
  seed <- as.integer(611050000L+i*7919L)
  if(file.exists(path)){
    z <- readRDS(path)
    stopifnot(z$signature==signature,z$seed==seed,identical(z$ids,input$sourceIDs))
    return(TRUE)
  }
  z <- prag_impute_envelope(input,m=1L,numit=40L,seed=seed,rjlimit=100000L,maxit=1000L)
  a <- z$exact_envelope_diagnostics
  stopifnot(isTRUE(a$all_returned_glms_converged),a$max_scaled_probability<=1+1e-8,
    a$max_offset_difference<1e-9,a$max_inactive_matrix_error<1e-9,
    length(z$rejection_failure_context)==0L)
  fatal <- grepl('not converge|not converged|failed to converge|rejection sampling',z$warnings,ignore.case=TRUE)
  expected <- length(z$warning_context)>0L && all(vapply(z$warning_context,function(w)
    identical(w$variable,'plow') && grepl('brglmFit: algorithm did not converge|brglmFit: fitted probabilities numerically 0 or 1',w$message),logical(1)))
  checked <- isTRUE(a$all_returned_glms_converged) && !is.null(a$constant_binary_intercept_checks) &&
    a$constant_binary_intercept_checks>0 && a$constant_binary_max_analytic_error<1e-6
  stopifnot(!any(fatal) || (expected && checked))
  saved <- list(d=z$impDatasets[[1]],trace=z$smCoefIter,warnings=z$warnings,
    warning_context=z$warning_context,exact=a,seed=seed,signature=signature,
    ids=input$sourceIDs,formula=contract$primary_formula,
    component_initialization_warnings_adjudicated=isTRUE(expected&&checked))
  saveRDS(saved,paste0(path,'.tmp'))
  stopifnot(file.rename(paste0(path,'.tmp'),path))
  cat('COMPLETE',i,round(z$seconds,1),'seconds; warnings',length(z$warnings),'\n');flush.console()
  TRUE
}
M <- M_PRIMARY
workers <- WORKERS
worker_log <- tempfile('primary_workers_',fileext='.txt')
writeLines(worker_log,file.path(DEST,'logs/workers_log_path.txt'))
cl <- parallel::makePSOCKcluster(workers,outfile=worker_log)
parallel::clusterExport(cl,c('CODE','DATA','OUT','INPUT','DEST','FORM_PRIMARY','FUN','KN','setup','make_input','signature','job','contract'))
invisible(parallel::clusterEvalQ(cl,{
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
  requireNamespace('smcfcs',quietly=TRUE)
  input <- make_input()
  NULL
}))
cat('START',M,'chains with',workers,'workers\n');flush.console()
res <- tryCatch(parallel::parLapplyLB(cl,seq_len(M),job,chunk.size=1L),finally=parallel::stopCluster(cl))
stopifnot(all(unlist(res)))
file.copy(worker_log,file.path(DEST,'logs/MI_workers.txt'),overwrite=TRUE)
writeLines(capture.output(sessionInfo()),file.path(DEST,'logs/R_session_info.txt'))
cat('ALL MI COMPLETE\n');flush.console()
source(file.path(CODE,'R/02_primary_models.R'))
