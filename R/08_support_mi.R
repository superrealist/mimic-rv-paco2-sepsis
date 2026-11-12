source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
O<-file.path(OUT,'support')
make_dirs(O,c('completed','quality','tables'))
args<-commandArgs(TRUE);endpoints<-if(length(args))args else SUPPORT_ENDPOINTS
stopifnot(all(endpoints%in%SUPPORT_ENDPOINTS))
workers<-WORKERS
sources<-setNames(file.path(DATA,'support',paste0(SUPPORT_ENDPOINTS,'_input.csv')),SUPPORT_ENDPOINTS)
formula_string<-FORM_CATEGORY
setup<-function(ep){
 source(FUN('imputation.R'))
 source(FUN('numerical_envelope.R'))
 source(FUN('category_setup.R'))
 source(FUN('support_imputation.R'))
 inp<-prag_read_input(sources[[ep]],KN)
 inp$data$y<-as.numeric(inp$source$analysis_y)
 inp$source$death28<-inp$source[[if(grepl('28$',ep))'death28_echo' else 'death7_echo']]
 if(!grepl('28$',ep)) {
   altered<-paste(deparse(prepare_vfd_mice),collapse='\n')
   altered<-gsub('c(0, 14, 28)','c(0, 3.5, 7)',altered,fixed=TRUE)
   prepare_vfd_mice<-eval(parse(text=altered))
 }
 inp$data$co2lo<-as.integer(inp$data$x<35);inp$data$co2hi<-as.integer(inp$data$x>45)
 prag_generic<<-unique(c(prag_generic,'co2lo','co2hi'));prag_formula<<-formula_string;prag_knots<<-KN
 inp$data$lcount<-inp$source$recorded_L_drug_count;prag_generic<<-unique(c(prag_generic,'lcount'))
 p<-prepare_vfd_mice(inp,'binary',FALSE)
 # Add observed category indicators and RV interactions to the
 # outcome-enriched bounded donor imputation specification.
 for(t in c('co2lo','co2hi')){
   nm<-paste0('rv__',t);p$data[[nm]]<-p$data$rv*p$data[[t]]
   p$method[nm]<-'';p$where<-cbind(p$where,FALSE);colnames(p$where)[ncol(p$where)]<-nm
   p$predictorMatrix<-rbind(cbind(p$predictorMatrix,0L),0L);dimnames(p$predictorMatrix)<-list(names(p$data),names(p$data))
   active<-names(p$method)[p$method=='pmm'];p$predictorMatrix[active,c(t,nm)]<-1L
 }
 inp$prepared<-p;inp
}
job<-function(i){
 ep<-EP;fn<-file.path(O,'completed',sprintf('%s_cat_imp_%03d.rds',ep,i))
 if(file.exists(fn))return(TRUE)
 seed<-as.integer(630000000L+match(ep,names(sources))*1000000L+7919L*i)
 warn<-character();logged<-NULL
   p<-input$prepared
   z<-withCallingHandlers(mice::mice(p$data,m=1L,maxit=40L,method=p$method,predictorMatrix=p$predictorMatrix,visitSequence=p$visitSequence,
        where=p$where,donors=p$donors,matchtype=p$matchtype,ridge=p$ridge,remove.constant=FALSE,remove.collinear=FALSE,printFlag=FALSE,seed=seed),
        warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')})
   stopifnot(identical(z$method,p$method),all(z$predictorMatrix==p$predictorMatrix))
   d<-mice::complete(z);logged<-z$loggedEvents;trace<-z$chainMean
 prag_validate(d,input)
 stopifnot(all(d$y==input$source$analysis_y),all(d$co2lo==(d$x<35)),all(d$co2hi==(d$x>45)),is.null(logged)||!nrow(logged))
 saveRDS(list(d=d,trace=trace,warning=warn,logged=logged,seed=seed,ids=input$sourceIDs),fn);TRUE
}
for(EP in endpoints){
 logfile<-tempfile('secondary_cat_workers_',fileext='.log')
 cl<-parallel::makePSOCKcluster(workers,outfile=logfile)
 parallel::clusterExport(cl,c('CODE','O','EP','FUN','KN','sources','setup','job','formula_string'))
 invisible(parallel::clusterEvalQ(cl,{input<-setup(EP);NULL}))
 results<-tryCatch(parallel::parLapplyLB(cl,1:40,job,chunk.size=1),finally=parallel::stopCluster(cl))
 pieces<-lapply(1:40,function(i)readRDS(file.path(O,'completed',sprintf('%s_cat_imp_%03d.rds',EP,i))))
 z<-list(impDatasets=lapply(pieces,'[[','d'),sourceIDs=pieces[[1]]$ids,formula=formula_string)
 stopifnot(all(vapply(pieces,function(s)identical(s$ids,z$sourceIDs),logical(1))))
 saveRDS(z,file.path(O,'completed',paste0(EP,'_imputations.rds')))
 warn<-unlist(lapply(pieces,'[[','warning'))
 # Examine individual tracked variable/coefficient traces across chains from the last 20 iterations.
 mats<-lapply(pieces,function(v)t(v$trace[,21:40,1,drop=FALSE][,,1]))
 rh<-vapply(seq_len(ncol(mats[[1]])),function(k){
   v<-sapply(mats,function(m)m[,k]);if(any(!is.finite(v))||sd(as.numeric(v))<1e-12)NA_real_ else posterior::rhat(v)
 },numeric(1))
 jsonlite::write_json(list(endpoint=EP,M=40,iterations=40,N=length(z$sourceIDs),source_sha256=digest::digest(file=sources[[EP]],algo='sha256'),
             formula=formula_string,outcome_unchanged=TRUE,observed_values_validated=TRUE,max_Rhat=max(rh,na.rm=TRUE),
             warning_count=length(warn),warnings=unique(warn)),file.path(O,'quality',paste0(EP,'_mi.json')),pretty=TRUE,auto_unbox=TRUE)
}
