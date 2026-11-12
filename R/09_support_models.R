source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
O<-file.path(OUT,'support')
make_dirs(O,c('tables','fits','quality'))
source(FUN('pooling.R'))
spec<-as.formula(FORM_CATEGORY)
allrows<-alltests<-list()
for(ep in c('vpfd7','vfd7','crrtfd7','vpfd28','vfd28','crrtfd28')){
  path<-file.path(O,'completed',paste0(ep,'_imputations.rds'))
  z<-readRDS(path);raw<-read.csv(file.path(DATA,'support',paste0(ep,'_input.csv')))
  stopifnot(identical(as.numeric(z$sourceIDs),as.numeric(raw$subject_id)))
  fits<-lapply(z$impDatasets,function(d){
    stopifnot(all(d$y==raw$analysis_y),all(d$rv==raw$rv_binary),all(d$x==raw$pco2))
    f<-lm(spec,d)
    stopifnot(f$rank==length(coef(f)),all(is.finite(coef(f))),nobs(f)==nrow(raw))
    f
  })
  Q<-sapply(fits,coef);nm<-rownames(Q);m<-ncol(Q);p<-nrow(Q)
  U<-array(NA_real_,c(p,p,m),dimnames=list(nm,nm,NULL))
  for(i in 1:m)U[,,i]<-sandwich::vcovHC(fits[[i]],type='HC3')
  stopifnot(all(is.finite(U)),all(vapply(1:m,function(i)min(eigen(U[,,i],symmetric=TRUE,only.values=TRUE)$values)>0,logical(1))))
  sc<-expand.grid(category=c('low','high'),rv=0:1)
  H<-matrix(0,4,p,dimnames=list(NULL,nm))
  for(k in 1:4){
    term<-if(sc$category[k]=='low')'co2lo' else 'co2hi'
    H[k,term]<-1
    if(sc$rv[k])H[k,paste0('rv:',term)]<-1
  }
  pp<-do.call(rbind,lapply(1:4,function(k){
    qq<-drop(H[k,,drop=FALSE]%*%Q)
    uu<-vapply(1:m,function(i)drop(H[k,,drop=FALSE]%*%U[,,i]%*%H[k,]),numeric(1))
    pi_pool_scalar(qq,uu,df.residual(fits[[1]]))
  }))
  idx<-match(c('rv:co2lo','rv:co2hi'),nm)
  test<-pi_d1_q_u(Q[idx,,drop=FALSE],U[idx,idx,,drop=FALSE],df.residual(fits[[1]]))
  # Independent explicit standardization verifies the linear coefficient contrasts.
  for(i in c(1,m))for(k in 1:4){
    d<-z$impDatasets[[i]];d$rv<-sc$rv[k];d$co2lo<-0;d$co2hi<-0
    ref<-predict(fits[[i]],newdata=d)
    d$co2lo<-as.integer(sc$category[k]=='low');d$co2hi<-as.integer(sc$category[k]=='high')
    diff<-mean(predict(fits[[i]],newdata=d)-ref)
    stopifnot(abs(diff-drop(H[k,,drop=FALSE]%*%Q[,i]))<1e-9)
  }
  result<-cbind(outcome=ep,sc,n=nrow(raw),pp[,c('estimate','lower','upper','p.value')],interaction_p=test$p.value)
  allrows[[ep]]<-result
  alltests[[ep]]<-data.frame(outcome=ep,N=nrow(raw),M=m,parameters=p,df1=test$df1,df2=test$df2,statistic=test$statistic,interaction_p=test$p.value,max_MCSE_SE=max(pp$mi_mcse/pp$std.error))
  write.csv(result,file.path(O,'tables',paste0('forest_',ep,'.csv')),row.names=FALSE)
  saveRDS(list(Q=Q,U=U,H=H,dfcom=df.residual(fits[[1]]),source_ids=raw$subject_id),file.path(O,'fits',paste0(ep,'_fit.rds')))
}
write.csv(do.call(rbind,allrows),file.path(O,'tables/forest.csv'),row.names=FALSE)
write.csv(do.call(rbind,alltests),file.path(O,'quality/model_checks.csv'),row.names=FALSE)
