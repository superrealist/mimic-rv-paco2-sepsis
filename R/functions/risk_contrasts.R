# Marginal risk and contrast calculations from frozen fits; no model refitting.
fixed_rcs<-function(x,k=c(30,40,55)){
 (pmax(x-k[1],0)^3-(k[3]-k[1])/(k[3]-k[2])*pmax(x-k[2],0)^3+(k[2]-k[1])/(k[3]-k[2])*pmax(x-k[3],0)^3)/(k[3]-k[1])^2
}
risk_one_completion<-function(job){
 d<-job$data;b<-job$b;V<-job$V;sc<-job$sc;nm<-names(b);n<-nrow(d);p<-length(b)
 X<-model.matrix(as.formula(job$formula),d)[,nm,drop=FALSE]
 stopifnot(job$encoding=='mimic_raw')
 exposure<-c('rv','x','xnl','rv:x','rv:xnl')
 j<-match(exposure,nm);stopifnot(!anyNA(j))
 nd<-d[rep(1,nrow(sc)),,drop=FALSE];nd$rv<-sc$rv
 nd$x<-sc$pco2;nd$xnl<-fixed_rcs(sc$pco2)
 A<-model.matrix(as.formula(job$formula),nd)[,nm,drop=FALSE];A[,-j]<-0
 X[,j]<-0
 pr<-plogis(outer(drop(X%*%b),drop(A%*%b),'+'));w<-pr*(1-pr)
 q<-colMeans(pr)
 G<-t(crossprod(X,w))/n+A*rep(colMeans(w),times=p)
 u<-rowSums((G%*%V)*G)
 stopifnot(all(is.finite(c(q,G,u))),all(q>0&q<1),all(u>=0))
 list(q=q,u=u,G=G)
}
scalar_pool_table<-function(Q,U,dfcom){
 out<-do.call(rbind,lapply(seq_len(nrow(Q)),function(k)pi_pool_scalar(Q[k,],U[k,],dfcom)))
 data.frame(estimate=out$estimate,se=out$std.error,lower=out$lower,upper=out$upper,df=out$df,p_value=out$p.value,
            within_variance=out$ubar,between_variance=out$b,total_variance=out$total_variance)
}
contrast_design<-function(sc){
 n<-nrow(sc);at<-function(g,x){k<-which(sc$rv==g & abs(sc$pco2-x)<1e-8);stopifnot(length(k)==1);k}
 one<-function(g,a,b){z<-numeric(n);z[at(g,b)]<-z[at(g,b)]+1;z[at(g,a)]<-z[at(g,a)]-1;z}
 relative<-do.call(rbind,lapply(seq_len(nrow(sc)),function(k)one(sc$rv[k],40,sc$pco2[k])))
 starts<-seq(25,55,by=.1);slide<-list();sm<-list();fixed<-list();fm<-list()
 add<-function(a,b){z0<-one(0,a,b);z1<-one(1,a,b);rbind(change_rv0=z0,change_rv1=z1,interaction=z1-z0)}
 for(a in starts){z<-add(a,a+10);slide[[length(slide)+1]]<-z;sm[[length(sm)+1]]<-data.frame(x_start=a,x_end=a+10,metric=rownames(z))}
 for(v in list(c(30,40),c(40,50),c(50,60),c(40,60))){z<-add(v[1],v[2]);fixed[[length(fixed)+1]]<-z;fm[[length(fm)+1]]<-data.frame(start=v[1],end=v[2],metric=rownames(z))}
 list(relative_L=relative,relative_meta=sc,sliding_L=do.call(rbind,slide),sliding_meta=do.call(rbind,sm),fixed_L=do.call(rbind,fixed),fixed_meta=do.call(rbind,fm))
}
calculate_contrasts<-function(cohort,datasets,Q,U,formula,encoding,tag,dfcom,outdir,source_metadata=list(),existing_curves=NULL,workers=8L,model_refit=FALSE){
 for(sub in c('tables','private','quality'))dir.create(file.path(outdir,sub),recursive=TRUE,showWarnings=FALSE)
 M<-length(datasets);n<-nrow(datasets[[1]]);p<-nrow(Q)
 stopifnot(ncol(Q)==M,identical(dim(U),c(p,p,M)),all(vapply(datasets,nrow,integer(1))==n))
 sc<-expand.grid(pco2=seq(25,65,by=.1),rv=0:1);design<-contrast_design(sc)
 jobs<-lapply(seq_len(M),function(i)list(data=datasets[[i]],b=Q[,i],V=U[,,i],sc=sc,formula=formula,encoding=encoding))
 workerlog<-tempfile(paste0('rv_',tag,'_'),fileext='.log')
 cl<-parallel::makePSOCKcluster(min(workers,M),outfile=workerlog)
 parallel::clusterExport(cl,c('fixed_rcs','risk_one_completion'),envir=environment(risk_one_completion))
 parallel::clusterEvalQ(cl,{Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1);NULL})
 ans<-tryCatch(parallel::parLapplyLB(cl,jobs,risk_one_completion,chunk.size=1L),finally=parallel::stopCluster(cl))
 if(file.exists(workerlog))file.copy(workerlog,file.path(outdir,'quality',paste0(tag,'_workers.log')),overwrite=TRUE)
 RQ<-do.call(cbind,lapply(ans,'[[','q'));RU<-do.call(cbind,lapply(ans,'[[','u'))
 G<-simplify2array(lapply(ans,'[[','G'))
 risk<-cbind(cohort=cohort,sc,scalar_pool_table(RQ,RU,dfcom))
 sets<-list(relative40=list(L=design$relative_L,meta=design$relative_meta),sliding=list(L=design$sliding_L,meta=design$sliding_meta),fixed_contrasts=list(L=design$fixed_L,meta=design$fixed_meta))
 stored<-list();tables<-list(risks=risk)
 for(key in names(sets)){
  L<-sets[[key]]$L;CQ<-L%*%RQ;CU<-matrix(NA_real_,nrow(L),M)
  for(i in seq_len(M)){gg<-L%*%G[,,i];CU[,i]<-rowSums((gg%*%U[,,i])*gg)}
  tables[[key]]<-cbind(cohort=cohort,sets[[key]]$meta,scalar_pool_table(CQ,CU,dfcom))
  stored[[key]]<-list(Q=CQ,U_diag=CU,L=L,metadata=sets[[key]]$meta)
 }
 # The complete cross-grid covariance of I(x), including shared endpoints.
 ii<-which(design$sliding_meta$metric=='interaction');LI<-design$sliding_L[ii,,drop=FALSE]
 IQ<-LI%*%RQ;Ubar<-matrix(0,nrow(LI),nrow(LI))
 for(i in seq_len(M)){gg<-LI%*%G[,,i];Ubar<-Ubar+(gg%*%U[,,i]%*%t(gg))/M}
 B<-cov(t(IQ));T<-Ubar+(1+1/M)*B
 # Fixed contrasts retain their entire joint covariance as a separate family.
 LF<-design$fixed_L;FUbar<-matrix(0,nrow(LF),nrow(LF))
 for(i in seq_len(M)){gg<-LF%*%G[,,i];FUbar<-FUbar+(gg%*%U[,,i]%*%t(gg))/M}
 FB<-cov(t(stored$fixed_contrasts$Q));FT<-FUbar+(1+1/M)*FB
 frozen_error<-NA_real_
 if(!is.null(existing_curves)){
  keys<-function(z)paste(z$rv,sprintf('%.1f',z$pco2))
  z<-existing_curves[match(keys(risk),keys(existing_curves)),]
  stopifnot(!anyNA(z$risk))
  frozen_error<-max(abs(c(risk$estimate-z$risk,risk$lower-z$risk_lower,risk$upper-z$risk_upper)))
  stopifnot(frozen_error<1e-8)
 }
 slideI<-tables$sliding[ii,]
 stopifnot(max(abs(diag(T)-slideI$total_variance))<1e-10,max(abs(T-t(T)))<1e-12)
 zero<-tables$relative40[abs(tables$relative40$pco2-40)<1e-8,]
 stopifnot(all(abs(c(zero$estimate,zero$se,zero$lower,zero$upper))<1e-10))
 # I(40) is exactly the predeclared 40-to-50 difference in changes.
 ci<-which(abs(slideI$x_start-40)<1e-8);fi<-which(design$fixed_meta$start==40 & design$fixed_meta$end==50 & design$fixed_meta$metric=='interaction')
 stopifnot(max(abs(IQ[ci,]-stored$fixed_contrasts$Q[fi,]))<1e-10)
 for(key in names(tables))write.csv(tables[[key]],file.path(outdir,'tables',paste0(tag,'_',key,'.csv')),row.names=FALSE)
 y<-datasets[[1]]$y;rv<-datasets[[1]]$rv
 summary<-list(cohort=cohort,N=n,RV0=sum(rv==0),RV1=sum(rv==1),events=sum(y),M=M,parameters=p,dfcom=dfcom,formula=formula,
  exposure_encoding=encoding,exposure_knots=c(30,40,55),reference_PaCO2=40,slide_increment=10,risk_grid=c(25,65,.1),interaction_grid=c(25,55,.1),
  frozen_curve_max_absolute_error=frozen_error,interaction_T_diag_max_error=max(abs(diag(T)-slideI$total_variance)),
  interaction_T_max_asymmetry=max(abs(T-t(T))),no_new_imputation=TRUE,no_model_refit=!model_refit,
  variance_scope='Conditional on each completed empirical target covariate distribution; full coefficient covariance delta method, then Rubin pooling.',source_metadata=source_metadata)
 saveRDS(list(cohort=cohort,scenarios=sc,risk_Q=RQ,risk_U=RU,risk_gradients=G,coef_Q=Q,coef_V=U,dfcom=dfcom,
  contrasts=stored,interaction=list(metadata=design$sliding_meta[ii,],Q=IQ,L=LI,Ubar=Ubar,B=B,T=T),
  fixed_joint=list(metadata=design$fixed_meta,Q=stored$fixed_contrasts$Q,L=LF,Ubar=FUbar,B=FB,T=FT),summary=summary),
  file.path(outdir,'private',paste0(tag,'_inference.rds')),compress=FALSE)
 jsonlite::write_json(summary,file.path(outdir,'quality',paste0(tag,'_source_and_calculation_audit.json')),auto_unbox=TRUE,pretty=TRUE,digits=15)
 cat(cohort,'completed:',n,'patients,',M,'completions, frozen curve error',format(frozen_error,scientific=TRUE),'\n')
 invisible(summary)
}
