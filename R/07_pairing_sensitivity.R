source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
O<-file.path(OUT,'sensitivity_6h')
make_dirs(O,c('tables','quality','private'))
source(FUN('pooling.R'));source(FUN('risk_contrasts.R'))
rawfile<-INPUT
fitfile<-file.path(OUT,'primary/fits/model.rds')
raw<-read.csv(rawfile,check.names=FALSE)
frozen<-readRDS(fitfile)
stopifnot(nrow(raw)==2051L,length(frozen$data)==100L,length(frozen$fits)==100L,
 identical(as.character(frozen$sourceIDs),as.character(raw$subject_id)),!anyDuplicated(raw$subject_id))
stamp<-function(x)as.POSIXct(x,format='%Y-%m-%d %H:%M:%S',tz='UTC')
gap<-abs(as.numeric(difftime(stamp(raw$bg_time),stamp(raw$echo_time),units='hours')))
keep<-which(gap<=6)
stopifnot(length(keep)==1812L,sum(raw$death28_echo[keep])==885L,
 sum(raw$rv_binary[keep]==0)==1102L,sum(raw$rv_binary[keep]==1)==710L,
 all(stamp(raw$t0[keep])==stamp(raw$echo_time[keep])))
for(d in frozen$data)stopifnot(nrow(d)==2051L,all(d$y==raw$death28_echo),
 all(d$rv==raw$rv_binary),max(abs(d$x-raw$pco2))<1e-10)
form<-formula(frozen$fits[[1]])
formstr<-paste(deparse(form,width.cutoff=500L),collapse=' ')
stopifnot(!'hist'%in%all.vars(form),length(coef(frozen$fits[[1]]))==25L)
datasets<-lapply(frozen$data,function(d)d[keep,,drop=FALSE])
warnings<-vector('list',100L)
checks<-list()
fits<-lapply(seq_len(100),function(i){
 w<-character()
 f<-withCallingHandlers(glm(form,family=binomial(),data=datasets[[i]],control=glm.control(maxit=500)),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')})
 stopifnot(f$converged,nobs(f)==1812L,f$rank==25L,all(is.finite(coef(f))),all(is.finite(vcov(f))))
 if(i%in%c(1L,100L)){
   sep<-detectseparation::detect_separation(model.matrix(f),f$y,family=binomial(),intercept=TRUE)
   stopifnot(!any(is.infinite(sep$coefficients)))
 }
 warnings[[i]]<<-w
 checks[[i]]<<-data.frame(imputation=i,n=nobs(f),events=sum(f$y),parameters=length(coef(f)),iterations=f$iter,converged=f$converged,min_fitted=min(fitted(f)),max_fitted=max(fitted(f)),warnings=length(w))
 f
})
ex<-pi_extract(fits)
nm<-rownames(ex$Q)
unit<-function(term){a<-numeric(ex$p);a[match(term,nm)]<-1;a}
CT<-rbind(unit('x'),unit('xnl'));IT<-rbind(unit('rv:x'),unit('rv:xnl'))
contrasts<-list(interaction=IT,rv0_overall=CT,rv1_overall=CT+IT,
 overall_association=rbind(CT,IT),rv0_nonlinear=CT[2,,drop=FALSE],
 rv1_nonlinear=(CT+IT)[2,,drop=FALSE],nonlinear_interaction=IT[2,,drop=FALSE],
 any_nonlinearity=rbind(CT[2,,drop=FALSE],IT[2,,drop=FALSE]))
independent<-numeric()
tests<-do.call(rbind,lapply(names(contrasts),function(name){
 H<-contrasts[[name]];q<-H%*%ex$Q;u<-array(0,c(nrow(H),nrow(H),100L))
 for(i in seq_len(100L))u[,,i]<-H%*%ex$U[,,i]%*%t(H)
 p<-pi_d1_q_u(q,u,ex$dfcom)
 rownames(q)<-paste0('C',seq_len(nrow(H)))
 dimnames(u)<-list(rownames(q),rownames(q),NULL)
 qlist<-lapply(seq_len(100L),function(i)setNames(as.numeric(q[,i]),rownames(q)))
 ulist<-lapply(seq_len(100L),function(i)matrix(u[,,i],nrow=nrow(H),dimnames=list(rownames(q),rownames(q))))
 alt<-mitml::testConstraints(qhat=qlist,uhat=ulist,constraints=rownames(q),df.com=ex$dfcom)
 independent<<-c(independent,abs(p$p.value-alt$test[1,4]))
 data.frame(test=name,df=p$df1,p_value=p$p.value,statistic=p$statistic,df_den=p$df2)
}))
stopifnot(max(independent)<1e-8)
write.csv(tests,file.path(O,'tables/mimic6h_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(O,'quality/mimic6h_refit_checks.csv'),row.names=FALSE)
saveRDS(list(data=datasets,Q=ex$Q,U=ex$U,dfcom=ex$dfcom,formula=formstr,
 sourceIDs=frozen$sourceIDs[keep],main_source_positions=keep),file.path(O,'private/mimic6h_refits.rds'))
metadata<-list(analysis='fixed-original-pair bilateral <=6-hour sensitivity',main_N=2051L,
 N=1812L,deaths=885L,rv0=1102L,rv1=710L,M=100L,parameters=25L,
 source_ID_order_exact=TRUE,outcome='death28_echo',all_t0_equals_echo=TRUE,
 same_original_echo_gas_pair=TRUE,no_new_gas_selected=TRUE,
 original_IMV_episode_and_pressor_course_eligibility_inherited=TRUE,
 covariance_policy='current no-hist formula and fixed knots; refit after subset',
 no_new_imputation=TRUE,completed_data_source='current main model 100 completed datasets',
 source_fit_sha256=digest::digest(file=fitfile,algo='sha256'),
 source_data_sha256=digest::digest(file=rawfile,algo='sha256'),
 fit_warnings=sum(lengths(warnings)),independent_D1_checks=length(independent),
 independent_D1_max_difference=max(independent))
jsonlite::write_json(metadata,file.path(O,'quality/mimic6h_refit_audit.json'),auto_unbox=TRUE,pretty=TRUE,digits=15)
calculate_contrasts(cohort='MIMIC-IV bilateral 6 h',datasets=datasets,Q=ex$Q,U=ex$U,
 formula=formstr,encoding='mimic_raw',tag='mimic6h',dfcom=ex$dfcom,outdir=O,
 source_metadata=metadata,workers=WORKERS,model_refit=TRUE)
writeLines(capture.output(sessionInfo()),file.path(O,'quality/sessionInfo.txt'))
print(tests)
