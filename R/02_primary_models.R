source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
O<-file.path(OUT,'primary')
make_dirs(O,c('tables','quality','fits'))
source(FUN('imputation.R'));source(FUN('pooling.R'))
con<-jsonlite::read_json(file.path(O,'quality/model_contract.json'),simplifyVector=TRUE)
raw<-read.csv(INPUT,check.names=FALSE)
prag_knots<-KN
inp <- prag_read_input(INPUT,KN)
imp <- lapply(seq_len(con$M),function(i)readRDS(file.path(O,'checkpoints',sprintf('imp_%03d.rds',i))))
ds <- lapply(imp,'[[','d')
stopifnot(all(vapply(imp,function(z)z$signature==con$signature && z$formula==con$formula && identical(z$ids,inp$sourceIDs),logical(1))))
fits <- lapply(ds,function(d){
  prag_validate(d,inp)
  glm(as.formula(con$formula),data=d,family=binomial(),control=glm.control(maxit=1000,epsilon=1e-9))
})
stopifnot(length(fits)==100L,all(vapply(fits,function(f)nobs(f)==2051L&&f$converged&&f$rank==length(coef(f))&&all(is.finite(vcov(f))),logical(1))),
  all(vapply(fits,function(f)!'hist'%in%names(coef(f)),logical(1))))
write_out <- function(x,label)write.csv(x,file.path(O,'tables',paste0(label,'.csv')),row.names=FALSE)
e <- pi_extract(fits);M <- e$m
I <- diag(e$p);dimnames(I)<-list(e$terms,e$terms)
H <- list(interaction=I[c('rv:x','rv:xnl'),,drop=FALSE],overall_PaCO2=I[c('x','xnl','rv:x','rv:xnl'),],
  rv0_overall=I[c('x','xnl'),],rv1_overall=I[c('x','xnl'),]+I[c('rv:x','rv:xnl'),],
  rv0_nonlinear=I['xnl',,drop=FALSE],rv1_nonlinear=I['xnl',,drop=FALSE]+I['rv:xnl',,drop=FALSE],
  any_nonlinearity=I[c('xnl','rv:xnl'),],nonlinear_interaction=I['rv:xnl',,drop=FALSE])
tests <- do.call(rbind,lapply(names(H),function(name){
  L<-H[[name]];U<-array(0,c(nrow(L),nrow(L),M))
  for(i in seq_len(M))U[,,i]<-L%*%e$U[,,i]%*%t(L)
  cbind(test=name,pi_d1_q_u(L%*%e$Q,U,e$dfcom))
}))
write_out(tests,'tests')
design <- function(d){d$xnl<-prag_rcs(d$x,c(30,40,55));d}
risks <- do.call(rbind,lapply(list(c(30,40),c(40,60)),function(xx){
  z<-pi_pool_standardized(fits,ds,design_fun=design,levels=xx)$table
  cbind(start=xx[1],end=xx[2],z)
}))
write_out(risks,'risks')
write_out(pi_or_curve(fits,ds,c(30,60),design_fun=design,reference=40),'withinOR')
coefs <- pi_pool_fits(fits)$table
write_out(coefs,'coefficients')
ratios <- do.call(rbind,lapply(list(c(30,40),c(40,60)),function(xx){
  L<-setNames(rep(0,e$p),e$terms)
  L['rv:x']<-diff(xx);L['rv:xnl']<-diff(prag_rcs(xx,c(30,40,55)))
  z<-pi_pool_contrast(fits,L)
  cbind(start=xx[1],end=xx[2],OR_ratio=exp(z$estimate),OR_lower=exp(z$lower),OR_upper=exp(z$upper),z)
}))
write_out(ratios,'interaction_OR_ratio')
trace <- abind::abind(lapply(imp,'[[','trace'),along=1)
draws <- aperm(trace[,,21:40,drop=FALSE],c(3,1,2));dimnames(draws)<-list(NULL,NULL,e$terms)
rh <- posterior::summarise_draws(posterior::as_draws_array(draws),rhat=posterior::rhat)
write_out(rh,'Rhat')
qa <- list(N=2051,deaths=998,M=M,formula=con$formula,parameters=e$p,
  history_flag_removed_from_outcome_model=!'hist'%in%e$terms,
  history_retained_as_imputation_auxiliary=TRUE,all_fits_converged_full_rank=TRUE,
  all_observed_values_validated=TRUE,max_Rhat=max(rh$rhat),
  max_interaction_Rhat=max(rh$rhat[rh$variable%in%c('rv:x','rv:xnl')]),
  max_interaction_MCSE_SE=max(with(coefs[coefs$term%in%c('rv:x','rv:xnl'),],mi_mcse/std.error)),
  max_delta_MCSE_SE=max(with(risks[risks$term=='delta',],mi_mcse/std.error)),
  imputation_warning_count=sum(vapply(imp,function(z)length(z$warnings),integer(1))),
  warnings=unique(unlist(lapply(imp,'[[','warnings'))),
  adjudicated_component_initialization_chains=sum(vapply(imp,'[[',logical(1),'component_initialization_warnings_adjudicated')))
jsonlite::write_json(qa,file.path(O,'quality/nohist_QA.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
saveRDS(list(fits=fits,data=ds,formula=con$formula,sourceIDs=inp$sourceIDs,signature=con$signature),file.path(O,'fits/model.rds'))
print(tests[,c('test','p.value')],row.names=FALSE)
print(risks[risks$term=='delta',c('start','end','estimate','lower','upper')],row.names=FALSE)
