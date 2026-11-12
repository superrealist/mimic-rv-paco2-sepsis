# Exact rejection-envelope normalization for the unchanged logistic SMC target.
# The installed smcfcs source is not modified; the wrapper preserves the RNG sequence.
vent_rcs_scalar<-function(x,k){
 (pmax(x-k[1],0)^3-(k[3]-k[1])/(k[3]-k[2])*pmax(x-k[2],0)^3+(k[2]-k[1])/(k[3]-k[2])*pmax(x-k[3],0)^3)/(k[3]-k[1])^2
}
vent_bmi_extrema<-function(lo,hi,b_linear,b_nonlinear,knots){
 k<-as.numeric(knots);stopifnot(length(k)==3L,all(diff(k)>0),length(lo)==length(hi),all(is.finite(c(lo,hi,b_linear,b_nonlinear))),all(lo<=hi))
 D<-k[3]-k[1];A<-(k[3]-k[1])/(k[3]-k[2]);roots<-numeric()
 # On each inner knot interval the derivative is a quadratic polynomial.
 # Below k1 and above k3 it is constant, so endpoints suffice there.
 for(seg in 1:2){
  w<-if(seg==1)1 else c(1,-A);kk<-k[seq_along(w)]
  pol<-c(b_linear+3*b_nonlinear/D^2*sum(w*kk^2),-6*b_nonlinear/D^2*sum(w*kk),3*b_nonlinear/D^2*sum(w))
  if(any(pol!=0)){
   zz<-polyroot(pol);rr<-Re(zz[abs(Im(zz))<=1e-7*pmax(1,abs(Re(zz)))])
   rr<-rr[rr>=k[seg]-1e-8&rr<=k[seg+1]+1e-8];roots<-c(roots,pmin(k[seg+1],pmax(k[seg],rr)))
  }
 }
 candidates<-sort(unique(c(k,roots)));f<-function(x)b_linear*x+b_nonlinear*vent_rcs_scalar(x,k)
 v0<-f(lo);v1<-f(hi);mn<-pmin(v0,v1);mx<-pmax(v0,v1)
 for(ca in candidates){ix<-lo<=ca&ca<=hi;if(any(ix)){v<-f(ca);mn[ix]<-pmin(mn[ix],v);mx[ix]<-pmax(mx[ix],v)}}
 list(minimum=mn,maximum=mx,candidates=candidates,lower=lo,upper=hi)
}

prag_impute_exact_envelope<-function(inp,m=1L,numit=40L,seed,rjlimit=100000L,maxit=1000L){
 ns<-asNamespace('smcfcs');hi<-new.env(parent=ns);source_ids<-inp$sourceIDs;bmi_k<-as.numeric(inp$knots$bmi);pf_k<-as.numeric(inp$knots$pf)
 hi$glm<-local({cap<-as.integer(maxit);function(...,control=list()){
  control$maxit<-max(cap,control$maxit)
  z<-stats::glm(...,control=control)
  if(!isTRUE(z$converged)&&inherits(z,'brglmFit')){
   rng_before<-get('.Random.seed',envir=.GlobalEnv)
   args<-list(...);args$control<-control;args$control$maxit<-10000L;args$control$slowit<-.25
   args$start<-coef(z);z<-do.call(stats::glm,args)
   stopifnot(identical(rng_before,get('.Random.seed',envir=.GlobalEnv)))
   previous<-audit$bias_reduced_numeric_retries;if(is.null(previous))previous<-0L
   audit$bias_reduced_numeric_retries<-previous+1L
  }
  stopifnot(isTRUE(z$converged),all(is.finite(coef(z))))
  if(inherits(z,'brglmFit') && ncol(model.matrix(z))==1L && length(unique(z$y))==1L){
   expected<-qlogis((sum(z$y)+.5)/(length(z$y)+1));err<-max(abs(coef(z)-expected))
   stopifnot(err<1e-6,all(is.finite(vcov(z))),all(diag(vcov(z))>0))
   oldn<-audit$constant_binary_intercept_checks;if(is.null(oldn))oldn<-0L
   audit$constant_binary_intercept_checks<-oldn+1L
   olderr<-audit$constant_binary_max_analytic_error;if(is.null(olderr))olderr<-0
   audit$constant_binary_max_analytic_error<-max(olderr,err)
  }
  audit$all_returned_glms_converged<-TRUE
  z
 }})
 audit<-new.env(parent=emptyenv());audit$counts<-list();audit$inactive_counts<-list();audit$active_bmi_counts<-list();audit$active_pf_counts<-list();audit$max_inactive_matrix_error<-0;audit$max_offset_difference<-0;audit$max_scaled_probability<-0;audit$minimum_log_envelope<-0;audit$callbacks<-0L;audit$failures<-list()
 hi$.vent_exact_prob<-function(prob,fr,batch){
  getf<-function(n)get(n,envir=fr,inherits=FALSE)
  target<-names(getf('originaldata'))[getf('targetCol')]
  if(!target%in%c('hlat','wlat','zfio','zpeep','zmv'))return(prob)
  if(!getf('method')[getf('targetCol')]%in%c('norm','latnorm'))return(prob)
  rng0<-get('.Random.seed',envir=.GlobalEnv)
  row_ids<-if(batch)getf('i') else getf('imputationNeeded')
  d<-if(batch)getf('tempData') else getf('imputations')[[getf('imp')]][row_ids,,drop=FALSE]
  eta<-as.numeric(if(batch)getf('outmodxb') else getf('outmodxb')[row_ids]);y<-as.numeric(d[[getf('outcomeCol')]])
  stopifnot(length(prob)==nrow(d),length(eta)==nrow(d),all(is.finite(eta)),all(y%in%c(0,1)))
  yes<-function(n)if(n%in%names(d))!is.na(d[[n]])&d[[n]]==1 else rep(FALSE,nrow(d))
  inactive<-switch(target,hlat=yes('bmiknown'),wlat=yes('bmiknown'),zfio=yes('pfknown')|yes('plow')|yes('phigh'),zpeep=yes('peepknown')|yes('pzero'),zmv=yes('mvknown'))
  logL<-plogis((2*y-1)*eta,log.p=TRUE);logM<-rep(0,length(logL))
  # Algebraic inactivity is also verified by finite, broad latent perturbation.
  if(any(inactive)){
   check<-d[which(inactive),,drop=FALSE];if(batch)check<-check[1,,drop=FALSE]
   du<-dv<-check;du[[target]]<- -4;dv[[target]]<-4
   update<-get('updatePassiveVars',envir=ns)
   du<-update(du,getf('method'),getf('passiveVars'));dv<-update(dv,getf('method'),getf('passiveVars'))
   Xu<-model.matrix(as.formula(getf('smformula')),du);Xv<-model.matrix(as.formula(getf('smformula')),dv)
   dif<-max(abs(Xu-Xv));stopifnot(is.finite(dif),dif==0)
   audit$max_inactive_matrix_error<-max(audit$max_inactive_matrix_error,dif)
   logM[inactive]<-logL[inactive]
  }
  active_bmi<-!inactive&target%in%c('hlat','wlat')
  active_pf<-!inactive&target=='zfio';active_spline<-active_bmi|active_pf
  if(any(active_spline)){
   # A target-dependent BMI/PF is the only substantive dependency of these
   # bounded coordinates; interactions or direct height/weight/PaO2 terms would
   # invalidate this envelope and are rejected explicitly.
   check<-if(batch)d[1,,drop=FALSE] else d
   X<-model.matrix(as.formula(getf('smformula')),check);b<-as.numeric(getf('outcomeModBeta'));names(b)<-colnames(X)
   pair<-if(target=='zfio')c('pf','pfnl') else c('bmi','bminl');kw<-if(target=='zfio')pf_k else bmi_k
   stopifnot(all(pair%in%names(b)),!any(grepl(pair[1],names(b),fixed=TRUE)&!names(b)%in%pair),
    !any(c('height','weight','hlat','wlat','zfio','po2','flim')%in%all.vars(as.formula(getf('smformula')))))
   X0<-X;X0[,pair]<-0;offset<-drop(X0%*%b)
   if(target=='hlat'){lo<-check$weight/6.25;up<-check$weight}else if(target=='wlat'){lo<-20/(check$height/100)^2;up<-400/(check$height/100)^2}else{stopifnot(all(check$flim>=.21&check$flim<=1));lo<-check$po2;up<-check$po2/check$flim}
   ext<-vent_bmi_extrema(lo,up,b[[pair[1]]],b[[pair[2]]],kw)
   ymax<-offset+ifelse(check[[getf('outcomeCol')]]==1,ext$maximum,ext$minimum)
   lm<-plogis((2*check[[getf('outcomeCol')]]-1)*ymax,log.p=TRUE)
   # Conservative positive padding may only lower efficiency, not alter target.
   lm<-pmin(0,lm+1e-10*(1+abs(ymax)))
   if(batch){stopifnot(length(lm)==1L);logM[active_spline]<-lm}else logM[active_spline]<-lm[active_spline]
   # Independent zeroing of the two spline columns makes offset independent.
   dd<-check;dd[[target]]<-0
   dd<-get('updatePassiveVars',envir=ns)(dd,getf('method'),getf('passiveVars'))
   XX<-model.matrix(as.formula(getf('smformula')),dd);XX[,pair]<-0
   dif<-max(abs(drop(XX%*%b)-offset));stopifnot(is.finite(dif),dif==0)
   audit$max_offset_difference<-max(audit$max_offset_difference,dif)
  }
  p<-exp(logL-logM);stopifnot(all(is.finite(p)),all(p>=0),max(p)<=1+1e-8)
  p<-pmin(1,p);out<-prob;out[]<-p
  bump<-function(slot,n){old<-audit[[slot]][[target]];if(is.null(old))old<-0;audit[[slot]][[target]]<-old+n}
  bump('counts',length(p));bump('inactive_counts',sum(inactive));bump('active_bmi_counts',sum(active_bmi));bump('active_pf_counts',sum(active_pf))
  audit$callbacks<-audit$callbacks+1L;audit$max_scaled_probability<-max(audit$max_scaled_probability,max(p));audit$minimum_log_envelope<-min(audit$minimum_log_envelope,min(logM))
  stopifnot(identical(rng0,get('.Random.seed',envir=.GlobalEnv)))
  out
 }
 hi$.vent_record_failure<-function(fr){
  val<-function(n)if(exists(n,envir=fr,inherits=FALSE))get(n,envir=fr,inherits=FALSE) else NULL
  t<-val('targetCol');i<-val('i');pr<-as.numeric(val('prob'));lp<-as.numeric(val('outmodxb'))
  audit$failures[[length(audit$failures)+1L]]<-list(variable=names(val('originaldata'))[t],iteration=val('cyclenum'),row=i,source_subject_id=source_ids[i],rjlimit=val('rjlimit'),probability_range=range(pr),linear_predictor_range=range(lp),nonfinite_probability=sum(!is.finite(pr)))
  invisible(NULL)
 }
 core<-get('smcfcs.core',envir=ns);old<-paste(deparse(core,width.cutoff=500L),collapse='\n');txt<-old;patches<-list()
 patch<-function(needle,replacement){
  hit<-gregexpr(needle,txt,fixed=TRUE)[[1]];stopifnot(length(hit)==1L,hit[1]>0)
  txt<<-sub(needle,replacement,txt,fixed=TRUE);patches[[length(patches)+1L]]<<-list(original=needle,replacement=replacement)
 }
 a<-'prob <- prob * imputations[[imp]][imputationNeeded, outcomeCol] + (1 - prob) * (1 - imputations[[imp]][imputationNeeded, outcomeCol])'
 patch(a,paste0(a,'\n                        prob <- .vent_exact_prob(prob, environment(), FALSE)'))
 a<-'prob <- prob * tempData[, outcomeCol] + (1 - prob) * (1 - tempData[, outcomeCol])'
 patch(a,paste0(a,'\n                        prob <- .vent_exact_prob(prob, environment(), TRUE)'))
 a<-'rjFailCount <- rjFailCount + 1';patch(a,paste0('{ ',a,'; .vent_record_failure(environment()) }'))
 inv<-txt;for(p in rev(patches))inv<-sub(p$replacement,p$original,inv,fixed=TRUE);stopifnot(identical(inv,old))
 hi$smcfcs.core<-eval(parse(text=txt),envir=hi)
 env<-new.env(parent=.GlobalEnv);env$prag_smcfcs_exact<-get('smcfcs',envir=ns);environment(env$prag_smcfcs_exact)<-hi
 runner<-eval(parse(text=gsub('smcfcs::smcfcs','prag_smcfcs_exact',paste(deparse(prag_impute),collapse='\n'),fixed=TRUE)),envir=env)
 ans<-runner(inp,m=m,numit=numit,seed=seed,rjlimit=rjlimit)
 ans$rejection_failure_context<-audit$failures;ans$rejection_events<-audit$failures
 ans$exact_envelope_diagnostics<-as.list(audit)
 ans$metadata$numerical_glm_maxit<-as.integer(maxit);ans$metadata$rjlimit<-as.integer(rjlimit)
 ans$metadata$exact_envelope<-list(method='Exact global BMI/PF-likelihood envelope and inactive-coordinate factor cancellation; normal proposal q and substantive logistic formula unchanged.',installed_core_sha256=digest::digest(old,algo='sha256',serialize=FALSE),patched_core_sha256=digest::digest(txt,algo='sha256',serialize=FALSE),patches=patches,inverse_patch_identity_verified=TRUE,no_extra_RNG_calls=TRUE,BMI_domain='hlat: weight/6.25..weight; wlat:20/(height/100)^2..400/(height/100)^2',PF_domain='zfio with pfknown=0,plow=0,phigh=0: po2..po2/flim; flim=max(.21,po2/713)',global_candidates='Both row-domain endpoints, all contained fixed knots, real roots of the two interior quadratic derivative polynomials; outer derivative is constant.',scaled_probability_bound_assertion=TRUE)
 ans
}
prag_impute_envelope<-prag_impute_exact_envelope
