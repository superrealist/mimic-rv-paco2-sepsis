# Preparation helpers only: sourcing this file neither imputes nor fits a model.
# Outcome-enriched FCS is the approximate main MI strategy, not exact SMC-FCS.
# Requires imputation.R to have been sourced in this R process.
VFD_MI_DIAGNOSTICS<-new.env(parent=emptyenv())
VFD_MI_DIAGNOSTICS$logreg_calls<-0L;VFD_MI_DIAGNOSTICS$constant_binary_calls<-0L
VFD_MI_DIAGNOSTICS$max_glm_iterations<-0L;VFD_MI_DIAGNOSTICS$warnings<-character()

mice.impute.vent_binary<-function(y,ry,x,wy=NULL,...) {
 if(is.null(wy))wy<-!ry
 obs<-as.numeric(y[ry]);stopifnot(length(obs)>0L,all(obs%in%c(0,1)))
 if(length(unique(obs))==1L) {
  # A constant observed endpoint cannot support a covariate regression. The
  # prespecified intercept-only Jeffreys posterior retains probability uncertainty.
  prob<-stats::rbeta(1L,sum(obs==1)+.5,sum(obs==0)+.5)
  VFD_MI_DIAGNOSTICS$constant_binary_calls<-VFD_MI_DIAGNOSTICS$constant_binary_calls+1L
  return(stats::rbinom(sum(wy),1L,prob))
 }
 # Preserve mice's data augmentation and coefficient/draw algorithm; only
 # raise the iterative GLM cap and expose numerical diagnostics. No RNG call
 # is added by the cap wrapper itself.
 hi<-new.env(parent=asNamespace('mice'))
 hi$glm.fit<-function(...,control=list()){
  control$maxit<-max(1000L,control$maxit)
  f<-withCallingHandlers(stats::glm.fit(...,control=control),warning=function(w){
   VFD_MI_DIAGNOSTICS$warnings<-c(VFD_MI_DIAGNOSTICS$warnings,conditionMessage(w))
  })
  VFD_MI_DIAGNOSTICS$logreg_calls<-VFD_MI_DIAGNOSTICS$logreg_calls+1L
  VFD_MI_DIAGNOSTICS$max_glm_iterations<-max(VFD_MI_DIAGNOSTICS$max_glm_iterations,f$iter)
  if(!isTRUE(f$converged))stop('VFD binary component GLM did not converge at cap1000; numerical adjudication required.')
  f
 }
 fun<-get('mice.impute.logreg',envir=asNamespace('mice'));environment(fun)<-hi
 fun(y=y,ry=ry,x=x,wy=wy,...)
}

prepare_vfd_mice<-function(inp,exposure=c('binary','graded'),adjust_dose=TRUE) {
 exposure<-match.arg(exposure)
 stopifnot(exists('prag_prepare',mode='function'),exists('prag_rcs',mode='function'))
 stopifnot(is.data.frame(inp$data),is.data.frame(inp$source),nrow(inp$data)==nrow(inp$source))
 stopifnot(all(is.finite(inp$data$y)),all(inp$data$y>=0&inp$data$y<=28),all(inp$source$death28%in%c(0,1)))
 stopifnot(all(inp$data$y[inp$source$death28==1]==0))
 # The caller must have set y to the new observed VFD28 before this function.
 p<-prag_prepare(inp);d<-p$data;meth<-p$methods_named
 d$sofa_imp<-as.numeric(d$sofa)
 d$sofa_ord<-NULL;meth<-meth[names(d)[names(d)%in%names(meth)]]
 meth['sofa_imp']<-if(anyNA(d$sofa_imp))'pmm' else ''
 d$death28_aux<-as.integer(inp$source$death28)
 d$ynl<-prag_rcs(d$y,c(0,14,28))
 stopifnot(all(is.finite(d$ynl)))
 extras<-'lcount'
 stopifnot('lcount'%in%names(d),all(is.finite(d$lcount)))
 if(adjust_dose){extras<-c(extras,'lnee');stopifnot('lnee'%in%names(d),all(is.finite(d$lnee)))}
 bases<-if(exposure=='binary')'rv' else c('rv_g_mild','rv_g_modsev','rv_g_ungraded')
 stopifnot(all(bases%in%names(d)),all(vapply(d[bases],function(z)all(z%in%c(0,1)),logical(1))))
 if(exposure=='graded')stopifnot(all(rowSums(d[bases])<=1L))
 fullobs<-c('x','xnl','age','agenl','sex','copd','chf','ph','cad','ckd','hist','lvent',extras,bases,'y','ynl','death28_aux')
 products<-character()
 for(g in bases)for(v in c('x','xnl','y','ynl')){
  nm<-paste0(g,'__',v);d[[nm]]<-d[[g]]*d[[v]];products<-c(products,nm)
 }
 fullobs<-unique(c(fullobs,products))
 stopifnot(!anyNA(d[fullobs]),all(is.finite(as.matrix(d[fullobs]))))
 primitive<-c('hlat','wlat','lmap','llac','plow','phigh','zfio','pzero','zpeep','zmv','sofa_imp')
 passive<-c(height='100+150*plogis(hlat)',weight='20+380*plogis(wlat)',
  bmi='ifelse(bmiknown==1,bmi0,weight/(height/100)^2)',map='exp(lmap)',
  pf='ifelse(pfknown==1,pf0,ifelse(plow==1,po2/.21,ifelse(phigh==1,po2,po2/(flim+(1-flim)*plogis(zfio)))))',
  peep='ifelse(peepknown==1,peep0,ifelse(pzero==1,0,40*plogis(zpeep)))',
  lmv='ifelse(mvknown==1,lmv0,log2(50*plogis(zmv)))',sofa='sofa_imp',
  bminl='prag_nl_bmi(bmi)',mapnl='prag_nl_map(map)',pfnl='prag_nl_pf(pf)')
 # All listed variables stay in the data so passive expressions retain source
 # flags. Only whitelisted primitive or fully observed columns predict targets.
 method<-setNames(rep('',ncol(d)),names(d))
 for(v in primitive)if(anyNA(d[[v]]))method[v]<-if(v%in%c('plow','phigh','pzero'))'vent_binary' else 'pmm'
 for(v in names(passive))if(anyNA(d[[v]]))method[v]<-paste0('~I(',passive[[v]],')')
 pred<-matrix(0L,ncol(d),ncol(d),dimnames=list(names(d),names(d)))
 active<-names(method)[method%in%c('pmm','vent_binary')]
 for(v in active)pred[v,setdiff(unique(c(fullobs,primitive)),v)]<-1L
 # Keep the native mixture-coordinate independence restriction. The outcome
 # predictors enrich FCS directly; there is no substantive rejection step.
 pred[c('plow','phigh','zfio'),c('plow','phigh','zfio')]<-0L
 pred[c('pzero','zpeep'),c('pzero','zpeep')]<-0L
 # Rare endpoint binary regressions use a compact, fixed clinical set.
 compact<-unique(c(bases,'x','age','sex','lvent','y','death28_aux',paste0(bases,'__y')))
 for(v in intersect(active,c('plow','pzero'))){pred[v,]<-0L;pred[v,compact]<-1L}
 if(exposure=='graded')pred[,'rv']<-0L
 constant<-names(d)[vapply(d,function(z)!anyNA(z)&&length(unique(z))<2L,logical(1))]
 pred[,constant]<-0L
 # With only primitive predictors, passive variables can be updated once at
 # the end of each cycle without entering any target with stale derivatives.
 visit<-c(intersect(primitive,active),names(passive)[grepl('^~',method[names(passive)])])
 stopifnot(!anyDuplicated(visit),all(method[c('y','ynl','death28_aux',bases)]==''),all(rowSums(pred[setdiff(names(d),active),,drop=FALSE])==0L))
 for(v in active[method[active]=='pmm'])stopifnot(sum(!is.na(d[[v]]))>=5L)
 list(data=d,method=method,predictorMatrix=pred,visitSequence=visit,
  where=is.na(d),donors=5L,matchtype=1L,ridge=1e-5,
  mice_controls=list(remove.constant=FALSE,remove.collinear=FALSE),
  policy=list(outcome='Observed VFD28, never imputed; all 28d decedents retain 0',
   exposure=exposure,observed_values_unchanged=TRUE,ordinary_covariates_only=TRUE,
   label='Main outcome-enriched approximate FCS; not exact substantive-model-compatible MI. Complete-case models are selection sensitivities.',
   outcome_features='VFD RCS knots0/14/28, death28 component, RV-by-VFD and RV-by-PaCO2 bases',
   bounded_components='Frozen native height/weight, FiO2 mixture, PEEP mixture, MV coordinates; donor PMM for continuous primitives and numeric integer SOFA',
   rare_binary='mice augmented logreg with glm maxit1000 and convergence audit; if observed binary constant, prespecified Jeffreys Beta(0.5,0.5) intercept posterior',
   controls='Do not let mice silently delete predictors; archive loggedEvents, warnings, full predictor matrix and package versions',
   completeness='Source outcome/pressor/IMV states and unknown NEE are excluded from imputation',
   validation='After completing, apply frozen physical/observed validator plus y/death/exposure invariance and derived-equation checks',
   limitation='Predictor enrichment does not establish exact congeniality, MAR, or correct conditional distributions.'))
}
