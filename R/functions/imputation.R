# Native SMC-FCS with bounded component/mixture parametrizations.
prag_rcs <- function(x,k) {
  (pmax(x-k[1],0)^3-(k[3]-k[1])/(k[3]-k[2])*pmax(x-k[2],0)^3+
     (k[2]-k[1])/(k[3]-k[2])*pmax(x-k[3],0)^3)/(k[3]-k[1])^2
}
prag_knots <- NULL
prag_nl_age<-function(x) prag_rcs(x,prag_knots$age)
prag_nl_bmi<-function(x) prag_rcs(x,prag_knots$bmi)
prag_nl_map<-function(x) prag_rcs(x,prag_knots$map)
prag_nl_pf<-function(x) prag_rcs(x,prag_knots$pf)
prag_formula<-paste('y~rv+x+xnl+rv:x+rv:xnl+age+agenl+bmi+bminl+map+mapnl+pf+pfnl+',
  'llac+lmv+sofa+peep+lvent+sex+copd+chf+ph+cad+ckd')
prag_generic<-c('y','rv','x','age','bmi','map','pf','llac','lmv','sofa','peep','lvent',
  'sex','copd','chf','ph','cad','ckd','hist','sdoc')

prag_read_input<-function(path,knots) {
  s<-read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(s$subject_id),all(s$death28%in%c(0,1)))
  d<-data.frame(y=s$death28,rv=s$rv_binary,x=s$pco2,age=s$age,bmi=s$bmi,map=s$map,pf=s$pf,
    llac=log2(s$lactate),lmv=log2(s$minute_volume),sofa=s$noncardiopulmonary_sofa,
    peep=s$peep,lvent=log1p(s$vent_hours),sex=as.integer(s$sex=='M'),
    copd=s$copd,chf=s$chf,ph=s$pulmonary_hypertension,cad=s$cad,ckd=s$ckd,
    hist=s$prior_history_observed,sdoc=s$s_doc)
  stopifnot(!anyNA(d[c('y','rv','x','age','lvent','sex','copd','chf','ph','cad','ckd','hist','sdoc')]))
  list(data=d,source=s,sourceIDs=s$subject_id,knots=knots,path=path)
}

prag_prepare<-function(inp) {
  d<-inp$data;s<-inp$source;prag_knots<<-inp$knots
  # Component bounds match the predefined valid source domains; all observed
  # component extrema lie strictly inside these domains, including extreme BMI.
  stopifnot(all(s$height>100&s$height<250,na.rm=TRUE),
    all(s$weight>20&s$weight<400,na.rm=TRUE),all(s$po2>0&s$po2<713))
  d$hlat<-qlogis((s$height-100)/150)
  d$wlat<-qlogis((s$weight-20)/380)
  d$height<-s$height;d$weight<-s$weight
  d$bmiknown<-as.integer(!is.na(d$bmi));d$bmi0<-ifelse(d$bmiknown==1,d$bmi,0)
  d$lmap<-log(d$map)
  # PF model: observed PF is never altered. At missing PF, FiO2 is a bounded
  # three-state mixture (room-air endpoint, FiO2=1 endpoint, interior).
  # A flagged measured FiO2 is not trusted or 'repaired'; its analysis PF is
  # missing and its imputation remains explicitly model-based.
  d$po2<-s$po2;d$flim<-pmax(.21,d$po2/713)
  d$pfknown<-as.integer(!is.na(d$pf));d$pf0<-ifelse(d$pfknown==1,d$pf,0)
  fi<-d$po2/d$pf
  low<-!is.na(fi)&abs(fi-.21)<1e-9
  high<-!is.na(fi)&abs(fi-1)<1e-9
  interior<-!is.na(fi)&!low&!high
  stopifnot(all(fi[interior]>d$flim[interior]),all(fi[interior]<1))
  d$plow<-ifelse(d$pfknown==1,as.integer(low),NA_integer_)
  # At PaO2>713*.21 the room-air endpoint is logically impossible in the
  # quality-screened analysis definition, so the low state is known zero.
  d$plow[d$po2>713*.21]<-0L
  d$phigh<-ifelse(d$pfknown==1&!low,as.integer(high),NA_integer_)
  d$zfio<-NA_real_
  d$zfio[interior]<-qlogis((fi[interior]-d$flim[interior])/(1-d$flim[interior]))
  d$peepknown<-as.integer(!is.na(d$peep));d$peep0<-ifelse(d$peepknown==1,d$peep,0)
  d$pzero<-ifelse(is.na(d$peep),NA_integer_,as.integer(d$peep==0))
  d$zpeep<-ifelse(d$peep>0,qlogis(d$peep/40),NA_real_)
  d$mvknown<-as.integer(!is.na(d$lmv));d$lmv0<-ifelse(d$mvknown==1,d$lmv,0)
  d$zmv<-qlogis(s$minute_volume/50)
  d$sofa_ord<-ordered(d$sofa,levels=sort(unique(na.omit(d$sofa))))
  d$xnl<-prag_rcs(d$x,c(30,40,55));d$agenl<-prag_nl_age(d$age)
  d$bminl<-prag_nl_bmi(d$bmi);d$mapnl<-prag_nl_map(d$map);d$pfnl<-prag_nl_pf(d$pf)
  meth<-setNames(rep('',ncol(d)),names(d))
  for(v in c('hlat','wlat','lmap','llac','zfio','zpeep','zmv')) if(anyNA(d[[v]]))meth[v]<-'norm'
  for(v in c('plow','phigh','pzero')) if(anyNA(d[[v]]))meth[v]<-'brlogreg'
  if(anyNA(d$sofa_ord))meth['sofa_ord']<-'podds'
  passive<-c(height='100+150*plogis(hlat)',weight='20+380*plogis(wlat)',
    bmi='ifelse(bmiknown==1,bmi0,weight/(height/100)^2)',map='exp(lmap)',
    pf='ifelse(pfknown==1,pf0,ifelse(plow==1,po2/.21,ifelse(phigh==1,po2,po2/(flim+(1-flim)*plogis(zfio)))))',
    peep='ifelse(peepknown==1,peep0,ifelse(pzero==1,0,40*plogis(zpeep)))',
    lmv='ifelse(mvknown==1,lmv0,log2(50*plogis(zmv)))',
    sofa='as.numeric(as.character(sofa_ord))',
    bminl='prag_nl_bmi(bmi)',mapnl='prag_nl_map(map)',pfnl='prag_nl_pf(pf)')
  for(v in names(passive))if(anyNA(d[[v]]))meth[v]<-passive[v]
  primitive<-c('rv','x','xnl','age','agenl','sex','copd','chf','ph','cad','ckd','hist','sdoc','lvent',
    'hlat','wlat','lmap','llac','zmv','sofa','plow','phigh','zfio','pzero','zpeep')
  pred<-matrix(0,ncol(d),ncol(d),dimnames=list(names(d),names(d)))
  active<-names(meth)[meth%in%c('norm','brlogreg','podds')]
  for(v in active)pred[v,setdiff(primitive,v)]<-1L
  if('sofa_ord'%in%active)pred['sofa_ord','sofa']<-0L
  # Independent latent components given the remaining clinical variables.
  # Outcome likelihood creates their substantive-model dependence.
  pred[c('plow','phigh','zfio'),c('plow','phigh','zfio')]<-0L
  pred[c('pzero','zpeep'),c('pzero','zpeep')]<-0L
  # One observed room-air endpoint cannot support a high-dimensional
  # conditional regression. A parsimonious bias-reduced model is specified
  # for this latent mixture probability; the outcome model is unchanged.
  pred['plow',]<-0L
  pred['plow',c('rv','x','age','sex','sdoc','lvent')]<-1L
  # Forty-one observed zero PEEP values cannot support every rare latent
  # mixture indicator. Use clinical ventilation/oxygenation summaries instead.
  pred['pzero',]<-0L
  pred['pzero',c('rv','x','age','sex','sdoc','lvent','sofa','lmap','pf','lmv')]<-1L
  list(data=d,method=unname(meth),methods_named=meth,predictorMatrix=pred,
    pf_counts=c(low=sum(low),high=sum(high),interior=sum(interior),missing=sum(is.na(d$pf))))
}

prag_validate<-function(z,inp) {
  d<-inp$data;s<-inp$source
  stopifnot(!anyNA(z[prag_generic]),all(is.finite(as.matrix(z[prag_generic]))))
  for(v in prag_generic) {
    ok<-!is.na(d[[v]])
    if(!isTRUE(all.equal(as.numeric(z[[v]][ok]),as.numeric(d[[v]][ok]),tolerance=1e-10)))
      stop(paste('Observed analysis value changed:',v))
  }
  stopifnot(all(z$height>100&z$height<250),all(z$weight>20&z$weight<400),
    all(z$bmi>0),all(z$map>0),all(z$peep>=0&z$peep<40),
    all(z$sofa>=0&z$sofa<=16&z$sofa==round(z$sofa)),all(2^z$lmv>0&2^z$lmv<50))
  mis<-is.na(d$pf)
  stopifnot(all(z$pf[mis]>=s$po2[mis]-1e-8),
    all(z$pf[mis]<=s$po2[mis]/.21+1e-8),all(z$pf[mis]<=713+1e-8))
  bmi_mis<-is.na(d$bmi)
  stopifnot(max(abs(z$bmi[bmi_mis]-z$weight[bmi_mis]/(z$height[bmi_mis]/100)^2))<1e-8)
  for(v in c('height','weight')) {
    obs<-!is.na(s[[v]])
    stopifnot(max(abs(z[[v]][obs]-s[[v]][obs]))<1e-8)
  }
  invisible(TRUE)
}

prag_impute<-function(inp,m=1L,numit=10L,seed=2026102701L,rjlimit=10000L) {
  prag_knots<<-inp$knots;p<-prag_prepare(inp);set.seed(seed)
  warns<-character();warn_context<-list();started<-proc.time()[3]
  ans<-withCallingHandlers(smcfcs::smcfcs(originaldata=p$data,smtype='logistic',
    smformula=prag_formula,method=p$method,predictorMatrix=p$predictorMatrix,
    m=m,numit=numit,rjlimit=rjlimit,noisy=FALSE),
    warning=function(w){
      ctx<-list(variable='unknown',iteration=NA_integer_,imputation=NA_integer_,message=conditionMessage(w))
      for(fr in rev(sys.frames()))if(exists('targetCol',envir=fr,inherits=FALSE)&&exists('originaldata',envir=fr,inherits=FALSE)) {
        ctx$variable<-names(get('originaldata',envir=fr))[get('targetCol',envir=fr)]
        if(exists('cyclenum',envir=fr,inherits=FALSE))ctx$iteration<-get('cyclenum',envir=fr)
        if(exists('imp',envir=fr,inherits=FALSE))ctx$imputation<-get('imp',envir=fr)
        break
      }
      warn_context[[length(warn_context)+1L]]<<-ctx
      warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
  invisible(lapply(ans$impDatasets,prag_validate,inp=inp))
  keep<-unique(c(prag_generic,'xnl','agenl','bminl','mapnl','pfnl','height','weight','po2','plow','phigh','zfio','pzero','zpeep'))
  ans$impDatasets<-lapply(ans$impDatasets,function(x)x[keep])
  ans$warnings<-warns;ans$warning_context<-warn_context;ans$seed<-seed;ans$numit<-numit;ans$seconds<-unname(proc.time()[3]-started)
  ans$methods<-p$methods_named;ans$predictorMatrix<-p$predictorMatrix;ans$sourceIDs<-inp$sourceIDs
  ans$metadata<-list(formula=prag_formula,knots=inp$knots,source=inp$path,
    component_model='Native SMC-FCS: bounded height/weight, three-state FiO2 with two bias-reduced hurdle models and logit interior, zero/positive bounded PEEP, bounded MV, ordinal SOFA.',
    pf_latent_counts=p$pf_counts,observed_values_unchanged=TRUE,
    caveat='Imputed PF is model-based, including flagged-source cases. The model assumes the PaO2 component is reliable and conditional latent-component independence; no imputed value is confirmed FiO2.')
  ans
}
