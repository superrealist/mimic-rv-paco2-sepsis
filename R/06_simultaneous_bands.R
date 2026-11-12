source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
O<-file.path(OUT,'contrasts')
nsim<-200000L;chunk<-5000L;confidence<-.95
summaries<-list()
for(tag in 'mimic'){
 z<-readRDS(file.path(O,'private',paste0(tag,'_inference.rds')))
 T<-z$interaction$T;T<-(T+t(T))/2;se<-sqrt(diag(T));ok<-se>1e-14
 stopifnot(all(ok),nrow(T)==301)
 R<-cov2cor(T[ok,ok,drop=FALSE]);ee<-eigen((R+t(R))/2,symmetric=TRUE)
 negtol<-1e-10*max(ee$values);stopifnot(min(ee$values)>=-negtol)
 positive_threshold<-1e-12*max(ee$values);keep<-ee$values>positive_threshold
 F<-sweep(ee$vectors[,keep,drop=FALSE],2,sqrt(ee$values[keep]),'*')
 approxR<-tcrossprod(F);approxerror<-max(abs(approxR-R));stopifnot(approxerror<1e-8)
 seed<-202611091L;set.seed(seed)
 maxima<-numeric(nsim)
 for(a in seq.int(1L,nsim,by=chunk)){
  ii<-a:min(nsim,a+chunk-1L)
  draws<-matrix(rnorm(length(ii)*ncol(F)),nrow=length(ii))%*%t(F)
  maxima[ii]<-matrixStats::rowMaxs(abs(draws))
 }
 critical<-unname(quantile(maxima,confidence,type=8))
 orderstat<-sort(maxima);ranks<-pmin(nsim,pmax(1L,qbinom(c(.025,.975),nsim,confidence)+1L))
 critical_mc95<-orderstat[ranks]
 s<-read.csv(file.path(O,'tables',paste0(tag,'_sliding.csv')))
 ii<-which(s$metric=='interaction');stopifnot(max(abs(s$se[ii]-se))<1e-10)
 s$simultaneous_lower<-NA_real_;s$simultaneous_upper<-NA_real_
 s$simultaneous_lower[ii]<-s$estimate[ii]-critical*se
 s$simultaneous_upper[ii]<-s$estimate[ii]+critical*se
 write.csv(s,file.path(O,'tables',paste0(tag,'_sliding.csv')),row.names=FALSE)
 band<-s[ii,]
 write.csv(band,file.path(O,'tables',paste0(tag,'_interaction_band.csv')),row.names=FALSE)
 saveRDS(list(maxima=maxima,critical=critical,seed=seed,correlation=R,factor=F,metadata=z$interaction$metadata),file.path(O,'private',paste0(tag,'_band_simulation.rds')),compress=FALSE)
 summaries[[tag]]<-list(cohort=z$cohort,confidence=confidence,critical=critical,critical_simulation_MC95=critical_mc95,nsim=nsim,seed=seed,
  grid_start=25,grid_end=55,grid_step=.1,grid_points=301,increment=10,
  zero_variance_grid_points=sum(!ok),minimum_correlation_eigenvalue=min(ee$values),negative_eigenvalues=sum(ee$values<0),
  negative_eigenvalue_tolerance=negtol,positive_eigenvalue_retention_threshold=positive_threshold,retained_rank=sum(keep),
  correlation_reconstruction_max_error=approxerror,method='Gaussian max-absolute-standardized-error simulation from full Rubin total covariance',
  asymptotic=TRUE,finite_sample_or_finite_MI_exact_coverage=FALSE,
  coverage_family='The 301 prespecified I(x) points within this cohort only; not a joint band across cohorts or other contrast families.',
  continuum_coverage='No claim of exact coverage between grid points.',
  pointwise_positive_points=sum(band$lower>0),pointwise_negative_points=sum(band$upper<0),
  simultaneous_positive_points=sum(band$simultaneous_lower>0),simultaneous_negative_points=sum(band$simultaneous_upper<0))
 cat(tag,'critical=',critical,'rank=',sum(keep),'simultaneous excludes zero=',sum(band$simultaneous_lower>0|band$simultaneous_upper<0),'\n')
}
jsonlite::write_json(summaries,file.path(O,'quality/band_simulation.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
for(kind in c('risks','relative40','sliding','fixed_contrasts','tests','interaction_band')){
 d<-do.call(rbind,lapply('mimic',function(tag)read.csv(file.path(O,'tables',paste0(tag,'_',kind,'.csv')))))
 write.csv(d,file.path(O,'tables',paste0(kind,'.csv')),row.names=FALSE)
}
