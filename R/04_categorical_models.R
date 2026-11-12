source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
O<-file.path(OUT,'category')
make_dirs(O,c('tables'))
source(FUN('pooling.R'))
base<-read.csv(INPUT,check.names=FALSE)
spec<-as.formula(FORM_CATEGORY)
sc <- data.frame(rv=rep(0:1,each=3),category=rep(c('low','middle','high'),2))
ci <- c(1,3,4,6); ref <- c(2,2,5,5)
L <- diag(6)[ci,,drop=FALSE]-diag(6)[ref,,drop=FALSE]
for(endpoint in 'death28') {
  path<-file.path(OUT,'category/imputations.rds')
  z <- readRDS(path)
  raw <- base[match(as.character(z$sourceIDs),as.character(base$subject_id)),]
  stopifnot(!anyNA(raw$subject_id),!anyDuplicated(raw$subject_id))
  truth<-raw$death28
  binary <- endpoint %in% c('death28','lib7')
  res<-lapply(z$impDatasets,function(d){
    stopifnot(isTRUE(all.equal(as.numeric(d$x),as.numeric(raw$pco2))),all(d$rv==raw$rv_binary),max(abs(d$y-truth))<1e-8)
    d$co2lo<-as.integer(d$x<35); d$co2hi<-as.integer(d$x>45)
    f<-if(binary) glm(spec,binomial(),d,control=glm.control(maxit=500)) else lm(spec,d)
    stopifnot(nobs(f)==nrow(raw),f$rank==length(coef(f)),all(is.finite(coef(f))), !binary||f$converged)
    V<-if(binary) vcov(f) else sandwich::vcovHC(f,type='HC3')
    b<-coef(f); mu<-numeric(6);G<-H<-matrix(0,6,length(b))
    for(k in 1:6){
      nd<-d;nd$rv<-sc$rv[k];nd$co2lo<-as.integer(sc$category[k]=='low');nd$co2hi<-as.integer(sc$category[k]=='high')
      X<-pi_model_matrix(f,nd);prob<-if(binary)plogis(X%*%b) else drop(X%*%b)
      mu[k]<-mean(prob);G[k,]<-if(binary)colMeans(X*as.numeric(prob*(1-prob))) else colMeans(X)
      H[k,]<-X[1,]
    }
    D<-L%*%G;C<-L%*%H
    list(q=drop(L%*%mu),u=diag(D%*%V%*%t(D)),conditional=drop(C%*%b),
         cond_var=diag(C%*%V%*%t(C)),df=df.residual(f))
  })
  pool<-function(a,b){do.call(rbind,lapply(1:4,function(i)pi_pool_scalar(
    vapply(res,function(v)v[[a]][i],numeric(1)),vapply(res,function(v)v[[b]][i],numeric(1)),res[[1]]$df)))}
  pp<-pool('q','u')
  public_endpoint<-endpoint
  ftab<-cbind(outcome=public_endpoint,sc[ci,],n=nrow(raw),pp[c('estimate','lower','upper','p.value')])
  target<-file.path(O,'tables',paste0('risk_contrasts_',public_endpoint,'.csv'))
  write.csv(ftab,target,row.names=FALSE)
  if(endpoint=='death28'){
    cp<-pool('conditional','cond_var')
    tab<-sc;cc<-ifelse(raw$pco2<35,'low',ifelse(raw$pco2>45,'high','middle'))
    tab$n<-sapply(1:6,function(i)sum(raw$rv_binary==sc$rv[i]&cc==sc$category[i]))
    tab$events<-sapply(1:6,function(i)sum(raw$death28[raw$rv_binary==sc$rv[i]&cc==sc$category[i]]))
    tab$or<-1;tab$lower<-tab$upper<-tab$p<-NA_real_
    tab$or[ci]<-exp(cp$estimate);tab$lower[ci]<-exp(cp$lower);tab$upper[ci]<-exp(cp$upper);tab$p[ci]<-cp$p.value
    stopifnot(sum(tab$n)==2051,sum(tab$events)==998)
    write.csv(tab,file.path(O,'tables/categories_within.csv'),row.names=FALSE)
  }
}
