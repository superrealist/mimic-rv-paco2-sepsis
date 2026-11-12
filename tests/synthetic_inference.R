# Mathematical checks using synthetic values only; no patient input is read.
root<-Sys.getenv('MIMIC_CODE_ROOT','.')
source(file.path(root,'R/functions/pooling.R'))
source(file.path(root,'R/functions/risk_contrasts.R'))
q<-c(-.03,.01,.02);u<-c(.001,.002,.0015)
p<-pi_pool_scalar(q,u,dfcom=100)
stopifnot(abs(p$estimate-mean(q))<1e-14,
 abs(p$total_variance-(mean(u)+(1+1/length(q))*var(q)))<1e-14)
d<-data.frame(y=rep(0:1,4),rv=rep(0:1,4),x=seq(25,60,5),age=seq(40,75,5))
d$xnl<-fixed_rcs(d$x)
form<-'y~rv*(x+xnl)+age'
nm<-colnames(model.matrix(as.formula(form),d))
b<-setNames(c(-2,.25,.01,.02,.015,-.003,-.001),nm)
V<-diag(length(b))*.0001;dimnames(V)<-list(nm,nm)
sc<-expand.grid(pco2=seq(25,65,.1),rv=0:1)
job<-list(data=d,b=b,V=V,sc=sc,formula=form,encoding='mimic_raw')
z<-risk_one_completion(job)
step<-1e-6
for(k in seq_along(b)){
 plus<-minus<-job;plus$b[k]<-b[k]+step;minus$b[k]<-b[k]-step
 numeric_gradient<-(risk_one_completion(plus)$q-risk_one_completion(minus)$q)/(2*step)
 stopifnot(max(abs(numeric_gradient-z$G[,k]))<1e-6)
}
co<-contrast_design(sc)
stopifnot(all(co$relative_L[sc$pco2==40,]==0),sum(co$sliding_meta$metric=='interaction')==301)
g<-co$fixed_L%*%z$G
v<-rowSums((g%*%V)*g)
full<-co$fixed_L%*%(z$G%*%V%*%t(z$G))%*%t(co$fixed_L)
stopifnot(max(abs(v-diag(full)))<1e-12)
cat('PASS: scalar pooling, finite-difference gradients, shared-reference and joint covariance checks.\n')
