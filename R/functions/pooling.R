# Pooled inference helpers; this module reads no patient data.
# Rubin (1987); Barnard-Rubin (1999); D1 as in mitml 0.4-5 / mice 3.19.
# Sources: https://amices.org/mice/reference/pool.html
#          https://amices.org/mice/reference/D1.html

pi_extract <- function(fits, dfcom=NULL) {
  if(!is.list(fits) || !length(fits)) stop('fits must be a nonempty list.')
  nm <- names(coef(fits[[1L]])); p <- length(nm); m <- length(fits)
  if(!p || anyDuplicated(nm)) stop('Invalid coefficient names.')
  Q <- matrix(NA_real_,p,m,dimnames=list(nm,paste0('MI',seq_len(m))))
  U <- array(NA_real_,c(p,p,m),dimnames=list(nm,nm,paste0('MI',seq_len(m))))
  for(i in seq_len(m)) {
    qi<-coef(fits[[i]]); ui<-vcov(fits[[i]])
    if(!identical(names(qi),nm) || !identical(rownames(ui),nm) || !identical(colnames(ui),nm))
      stop('Coefficient/covariance names or order differ across imputations.')
    if(any(!is.finite(qi)) || any(!is.finite(ui))) stop('Nonfinite/aliased fit cannot be pooled.')
    Q[,i]<-qi;U[,,i]<-ui
  }
  if(is.null(dfcom)) {
    dfs<-vapply(fits,function(z) {v<-df.residual(z);if(length(v)!=1L || !is.finite(v)) Inf else v},numeric(1))
    dfcom<-min(dfs)
  }
  if(length(dfcom)!=1L || is.na(dfcom) || dfcom<=0) stop('Invalid complete-data df.')
  list(Q=Q,U=U,m=m,p=p,terms=nm,dfcom=dfcom)
}

pi_pool_scalar <- function(q,u,dfcom=Inf,conf.level=.95,null=0) {
  q<-as.numeric(q);u<-as.numeric(u);m<-length(q)
  if(!m || length(u)!=m || any(!is.finite(c(q,u))) || any(u< -1e-12)) stop('Invalid scalar estimates/variances.')
  u<-pmax(u,0);qb<-mean(q);ub<-mean(u);b<-if(m>1L)var(q) else 0
  tt<-ub+(1+1/m)*b;lambda<-if(tt>0)(1+1/m)*b/tt else 0
  if(m==1L) df<-dfcom else if(is.infinite(dfcom)) df<-if(lambda>0)(m-1)/lambda^2 else Inf else {
    tmp<-(1-lambda)*(1+dfcom)*dfcom
    df<-(m-1)*tmp/((dfcom+3)*(m-1)+lambda^2*tmp)
  }
  if(is.na(df) || df<=0) stop('Nonpositive pooled df.')
  se<-sqrt(tt);stat<-if(se>0)(qb-null)/se else if(qb==null)0 else sign(qb-null)*Inf
  p<-if(is.infinite(df))2*pnorm(-abs(stat)) else 2*pt(-abs(stat),df)
  crit<-if(is.infinite(df))qnorm(1-(1-conf.level)/2) else qt(1-(1-conf.level)/2,df)
  r<-if(ub>0)(1+1/m)*b/ub else if(b>0)Inf else 0
  fmi<-if(is.finite(r))(r+2/(df+3))/(r+1) else 1
  data.frame(m=m,estimate=qb,std.error=se,statistic=stat,df=df,p.value=p,
             lower=qb-crit*se,upper=qb+crit*se,ubar=ub,b=b,total_variance=tt,
             riv=r,lambda=lambda,fmi=fmi,mi_mcse=sqrt(b/m),row.names=NULL)
}

pi_pool_q_u <- function(Q,U,dfcom=Inf,conf.level=.95) {
  Q<-as.matrix(Q);p<-nrow(Q);m<-ncol(Q)
  if(!identical(dim(U),c(p,p,m))) stop('U must be p x p x m, matching Q.')
  if(any(!is.finite(Q)) || any(!is.finite(U))) stop('Nonfinite pooling input.')
  nm<-rownames(Q);if(is.null(nm))nm<-paste0('Q',seq_len(p))
  ub<-apply(U,c(1,2),mean);bb<-if(m>1L)cov(t(Q)) else matrix(0,p,p)
  dimnames(ub)<-dimnames(bb)<-list(nm,nm);tt<-ub+(1+1/m)*bb
  tab<-do.call(rbind,lapply(seq_len(p),function(j)cbind(term=nm[j],pi_pool_scalar(Q[j,],U[j,j,],dfcom,conf.level))))
  rownames(tab)<-NULL
  list(table=tab,estimate=setNames(rowMeans(Q),nm),Q=Q,Ubar=ub,B=bb,T=tt,m=m,dfcom=dfcom)
}

pi_pool_fits <- function(fits,dfcom=NULL,conf.level=.95) {
  e<-pi_extract(fits,dfcom);pi_pool_q_u(e$Q,e$U,e$dfcom,conf.level)
}

pi_d1_q_u <- function(Q,U,dfcom=Inf,null=0) {
  Q<-as.matrix(Q);k<-nrow(Q);m<-ncol(Q)
  if(!identical(dim(U),c(k,k,m))) stop('D1 array dimensions differ.')
  if(length(null)==1L)null<-rep(null,k)
  if(length(null)!=k)stop('D1 null dimension differs.')
  qb<-rowMeans(Q)-null;ub<-apply(U,c(1,2),mean)
  if(any(!is.finite(c(Q,U))))stop('D1 inputs are nonfinite.')
  invu<-chol2inv(chol(ub));bb<-if(m>1L)cov(t(Q)) else matrix(0,k,k)
  r<-if(m>1L)(1+1/m)*sum(diag(bb%*%invu))/k else 0
  if(r < -1e-10)stop('Invalid negative relative variance.')
  r<-max(r,0);Fval<-as.numeric(crossprod(qb,invu%*%qb))/(k*(1+r))
  if(m==1L) v<-dfcom else {
    tt<-k*(m-1)
    if(is.finite(dfcom) && tt>4) {
      # Reiter's finite complete-data adjustment, matching mitml:::.D1.
      a<-r*tt/(tt-2);vstar<-(dfcom+1)/(dfcom+3)*dfcom
      c0<-1/(tt-4);c1<-vstar-2*(1+a);c2<-vstar-4*(1+a)
      z<-1/c2+c0*a^2*c1/((1+a)^2*c2)+
        c0*(8*a^2*c1/((1+a)*c2^2)+4*a^2/((1+a)*c2))+
        c0*(4*a^2/(c2*c1)+16*a^2*c1/c2^3)+c0*8*a^2/c2^2
      v<-4+1/z
    } else if(r==0) v<-Inf else if(tt>4) v<-4+(tt-4)*(1+(1-2/tt)/r)^2 else
      v<-tt*(1+1/k)*(1+1/r)^2/2
  }
  if(!is.finite(Fval) || is.na(v) || v<=0)stop('D1 computation invalid; inspect inputs.')
  data.frame(statistic=Fval,df1=k,df2=v,dfcom=dfcom,p.value=pf(Fval,k,v,lower.tail=FALSE),riv=r,m=m,
             method='D1 multivariate Wald',row.names=NULL)
}

pi_d1 <- function(fits,terms=c('rv:x','rv:xnl'),dfcom=NULL) {
  e<-pi_extract(fits,dfcom);j<-match(terms,e$terms)
  if(anyNA(j) || anyDuplicated(j))stop('D1 terms missing or duplicate.')
  z<-pi_d1_q_u(e$Q[j,,drop=FALSE],e$U[j,j,,drop=FALSE],e$dfcom)
  z$terms<-paste(terms,collapse=' | ');z
}

pi_model_matrix <- function(fit,newdata,design_fun=identity) {
  if(length(attr(terms(fit),'offset')))stop('Offset models are not implemented.')
  z<-design_fun(newdata)
  X<-if(is.matrix(z))z else model.matrix(delete.response(terms(fit)),data=z,contrasts.arg=fit$contrasts,xlev=fit$xlevels)
  nm<-names(coef(fit))
  if(!all(nm%in%colnames(X)))stop('New design lacks fitted terms.')
  X<-X[,nm,drop=FALSE]
  if(any(!is.finite(X)))stop('Nonfinite new design.')
  X
}

pi_standardized <- function(fit,data,design_fun=identity,rv_var='rv',x_var='x',levels=c(40,60),weights=NULL) {
  if(!identical(fit$family$link,'logit'))stop('Standardization currently supports logit fits.')
  if(length(levels)!=2L)stop('Specify two PaCO2 levels.')
  n<-nrow(data);if(!n)stop('Empty standardization target.')
  if(is.null(weights))weights<-rep(1/n,n) else {
    if(length(weights)!=n || any(!is.finite(weights)) || any(weights<0) || sum(weights)<=0)stop('Invalid target weights.')
    weights<-weights/sum(weights)
  }
  mu<-numeric(4);G<-matrix(NA_real_,4,length(coef(fit)))
  nm<-c('r0x40','r0x60','r1x40','r1x60')
  nm<-paste0('r',rep(0:1,each=2),'x',rep(levels,2))
  scenarios<-expand.grid(x=levels,rv=0:1)
  for(j in 1:4) {
    z<-data;z[[rv_var]]<-scenarios$rv[j];z[[x_var]]<-scenarios$x[j]
    X<-pi_model_matrix(fit,z,design_fun);pr<-plogis(drop(X%*%coef(fit)))
    mu[j]<-sum(weights*pr);G[j,]<-colSums(X*(weights*pr*(1-pr)))
  }
  L<-rbind(diag(4),c(-1,1,0,0),c(0,0,-1,1),c(-1,0,1,0),c(0,-1,0,1),c(1,-1,-1,1))
  rn<-c(nm,'change_r0','change_r1',paste0('rv_rd_x',levels),'delta')
  gg<-L%*%G;rownames(gg)<-rn;colnames(gg)<-names(coef(fit))
  v<-gg%*%vcov(fit)%*%t(gg);q<-setNames(drop(L%*%mu),rn)
  list(q=q,U=v,gradient=gg,risks=setNames(mu,nm),target_n=n,
       variance_scope='Conditional on the completed empirical target covariate distribution; beta covariance delta method.')
}

pi_pool_standardized <- function(fits,datasets,design_fun=identity,rv_var='rv',x_var='x',levels=c(40,60),dfcom=NULL,conf.level=.95) {
  if(length(fits)!=length(datasets))stop('Fit/data MI counts differ.')
  a<-lapply(seq_along(fits),function(i)pi_standardized(fits[[i]],datasets[[i]],design_fun,rv_var,x_var,levels))
  Q<-do.call(cbind,lapply(a,'[[','q'));p<-nrow(Q);m<-ncol(Q)
  U<-array(unlist(lapply(a,'[[','U')),c(p,p,m),dimnames=list(rownames(Q),rownames(Q),NULL))
  e<-pi_extract(fits,dfcom);out<-pi_pool_q_u(Q,U,e$dfcom,conf.level)
  out$per_imputation<-a;out$variance_scope<-a[[1]]$variance_scope;out
}

pi_or_curve <- function(fits,datasets,xgrid,design_fun=identity,reference=40,rv_var='rv',x_var='x',dfcom=NULL,conf.level=.95) {
  if(length(fits)!=length(datasets))stop('Fit/data MI counts differ.')
  e<-pi_extract(fits,dfcom);rows<-list();jj<-0L
  for(r in 0:1) for(x in xgrid) {
    qi<-ui<-numeric(length(fits))
    for(i in seq_along(fits)) {
      z<-datasets[[i]][rep(1L,2),,drop=FALSE];z[[rv_var]]<-r;z[[x_var]]<-c(x,reference)
      X<-pi_model_matrix(fits[[i]],z,design_fun);cc<-X[1,]-X[2,]
      qi[i]<-sum(cc*coef(fits[[i]]));ui[i]<-drop(crossprod(cc,vcov(fits[[i]])%*%cc))
    }
    jj<-jj+1L;po<-pi_pool_scalar(qi,ui,e$dfcom,conf.level)
    rows[[jj]]<-cbind(rv=r,pco2=x,reference=reference,po,OR=exp(po$estimate),OR_lower=exp(po$lower),OR_upper=exp(po$upper))
  }
  out<-do.call(rbind,rows);rownames(out)<-NULL;out
}

pi_pool_contrast <- function(fits,contrast,dfcom=NULL,conf.level=.95) {
  e<-pi_extract(fits,dfcom);m<-length(fits);q<-u<-numeric(m)
  for(i in seq_len(m)) {
    cc<-if(is.list(contrast))contrast[[i]] else contrast
    if(is.null(names(cc))) {if(length(cc)!=e$p)stop('Unnamed contrast size differs.')} else {
      if(!all(e$terms%in%names(cc)))stop('Named contrast lacks model coefficients.')
      cc<-cc[e$terms]
    }
    q[i]<-sum(cc*e$Q[,i]);u[i]<-drop(crossprod(cc,e$U[,,i]%*%cc))
  }
  pi_pool_scalar(q,u,e$dfcom,conf.level)
}
