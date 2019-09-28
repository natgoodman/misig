#################################################################################
##
## Author:  Nat Goodman
## Created: 19-01-01
##          from repwr/stats.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Statistical functions for repwr.R
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- Statistical Functions for d2t Distribution ----
## t.test related functions for one and two samples
## two sample functions all assume equal samples size and variance

## My formula for pooled_sd, though independent of n, is correct. It's a simplification
##   of the standard formula for n0=n1
##   standard formula: sqrt(((n-1)*sd0^2+(n-1)*sd1^2)/(n+n-2));
pooled_sd=function(sd0,sd1) sqrt((sd0^2+sd1^2)/2);

## t statistic (not used)
tstat=function(n,mean0,mean1,sd0,sd1) (mean0-mean1)/sqrt((sd0^2+sd1^2)/n);
## NG 19-01-21: Geez Louise - how did I not simplify t2d & d2t from the gitgo???
## t2d=function(n,t) t*sqrt((2*n)/n^2)
## d2t=function(n,d) d*sqrt(n^2/(2*n))

## t statistic to Cohen's d & pval
t2d=function(n,t) t*sqrt(2/n);
t2pval=function(n,t) 2*pt(-abs(t),df=2*(n-1))
## Cohen's d to t statistic & pval
d2t=function(n,d) d*sqrt(n/2);
d2pval=function(n,d) t2pval(n,d2t(n,d))
## pval to t statistic & Cohen's d
pval2t=function(n,pval) qt(pval/2,df=2*(n-1),lower.tail=F)
pval2d=function(n,pval) q_d2t(n,p=pval/2,lower.tail=F)
## confidence interval of d.raw
ci_draw=function(n,d.raw,sd,conf.level=0.95) {
  p0=(1-conf.level)/2; p1=1-p0;
  tstar=sd/sqrt(n/2)*qt(p1,df=2*n-2);
  setNames(c(d.raw-tstar,d.raw+tstar),paste(sep='',100*c(p0,p1),'%'));
}
## significance boundary for Cohen's d
d_crit=d_sig=function(n,sig.level=param(sig.level)) pval2d(n,pval=sig.level)

## probability functions for t-distribution of d
ncp=function(n,d) sqrt(n/2)*d
## NG 19-02-14: I messed up density here and in repwr
##   need to scale by d2t (sqrt(n/2)) to account for difference in x-density
d_d2t=function(n,d,d0=NULL) {
  df=2*(n-1);
  t=d2t(n,d);
  sqrt(n/2)*
    if (!is.null(d0)) suppressWarnings(dt(t,df=df,ncp=ncp(n,d0)))
    else dt(t,df=df)
}
p_d2t=function(n,d,d0=NULL,lower.tail=TRUE) {
  df=2*(n-1);
  t=d2t(n,d);
  if (!is.null(d0)) suppressWarnings(pt(t,df=df,ncp=ncp(n,d0),lower.tail=lower.tail))
    else pt(t,df=df,lower.tail=lower.tail)
}
q_d2t=function(n,p,d0=NULL,lower.tail=TRUE) {
  df=2*(n-1);
  if (!is.null(d0)) t=suppressWarnings(qt(p,df=df,ncp=ncp(n,d0),lower.tail=lower.tail))
    else t=qt(p,df=df,lower.tail=lower.tail)
  t2d(n,t);
}
r_d2t=function(m,n,d0=NULL) {
  df=2*(n-1);
  if (!is.null(d0)) t=suppressWarnings(rt(m,df=df,ncp=ncp(n,d0))) else t=rt(m,df=df);
  t2d(n,t)
}
## mean and sd for t-distribution of d 
mean_d2t=function(n,d0=NULL) {
  df=2*(n-1);
  if (!is.null(d0)) {
    ncp=ncp(n,d0);
    ## NG 18-02-07. gamma blows up when n>100 or so. use lgamma instead
    ## theo.mean=sqrt(df/2)*ncp*gamma((df-1)/2)/gamma(df/2)
    theo.mean=sqrt(df/2)*ncp*exp(lgamma((df-1)/2)-lgamma(df/2))
    t2d(n,theo.mean)
  } else 0;
}
sd_d2t=function(n,d0=NULL) {
  df=2*(n-1);
  if (!is.null(d0)) {
    ncp=ncp(n,d0);
    ## NG 18-02-07. gamma blows up when n>100 or so. use lgamma instead
    ## theo.mean=sqrt(df/2)*ncp*gamma((df-1)/2)/gamma(df/2)
    theo.mean=sqrt(df/2)*ncp*exp(lgamma((df-1)/2)-lgamma(df/2))
    theo.var=(1+ncp^2)*df/(df-2)-(theo.mean^2)
    theo.sd=sqrt(theo.var)
    t2d(n,theo.sd)
  } else
    (sqrt(2*n)/n)*sdt(2*(n-1));
}
## sd of (central) t distribution
sdt=function(df) sqrt(df/(df-2))

## confidence and prediction intervals for t-distribution of d
## my adaptation of confidence interval function from
## http://urisohn.com/sohn_files/BlogAppendix/Colada20.ConfidenceIntervalsForD.R
## can possibly make it a bit faster by unwrapping p_d2t
ci_d2t=function(n,d,simplify=T,conf.level=param(conf.level)) {
  ci=ci_d2t_(n,d,conf.level);
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
ci_d2t_=Vectorize(function(n,d,conf.level) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=suppressWarnings(
    uniroot(function(d0) p_d2t(n,d,d0,lower.tail=F)-p0,interval=c(-10,10))$root);
  hi=suppressWarnings(
    uniroot(function(d0) p_d2t(n,d,d0,lower.tail=F)-p1,interval=c(-10,10))$root);
  c(lo,hi);
})
## my adaptation of prediction interval function from
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5028066/ and predictionInterval pacakge
pi_d2t=function(n1,n2,d,ci1=NULL,ci2=NULL,pred.level=0.95) {
  if (is.null(ci1)) ci1=ci_d2t(n1,d,pred.level);
  if (is.null(ci2)) ci2=ci_d2t(n2,d,pred.level);
  l1=ci1[1]; u1=ci1[2];
  l2=ci2[1]; u2=ci2[2];
  c(d-sqrt((d-l1)^2+(u2-d)^2),d+sqrt((d-l2)^2+(u1-d)^2));
}
## mean signficant effect size
## TODO: not yet used in running code
meand_d2t=Vectorize(
  function(n,d0,d.crit) {
    if (missing (d.crit)) d.crit=d_crit(n);
    meand=integrate(function(d) d_d2t(n=n,d=d,d0=d0)*d,lower=d.crit,upper=Inf)$value / p_d2t(n=n,d0=d0,d=d.crit,lower.tail=F);
  })
## power - wrapper for power.t.test for stylistic consistency with power_d2ht 
power_d2t=function(n,d,sig.level=param(sig.level),strict=F)
  power.t.test(n=n,delta=d,sig.level=sig.level,strict=strict)$power

## ---- Statistical Functions for d2ht Distribution ----
## d2ht is t-distribution with effect size heterogeneity
## I have two implementationss:
##  1) scaled d2t based on Kenny & Judd paper
##  2) my own starting from first principles
## running code uses KJ. I used mine during development to check KJ

## Kenny & Judd implementation scales t by q=sqrt((2/n)/((2/n)+sd.het^2))
q_=function(n,sd.het) sqrt((2/n)/((2/n)+sd.het^2));

t2ht=function(n,sd.het,t) t*q_(n,sd.het);
ht2t=function(n,sd.het,ht) ht/q_(n,sd.het);
d2ht=function(n,sd.het,d) t2ht(n,sd.het,d2t(n,d));
ht2d=function(n,sd.het,ht) t2d(n,ht2t(n,sd.het,ht));

d2htpval=function(n,sd.het,d) t2pval(n,d2ht(n,sd.het,d))
htpval2d=function(n,sd.het,pval) q_d2ht(n,sd.het=sd.het,p=pval/2,lower.tail=F)

## significance boundary for Cohen's d
d_htcrit=function(n,sd.het,sig.level=param(sig.level)) htpval2d(n,sd.het,pval=sig.level)

## need to scale by q_*d2t (sqrt(n/2)) to account for t scaling and difference in x-density
ncp_ht=function(n,d.het,sd.het) sqrt(n/(2+(n*sd.het^2)))*d.het;
d_d2ht=function(n,d.het=NULL,sd.het,d) {
  df=2*(n-1);
  ht=d2ht(n,sd.het,d);
  sqrt(n/(2+(n*sd.het^2)))*              # simplification of q_(n,sd.het)*sqrt(n/2)*
    if (!is.null(d.het)) suppressWarnings(dt(ht,df=df,ncp=ncp_ht(n,d.het,sd.het)))
    else dt(ht,df=df)
}
p_d2ht=function(n,d.het=NULL,sd.het,d,lower.tail=TRUE) {
  df=2*(n-1);
  ht=d2ht(n,sd.het,d);
  if (!is.null(d.het))
    suppressWarnings(pt(ht,df=df,ncp=ncp_ht(n,d.het,sd.het),lower.tail=lower.tail))
  else pt(ht,df=df,lower.tail=lower.tail)
}
q_d2ht=function(n,d.het=NULL,sd.het,p,lower.tail=TRUE) {
  df=2*(n-1);
  if (!is.null(d.het))
    ht=suppressWarnings(qt(p,df=df,ncp=ncp_ht(n,d.het,sd.het),lower.tail=lower.tail))
  else ht=qt(p,df=df,lower.tail=lower.tail)
  ht2d(n,sd.het,ht);
}
r_d2ht=function(n,d.het=NULL,sd.het,m) {
  df=2*(n-1);
  if (!is.null(d.het)) ht=suppressWarnings(rt(m,df=df,ncp=ncp_ht(n,d.het,sd.het)))
  else ht=rt(m,df=df);
  ht2d(n,sd.het,ht)
}

## my adaptation of confidence interval function from
## http://urisohn.com/sohn_files/BlogAppendix/Colada20.ConfidenceIntervalsForD.R
## can possibly make it a bit faster by unwrapping p_d2ht
ci_d2ht=function(n,sd.het=0,d,simplify=T,conf.level=param(conf.level)) {
  ci=ci_d2ht_(n,sd.het,d,conf.level);
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
ci_d2ht_=Vectorize(function(n,sd.het,d,conf.level) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=suppressWarnings(
    uniroot(function(d.het)
      p_d2ht(n,d.het=d.het,sd.het=sd.het,d=d,lower.tail=F)-p0,interval=c(-10,10))$root);
  hi=suppressWarnings(
    uniroot(function(d.het)
      p_d2ht(n,d.het=d.het,sd.het=sd.het,d=d,lower.tail=F)-p1,interval=c(-10,10))$root);
  c(lo,hi);
})
## mean signficant effect size
meand_d2ht=Vectorize(
  function(n,d.het,sd.het,d.crit) {
    if (missing(d.crit)) d.crit=d_htcrit(n=n,sd.het=sd.het);
    integrate(function(d)
      d_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=d)*d,lower=d.crit,upper=Inf)$value / p_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=d.crit,lower.tail=F);
  })
## power. strict=F only considers positive tail. default for power.t.test
power_d2ht=function(n,d,sd.het,d.crit,sig.level=param(sig.level),strict=F) {
  if (missing(d.crit)) d.crit=d_htcrit(n=n,sd.het=sd.het,sig.level=sig.level);
  p_d2ht(n=n,d.het=d,sd.het=sd.het,d=d.crit,lower.tail=F)+
    ifelse(strict,p_d2ht(n=n,d.het=d,sd.het=sd.het,d=-d.crit,lower.tail=T),0)
}

########################################
## my from-scratch implementations of d2ht
## probability density
d_d2ht_ng=
  Vectorize(function(n,d.het=0,sd.het=0.2,d) {
    integrate(function(d0,d.het,sd.het,n,d) dnorm(d0,mean=d.het,sd=sd.het)*d_d2t(n=n,d=d,d0=d0),
              d.het=d.het,sd.het=sd.het,n=n,d=d,lower=Inf,upper=Inf)$value;
  },vectorize.args='d')
## remaining functions are quick hacks
## TODO: do it right! when I figure out how...
## CAUTION: interpolation code duplicated in dat.R. FIX this!!
## interpolation function used by p_d2ht,q_d2ht,r_d2ht
## SLOW!! calls integrate dlen times to get successive values of integral
## TODO: there must be a direct way to generate successive integrals
ht_interp=function(n,d.het=NULL,sd.het=0.2,dlim=c(-3,3),dlen=1e2) {
  if (is.null(d.het)) d.het=0;
  dlim=dlim+d.het;
  dx=seq(min(dlim),max(dlim),len=dlen);
  cum=ht_interp_(dx,n=n,d.het=d.het,sd.het=sd.het);
  approxfun(dx,cum);
}
ht_interp_=Vectorize(function(d,...)
  integrate(function(d) d_d2ht_ng(d=d,...),lower=-Inf,upper=d)$value);
  
p_d2ht_ng=function(ht.interp,d,...) {
  if (missing(ht.interp)) ht.interp=ht_interp(...);
  ht.interp(d);
}
q_d2ht_ng=function(ht.interp,p,...) {
  if (missing(ht.interp)) ht.interp=ht_interp(...);
  sapply(p,function(p) uniroot(function(d) ht.interp(d)-p,interval=c(-2,2))$root);
}
## random number generator.
## basic method from https://en.wikibooks.org/wiki/R_Programming/Random_Number_Generation
r_d2ht_ng=function(ht.interp,m,...) q_d2ht_ng(ht.interp,p=runif(m),...)

## significance boundary
d_htcrit_ng=function(sig.level=param(sig.level),...) q_d2ht_ng(ht.interp,q=1-(sig.level/2),...)

