####################
## functions extracted from R/stats.R.
## wrapper for R's t distribution specialized for 
##   two group difference-of-mean. Cohen's d.
##   same sample size each group. normal populations w/ unit variance
## power_d2ht (the last function in the file) computes power
##
## Kenny & Judd scale factor
q_=function(n,sd.het) sqrt((2/n)/((2/n)+sd.het^2));

## significance boundary for Cohen's d
d_htcrit=function(n,sd.het,sig.level=0.05) htpval2d(n,sd.het,pval=sig.level)

## pval to Cohen's d
htpval2d=function(n,sd.het,pval) q_d2ht(n,sd.het=sd.het,p=pval/2,lower.tail=F)

## noncentrality parameter
ncp_ht=function(n,d.het,sd.het) sqrt(n/(2+(n*sd.het^2)))*d.het;

## het t statistic to/from Cohen's d, conventional t, pval
t2ht=function(n,sd.het,t) t*q_(n,sd.het);
ht2t=function(n,sd.het,ht) ht/q_(n,sd.het);
d2ht=function(n,sd.het,d) t2ht(n,sd.het,d2t(n,d));
ht2d=function(n,sd.het,ht) t2d(n,ht2t(n,sd.het,ht));
d2htpval=function(n,sd.het,d) t2pval(n,d2ht(n,sd.het,d))
htpval2d=function(n,sd.het,pval) q_d2ht(n,sd.het=sd.het,p=pval/2,lower.tail=F)

## conventional t statistic to/from Cohen's d
t2d=function(n,t) t*sqrt(2/n);
d2t=function(n,d) d*sqrt(n/2);

## quantile function
q_d2ht=function(n,d.het=NULL,sd.het,p,lower.tail=TRUE) {
  df=2*(n-1);
  if (!is.null(d.het))
    ht=suppressWarnings(qt(p,df=df,ncp=ncp_ht(n,d.het,sd.het),lower.tail=lower.tail))
  else ht=qt(p,df=df,lower.tail=lower.tail)
  ht2d(n,sd.het,ht);
}
## cumulative distribution function
p_d2ht=function(n,d.het=NULL,sd.het,d,lower.tail=TRUE) {
  df=2*(n-1);
  ht=d2ht(n,sd.het,d);
  if (!is.null(d.het))
    suppressWarnings(pt(ht,df=df,ncp=ncp_ht(n,d.het,sd.het),lower.tail=lower.tail))
  else pt(ht,df=df,lower.tail=lower.tail)
}
## power. strict=F only considers positive tail. default for power.t.test
power_d2ht=function(n,d,sd.het,d.crit,sig.level=0.05,strict=F) {
  if (missing(d.crit)) d.crit=d_htcrit(n=n,sd.het=sd.het,sig.level=sig.level);
  p_d2ht(n=n,d.het=d,sd.het=sd.het,d=d.crit,lower.tail=F)+
    ifelse(strict,p_d2ht(n=n,d.het=d,sd.het=sd.het,d=-d.crit,lower.tail=T),0)
}
